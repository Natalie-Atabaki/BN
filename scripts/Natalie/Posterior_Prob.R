# ================================ #
# Bayesian Network Analysis (Simulated)
# Trait1 ~ 34 Proteins
# ================================ #

# -------- Load Required Libraries --------
library(bnlearn)
library(gRain)
library(dplyr)
library(ggplot2)

# -------- Simulate Data --------
set.seed(123)

n_samples <- 964
n_traits <- 19
n_proteins <- 34

# Simulated clinical data: ID + 19 traits
Clinical <- data.frame(
  ID = sprintf("ID_%04d", 1:n_samples),
  matrix(rnorm(n_samples * n_traits), ncol = n_traits)
)
colnames(Clinical)[-1] <- paste0("Trait", 1:n_traits)

# Simulated protein data: ID + 34 proteins
Protein <- data.frame(
  ID = Clinical$ID,
  matrix(rnorm(n_samples * n_proteins), ncol = n_proteins)
)
colnames(Protein)[-1] <- paste0("Prot", 1:n_proteins)

# -------- Merge & Discretize --------
ClinicalProtein <- merge(Clinical, Protein, by = "ID")

# Discretize all numeric variables (excluding ID)
Discretized <- discretize(ClinicalProtein[,-1],
                          method = "hartemink",
                          breaks = 3, ibreaks = 9, idisc = "quantile")

# Label levels: LOW / AVG / HIGH
for (col in names(Discretized)) {
  levels(Discretized[[col]]) <- c("LOW", "AVG", "HIGH")
}

# -------- Fit Bayesian Network --------
dag_model <- hc(Discretized)
fitted_bn <- bn.fit(dag_model, data = Discretized, method = "bayes", iss = 5)

# Convert to junction tree
jtree <- compile(as.grain(fitted_bn))

# Prior for Trait1
cat("Unconditional probability of Trait1:\n")
print(querygrain(jtree, nodes = "Trait1")$Trait1)

# -------- Conditional Inference --------
# Assume: first 27 proteins positively associated, last 7 negatively

Proteins_High <- paste0("Prot", 1:27)
Proteins_Low  <- paste0("Prot", 28:34)

# P(Trait1 = HIGH | ProtX = HIGH)
results_high <- sapply(Proteins_High, function(protein) {
  particles <- cpdist(fitted_bn, nodes = "Trait1", evidence = (Discretized[[protein]] == "HIGH"))
  prop.table(table(particles))[["HIGH"]]
})

# P(Trait1 = HIGH | ProtX = LOW)
results_low <- sapply(Proteins_Low, function(protein) {
  particles <- cpdist(fitted_bn, nodes = "Trait1", evidence = (Discretized[[protein]] == "LOW"))
  prop.table(table(particles))[["HIGH"]]
})

# -------- Format Results --------
df_high <- data.frame(
  Protein = names(results_high),
  Probability = unname(results_high),
  BorderColor = "lightcoral"
)

df_low <- data.frame(
  Protein = names(results_low),
  Probability = unname(results_low),
  BorderColor = "deepskyblue3"
)

ranked_results <- bind_rows(df_high, df_low) %>%
  arrange(desc(Probability)) %>%
  mutate(Protein = factor(Protein, levels = Protein))

# -------- Visualization --------
ggplot(ranked_results, aes(x = Protein, y = Probability, color = BorderColor)) +
  geom_bar(stat = "identity", fill = "white", size = 0.3, position = position_dodge(width = 1)) +
  geom_text(aes(label = round(Probability, 2)), vjust = -0.3, size = 2.5) +
  coord_polar() +
  theme_minimal() +
  labs(title = "P(Trait1 = HIGH) by Protein Level") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 90, size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(size = 10, hjust = 0.5)
  ) +
  scale_color_manual(values = c("lightcoral", "deepskyblue3"))
