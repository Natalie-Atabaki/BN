# ----------------------------- Load Required Libraries -----------------------------
library(dplyr)
library(bnlearn)
library(parallel)
library(openxlsx)

set.seed(123)  # For reproducibility

# ----------------------------- Generate Synthetic Data -----------------------------

# Simulate clinical data: 964 samples, 19 traits + ID
n_samples <- 964
n_traits <- 19
n_proteins <- 446

Clinical <- data.frame(
  ID = sprintf("ID_%04d", 1:n_samples),
  matrix(rnorm(n_samples * n_traits), ncol = n_traits)
)
colnames(Clinical)[-1] <- paste0("Trait", 1:n_traits)

# Simulate protein data: 964 samples, 446 proteins + matching ID
Protein <- data.frame(
  ID = Clinical$ID,
  matrix(rnorm(n_samples * n_proteins), ncol = n_proteins)
)
colnames(Protein)[-1] <- paste0("Prot", 1:n_proteins)

# ----------------------------- Merge Clinical and Protein Data -----------------------------
ClinicalProtein <- merge(Clinical, Protein, by = "ID")

# ----------------------------- Define Identifiers -----------------------------
ids <- names(ClinicalProtein)[1]
traits <- names(ClinicalProtein)[2:(n_traits + 1)]
proteins <- names(ClinicalProtein)[(n_traits + 2):ncol(ClinicalProtein)]

# ----------------------------- Bayesian Network Model Function -----------------------------
fit.the.model <- function(data, alpha) {
  cpc <- vector("list", length(traits))
  names(cpc) <- traits
  for (t in seq_along(traits)) {
    cpc[[t]] <- learn.nbr(data[, c(traits, proteins)], node = traits[t],
                          method = "si.hiton.pc", test = "cor", alpha = alpha)
  }
  nodes <- unique(c(traits, unlist(cpc)))
  hc(data[, nodes])
}

# ----------------------------- k-Fold Cross-Validation Function -----------------------------
xval.the.model <- function(data, k = 10, cluster, alpha) {
  n <- nrow(data)
  predcor <- numeric(length(traits))
  names(predcor) <- traits
  postcor <- numeric(length(traits))
  names(postcor) <- traits
  
  kcv <- split(sample(n), rep(1:k, length.out = n))
  
  predicted <- parLapply(kcv, cl = cluster, function(test) {
    dtraining <- data[-test, ]
    dtest <- data[test, ]
    model <- fit.the.model(dtraining, alpha)
    fitted <- bn.fit(model, dtraining[, nodes(model)])
    dtest <- dtest[, nodes(model), drop = FALSE]
    
    available_traits <- traits[traits %in% nodes(model)]
    pred_partial <- sapply(available_traits, function(t) predict(fitted, node = t, data = dtest))
    pred_full <- matrix(NA, nrow = nrow(dtest), ncol = length(traits))
    colnames(pred_full) <- traits
    pred_full[, available_traits] <- pred_partial
    
    post_list <- lapply(seq(nrow(dtest)), function(i) {
      evidence <- as.list(dtest[i, names(dtest) %in% proteins])
      result <- rep(NA, length(traits))
      names(result) <- traits
      available <- traits[traits %in% nodes(model)]
      result[available] <- colMeans(cpdist(fitted, nodes = available, evidence = evidence,
                                           method = "lw", n = 1000))
      result
    })
    post_full <- do.call(rbind, post_list)
    
    list(model = fitted, pred = pred_full, post = post_full)
  })
  
  causal <- do.call(rbind, lapply(predicted, `[[`, "pred"))
  posterior <- do.call(rbind, lapply(predicted, `[[`, "post"))
  
  for (t in traits) {
    predcor[t] <- cor(causal[, t], data[unlist(kcv), t], use = "pairwise.complete.obs")
    postcor[t] <- cor(posterior[, t], data[unlist(kcv), t], use = "pairwise.complete.obs")
  }
  
  list(predicted = causal, posterior = posterior,
       observed = data[unlist(kcv), traits],
       predcor = predcor, postcor = postcor,
       models = lapply(predicted, `[[`, "model"))
}

# ----------------------------- Train and Evaluate Models -----------------------------
cl <- makeCluster(3)
clusterEvalQ(cl, library(bnlearn))
clusterExport(cl = cl, c("traits", "proteins", "ids", "fit.the.model"))

pr001 <- lapply(1:3, function(i) xval.the.model(ClinicalProtein, cluster = cl, alpha = 0.01))
stopCluster(cl)

pred.summary <- sapply(pr001, `[[`, "predcor")
print(rowMeans(pred.summary, na.rm = TRUE))

post.summary <- sapply(pr001, `[[`, "postcor")
print(rowMeans(post.summary, na.rm = TRUE))

# ----------------------------- Construct Averaged Network -----------------------------
arclist <- unlist(lapply(pr001, function(x) lapply(x$models, arcs)), recursive = FALSE)
nodes <- unique(unlist(arclist))
strength <- custom.strength(arclist, nodes = nodes)
averaged <- averaged.network(strength)

relevant.nodes <- nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
averaged_relevant <- subgraph(averaged, relevant.nodes)
strength_relevant <- strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]

saveRDS(list(network = averaged_relevant, strengths = strength_relevant), "~/network_simulated.rds")

# ----------------------------- Filter by Strength & Export -----------------------------
strength80_50 <- strength_relevant[strength_relevant$strength >= 0.8 & strength_relevant$direction >= 0.5, ]
averaged80_50 <- averaged.network(strength80_50)

strength80_80 <- strength_relevant[strength_relevant$strength >= 0.8 & strength_relevant$direction >= 0.8, ]
averaged80_80 <- averaged.network(strength80_80)

wb <- createWorkbook()
addWorksheet(wb, "Arcs_80_50")
addWorksheet(wb, "Arcs_80_80")
writeData(wb, sheet = "Arcs_80_50", strength80_50)
writeData(wb, sheet = "Arcs_80_80", strength80_80)
saveWorkbook(wb, file = "~/Simulated_ArcStrengthDirection.xlsx", overwrite = TRUE)

# ----------------------------- Compare BIC Scores to Random -----------------------------

# We need a learned network without any undirected arcs since in this random data averaged80_50 had undirected arcs, averaged 80_60 was generated

# Filter for arcs with strength ≥ 0.8 and direction ≥ 0.6
strength80_60 <- strength_relevant[strength_relevant$strength >= 0.8 & strength_relevant$direction >= 0.6, ]
averaged80_60 <- averaged.network(strength80_60)

# Extract nodes with at least one connection
relevant.nodes80_60 <- nodes(averaged80_60)[sapply(nodes(averaged80_60), degree, object = averaged80_60) > 0]

# Subset the averaged network and data
averaged80_60 <- subgraph(averaged80_60, relevant.nodes80_60)
DataSubset_80_60 <- ClinicalProtein[, relevant.nodes80_60]

# Fit the learned network
fitted_bn <- bn.fit(averaged80_60, data = DataSubset_80_60)
bic_learned <- BIC(fitted_bn, DataSubset_80_60)

# Generate 1000 random networks and compute BICs
random_bics <- numeric(1000)
for (i in 1:1000) {
  random_bn <- random.graph(nodes = relevant.nodes80_60)
  fitted_random <- bn.fit(random_bn, DataSubset_80_60)
  random_bics[i] <- BIC(fitted_random, DataSubset_80_60)
}

# Set x-axis range for consistent scaling
x_range <- range(c(random_bics, bic_learned + 10))

# ----------------------------- Plot Histogram and Density -----------------------------
hist(random_bics, breaks = 50, col = "lightgray", main = "BIC Comparison: Learned vs. Random BNs",
     xlab = "BIC Score", ylab = "Frequency", border = "black", xlim = x_range, probability = TRUE)

lines(density(random_bics), col = "blue", lwd = 2)
abline(v = bic_learned, col = "red", lwd = 3, lty = 2)
text(x = bic_learned, y = max(density(random_bics)$y) * 0.7,
     labels = paste("BIC =", round(bic_learned, 2)), col = "red", font = 2, pos = 4, cex = 0.8)

legend("topright", legend = c("Learned BN BIC", "Density of Random BICs"),
       col = c("red", "blue"), lwd = 2, lty = c(2, 1), box.lty = 0, cex = 0.8)
