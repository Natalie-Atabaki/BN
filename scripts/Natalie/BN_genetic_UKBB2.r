
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(purrr)
library(rbgen)
library(ggplot2)


tab_location1 <- "/ludc/Active_Projects/UKBB_18274/Private/ukb26608.tab"

tab_vars1 <- tribble(
    ~htmlpos, ~colnam,
    0, "eid",
    4, "sex",
    15, "waist_c",
    18, "hip_c",
    200,"DiabetesDiagnosed",
    338, "DBP",
    344, "SBP",
    825, "AlcoholStatus",
    837, "Ethnicity",
    840, "BMI",
    843, "Weight",
    846, "age",
    849, "ageRecruit",
    851, "genetic_sex",
    855, "genetic_grouping",
    ## We'll use here only 5 principal components for simplicity, but generally 10 are used
    858, "gpc1",
    859, "gpc2",
    860, "gpc3",
    861, "gpc4",
    862, "gpc5",
    863, "gpc6",
    864, "gpc7",
    865, "gpc8",
    866, "gpc9",
    867, "gpc10",
    898, "sex_chromosome_aneuploidy",
    899, "used_in_gpc_calculation",
    914, "VAT",
    915, "SAT",
    917, "TrunkFat",
    951, "BodyFatMass",
    961, "ImpedanceWholeBody",
    1008, "TrunkFatPercentage",
    1076, "OverallAccelerationAverage") 

tab_vars1 <- tab_vars1 %>%
    mutate(colpos = htmlpos + 1)

field_indices1 <- tab_vars1 %>% pull(colpos) %>% paste(collapse = ",")

pheno_dat_cmd1 <- paste0("cut -f", field_indices1, " ", tab_location1)

pheno_dat1 <- pheno_dat_cmd1 %>% pipe %>% read_tsv(skip = 1, col_names = tab_vars1$colnam, show_col_types = FALSE)

pheno_dat1_sub <- pheno_dat1 %>% filter(sex == genetic_sex, genetic_grouping == 1, is.na(sex_chromosome_aneuploidy), used_in_gpc_calculation == 1) 
    
pheno_dat1_sub <- subset(pheno_dat1_sub,  select=-c(genetic_sex, genetic_grouping, sex_chromosome_aneuploidy, used_in_gpc_calculation,  hip_c))


tab_location2 <- "/ludc/Active_Projects/UKBB_18274/Private/ukb26876.tab"

tab_vars2 <- tribble(
    ~htmlpos, ~colnam,
    0, "eid",
    1, "LiverFat"
)
 
tab_vars2 <- tab_vars2 %>%
    mutate(colpos = htmlpos + 1)


field_indices2 <- tab_vars2 %>% pull(colpos) %>% paste(collapse = ",")

pheno_dat_cmd2 <- paste0("cut -f", field_indices2, " ", tab_location2)

pheno_dat2 <- pheno_dat_cmd2 %>% pipe %>% read_tsv(skip = 1, col_names = tab_vars2$colnam, show_col_types = FALSE)

Pheno12 <- merge(pheno_dat2,pheno_dat1_sub, by="eid")



tab_location3 <- "/ludc/Active_Projects/UKBB_18274/Private/ukb27539.tab"

tab_vars3 <- tribble(
    ~htmlpos, ~colnam,
    0, "eid",
    1, "LiverIron",
    2,"LIF")

tab_vars3 <- tab_vars3 %>%
    mutate(colpos = htmlpos + 1)

field_indices3 <- tab_vars3 %>% pull(colpos) %>% paste(collapse = ",")

pheno_dat_cmd3 <- paste0("cut -f", field_indices3, " ", tab_location3)


pheno_dat3 <- pheno_dat_cmd3 %>% pipe %>% read_tsv(skip = 1, col_names = tab_vars3$colnam, show_col_types = FALSE)

Pheno123 <- merge(Pheno12,pheno_dat3, by="eid")

tab_location4 <- "/ludc/Active_Projects/UKBB_18274/Private/ukb28686.tab"

tab_vars4 <- tribble(
    ~htmlpos, ~colnam,
    0, "eid",
    1, "Albumin",
    3, "ALT",
    5, "AST",
    7, "DirectBilirubin",
    9, "Cholesterol",
    11, "GGTP",
    13, "Glucose",
    15, "HbA1c",
    17, "HDL",
    19, "LDL",
    21, "TotalBilirubin",
    23, "Triglycerides"
) 
tab_vars4 <- tab_vars4 %>%
    mutate(colpos = htmlpos + 1)

field_indices4 <- tab_vars4 %>% pull(colpos) %>% paste(collapse = ",")

pheno_dat_cmd4 <- paste0("cut -f", field_indices4, " ", tab_location4)

pheno_dat4 <- pheno_dat_cmd4 %>% pipe %>% read_tsv(skip = 1, col_names = tab_vars4$colnam, show_col_types = FALSE)

Pheno1234 <- merge(Pheno123,pheno_dat4, by="eid")

Pheno1234$sex <- factor(Pheno1234$sex)
Pheno1234$DiabetesDiagnosed <- factor(Pheno1234$DiabetesDiagnosed)
Pheno1234$AlcoholStatus <- factor(Pheno1234$AlcoholStatus)

Pheno1234_sub <- subset(Pheno1234, select=-c(Ethnicity,LiverIron,LIF,AlcoholStatus,TrunkFat,ImpedanceWholeBody,
          TrunkFatPercentage,OverallAccelerationAverage,DirectBilirubin,Cholesterol,Glucose,LDL,ageRecruit))

Pheno1234_sub_LiverComp <- Pheno1234_sub[complete.cases(Pheno1234_sub[,c("LiverFat")]),]#1011  577

Pheno1234_sub_CompAll <- Pheno1234_sub_LiverComp[complete.cases(Pheno1234_sub_LiverComp),]

F5 <- apply(Pheno1234_sub_CompAll[c(2,4,6:9,21:31)], 2, FUN=function(x){lm( x~Pheno1234_sub_CompAll[,3]+
Pheno1234_sub_CompAll[,10]+Pheno1234_sub_CompAll[,11]+Pheno1234_sub_CompAll[,12]+Pheno1234_sub_CompAll[,13]
+Pheno1234_sub_CompAll[,14]+Pheno1234_sub_CompAll[,15]+Pheno1234_sub_CompAll[,16]+Pheno1234_sub_CompAll[,17]
+Pheno1234_sub_CompAll[,18]+Pheno1234_sub_CompAll[,19]+Pheno1234_sub_CompAll[,20])}) #Adjust for Sex 'age''gpc1''gpc2''gpc3''gpc4''gpc5''gpc6''gpc7''gpc8''gpc9''gpc10'
residuals5 <- sapply(F5, function(x){resid(x)})
residuals5 <- as.data.frame(residuals5)
residuals5 <- droplevels(residuals5) 

RankedComp <- sapply(residuals5, function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))})
RankedComp <- as.data.frame(RankedComp) 

RankedComp$eid <-Pheno1234_sub_CompAll$eid #2620   18
#RankedComp$eid < as.character(RankedComp$eid)############ Check!!!

library("reshape")
library("MASS")
library("stats")
#library("bnlearn")
library("GGally")
library("gplots")
library("tidyverse")
#library("Rgraphviz")
library("gRain")
#library("parallel")
library("devtools")


ALT_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Alanine aminotransferase.csv", show_col_types = FALSE)
ASTlevels_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Aspartate aminotransferase levels.csv", show_col_types = FALSE)
AST_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Aspartate aminotransferase.csv", show_col_types = FALSE)
BodyfatPercent_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Body fat percentage.csv", show_col_types = FALSE)
BMI_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Body mass index.csv", show_col_types = FALSE)
DBP_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Diastolic blood pressure  automated reading.csv", show_col_types = FALSE)
fastingGlucose_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Fasting blood glucose.csv", show_col_types = FALSE)
fastingInsulin_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Fasting blood insulin.csv", show_col_types = FALSE)
GGTPlevels_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Gamma glutamyl transferase levels.csv", show_col_types = FALSE)
GGTP_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Gamma glutamyltransferase.csv", show_col_types = FALSE)
GlycatedHaemoglobin_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Glycated haemoglobin.csv", show_col_types = FALSE)
HDL_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/HDL cholesterol.csv", show_col_types = FALSE)
HbA1c_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/HbA1C.csv", show_col_types = FALSE)
StumvollBMIadj_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Modified Stumvoll Insulin Sensitivity Index (model adjusted for BMI).csv", show_col_types = FALSE)
NAFLD_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Nonalcoholic fatty liver disease.csv", show_col_types = FALSE)
Albumin_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Serum albumin level.csv", show_col_types = FALSE)
SBP_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Systolic blood pressure  automated reading.csv", show_col_types = FALSE)
bilirubin_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Total bilirubin.csv", show_col_types = FALSE)
TG_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Triglycerides.csv", show_col_types = FALSE)
twoHourGlucose_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Two-hour glucose challenge.csv", show_col_types = FALSE)
twoDiabetesBMIadj_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Type 2 diabetes (adjusted for BMI).csv", show_col_types = FALSE)
twoDiabetes_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Type 2 diabetes.csv", show_col_types = FALSE)
Waist_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Waist circumference.csv", show_col_types = FALSE)
bodyFatMass_snps <- read_csv("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/GWAS_hits/Whole body fat mass.csv", show_col_types = FALSE)


ClumpSome_snps <- rbind(ASTlevels_snps,AST_snps,BodyfatPercent_snps,BMI_snps,DBP_snps,fastingGlucose_snps,
fastingInsulin_snps,GGTPlevels_snps, GGTP_snps, GlycatedHaemoglobin_snps, HDL_snps, HbA1c_snps,
StumvollBMIadj_snps,NAFLD_snps,Albumin_snps, SBP_snps, bilirubin_snps, TG_snps, twoHourGlucose_snps,
twoDiabetesBMIadj_snps, twoDiabetes_snps,Waist_snps,bodyFatMass_snps)

ClumpSome_snps_ColSub <- subset(ClumpSome_snps, select=c("chr","rsid","ea","nea"))


ClumpSome_snps_ColSub_ord <- ClumpSome_snps_ColSub[order(ClumpSome_snps_ColSub$rsid),]


ClumpSome_snps_ColSub_ord_Uniq <- ClumpSome_snps_ColSub_ord %>% distinct(rsid, .keep_all = TRUE)


ClumpSome_snps_ColSub_ord_Uniq <- ClumpSome_snps_ColSub_ord_Uniq %>%
    filter(chr != "X")

ClumpSome_snps_ColSub_ord_Uniq$chr <- as.numeric(ClumpSome_snps_ColSub_ord_Uniq$chr)

clumpSNP <- read.table("/ludc/Home/naeimeh_a/clumpSNP", header=TRUE, quote = "\"")
Sub <- ClumpSome_snps_ColSub_ord_Uniq[ClumpSome_snps_ColSub_ord_Uniq$rsid %in% clumpSNP$rsid,]#745   4

ClumpSome_genotypes <- Sub %>%
    group_split(chr) %>%
    map(
        ~bgen.load(
            filename = paste0("/ludc/Raw_Data_Archive/UKBB/imp/ukb_imp_chr", unique(.x$chr), "_v3.bgen"),
            rsids = .x$rsid,
            index.filename = paste0("/ludc/Active_Projects/UKBB_18274/Private/Data/ukb_bgi/ukb_imp_chr", unique(.x$chr), "_v3.bgen.bgi")
        )
)


sample_ids <- read_table("/ludc/Active_Projects/UKBB_18274/Private/ukb18274_imp_chr1_v3_s487314.sample", show_col_types = FALSE) %>% slice(-1) ##487409      4

geno_dat <- ClumpSome_genotypes %>%
    ## One chromosome at a time
    map_dfc(~{
        ## SNP data
        snpdat <- .x$variants %>%
            tibble %>%
            ## Are alleles aligned?
            inner_join(Sub, by = "rsid") %>%
            mutate(
                aligned = case_when(
                    ea == allele1 & nea == allele0 ~ 0,  ## Aligned
                    ea == allele0 & nea == allele1 ~ -2, ## Unaligned
                    TRUE ~ NaN                           ## Do not correspond (needs to be dropped)
                )
            ) %>%
            filter(!is.na(aligned))
        ## One SNP at a time
        map2_dfc(
            snpdat$rsid, snpdat$aligned,
            function(snp, align){
                ## Extracting and simplifying genotype matrix
                gmat <- .x$data[snp,,] %*% c(0, 1, 2)
                ## Aligning and converting to data frame
                gdat <- data.frame(abs(gmat + align))
                ## Naming columns with rsids
                colnames(gdat) <- snp
                tibble(gdat) ## To avoid rownames
            }
        )
    }) %>%
    mutate(eid = sample_ids$ID_1) %>% ## Adding individual IDs
    inner_join(RankedComp) ## Joining with the phenotype data


library(lme4)
library(bnlearn)
library(parallel)

# split the variables.
ids = names(geno_dat)[591] ##608
traits = names(geno_dat)[592:ncol(geno_dat)]
genes = names(geno_dat)[1:590]

#ids = names(geno_dat)[2704]
#traits = names(geno_dat)[2705:ncol(geno_dat)]
#genes = names(geno_dat)[1:2703]

#Performing cross-validation
#The Bayesian networks model is fitted by the fit.the.model() function below,
# which takes the data and the type I error threshold alpha to use for structure learning as arguments.


fit.the.model = function(data, alpha) {

cpc = vector(length(traits), mode = "list")
names(cpc) = traits

# find the parents of each trait (may be genes or other traits).
for (t in seq_along(traits)) {

# BLUP away the family structure.
#m = lmer(as.formula(paste(traits[t], "~ (1|FUNNEL:PLANT)")), data = data)
#data[!is.na(data[, traits[t]]), traits[t]] =
#data[, traits[t]] - ranef(m)[[1]][paste(data$FUNNEL, data$PLANT, sep = ":"), 1]
# find out the parents.
cpc[[t]] = learn.nbr(data[, c(traits, genes)], node = traits[t], debug = FALSE,
method = "si.hiton.pc", test = "cor", alpha = alpha)
}#FOR

# merge the relevant variables to use for learning.
nodes = unique(c(traits, unlist(cpc)))
#blacklist = tiers2blacklist(list(nodes[!(nodes %in% traits)],
#traits[traits != "YLD"], "YLD"))
blacklist = tiers2blacklist(list(nodes[!(nodes %in% traits)],traits))

# build the Bayesian network.
bn = hc(data[, nodes], blacklist = blacklist)

return(bn)
}#FIT.THE.MODEL


#In order to have multiple models to average and to assess the predictive power of the Bayesian network
 # model, we run fit.the.model() under 10-fold cross-validation for 10 times. 
 #Since structure learning has been customised for the analysis, we cannot use bn.cv() from bnlearn. 
 #Instead, we implement a custom xval.the.model() function, shown below.

## k=10
xval.the.model = function(data, k = 10, cluster, alpha, ridge) {
n = nrow(data)
predcor = numeric(length(traits))
names(predcor) = traits
postcor = numeric(length(traits))
names(postcor) = traits

# shuffle the data to get unbiased splits.
kcv = split(sample(n), seq_len(k))
# store the length of each test set.
kcv.length = sapply(kcv, length)

predicted = parLapply(kcv, cl = cluster, function(test) {

# create a matrix to store the predicted values.
pred = matrix(0, nrow = length(test), ncol = length(traits))
colnames(pred) = traits
# create a matrix to store posterior estimates.
post = matrix(0, nrow = length(test), ncol = length(traits))
colnames(post) = traits

cat("* beginning cross-validation fold.\n")

# split training and test.
dtraining = data[-test, ]
dtest = data[test, ]
# fit the model on the training data.
model = fit.the.model(dtraining, alpha = alpha)
fitted = bn.fit(model, dtraining[, nodes(model)])
# maybe re-fit with ridge regression.
if (ridge) {

library(penalized)

for (no in nodes(fitted)) {

node.parents = parents(fitted, no)

if (length(node.parents) < 3)
next

opt.lambda = optL2(response = dtraining[, no],
penalized = dtraining[, node.parents],
model = "linear", trace = FALSE,
minlambda2 = 10e-5, maxlambda = 500)$lambda
fitted[[no]] = penalized(response = dtraining[, no],
penalized = dtraining[, node.parents],
model = "linear", trace = FALSE,
lambda1 = 0, lambda2 = opt.lambda)

}#FOR

}#THEN
# subset the test data.
dtest = dtest[, nodes(model)]

cat("  > model has", length(nodes(model)), "nodes.\n")

# predict each trait in turn, given all the parents.
for (t in traits)
pred[, t] = predict(fitted, node = t, data = dtest[, nodes(model)])

for (i in seq(nrow(dtest)))
post[i, traits] = colMeans(cpdist(fitted, nodes = traits,
evidence = as.list(dtest[i, names(dtest) %in% genes]),
method = "lw", n = 1000))

return(list(model = fitted, pred = pred, post = post))

})

# merge all the predicted values.
posterior = do.call(rbind, lapply(predicted, `[[`, "post"))
causal = do.call(rbind, lapply(predicted, `[[`, "pred"))

cat("* overall cross-validated correlations:\n")
for (t in traits) {

predcor[t] = cor(causal[, t], data[unlist(kcv), t])
cat("  > PREDCOR(", t, "):", predcor[t], "\n")
postcor[t] = cor(posterior[, t], data[unlist(kcv), t])
cat("  > POSTCOR(", t, "):", postcor[t], "\n")

}#FOR

return(list(predicted = causal, posterior = posterior,
observed = data[unlist(kcv), t], predcor = predcor, postcor = postcor,
models = lapply(predicted, `[[`, "model")))

}#XVAL.THE.MODEL


cl = makeCluster(10)##10
invisible(clusterEvalQ(cl, library(bnlearn)))
invisible(clusterEvalQ(cl, library(lme4)))
clusterExport(cl = cl, c("traits", "genes", "ids", "fit.the.model"))
pr001 = vector(10, mode = "list")
for (i in 1:10)
pr001[[i]] = xval.the.model(geno_dat, cluster = cl, alpha = 0.01, ridge = FALSE)
stopCluster(cl)



pred.summary = sapply(pr001, `[[`, "predcor")
print(rowMeans(pred.summary))
 
post.summary = sapply(pr001, `[[`, "postcor")
print(rowMeans(post.summary))



# gather all the arc lists.
arclist = list()
 
for (i in seq_along(pr001)) {

# extract the models.
run = pr001[[i]]$models

for (j in seq_along(run))
arclist[[length(arclist) + 1]] = arcs(run[[j]])

}#FOR


# compute the arc strengths.
nodes = unique(unlist(arclist))
strength = custom.strength(arclist, nodes = nodes)
# estimate the threshold and average the networks.
averaged = averaged.network(strength)

# subset the network to remove isolated nodes.
relevant.nodes = nodes(averaged)[sapply(nodes, degree, object = averaged) > 0]
averaged2 = subgraph(averaged, relevant.nodes)
strength2 = strength[(strength$from %in% relevant.nodes) & (strength$to %in% relevant.nodes), ]

library(Rgraphviz)

pdf("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/ClumpSNPsClinicalBNsSome.pdf")
par(cex=0.7)
gR = strength.plot(averaged, strength, shape = "rectangle", layout = "fdp")

nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$col = "darkblue"
nodeRenderInfo(gR)$fill[traits] = "limegreen"
nodeRenderInfo(gR)$col[traits] = "darkgreen"
a = arcs(subgraph(averaged, traits))
a = as.character(interaction(a[, "from"], a[, "to"], sep = "~"))
edgeRenderInfo(gR)$col = "grey"
edgeRenderInfo(gR)$col[a] = "darkgreen"
renderGraph(gR)
dev.off()


sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/relevantNodes.txt")
relevant.nodes
sink()
sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/ArcStrengthSub.txt")
strength
sink()
sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/AvgSub.txt")
averaged
sink()
sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/AvgmbSubMB.txt")
mb(averaged, "LiverFat")
sink()


averaged7 = averaged.network(strength[strength$strength > 0.7 & strength$direction > 0.7,])#80
strength7 = strength[strength$strength > 0.7 & strength$direction > 0.7,]
pdf("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/SNPsClinicalBN70BothSubSome.pdf")
par(cex=0.6)
gR = strength.plot(averaged7, strength7, shape = "rectangle", layout = "fdp")
nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$col = "darkblue"
nodeRenderInfo(gR)$fill[traits] = "limegreen"
nodeRenderInfo(gR)$col[traits] = "darkgreen"
a = arcs(subgraph(averaged, traits))
a = as.character(interaction(a[, "from"], a[, "to"], sep = "~"))
edgeRenderInfo(gR)$col = "grey"
edgeRenderInfo(gR)$col[a] = "darkgreen"
renderGraph(gR)
dev.off()


sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/ArcStrengthSub_Clean.txt")
strength2
sink()
sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/AvgSub_Clean.txt")
averaged2
sink()
sink("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/AvgmbSubMB_Clean.txt")
mb(averaged2, "LiverFat")
sink()


pdf("/ludc/Active_Projects/UKBB_18274/Private/UKBB_18274/Private/Network_project/BN_Genetic_UKBB/SNPsCleanClinicalBNsubSome2.pdf") 
par(cex=0.7)
gR = strength.plot(averaged2, strength2, shape = "rectangle", layout = "fdp")

nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$fill = "lightblue"
nodeRenderInfo(gR)$col = "darkblue"
nodeRenderInfo(gR)$fill[traits] = "limegreen"
nodeRenderInfo(gR)$col[traits] = "darkgreen"
a = arcs(subgraph(averaged2, traits))
a = as.character(interaction(a[, "from"], a[, "to"], sep = "~"))
edgeRenderInfo(gR)$col = "grey"
edgeRenderInfo(gR)$col[a] = "darkgreen"
renderGraph(gR)
dev.off()

