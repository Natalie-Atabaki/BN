############################################  Loading libraries #########################################################
library(reshape)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(lme4)
library(bnlearn)
library(parallel)
library(GGally)
library(gplots)
library(Rgraphviz)
library(openxlsx)

##################################################### Reading WP2.1 #####################################################

Data2_1 <- read_tsv('/users/home/Data/Repository/MultiOmics/WP2/2019_07_18/ClinicalVariables/Matrices/_10/merged_pheno_matrix_short_raw_2-1_direct_19_03_2019_10.txt' ,show_col_types=FALSE)#2234   87
SubsetD.1 <- subset(Data2_1, select=c(  
  "studyid","center","CRF_1_100_section_a.age_at_visit","NewGender",
  "CRF_1_101_section_a.height" ,"CRF_1_101_section_a.weight" ,"CRF_1_101_section_a.waist","bmi","BP_s_mean","BP_d_mean", "Hx_yes",                          
  "fasting.HbA1c","fasting.Glucose","fasting.Insulin","ogtt.120.Glucose","ogtt.120.Insulin","fasting.HDL","fasting.LDL",                        
  "fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr",                           
  "glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","active_glp1_conc_0min",               
  "total_glp1_conc_0min","Glucagon_conc_0min_pg_ml","glp1.inc.60",                         
  "glucagon.inc.60","liver.iron","panc.iron","liver.fat" ,"panc.fat", "iaat","taat" ,"asat" ))
Data2.1full <- read.table("/users/home/Data/Repository/MultiOmics/WP2/2019_07_18/ClinicalVariables/Matrices/_10/merged_pheno_matrix_full_raw_2-1_direct_19_03_2019_10.txt", sep="\t", header=TRUE, na.string=c("nk", "", "NA", "not known","N/A"))#795 372
Data2.1full <- Data2.1full[,c("studyid","CRF_1_101_section_a.waist","CRF_1_100_section_b.Hip"  
                              ,"CRF_1_100_section_d.hypercholesterolaemia","CRF_1_100_section_f.family_history_parent","CRF_1_100_section_f.family_history_sibling",
                              "CRF_1_100_section_f.family_history_child","Clinsb","Clins")]##1011   11
Data2.1full$WaistToHip <- Data2.1full$CRF_1_101_section_a.waist/Data2.1full$CRF_1_100_section_b.Hip 
Data2.1full$Hx <- "no"
Data2.1full$Hx[Data2.1full$CRF_1_100_section_f.family_history_parent=="yes"| Data2.1full$CRF_1_100_section_f.family_history_sibling=="yes"| Data2.1full$CRF_1_100_section_f.family_history_child=="yes"] <- "yes"
Data2.1full$Hx <- as.factor(Data2.1full$Hx)
Data2.1full <- subset(Data2.1full, select=c("studyid", "WaistToHip","CRF_1_100_section_d.hypercholesterolaemia","Hx","Clinsb","Clins"))#2234    6
SubsetD.1$ASTtoALT <- SubsetD.1$fasting.AST/SubsetD.1$fasting.ALT ##2282   68
SubsetD.1 <- merge(SubsetD.1, Data2.1full, by="studyid") ##2234   73
SubsetD.1 <- SubsetD.1[complete.cases(SubsetD.1[,c("liver.fat")]),] ##1008  52 ###1025  50
#Fasting insulin (pmol/L)  but we need to convert
#1 μIU/mL = 6.00 pmol/L
SubsetD.1$fasting.Insulin <- (SubsetD.1$fasting.Insulin)*6
SubsetD.1$ogtt.120.Insulin <- (SubsetD.1$ogtt.120.Insulin)*6
SubsetD.1$HOMA_IR <- (SubsetD.1$fasting.Insulin*SubsetD.1$fasting.Glucose)/22.5 #1025 51
SubsetD.1$DiabetesStatus <- factor(as.numeric(gsub("([0-9]+)[-]([0-9]+)[-]([0-9]+).*", "\\2",SubsetD.1$studyid))) 

Sub2.1_small <- SubsetD.1[,c("studyid", "center","CRF_1_100_section_a.age_at_visit","NewGender",
                             "fasting.Glucose" ,"fasting.HbA1c","fasting.Insulin" ,"ogtt.120.Glucose","ogtt.120.Insulin",
                             "fasting.HDL","fasting.TG",            
                             "basal.isr" ,"glu.sens","X2.h.OGIS" ,
                             "total_glp1_conc_0min","Glucagon_conc_0min_pg_ml",
                             "liver.fat" ,"panc.fat" , "iaat","asat","Clins" ,"Clinsb", "DiabetesStatus","HOMA_IR")]#1025  24
Sub2.1_comp_small <- Sub2.1_small[complete.cases(Sub2.1_small),] #964  24

##################################################### Merging WP2.1 and WP2.2 #####################################################

colnames(Sub2.1_comp_small) <- c("studyid", "center" ,"Age","Sex" ,                      
                                 "Glucose" ,"HbA1c","Insulin" ,"TwoGlucose","TwoInsulin",
                                 "HDL","TG","BasalISR" ,"GlucoseSens","OGIS" ,
                                 "TotGLP1min0","Glucagonmin0","LiverFat" ,"PancFat" ,               
                                 "VAT","SAT","Clins" , "Clinsb","DiabetesStatus","HOMA_IR")

Sub2.1_comp_small <- droplevels(Sub2.1_comp_small) #
Sub2.1_comp_small$center <- as.factor(Sub2.1_comp_small$center)
Sub2.1_comp_small$Sex <- as.factor(Sub2.1_comp_small$Sex)
write.table(Sub2.1_comp_small,"/users/home/pasdar/WorkingDir/BN/BN_11_11/Sub_1_NoAdje_24_964.tab", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
