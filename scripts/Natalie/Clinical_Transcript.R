############################################  Loading libraries #########################################################
library(reshape)
library(readr)
library(dplyr, warn.conflicts = FALSE)
library(lme4)
library(bnlearn)
library(parallel)


############################################ Reading Transcriptomic #####################################################

reads <- read.table("/home/Teams/teamVIP/Analysis/MultiOmics/Data/Transcriptomics/RPKM_WP2_DIRECT_TechAgeGenderResiduals_24082017.rank_TSScorrected.bedv6.gz",sep="\t", header=T, comment.char="", check.names=FALSE)##16209  3035
TransposedReads <- t(reads[,7:ncol(reads)])
TransposedReads <- as.data.frame(TransposedReads)#(3029 16209)
ReadsNames <- reads[,1:6]##(16209     6)
colnames(TransposedReads) <- sub("\\.\\d+$", "", ReadsNames$GeneID)
names <- rownames(TransposedReads)
TransposedReads <- cbind(names, TransposedReads)
TransposedReads <- TransposedReads %>%
  rename( studyid= names)
row.names(TransposedReads) <-  NULL

##################################################### Reading WP2.2 #####################################################

list.files("/users/home/Data/Repository/MultiOmics/WP2/2019_07_18/ClinicalVariables/Matrices/_10")
Data2_2 <- read_tsv('/users/home/Data/Repository/MultiOmics/WP2/2019_07_18/ClinicalVariables/Matrices/_10/merged_pheno_matrix_short_raw_2-2_direct_19_03_2019_10.txt',show_col_types=FALSE) #795  85
SubsetD.2 <- subset(Data2_2, select=c(  
"studyid","center","CRF_2_100_eligibility.age_at_visit","NewGender","CRF_2_102_general.subject_height",
"CRF_2_102_general.subject_weight","CRF_2_102_general.subject_waist","bmi","BP_s_mean","BP_d_mean","treat.cat","metformin",                           
"fasting.HbA1c","fasting.Glucose","fasting.Insulin","mmtt.120.Glucose","mmtt.120.Insulin","fasting.HDL","fasting.LDL",                        
"fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr",                           
"glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","active_glp1_conc_0min",               
"total_glp1_conc_0min","Glucagon_conc_0min_pg_ml","glp1.inc.60",                         
"glucagon.inc.60","liver.iron","panc.iron","liver.fat" ,"panc.fat", "iaat","taat" ,"asat" ))
## All Diebetic , "Diab_baseline"
##"glucagon.inc.120" with 600 NA's
##"proinsulin_conc_60min_pmol_L" with 400 NA's
##"Pass_InclusionExclusion" No :  6   and Yes:789
colnames(SubsetD.2) <- c("studyid","center","Age","NewGender","Height",
"Weight","Waist","bmi","BP_s_mean","BP_d_mean","treat.cat","metformin", 
"fasting.HbA1c","fasting.Glucose","fasting.Insulin","TwoH.Glucose","TwoH.Insulin","fasting.HDL","fasting.LDL",  
"fasting.TG","fasting.ALT","fasting.AST","fasting.Chol","mean.glu","mean.ins","basal.isr", 
"glu.sens","rate.sens","PFR1","total.isr","X2.h.OGIS","Stumvoll","Matsuda","active_glp1_conc_0min",
"total_glp1_conc_0min","Glucagon_conc_0min_pg_ml","glp1.inc.60", 
"glucagon.inc.60","liver.iron","panc.iron","liver.fat" ,"panc.fat", "iaat","taat" ,"asat")
#SubsetD.2$DiabetesStatus <- factor(as.numeric(gsub("([0-9]+)[-]([0-9]+)[-]([0-9]+).*", "\\2",SubsetD.2$studyid))) ##795,46
GGTP2.2 <- read.table("/users/home/pasdar/2019_01_11_Exeter_GGTP_WP2.2.tab", sep="\t", header=TRUE)##(818   8)
GGTP2.2 <- GGTP2.2[,c("studyid","wpid","result_value")]  
GGTP2.2$wpid <- as.factor(GGTP2.2$wpid)
colnames(GGTP2.2) <- c("studyid","wpid","GGTP")
Data2.2full <- read.table("/users/home/Data/Repository/MultiOmics/WP2/2019_07_18/ClinicalVariables/Matrices/_10/merged_pheno_matrix_full_raw_2-2_direct_19_03_2019_10.txt", sep="\t", header=TRUE, na.string=c("nk", "", "NA", "not known","N/A"))#795 372
Data2.2full <- Data2.2full[,c("studyid", "CRF_2_102_general.subject_waist","CRF_2_102_general.subject_hip",
"CRF_2_100_medical_history_1.hypercholesterolaemia","CRF_2_100_diabetes_history.family_history_parent",
"CRF_2_100_diabetes_history.family_history_sibling",
"CRF_2_100_diabetes_history.family_history_child","Clinsb","Clins")]##795  9
Data2.2full$WaistToHip <- Data2.2full$CRF_2_102_general.subject_waist/Data2.2full$CRF_2_102_general.subject_hip  
Data2.2full$Hx <- "no"
Data2.2full$Hx[Data2.2full$CRF_2_100_diabetes_history.family_history_parent=="yes"| Data2.2full$CRF_2_100_diabetes_history.family_history_sibling=="yes"| Data2.2full$CRF_2_100_diabetes_history.family_history_child=="yes"] <- "yes"
Data2.2full$Hx <- as.factor(Data2.2full$Hx)
Data2.2full <- subset(Data2.2full, select=c("studyid", "WaistToHip","CRF_2_100_medical_history_1.hypercholesterolaemia","Hx","Clinsb","Clins"))#503   8
colnames(Data2.2full) <- c("studyid","WaistToHip","hypercholesterolaemia","Hx","Clinsb","Clins")
SubsetD.2$ASTtoALT <- SubsetD.2$fasting.AST/SubsetD.2$fasting.ALT ##795  46
SubsetD.2 <- merge(SubsetD.2, Data2.2full, by="studyid") ##795  51
SubsetD.2 <- merge(SubsetD.2,GGTP2.2, by="studyid")##793  53
SubsetD.2 <- SubsetD.2[complete.cases(SubsetD.2[,c("liver.fat")]),] ##518  53
Sub2.2 <- SubsetD.2[,c("studyid", "center","Age","NewGender","bmi","Waist","BP_s_mean" ,"BP_d_mean" , 
	"fasting.Glucose" ,"fasting.HbA1c","fasting.Insulin" ,"TwoH.Glucose","TwoH.Insulin",
	"fasting.HDL","fasting.TG","fasting.ALT","fasting.AST" ,            
	"basal.isr" ,"glu.sens","X2.h.OGIS" ,
	"total_glp1_conc_0min","Glucagon_conc_0min_pg_ml",
     "liver.iron", "panc.iron","liver.fat" ,"panc.fat" ,               
	"iaat","asat","Clins" , "wpid","GGTP")]#518  31
Sub2.2_comp <- Sub2.2[complete.cases(Sub2.2),] #281  31


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
GGTP2.1 <- read.table("/users/home/pasdar/2019_01_11_Exeter_GGTP_WP2.1.tab", sep="\t", header=TRUE)##(818   8)
GGTP2.1 <- GGTP2.1[,c("studyid","wpid","result_value")]  
GGTP2.1$wpid <- as.factor(GGTP2.1$wpid)
colnames(GGTP2.1) <- c("studyid","wpid","GGTP")
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
SubsetD.1 <- merge(SubsetD.1,GGTP2.1, by="studyid")##2200   52
SubsetD.1 <- SubsetD.1[complete.cases(SubsetD.1[,c("liver.fat")]),] ##1008  52
SubsetD.1$fasting.Insulin <- (SubsetD.1$fasting.Insulin)*6
SubsetD.1$ogtt.120.Insulin <- (SubsetD.1$ogtt.120.Insulin)*6
Sub2.1 <- SubsetD.1[,c("studyid", "center","CRF_1_100_section_a.age_at_visit","NewGender","bmi","CRF_1_101_section_a.waist","BP_s_mean" ,"BP_d_mean" , 
	"fasting.Glucose" ,"fasting.HbA1c","fasting.Insulin" ,"ogtt.120.Glucose","ogtt.120.Insulin",
	"fasting.HDL","fasting.TG","fasting.ALT","fasting.AST" ,            
	"basal.isr" ,"glu.sens","X2.h.OGIS" ,
	"total_glp1_conc_0min","Glucagon_conc_0min_pg_ml",
     "liver.iron", "panc.iron","liver.fat" ,"panc.fat" ,               
	"iaat","asat","Clins" , "wpid","GGTP")]#1008  31
 Sub2.1_comp <- Sub2.1[complete.cases(Sub2.1),] #773  31

##################################################### Merging WP2.1 and WP2.2 #####################################################

colnames(Sub2.1_comp) <- c("studyid", "center" ,"Age","Sex" ,                      
"bmi" , "Waist", "SBP","DBP",                  
"Glucose" ,"HbA1c","Insulin" ,"TwoGlucose","TwoInsulin",
	"HDL","TG","ALT","AST" ,            
	"BasalISR" ,"GlucoseSens","OGIS" ,
	"TotGLP1min0","Glucagonmin0",
     "LiverIron", "PancIron","LiverFat" ,"PancFat" ,               
	"VAT","SAT","Clins" , "wpid","GGTP")
colnames(Sub2.2_comp) <- c("studyid", "center" ,"Age","Sex" ,                      
"bmi" , "Waist", "SBP","DBP",                  
"Glucose" ,"HbA1c","Insulin" ,"TwoGlucose","TwoInsulin",
	"HDL","TG","ALT","AST" ,            
	"BasalISR" ,"GlucoseSens","OGIS" ,
	"TotGLP1min0","Glucagonmin0",
     "LiverIron", "PancIron","LiverFat" ,"PancFat" ,               
	"VAT","SAT","Clins" , "wpid","GGTP")
MergedAll <- rbind(Sub2.1_comp,Sub2.2_comp) ##1054   31
MergedAll <- droplevels(MergedAll) #1054   31
MergedAll$center <- as.factor(MergedAll$center)
MergedAll$Sex <- as.factor(MergedAll$Sex)
MergedAll$GGTP <- as.numeric(MergedAll$GGTP)

##################################################### Normalizing Clinical data #####################################################

#Since we will be using Gaussian BNs for the analysis, it also interesting to check whether the variables are normally distributed, at least marginally; and from the plots below that does not seem to be the case for all of them.
F5 <- apply(MergedAll[c(5:29,31)], 2, FUN=function(x){lm( x~MergedAll[,2]+MergedAll[,3]+MergedAll[,4])}) #Adjust for Center, Age, Sex
residuals5 <- sapply(F5, function(x){resid(x)})
residuals5 <- as.data.frame(residuals5) 
RankedCompMergedAll <- sapply(residuals5, function(x){qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))})
RankedCompMergedAll <- as.data.frame(RankedCompMergedAll) #1054 26
RankedCompMergedAll$studyid <-MergedAll$studyid 
#RankedCompMergedAll$wpid <-MergedAll$wpid  #1054 28

##################################################### Merging Clinical and Transcript #####################################################

TraitsTransposedReadsSub <- merge(RankedCompMergedAll,TransposedReads, by="studyid") #1054 16237

##################################################### Preparing the BN #####################################################
