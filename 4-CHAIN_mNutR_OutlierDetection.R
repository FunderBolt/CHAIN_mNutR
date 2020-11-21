#############################
###  CHAIN_NutR
###  Outlier Detection: Samples or variables
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","FactoMineR", "caret","outliers")


# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)



# load full data
idata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_fulldata_2020-11-20.csv"), row.names = 1) # ND (not done) equivalent to NA
colnames(idata)


### list factors to convert
str(idata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)




##############################################
################  Sample outlier
colnames(idata)[17:ncol(idata)]


### pull out raw data
rawdata<-idata[ , c(7, 17:ncol(idata))]


### create log dataset and add in group annotation
logdata<-log(rawdata[,-1] + 1)
logdata<-cbind(idata$group_adm, logdata)



colnames(data)
## PCA on all metabolite 
PCA(rawdata, scale.unit = TRUE, quali.sup=1)
out_PCA<-PCA(logdata, scale.unit = TRUE, quali.sup=1)

plot(out_PCA, axes = c(1, 2), choix = c("ind"), ellipse = NULL, habillage=1,
     col.quali=c("red"))


dput(colnames(idata)[17:ncol(idata)])
### PCA only vitamines
out_PCA<-PCA(rawdata[,c("group_adm","albumin_adm","Retinol", "OH25_VitD3", "aTocopherol", "bg_Tocopherol", "B1", 
               "B2","B2_corr", "B3_amide", "B5", "B6","B6_corr", "B7", "C")], scale.unit = TRUE, quali.sup=1)

PCA(logdata[,c("albumin_adm","Retinol", "OH25_VitD3", "aTocopherol", "bg_Tocopherol", "B1", 
               "B2","B2_corr", "B3_amide", "B5", "B6","B6_corr", "B7", "C")], scale.unit = TRUE)


plot(out_PCA, axes = c(1, 2), choix = c("ind"), ellipse = NULL, habillage=1,
     col.quali=c("red"))


idata$a_aminobutyric_ac
dput(colnames(idata)[17:ncol(idata)])
### PCA only amino acids 
PCA(rawdata[,c("Creatinine", "Glycine", 
               "Alanine", "Serine", "Proline", "Valine", "Threonine", "Taurine", 
               "Putrescine", "trans.Hydroxyproline", "Leucine", "Isoleucine", 
               "Asparagine", "Aspartic_ac", "Glutamine", "Glutamic_ac", "Methionine", 
               "Histidine", "a_Aminoadipic_ac", "Phenylalanine", "Methionine.sulfoxide", 
               "Arginine", "Acetyl.ornithine", "Citrulline", "Serotonin", "Tyrosine", 
               "Asymmetric_dimethylarginine", "total_dimethylarginine", "Tryptophan", 
               "Kynurenine", "Carnosine", "Ornithine", "Lysine", "Spermidine", 
               "Spermine", "Sarcosine", "Creatine", "Betaine", "Choline", "Trimethylamine_N.oxide", 
               "Methylhistidine", "Homocysteine", "Hypoxanthine", "Ethanolamine", 
               "Cystathionine", "N.acetylpetrescine", "Homocitrulline", "Urea", 
               "Hydroxy.lysine", "a_aminobutyric_ac", "g_aminobutyric_ac", "Lactic_ac", 
               "b_Hydroxybutyric_ac", "a_Ketoglutaric_ac", "Citric_ac", "Butyric_ac", 
               "Propionic_ac", "HPHPA", "p.Hydroxyhippuric_ac", "Succinic_ac", 
               "Fumaric_ac", "Pyruvic_ac", "Isobutyric_ac", "Hippuric_ac", "Methylmalonic_ac", 
               "Homovanillic_ac", "Indole_acetic_ac", "Uric_ac", "OH5_Indoleacetic_ac", 
               "p.Hydroxyphenylacetic_ac", "Valeric_ac", "Uridine", "Xanthine", 
               "Guanidineacetic_ac", "Acetyl.lysine", "Glucose")], scale.unit = TRUE)

PCA(logdata[,c("Creatinine", "Glycine", 
               "Alanine", "Serine", "Proline", "Valine", "Threonine", "Taurine", 
               "Putrescine", "trans.Hydroxyproline", "Leucine", "Isoleucine", 
               "Asparagine", "Aspartic_ac", "Glutamine", "Glutamic_ac", "Methionine", 
               "Histidine", "a_Aminoadipic_ac", "Phenylalanine", "Methionine.sulfoxide", 
               "Arginine", "Acetyl.ornithine", "Citrulline", "Serotonin", "Tyrosine", 
               "Asymmetric_dimethylarginine", "total_dimethylarginine", "Tryptophan", 
               "Kynurenine", "Carnosine", "Ornithine", "Lysine", "Spermidine", 
               "Spermine", "Sarcosine", "Creatine", "Betaine", "Choline", "Trimethylamine_N.oxide", 
               "Methylhistidine", "Homocysteine", "Hypoxanthine", "Ethanolamine", 
               "Cystathionine", "N.acetylpetrescine", "Homocitrulline", "Urea", 
               "Hydroxy.lysine", "a_aminobutyric_ac", "g_aminobutyric_ac", "Lactic_ac", 
               "b_Hydroxybutyric_ac", "a_Ketoglutaric_ac", "Citric_ac", "Butyric_ac", 
               "Propionic_ac", "HPHPA", "p.Hydroxyhippuric_ac", "Succinic_ac", 
               "Fumaric_ac", "Pyruvic_ac", "Isobutyric_ac", "Hippuric_ac", "Methylmalonic_ac", 
               "Homovanillic_ac", "Indole_acetic_ac", "Uric_ac", "OH5_Indoleacetic_ac", 
               "p.Hydroxyphenylacetic_ac", "Valeric_ac", "Uridine", "Xanthine", 
               "Guanidineacetic_ac", "Acetyl.lysine", "Glucose")], scale.unit = TRUE)



dput(colnames(idata)[17:ncol(idata)])
### PCA only lipids
PCA(rawdata[,c("LYSOC14.0", 
               "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", "LYSOC18.1", 
               "LYSOC18.0", "LYSOC20.4", "LYSOC20.3", "LYSOC24.0", "LYSOC26.1", 
               "LYSOC26.0", "LYSOC28.1", "LYSOC28.0", "Sphg_14.1SMOH", "Sphg_16.1SM", 
               "Sphg_16.0SM", "Sphg_16.1SMOH", "Sphg_18.1SM", "PC32.2_AA", "Sphg_18.0SM", 
               "Sphg_20.2SM", "PC36.0_AE", "PC36.6_AA", "PC36.0_AA", "Sphg_22.2SMOH", 
               "Sphg_22.1SMOH", "PC38.6_AA", "PC38.0_AA", "PC40.6_AE", "Sphg_24.1SMOH", 
               "PC40.6_AA", "PC40.2_AA", "PC40.1_AA")], scale.unit = TRUE)


PCA(logdata[,c("LYSOC14.0", 
               "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", "LYSOC18.1", 
               "LYSOC18.0", "LYSOC20.4", "LYSOC20.3", "LYSOC24.0", "LYSOC26.1", 
               "LYSOC26.0", "LYSOC28.1", "LYSOC28.0", "Sphg_14.1SMOH", "Sphg_16.1SM", 
               "Sphg_16.0SM", "Sphg_16.1SMOH", "Sphg_18.1SM", "PC32.2_AA", "Sphg_18.0SM", 
               "Sphg_20.2SM", "PC36.0_AE", "PC36.6_AA", "PC36.0_AA", "Sphg_22.2SMOH", 
               "Sphg_22.1SMOH", "PC38.6_AA", "PC38.0_AA", "PC40.6_AE", "Sphg_24.1SMOH", 
               "PC40.6_AA", "PC40.2_AA", "PC40.1_AA")], scale.unit = TRUE)





dput(colnames(idata)[17:ncol(idata)])
### PCA only carnitines
PCA(rawdata[,c("C0", "C2", "C3.1", "C3", 
               "C4.1", "C4", "C3OH", "C5.1", "C5", "C4OH", "C6.1", "C6", "C5OH", 
               "C5.1DC", "C5DC", "C8", "C5MDC", "C9", "C7DC", "C10.2", "C10.1", 
               "C10", "C12.1", "C12", "C14.2", "C14.1", "C14", "C12DC", "C14.2OH", 
               "C14.1OH", "C16.2", "C16.1", "C16", "C16.2OH", "C16.1OH", "C16OH", 
               "C18.2", "C18.1", "C18", "C18.1OH")], scale.unit = TRUE)


PCA(logdata[,c("C0", "C2", "C3.1", "C3", 
               "C4.1", "C4", "C3OH", "C5.1", "C5", "C4OH", "C6.1", "C6", "C5OH", 
               "C5.1DC", "C5DC", "C8", "C5MDC", "C9", "C7DC", "C10.2", "C10.1", 
               "C10", "C12.1", "C12", "C14.2", "C14.1", "C14", "C12DC", "C14.2OH", 
               "C14.1OH", "C16.2", "C16.1", "C16", "C16.2OH", "C16.1OH", "C16OH", 
               "C18.2", "C18.1", "C18", "C18.1OH")], scale.unit = TRUE)




###################################
dput(colnames(idata)[17:ncol(idata)])
### PCA only metals
PCA(rawdata[,c("Barium", "Sodium", "Magnesium", 
               "Phosphorous", "Calcium", "Titanium", "Vanadium", "Manganese", 
               "Iron", "Cobalt", "Copper", "Zinc", "Gallium", "Arsenic", "Selenium", 
               "Rubidium", "Strontium", "Antimony", "Tellurium")], scale.unit = TRUE)


PCA(logdata[,c("Barium", "Sodium", "Magnesium", 
               "Phosphorous", "Calcium", "Titanium", "Vanadium", "Manganese", 
               "Iron", "Cobalt", "Copper", "Zinc", "Gallium", "Arsenic", "Selenium", 
               "Rubidium", "Strontium", "Antimony", "Tellurium")], scale.unit = TRUE)




#### NO outlier to remove --> sample 52, strong outlier to remove
#### almost removed 78 (--- > but before logtransform) to consider





colnames(idata)
list.plots<-colnames(idata)[17:ncol(idata)]
###################### Histograms for all metabolites

###create the histogram using a loop
for(i in 2:ncol(rawdata)){
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Histograms/",colnames(rawdata)[i],"_Histograms_raw.jpeg")
  jpeg(filename=file.name)
  hist(rawdata[,i], xlab = colnames(rawdata)[i], main = paste0("Histogram",colnames(rawdata)[i]))
  dev.off()
}


###create the histogram using a loop
for(i in 1:ncol(logdata)){
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Histograms/",colnames(logdata)[i],"_Histograms_log.jpeg")
  jpeg(filename=file.name)
  hist(logdata[,i], xlab = colnames(logdata)[i], main = paste0("Histogram",colnames(logdata)[i]))
  dev.off()
}




###################### Boxplots for all metabolites

###create the boxplots using a loop
for(i in 17:ncol(rawdata)){
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Boxplots/",colnames(rawdata)[i],"_Boxplots_raw.jpeg")
  jpeg(filename=file.name)
  boxplot(rawdata[,i]~idata$group_adm, xlab = colnames(rawdata)[i], main = paste0("Boxplots",colnames(rawdata)[i]))
  dev.off()
}



###create the boxplots using a loop
for(i in 17:ncol(logdata)){
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Boxplots/",colnames(logdata)[i],"_Boxplots_log.jpeg")
  jpeg(filename=file.name)
  boxplot(logdata[,i]~idata$group_adm, xlab = colnames(logdata)[i], main = paste0("Boxplots",colnames(logdata)[i]))
  dev.off()
}








#########################################
#Boxcox Transform
colnames(idata)
trans_setup<-caret::preProcess(idata[,17:ncol(idata)], method=c("BoxCox"))
idata_trans<-predict(trans_setup,idata[,17:ncol(idata)])

idata_trans<-cbind(idata[,c(1:16)],idata_trans)
colnames(idata_trans)



###create the boxplots using a loop
for(i in 17:ncol(idata_trans)){
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Histograms/",colnames(idata_trans)[i],"_Histograms_BoxCox.jpeg")
  jpeg(filename=file.name)
  hist(idata_trans[,i], xlab = colnames(idata_trans)[i], main = paste0("Histogram",colnames(idata_trans)[i]))
  dev.off()
}



outlier_all_SD<-c()
##### Hampel filter
for(i in 17:ncol(idata)){
  
  lower_bound <- signif( median(idata[,i], na.rm=T) - 15 * mad(idata[,i], na.rm=T) , 4)
  upper_bound <- signif( median(idata[,i], na.rm=T) + 15 * mad(idata[,i], na.rm=T) , 4)
  outlier_under <- length(which(idata[,i] < lower_bound))
  outlier_over <- length(which(idata[,i] > upper_bound))
  
  outlier__SD <- cbind("metabolite"= colnames(idata)[i], 
                                 "Lower"=lower_bound, "Upper"=upper_bound, 
                                 "Outlier_under"=outlier_under, "Outlier_over"=outlier_over)
         
  outlier_all_SD <-rbind(outlier_all_SD, outlier__SD)
}


outlier_all_SD<-data.frame(outlier_all_SD)
outlier_all_SD


write.csv(outlier_all_SD, file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\QC_cleaning\\CHAIN_mNutR_outliers_Hampel_all_15SD_",Sys.Date(),".csv"))




###################################
### [Note] Grubbs test depends of having some form of normal distibution, so using BoxCox transformed data
# test <- grubbs.test(idata$albumin_adm)
# test
# 
# 
# ### apply the function to all metabolites
# outlier_out <- apply(idata_trans[c(17:ncol(idata_trans))], 2,grubbs.test )
# outlier_out
# 
# outlier_compile<-c()
# for (i in 1:183) {
#   outlier_compile<-rbind(outlier_compile,
#                          c( names(outlier_out[i]), outlier_out[[i]]$p.value, outlier_out[[i]]$alternative))
# }
# 
# outlier_compile
# 
# 
# write.csv(outlier_compile, file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\QC_cleaning\\CHAIN_mNutR_outliers_GrubsTest_all_",Sys.Date(),".csv"))




########################################
#######outliers  Removed

colnames(idata)
### a-aminoadipic acid, n=2
max(idata$a_Aminoadipic_ac, na.rm=T)
idata$a_Aminoadipic_ac [idata$a_Aminoadipic_ac == max(idata$a_Aminoadipic_ac, na.rm=T)]<-NA
max(idata$a_Aminoadipic_ac, na.rm=T)
idata$a_Aminoadipic_ac [idata$a_Aminoadipic_ac == max(idata$a_Aminoadipic_ac, na.rm=T)]<-NA
max(idata$a_Aminoadipic_ac, na.rm=T)
idata$a_Aminoadipic_ac

### butyric acid, n=1
idata$Butyric_ac[idata$Butyric_ac == max(idata$Butyric_ac, na.rm=T)]<-NA
idata$Butyric_ac

### C5, n=1
idata$C5 [idata$C5 == max(idata$C5, na.rm=T)]<-NA
idata$C5

### g-aminobutyric acid, n=1
idata$g_aminobutyric_ac [idata$g_aminobutyric_ac == max(idata$g_aminobutyric_ac, na.rm=T)]<-NA
idata$g_aminobutyric_ac

### Hippuric acid, n=1
idata$Hippuric_ac [idata$Hippuric_ac == max(idata$Hippuric_ac, na.rm=T)]<-NA
#idata$Hippuric_ac [idata$Hippuric_ac == max(idata$Hippuric_ac, na.rm=T)]<-NA
idata$Hippuric_ac


### Homovanillic acid, n=2
idata$Homovanillic_ac [idata$Homovanillic_ac == max(idata$Homovanillic_ac, na.rm=T)]<-NA
idata$Homovanillic_ac [idata$Homovanillic_ac == max(idata$Homovanillic_ac, na.rm=T)]<-NA
idata$Homovanillic_ac


### Kynurenine_Trp
idata$Kynurenine_Trp [idata$Kynurenine_Trp == max(idata$Kynurenine_Trp, na.rm=T)]<-NA
idata$Kynurenine_Trp



### LysoC26.1, n=1
idata$LYSOC26.1 [idata$LYSOC26.1 == max(idata$LYSOC26.1, na.rm=T)]<-NA
idata$LYSOC26.1

### LysoC24.0, n=1
idata$LYSOC24.0 [idata$LYSOC24.0 == max(idata$LYSOC24.0, na.rm=T)]<-NA
idata$LYSOC24.0


### LysoC28.1, n=1
idata$LYSOC28.1 [idata$LYSOC28.1 == max(idata$LYSOC28.1, na.rm=T)]<-NA
idata$LYSOC28.1


###p.Hydroxyphenylacetic_ac
idata$p.Hydroxyphenylacetic_ac [idata$p.Hydroxyphenylacetic_ac == max(idata$p.Hydroxyphenylacetic_ac, na.rm=T)]<-NA
idata$p.Hydroxyphenylacetic_ac



### PC40.2_AA, n=1
idata$PC40.2_AA [idata$PC40.2_AA == max(idata$PC40.2_AA, na.rm=T)]<-NA
idata$PC40.2_AA


### PC40.1, n=1
#idata$PC40.1 [idata$PC40.1 == max(idata$PC40.1, na.rm=T)]<-NA
#idata$PC40.1 [idata$PC40.1 == max(idata$PC40.1, na.rm=T)]<-NA
#idata$PC40.1


### Succinic_ac
idata$Succinic_ac [idata$Succinic_ac == max(idata$Succinic_ac, na.rm=T)]<-NA
idata$Succinic_ac

### Valeric acid, n=1
idata$Valeric_ac [idata$Valeric_ac == max(idata$Valeric_ac, na.rm=T)]<-NA
#idata$Valeric_ac [idata$Valeric_ac == max(idata$Valeric_ac, na.rm=T)]<-NA
idata$Valeric_ac


### metals checked
colnames(idata) 
### removing outlier row 52, based on PCA
idata[52:54,178:200]
idata[52,178:200]<-NA
idata[52:54,178:200]



PCA(idata[,c("albumin_adm","Retinol", "OH25_VitD3", "aTocopherol", "bg_Tocopherol", "B1", 
               "B2","B2_corr", "B3_amide", "B5", "B6","B6_corr", "B7", "C")], scale.unit = TRUE)


PCA(idata[,c("Creatinine", "Glycine", 
               "Alanine", "Serine", "Proline", "Valine", "Threonine", "Taurine", 
               "Putrescine", "trans.Hydroxyproline", "Leucine", "Isoleucine", 
               "Asparagine", "Aspartic_ac", "Glutamine", "Glutamic_ac", "Methionine", 
               "Histidine", "a_Aminoadipic_ac", "Phenylalanine", "Methionine.sulfoxide", 
               "Arginine", "Acetyl.ornithine", "Citrulline", "Serotonin", "Tyrosine", 
               "Asymmetric_dimethylarginine", "total_dimethylarginine", "Tryptophan", 
               "Kynurenine", "Carnosine", "Ornithine", "Lysine", "Spermidine", 
               "Spermine", "Sarcosine", "Creatine", "Betaine", "Choline", "Trimethylamine_N.oxide", 
               "Methylhistidine", "Homocysteine", "Hypoxanthine", "Ethanolamine", 
               "Cystathionine", "N.acetylpetrescine", "Homocitrulline", "Urea", 
               "Hydroxy.lysine", "a_aminobutyric_ac", "g_aminobutyric_ac", "Lactic_ac", 
               "b_Hydroxybutyric_ac", "a_Ketoglutaric_ac", "Citric_ac", "Butyric_ac", 
               "Propionic_ac", "HPHPA", "p.Hydroxyhippuric_ac", "Succinic_ac", 
               "Fumaric_ac", "Pyruvic_ac", "Isobutyric_ac", "Hippuric_ac", "Methylmalonic_ac", 
               "Homovanillic_ac", "Indole_acetic_ac", "Uric_ac", "OH5_Indoleacetic_ac", 
               "p.Hydroxyphenylacetic_ac", "Valeric_ac", "Uridine", "Xanthine", 
               "Guanidineacetic_ac", "Acetyl.lysine", "Glucose")], scale.unit = TRUE)





PCA(idata[,c("LYSOC14.0", 
             "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", "LYSOC18.1", 
             "LYSOC18.0", "LYSOC20.4", "LYSOC20.3", "LYSOC24.0", "LYSOC26.1", 
             "LYSOC26.0", "LYSOC28.1", "LYSOC28.0", "Sphg_14.1SMOH", "Sphg_16.1SM", 
             "Sphg_16.0SM", "Sphg_16.1SMOH", "Sphg_18.1SM", "PC32.2_AA", "Sphg_18.0SM", 
             "Sphg_20.2SM", "PC36.0_AE", "PC36.6_AA", "PC36.0_AA", "Sphg_22.2SMOH", 
             "Sphg_22.1SMOH", "PC38.6_AA", "PC38.0_AA", "PC40.6_AE", "Sphg_24.1SMOH", 
             "PC40.6_AA", "PC40.2_AA", "PC40.1_AA")], scale.unit = TRUE)



PCA(idata[,c("C0", "C2", "C3.1", "C3", 
               "C4.1", "C4", "C3OH", "C5.1", "C5", "C4OH", "C6.1", "C6", "C5OH", 
               "C5.1DC", "C5DC", "C8", "C5MDC", "C9", "C7DC", "C10.2", "C10.1", 
               "C10", "C12.1", "C12", "C14.2", "C14.1", "C14", "C12DC", "C14.2OH", 
               "C14.1OH", "C16.2", "C16.1", "C16", "C16.2OH", "C16.1OH", "C16OH", 
               "C18.2", "C18.1", "C18", "C18.1OH")], scale.unit = TRUE)




PCA(idata[,c("Barium", "Sodium", "Magnesium", 
               "Phosphorous", "Calcium", "Titanium", "Vanadium", "Manganese", 
               "Iron", "Cobalt", "Copper", "Zinc", "Gallium", "Arsenic", "Selenium", 
               "Rubidium", "Strontium", "Antimony", "Tellurium")], scale.unit = TRUE)


colnames(idata)
# final clean data
write.csv(idata, file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\5-Data\\CHAIN_mNutR_fulldata_X_",Sys.Date(),".csv"))




colnames(idata)
trans_setup<-caret::preProcess(idata[,17:ncol(idata)], method=c("BoxCox"))
idata_trans<-predict(trans_setup,idata[,17:ncol(idata)])

idata_trans<-cbind(idata[,c(1:16)],idata_trans)
colnames(idata_trans)

# final clean data
write.csv(idata, file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\5-Data\\CHAIN_mNutR_fulldata_X_BoxCox_",Sys.Date(),".csv"))

