# ====  CHAIN_NutR | Making metabolite Table, Median IQR and GLM results ===================================


# 0. Packages ============================================================
# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","arsenal","gt","multcomp")
install.packages("tidymodels")

# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)

##set overall theme
#theme_set(theme_minimal())


# 1. Load data ============================================================

# load full data
idata <- read.csv(paste0(dirname(dirname(here())), "/5-Data","/CHAIN_mNutR_fulldata_X_2020-11-21.csv"), row.names = 1) 


# 2. Set Factor ============================================================
### list factors to convert
str(idata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)


# 3. Function to calculate IQR across columns =================================================

cal.iqr<-function(x){quantile(table.data[,x], na.rm=TRUE)}



# 4. Create tables median IQR for all metabolites of interest within sub-groups =======================================

# * 1. Pull metabolites of interest and subset tables for sub-groups =======================================

colnames(idata)
list_met <- colnames(idata)[-c(1:16)]
list_met

#####split dataset by diagnositic, creating 3 datasets for each
summary(idata$group_adm)
SM<-subset(idata,idata$group_adm =="SM")
CP<-subset(idata,idata$group_adm =="CP")



# * 2. SM group: Calculate IQRS  ============================================================

### Subset SM table to include metabolites of interest 
table.data<-SM[list_met]
str(table.data)

### apply function to each column of table
iqr.list<-sapply(1:ncol(table.data),cal.iqr)
### round off to 3 significant values
iqr.list<-signif(iqr.list,3)
iqr.list

## transpose table
iqr.list<-t(iqr.list)
iqr.list

## create vector of median IQR values for SM group
medians.iqrs.SM<-paste0(iqr.list[,3]," [",iqr.list[,2],"-",iqr.list[,4],"]")
medians.iqrs.SM


# * 3. CP group: Calculate IQRS  ============================================================

### Subset CP table to include metabolites of interest 
table.data<-CP[list_met]
str(table.data)

### apply function to each column of table
iqr.list<-sapply(1:ncol(table.data),cal.iqr)
### round off to 3 significant values
iqr.list<-signif(iqr.list,3)
iqr.list

## transpose table
iqr.list<-t(iqr.list)
iqr.list

## create vector of median IQR values for SM group
medians.iqrs.CP<-paste0(iqr.list[,3]," [",iqr.list[,2],"-",iqr.list[,4],"]")
medians.iqrs.CP



# * 4. Collate into table  ============================================================

table<-data.frame(cbind(medians.iqrs.SM,medians.iqrs.CP))
rownames(table)<-list_met
head(table)

### set column names
colnames(table) <-c("SM", "CP")


### clean out column names
#dput(table[,1])
table %>% gt()


# * 5. Save table  ============================================================

write.csv(table, file=paste0("Table2_Median_IQR_Metabolites_",Sys.Date(),".csv"))




# 5. Functions: Set up GLMs for all metabolites  ============================================================



# * 1. Gaussian: Group adm only  ============================================================

# GLM.run1<-function(y) {
#   form <- as.formula(paste0(y,"~ group_adm"))
#   fit<-(glm(form, data=idata, na.action = na.exclude))
# }


# * 2. Binomial: Group adm only  ============================================================


GLM.run2<-function(y) {
  form <- as.formula(paste0("group_adm~", y))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}



# * 3. Gaussian: Group adm, adjusted sex, age ============================================================

# GLM.run3<-function(y) {
#   form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm"))
#   fit<-(glm(form, data=idata, na.action = na.exclude))
# }


# * 4. Binomial: Group adm, adjusted sex, age ============================================================

GLM.run4<-function(y) {
  form <- as.formula(paste0("group_adm~", y," + sex_adm + age_adm"))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}


# * 5. Gaussian: Group adm, adjusted sex, age, site ============================================================

# GLM.run5<-function(y) {
#   form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
#   fit<-(glm(form, data=idata, na.action = na.exclude))
# }


# * 6. Binomial: Group adm, adjusted sex, age, site ============================================================

GLM.run6<-function(y) {
  form <- as.formula(paste0("group_adm~", y," + sex_adm + age_adm + site"))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}


# * 7. Binomial: Group adm, site ============================================================

GLM.run7<-function(y) {
  form <- as.formula(paste0("group_adm~",y,"+ site"))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}



# 6. GLMs: Apply functions & collect results ================================================

colnames(idata)
### apply the function to all metabolites
GLMs.out <- lapply(colnames(idata)[c(17:ncol(idata))], GLM.run6)


#### print out results
results<-lapply(GLMs.out, function(x){summary(x)})
results


#Pull coefficients from models
estim.coef.results<-sapply(results, function(x){coef(x)})
estim.coef.results


# 7. GLMs: collate into dataframe ================================================

## transpose the table
estim.coef.results<-t(estim.coef.results)
estim.coef.results

#the output of above commands create a matrix, but dataframes are better to indicate column names and rownames
class(estim.coef.results)
estim.coef.results<-as.data.frame(estim.coef.results)

#add rownames for identifiers
rownames(estim.coef.results)<-colnames(idata)[c(17:ncol(idata))]
head(estim.coef.results)


#add columns names
results[1]
colnames(estim.coef.results)<-c(paste0("est.",rownames(coef(results[[1]]))),
                                paste0("se.",rownames(coef(results[[1]]))),
                                paste0("t.",rownames(coef(results[[1]]))),
                                paste0("p.",rownames(coef(results[[1]]))))
head(estim.coef.results)



### collect different models into a data frame
ALL_estim.coef.results<-estim.coef.results
ALL_estim.coef.results<-data.frame(cbind(ALL_estim.coef.results,  estim.coef.results))

### check expected structure
colnames(ALL_estim.coef.results)
row.names(ALL_estim.coef.results)



# 8. FDR correction: calculate FDR and add to table =====================================

### pull p-values for correction
head(ALL_estim.coef.results)
 p.vals<-c(ALL_estim.coef.results$p.group_admSM)
# p.vals<-c(ALL_estim.coef.results$p.group_admSM.1)
# p.vals<-c(ALL_estim.coef.results$p.group_admSM.2)


colnames(ALL_estim.coef.results)

#### add in FDR p-values
# apply correction -- > some correction option methods: "bonferroni", "fdr"
ALL_estim.coef.results$FDR_p.group_admSM_M1<-p.adjust(ALL_estim.coef.results$p.albumin_adm, method="fdr")
ALL_estim.coef.results$FDR_p.group_admSM_M2<-p.adjust(ALL_estim.coef.results$p.albumin_adm.1, method="fdr")
ALL_estim.coef.results$FDR_p.group_admSM_M3<-p.adjust(ALL_estim.coef.results$p.albumin_adm.2, method="fdr")
ALL_estim.coef.results$FDR_p.group_admSM_M4<-p.adjust(ALL_estim.coef.results$p.albumin_adm.3, method="fdr")
#dput(final)

head(ALL_estim.coef.results)
###


# 9. Save table: GLM results =====================================

### apply rounding to significant digits
ALL_estim.coef.results<-apply(ALL_estim.coef.results,2, signif, digits=3)
head(ALL_estim.coef.results)


#write.csv(ALL_estim.coef.results,file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\BoxCoxTransform_coefficients_all_M1M2&M3_Gamma_3digits_",Sys.Date(), ".csv"))
write.csv(ALL_estim.coef.results,file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\Coefficients_all_M1M2&M3M4_binomial_3digits_",Sys.Date(), ".csv"))



# 10. Collecting median IQR and GLM data =====================================

colnames(table)
head(table)

colnames(ALL_estim.coef.results)
table_2<-merge(table, ALL_estim.coef.results[,c("est.albumin_adm.3","se.albumin_adm.3","p.albumin_adm.3","FDR_p.group_admSM_M3")], by="row.names",sort = FALSE)
head(table_2)

### add in colnames
colnames(table_2)<-c("Row.names", "Micronutrient / Metabolite", "SM", "CP", "b", 
                     "SE", "p", "q")

write.csv(table_2[,-1],file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\MedianIQR_GLM_results_all_M3_binomial_3digits_",Sys.Date(), ".csv"))



# 11. Model checks | Residual plots for GLMs  =====================================


# * 1. Residual plots: Create loop function: create and save each plot ================================

### Residual Plots: Family Gaussian
Residual_plot.run <-function(y) {
  form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
  fit<-glm(form, data=idata, na.action = na.exclude)
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Residuals/", y,".GLM_BoxCox_M3_residual_plot.jpeg")
  jpeg(filename=file.name)
  print(plot(fit,which=1))
  dev.off()
}


### Residual Plots: Family Gamma
Residual_plot.run <-function(y) {
  form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
  fit<-glm(form, data=idata, family="Gamma", na.action = na.exclude)
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Residuals_gamma/", y,".GLM_BoxCox_M3_residual_plot.jpeg")
  jpeg(filename=file.name)
  print(plot(fit, which=1))
  dev.off()
}


### Apply function | create and save each residual plot ===============================
lapply(colnames(idata)[17:ncol(idata)],Residual_plot.run )



# * 2. QQ plots: Create loop function: create and save each plot ================================

###QQ-plots
qqplots.run <-function(y) {
  form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
  fit<-glm(form, data=idata, na.action = na.exclude)
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/qqPlot/", y,".GLM_BoxCox_M3_qqplot.jpeg")
  jpeg(filename=file.name)
  print(plot(fit,which=2))
  dev.off()
}


### Apply function | create and save each QQ plot ===============================
lapply(colnames(idata)[c(17:ncol(idata))],qqplots.run )


