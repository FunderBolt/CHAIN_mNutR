
#############################
###  CHAIN_NutR
###  Making Tables
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","arsenal","gt","multcomp")



# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)

##set overall theme
#theme_set(theme_minimal())



# load full data
idata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_fulldata_X_2020-11-21.csv"), row.names = 1) 



### list factors to convert
str(idata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)




colnames(idata)
# How many ABC groups in the datase
summary(idata$group_adm) # equal proportion of community and SAM
###### Making Participant characteristic table
labels(idata)  <- c(site = 'Site', group_adm = "Group")


mycontrols  <- tableby.control(test=TRUE, digits.p = 3,total=FALSE, digits.pct = 0, cat.simplify =F,
                               numeric.test="anova", cat.test="fe",
                               numeric.stats=c("meansd"), digits = 1
                               #cat.stats=c("countpct"),
                               #stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3')
)


tab1 <- tableby(group_adm ~ sex_adm + age_adm + oedema_adm + muac_adm + haz_adm + waz_adm + whz_adm , data=idata, control=mycontrols)
summary(tab1, text=TRUE) 

### write table
write2html(tab1, file=paste0("Table1_Participant_Characteristic_", Sys.Date(),".html"))







#######################################################
###### Making metabolite tables


####pull the data you want to include in table
colnames(idata)
list_met <- colnames(idata)[-c(1:16)]
list_met

#####split dataset by diagnositic, creating 3 datasets for each
summary(idata$group_adm)
SM<-subset(idata,idata$group_adm =="SM")
CP<-subset(idata,idata$group_adm =="CP")


table.data<-SM[list_met]
str(table.data)

cal.iqr<-function(x){quantile(table.data[,x], na.rm=TRUE)}
iqr.list<-sapply(1:ncol(table.data),cal.iqr)
iqr.list
iqr.list<-signif(iqr.list,3)
iqr.list

iqr.list<-t(iqr.list)
iqr.list

medians.iqrs<-paste0(iqr.list[,3]," [",iqr.list[,2],"-",iqr.list[,4],"]")
medians.iqrs


medians.iqrs.SM<-medians.iqrs




########################### calc table for CP
table.data<-CP[list_met]
str(table.data)

cal.iqr<-function(x){quantile(table.data[,x], na.rm=TRUE)}
iqr.list<-sapply(1:ncol(table.data),cal.iqr)
iqr.list
iqr.list<-signif(iqr.list,3)
iqr.list

iqr.list<-t(iqr.list)
iqr.list

medians.iqrs<-paste0(iqr.list[,3]," [",iqr.list[,2],"-",iqr.list[,4],"]")
medians.iqrs


medians.iqrs.CP<-medians.iqrs


table<-cbind(medians.iqrs.SM,medians.iqrs.CP)
rownames(table)<-colnames(table.data)
head(table)



table<-data.frame(cbind(row.names(table),table))
colnames(table) <-c("Micronutrient / Metabolite" , "SM", "CP")


### clean out column names
#dput(table[,1])
table %>% gt()



#write.csv(table, file=paste0("Table2_Median_IQR_Metabolites_",Sys.Date(),".csv"))









######################################
#### adding GLMs

# load full data
#idata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_fulldata_X_BoxCox_2020-11-08.csv"), row.names = 1) 
idata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_fulldata_X_2020-11-21.csv"), row.names = 1) 
colnames(idata)


### list factors to convert
str(idata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)




colnames(idata)
# How many ABC groups in the datase
summary(idata$group_adm) # equal proportion of community and SAM
###### Making Participant characteristic table
labels(idata)  <- c(site = 'Site', group_adm = "Group")

# idata$Group <- relevel(idata$group_adm, ref = "CP")



#idata[,17:ncol(idata)]<-idata[,17:ncol(idata)]+1
#################################################
###GLMs for all metabolites --- create a function

#### exclude lead (giving issues, to many zeros)


### M1
# GLM.run<-function(y) {
#   form <- as.formula(paste0(y,"~ group_adm"))
#   fit<-(glm(form, data=idata, na.action = na.exclude))
# }

GLM.run<-function(y) {
  form <- as.formula(paste0("group_adm~", y))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}




#### M2
# GLM.run<-function(y) {
#   form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm"))
#   fit<-(glm(form, data=idata, na.action = na.exclude))
# }

GLM.run<-function(y) {
  form <- as.formula(paste0("group_adm~", y," + sex_adm + age_adm"))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}



#### M3
# GLM.run<-function(y) {
#   form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
#   fit<-(glm(form, data=idata, na.action = na.exclude))
# }


GLM.run<-function(y) {
  form <- as.formula(paste0("group_adm~", y," + sex_adm + age_adm + site"))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}


#### M4
GLM.run<-function(y) {
  form <- as.formula(paste0("group_adm~",y,"+ site"))
  fit<-(glm(form, data=idata, family=binomial, na.action = na.exclude))
}






####################################################################
colnames(idata)
### apply the function to all metabolites
GLMs.out <- lapply(colnames(idata)[c(17:ncol(idata))],GLM.run )


#### print out results
results<-lapply(GLMs.out, function(x){summary(x)})
results


#Pull coefficients from models
estim.coef.results<-sapply(results, function(x){coef(x)})
estim.coef.results


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



#ALL_estim.coef.results<-estim.coef.results
ALL_estim.coef.results<-data.frame(cbind(ALL_estim.coef.results,  estim.coef.results))




ALL_estim.coef.results<-data.frame(ALL_estim.coef.results)


colnames(ALL_estim.coef.results)
row.names(ALL_estim.coef.results)
### add in FDR corrected p-values
head(ALL_estim.coef.results)
# p.vals<-c(ALL_estim.coef.results$p.group_admSM)
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




###### prep for saving
ALL_estim.coef.results<-apply(ALL_estim.coef.results,2, signif, digits=3)
head(ALL_estim.coef.results)

#write.csv(ALL_estim.coef.results,file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\BoxCoxTransform_coefficients_all_M1M2&M3_Gamma_3digits_",Sys.Date(), ".csv"))
write.csv(ALL_estim.coef.results,file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\Coefficients_all_M1M2&M3M4_binomial_3digits_",Sys.Date(), ".csv"))



################################################
### collecting median IQR and GLM data
colnames(table)
head(table)

colnames(ALL_estim.coef.results)
table_2<-merge(table, ALL_estim.coef.results[,c("est.albumin_adm.3","se.albumin_adm.3","p.albumin_adm.3","FDR_p.group_admSM_M3")], by="row.names",sort = FALSE)
head(table_2)
colnames(table_2)<-c("Row.names", "Micronutrient / Metabolite", "SM", "CP", "b", 
                     "SE", "p", "q")


write.csv(table_2[,-1],file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\MedianIQR_GLM_results_all_M3_binomial_3digits_",Sys.Date(), ".csv"))

#################################################
#################################################


###Residual Plots
Residual_plot.run <-function(y) {
  form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
  fit<-glm(form, data=idata, na.action = na.exclude)
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Residuals/", y,".GLM_BoxCox_M3_residual_plot.jpeg")
  jpeg(filename=file.name)
  print(plot(fit,which=1))
  dev.off()
}


###Residual Plots
Residual_plot.run <-function(y) {
  form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
  fit<-glm(form, data=idata, family="Gamma", na.action = na.exclude)
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/Residuals_gamma/", y,".GLM_BoxCox_M3_residual_plot.jpeg")
  jpeg(filename=file.name)
  print(plot(fit, which=1))
  dev.off()
}



lapply(colnames(idata)[17:ncol(idata)],Residual_plot.run )


###QQ-plots
qqplots.run <-function(y) {
  form <- as.formula(paste0(y,"~ group_adm + sex_adm + age_adm + site"))
  fit<-glm(form, data=idata, na.action = na.exclude)
  file.name<-paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/qqPlot/", y,".GLM_BoxCox_M3_qqplot.jpeg")
  jpeg(filename=file.name)
  print(plot(fit,which=2))
  dev.off()
}

 lapply(colnames(idata)[c(17:ncol(idata))],qqplots.run )





 ############################# Kruskall wallis test
 
 # ####
 # GLM.run<-function(y) {
 #   form <- as.formula(paste0(y,"~ group_adm"))
 #   fit<-(kruskal.test(form, data=idata, na.action = na.exclude))
 # }
 # 
 # 
 # 
 # for(i in 1:length(GLMs.out))
 #   results[1]<-GLMs.out[[1]]$p.value
 # 

