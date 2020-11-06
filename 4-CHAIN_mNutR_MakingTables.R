
#############################
###  CHAIN_NutR
###  Making Tables
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","arsenal")


#"summarytools","ggthemes","tidyselect","lubridate","ggplotify",
# "vcd","reshape2","pracma",
# "naniar","visdat", "UpSetR", # exploring missing data
# "kableExtra","RColorBrewer", "htmltools", 
# "flextable","DT","zscorer","janitor")

# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)

##set overall theme
#theme_set(theme_minimal())



colnames(fulldata)
# How many ABC groups in the datase
summary(fulldata$group_adm) # equal proportion of community and SAM




###### Making Participant characteristic table
labels(fulldata)  <- c(site = 'Site', group_adm = "Group")


mycontrols  <- tableby.control(test=TRUE, digits.p = 3,total=FALSE, digits.pct = 0, cat.simplify =F,
                               numeric.test="anova", cat.test="fe",
                               numeric.stats=c("meansd"), digits = 1
                               #cat.stats=c("countpct"),
                               #stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3')
)


tab1 <- tableby(group_adm ~ sex_adm + age_adm + oedema_adm + muac_adm + haz_adm + waz_adm + whz_adm , data=fulldata, control=mycontrols)
summary(tab1, text=TRUE) 

### write table
write2html(tab1, file=paste0("Table1_Participant_Characteristic_",Sys.Date(),".html"))







#######################################################
###### Making metabolite tables


####pull the data you want to include in table
colnames(fulldata)
list_met <- colnames(fulldata)[-c(1:16)]

#####split dataset by diagnositic, creating 3 datasets for each
summary(fulldata$group_adm)
SM<-subset(fulldata,fulldata$group_adm =="SM")
CP<-subset(fulldata,fulldata$group_adm =="CP")



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



dput(table[,1])

table %>% gt()



write.csv(table, file=paste0("Table2_Median_IQR_Metabolites_",Sys.Date(),".csv"))





