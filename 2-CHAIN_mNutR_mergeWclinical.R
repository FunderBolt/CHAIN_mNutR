#############################
###  CHAIN_NutR
###  Merge with clinical data
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse")


# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)



# load met data
met_data <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_metabolite_data_2020-11-20.csv"), row.names = 1) 

# load clinical data
#choose.files()
anthro <-read.csv("C:\\Users\\cbour\\OneDrive\\Desktop\\2019_DataCuration\\CHAIN_Data\\CHAIN_Data_v9\\demographics\\anthropometry.csv",na.string=".")
dput(colnames(anthro))
anthro_data<-c("record_id","oedema_adm","muac_adm", "haz_adm", "waz_adm", "whz_adm")
anthro<-anthro[anthro_data]


demograph <-read.csv("C:\\Users\\cbour\\OneDrive\\Desktop\\2019_DataCuration\\CHAIN_Data\\CHAIN_Data_v9\\demographics\\demographics.csv",na.string=".")
dput(colnames(demograph))
demograph_data<-c("record_id","sex_adm", "age_adm", "parttype_adm", "site","africa", 
                  "group_adm", "age_group_adm", "urban")
demograph<-demograph[demograph_data]


outcome <- read.csv("C:\\Users\\cbour\\OneDrive\\Desktop\\2019_DataCuration\\CHAIN_Data\\CHAIN_Data_v9\\demographics\\outcome.csv", na.string=".")
dput(colnames(outcome))
outcome_data<-c("record_id","adm_dead", "dead")
outcome<-outcome[outcome_data]


adm_blood<- read.csv("C:\\Users\\cbour\\OneDrive\\Desktop\\2019_DataCuration\\CHAIN_Data\\CHAIN_Data_v9\\laboratory\\Biochemistry\\clean_adm_clinical_chem_data.csv", na.string=".")
head(adm_blood)
colnames(adm_blood)[1] <- "record_id"
blood_data<-c("record_id","albumin_adm")
adm_blood<-adm_blood[blood_data]
#adm_cbc <- read.csv("C:\\Users\\cbour\\OneDrive\\Desktop\\2019_DataCuration\\CHAIN_Data\\CHAIN_Data_v9\\laboratory\\Complete Blood Count\\clean_adm_cbc_data.csv")



# left join for clinical data 
idata<-merge(demograph, anthro, by="record_id", all=TRUE)
idata<-merge(idata, outcome, by="record_id", all=TRUE)
idata<-merge(idata, adm_blood, by="record_id", all=TRUE)
### merge and restrict to metabolite dataset
idata<-merge(idata, met_data, by="record_id", all.y=TRUE)



dput(colnames(idata))
# final clean data
write.csv(idata, file=paste0("CHAIN_mNutR_idata_", Sys.Date(),".csv"))


