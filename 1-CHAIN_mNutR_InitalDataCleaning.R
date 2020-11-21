#############################
###  CHAIN_NutR
###  Raw data load and inital clean
#############################

# make list of needed packages
list_packages <- c("readxl","here","mgsub","tidyverse")



# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)





######################## Fat soluble vitamins
##load raw data
#choose.files()
fatsol_vits<-read_excel(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/TMIC016R fat-soluble vitamin data.xlsx"),skip=12)

### clean up names
#dput(colnames(fatsol_vits))
colnames(fatsol_vits)<-c("record_id", "Retinol", "OH25_VitD3", "aTocopherol", "bg_Tocopherol")


### Remove "-S) record_id names
fatsol_vits$record_id <- mgsub(fatsol_vits$record_id, c("-S", "-A0-S"), c("",""))

## check if duplicates 
fatsol_vits$record_id[duplicated(fatsol_vits$record_id)]
### check general format
head(fatsol_vits)






######################## Metabolites TMIC PRIME DI_LC-MS_MS
#choose.files()
met<-read_excel(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/TMIC016R TMIC PRIME DI_LC-MS_MS data.xlsx"),skip=12)

### clean up names
#dput(colnames(met))
colnames(met)[1]<-c("record_id")
colnames(met)<-mgsub(colnames(met), c("gamma-","alpha-","beta-", " ","acid","5-Hydroxy","AA","AE"), c('g_','a_',"b_","_","ac", "OH5","_AA","_AE"))

### suffix infront of any names with numbers related to sphingomyelins
SMnames<-c(colnames(met)[grep(c("^1"),colnames(met))], colnames(met)[grep(c("^2"),colnames(met))] )
SMnames

### create new names string
SMnames_new <- paste0("Sphg_", SMnames)
SMnames_new

## replace in dataframe
colnames(met)[colnames(met) %in% SMnames] <- SMnames_new
colnames(met)


### Remove "-S) record_id names
met$record_id <- mgsub(met$record_id, c("-S", "-A0-S"), c("",""))
met$record_id





######################## water vits
#choose.files()
Water_vits<-read_excel(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/TMIC016R water-soluble vitamin data.xlsx"),skip=11)

### clean up names
#dput(colnames(Water_vits))
colnames(Water_vits)<-c("record_id", "B1", "B2", "B3_amide", "B5", "B6", "B7", "C")


### Remove "-S) record_id names
Water_vits$record_id <- mgsub(Water_vits$record_id, c("-S", "-A0-S"), c("",""))

## check if duplicates 
Water_vits$record_id[duplicated(Water_vits$record_id)]
### check general format
head(Water_vits)




######################## Metals
#choose.files()
metals<-read_excel(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/TMIC016R metal data.xlsx"),skip=13)

### clean up names
#dput(colnames(metals))
colnames(metals)[1]<-c("record_id")


### Remove "Plasma_  record_id names
metals$record_id <- mgsub(metals$record_id, c("Plasma_ "), c(""))

## check if duplicates 
metals$record_id[duplicated(metals$record_id)]
### check general format
head(metals)



####################### Merge files
data<-merge(fatsol_vits, Water_vits, by="record_id", all=TRUE)
data<-merge(data, met, by="record_id", all=TRUE)



####################### Correct ids
### List of ids that need to be corrected
list_ids<-c(30002751, 30002802,30002818,30002821,30002837,30002842)
list_ids_new<-gsub("30002", "30001", list_ids)

### set to new ids
data$record_id[data$record_id %in% list_ids] <- list_ids_new
data$record_id


data<-merge(data, metals, by="record_id", all=TRUE)

write.csv(data, file=paste0("CHAIN_mNutR_metabolite_data_",Sys.Date(),".csv"))
