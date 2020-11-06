
#############################
###  CHAIN_NutR
###  Cleaning and Creating new variables
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse")


# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)


### list factors to convert
str(fulldata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

fulldata[, factor_list] <- lapply(fulldata[, factor_list], as.factor)
str(fulldata)




#### replace LOD by NA
fulldata$Carnosine<-gsub("< LOD", 0, fulldata$Carnosine)
fulldata$Carnosine<-as.numeric(fulldata$Carnosine)


