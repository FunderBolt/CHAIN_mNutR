
# ====  CHAIN_NutR | Participant Characteristic Table =======================================================


# 0. Packages ============================================================
# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","arsenal")


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



# 3. Arsenal | Participant Table ============================================================

colnames(idata)
summary(idata$group_adm) # equal proportion of community and SAM


# * 1. Set labels for column names  ============================================================

labels(idata)  <- c(site = 'Site', group_adm = "Enrolment group", 
                    sex_adm="Sex", age_adm="Age, mo", 
                    oedema_adm="Oedema status",
                    muac_adm="MUAC, mm", haz_adm="Height-for-age, z-score",
                    waz_adm="Weight-for-age, z-score",
                    whz_adm="Weight-for-height, z-score")

# * 2. Set control parameters for table  ============================================================

mycontrols <- tableby.control(test=TRUE, digits.p = 3, 
                              total=FALSE, 
                              digits.pct = 0, cat.simplify =F, cat.test="fe",  
                              #cat.stats=c("countpct"),
                              numeric.test="anova", numeric.stats=c("meansd"), digits = 1
                              #stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3')
)


# * 3. Run table  ============================================================

tab1 <- tableby(group_adm ~ sex_adm + age_adm + oedema_adm + muac_adm + haz_adm + 
                  waz_adm + whz_adm , data=idata, control=mycontrols)
summary(tab1, text=TRUE) 


# * 4. Write Table | Participant Characteristics =================================================

### write table
write2html(tab1, file=paste0(dirname(dirname(here())), "/7-Analysis-Results/ParticipantCharacteristic_Tables/Table1_Participant_Characteristic_", Sys.Date(),".html"))
write2word(tab1, file=paste0(dirname(dirname(here())), "/7-Analysis-Results/ParticipantCharacteristic_Tables/Table1_Participant_Characteristic_", Sys.Date(),".doc"))




