
#############################
###  CHAIN_NutR
###  Cleaning and Creating new variables
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","mosaic","arsenal")


# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)



# load full data
fulldata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_fulldata_2020-11-06.csv"), row.names = 1) # ND (not done) equivalent to NA



### list factors to convert
str(fulldata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

fulldata[, factor_list] <- lapply(fulldata[, factor_list], as.factor)
str(fulldata)



### make site human readible
fulldata <- fulldata %>% mutate(site = derivedFactor(
  "Banfora" = site=="60001",
  "Blantyre" = site=="30001",
  "Dhaka" = site=="50001",
  "Kampala" = site=="20001",
  "Karachi" = site=="40001",
  "Kilifi" = site=="10001",
  "Matlab" = site=="50002",
  "Migori" = site=="10003",
  "Nairobi" = site=="10002",
  .method = "first",
  .default = "NA"))

summary(fulldata$site)
fulldata$site<-droplevels(fulldata$site)



###############################
### make groups human readible
fulldata <- fulldata %>% mutate(group_adm = derivedFactor(
  "SM" = group_adm=="1",
  "MAM" = group_adm=="2",
  "NAM" = group_adm=="3",
  "CP" = group_adm=="4",
  .method = "first",
  .default = "NA"))

summary(fulldata$group_adm)
fulldata$group_adm<-droplevels(fulldata$group_adm)



###### check summary
mycontrols  <- tableby.control(test=FALSE, total=TRUE, digits.pct = 0
                               #numeric.test="kwt", cat.test="chisq",
                               #numeric.stats=c("N", "median", "q1q3"),
                               #cat.stats=c("countpct"),
                               #stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3')
                               )

labels(fulldata)  <- c(site = 'Site', group_adm = "Group")
tab1 <- tableby(site ~ group_adm, data=fulldata, control=mycontrols)
summary(tab1, text=TRUE) 

### write table
write2html(tab1, file=paste0("Summary_Tally_SamplesPerSite_",Sys.Date(),".html"))





### make sex human readible
fulldata <- fulldata %>% mutate(sex_adm = derivedFactor(
  "Male" = sex_adm=="1",
  "Female" = sex_adm=="2",
  .method = "first",
  .default = "NA"))

summary(fulldata$sex_adm)




#### replace LOD by NA
fulldata$Carnosine<-gsub("< LOD", 0, fulldata$Carnosine)
fulldata$Carnosine<-as.numeric(fulldata$Carnosine)




head(fulldata)[1:20]


