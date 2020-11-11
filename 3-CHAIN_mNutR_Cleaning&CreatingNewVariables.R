
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
idata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_idata_2020-11-10.csv"), row.names = 1) # ND (not done) equivalent to NA


### list factors to convert
str(idata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)



### make site human readible
idata <- idata %>% mutate(site = derivedFactor(
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

summary(idata$site)
idata$site<-droplevels(idata$site)



###############################
### make groups human readible
idata <- idata %>% mutate(group_adm = derivedFactor(
  "SM" = group_adm=="1",
  "MAM" = group_adm=="2",
  "NAM" = group_adm=="3",
  "CP" = group_adm=="4",
  .method = "first",
  .default = "NA"))

summary(idata$group_adm)
idata$group_adm<-droplevels(idata$group_adm)



###### check summary
mycontrols  <- tableby.control(test=FALSE, total=TRUE, digits.pct = 0
                               #numeric.test="kwt", cat.test="chisq",
                               #numeric.stats=c("N", "median", "q1q3"),
                               #cat.stats=c("countpct"),
                               #stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3')
                               )

labels(idata)  <- c(site = 'Site', group_adm = "Group")
tab1 <- tableby(site ~ group_adm, data=idata, control=mycontrols)
summary(tab1, text=TRUE) 

### write table
write2html(tab1, file=paste0("Summary_Tally_SamplesPerSite_",Sys.Date(),".html"))





### make sex human readible
idata <- idata %>% mutate(sex_adm = derivedFactor(
  "Male" = sex_adm=="1",
  "Female" = sex_adm=="2",
  .method = "first",
  .default = "NA"))

summary(idata$sex_adm)



#### replace LOD by 0.0104/2  half lowest detection range
idata$Carnosine<-gsub("< LOD", 0.0052, idata$Carnosine)
idata$Carnosine<-as.numeric(idata$Carnosine)



#### correcting units for albumin
hist(idata$albumin_adm)
idata[which(idata$albumin_adm < 10),]

#Blantyre<-idata[idata$site =="Blantyre", ]
#write.csv(Blantyre, file="CHAIN_admBlood_albumin_unit_check.csv")



###### replace after changing units
idata[which(idata$albumin_adm < 10),]$albumin_adm <- idata[which(idata$albumin_adm < 10),]$albumin_adm * 10
idata[which(idata$albumin_adm < 10),]$albumin_adm

#### correction B2
idata$B2_corr<-signif(idata$B2/idata$albumin_adm, 4)
idata$B2_corr


#### correction B6
idata$B6_corr<-signif(idata$B6/idata$albumin_adm, 4)
idata$B6_corr


head(idata)[1:20]



######################### calculating other summary ratios


###create AAA aromatic amino acids with aromatic ring
### tyrosine,tryptophane, phenylalanine, histidine, + 5-hydroxytryptophan, thyroxine, L-DOPA
colnames(idata)
idata$AAA <- idata$Tyrosine + idata$Tryptophan + idata$Phenylalanine + idata$Histidine
hist(idata$AAA)


### BCAA branched amino acids 
# essential (isoleucine, leucine, and valine)  + 2-aminoisobutyric acid
colnames(idata)
idata$BCAA<- idata$Isoleucine + idata$Leucine + idata$Valine
hist(idata$BCAA)


### BCAA/AAA ratio
colnames(idata)
idata$BCAA_AAA<- idata$BCAA / idata$AAA
hist(idata$BCAA_AAA)


### Fischer ratio 
### essential (isoleucine, leucine, and valine) to the AAA tyrosine and phenylalanine
colnames(idata)
idata$Fischer_ratio <- idata$BCAA / (idata$Tyrosine + idata$Phenylalanine)
hist(idata$Fischer_ratio)


#### Essential amino acids 
### phenylalanine, valine, threonine, tryptophan, methionine, leucine, isoleucine, lysine, and histidine
colnames(idata)
idata$Essential_AA <- idata$Threonine + idata$Tryptophan + idata$Methionine + idata$Phenylalanine + 
  idata$Histidine + idata$Isoleucine + idata$Leucine + idata$Valine + idata$Lysine
hist(idata$Essential_AA)


#### glucogenic amino acids
### Alanine, Arginine, Asparagine, Aspartic,Glutamic,Glutamine,Glycine,Histidine,Methionine,Proline,Serine,Valine
### missing --> Cysteine
colnames(idata)
idata %>% select(contains("Aspartic")) %>% colnames()

idata$Glucogenic_AA<-idata$Alanine + idata$Arginine + idata$Asparagine + idata$Aspartic_ac + 
  idata$Glutamic_ac + idata$Glutamine + idata$Glycine + idata$Histidine + 
  idata$Methionine + idata$Proline + idata$Serine + idata$Valine 
hist(idata$Glucogenic_AA)

### ketogenic 
## leucine lysine
colnames(idata)
idata$Keto_AA<- idata$Leucine + idata$Lysine
hist(idata$Keto_AA)


##### Kynurenine_Trp
colnames(idata)
idata$Kynurenine_Trp<-idata$Kynurenine/idata$Tryptophan
hist(idata$Kynurenine_Trp)


###
colnames(idata)
idata$Orn_Arg<- idata$Ornithine / idata$Arginine
hist(idata$Orn_Arg)


###Putrescine_Orn
colnames(idata)
idata$Putrescine_Orn<- idata$Ornithine / idata$Putrescine
hist(idata$Putrescine_Orn)


##### Serotonin_Trp
colnames(idata)
idata$Serotonin_Trp<- idata$Serotonin / idata$Tryptophan
hist(idata$Serotonin_Trp)

### Spermidine_Putrescine
colnames(idata)
idata$Spermidine_Putrescine<- idata$Spermidine / idata$Putrescine
hist(idata$Spermidine_Putrescine)


#### Spermine_Spermidine
colnames(idata)
idata$Spermine_Spermidine<- idata$Spermine/idata$Spermidine
hist(idata$Spermine_Spermidine)
hist(idata$Spermine)
hist(idata$Spermidine)

### Tyr_Phe
colnames(idata)
idata$Tyr_Phe<- idata$Tyrosine / idata$Phenylalanine
hist(idata$Tyr_Phe)


#### Total_DMA_Arg
colnames(idata)
idata$Total_DMA_Arg<- idata$total_dimethylarginine / idata$Arginine
hist(idata$Total_DMA_Arg)


#### Cit_Orn
colnames(idata)
idata$Cit_Orn<- idata$Citrulline / idata$Ornithine
hist(idata$Cit_Orn)

### Cit_Arg
colnames(idata)
idata$Cit_Arg<- idata$Citrulline / idata$Arginine
hist(idata$Cit_Arg)


#### C2_C0
colnames(idata)
idata$C2_C0<- idata$C2 / idata$C0
hist(idata$C2_C0)


#### total amino acids
colnames(idata)
idata$Total_AAs<- idata$Arginine + idata$Methionine + idata$Valine + idata$Alanine + idata$Asparagine +
  idata$Aspartic_ac + idata$Glutamine + idata$Glutamic_ac + idata$Glycine + idata$Histidine +
  idata$Isoleucine + idata$Leucine + idata$Lysine + idata$Phenylalanine + idata$Proline + 
  idata$Serine + idata$Threonine + idata$Tryptophan + idata$Tyrosine
#missing: cysteine
hist(idata$Total_AAs)


### urea cycle amino acids
### Citrulline, l-ornitine, L-argining, urea, aspartic acid
colnames(idata)
idata$urea_cycle<-idata$Ornithine + idata$Arginine + idata$Citrulline + idata$Aspartic_ac
hist(idata$urea_cycle)
hist(idata$Ornithine)
hist(idata$Arginine)
hist(idata$Citrulline)
hist(idata$Aspartic_ac)




choose.dir()
write.csv(idata, file=paste0("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\5-Data\\CHAIN_mNutR_fulldata_", Sys.Date(),".csv"))



