# ====  CHAIN_NutR | Cleaning and Creating new variables ===============================

### author: CBD
### date: format(Sys.Date(), "%B %d, %Y")



# 0. Packages ============================================================
# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","mosaic","arsenal")


# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)




# 1. Load data ============================================================

# load full data
idata <- read.csv(paste0(dirname(dirname(here())), "/5-Data","/CHAIN_mNutR_idata_2020-11-20.csv"), row.names = 1) 



# 2. Set Factor ============================================================
### list factors to convert

str(idata)
dput(colnames(idata))
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)


# 3. Make human readible ============================================================

# * 1. Site ============================================================

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
summary(idata$site)


# * 2. Group ============================================================

idata <- idata %>% mutate(group_adm = derivedFactor(
  "SM" = group_adm=="1",
  "MAM" = group_adm=="2",
  "NAM" = group_adm=="3",
  "CP" = group_adm=="4",
  .method = "first",
  .default = "NA"))

summary(idata$group_adm)
idata$group_adm<-droplevels(idata$group_adm)




# * 3. Sex ============================================================
### make sex human readible
idata <- idata %>% mutate(sex_adm = derivedFactor(
  "Male" = sex_adm=="1",
  "Female" = sex_adm=="2",
  .method = "first",
  .default = "NA"))

summary(idata$sex_adm)



# * 4. Check summary: for human readible ============================================================


### set labels for table column names for output
labels(idata)  <- c(site = 'Site', group_adm = "Enrolment group", sex_adm="Sex")


mycontrols  <- tableby.control(test=FALSE, total=TRUE, digits.pct = 0
                               #numeric.test="kwt", cat.test="chisq",
                               #numeric.stats=c("N", "median", "q1q3"),
                               #cat.stats=c("countpct"),
                               #stats.labels=list(N='Count', median='Median', q1q3='Q1,Q3')
)



tab1 <- tableby(group_adm ~ site + sex_adm , data=idata, control=mycontrols)
summary(tab1, text=TRUE) 


### write table
write2html(tab1, file=paste0("Summary_Tally_SamplesPerSite_",Sys.Date(),".html"))





# 4. Calculate albumin correction ratios  ============================================================

# * 1. check for albumin outliers ============================================================

#### correcting units for albumin
hist(idata$albumin_adm)
idata[which(idata$albumin_adm < 10),]

#Blantyre<-idata[idata$site =="Blantyre", ]

### write list of albumin outlier cases
#write.csv(Blantyre, file="CHAIN_admBlood_albumin_unit_check.csv")


# * 2. change albumin units ============================================================

###### replace after changing units
idata[which(idata$albumin_adm < 10),]$albumin_adm <- idata[which(idata$albumin_adm < 10),]$albumin_adm * 10
idata[which(idata$albumin_adm < 10),]$albumin_adm


# * 3. calculate albumine correciton ratio ============================================================

#### correction B2
idata$B2_corr<-signif(idata$B2/idata$albumin_adm, 4)
idata$B2_corr


#### correction B6
idata$B6_corr<-signif(idata$B6/idata$albumin_adm, 4)
idata$B6_corr


head(idata)[1:20]



# 5. Calculate other ratios and summary variables ============================================================


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





# 6. Removing variables that have too few values ============================================================

colnames(idata)


# * 1. Function to count NAs across columns
na.counts<-c()
for (i in 17:ncol(idata)){
  na.counts[i-16]<- length(which(is.na(idata[,i]) | idata[i] == "ND" )) / length(idata[,i]) *100
}
na.counts

# * 2. Organise dataframe and identify those with more than 70% missing

### [Note] This criteria is usually applied to each group of interest seperatly

na.counts<-data.frame(cbind(colnames(idata)[17:ncol(idata)], na.counts))
na.counts
to_Xclude<-na.counts[which(na.counts$na.counts > 70),]$V1
to_Xclude

# * 3. Create new dataframe excluding variables with more than 70% missing

idata[, colnames(idata) %in% to_Xclude]
idata<-idata[, !colnames(idata) %in% to_Xclude]





# 7. Replace lower than LOD  ============================================================

# * 1. Example case | Replace lower than LOD  ============================================================

#### replace Limit of detection (LOD) by half lowest detection range
idata$Carnosine<-gsub("< LOD", min(as.numeric(idata$Carnosine), na.rm=T)/10 , idata$Carnosine)
idata$Carnosine<-as.numeric(idata$Carnosine)
idata$Carnosine


colnames(idata)

# * 2. Function | Identify those that are ND (not detected) or lower than LOD =========================================
nd.counts<-c()
for (i in 17:ncol(idata)){
  nd.counts[i-16]<- length(which(idata[i] == "ND" )) / length(idata[,i]) *100
}
nd.counts


# * 3. Apply Function | Identify those that are ND (not detected) or lower than LOD =========================================
nd.counts<-data.frame(cbind(colnames(idata)[17:ncol(idata)], nd.counts))
nd.counts
to_replace<-nd.counts[which(nd.counts$nd.counts > 0),]$V1
to_replace


# * 4. Function | Replace by minimim/10 those that are ND (not detected) or lower than LOD =========================================

for(i in 1:length(to_replace)){
  idata[ , to_replace[i]]<- ifelse(idata[ , to_replace[i]] =="ND", min(as.numeric(idata[,to_replace[i]]), na.rm=T)/10 ,idata[,to_replace[i]])
}

idata[, colnames(idata) %in% to_replace]




# 8. Save new dataset =========================================

#choose.dir()
write.csv(idata, file=paste0(dirname(dirname(here())), "/5-Data/CHAIN_mNutR_idata_", Sys.Date(),".csv")) 


