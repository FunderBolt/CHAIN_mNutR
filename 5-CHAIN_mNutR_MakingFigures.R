#############################
###  CHAIN_NutR
###  Making Figures
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","ggbeeswarm","ggpubr")


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



