#############################
###  CHAIN_NutR
###  Making Figures
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","ggbeeswarm","ggpubr","gridExtra")


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
theme_set(theme_minimal())




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



######################################################
color.set<-c(c('grey80', 'springgreen2','grey60', 'grey40'))
barplot(1:4, col = color.set)


treshold_data<-read.csv("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\GLM_results_Coefficient_tables\\Coefficients_all_M1M2&M3M4_binomial_3digits_2020-11-21.csv")
colnames(treshold_data)
treshold_data<-treshold_data[,c("X","FDR_p.group_admSM_M3")]
colnames(treshold_data) <- c("metabolites", "pvalue")
head(treshold_data)

treshold_data$pvalue<-as.numeric(format(treshold_data$pvalue, scientific=FALSE))
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue>0.05, "n.s", treshold_data$pvalue)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.05, "*", treshold_data$pvalue_sig)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.001, "**", treshold_data$pvalue_sig)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.0001, "***", treshold_data$pvalue_sig)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.00001, "****", treshold_data$pvalue_sig)
treshold_data$pvalue_sig


colnames(idata)
################

plot_list <- c()

for (i in 1:length(list_met)){
  y<-i+16
  
  stat.test<-data.frame("group1"="CP", "group2"="SM","p.adj"=treshold_data$pvalue_sig[i],
                        "y.position"= max(log(idata[,y]+1), na.rm=T)+ 0.1* max(log(idata[,y]+1), na.rm=T), 
                        "y_scale_lim"= max(log(idata[,y]+1), na.rm=T) + 0.2* max(log(idata[,y]+1), na.rm=T),
                        "x_scale_lim"= ifelse(min(log(idata[,y]+1), na.rm=T) - 0.1* min(log(idata[,y]+1), na.rm=T) > 0,
                                              min(log(idata[,y]+1), na.rm=T) - 0.1* min(log(idata[,y]+1), na.rm=T), 0) ) 
  
  plot_list[[i]]<-ggplot(idata, aes(y=log(idata[,y]+1), x=group_adm)) +
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(fill=color.set[1:2])+
    geom_quasirandom(method = "quasirandom", size=1.5, col ="black", width= 0.2)+
    stat_pvalue_manual(stat.test, y.position= "y.position", label= 'p.adj', tip.length = 0.01, label.size = 8)+
    #geom_bracket()
    scale_y_continuous(name = "Log(concentration, \u03bcM)",limits=c(stat.test$x_scale_lim, stat.test$y_scale_lim)) +
    scale_x_discrete(name = "") +
    ggtitle(colnames(idata)[y])+
    theme_bw() +
    theme(plot.title = element_text(size = 16, face = "plain"),
          text = element_text(size = 14),
          axis.title = element_text(face="plain"),
          axis.text.x=element_text(size = 12),
          axis.text.y=element_text(size = 12))
  
  
    # png(sprintf("p%s.png", y),width = 1600, height = 1400)
    # plot(plot_list[[i]])
    # dev.off()
    
    pdf(sprintf("p_log_%s.pdf", y),width = 3, height = 4)
    plot(plot_list[[i]])
    dev.off()
    
    svg(sprintf("p_log_%s.svg", y),width = 3, height = 4)
    plot(plot_list[[i]])
    dev.off()
    
}




## subset to final plot list
# colnames(data)[ colnames(data) %in% list.plots]
# 
# idata<-data[, colnames(data) %in% list.vars |  colnames(data) %in% list.plots]
# head(idata)
# colnames(idata)
# 
# 
# 
# #?cowplot
# plot_grid(p12_LYSOPC16.0,p13_LYSOPC20.3,p15_SMOH22.1,p14_SMOH22.2,p16_SMOH24.1,
#           p17_PC36.0AE,p18_PC36.6AA,p19_PC36.0AA,p20_C3.1,p21_C5.1DC,
#           p9_Leucine, p10_Aspartic.acid,p11_BCAA,p1_BCAA_AAA, p2_kyn_trp,p3_urea, p4_choline, p5_betaine, p6_glutamic.ac,
#           p7_C2_C0, p8_bHydroxybutyric.ac,
#           ncol=5)
# 
# 









