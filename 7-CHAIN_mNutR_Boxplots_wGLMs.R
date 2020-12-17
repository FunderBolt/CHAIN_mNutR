# ====  CHAIN_NutR | Boxplots with GLM results ===================================


# 0. Packages ============================================================
# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","ggbeeswarm","ggpubr","gridExtra")
install.packages("tidymodels")

# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)

##set overall theme
#theme_set(theme_minimal())



# 1. Load data ============================================================

### load full data
idata <- read.csv(paste0(dirname(dirname(here())), "/5-Data","/CHAIN_mNutR_fulldata_X_2020-11-21.csv"), row.names = 1) 


### load threshodl data
treshold_data<-read.csv("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\GLM_results_Coefficient_tables\\Coefficients_all_M1M2&M3M4_binomial_3digits_2020-11-21.csv")
colnames(treshold_data)
treshold_data<-treshold_data[,c("X","FDR_p.group_admSM_M3")]
colnames(treshold_data) <- c("metabolites", "pvalue")
head(treshold_data)



# 2. Set Factor ============================================================

### list factors to convert
str(idata)
factor_list<-c("record_id", "sex_adm","parttype_adm","site","africa","group_adm","age_group_adm",
               "urban","oedema_adm","adm_dead","dead")

idata[, factor_list] <- lapply(idata[, factor_list], as.factor)
str(idata)




# 3. Set colors ============================================================

color.set<-c(c('grey80', 'springgreen2','grey60', 'grey40'))
barplot(1:4, col = color.set)



# 4. Recode p-values ============================================================

treshold_data$pvalue<-as.numeric(format(treshold_data$pvalue, scientific=FALSE))
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue > 0.05, "n.s", treshold_data$pvalue)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.05, "*", treshold_data$pvalue_sig)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.001, "**", treshold_data$pvalue_sig)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.0001, "***", treshold_data$pvalue_sig)
treshold_data$pvalue_sig<-ifelse(treshold_data$pvalue < 0.00001, "****", treshold_data$pvalue_sig)
treshold_data$pvalue_sig



# 5. Function | Boxplots loop ============================================================

colnames(idata)
### This creates and saves each individual Boxplot into folder


plot_list <- c()


for (i in 1:length(list_met)){
  y<-i+16
  
  ### this section is to indicate the p-value to plot and the location of the "bar" 
  stat.test<-data.frame("group1"="CP", "group2"="SM","p.adj"=treshold_data$pvalue_sig[i],
                        ### this indicates the height of the "bar", i add ~10% in height from the highest value 
                        "y.position"= max(log(idata[,y]+1), na.rm=T)+ 0.1* max(log(idata[,y]+1), na.rm=T), 
                        ### to create enough room for the bar in set the y limits to + ~20% in height from the highest value 
                        "y_scale_lim"= max(log(idata[,y]+1), na.rm=T) + 0.2* max(log(idata[,y]+1), na.rm=T),
                        "x_scale_lim"= ifelse(min(log(idata[,y]+1), na.rm=T) - 0.1* min(log(idata[,y]+1), na.rm=T) > 0,
                                              min(log(idata[,y]+1), na.rm=T) - 0.1* min(log(idata[,y]+1), na.rm=T), 0) ) 
  
  plot_list[[i]]<-ggplot(idata, aes(y=log(idata[,y]+1), x=group_adm)) +
    ### to add whiskers
    stat_boxplot(geom = "errorbar", width = 0.2) + 
    geom_boxplot(fill=color.set[1:2])+
    ### to add the individual data points over the boxplots
    geom_quasirandom(method = "quasirandom", size=1.5, col ="black", width= 0.2)+
    ### adding the p-value bars
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
    
    pdf(file = paste0(dirname(dirname(here())), "/7-Analysis-Results/medianIQR_&_GLMs/Boxplots_with_GLM_results/", sprintf("p_log_%s.pdf", y)),width = 3, height = 4)
    plot(plot_list[[i]])
    dev.off()
    
    svg(file = paste0(dirname(dirname(here())), "/7-Analysis-Results/medianIQR_&_GLMs/Boxplots_with_GLM_results/", sprintf("p_log_%s.svg", y)),width = 3, height = 4)
    plot(plot_list[[i]])
    dev.off()
    
    
}



### I save both as .svg and pdf but I usually find the PDF more helpful

### [Note] as mentionned this code needs to be updated to save the plots within a list
###        so better able to control grid plotting using ggarrange


#ggarrange(plot_list[[1]], plot_list[[1]])



