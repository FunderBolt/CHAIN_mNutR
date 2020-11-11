#############################
###  CHAIN_NutR
###  Volcano Plot from GLM results
#############################

## install packages
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')


# make list of needed packages
list_packages <- c("here","mgsub","tidyverse","EnhancedVolcano","gtools","ggthemes")



# list any missing packages
new.packages <- list_packages[!(list_packages %in% installed.packages()[,"Package"])]
### if packages missing --> install
if(length(new.packages) > 0) {install.packages(new.packages,dependencies = TRUE)}
## load all packages
lapply(list_packages, require, character.only = TRUE)



##set overall theme
theme_set(theme_minimal())



# load full data
idata <- read.csv(paste0(dirname(dirname(dirname(here("6-Data")))), "/5-Data","/CHAIN_mNutR_fulldata_X_2020-11-08.csv"), row.names = 1) 



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



## Prepare the dataset
#calculate median value for each metabolite per case
case<- subset(idata, idata$group_adm=="SM")
head(case)
#case_m <- apply(case[,17:ncol(case)], 2, mean, na.rm=T)
case_m <- apply(case[,17:ncol(case)], 2, median, na.rm=T)


#calculate mean value for each metabolite per control
control<- subset(idata, idata$group_adm=="CP")
head(control)
#control_m <- apply(control[,17:ncol(control)], 2, mean, na.rm=T)
control_m <- apply(control[,17:ncol(control)], 2, median, na.rm=T)


# calculate log2 fold change for each metabolite 

# FC_data<-as.data.frame(cbind(case_m, control_m, case_m-control_m, foldchange(case_m,control_m)))
# colnames(FC_data) <- c("mean_case", "mean_control","group_diff",  "FC")
# format(FC_data, scientific=FALSE)

FC_data<-as.data.frame(cbind(case_m, control_m, case_m-control_m, foldchange(case_m,control_m)))
colnames(FC_data) <- c("median_case", "median_control","group_diff_median","FC_median")
format(FC_data, scientific=FALSE)

##add column for log2 base 
FC_data$log2FC<-log2(FC_data$median_case) - log2(FC_data$median_control)
hist(FC_data_compiled$FC_median)

## add in metabolites column for merge
FC_data$metabolites<-row.names(FC_data)
head(FC_data)
str(FC_data)





# add threshold 
### pulling the p-values from the glm results
#choose.files()
treshold_data<-read.csv("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\BoxCoxTransform_coefficients_all_M1M2&M3_3digits_2020-11-08.csv")
colnames(treshold_data)
treshold_data<-treshold_data[,c("X","FDR_p.group_admSM_M3")]
colnames(treshold_data) <- c("metabolites", "pvalue")
head(treshold_data)
str(treshold_data)
treshold_data$metabolites<-as.character(treshold_data$metabolites)

plot_data<-left_join(FC_data, treshold_data, by="metabolites" )



########################
write.csv(plot_data, file=paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/GLMs/List_Metabolites_FC&FDR_sig_", Sys.Date(), ".csv"))


### select label to plot
dput(plot_data$metabolites[plot_data$pvalue<0.05 & (plot_data$log2FC > 0.5849625 | plot_data$log2FC < -0.5849625)])

to_label<-c("B1", "Alanine", "Threonine", "trans.Hydroxyproline", "Glutamine", 
            "Arginine", "Citrulline", "Tyrosine", "Tryptophan", "Ornithine", 
            "Ethanolamine", "Hydroxy.lysine", "b_Hydroxybutyric_ac", "Butyric_ac", 
            "Isobutyric_ac", "Homovanillic_ac", "Valeric_ac", "LYSOC14.0", 
            "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", "LYSOC18.0", 
            "LYSOC20.3", "PC36.6_AA", "PC38.6_AA", "PC40.6_AE", "PC40.6_AA", 
            "C0", "C2", "C4OH", "C9", "C10.2", "C16.1")


logratio2foldchange(0.5, base=2)
foldchange2logratio(1.5, base=2)

head(plot_data)
volcano_plot<-EnhancedVolcano(plot_data,
                lab = plot_data$metabolites,
                x = 'log2FC',
                y = 'pvalue',
                xlim = c(-3, 3),
                ylim = c(0, 12.5),
                title = 'SM versus Community',
                subtitle = '',
                subtitleLabSize = 12,
                xlab = bquote(~Log[2]~ 'FC'),
                ylab = bquote(~-Log[10]~FDR-italic(p)),
                axisLabSize = 16,
                caption = paste0('Total = ', nrow(plot_data), ' vitamins and metabolites'),
                captionLabSize = 12,
                pCutoff = 0.05,
                FCcutoff = 0.5849625,
                pointSize = 3.0,
                labSize = 3.5,
                labCol = 'grey30',
                col=c('grey80', 'grey60', 'grey40', 'springgreen4'),
                colAlpha = 0.70,
              legendLabels=c(' n.s.', bquote(~Log[2]~FC), bquote(~FDR-italic(p)), bquote(~Log[2]~'FC &'~FDR-italic(p))),
                legendPosition = 'right',
             # legendPosition = 'none',
                #legendLabSize = 16,
               # legendIconSize = 5.0,
              drawConnectors = TRUE,
              widthConnectors = 0.2,
              colConnectors = 'grey30',
              selectLab = to_label)
                
              

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\GLMs\\results_Volcano_plot.svg", width = 9, height = 8, bg = "white")
volcano_plot
dev.off()


# ## Volcano plot with ggplot
# ggplot(final) +
#   geom_point(aes(x=counts, y=pvalue, colour=threshold.y)) +
#   ggtitle("Mov10 overexpression") +
#   xlab("log2 fold change") + 
#   ylab("-log10 adjusted p-value") +
#   #scale_y_continuous(limits = c(0,50)) +
#   theme(legend.position = "none",
#         plot.title = element_text(size = rel(1.5), hjust = 0.5),
#         axis.title = element_text(size = rel(1.25)))  












