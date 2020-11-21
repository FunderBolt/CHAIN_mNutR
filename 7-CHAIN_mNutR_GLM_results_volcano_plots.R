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
hist(FC_data$FC_median)

## add in metabolites column for merge
FC_data$metabolites<-row.names(FC_data)
head(FC_data)
str(FC_data)





# add threshold 
### pulling the p-values from the glm results
#choose.files()
treshold_data<-read.csv("D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\GLM_results_Coefficient_tables\\Coefficients_all_M1M2&M3M4_binomial_3digits_2020-11-21.csv")
colnames(treshold_data)
treshold_data<-treshold_data[,c("metabolites","FDR_p.group_admSM_M3")]
head(treshold_data)
str(treshold_data)


head(FC_data)
plot_data<-left_join(FC_data, treshold_data, by="metabolites" )
colnames(plot_data)<-c("median_case", "median_control", "group_diff_median", "FC_median", 
                       "log2FC", "metabolites", "pvalue")


########################
write.csv(plot_data, file=paste0("D:/Dropbox/Bandsma.Lab/1.Projects/1.CHAIN_NETWORK/2019_CHAIN_Micronutrient/7-Analysis-Results/medianIQR_&_GLMs/List_Metabolites_FC&FDR_sig_", Sys.Date(), ".csv"))




dput(plot_data$metabolites)
vitamins <-c("Retinol", "OH25_VitD3", "aTocopherol", "bg_Tocopherol", 
                    "B1", "B2", "B3_amide", "B5", "B6", "B7", "C","B2_corr", "B6_corr")

metals<-c("Barium", "Sodium", "Magnesium", 
            "Phosphorous", "Calcium", "Titanium", "Vanadium", "Manganese", 
            "Iron", "Cobalt", "Copper", "Zinc", "Gallium", "Arsenic", "Selenium", 
            "Rubidium", "Strontium", "Antimony", "Tellurium")

just_metabolites<-c("Creatinine", 
            "Glycine", "Alanine", "Serine", "Proline", "Valine", "Threonine", 
            "Taurine", "Putrescine", "trans.Hydroxyproline", "Leucine", "Isoleucine", 
            "Asparagine", "Aspartic_ac", "Glutamine", "Glutamic_ac", "Methionine", 
            "Histidine", "a_Aminoadipic_ac", "Phenylalanine", "Methionine.sulfoxide", 
            "Arginine", "Acetyl.ornithine", "Citrulline", "Serotonin", "Tyrosine", 
            "Asymmetric_dimethylarginine", "total_dimethylarginine", "Tryptophan", 
            "Kynurenine", "Carnosine", "Ornithine", "Lysine", "Spermidine", 
            "Spermine", "Sarcosine", "Creatine", "Betaine", "Choline", "Trimethylamine_N.oxide", 
            "Methylhistidine", "Homocysteine", "Hypoxanthine", "Ethanolamine", 
            "Cystathionine", "N.acetylpetrescine", "Homocitrulline", "Urea", 
            "Hydroxy.lysine", "a_aminobutyric_ac", "g_aminobutyric_ac", "Lactic_ac", 
            "b_Hydroxybutyric_ac", "a_Ketoglutaric_ac", "Citric_ac", "Butyric_ac", 
            "Propionic_ac", "HPHPA", "p.Hydroxyhippuric_ac", "Succinic_ac", 
            "Fumaric_ac", "Pyruvic_ac", "Isobutyric_ac", "Hippuric_ac", "Methylmalonic_ac", 
            "Homovanillic_ac", "Indole_acetic_ac", "Uric_ac", "OH5_Indoleacetic_ac", 
            "p.Hydroxyphenylacetic_ac", "Valeric_ac", "Uridine", "Xanthine", 
            "Guanidineacetic_ac", "Acetyl.lysine", "Glucose", "AAA", 
            "BCAA", "BCAA_AAA", "Fischer_ratio", "Essential_AA", "Glucogenic_AA", 
            "Keto_AA", "Kynurenine_Trp", "Orn_Arg", "Putrescine_Orn", "Serotonin_Trp", 
            "Spermidine_Putrescine", "Spermine_Spermidine", "Tyr_Phe", "Total_DMA_Arg", 
            "Cit_Arg", "Total_AAs", "urea_cycle")
            
            
just_lipids<-c("LYSOC14.0", 
            "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", "LYSOC18.1", 
            "LYSOC18.0", "LYSOC20.4", "LYSOC20.3", "LYSOC24.0", "LYSOC26.1", 
            "LYSOC26.0", "LYSOC28.1", "LYSOC28.0", "Sphg_14.1SMOH", "Sphg_16.1SM", 
            "Sphg_16.0SM", "Sphg_16.1SMOH", "Sphg_18.1SM", "PC32.2_AA", "Sphg_18.0SM", 
            "Sphg_20.2SM", "PC36.0_AE", "PC36.6_AA", "PC36.0_AA", "Sphg_22.2SMOH", 
            "Sphg_22.1SMOH", "PC38.6_AA", "PC38.0_AA", "PC40.6_AE", "Sphg_24.1SMOH", 
            "PC40.6_AA", "PC40.2_AA", "PC40.1_AA","PC40.2") 
            
            
just_carnitines<-c( "C0", "C2", "C3.1", "C3", 
            "C4.1", "C4", "C3OH", "C5.1", "C5", "C4OH", "C6.1", "C6", "C5OH", 
            "C5.1DC", "C5DC", "C8", "C5MDC", "C9", "C7DC", "C10.2", "C10.1", 
            "C10", "C12.1", "C12", "C14.2", "C14.1", "C14", "C12DC", "C14.2OH", 
            "C14.1OH", "C16.2", "C16.1", "C16", "C16.2OH", "C16.1OH", "C16OH", 
            "C18.2", "C18.1", "C18", "C18.1OH","C2_C0")


### select label to plot
all_results<-plot_data
plot_data<-all_results

plot_data<-plot_data[plot_data$metabolites %in% vitamins, ]
plot_data<-plot_data[plot_data$metabolites %in% metals, ]
plot_data<-plot_data[plot_data$metabolites %in% just_metabolites, ]
plot_data<-plot_data[plot_data$metabolites %in% just_lipids, ]
plot_data<-plot_data[plot_data$metabolites %in% just_carnitines, ]


dput(plot_data$metabolites[plot_data$pvalue<0.05 & (plot_data$log2FC > 0.5849625 | plot_data$log2FC < -0.5849625)])

to_label<-c("B1", "Alanine", "Threonine", "trans.Hydroxyproline", "Glutamine", 
            "Arginine", "Citrulline", "Tyrosine", "Tryptophan", "Ornithine", 
            "Ethanolamine", "Hydroxy.lysine", "b_Hydroxybutyric_ac", "Butyric_ac", 
            "Isobutyric_ac", "Homovanillic_ac", "Valeric_ac", "LYSOC14.0", 
            "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", "LYSOC18.0", 
            "LYSOC20.3", "PC36.6_AA", "PC38.6_AA", "PC40.6_AE", "PC40.6_AA", 
            "C0", "C2", "C4OH", "C9", "C10.2", "C16.1", "Barium", "Titanium", 
            "Manganese", "Arsenic", "Tellurium", "B2_corr", "Kynurenine_Trp", 
            "Putrescine_Orn", "Serotonin_Trp", "Tyr_Phe", "Total_DMA_Arg", 
            "C2_C0", "urea_cycle")

to_label<-c("B1", "B2_corr")

to_label<-c("Barium", "Titanium", "Manganese", "Arsenic", "Tellurium")

to_label<-c("Alanine", "Threonine", "trans.Hydroxyproline", "Glutamine", 
            "Arginine", "Citrulline", "Tyrosine", "Tryptophan", "Ornithine", 
            "Ethanolamine", "Hydroxy.lysine", "b_Hydroxybutyric_ac", "Butyric_ac", 
            "Isobutyric_ac", "Homovanillic_ac", "Valeric_ac", "Kynurenine_Trp", 
            "Putrescine_Orn", "Serotonin_Trp", "Tyr_Phe", "Total_DMA_Arg", 
             "urea_cycle")

to_label<-c("LYSOC14.0", "LYSOC16.1", "LYSOC16.0", "LYSOC17.0", "LYSOC18.2", 
            "LYSOC18.0", "LYSOC20.3", "PC36.6_AA", "PC38.6_AA", "PC40.6_AE", 
            "PC40.6_AA")

to_label<-c("C0", "C2", "C4OH", "C9", "C10.2", "C16.1","C2_C0")


#### checking threshold to use
logratio2foldchange(0.5, base=2)
foldchange2logratio(1.5, base=2)


log10(range(plot_data$pvalue)[1])
range(plot_data$log2FC)

### ploting
head(plot_data)
volcano_plot<-EnhancedVolcano(plot_data,
                lab = plot_data$metabolites,
                x = 'log2FC',
                y = 'pvalue',
                xlim = c(-5, 5),
                ylim = c(0, 6),
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
                
              
volcano_plot

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\results_Volcano_plot_ALL_2020-11-21.svg", width = 12, height = 8, bg = "white")
volcano_plot
dev.off()

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\results_Volcano_plot_vitamines_2020-11-13.svg", width = 9, height = 8, bg = "white")
volcano_plot
dev.off()

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\results_Volcano_plot_metals_2020-11-21.svg", width = 9, height = 8, bg = "white")
volcano_plot
dev.off()

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\results_Volcano_plot_Metabs_2020-11-13.svg", width = 9, height = 8, bg = "white")
volcano_plot
dev.off()

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\results_Volcano_plot_Lipids_2020-11-13.svg", width = 9, height = 8, bg = "white")
volcano_plot
dev.off()

svg(file = "D:\\Dropbox\\Bandsma.Lab\\1.Projects\\1.CHAIN_NETWORK\\2019_CHAIN_Micronutrient\\7-Analysis-Results\\medianIQR_&_GLMs\\results_Volcano_plot_Carnitines_2020-11-13.svg", width = 9, height = 8, bg = "white")
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












