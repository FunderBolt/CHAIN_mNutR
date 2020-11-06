
#############################
###  CHAIN_NutR
###  Making Tables
#############################

# make list of needed packages
list_packages <- c("here","mgsub","tidyverse")


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




fulldata %>% 
  group_by(test_name) %>% 
  summarise(
    avg = mean(results, na.rm = T),
    medn = median(results, na.rm = T),
    minval = min(results, na.rm = T),
    maxval = max(results, na.rm = T)
  ) %>% 
  arrange(desc(medn, avg))


# dot plots

plot_lst <- c()

# number of test available
test_names <- unique(met_data_long$test_name)

for (testname in test_names){
  
  plot_lst[[testname]] <- met_data_long %>% 
    filter(test_name == testname) %>% 
    ggplot(aes(categ_enrol, results)) +
    geom_quasirandom(method = "tukeyDense", size=1.5, col ="darkblue", width= 0.15)+
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.2, color="Red")+
    stat_compare_means(label.x.npc = "centre", 
                       method.args = list(alternative = "two.sided"))+
    theme_pubclean()+
    labs(title = paste("Dotplot of", testname, ":SAM vs Community", sep = " "), 
         subtitle = "unpaired two-samples Wilcoxon test", x = "", y= "Test results")
}


pdf(file = "results_dotplots.pdf", width = 10, height = 7, onefile = TRUE)
plot_lst
dev.off()


# + Notched boxplots
# basic ggplot comparing means

box_plots <- c()


for (testname in test_names){
  
  box_plots[[testname]] <- met_data_long %>% 
    filter(test_name == testname) %>% 
    ggplot(aes(categ_enrol, results)) + geom_boxplot(width = 0.3, 
                                                     notch = T, fill = "light grey")+
    stat_compare_means(label.x.npc = "centre", 
                       method.args = list(alternative = "two.sided"))+
    theme_pubclean()+
    labs(title = paste("Notched boxplots of", testname, ":SAM vs Community", sep = " "), 
         subtitle = "unpaired two-samples Wilcoxon test", x = "", y= "Test results")
}


pdf(file = "results_boxplots.pdf", width = 10, height = 7, onefile = TRUE)
box_plots
dev.off()



# Vitamins datasets -------------------------------------------------------

# load vitamins raw datasets
water_soluble <- read_csv("water_soluble_vitamins.csv") %>% 
  select(record_id = Sample, everything())

fat_soluble <- read_csv("fat_solube_vitamins.csv") %>% 
  select(record_id = `Sample Name`, everything())

clin_chem <- read_csv("C:/Users/cmaronga/Kemri Wellcome Trust/Narshion Ngao - CHAIN Data Curation/CHAIN_Data/CHAIN_Data_v9/laboratory/Biochemistry/clean_adm_clinical_chem_data.csv",
                      na = ".")
# change to long format
# correct for wrong ID structure too

water_soluble <- water_soluble %>% 
  mutate(
    first_five = str_sub(record_id, 1, 5),
    last_three = str_sub(record_id, 6, 8),
    first_five = replace(first_five, which(first_five == 30002), 30001),
    record_id = paste(first_five, last_three, sep = "")
  ) %>% select(-c(first_five, last_three)) %>% mutate(across(record_id, as.numeric)) %>% 
  left_join(clin_chem %>% 
              select(record_id_adm, albumin_adm), by = c("record_id" = "record_id_adm")) %>% 
  mutate(
    B2_corrected = round(B2/albumin_adm, 4),
    B6_corrected = round(B6/albumin_adm, 4)
  )

water_soluble_long <- water_soluble %>% 
  pivot_longer(cols = -record_id, 
               names_to = "vitamin", 
               values_to = "test_results") %>% 
  mutate(across(.cols = vitamin, ~paste("vitamin", .))) %>%  # add prefix vitamin
  filter(vitamin != "vitamin albumin_adm")



# change fat soluble to long format
fat_soluble <- fat_soluble %>% 
  mutate(
    first_five = str_sub(record_id, 1, 5),
    last_three = str_sub(record_id, 6, 8),
    first_five = replace(first_five, which(first_five == 30002), 30001),
    record_id = paste(first_five, last_three, sep = "")
  ) %>% select(-c(first_five, last_three))


fat_soluble_long <- fat_soluble %>% 
  pivot_longer(cols = -record_id, 
               names_to = "vitamin", 
               values_to = "test_results") %>% 
  mutate(
    record_id = str_sub(record_id, 1, 8)
  ) %>% mutate(across(record_id, as.numeric))


# combine all vitamins data
vitamins_data <- bind_rows(water_soluble_long, 
                           fat_soluble_long)



# add site and ABC group
vitamins_data <- vitamins_data %>% 
  left_join(chain_data %>% mutate_at(vars(record_id), as.numeric) %>% 
              select(record_id, site, categ_enrol), by = "record_id") 

# Export dataset

write.csv(vitamins_data, "vitamins_data.csv", row.names = F)

# make dot plots
plot_lst <- c()

# number of test available
vitamins <- unique(vitamins_data$vitamin)

for (vitamin_name in vitamins){
  
  plot_lst[[vitamin_name]] <- vitamins_data %>% 
    filter(vitamin == vitamin_name) %>% 
    ggplot(aes(categ_enrol, test_results)) +
    geom_quasirandom(method = "tukeyDense", size=1.5, col ="darkblue", width= 0.15)+
    stat_summary(fun = median, fun.min = median, fun.max = median,
                 geom = "crossbar", width = 0.2, color="Red")+
    stat_compare_means(label.x.npc = "centre", 
                       method.args = list(alternative = "two.sided"))+
    theme_pubclean()+
    labs(title = paste("Dotplot of", vitamin_name, ":SAM vs Community", sep = " "), 
         subtitle = "unpaired two-samples Wilcoxon test", x = "", y= "Test results")
}


pdf(file = "vitamins_dotplots.pdf", width = 10, height = 7, onefile = TRUE)
plot_lst
dev.off()



# Make notched boxplots
box_plots <- c()


for (vitamin_name in vitamins){
  
  box_plots[[vitamin_name]] <- vitamins_data %>% 
    filter(vitamin == vitamin_name) %>% 
    ggplot(aes(categ_enrol, test_results)) + geom_boxplot(width = 0.3, 
                                                          notch = T, fill = "light grey")+
    stat_compare_means(label.x.npc = "centre", 
                       method.args = list(alternative = "two.sided"))+
    theme_pubclean()+
    labs(title = paste("Notched boxplots of", vitamin_name, ":SAM vs Community", sep = " "), 
         subtitle = "unpaired two-samples Wilcoxon test", x = "", y= "Test results")
}


pdf(file = "vitamins_boxplots.pdf", width = 10, height = 7, onefile = TRUE)
box_plots
dev.off()


# Summary table for vitamins data -----------------------------------------

SAM <- vitamins_data %>%
  filter(categ_enrol == "Acute A") %>% 
  group_by(vitamin) %>% 
  summarise(
    Mean = round(mean(test_results, na.rm = T),4),
    SD = round(sd(test_results, na.rm = T),4),
    Median = round(median(test_results, na.rm = T),4),
    Q1 = round(quantile(test_results, probs = c(0.25), na.rm = T),4),
    Q3 = round(quantile(test_results, probs = c(0.75), na.rm = T),4),
    IQR = paste(Median, "(", Q3," - ",Q1,")", sep = "")
  ) %>% 
  select(-c(Q1, Q3))



CP <- vitamins_data %>% 
  filter(categ_enrol == "Community") %>% 
  group_by(vitamin) %>% 
  summarise(
    Mean = round(mean(test_results, na.rm = T),4),
    SD = round(sd(test_results, na.rm = T),4),
    Median = round(median(test_results, na.rm = T),4),
    Q1 = round(quantile(test_results, probs = c(0.25), na.rm = T),4),
    Q3 = round(quantile(test_results, probs = c(0.75), na.rm = T),4),
    IQR = paste(Median, "(", Q3," - ",Q1,")", sep = "")
  ) %>% 
  select(-c(Q1, Q3))


summary_table <- SAM %>% 
  left_join(CP, by = "vitamin", suffix = c("_SAM","_CP")) 


write.csv(summary_table, 
          "vitamins_summary_table.csv", row.names = F)



