##Figure.S2, pie plot----
load("ukb_lung_proteomic_data_with_covar_clean.Rdata")

library(ggplot2)
library(ggpubr)

lung_data <- tibble(
  Subtype = c("Adenocarcinoma", 
              "Squamous Cell Carcinoma", 
              "Large Cell Carcinoma", 
              "Other Specified Carcinoma",
              "Unspecified Malignant Neoplasms",
              "Small Cell Lung Cancer"),
  Cases = c(sum(train$cc_lac == 1,na.rm = T),
            sum(train$cc_lscc == 1,na.rm = T),
            sum(train$cc_lcc == 1,na.rm = T),
            sum(train$cc_others == 1,na.rm = T),
            sum(train$cc_unknown == 1,na.rm = T),
            sum(train$cc_sclc == 1,na.rm = T)))

lung_data <- lung_data %>%
  mutate(Total = sum(Cases),                 
         Percentage = Cases / Total * 100)    


lung_data <- lung_data %>%
  mutate(Label = sprintf("%s\nCases: %d (%.1f%%)", Subtype, Cases, Percentage))


pie <- ggpie(lung_data, "Cases", label = "Label", 
             palette = "jco", 
             fill = "Subtype", 
             color = "white", 
             lab.pos = 'out',
             lab.font = c(4, 'black',"bold"), 
             title = "Distribution of Lung Cancer Subtypes in Discovery Cohort")


ggsave(pie, height = 7, width = 7, file = "Figure.S2.pdf")
