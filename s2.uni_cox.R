setwd("~/Library/CloudStorage/OneDrive-Personal")
library(pacman)
p_load(tidyverse,survival,openxlsx,ggrepel,data.table)

load("ukb_lung_proteomic_data_with_covar_clean.Rdata")


protein_position <- c(which(names(dat) == "p845"):which(names(dat) == "p648"))
vars <- names(dat)[protein_position]


Unicox_crude <- function(x,dt,cc){
  FML <- as.formula(paste0("Surv(py,",cc,") ~ ",x,"+age+age2+sex+smoking+race"))
  fit <- coxph(formula = FML,data = dt)
  Gsum <- summary(fit)
  output <- tibble("Characteristic" = names(fit$coefficients),
                   "Hazard Ratio" = round(Gsum$coefficients[,2],2),
                   "Confidence Interval" = paste0(round(Gsum$conf.int[,3],2),"-",
                                                  round(Gsum$conf.int[,4],2)),
                   "lower" = round(Gsum$conf.int[,3],2),
                   "upper" = round(Gsum$conf.int[,4],2),
                   "P" = Gsum$coefficients[,5] )
  return(output)
}

Unicox_adjusted <- function(x,dt,cc){
  FML <- as.formula(paste0("Surv(py,",cc,") ~ ",x,"+age+age2+sex+smoking+race+BMI+
                           fev_ratio+education_q+nd_2010+no2_2010+pm10_2010+
                           pm25_2010+pm_2010+COPD+vegetable+family_history"))
  fit <- coxph(formula = FML,data = dt)
  Gsum <- summary(fit)
  output <- tibble("Characteristic" = names(fit$coefficients),
                   "Hazard Ratio" = round(Gsum$coefficients[,2],2),
                   "Confidence Interval" = paste0(round(Gsum$conf.int[,3],2),"-",
                                                  round(Gsum$conf.int[,4],2)),
                   "lower" = round(Gsum$conf.int[,3],2),
                   "upper" = round(Gsum$conf.int[,4],2),
                   "P" = Gsum$coefficients[,5] )
  return(output)
}

##cox----
cox_model1 <- tibble()
for (i in 1:length(vars)) {
  tmp <- Unicox_crude(vars[[i]],train,"cc_nsclc")
  cox_model1 <- rbind(cox_model1,tmp)
}
cox_model1 <- cox_model1 %>% filter(!Characteristic %in% c("age","age2","sex",
                                                           "smoking","race"))  
  
cox_model2 <- tibble()
for (i in 1:length(vars)) {
  tmp <- Unicox_adjusted(vars[[i]],train,"cc_nsclc")
  cox_model2 <- rbind(cox_model2,tmp)
}
cox_model2 <- cox_model2 %>% 
  filter(!Characteristic %in% c("age","age2","sex","smoking","BMI","fev_ratio",
                                "education_q","nd_2010","no2_2010","pm10_2010",
                                "pm25_2010","pm_2010","COPD","vegetable",
                                "family_history","race"))    









