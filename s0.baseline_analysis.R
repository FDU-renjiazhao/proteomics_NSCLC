setwd("~/Library/CloudStorage/OneDrive-Personal")
library(pacman)
p_load(tidyverse,survival,glmnet,openxlsx,ggrepel,data.table)

load("ukb_lung_proteomic_data.Rdata")

#3.
dat <- dat %>%
  mutate(Country = case_when(
    location %in% c(11005,11004) ~ "Scotland",
    location %in% c(11003,11022,11023) ~ "Wales",
    TRUE ~ "England" 
  ))

table(dat$Country)

dat <- dat %>% mutate(class = case_when(
  Country %in% c("Wales","Scotland") ~"Test",
  Country %in% c("England") ~"Train",
))

#1.
dat$age2 <- dat$age*dat$age
## smoking
smoking <- fread("ukb_smoking.csv.gz") %>% select(ID,ever_smoked)
# covar
covar <- fread("ukb_lung_covar.csv.gz")
covar$smoking <- ifelse(covar$smoking == "",NA,covar$smoking)
covar$smoking <- covar$smoking %>% as.numeric()
covar <- covar %>% mutate(education_q = coalesce(education_eng_quintile, 
                                                 education_scot_quintile,
                                                 education_wales_quintile))
covar <- covar %>% select(!c("education_eng_quintile", 
                             "education_scot_quintile",
                             "education_wales_quintile",
                             "education_eng", 
                             "education_scot",
                             "education_wales"))
covar <- covar %>% select(!alcohol)
#COPD
copd <- fread("ukb_disease_long.csv.gz") %>% 
  filter(dis_before_entry == "Yes") %>% select(ID,ICD_code)
copd <- copd[grep("J44",copd$ICD_code),]
#lifestyle
lifestyle <- read_csv("ukb_lifestyle_score.csv") %>% select(ID,vegetable = vegetable2)
#family_history
family_history <- fread("ukb_family_history_long.csv.gz")
lung_history <- family_history %>% filter(Name == "Lung cancer") %>% select(ID,Name)
lung_history <- lung_history[!duplicated(lung_history$ID),]
lung_history <- lung_history %>% rename(family_history = Name)
lung_history$family_history <- 1
#combine
covar <- covar %>% full_join(copd,by = "ID") %>% full_join(lifestyle,by = "ID") %>% full_join(lung_history,by = "ID")
covar <- covar %>% rename(COPD = ICD_code)
covar$COPD <- ifelse(is.na(covar$COPD),0,1)
covar$family_history <- ifelse(is.na(covar$family_history),0,1)
covar <- covar[!duplicated(covar$ID),]
covar <- covar %>% left_join(smoking)

dat <- dat %>% left_join(covar)

dat <- dat %>% select(ID:cc_nsclc_above_10,Country:ever_smoked,p845:p648)
dat$class <- factor(dat$class,
                    levels = c("Train","Test"))


dat <- dat %>% filter(!is.na(cc_nsclc))
library(tableone)

dat <- dat %>% mutate(race = case_when(
  race %in% c(1,1001,1002,1003) ~ "White",
  race %in% c(2,2001,2002,2003,2004) ~ "Mixed",
  race %in% c(3,3001,3002,3003,3004,5) ~ "Asian",
  race %in% c(4,4001,4002,4003) ~ "Black",
  race %in% c(-3,-1,6) ~ "Others/Unknown",
  is.na(race) ~"Others/Unknown"
))

dat$race <- factor(dat$race,
                   levels = c("Asian","Black","Mixed","White","Others/Unknown"))

dat$smoking <- factor(dat$smoking,
                      levels = c(0,2,1),
                      labels = c("No","Occasionally","Yes"))

dat$sex <- factor(dat$sex,
                  levels = c(0,1),
                  labels = c("Female","Male"))

dat$COPD <- factor(dat$COPD,
                  levels = c(0,1),
                  labels = c("No","Yes"))

dat$family_history <- factor(dat$family_history,
                  levels = c(0,1),
                  labels = c("No","Yes"))

library(tableone)

catvars <- c("race","sex","smoking","COPD","family_history","education_q",
             "vegetable","ever_smoked")

train <- dat %>% filter(class == "Train")
test <- dat %>% filter(class == "Test")


tableOne <- CreateTableOne(vars = vars, 
                           data = dat,
                           strata = "class",
                           factorVars = catvars,
                           addOverall = T,
                           includeNA = T)

tableOne 
print(tableOne,showAllLevels = T)



tableOne <- CreateTableOne(vars = vars, 
                           data = train,
                           strata = "cc_nsclc",
                           factorVars = catvars,
                           addOverall = T,
                           includeNA = T)
tableOne 
print(tableOne,showAllLevels = T)

