setwd("~/Library/CloudStorage/OneDrive-Personal")
library(pacman)
p_load(tidyverse,survival,openxlsx,data.table,pROC,nricens,survival)

load("ukb_lung_proteomic_data_with_covar_clean.Rdata")


fit1 <- coxph(Surv(py, cc_nsclc) ~ age + sex + smoking, 
              x = TRUE,
              data = train)

fit2 <- coxph(Surv(py, cc_nsclc) ~ age + fev_ratio +smoking+
                nd_2010+pm25_2010+pm_2010+education_q+COPD+family_history+vegetable, 
              x = TRUE,
              data = train)

P.eallY.fit1 <- get.risk.coxph(fit1, time_point)
P.eallY.fit2 <- get.risk.coxph(fit2, time_point)

deciles <- quantile(P.eallY.fit1, 
                    probs = seq(0.1, 1, by = 0.1))


nri_result <- nricens(mdl.std = fit1, 
                      mdl.new = fit2, 
                      t0 = time_point, 
                      p.std=P.eallY.fit1,
                      p.new=P.eallY.fit2,
                      cut = deciles,
                      niter = 10)


fit3 <- coxph(Surv(py, cc_nsclc) ~ p746+p512+p2891, 
              x = TRUE,
              data = train)

P.eallY.fit3 <- get.risk.coxph(fit3, time_point)


fit4 <- coxph(Surv(py, cc_nsclc) ~ age + fev_ratio + smoking+
                COPD+p746+p512+p2891, 
              x = TRUE,
              data = train)
P.eallY.fit4 <- get.risk.coxph(fit4, time_point)

deciles <- quantile(P.eallY.fit3, 
                    probs = seq(0.1, 1, by = 0.1))

nri_result <- nricens(mdl.std = fit3, 
                      mdl.new = fit4, 
                      t0 = time_point, 
                      p.std=P.eallY.fit3,
                      p.new=P.eallY.fit4,
                      cut = deciles,
                      niter = 10)

##IDI
library(survIDINRI)

y <- train %>% select(py,cc_nsclc) %>% as.matrix()

model1 <- train %>% 
  select(age,sex,smoking) %>% 
  as.matrix()

model2 <- train %>% 
  select(age , fev_ratio,smoking,vegetable,
         nd_2010,pm25_2010,pm_2010,education_q,COPD,family_history) %>% 
  as.matrix()

model3 <- train %>% 
  select("p746","p512","p2891") %>% 
  as.matrix()

model4 <- train %>% 
  select(age , fev_ratio , smoking,
         COPD,p746,p512,p2891) %>% 
  as.matrix()

x <- IDI.INF(indata = y,
             covs0 = model1,
             covs1 = model2,
             t0=5,
             npert = 200)

IDI.INF.OUT(x)

x <- IDI.INF(indata = y,
             covs0 = model3,
             covs1 = model4,
             t0=5,
             npert = 200)

IDI.INF.OUT(x)




