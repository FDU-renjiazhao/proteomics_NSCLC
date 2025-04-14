setwd("~/Library/CloudStorage/OneDrive-Personal")
library(pacman)
p_load(tidyverse,survival,glmnet,openxlsx,data.table,pROC)

load("ukb_lung_proteomic_data_with_covar_clean.Rdata")

load("sig_protein_cox_brof.Rdata")

sig_proteins <- sig_protein_cox_brof %>% filter(Characteristic %in% 
                                                   c("p746","p512","p2891"))

dat <- dat %>% select(ID:family_history,all_of(sig_proteins$Characteristic)) %>% 
  filter(!is.na(cc_nsclc))

train <- train %>% select(ID:family_history,all_of(sig_proteins$Characteristic)) %>% 
  filter(!is.na(cc_nsclc))

test <- test %>% select(ID:family_history,all_of(sig_proteins$Characteristic)) %>% 
  filter(!is.na(cc_nsclc))


train$smoking <- factor(train$smoking)
test$smoking <- factor(test$smoking)

#3.Model----

train5 <- train %>% select(ID:location,py,cc_nsclc_5,Country:names(train)[ncol(train)]) %>% na.omit()
test5 <- test %>% select(ID:location,py,cc_nsclc_5,Country:names(train)[ncol(train)]) %>% na.omit()

fit1 <- coxph(Surv(py,cc_nsclc_5)~age+sex+smoking,
              x = TRUE,
              data = train5)
summary(fit1)

test5$prob1 <- predict(fit1,
                       newdata = test5,
                       type = "risk")
roc1 <- roc(test5$cc_nsclc_5,test5$prob1)
print(roc1)
auc1 <- ci.auc(roc1)

###3.1.2----
x <- data.matrix(train5[,c(3,4,11:22)])
y <- data.matrix(Surv(train5$py,train5$cc_nsclc_5))

cvfit_lasso <- cv.glmnet(x,y,
                         family = "cox",
                         nfold = 10,
                         alpha = 1)

lasso_coefs <- coef(cvfit_lasso, s = "lambda.min")
selected_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

formula_lifestyle <- as.formula(paste("Surv(py,cc_nsclc_5) ~", 
                                      paste(selected_vars, 
                                            collapse = " + ")))

fit2 <- coxph(formula_lifestyle,
              x = TRUE,
              data = train5)
summary(fit2)

test5$prob2 <- predict(fit2,
                       newdata = test5,
                       type = "risk")
roc2 <- roc(test5$cc_nsclc_5,test5$prob2)
print(roc2)
auc2 <- ci.auc(roc2)

###3.1.3----
x <- data.matrix(train5[,c(23:ncol(train5))])
y <- data.matrix(Surv(train5$py,train5$cc_nsclc_5))

cvfit_lasso <- cv.glmnet(x,y,
                         family = "cox",
                         nfold = 10,
                         alpha = 1)
lasso_coefs <- coef(cvfit_lasso, s = "lambda.min")
selected_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

formula_protein <- as.formula(paste("Surv(py,cc_nsclc_5) ~", 
                                    paste(selected_vars, 
                                          collapse = " + ")))

fit3 <- coxph(formula_protein,
              data = train5)
summary(fit3)

test5$prob3 <- predict(fit3,
                       newdata = test5,
                       type = "risk")
roc3 <- roc(test5$cc_nsclc_5,test5$prob3)
print(roc3)
auc3 <- ci.auc(roc3)

###3.1.4 ----
x <- data.matrix(train5[,c(3,4,11:ncol(train5))])
y <- data.matrix(Surv(train5$py,train5$cc_nsclc_5))

cvfit_lasso <- cv.glmnet(x,y,
                         family = "cox",
                         nfold = 10,
                         alpha = 1)
lasso_coefs <- coef(cvfit_lasso, s = "lambda.1se")
selected_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]

formula_full <- as.formula(paste("Surv(py,cc_nsclc_5) ~", 
                                 paste(selected_vars, 
                                       collapse = " + ")))

fit4 <- coxph(formula_full,
              data = train5)
summary(fit4)

test5$prob4 <- predict(fit4,
                       newdata = test5,
                       type = "risk")
roc4 <- roc(test5$cc_nsclc_5,test5$prob4)
print(roc4)
auc4 <- ci.auc(roc4)

#4 ROC----
library(plotROC)
best.1 <- pROC::coords(roc1, "best", ret = "all", transpose = FALSE)
best.2 <- pROC::coords(roc2, "best", ret = "all", transpose = FALSE)
best.3 <- pROC::coords(roc3, "best", ret = "all", transpose = FALSE)
best.4 <- pROC::coords(roc4, "best", ret = "all", transpose = FALSE)
plot_roc_5 <- ggplot(test5, aes(d = cc_nsclc_5, m = prob1)) +
  geom_roc(data = test5, aes(d = cc_nsclc_5, m = prob1), #age+sex
           n.cuts = 0,
           col = "#5FA67F") +
  geom_roc(data = test5, aes(d = cc_nsclc_5, m = prob2), #age+sex+smoking
           n.cuts = 0,
           col = "#403990") +
  geom_roc(data = test5, aes(d = cc_nsclc_5, m = prob3), #Full
           n.cuts = 0,
           col = "#F46F43") +
  geom_roc(data = test5, aes(d = cc_nsclc_5, m = prob4), #bio markers
           n.cuts = 0,
           col = "black") +
  geom_point(y = best.1$sensitivity,
             x = 1-best.1$specificity,
             color = "#5FA67F",
             size = 2)+
  geom_point(y = best.2$sensitivity,
             x = 1-best.2$specificity,
             color = "#403990",
             size = 2)+
  geom_point(y = best.3$sensitivity,
             x = 1-best.3$specificity,
             color = "#F46F43",
             size = 2)+
  geom_point(y = best.4$sensitivity,
             x = 1-best.4$specificity,
             color = "black", 
             size = 2)+
  coord_equal() +
  theme_bw()+
  style_roc(xlab = "1 - Specificity",
            ylab = "Sensitivity")+
  theme(axis.text.x = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill=NA, size=1.5),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(size = 15, face = "bold", hjust = 0.5))+
  geom_abline(intercept=seq(-100, 100, 25), 
              slope=1, 
              colour="darkgrey")+
  labs(tag = "(A)",size = 40,fontface = "bold")+
  ggtitle("5-Year Prediction")


