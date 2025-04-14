setwd("~/Library/CloudStorage/OneDrive-Personal/Doctor/Project/Leo/2024-03-20_蛋白组癌症预测/1.Data")
library(pacman)
p_load(tidyverse,openxlsx,data.table,patchwork)


load("sig_protein_cox_brof.Rdata")

cis_pQTL_PPP <- fread("Protein/pQTL/UKB_PPP/ukb_pqtl_sigbrof_annotated.csv.gz")
clumped_dat <- fread("Protein/pQTL/UKB_PPP/clumped_4_cis.txt")

cis_pQTL_PPP <- cis_pQTL_PPP %>% semi_join(clumped_dat,
                                           by = c("CHROM" = "chr.exposure",
                                                  "GENPOS_37" = "pos.exposure",
                                                  "PROTEIN" = "PROTEIN",
                                                  "rsid" = "SNP",
                                                  "BETA" = "beta.exposure"))

cis_pQTL_PPP <- cis_pQTL_PPP %>% filter(PROTEIN %in% sig_proteins$PID)

library(TwoSampleMR)

MR_outcome <- fread("MR/MR_outcomes.csv.gz")
MR_outcome <- MR_outcome %>%
  filter(grepl("lung", trait, ignore.case = TRUE))

MR_outcome <- MR_outcome[grepl("ca", MR_outcome$trait, ignore.case = TRUE), ]


valid <- tibble()
for (i in 1:nrow(sig_proteins)) {
  
exposure <- clumped_dat %>% filter(PROTEIN == sig_proteins$PID[[i]])
print(sig_proteins$PID[[i]])

exposure_data <- read_exposure_data(
  "exposure.txt",
  clump = FALSE,
  sep = " ",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  units_col = "units",
  gene_col = "id.exposure",
  id_col = "id.exposure",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr.exposure",
  pos_col = "pos.exposure"
)

for (j in 1:nrow(MR_outcome)) {
  outcomes <- MR_outcome$id[[j]]
outcome_data <- extract_outcome_data(
    snps = exposure_data$SNP,
    outcomes = outcomes)
if (is.null(outcome_data)) {
  next  
}

harmonised_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)

mr_results <- mr(harmonised_data)
print(mr_results)
kk <- print(mr_results)
if (any(mr_results$pval<0.05)) {
  valid <- valid %>% rbind(kk)
}
next
}
}

MR_result <- valid
MR_result$`OR (95%CI)` <- paste0(round(exp(MR_result$b),2)," (",
                                 round(exp(MR_result$b-1.96*MR_result$se),2),"-",
                                 round(exp(MR_result$b+1.96*MR_result$se),2),")")

exposure_data <- read_exposure_data(
  "exposure.txt",
  clump = FALSE,
  sep = " ",
  phenotype_col = "exposure",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  eaf_col = "eaf.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  pval_col = "pval.exposure",
  units_col = "units",
  gene_col = "id.exposure",
  id_col = "id.exposure",
  min_pval = 1e-200,
  log_pval = FALSE,
  chr_col = "chr.exposure",
  pos_col = "pos.exposure"
)

outcome_data <- extract_outcome_data(
    snps = exposure_data$SNP,
    outcomes = "ieu-a-988")

harmonised_data <- harmonise_data(
  exposure_dat = exposure_data, 
  outcome_dat = outcome_data
)
mr_results <- mr(harmonised_data)
print(mr_results)

pleiotropy_test <- mr_pleiotropy_test(harmonised_data)
print(pleiotropy_test)



