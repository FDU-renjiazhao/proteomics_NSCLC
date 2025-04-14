setwd("~/Library/CloudStorage/OneDrive-Personal/Doctor")
library(pacman)
p_load(tidyverse,glmnet,openxlsx,data.table,pROC,survival,patchwork,TwoSampleMR)

#1.目标蛋白----
load("sig_protein_cox_brof.Rdata")


Sys.setenv(OPENGWAS_JWT)

cis_pQTL_PPP <- fread("Protein/pQTL/UKB_PPP/ukb_pqtl_sig5e8_annotated.csv.gz")
clumped_dat <- fread("Protein/pQTL/UKB_PPP/ukb_cis_pqtl_clumped.txt")

cis_pQTL_PPP <- cis_pQTL_PPP %>% semi_join(clumped_dat,
                                           by = c("CHROM" = "chr.exposure",
                                                  "GENPOS_37" = "pos.exposure",
                                                  "PROTEIN" = "PROTEIN",
                                                  "rsid" = "SNP",
                                                  "BETA" = "beta.exposure"))

cis_pQTL_PPP <- cis_pQTL_PPP %>% filter(PROTEIN %in% sig_proteins$PID)
cis_pQTL_PPP <- cis_pQTL_PPP %>% select(CHR = CHROM,ID = rsid)


library(tidyverse)
library(openxlsx)
library(data.table)
library(rbgen)
snp.list <- read.xlsx("ppp_pqtl.xlsx")
output <- fread("ukb_output.csv")##ID 数据

##imputation
imputation <- function(CHR,ID){
  path_imp <- paste0("/disk/disk1/UKB/Gene/Imputation/ukb_imp_chr",CHR,"_v3.bgen")
  rsid <- as.character(ID)
  snp <- bgen.load(path_imp, rsids = c(rsid))
  A1 <- snp$variants$allele0
  A2 <- snp$variants$allele1
  A11 <- paste0(A1,A1)
  A12 <- paste0(A1,A2)
  A22 <- paste0(A2,A2)
  if(dim(snp$variants)[[1]] == 0){
    return(tibble(ID_2 = output$ID_2,!!rsid <- rep(0,487409)))
  }else{
    #analysis
    path_sample <- paste0("/disk/disk1/UKB/Gene/Imputation/sampleID/ukb63726_imp_chr",CHR,"_v3_s487296.sample")
    sample <- bigreadr::fread2(path_sample)
    sample <- sample %>% filter(sex != 'D')
    genotype <- as_tibble(snp$data[1,,])
    sample2 <- sample %>% bind_cols(genotype) %>% select(-c(ID_1,missing,sex))
    sample2 <- sample2 %>% mutate(!!rsid := ifelse(`g=0` == 1,A11,ifelse(`g=1`== 1,A12,A22))) 
    sample2 <- sample2 %>% select(ID_2,!!rsid)
    return(sample2)
  }}
snp_all <- output
for (i in 1:j) {
  CHR <- snp.list$CHR[[i]]
  ID <- snp.list$ID[[i]]
  kk <- imputation(CHR,ID)
  snp_all <- snp_all %>% left_join(kk,by = "ID_2")
}