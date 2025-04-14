setwd("~/Library/CloudStorage/OneDrive-Personal/Doctor")
library(pacman)
p_load(tidyverse,openxlsx,data.table)

#1.----
##1.1 UKB-PPP----
PPP <- fread("Protein/pQTL/UKB_PPP/ukb_protein_gwas_sig1e5_original.csv.gz")
PPP$PROTEIN <- substring(PPP$PROTEIN, 8)
PPP <- PPP %>% rename(GENPOS_38 = GENPOS)
PPP$GENPOS_37 <- sub(".*:(.*?):.*", "\\1", PPP$ID) %>% as.integer()
PPP <- PPP %>% select(CHROM,GENPOS_37,GENPOS_38,everything())


snp_annotation <- fread("Protein/pQTL/UKB_PPP/snp_annotation_ref.csv.gz")

PPP <- PPP %>% left_join(snp_annotation,
                         by = c("ID",
                                 "GENPOS_37" = "POS19",
                                 "GENPOS_38" = "POS38",
                                 "ALLELE0" = "REF",
                                 "ALLELE1" = "ALT"))

protein_map <- fread("Protein/pQTL/UKB_PPP/olink_protein_map_3k_v1.tsv")
protein_map <- protein_map %>% select(OlinkID,UniProt,UniProt2,Assay,HGNC.symbol,
                                      ensembl_id,
                                      chr,gene_start,gene_end,Strand)
protein_map$chr <- ifelse(protein_map$chr == "X",23,protein_map$chr) %>% as.integer()

PPP <- PPP %>% left_join(protein_map,by = c("PROTEIN" = "HGNC.symbol",
                                             "OlinkID",
                                             "UniProt"))

log_threshold <- -log10(5e-8/2911)
PPP <- fread("Protein/pQTL/UKB_PPP/ukb_protein_gwas_sig1e5_annotated.csv.gz")

PPP <- PPP %>% mutate(F_stat = (BETA^2) / (SE^2)) %>% 
  filter(LOG10P >= log_threshold,
         str_starts(rsid, "rs"))#只保留以 rs 的 SNP

PPP$P <- -PPP$LOG10P
PPP$P <- 10^PPP$P

cis_pQTL_PPP <- PPP %>% filter(CHROM == chr) %>% 
                        filter(GENPOS_38 >= (gene_start - 1e+6) & 
                               GENPOS_38 <= (gene_end + 1e+6))
table(cis_pQTL_PPP$PROTEIN)






















