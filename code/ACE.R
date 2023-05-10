# ACE inhibition and Alzheimer's disease tutorial  
# CHARGE 2023
# David Ryan, Neil Davies, Venexia Walker

# Setup ------------------------------------------------------------------------

# Uncomment the following lines to install libraries that you do not already have
# install.packages('data.table')
# install.packages("devtools")
# devtools::install_github("MRCIEU/TwoSampleMR")

# Import libraries 
library(data.table)
library(TwoSampleMR)

# Create empty results table 
results <- NULL

# Load data --------------------------------------------------------------------

# Download the exposure data - Evangelou et al, 2018 SBP GWAS 
sbp_gwas <- fread("data/sbp_subset.csv")

# Download the outcome data - de Rojas et al. 2020 Alzhemier's GWAS
alz_gwas <- fread("data/alz_subset.csv")

# Format data for TwoSampleMR package ------------------------------------------

sbp_gwas$phenotype <- "Systolic blood pressure"
sbp_gwas <- format_data(
  sbp_gwas,
  type="exposure",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "StdErr",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  eaf_col = "Freq1",
  pval_col = "P-value",
  chr_col = "chromosome",
  pos_col = "position",
  phenotype_col = "phenotype")

alz_gwas$phenotype <- "Alzheimer's disease"
alz_gwas <- format_data(
  alz_gwas,
  type="outcome",
  snp_col = "RS",
  beta_col = "Beta",
  se_col = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col = "P",
  chr_col = "CHR",
  pos_col = "BP",
  phenotype_col = "phenotype")

# Restrict the exposure to SNPs that are genome-wide significant ---------------

sbp_significant <- subset(sbp_gwas, pval.exposure < 5e-8)

# Specify target details for ACE (hg19, chr17:61,554,422-61,599,205) -----------

# Specify the chromosome the target gene is on
target_chr <- 17 

# Specify the start and end positions of the target gene
target_pos_start <- 61554422
target_pos_end <- 61599205 

# Specify how far outside of the gene you are willing to look (e.g., 20kb)
tolerance <- 20000 

# Analysis 1: ACE > Alzheimer's disease ----------------------------------------

# Filter to SNPs within the target gene
ace_unclumped <- subset(sbp_significant, 
               chr.exposure==target_chr & 
                 pos.exposure > (target_pos_start-tolerance) & 
                 pos.exposure < (target_pos_end+tolerance))

# Update exposure name
ace_unclumped$exposure <- "ACE"

# Clump instrument
ace <- clump_data(
  ace_unclumped,
  clump_kb = 10000,
  clump_r2 = 0.01,
  pop = "EUR")
data.table::fwrite(ace,"data/ace_clumped.csv", row.names = FALSE)

# If you have to wifi, please read in this already clumped data instead
# ace <- data.table::fread("data/ace_clumped.csv")

# How many SNPs are in our instrument? 
length(unique(ace$SNP))

# Harmonise exposure and outcome data
ace_alz_harmonised <- harmonise_data(ace, alz_gwas)

# Perform MR analysis 
ace_alz_results <- mr(ace_alz_harmonised, method_list=c("mr_wald_ratio"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
ace_alz_results$or <- exp(-1*ace_alz_results$b)

# Calculate confidence interval 
ace_alz_results$or_lci <- exp(-1*ace_alz_results$b-qnorm(0.975)*ace_alz_results$se)
ace_alz_results$or_uci <-exp(-1*ace_alz_results$b+qnorm(0.975)*ace_alz_results$se)

# Display results
ace_alz_results[,c("exposure","outcome","method","nsnp","or","or_lci","or_uci","pval")]

# Add results to table
results <- rbind(results, ace_alz_results)

# Analysis 2: Systolic blood pressure > Alzheimer's disease --------------------

# The dataframe 'sbp_significant' is already restricted to significant SNPs
# This means no further processing is required for systolic blood pressure

# Clump instrument
sbp <- clump_data(
  sbp_significant,
  clump_kb = 10000,
  clump_r2 = 0.01,
  pop = "EUR")
data.table::fwrite(sbp,"data/sbp_clumped.csv", row.names = FALSE)

# If you have to wifi, please read in this already clumped data instead
# sbp <- data.table::fread("data/sbp_clumped.csv")

# How many SNPs are in our instrument? 
length(unique(sbp$SNP))

# Harmonise exposure and outcome data
sbp_alz_harmonised <- harmonise_data(sbp, alz_gwas)

# Perform MR analysis 
sbp_alz_results <- mr(sbp_alz_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_ivw_mre"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
sbp_alz_results$or <- exp(-1*sbp_alz_results$b)

# Calculate confidence interval 
sbp_alz_results$or_lci <- exp(-1*sbp_alz_results$b - qnorm(0.975)*sbp_alz_results$se)
sbp_alz_results$or_uci <- exp(-1*sbp_alz_results$b + qnorm(0.975)*sbp_alz_results$se)

# Display results
sbp_alz_results[,c("exposure","outcome","method","nsnp","or","or_lci","or_uci","pval")]

# Add results to table
results <- rbind(results, sbp_alz_results)

# Analysis 3: Systolic blood pressure, without ACE > Alzheimer's disease -------

# Remove the effect at the ACE locus from the systolic blood pressure instrument
sbp_noace_unclumped <- subset(sbp_significant, 
                              !(chr.exposure==target_chr & 
                                  pos.exposure > (target_pos_start-tolerance) & 
                                  pos.exposure < (target_pos_end+tolerance)))

# Update exposure name
sbp_noace_unclumped$exposure <- "Systolic blood pressure (no ACE)"

# Clump instrument
sbp_noace <- clump_data(
  sbp_noace_unclumped,
  clump_r2 = 0.01,
  clump_kb = 10000,
  pop = "EUR")
data.table::fwrite(sbp_noace,"data/sbp_noace_clumped.csv", row.names = FALSE)

# If you have to wifi, please read in this already clumped data instead
# sbp_noace <- data.table::fread("data/sbp_noace_clumped.csv")

# How many SNPs are in our instrument? 
length(unique(sbp_noace$SNP))

# Harmonise exposure and outcome data
sbp_noace_alz_harmonised <- harmonise_data(sbp_noace, alz_gwas)

# Perform MR analysis 
sbp_noace_alz_results <- mr(sbp_noace_alz_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_ivw_mre"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
sbp_noace_alz_results$or <- exp(-1*sbp_noace_alz_results$b)

# Calculate confidence interval 
sbp_noace_alz_results$or_lci <- exp(-1*sbp_noace_alz_results$b - qnorm(0.975)*sbp_noace_alz_results$se)
sbp_noace_alz_results$or_uci <- exp(-1*sbp_noace_alz_results$b + qnorm(0.975)*sbp_noace_alz_results$se)

# Display results
sbp_noace_alz_results[,c("exposure","outcome","method","nsnp","or","or_lci","or_uci","pval")]

# Add results to table
results <- rbind(results, sbp_noace_alz_results)

# Display all results ----------------------------------------------------------

results <- results[,c("exposure","outcome","method","nsnp","or","or_lci","or_uci","pval")]
results 
fwrite(results,"output/ace_results.csv", row.names = FALSE)