# IL-6 tutorial  
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

# Download the exposure data - UK Biobank CRP GWAS 
crp_gwas <- fread("data/crp_subset.csv")

# Format exposure data for TwoSampleMR package ---------------------------------

crp_gwas$phenotype <- "C-reactive protein"
crp_gwas <- format_data(
  crp_gwas,
  type="exposure",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  eaf_col = "EAF",
  pval_col = "p",
  chr_col = "chr",
  pos_col = "pos",
  phenotype_col = "phenotype")

# Specify target details for IL-6R (hg19, chr1:154,377,669-154,441,926) --------

# Specify the chromosome the target gene is on
target_chr <- 1

# Specify the start and end positions of the target gene
target_pos_start <- 154377669
target_pos_end <- 154441926 

# Specify how far outside of the gene you are willing to look (e.g., 20kb)
tolerance <- 300000 

# Make IL-6R instrument --------------------------------------------------------

# Restrict to genome-wide significant SNPs
crp_significant <- subset(crp_gwas, pval.exposure < 5e-8)

# Select variants within 300kB of the IL-6R gene 
il6r_unclumped <- subset(crp_significant, 
                         chr.exposure==target_chr & 
                           pos.exposure > (target_pos_start-tolerance) & 
                           pos.exposure < (target_pos_end+tolerance))

# Clump instrument
il6r <- clump_data(
  il6r_unclumped,
  clump_kb = 10000,
  clump_r2 = 0.01,
  pop = "EUR")

# How many SNPs are in our instrument? 
length(unique(il6r$SNP))

# Analysis 1: IL-6R > CAD ------------------------------------------------------

# Load outcome data
cad_gwas <-  read_outcome_data(
  filename = "data/cad_subset.csv",
  sep = ",",
  snp_col = "markername",
  beta_col = "beta",
  se_col = "se_dgc",
  effect_allele_col = "effect_allele",
  eaf_col = 'effect_allele_freq',
  other_allele_col = 'noneffect_allele',
  pval_col = "p_dgc")
cad_gwas$outcome <- "Coronary artery disease"

# Harmonise exposure and outcome data
il6r_cad_harmonised <- harmonise_data(il6r, cad_gwas)

# Perform MR analysis
il6r_cad_results <- mr(il6r_cad_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_ivw_mre"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
il6r_cad_results$est <- exp(il6r_cad_results$b)

# Calculate confidence interval 
il6r_cad_results$est_lci <- exp(il6r_cad_results$b - qnorm(0.975)*il6r_cad_results$se)
il6r_cad_results$est_uci <- exp(il6r_cad_results$b + qnorm(0.975)*il6r_cad_results$se)

# Display results
il6r_cad_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_cad_results$est_type <- "OR"
results <- rbind(results, il6r_cad_results)

# Analysis 2: IL-6R > Acute ischaemic stroke -----------------------------------

# Load outcome data
ais_gwas <- read_outcome_data(
  filename = "data/ais_subset.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  pval_col = "P")

# Harmonise exposure and outcome data
il6r_ais_harmonised <- harmonise_data(il6r, ais_gwas)

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
il6r_ais_results$est <- exp(il6r_ais_results$b)

# Calculate confidence interval 
il6r_ais_results$est_lci <- exp(il6r_ais_results$b - qnorm(0.975)*il6r_ais_results$se)
il6r_ais_results$est_uci <- exp(il6r_ais_results$b + qnorm(0.975)*il6r_ais_results$se)

# Display results
il6r_ais_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_ais_results$est_type <- "OR"
results <- rbind(results, il6r_ais_results)

# Analysis 3: IL-6R > Cardioembolic stroke -------------------------------------

# Load outcome data
ces_gwas <- read_outcome_data(
  filename = "data/ces_subset.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  pval_col = "P")
ces_gwas$outcome <- "Cardioembolic stroke"

# Harmonise exposure and outcome data
il6r_ces_harmonised <- harmonise_data(il6r, ces_gwas)

# Perform MR analysis
il6r_ces_results <- mr(il6r_ces_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_ivw_mre"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
il6r_ces_results$est <- exp(il6r_ces_results$b)

# Calculate confidence interval 
il6r_ces_results$est_lci <- exp(il6r_ces_results$b - qnorm(0.975)*il6r_ces_results$se)
il6r_ces_results$est_uci <- exp(il6r_ces_results$b + qnorm(0.975)*il6r_ces_results$se)

# Display results
il6r_ces_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_ces_results$est_type <- "OR"
results <- rbind(results, il6r_ces_results)

# Analysis 4: IL-6R > Any stroke -----------------------------------

# Load outcome data
as_gwas <- read_outcome_data(
  filename = "data/any_stroke_subset.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  pval_col = "P")
as_gwas$outcome <- "Any stroke"

# Harmonise exposure and outcome data
il6r_as_harmonised <- harmonise_data(il6r, as_gwas)

# Perform MR analysis
il6r_as_results <- mr(il6r_as_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_ivw_mre"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
il6r_as_results$est <- exp(il6r_as_results$b)

# Calculate confidence interval 
il6r_as_results$est_lci <- exp(il6r_as_results$b - qnorm(0.975)*il6r_as_results$se)
il6r_as_results$est_uci <- exp(il6r_as_results$b + qnorm(0.975)*il6r_as_results$se)

# Display results
il6r_as_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_as_results$est_type <- "OR"
results <- rbind(results, il6r_as_results)

# Analysis 5: IL-6R > Chronic kidney disease -----------------------------------

# Load outcome data
ckd_gwas <-read_outcome_data(
  filename = "data/ckd_subset.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  eaf_col = "EAF",
  pval_col = "P")
ckd_gwas$outcome <- "Chronic kidney disease"

# Harmonise exposure and outcome data
il6r_ckd_harmonised <- harmonise_data(il6r, ckd_gwas)

# Perform MR analysis
il6r_ckd_results <- mr(il6r_ckd_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_ivw_mre"))

# Multiply by -1 to represent a decrease in SBP and exponentiate the MR estimate to get the odds ratio 
il6r_ckd_results$est <- exp(il6r_ckd_results$b)

# Calculate confidence interval 
il6r_ckd_results$est_lci <- exp(il6r_ckd_results$b - qnorm(0.975)*il6r_ckd_results$se)
il6r_ckd_results$est_uci <- exp(il6r_ckd_results$b + qnorm(0.975)*il6r_ckd_results$se)

# Display results
il6r_ckd_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_ckd_results$est_type <- "OR"
results <- rbind(results, il6r_ckd_results)

# Analysis 6: IL-6R > eGFR -----------------------------------------------------

# Load outcome data
egfr_gwas <-read_outcome_data(
  filename = "data/egfr_subset.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  eaf_col = "EAF",
  pval_col = "P")
egfr_gwas$outcome <- "eGFR"

# Harmonise exposure and outcome data
il6r_egfr_harmonised <- harmonise_data(il6r, egfr_gwas)

# Perform MR analysis
il6r_egfr_results <- mr(il6r_egfr_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_ivw_mre"))

# Continuous outcome so beta is our final estimate
il6r_egfr_results$est <- il6r_egfr_results$b

# Calculate confidence interval 
il6r_egfr_results$est_lci <- il6r_egfr_results$b - qnorm(0.975)*il6r_egfr_results$se
il6r_egfr_results$est_uci <- il6r_egfr_results$b + qnorm(0.975)*il6r_egfr_results$se

# Display results
il6r_egfr_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_egfr_results$est_type <- "Beta"
results <- rbind(results, il6r_egfr_results)

# Analysis 7: IL-6R > Blood urea nitrogen -----------------------------------

# Load outcome data
bun_gwas <-read_outcome_data(
  filename = "data/bun_subset.csv",
  sep = ",",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "EA",
  other_allele_col = 'NEA',
  eaf_col = "EAF",
  pval_col = "P")
bun_gwas$outcome <- "Blood urea nitrogen"

# Harmonise exposure and outcome data
il6r_bun_harmonised <- harmonise_data(il6r, bun_gwas)

# Perform MR analysis
il6r_bun_results <- mr(il6r_bun_harmonised, method_list=c("mr_simple_median", "mr_weighted_median", "mr_egger_regression", "mr_ivw_mre"))

# Continuous outcome so beta is our final estimate
il6r_bun_results$est <- il6r_bun_results$b

# Calculate confidence interval 
il6r_bun_results$est_lci <- il6r_bun_results$b - qnorm(0.975)*il6r_bun_results$se
il6r_bun_results$est_uci <- il6r_bun_results$b + qnorm(0.975)*il6r_bun_results$se

# Display results
il6r_bun_results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]

# Add results to table
il6r_bun_results$est_type <- "Beta"
results <- rbind(results, il6r_bun_results)

# Display all results ----------------------------------------------------------

results <- results[,c("exposure","outcome","method","nsnp","est","est_lci","est_uci","pval")]
results 
fwrite(results,"output/il-6_results.csv", row.names = FALSE)
