In this repository you will find all of the code, config files, and dependency files used to perform the Alzheimer's Disease (AD) and Cardiometabolic trait bivariate GWAS analyses to detect novel loci associated with AD and each Cardiometabolic trait

Instructions on running the bivariate scan pipeline:

* Create a dependency file directory and edit your desired bivariate GWAS script (Bivar_Script_generalized_GWAS_sep_2019.R or Bivar_Script_generalized_GWAS_trait1_focus.R) to point to the dependency file in this directory

* Create an analysis directory and add your bivariate GWAS script and config.R file corresponding the GWAS files you would like to perform a bivariate GWAS on. Then execute the bivariate GWAS script in that directory (eg Rscript ./Bivar_Script_generalized_GWAS_trait1_focus.R)


Code files:
- bivariate_scan.R => Code to perform the chi-squared statistical test given a bivariate normal distribution
- Bivar_Script_generalized_GWAS_sep_2019.R => Code to take in two hg37 GWAS summary statistic files with rsnumbers, chromosome, position, effect allele frequency, betas, and standard error, and perform a bivariate GWAS to detect novel loci associated with both traits.
- Bivar_Script_generalized_GWAS_trait1_focus.R => Code to take in two hg37 GWAS summary statistic files with rsnumbers, chromosome, position, effect allele frequency, betas, and standard error, and perform a bivariate GWAS followed by colocalization analysis of all bivariate significant loci that meet a minimum single-trait P-value for trait 1 (AD in these experiments).

Individual Run Config files:
These are the config files used for each bivariate scan we performed. NOTE: if you want to repeat our analyses you should update Bivar_Script_generalized_GWAS_sep_2019.R to include the trait string infront of the "config.R". example "BFP_AD_config.R" instead of "config.R"

- BFP_AD_config.R => config file for body fat percentage and AD bivariate GWAS
- BMI_AD_config.R => config file for body mass index  and AD bivariate GWAS
- CHD_AD_config.R => config file for congestive heart disease and AD bivariate GWAS
- DBP_AD_config.R => config file for diastolic blood pressure and AD bivariate GWAS
- HDL_AD_config.R => config file for high-density lipoproteins and AD bivariate GWAS
- LDL_AD_config.R => config file for low-density lipoproteins and AD bivariate GWAS
- SBP_AD_config.R => config file for systolic blood pressure and AD bivariate GWAS
- T2D_AD_config.R => config file for type 2 diabetes and AD bivariate GWAS
- TC_AD_config.R => config file for total cholesterol and AD bivariate GWAS
- TG_AD_config.R => config file for triglycerides and AD bivariate GWAS
- WHRadjBMI_AD_config.R => config file for waist-hip-ratio adjusted for Body mass index and AD bivariate GWAS

Dependency files:

- EUR.final => holds the 503 population ids associated with 1KG Phase III vcf NOTE: self-generated
- EUR.final.plink => if the plink file than contains the data for only the individuals from EUR.final
- GWAShg37.bed => BED file of the January 2019 GWAS Catalog flat file converted into BED file format
- genes_biomart_sorted.bed => bed file of gene positions from biomart

