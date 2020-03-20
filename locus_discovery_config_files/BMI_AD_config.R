trait1 = "AD"
trait2 = "BMI"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("BMI","body mass index","Body Mass Index")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_datasets/GWAS/04_giant/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"
trait1_BPcol = "BP"
trait2_BPcol = "POS"
trait1_CHRcol = "CHR"
trait2_CHRcol = "CHR"
trait1_Pcol = "P"
trait2_Pcol = "P"
exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_datasets/GWAS/04_giant/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt",sep="\t",snp_col="SNP",effect_allele_col="Tested_Allele",other_allele_col="Other_Allele",eaf_col="Freq_Tested_Allele_in_HRS",se_col="SE",pval_col=trait2_Pcol,beta_col="BETA")
