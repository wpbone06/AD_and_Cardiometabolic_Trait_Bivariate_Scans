trait1 = "AD"
trait2 = "WHRadjBMI"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("WHR","Waist-hip ratio","Waist-to-hip ratio", "Waist-to-hip ratio adjusted for body mass index","Waist-hip ratio adjusted for body mass index","Waist hip ratio","waist-hip ratio","waist-to-hip ratio","WHRadjBMI")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/WHRadjBMI_AD_bivar_scan_Pulit_Sept_2019/WHRadjBMI_Pulit_data_rsID_ready.txt"
trait1_BPcol = "BP"
trait2_BPcol = "POS"
trait1_CHRcol = "CHR"
trait2_CHRcol = "CHR"
trait1_Pcol = "P"
trait2_Pcol = "P"
exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/WHRadjBMI_AD_bivar_scan_Pulit_Sept_2019/WHRadjBMI_Pulit_data_rsID_ready.txt",sep="\t",snp_col="SNP",effect_allele_col="Tested_Allele",other_allele_col="Other_Allele",eaf_col="Freq_Tested_Allele",se_col="SE",pval_col=trait2_Pcol,beta_col="BETA")
