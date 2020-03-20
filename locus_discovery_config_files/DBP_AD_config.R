trait1 = "AD"
trait2 = "DBP"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("diastolic","Diastolic","Diastolic Blood Pressure","DBP")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/DBP_AD_bivar_scan/DBP_dbSNP150_onlySNPs_data.txt"
trait1_BPcol = "BP"
trait2_BPcol = "chromEnd"
trait1_CHRcol = "CHR"
trait2_CHRcol = "chrom"
trait1_Pcol = "P"
trait2_Pcol = "P"
exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/DBP_AD_bivar_scan/DBP_dbSNP150_onlySNPs_data.txt",sep="\t",snp_col="name",effect_allele_col="Allele1",other_allele_col="Allele2",eaf_col="Freq1",se_col="StdErr",pval_col=trait2_Pcol,beta_col="Effect")
