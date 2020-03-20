trait1 = "AD"
trait2 = "CHD"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("Heart Disease","Coronary","Artery disease")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_datasets/GWAS/21_CADMeta/CAD_META"
trait1_BPcol = "BP"
trait2_BPcol = "BP"
trait1_CHRcol = "CHR"
trait2_CHRcol = "CHR"
trait1_Pcol = "P"
trait2_Pcol = "P-value"
exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_datasets/GWAS/21_CADMeta/CAD_META",sep="\t",snp_col="oldID",effect_allele_col="Allele1",other_allele_col="Allele2",eaf_col="Freq1",se_col="StdErr",pval_col=trait2_Pcol,beta_col="Effect")
