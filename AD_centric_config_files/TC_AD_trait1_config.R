trait1 = "AD"
trait2 = "TC"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("Cholesterol","Total cho","LDL","HDL")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/TC_AD_bivar_scan_re_meta/TC_remeta_data_with_MAF.txt"
trait1_A1col = "A1"
trait1_A2col = "A2"
trait2_A1col = "A1"
trait2_A2col = "A2"
trait1_SNPcol = "SNP"
trait2_SNPcol = "SNP"
trait1_BPcol = "BP"
trait2_BPcol = "BP"
trait1_CHRcol = "CHR"
trait2_CHRcol = "CHR"
trait1_Pcol = "P"
trait2_Pcol = "P"
trait1_Ncol = "Nsum"
trait2_Ncol = "N"
trait1_MAFcol = "MAF"
trait2_MAFcol = "MAF"

exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/TC_AD_bivar_scan_re_meta/TC_remeta_data_with_MAF.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait2_Pcol,beta_col="Beta")
#trait info not in the input file
#traitType is set either to "cc" or "quant"

trait1Type = "cc"
trait1Prop = 0.157888494
trait2Type = "quant"
trait2Prop = ""

#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
