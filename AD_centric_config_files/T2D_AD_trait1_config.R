trait1 = "AD"
trait2 = "T2D"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("T2D","type II diabetes","Type 2 diabetes","type ii diabetes")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/T2D_AD_bivar_scan/T2D_meta_data_with_rsIDs.txt"
trait1_A1col = "A1"
trait1_A2col = "A2"
trait2_A1col = "EA"
trait2_A2col = "NEA"
trait1_SNPcol = "SNP"
trait2_SNPcol = "SNP"
trait1_BPcol = "BP"
trait2_BPcol = "Pos"
trait1_CHRcol = "CHR"
trait2_CHRcol = "Chr"
trait1_Pcol = "P"
trait2_Pcol = "Pvalue"
trait1_Ncol = "Nsum"
trait2_Ncol = "Neff"
trait1_MAFcol = "MAF"
trait2_MAFcol = "EAF"

exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/T2D_AD_bivar_scan/T2D_meta_data_with_rsIDs.txt",sep="\t",snp_col="SNP",effect_allele_col="EA",other_allele_col="NEA",eaf_col="EAF",se_col="SE",pval_col=trait2_Pcol,beta_col="Beta")
#trait info not in the input file
#traitType is set either to "cc" or "quant"

trait1Type = "cc"
trait1Prop = 0.157888494
trait2Type = "cc"
trait2Prop = 0.1800179

#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
