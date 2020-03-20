trait1 = "AD"
trait2 = "CHD"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("Heart Disease","Coronary","Artery disease")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/CHD_AD_centric_bivar_scan_Dec_2019/CAD_META_with_N"
trait1_A1col = "A1"
trait1_A2col = "A2"
trait2_A1col = "Allele1"
trait2_A2col = "Allele2"
trait1_SNPcol = "SNP"
trait2_SNPcol = "oldID"
trait1_BPcol = "BP"
trait2_BPcol = "BP"
trait1_CHRcol = "CHR"
trait2_CHRcol = "CHR"
trait1_Pcol = "P"
trait2_Pcol = "P-value"
trait1_Ncol = "Nsum"
trait2_Ncol = "N"
trait1_MAFcol = "MAF"
trait2_MAFcol = "Freq1"

exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/CHD_AD_centric_bivar_scan_Dec_2019/CAD_META_with_N",eaf_col="Freq1",effect_allele_col = "Allele1",other_allele_col="Allele2",pval_col=trait1_Pcol,snp_col="oldID",se_col="StdErr",beta_col="Effect",sep="\t")

#trait info not in the input file
#traitType is set either to "cc" or "quant"

trait1Type = "cc"
trait1Prop = 0.157888494
trait2Type = "cc"
trait2Prop = 0.2242678

#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
