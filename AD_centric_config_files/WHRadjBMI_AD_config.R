trait1 = "AD"
trait2 = "WHRadjBMI"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("WHR","Waist-hip ratio","Waist-to-hip ratio", "Waist-to-hip ratio adjusted for body mass index","Waist-hip ratio adjusted for body mass index","Waist hip ratio","waist-hip ratio","waist-to-hip ratio","WHRadjBMI")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/WHRadjBMI_AD_bivar_scan_Pulit_Sept_2019/WHRadjBMI_Pulit_data_rsID_ready.txt"
trait1_A1col = "A1"
trait1_A2col = "A2"
trait2_A1col = "Tested_Allele"
trait2_A2col = "Other_Allele"
trait1_SNPcol = "SNP"
trait2_SNPcol = "SNP"
trait1_BPcol = "BP"
trait2_BPcol = "POS"
trait1_CHRcol = "CHR"
trait2_CHRcol = "CHR"
trait1_Pcol = "P"
trait2_Pcol = "P"
trait1_Ncol = "Nsum"
trait2_Ncol = "N"
trait1_MAFcol = "MAF"
trait2_MAFcol = "Freq_Tested_Allele"

exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/WHRadjBMI_AD_bivar_scan_Pulit_Sept_2019/WHRadjBMI_Pulit_data_rsID_ready.txt",sep="\t",snp_col="SNP",effect_allele_col="Tested_Allele",other_allele_col="Other_Allele",eaf_col="Freq_Tested_Allele",se_col="SE",pval_col=trait2_Pcol,beta_col="BETA")
#trait info not in the input file
#traitType is set either to "cc" or "quant"

trait1Type = "cc"
trait1Prop = 0.157888494
trait2Type = "quant"
trait2Prop = ""

#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
