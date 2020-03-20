trait1 = "AD"
trait2 = "SBP"
trait1GWASStr = c("Alzheimer")
trait2GWASStr = c("Systolic","Systolic Blood Pressure","SBP")
expPath="/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt"
outPath="/project/voight_GWAS/wbone/Alzheimers_centric_bivariate_scans/SBP_AD_bivar_scan/SBP_dbSNP150_onlySNPs_data.txt"
trait1_A1col = "A1"
trait1_A2col = "A2"
trait2_A1col = "Allele1"
trait2_A2col = "Allele2"
trait1_SNPcol = "SNP"
trait2_SNPcol = "name"
trait1_BPcol = "BP"
trait2_BPcol = "chromEnd"
trait1_CHRcol = "CHR"
trait2_CHRcol = "chrom"
trait1_Pcol = "P"
trait2_Pcol = "P"
trait1_Ncol = "Nsum"
trait2_Ncol = "N_effective"
trait1_MAFcol = "MAF"
trait2_MAFcol = "Freq1"

exp_dat = read_exposure_data("/project/voight_datasets/GWAS/01_alzD/AD_sumstats_Jansenetal.txt",sep="\t",snp_col="SNP",effect_allele_col="A1",other_allele_col="A2",eaf_col="MAF",se_col="SE",pval_col=trait1_Pcol,beta_col="BETA")
out_dat = read_outcome_data("/project/voight_GWAS/wbone/Alzheimers_centric_bivariate_scans/SBP_AD_bivar_scan/SBP_dbSNP150_onlySNPs_data.txt",sep="\t",snp_col="name",effect_allele_col="Allele1",other_allele_col="Allele2",eaf_col="Freq1",se_col="StdErr",pval_col=trait2_Pcol,beta_col="Effect")
#trait info not in the input file
#traitType is set either to "cc" or "quant"

trait1Type = "cc"
trait1Prop = 0.157888494
trait2Type = "quant"
trait2Prop = ""

#This is the proportion of samples that are cases in a case control GWAS, if you are using a quantitative trait this should be set to "". traitProp = cases / case + controls
