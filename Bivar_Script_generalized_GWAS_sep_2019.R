#Load the following on an interactive session before using this script:
#module load JAGS
#module load OpenBUGS
#module load plink/1.90Beta4.5
#module load R


library(TwoSampleMR)
library(miscF)
library(data.table)
source("/project/voight_selscan/ksiewert/Migraine_Project/Scripts/bivariate_scan.R")
source("config.R")
ts = paste0(trait1,trait2,sep = "")


#Takes list of GWAS traits that overlap current SNP, returns if they overlap any in strList1 or strList2
noMatch <-function(currTraits,strList1,strList2){
	for (tStr in strList1){
		if (any(grep(tStr,currTraits,ignore.case=TRUE)))
		{
			return(FALSE)
		}
	}
	for (tStr in strList2){
		if (any(grep(tStr, currTraits,ignore.case=TRUE)))
		{
			return(FALSE)
		}}
	return(TRUE)
}

#r is row of data showing bivariate scan data for that SNP, this function finds all GWAS hits from catalog within 500kb, without considering LD
findNearGWASCatalogHits<-function(r){
	Chrom = trimws(r["CHR"])
	Coord=r["BP"]
	res = c("","","","")
	col = paste0("chr",Chrom)
	for (i in 1:dim(GWASlist[[col]])[1]){
		G = GWASlist[[col]][i,]
		if (abs(as.numeric(Coord)-as.numeric(G["coord"]))<500000){
			res[1]= paste(res[1],G["trait"])
			res[2] = paste(res[2],G["PMID"])
			res[3] = paste(res[3],G["rs"])
			res[4] = paste(res[4],abs(as.numeric(Coord)-as.numeric(G["coord"])))
			}}
	return (res)
}



#Takes in line with SNP proxies in, then looks at if any of those proxies are in the GWAS catalog
f = function(SNPProxies){
	leadSNP = SNPProxies["SNP"]
	SNPProxyList = unlist(strsplit(as.character(SNPProxies["TAGS"]),'|',fixed=TRUE))
	matches = which(GWAS[,"rs"] %in% SNPProxyList)
	if(length(matches)==0){
	return(c(leadSNP,"NONE","NONE","NONE"))
	}
	return(c(leadSNP,paste(GWAS[matches,"rs"],collapse = ","),paste(GWAS[matches,"trait"],collapse = ","),paste(GWAS[matches,"PMID"],collapse=",")))
}


#Find most significant SNP in peak using single trait GWAS p-values
findHighestSNP = function(proxySnps,GWASFile1,pvalcolName){
	proxyLines = GWASFile1[which(GWASFile1$SNP %in% unlist(proxySnps)),]
	maxLine = proxyLines[which.min(unlist(proxyLines[pvalcolName])),]
	return(maxLine)
}

getProxies = function(SNPLine){
	SNPProxyList = append(SNPLine["SNP"],strsplit(as.character(SNPLine["TAGS"]),'|',fixed=TRUE))
	return(SNPProxyList)
}



#Warning: The harmonise_data function harmonizes the outcome to the exposure. If the alleles in your exposure file are polarized so that the trait increasing allele is always the effect allele, instead of random, this will cause the bivariate normal to be skewed. For this reason, if one of your traits has the trait increaseing allele as the effect allele for all SNPs, use it as the outcome. If both are polarized, then you'll have to randomize them.
har = harmonise_data(exposure_dat=exp_dat, outcome_dat= out_dat) 
scoreOut = write.table(har, quote = FALSE,row.names=FALSE,file=paste0("AllSNPs_",ts,"_bivarpvalues_nomin.txt"),sep="\t")
#har$z_out = har$beta.outcome/har$se.outcome #Modify this to get GC corrected z scores
#har$z_exp = har$beta.exposure/har$se.exposure
har$pval.outcome = pmax(har$pval.outcome,1e-16) #Need to take pmax so don't get underflow with estimating p
har$pval.exposure = pmax(har$pval.exposure,1e-16)

har$z_out =  sqrt(qchisq(1 - har$pval.outcome,1))*sign(har$beta.outcome)
har$z_exp =  sqrt(qchisq(1 - har$pval.exposure,1))*sign(har$beta.exposure)

write.table(har[c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","beta.exposure","beta.outcome","pval.outcome","pval.exposure","se.outcome","se.exposure","z_out","z_exp")],paste0("AllRs",ts,".txt"),quote = FALSE,row.names=FALSE)


#BIVARIATE SCAN PARAM ESTIMATE & Scan
#Get LD-independent SNPs for estimation of bivariate params
for(ch in 1:22){
	system(paste0("plink --noweb --vcf /project/voight_datasets_01/1kg/phaseIII_2013/ALL.chr",toString(ch),".phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --keep /project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink  --extract AllRs",ts,".txt --out onlyGWASNPS_",ts,"_ch",toString(ch)," --indep-pairwise 1000 5 0.2 --chr ",toString(ch)))
}
system(paste0("cat onlyGWASNPS_",ts,"_ch*.prune.in > onlyGWASNPS_",ts,"_allch.prune.in"))

indRS = read.table(paste0("onlyGWASNPS_",ts,"_allch.prune.in"), header=FALSE,stringsAsFactors = FALSE,fill=TRUE,blank.lines.skip=TRUE, colClasses = "character")


#estimate params of bivariate normal distribution
res = har[har$SNP %in% indRS$V1,]
save(res,file="res_obj_onlyLDIndepInfo")
zs=res[,c("z_out","z_exp")]
p= mvn.ub(zs)
print(p)
inv = solve(p$hatSigma)
print(inv)


#Do bivariate scan of all SNPs
bi_neglogpvalue = apply(har[,c("z_out","z_exp")],c(1),df2Test,p$hatMu,inv)
scores = cbind(data.frame(har,bi_neglogpvalue))
scores$bipvalue = 10^(-1*scores$bi_neglogpvalue)
colnames(scores)[colnames(scores) == 'bipvalue'] <- 'P'
scoreOut =write.table(scores, quote = FALSE,row.names=FALSE,file=paste0("AllSNPs_",ts,"_bivarresults.txt"),sep="\t")


#At this point, we have a bivariate p-value for every SNP in the genome that's been assayed for both phenotypes

#Do LD-clumping so we get clumps/loci of significant SNPs.
system(paste0("plink --noweb --bfile /project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed  --keep /project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink --clump-p1 .000001 --clump-r2 0.2  --clump-kb 1000 --clump AllSNPs_",ts,"_bivarresults.txt --out allChrMergedClumped"))
cData = read.table("allChrMergedClumped.clumped",header=TRUE,stringsAsFactors = FALSE)



#Find nearest gene, Turn clumped file into bed format
system("awk ' {OFS=\"\t\"} (NR>1 && NF>0) {print \"chr\"$1,$4,$4+1,$3}' allChrMergedClumped.clumped | sort -k1,1 -k2,2n>  allChrMergedClumped.bed")  #Puts top SNP into bedfile format
system(paste0("closestBed -a allChrMergedClumped.bed -b /project/voight_selscan/ksiewert/BetaGenomeScan/genes_biomart_sorted.bed -d > ",ts,"_NearestGene.bed"))

nearGene = read.table(paste0(ts,"_NearestGene.bed"),header=FALSE,stringsAsFactors = FALSE)
colnames(nearGene) = c("chr","coord1","coord2","SNP","chr2","c1","c2","gene","dist")
nearGene = aggregate(nearGene[,c("gene","dist")],list(nearGene$SNP),toString)
colnames(nearGene)[1] = "SNP"
d = merge(cData,har[,c("SNP","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","eaf.exposure","eaf.outcome","beta.exposure","beta.outcome","se.exposure","se.outcome","pval.exposure","pval.outcome","z_exp","z_out")],by="SNP")
d = merge(d,nearGene[c("SNP","gene","dist")],by="SNP")
d$sameDir =  sign(d$beta.exposure)==sign(d$beta.outcome)


#Filter for loci that are at least nominally significant
o = order(d$P,decreasing=FALSE)
d = d[o,]
d = d[which(d$P<1e-7),] #contains all clustered SNP with p<1e-7 whether or not novel
write.table(d, quote = FALSE,row.names=FALSE,file=paste0("AllSigSNPs_",ts,"_r2Indepdent.txt"),sep="\t") #NOTE: May not be 1MB from each other, but will be r2<.2


#First, return hits from the GWAS catalog within 500kb, no matter what trait they correspond to
GWAS = read.table("/project/voight_datasets/GWAS/15_GWASCatalog/Jan_2019_version/GWAShg37.bed",header=FALSE,sep="\t",colClasses = "character",fill=TRUE,quote = "",stringsAsFactors = FALSE)
colnames(GWAS)=c("chr","coord","coord+1","PMID","rs","trait","P")
GWAS = GWAS[which(as.numeric(GWAS$P)<=5e-8),] #only want significant SNPs from GWAS catalog
GWASlist = split(GWAS,GWAS$chr)
GWASoverlap =t(apply(d,1,findNearGWASCatalogHits)) #finds gwas hits within 500kb
colnames(GWASoverlap) = c("nearby_trait","nearby_PMID","nearby_rs","nearby_coord")


#Next, want to filter out any SNPs that have a relevant GWAS result within 500KB bases, because this means they are not novel
GWASannot = cbind(d,GWASoverlap)
g = unlist(lapply(GWASoverlap[,1],noMatch,trait1GWASStr,trait2GWASStr)) #Checks if any of the GWAS catalog traits are for the traits of interest (using strings from config.R file)
novel = GWASannot[which(g),]
novelSNPs = novel[which(novel$pval.outcome>5e-8 & novel$pval.exposure>5e-8 & novel$pval.outcome<.005 & novel$pval.exposure<.005),] #pull out ones with indiv p-values that are nominally, but not genome-wide significant, because these are the ones that are novel and that there's compelling evidence for it being associated with both traits

write.table(novelSNPs, quote = FALSE,row.names=FALSE,file=paste0("Novelhits_",ts,"_bivarresults.txt"),sep="\t")
write.table(novelSNPs$SNP,row.names=FALSE,col.names=FALSE,quote=FALSE,file="novelSNPrsNums.txt",sep="\t")


#See what other traits are within r^2=.2 of novel GWAS hits
system(paste0("plink --bfile /project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed --keep /project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink --show-tags novelSNPrsNums.txt --out NovelSNPTags  --tag-kb 1000 --list-all --tag-r2 0.2"))
SNPs = read.table("NovelSNPTags.tags.list", header=TRUE,stringsAsFactors = FALSE)
overlappingTraits = t(apply(SNPs,1,f))
colnames(overlappingTraits)=c("SNP","Linked_SNPs","Linked_Trait","Linked_Link")



#Will need to make sure none of the proxies are relevant GWAS hits-similar process to above but with r^2 instead of distance
notKnowns = lapply(overlappingTraits[,"Linked_Trait"],noMatch,trait1GWASStr,trait2GWASStr)
novelSNPs = novelSNPs[which(unlist(notKnowns)),]
novelSNPs = merge(novelSNPs,overlappingTraits,by="SNP")
save(novelSNPs,file="novelSNPs")


#Get the highest single-trait SNP and all its info in haplotype
topOut = t(sapply(apply(SNPs,1,getProxies),findHighestSNP,out_dat,"pval.outcome"))
colnames(topOut)[colnames(topOut) == 'SNP'] <- 'Top_Outcome_SNP'
topExp = t(sapply(apply(SNPs,1,getProxies),findHighestSNP,exp_dat,"pval.exposure"))
colnames(topExp)[colnames(topExp) == 'SNP'] <- 'Top_Exposure_SNP'


topOut = cbind(topOut,SNP=SNPs$SNP)
topExp = cbind(topExp,SNP=SNPs$SNP)
results = merge(novelSNPs[,c("SNP","CHR","BP","P","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","eaf.exposure","eaf.outcome","beta.exposure","beta.outcome","se.exposure","se.outcome","pval.exposure","pval.outcome","z_exp","z_out","gene","dist","sameDir","Linked_SNPs","Linked_Trait","Linked_Link","nearby_trait","nearby_PMID","nearby_rs")],topOut[,c("SNP","Top_Outcome_SNP","beta.outcome","se.outcome","pval.outcome")],by="SNP")
results = merge(results,topExp[,c("SNP","Top_Exposure_SNP","beta.exposure","se.exposure","pval.exposure")],by="SNP")
fwrite(results,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("Novel_SNPs_",ts,"_FullyAnnot.txt"),sep="\t")
