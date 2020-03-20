#Load the following on an interactive session before using this script:
#module load JAGS
#module load OpenBUGS
#module load plink/1.90Beta4.5
#module load bedtools2/2.25.0
#module load R


library(coloc)
library(data.table)
library(TwoSampleMR)
library(miscF)
library(rlist)
source("/project/voight_selscan/ksiewert/Migraine_Project/Scripts/bivariate_scan.R")
source("trait1_config.R")

ts = paste0(trait1,trait2,sep = "")

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

har$pval.outcome = pmax(har$pval.outcome,1e-16) #Need to take pmax so don't get underflow with estimating p
har$pval.exposure = pmax(har$pval.exposure,1e-16)

har$z_out =  sqrt(qchisq(1 - har$pval.outcome,1))*sign(har$beta.outcome)
har$z_exp =  sqrt(qchisq(1 - har$pval.exposure,1))*sign(har$beta.exposure)

write.table(har[c("SNP","effect_allele.exposure","other_allele.exposure","eaf.exposure","beta.exposure","beta.exposure","beta.outcome","pval.outcome","pval.exposure","se.outcome","se.exposure","z_out","z_exp")],paste0("AllRs",ts,".txt"),quote = FALSE,row.names=FALSE)


#The har_chompos strategy should work, but I am worrier it will eat up a lot of time debugging. So just going to read the gwas summary stats into R after the bivar scans are done
#add chrom and Pos for hg19 to the har dataframe. Need this for coloc steps later


#read the chr pos rsID file into R
#rsID2chrPos = fread(file="/project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/at_least_AD/rsID2ChrPos_table_dbSNP150_hg19.txt", sep="\t", header=TRUE)


#inner join the rsID2chrPos and the har data frames
#har_chrpos <- merge(rsID2chrPos,har, by = "SNP")

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


##Filter for loci that are at least nominally significant
o = order(d$P,decreasing=FALSE)
d = d[o,]
d = d[which(d$P<1e-7),] #contains all clustered SNP with p<1e-7 whether or not novel
write.table(d, quote = FALSE,row.names=FALSE,file=paste0("AllSigSNPs_",ts,"_r2Indepdent.txt"),sep="\t") #NOTE: May not be 1MB from each other, but will be r2<.2

#read in the AllSigSNPs_",ts,"_r2Indepdent.txt" file from the config file
#d <- AllIndSigSNPs


#First, return hits from the GWAS catalog within 500kb, no matter what trait they correspond to
GWAS = read.table("/project/voight_datasets/GWAS/15_GWASCatalog/Jan_2019_version/GWAShg37.bed",header=FALSE,sep="\t",colClasses = "character",fill=TRUE,quote = "",stringsAsFactors = FALSE)
colnames(GWAS)=c("chr","coord","coord+1","PMID","rs","trait","P")
GWAS = GWAS[which(as.numeric(GWAS$P)<=5e-8),] #only want significant SNPs from GWAS catalog
GWASlist = split(GWAS,GWAS$chr)
GWASoverlap =t(apply(d,1,findNearGWASCatalogHits)) #finds gwas hits within 500kb
colnames(GWASoverlap) = c("nearby_trait","nearby_PMID","nearby_rs","nearby_coord")


#We want to remove the GWASCatalog filters for this version of the bivar scan
#Next, want to filter out any SNPs that have a relevant GWAS result within 500KB bases, because this means they are not novel
GWASannot = cbind(d,GWASoverlap)
#g = unlist(lapply(GWASoverlap[,1],noDoubleMatch,trait1GWASStr,trait2GWASStr)) #Checks if any of the GWAS catalog traits are for the traits of interest (using strings from config.R file)
#novel = GWASannot[which(g),]

#NOTE we are expecting the trait1 to be the exposure variable MAKE SURE THIS IS THE CASE IN YOUR config.R file
Trait1FocusSNPs = GWASannot[which(GWASannot$pval.outcome<.005 & GWASannot$pval.exposure<1e-6),]
#pull out ones with indiv p-values that are nominally, but not genome-wide significant for trait 2, and near significant for trait 1

write.table(Trait1FocusSNPs, quote = FALSE,row.names=FALSE,file=paste0("Trait1Focusedhits_",ts,"_bivarresults.txt"),sep="\t")
write.table(Trait1FocusSNPs$SNP,row.names=FALSE,col.names=FALSE,quote=FALSE,file="Trait1FocusSNPrsNums.txt",sep="\t")


#See what other traits are within r^2=.2 of novel GWAS hits
system(paste0("plink --bfile /project/voight_selscan/ksiewert/CardioMetaAnalysis/LDL_CHD_Bivar/LDClump/PlinkFilesOnlyRs/mergedBed --keep /project/voight_GWAS/wbone/neuro_degenerative_and_cardiometabolic_Bivariate_Scans/AD_bivariate_scan_code/EUR.final.plink --show-tags Trait1FocusSNPrsNums.txt --out Trait1SNPTags  --tag-kb 1000 --list-all --tag-r2 0.2"))
SNPs = read.table("Trait1SNPTags.tags.list", header=TRUE,stringsAsFactors = FALSE)
overlappingTraits = t(apply(SNPs,1,f))
colnames(overlappingTraits)=c("SNP","Linked_SNPs","Linked_Trait","Linked_Link")


#Remove this filter
#notKnowns = lapply(overlappingTraits[,"Linked_Trait"],noDoubleMatch,trait1GWASStr,trait2GWASStr)
#novelSNPs = novelSNPs[which(unlist(notKnowns)),]
Trait1FocusSNPs = merge(Trait1FocusSNPs,overlappingTraits,by="SNP")
save(Trait1FocusSNPs,file="Trait1FocusSNPs")


#Get the highest single-trait SNP and all its info in haplotype
topOut = t(sapply(apply(SNPs,1,getProxies),findHighestSNP,out_dat,"pval.outcome"))
colnames(topOut)[colnames(topOut) == 'SNP'] <- 'Top_Outcome_SNP'
topExp = t(sapply(apply(SNPs,1,getProxies),findHighestSNP,exp_dat,"pval.exposure"))
colnames(topExp)[colnames(topExp) == 'SNP'] <- 'Top_Exposure_SNP'


topOut = cbind(topOut,SNP=SNPs$SNP)
topExp = cbind(topExp,SNP=SNPs$SNP)
results = merge(Trait1FocusSNPs[,c("SNP","CHR","BP","P","effect_allele.exposure","other_allele.exposure","effect_allele.outcome","other_allele.outcome","eaf.exposure","eaf.outcome","beta.exposure","beta.outcome","se.exposure","se.outcome","pval.exposure","pval.outcome","z_exp","z_out","gene","dist","sameDir","Linked_SNPs","Linked_Trait","Linked_Link","nearby_trait","nearby_PMID","nearby_rs")],topOut[,c("SNP","Top_Outcome_SNP","beta.outcome","se.outcome","pval.outcome")],by="SNP")
results = merge(results,topExp[,c("SNP","Top_Exposure_SNP","beta.exposure","se.exposure","pval.exposure")],by="SNP")
fwrite(results,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("Trait1Focus_SNPs_",ts,"_FullyAnnot.txt"),sep="\t")

#Time to run coloc
#filter the input file to only include SNPs in the area we are interested in colocing

#read the gwas files into R
trait1_gwas_data = fread(file=expPath, sep="\t", header=TRUE)
trait2_gwas_data = fread(file=outPath, sep="\t", header=TRUE)

# initialize a list of the summary results for each coloc run
results_coloc_summary_list <- list()


for (hit in 1:nrow(results)){ 
    #go through each locus and coloc

    leadSNP <- results[hit,"SNP"] #lead SNP that we will base the names of files and the analysis will be centered around
    print(leadSNP)

    chrom = results[hit,"CHR"]
    colocStart = results[hit,"BP"] - 250000
    colocStop = results[hit,"BP"] + 250000

    colocStart_str = as.character(colocStart)
    colocStop_str = as.character(colocStop)
    chrom_str = as.character(chrom)

    print(colocStart_str)
    print(colocStop_str)
    print(chrom_str)

    out_prefix = paste(leadSNP,chrom_str,colocStart_str,colocStop_str, sep="_")
    ########################### NOT USING THIS ######################################
    #need to grab it all from har_chrpos
    #grab the snps that are within the start stop and on the correct chromosome from the the har_chrompos dataframe

    #coloc_region = har_chrpos$chrom == chrom & har_chrpos$BP >= colocStart & har_chrpos$BP <= colocStop

    ########################### NOT USING THIS ######################################

    #grab the snps that are within the start stop and on the correct chromosome from the trait files

    trait1_region = trait1_gwas_data[trait1_gwas_data[[trait1_CHRcol]] == chrom & trait1_gwas_data[[trait1_BPcol]] >= colocStart & trait1_gwas_data[[trait1_BPcol]] <= colocStop,]

    trait2_region = trait2_gwas_data[trait2_gwas_data[[trait2_CHRcol]] == chrom & trait2_gwas_data[[trait2_BPcol]] >= colocStart & trait2_gwas_data[[trait2_BPcol]] <= colocStop,]

    
    #Change column headers to avoid overlaps in coloc input file
    #trait1_region column header changes
    names(trait1_region)[names(trait1_region) == trait1_Pcol] <- "trait1_Pcol"
    names(trait1_region)[names(trait1_region) == trait1_Ncol] <- "trait1_Ncol"
    names(trait1_region)[names(trait1_region) == trait1_MAFcol] <- "trait1_MAFcol"

    #Change column headers to avoid overlaps in coloc input file
    #trait1_region column header changes
    names(trait2_region)[names(trait2_region) == trait2_Pcol] <- "trait2_Pcol"
    names(trait2_region)[names(trait2_region) == trait2_Ncol] <- "trait2_Ncol"
    names(trait2_region)[names(trait2_region) == trait2_MAFcol] <- "trait2_MAFcol"

    print("merging the trait and eqtl data on rsID")
    #merge the trait files on the rsID or SNP name
    colocInputFile = merge(trait1_region, trait2_region, by.x=trait1_SNPcol, by.y=trait2_SNPcol)
    #print(colnames(colocInputFile))

    ######################DON'T THINK I actually need this######################
    #filter rows that dont have either A1_trait1 == A1_trait2 and A2_trait1 ==A2_trait2 OR A1_trait1 == A2_trait2 and A2_trait1 == A1_trait2
    #colocInputFile = colocInputFile[eGeneTissue_region$ref == eGeneTissue_region$A1_eqtl & eGeneTissue_region$alt == eGeneTissue_region$A2_eqtl | eGeneTissue_region$alt == eGeneTissue_region$A1_eqtl & eGeneTissue_region$ref == eGeneTissue_region$A2_eqtl,]
    ######################DON'T THINK I actually need this######################

    #remove any rows with NAs
    colocInputFile = colocInputFile[complete.cases(colocInputFile), ]

    #removed any rows with odd values for the trait P-values or the MAF
    colocInputFile <- colocInputFile[colocInputFile$trait1_Pcol < 1 & colocInputFile$trait1_Pcol > 0 & colocInputFile$trait2_Pcol < 1 & colocInputFile$trait2_Pcol > 0 & colocInputFile$trait1_MAFcol < 1 & colocInputFile$trait1_MAFcol > 0,]

    #write colocInputFile to file for making locus zoom plots
    colocInputFile_outputStr = paste(out_prefix,"coloc_input_data.txt",sep="_")
    write.table(colocInputFile, file= colocInputFile_outputStr, sep="\t", row.names=FALSE, quote=FALSE)

    print("Running coloc")
    #run coloc taking into account the types of traits
    if (trait1Type == "cc" & trait2Type == "quant" ){

        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile$trait1_Pcol, N=colocInputFile$trait1_Ncol, type=trait1Type, s=trait1Prop), dataset2=list(pvalues=colocInputFile$trait2_Pcol, N=colocInputFile$trait2_Ncol, type=trait2Type),MAF=colocInputFile$trait1_MAFcol)

    } else if (trait1Type == "cc" & trait2Type == "cc" ){

       coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile$trait1_Pcol, N=colocInputFile$trait1_Ncol, type=trait1Type, s=trait1Prop), dataset2=list(pvalues=colocInputFile$trait2_Pcol, N=colocInputFile$trait2_Ncol, type=trait2Type, s=trait2Prop),MAF=colocInputFile$trait1_MAFcol) 

    } else if (trait1Type == "quant" & trait2Type == "cc"){

        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile$trait1_Pcol, N=colocInputFile$trait1_Ncol, type=trait1Type), dataset2=list(pvalues=colocInputFile$trait2_Pcol, N=colocInputFile$trait2_Ncol, type=trait2Type, s=trait2Prop),MAF=colocInputFile$trait1_MAFcol)

    }else {
        #for when both traits are "quant"
        coloc_results <- coloc.abf(dataset1=list(pvalues=colocInputFile$trait_Pcol, N=colocInputFile$trait_Ncol, type=traitType), dataset2=list(pvalues=colocInputFile$trait2_Pcol, N=colocInputFile$trait2_Ncol, type=trait2Type),MAF=colocInputFile$trait_MAFcol_str)

    }

    #prepare useful outputs
    coloc_results_summary = coloc_results$summary
    coloc_results_full = coloc_results$results

    #add leadSNP to coloc results and prep to add to the results dataframe
    #coloc_results_summary <- append(coloc_results_summary, leadSNP)
    #names(coloc_results_summary)[names(coloc_results_summary) == ""] <- "SNP"

    #add this locus' coloc results to the list of coloc results    
    #results_coloc_summary_list <- list.append(results_coloc_summary_list, coloc_results_summary )

    #calculate pp4 / pp3 + pp4
    PP3andPP4 = coloc_results_summary[5] + coloc_results_summary[6]

    pp4_conditional = coloc_results_summary[6] / PP3andPP4

    #prep coloc output strings
    coloc_results_summary_outputStr = paste(out_prefix,"coloc_results_summary.txt",sep="_")
    coloc_results_full_outputStr = paste(out_prefix,"coloc_results_full.txt",sep="_")
    coloc_results_pp4_cond_outputStr = paste(out_prefix,"coloc_results_pp4_cond.txt",sep="_")

    #add leadSNP to coloc results and prep to add to the results dataframe
    leadSNP <- as.character(leadSNP)
    coloc_results_summary <- append(coloc_results_summary, leadSNP)
    names(coloc_results_summary)[names(coloc_results_summary) == ""] <- "SNP"

    #add this locus' coloc results to the list of coloc results    
    results_coloc_summary_list <- list.append(results_coloc_summary_list, coloc_results_summary )

    #write to file
    write.table(coloc_results_summary, file=coloc_results_summary_outputStr, sep="\t", row.names=TRUE, quote=FALSE)
    write.table(coloc_results_full, file=coloc_results_full_outputStr, sep="\t", row.names=FALSE, quote=FALSE)
    write.table(pp4_conditional, file=coloc_results_pp4_cond_outputStr, sep="\t", row.names=FALSE, quote=FALSE)

}

#once all the lines have been coloc"ed" make a dataframe of the coloc result then add it to the bivar results file
results_coloc_summary_df <- as.data.frame(do.call(rbind,results_coloc_summary_list))

#merge coloc summnary results with bivariate scan results

results<- merge(results,results_coloc_summary_df,by="SNP")

#write to file
fwrite(results,row.names=FALSE,col.names=TRUE,quote=FALSE,file=paste0("Trait1Focus_SNPs_",ts,"_FullyAnnot_with_coloc_summary.txt"),sep="\t")
