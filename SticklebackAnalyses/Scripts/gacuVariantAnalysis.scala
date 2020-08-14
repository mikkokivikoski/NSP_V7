#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 16

object variantFiltering {
	
	val gacuVCF=INPUT("result_GITgacuMappingVariantCalling/gacuMpileup/folder1/Gacu_SRR5626529.vcf.gz")
	val gacuNewVCF=INPUT("result_GITgacuMappingVariantCalling/gacuNewMpileup/folder1/GacuNew_SRR5626529.vcf.gz")
	val VCF = Map("gacu" -> gacuVCF, "gacuNew" -> gacuNewVCF)
	
	
	val liftover = INPUT("liftover_fix.awk")
	
	val gasAcuRefIndex = INPUT("Gac-HiC_revised_genome_assembly.fa.fai")
    val gasAcuRefLA10XIndex = INPUT("gacuNEW/3SP2.fasta.fai")
    val gasAcuRefAgp = INPUT("gacuNEW/contigs.agp")
    val gasAcuRefLA10XAgp = INPUT("gacuNEW/3SP2.agp")

	val haplotypes = INPUT("result_gacuAnalyses/haplotypeContigs/table.csv")
	val haplotypeEnds = INPUT("result_gacuAnalyses/searchHaplotypes/table.csv")
	
	val gacuPositiveMask = INPUT("gacuAnalyses/masksInFolder/out/gacu_out1")
	val gacuNewPositiveMask = INPUT("gacuAnalyses/masksInFolder/out/gasAcuLA10X_out1")
	
	val gacuAgp = BashEvaluate(
	    var1=gasAcuRefAgp,
	    command="""
	      awk '{print $6,$7,$8,$4,$5,$1,$2,$3,$9}' OFS="\t" @var1@ > @out1@
	    """
	 )
	
	
	val liftGacuSNPs=BashEvaluate(
	    var1 = liftover,
	    var2 = gacuPositiveMask,
	    var3 = gacuAgp.out1,
	    var4 = gasAcuRefLA10XAgp,
	    var5 = haplotypes,
   	    var6 = haplotypeEnds,
   	    array1 = VCF,
        command = """
		    keys=( $( getarraykeys array1 ) )
	        files=( $( getarrayfiles array1 ) )    
	        gacuInd=$( getarraykeyindex array1 gacu )
	    
		    cd @folder1@
		    bcftools view -H -v snps -e '%QUAL<20 | DP<7 | DP>23' ${files[$gacuInd]} -T @var2@ > gacuvariants  
		    
		    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' gacuvariants > gacuvariantsID
		    
		    grep -v "CONTIG" @var5@ |sed 's/\"//g' > haplotypes.tmp
		    awk '{FS="\t"; row=$1 FS $2 FS $3; print row}' haplotypes.tmp > haplotypes.bed
		    awk -vinverse=1 -f @var1@ @var3@ gacuvariantsID |cut -f1,2,11| awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' > gacuvariantsInContigs
		    bedtools intersect -a haplotypes.bed -b gacuvariantsInContigs -wb > gacuhaplotypeVariants
		    sed -i '1i CONTIG\tSTART\tEND\tVARIANTGONTIG\tVARIANTSTART\tVARIANTEND\tORIGLOCUS' gacuhaplotypeVariants
		    cp gacuhaplotypeVariants @out3@
		    
		    ## PARTIAL HAPLOTYPES, i.e. where variant is in a region whose variant haplotype has been trimmed
		    grep -v "CONTIG" @var6@ |sed 's/\"//g' > gacuremainingHaplotypeRegions.bed
	 	    bedtools intersect -a gacuremainingHaplotypeRegions.bed -b gacuvariantsInContigs -wb > gacuhaplotypeVariants2
	 	    
		    
		    awk -vinverse=1 -f @var1@ @var3@ gacuvariantsID|awk -f @var1@ @var4@ - > gacuvariantsIngacuNewtmp
	 	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS' gacuvariantsID
		    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' gacuvariantsIngacuNewtmp > gacuvariantsIngacuNew
		    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS\tLIFTLOCUS' gacuvariantsIngacuNew
		    cut -f1,2,11 gacuvariantsID > @out1@
		    cut -f1,2,11,12 gacuvariantsIngacuNew > @out2@
	    """
	)
	liftGacuSNPs._keep=false
	//~ liftGacuSNPs._execute="once"
	
  val liftGacuNewSNPs=BashEvaluate(
	    var1 = liftover,
	    var2 = gacuNewPositiveMask,
	    var3 = gacuAgp.out1,
	    var4 = gasAcuRefLA10XAgp,
	    var5 = haplotypes,
	    var6 = haplotypeEnds,
	    array1 = VCF,
	    command = """
	        keys=( $( getarraykeys array1 ) )
	        files=( $( getarrayfiles array1 ) )    
	        gacuLA10XInd=$( getarraykeyindex array1 gacuNew )
	    
	    cd @folder1@
	    bcftools view -H -v snps -e '%QUAL<20 | DP<7 | DP>23' ${files[$gacuLA10XInd]} -T @var2@ > gacuNewvariants 
	    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' gacuNewvariants > gacuNewvariantsID
	    
	    ## FULL HAPLOTYPES, i.e. where one contig has been removed entirely
	    grep -v "CONTIG" @var5@ |sed 's/\"//g' > haplotypes.tmp
	    awk '{FS="\t"; row=$1 FS $2 FS $3; print row}' haplotypes.tmp > haplotypes.bed
	    awk -vinverse=1 -f @var1@ @var4@ gacuNewvariantsID |cut -f1,2,11| awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' > gacuNewvariantsInContigs
	    bedtools intersect -a haplotypes.bed -b gacuNewvariantsInContigs -wb > gacuNewhaplotypeVariants
	    
	    ## PARTIAL HAPLOTYPES, i.e. where variant is in a region whose variant haplotype has been trimmed
	    grep -v "CONTIG" @var6@ |sed 's/\"//g' > gacuNewremainingHaplotypeRegions.bed
 	    bedtools intersect -a gacuNewremainingHaplotypeRegions.bed -b gacuNewvariantsInContigs -wb > gacuNewhaplotypeVariants2
 	    
	    awk -vinverse=1 -f @var1@ @var4@ gacuNewvariantsID|awk -f @var1@ @var3@ - > gacuNewvariantsIngacutmp
	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS' gacuNewvariantsID
	    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' gacuNewvariantsIngacutmp > gacuNewvariantsIngacu
	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS\tLIFTLOCUS' gacuNewvariantsIngacu
	    cut -f1,2,11 gacuNewvariantsID > @out1@
	    cut -f1,2,11,12 gacuNewvariantsIngacu > @out2@
	    
	"""
	) 
	liftGacuNewSNPs._keep=false
	//~ liftGacuNewSNPs._execute="once"
   
   
  val forcedVariantcalling=BashEvaluate(
	    var1 =liftGacuNewSNPs.out2,
	    var3 = liftover,
	    var4 = gacuAgp.out1,
	    var5 = gasAcuRefLA10XAgp,
	    array1 = VCF,
	    command="""
		    keys=( $( getarraykeys array1 ) )
	        files=( $( getarrayfiles array1 ) )    
	        gacuLA10XInd=$( getarraykeyindex array1 gacuNew )
	    
		    cd @folder1@
		    grep -v "CHROM" @var1@ |cut -f1,2|sort -k1,1 -k2,2n > gacuNewingacu
            grep -v "CHROM" @var1@ |cut -f3|sed 's/:/\t/g' > gacuNewingacuNew
            bcftools view -v snps -e '%QUAL<20 | DP<7 | DP>23' ${files[$gacuLA10XInd]} -R gacuNewingacuNew -Oz -o gacuNew.vcf.gz #-e '%QUAL<20 | DP<30 | DP>70'
		    bcftools mpileup --threads 3 -Oz -o tmp.vcf.gz -f Gac-HiC_revised_genome_assembly.fa result_gacuMappingVariantCalling/WGSbwa/folder1/Gacu_SRR5626529_dedup.bam -R gacuNewingacu
		    bcftools index tmp.vcf.gz
		    bcftools call tmp.vcf.gz -R gacuNewingacu -m -Oz -o gacu.vcf.gz
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%MQ0F]\t[%MQ]\t[%DP4]\n' gacu.vcf.gz |sort -k1,1 -k2,2n > gacu.stats
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%MQ0F]\t[%MQ]\t[%DP4]\n' gacuNew.vcf.gz > gacuNew.stats
			awk -F'\t' '{$0=$0"\t"$1":"$2} 1' gacuNew.stats > gacuNew.stats.ID
			awk -vinverse=1 -f @var3@ @var5@ gacuNew.stats.ID | awk -f @var3@ @var4@ - |sort -k1,1 -k2,2n > gacuNew.stats.IDgacu
		    sed -i '1i CHROM\tPOS\tREF\tALT\tQUAL\tDP\tMQ0F\tMQ\tDP4' gacu.stats
		    sed -i '1i CHROM\tPOS\tREF\tALT\tQUAL\tDP\tMQ0F\tMQ\tDP4\tORIGLOCUS' gacuNew.stats.IDgacu
		    cp gacu.stats @out1@
		    cp gacuNew.stats.IDgacu @out2@
	    """
	)
    //~ forcedVariantcalling._execute="once"

	val variantFiles = Map("gacuAll" -> liftGacuSNPs.out1,"gacuLift" -> liftGacuSNPs.out2,"gacuNewAll" -> liftGacuNewSNPs.out1,"gacuNewLift" -> liftGacuNewSNPs.out2, "gacuStats" -> forcedVariantcalling.out1, "gacuNewStats" -> forcedVariantcalling.out2)
	
	
	val combine=REvaluate(
	    table3=liftGacuSNPs.out3,
	    var1=gasAcuRefIndex,
	    var2=gasAcuRefLA10XIndex,
	    var3=liftGacuNewSNPs.folder1,
	    inArray=variantFiles,
	    script="""
	    table.out=data.frame()
		rm('optOut1')
        fig.dir <- get.output(cf, 'optOut1')
		dir.create(fig.dir, recursive=TRUE)
		setwd(fig.dir)
		library(ggplot2)
		library(ggforce)
		OUT=list()
		gacuIndex=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
		gacuNewIndex=read.table(var2,header=F,stringsAsFactors=F,sep="\t")
	    gacuNewhaplotypeVariants=read.table(paste0(var3,"/gacuNewhaplotypeVariants"),header=F,stringsAsFactors=F,sep="\t")
	    gacuNewhaplotypeVariantsContigEnds=read.table(paste0(var3,"/gacuNewhaplotypeVariants2"),header=F,stringsAsFactors=F,sep="\t")
	    
	    gacu=subset(array[["gacuLift"]], select=c(ORIGLOCUS, LIFTLOCUS))
	    gacu=data.frame(gacu,stringsAsFactors=F)
	    gacu.2=subset(array[["gacuAll"]], !(array[["gacuAll"]]$ORIGLOCUS %in% array[["gacuLift"]]$ORIGLOCUS), select=c(ORIGLOCUS))
	    gacu.2=data.frame(gacu.2,LIFTLOCUS="Unassigned",stringsAsFactors=F)
	    gacu=rbind(gacu,gacu.2)
	  
	    gacuNew=subset(array[["gacuNewLift"]], select=c(ORIGLOCUS, LIFTLOCUS))
	    gacuNew=data.frame(gacuNew,stringsAsFactors=F)
	    gacuNew.2=subset(array[["gacuNewAll"]], !(array[["gacuNewAll"]]$ORIGLOCUS %in% array[["gacuNewLift"]]$ORIGLOCUS),select=c(ORIGLOCUS))
	    gacuNew.2=data.frame(gacuNew.2,LIFTLOCUS="Unassigned",stringsAsFactors=F)
	    gacuNew=rbind(gacuNew,gacuNew.2)
	    gacuNew$CALLED="NEW"
	    gacuNew=data.frame(gacuNew, stringsAsFactors=F)
	    print(head(gacuNew))
	    
	    gacuNew[gacuNew$ORIGLOCUS %in% gacu$LIFTLOCUS, "CALLED"]="OLD-NEW"
	    gacu.not.ingacuNew=gacu[gacu$LIFTLOCUS=="Unassigned",]
	    gacu.not.ingacuNew[which(!(unlist(sapply(strsplit(gacu.not.ingacuNew$LIFTLOCUS,split=":"),function(x) x[1])) %in% gacuNewIndex$V1)),"LIFTLOCUS"]="Rem. haplot."
	    gacu.not.ingacuNew=data.frame(gacu.not.ingacuNew,CALLED="OLD",stringsAsFactors=F)
	    colnames(gacu.not.ingacuNew)=c("LIFTLOCUS","ORIGLOCUS","CALLED")
	    gacu.not.called.ingacuNew=gacu[!(gacu$LIFTLOCUS %in% gacuNew$ORIGLOCUS) & gacu$LIFTLOCUS!="Unassigned",] 
	    gacu.not.called.ingacuNew=data.frame(gacu.not.called.ingacuNew,CALLED="OLD",stringsAsFactors=F)
	    print(head(gacuNew))
	    print(head(gacu.not.ingacuNew))
	    print(head(gacu.not.called.ingacuNew))
	    colnames(gacu.not.called.ingacuNew)=c("LIFTLOCUS","ORIGLOCUS","CALLED")
	    
	    gacuNew=rbind(gacuNew,gacu.not.ingacuNew,gacu.not.called.ingacuNew)
	    gacuNew[gacuNew$ORIGLOCUS %in% gacuNewhaplotypeVariants$V7,"CALLED"]=paste0(gacuNew[gacuNew$ORIGLOCUS %in% gacuNewhaplotypeVariants$V7,"CALLED"]," h")
	    gacuNew[gacuNew$ORIGLOCUS %in% gacuNewhaplotypeVariantsContigEnds$V7,"CALLED"]=paste0(gacuNew[gacuNew$ORIGLOCUS %in% gacuNewhaplotypeVariantsContigEnds$V7,"CALLED"]," hc")
	    gacuNew[gacuNew$LIFTLOCUS %in% table3$ORIGLOCUS & gacuNew$CALLED != "OLD-NEW h","CALLED"]=paste0(gacuNew[gacuNew$LIFTLOCUS %in% table3$ORIGLOCUS & gacuNew$CALLED != "OLD-NEW h","CALLED"]," h")

	    
	    gacuNew$gacuNewGROUP=sapply(strsplit(gacuNew$ORIGLOCUS,split=":"),function(x) x[1])
	    gacuNew$gacuGROUP=sapply(strsplit(gacuNew$LIFTLOCUS,split=":"),function(x) x[1])
	    gacuNew$CHANGE=paste0(gacuNew$gacuGROUP,":",gacuNew$gacuNewGROUP,":",gacuNew$CALLED)
	    gacuNew$LENGTH=1
	    OUT[["ALLSNPS"]]=gacuNew
	    total=aggregate(gacuNew$LENGTH,by=list(gacuNew$CHANGE),sum)
		colnames(total)=c("CHANGE","LENGTH")
		total$V6=sapply(strsplit(total$CHANGE,split=":"), function(x) x[1])
		total$V7=sapply(strsplit(total$CHANGE,split=":"), function(x) x[2])
		total$CALL=sapply(strsplit(total$CHANGE,split=":"), function(x) x[3])
		
		gacuNew.3=gacuNew
		gacuNew.3$gacuPOS=sapply(strsplit(gacuNew.3$LIFTLOCUS,split=":"), function(x) as.numeric(x[2]))
		gacuNew.3$gacuNewPOS=sapply(strsplit(gacuNew.3$ORIGLOCUS,split=":"), function(x) as.numeric(x[2]))
		gacuNew.3tmp=gacuNew.3
		gacuNew.3tmp$VARIANTGROUP="--"
		gacuNew.3=gacuNew.3tmp
		
		gacuNew.3[!(grepl("chrM",gacuNew.3$gacuGROUP)) & !(grepl("Un",gacuNew.3$gacuGROUP)) & grepl("CHR",gacuNew.3$gacuNewGROUP) & grepl("OLD-NEW",gacuNew.3$CALLED),"VARIANTGROUP"]="Remained same"
		gacuNew.3[grepl("chrUn:chrUn",gacuNew.3$CHANGE) & grepl("OLD-NEW",gacuNew.3$CALLED),"VARIANTGROUP"]="Remained same" 
		gacuNew.3[grepl("chrUn:CHR", gacuNew.3$CHANGE),"VARIANTGROUP"]="New contig" # & !(grepl("Gacu",gacuNew.3$CALLED))
		gacuNew.3[grepl("Unassigned:CHR", gacuNew.3$CHANGE),"VARIANTGROUP"]="New contig" # & !(grepl("Gacu",gacuNew.3$CALLED))
		gacuNew.3[!(grepl("chrM",gacuNew.3$gacuGROUP)) & !(grepl("Un",gacuNew.3$gacuGROUP)) & grepl("CHR",gacuNew.3$gacuNewGROUP) & !(grepl("OLD",gacuNew.3$CALLED)) & grepl("h",gacuNew.3$CALLED) & !(grepl("hc",gacuNew.3$CALLED)),"VARIANTGROUP"]="Haplotype removed"
		gacuNew.3[!(grepl("chrM",gacuNew.3$gacuGROUP)) & !(grepl("Un",gacuNew.3$gacuGROUP)) & grepl("CHR",gacuNew.3$gacuNewGROUP) & !(grepl("OLD",gacuNew.3$CALLED)) & grepl("hc",gacuNew.3$CALLED),"VARIANTGROUP"]="Haplotype trimmed"
		gacuNew.3[!(grepl("chrM",gacuNew.3$gacuGROUP)) & !(grepl("Un",gacuNew.3$gacuGROUP)) & !(grepl("CHR",gacuNew.3$gacuNewGROUP)),"VARIANTGROUP"]="Locus removed from Gacu_New" 
		
		gacuNew.3[!(grepl("chrM",gacuNew.3$gacuGROUP)) & !(grepl("Un",gacuNew.3$gacuGROUP)) & grepl("CHR",gacuNew.3$gacuNewGROUP) & !(grepl("NEW",gacuNew.3$CALLED)),"VARIANTGROUP"]="Variant not found anymore"
		gacuNew.3[grepl("Un",gacuNew.3$gacuGROUP) & grepl("Un",gacuNew.3$gacuNewGROUP) & !(grepl("OLD",gacuNew.3$CALLED)),"VARIANTGROUP"]="Unassigned in both Gacu_New"
		gacuNew.3[grepl("Un",gacuNew.3$gacuGROUP) & grepl("Un",gacuNew.3$gacuNewGROUP) & !(grepl("NEW",gacuNew.3$CALLED)),"VARIANTGROUP"]="Unassigned in both Gacu"
		gacuNew.3[!(grepl("OLD",gacuNew.3$CALLED)) & grepl("NEW",gacuNew.3$CALLED) & gacuNew.3$VARIANTGROUP=="--","VARIANTGROUP"] = "New variant for some reason"
		gacuNew.3[grepl("Un",gacuNew.3$gacuGROUP) & grepl("Rem. haplot",gacuNew.3$gacuNewGROUP),"VARIANTGROUP"]="Unassinged and removed"
		gacuNew.3[grepl("chrM",gacuNew.3$gacuGROUP),"VARIANTGROUP"]="Mitochondrial"
		
		total.3=aggregate(gacuNew.3$LENGTH,by=list(gacuNew.3$VARIANTGROUP),sum)
		total.3=aggregate(subset(gacuNew.3,!(grepl("chrXIX",gacuGROUP)) & !(grepl("CHR19",gacuNewGROUP)))$LENGTH,by=list(subset(gacuNew.3,!(grepl("chrXIX",gacuGROUP)) & !(grepl("CHR19",gacuNewGROUP)))$VARIANTGROUP),sum)
		OUT[["TOTAL3"]]=total.3
		
		
	    array.out=OUT
		
		"""
	)
		

}
