#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 16


object variantFiltering {
	val liftover = INPUT("liftover_fix.awk")
	val v6agp = INPUT("NSP_V6.agp")
	val v6ref = INPUT("NSP_V6.fasta")
	val v6Index = INPUT("NSP_V6.fasta.fai")
	val v7agp = INPUT("V7.agp")
	val v7Index = INPUT("V7.fasta.fai")
	val haplotypes = INPUT("result_ppunAnalyses/haplotypeContigs/table.csv")
	val haplotypeEnds = INPUT("result_ppunAnalyses/searchHaplotypes/table.csv")
    val referenceGenomes = INPUT("folder with reference sequences")
    val individuals = Map("FINPYO33"->"FIN-PYO-33","FINPYO22"->"FIN-PYO-22","FINPYO20"->"FIN-PYO-20","FINPYO21"->"FIN-PYO-21","FINPYO0"->"FIN-PYO-0")
    val vcfFolder = INPUT("result_mapping/vcfInFolder/out")
    val filterVariants = NamedMap[BinaryFolder]("filterVariants")
    val mask = INPUT("result_ppunAnalyses/masksInFolder/out") 
  

	for((ind, indName) <- individuals){
	
	    filterVariants(ind) = BashEvaluate(
	        var2=vcfFolder,
	        var3=referenceGenomes,
	        var4=mask,
	        param1=indName,
		    command="""
		     for i in v6 v7; do
		         bcftools filter -g 20 @var2@/"$i"_@param1@_filtered.vcf| bcftools view -M2 -V mnps,indels -T @var4@/"$i"_out1 -Oz -o @folder1@/"$i"_@param1@_clean.vcf.gz
		         bcftools index @folder1@/"$i"_@param1@_clean.vcf.gz
		     done
		    """
		).folder1
	}
	val vcfInFolder=Array2Folder(
	  in = filterVariants,
	  fileMode = ".",
	  keyCol = "Key",
	  link = false
	)  



	

	val liftv6SNPs=BashEvaluate(
	    var1 = liftover,
	    var2 = vcfInFolder.out,
	    var3 = v6agp,
	    var4 = v7agp,
	    var5 = haplotypes,
   	    var6 = haplotypeEnds,
        command = """
	    cd @folder1@
	    bcftools view -H -v snps -e '%QUAL<20 | DP<30 | DP>70' @var2@/v6_FIN-PYO-0_clean.vcf.gz > v6variants 
	    
	    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' v6variants > v6variantsID
	    
	    grep -v "CONTIG" @var5@ |sed 's/\"//g' > haplotypes.tmp
	    awk '{FS="\t"; row=$1"|quiver|pilon" FS $2 FS $3; print row}' haplotypes.tmp > haplotypes.bed
	    awk -vinverse=1 -f @var1@ @var3@ v6variantsID |cut -f1,2,11| awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' > v6variantsInContigs
	    bedtools intersect -a haplotypes.bed -b v6variantsInContigs -wb > v6haplotypeVariants
	    sed -i '1i CONTIG\tSTART\tEND\tVARIANTGONTIG\tVARIANTSTART\tVARIANTEND\tORIGLOCUS' v6haplotypeVariants
	    cp v6haplotypeVariants @out3@
	    
	    ## PARTIAL HAPLOTYPES, i.e. where variant is in a region whose variant haplotype has been trimmed
	    grep -v "CONTIG" @var6@ |sed 's/\"//g' > v6remainingHaplotypeRegions.bed
 	    bedtools intersect -a v6remainingHaplotypeRegions.bed -b v6variantsInContigs -wb > v6haplotypeVariants2
 	    
	    
	    awk -vinverse=1 -f @var1@ @var3@ v6variantsID|awk -f @var1@ @var4@ - > v6variantsInv7tmp
 	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS' v6variantsID
	    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' v6variantsInv7tmp > v6variantsInv7
	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS\tLIFTLOCUS' v6variantsInv7
	    cut -f1,2,11 v6variantsID > @out1@
	    cut -f1,2,11,12 v6variantsInv7 > @out2@
	    
	    """
	)
	liftv6SNPs._keep=false
	
	val liftv7SNPs=BashEvaluate(
	    var1 = liftover,
	    var2 = vcfInFolder.out,
	    var3 = v6agp,
	    var4 = v7agp,
	    var5 = haplotypes,
	    var6 = haplotypeEnds,
	    command = """
	    cd @folder1@
	    bcftools view -H -v snps -e '%QUAL<20 | DP<30 | DP>70' @var2@/v7_FIN-PYO-0_clean.vcf.gz > v7variants 
	    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' v7variants > v7variantsID
	    
	    ## FULL HAPLOTYPES, i.e. where one contig has been removed entirely
	    grep -v "CONTIG" @var5@ |sed 's/\"//g' > haplotypes.tmp
	    awk '{FS="\t"; row=$1"|quiver|pilon" FS $2 FS $3; print row}' haplotypes.tmp > haplotypes.bed
	    awk -vinverse=1 -f @var1@ @var4@ v7variantsID |cut -f1,2,11| awk 'BEGIN{FS=OFS="\t"} {$2 = $2 OFS $2} 1' > v7variantsInContigs
	    bedtools intersect -a haplotypes.bed -b v7variantsInContigs -wb > v7haplotypeVariants
	    
	    ## PARTIAL HAPLOTYPES, i.e. where variant is in a region whose variant haplotype has been trimmed
	    grep -v "CONTIG" @var6@ |sed 's/\"//g' > v7remainingHaplotypeRegions.bed
 	    bedtools intersect -a v7remainingHaplotypeRegions.bed -b v7variantsInContigs -wb > v7haplotypeVariants2
 	    
	    awk -vinverse=1 -f @var1@ @var4@ v7variantsID|awk -f @var1@ @var3@ - > v7variantsInv6tmp
	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS' v7variantsID
	    awk -F'\t' '{$0=$0"\t"$1":"$2} 1' v7variantsInv6tmp > v7variantsInv6
	    sed -i '1i CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\tORIGLOCUS\tLIFTLOCUS' v7variantsInv6
	    cut -f1,2,11 v7variantsID > @out1@
	    cut -f1,2,11,12 v7variantsInv6 > @out2@
	    
	"""
	) 
	liftv7SNPs._keep=false
	val forcedVariantcalling=BashEvaluate(
	    var1 =liftv7SNPs.out2,
	    var2 = vcfInFolder.out,
	    var3 = liftover,
	    var4 = v7agp,
	    var5 = v6agp,
	    var6 = v6ref,
	    command="""
		    cd @folder1@
		    grep -v "CHROM" @var1@ |cut -f1,2|sort -k1,1 -k2,2n > v7inv6
            grep -v "CHROM" @var1@ |cut -f3|sed 's/:/\t/g' > v7inv7
            bcftools view -v snps -e '%QUAL<20 | DP<30 | DP>70' @var2@/v7_FIN-PYO-0_clean.vcf.gz -R v7inv7 -Oz -o v7.vcf.gz 
		    bcftools mpileup --threads 3 -Oz -o tmp.vcf.gz -f @var6@ result_mapping/rmdup_FINPYO0/folder1/v6_FIN-PYO-0_dedup.bam -R v7inv6
		    bcftools index tmp.vcf.gz
		    bcftools call tmp.vcf.gz -R v7inv6 -m -Oz -o v6.vcf.gz
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%MQ0F]\t[%MQ]\t[%DP4]\n' v6.vcf.gz |sort -k1,1 -k2,2n > v6.stats
			bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%MQ0F]\t[%MQ]\t[%DP4]\n' v7.vcf.gz > v7.stats
			awk -F'\t' '{$0=$0"\t"$1":"$2} 1' v7.stats > v7.stats.ID
			awk -vinverse=1 -f @var3@ @var4@ v7.stats.ID | awk -f @var3@ @var5@ - |sort -k1,1 -k2,2n > v7.stats.IDv6
		    sed -i '1i CHROM\tPOS\tREF\tALT\tQUAL\tDP\tMQ0F\tMQ\tDP4' v6.stats
		    sed -i '1i CHROM\tPOS\tREF\tALT\tQUAL\tDP\tMQ0F\tMQ\tDP4\tORIGLOCUS' v7.stats.IDv6
		    cp v6.stats @out1@
		    cp v7.stats.IDv6 @out2@
	    """
	)
	forcedVariantcalling._execute="once"
	val variantFiles = Map("v6All" -> liftv6SNPs.out1,"v6Lift" -> liftv6SNPs.out2,"v7All" -> liftv7SNPs.out1,"v7Lift" -> liftv7SNPs.out2, "v6Stats" -> forcedVariantcalling.out1, "v7Stats" -> forcedVariantcalling.out2)
	
	
	/////////
	val combine=REvaluate(
	    table3=liftv6SNPs.out3,
	    var1=v6Index,
	    var2=v7Index,
	    var3=liftv7SNPs.folder1,
	    inArray=variantFiles,
	    script="""
	    table.out=data.frame()
		rm('optOut1')
        fig.dir <- get.output(cf, 'optOut1')
		dir.create(fig.dir, recursive=TRUE)
		setwd(fig.dir)
		library(ggplot2)
		OUT=list()
		v6Index=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
		v7Index=read.table(var2,header=F,stringsAsFactors=F,sep="\t")
	    v7haplotypeVariants=read.table(paste0(var3,"/v7haplotypeVariants"),header=F,stringsAsFactors=F,sep="\t")
	    v7haplotypeVariantsContigEnds=read.table(paste0(var3,"/v7haplotypeVariants2"),header=F,stringsAsFactors=F,sep="\t")
	    
	    v6=subset(array[["v6Lift"]], select=c(ORIGLOCUS, LIFTLOCUS))
	    v6=data.frame(v6,stringsAsFactors=F)
	    v6.2=subset(array[["v6All"]], !(array[["v6All"]]$ORIGLOCUS %in% array[["v6Lift"]]$ORIGLOCUS), select=c(ORIGLOCUS))
	    v6.2=data.frame(v6.2,LIFTLOCUS="Unassigned",stringsAsFactors=F)
	    v6=rbind(v6,v6.2)
	    
	    
	    v7=subset(array[["v7Lift"]], select=c(ORIGLOCUS, LIFTLOCUS))
	    v7=data.frame(v7,stringsAsFactors=F)
	    v7.2=subset(array[["v7All"]], !(array[["v7All"]]$ORIGLOCUS %in% array[["v7Lift"]]$ORIGLOCUS),select=c(ORIGLOCUS))
	    v7.2=data.frame(v7.2,LIFTLOCUS="Unassigned",stringsAsFactors=F)
	    v7.2[which(unlist(sapply(strsplit(v7.2$ORIGLOCUS,split=":"),function(x) x[1]))=="LG19" & unlist(sapply(strsplit(v7.2$ORIGLOCUS,split=":"),function(x) as.numeric(x[2])))>=5595611 & unlist(sapply(strsplit(v7.2$ORIGLOCUS,split=":"),function(x) as.numeric(x[2]))) <=7233324),"LIFTLOCUS"]="LG19 re-assembly"
	    v7=rbind(v7,v7.2)
	    v7$CALLED="ver. 7"
	    v7=data.frame(v7, stringsAsFactors=F)
	    
	    v7[v7$ORIGLOCUS %in% v6$LIFTLOCUS, "CALLED"]="ver. 6-ver. 7"
	    v6.not.inv7=v6[v6$LIFTLOCUS=="Unassigned",]
	    v6.not.inv7[which(!(unlist(sapply(strsplit(v6.not.inv7$LIFTLOCUS,split=":"),function(x) x[1])) %in% v7Index$V1)),"LIFTLOCUS"]="Rem. haplot."
	    v6.not.inv7=data.frame(v6.not.inv7,CALLED="ver. 6",stringsAsFactors=F)
	    colnames(v6.not.inv7)=c("LIFTLOCUS","ORIGLOCUS","CALLED")
	    v6.not.called.inv7=v6[!(v6$LIFTLOCUS %in% v7$ORIGLOCUS) & v6$LIFTLOCUS!="Unassigned",] 
	    v6.not.called.inv7=data.frame(v6.not.called.inv7,CALLED="ver. 6",stringsAsFactors=F)
	    colnames(v6.not.called.inv7)=c("LIFTLOCUS","ORIGLOCUS","CALLED")
	    
	    v7=rbind(v7,v6.not.inv7,v6.not.called.inv7)
	    v7[v7$ORIGLOCUS %in% v7haplotypeVariants$V7,"CALLED"]=paste0(v7[v7$ORIGLOCUS %in% v7haplotypeVariants$V7,"CALLED"]," h")
	    v7[v7$ORIGLOCUS %in% v7haplotypeVariantsContigEnds$V7,"CALLED"]=paste0(v7[v7$ORIGLOCUS %in% v7haplotypeVariantsContigEnds$V7,"CALLED"]," hc")
	    v7[v7$LIFTLOCUS %in% table3$ORIGLOCUS & v7$CALLED != "ver. 6-ver. 7 h","CALLED"]=paste0(v7[v7$LIFTLOCUS %in% table3$ORIGLOCUS & v7$CALLED != "ver. 6-ver. 7 h","CALLED"]," h")
	    
	    #LG19 inversion: LG19:5595611-7233324 
	    
	    v7$v7GROUP=sapply(strsplit(v7$ORIGLOCUS,split=":"),function(x) x[1])
	    v7$v6GROUP=sapply(strsplit(v7$LIFTLOCUS,split=":"),function(x) x[1])
	    v7$CHANGE=paste0(v7$v6GROUP,":",v7$v7GROUP,":",v7$CALLED)
	    v7$LENGTH=1
  	    OUT[["ALLSNPS"]]=v7


		v7.2=v7
		v7.2$v6POS=sapply(strsplit(v7.2$LIFTLOCUS,split=":"), function(x) as.numeric(x[2]))
		v7.2$v7POS=sapply(strsplit(v7.2$ORIGLOCUS,split=":"), function(x) as.numeric(x[2]))
		v7.4=subset(v7.2, v6GROUP == "LG12" & v6POS <= 25*(10^6) |v7GROUP == "LG12" & v7POS<=20*(10^6))
		v7.2tmp=subset(v7.2, v6GROUP != "LG12" & v7GROUP != "LG12" & v6GROUP != "LG19" & v7GROUP != "LG19" & grepl("LG",CHANGE))
		v7.2tmp=rbind(v7.2tmp,subset(v7.2, v6GROUP == "LG12" & v6POS > 25*(10^6) |v7GROUP == "LG12" & v7POS>20*(10^6)))
		v7.2tmp=rbind(v7.2tmp,subset(v7.2, v6GROUP == "LG19" & !(v6POS %in% seq(5500000,7500000)) |v7GROUP == "LG19" & !(v7POS %in% seq(5500000,7500000))))
		v7.2tmp$VARIANTGROUP="--"
		v7.2=v7.2tmp
		v7.2[grepl("LG",v7.2$v6GROUP) & grepl("LG",v7.2$v7GROUP) & grepl("ver. 6-ver. 7",v7.2$CALLED),"VARIANTGROUP"]="Remained same"
		v7.2[grepl("Unassigned:LG", v7.2$CHANGE),"VARIANTGROUP"]="New contig"
		v7.2[grepl("LG",v7.2$v6GROUP) & grepl("LG",v7.2$v7GROUP) & !(grepl("ver. 6",v7.2$CALLED)) & grepl("h",v7.2$CALLED) & !(grepl("hc",v7.2$CALLED)),"VARIANTGROUP"]="Haplotype removed"
		v7.2[grepl("LG",v7.2$v6GROUP) & grepl("LG",v7.2$v7GROUP) & !(grepl("ver. 6",v7.2$CALLED)) & grepl("hc",v7.2$CALLED),"VARIANTGROUP"]="Haplotype trimmed"
		v7.2[grepl("LG",v7.2$v6GROUP) & !(grepl("LG",v7.2$v7GROUP)),"VARIANTGROUP"]="Locus removed from ver. 7"
		v7.2[grepl("LG",v7.2$v6GROUP) & grepl("LG",v7.2$v7GROUP) & !(grepl("ver. 7",v7.2$CALLED)),"VARIANTGROUP"]="Variant not found anymore"
		v7.2[!(grepl("ver. 6",v7.2$CALLED)) & grepl("ver. 7",v7.2$CALLED) & v7.2$VARIANTGROUP=="--","VARIANTGROUP"] = "New variant for some reason"
		
		
		total.2=aggregate(v7.2$LENGTH,by=list(v7.2$VARIANTGROUP),sum)
		OUT[["TOTAL2"]]=total.2
		
		
	    v6.stats=array[["v6Stats"]]
        v6.stats$LOCUS=paste0(v6.stats$CHROM,":",v6.stats$POS)
        v6.stats$REFD=sapply(strsplit(v6.stats$DP4,split=","), function(x) sum(as.numeric(x[1]),as.numeric(x[2])))
		v6.stats$ALTD=sapply(strsplit(v6.stats$DP4,split=","), function(x) sum(as.numeric(x[3]),as.numeric(x[4])))
		v6.stats$vAD=v6.stats$ALTD/v6.stats$DP
		v6.stats=v6.stats[!(duplicated(v6.stats$LOCUS)),]
        v6.stats.a=subset(v6.stats, !(CHROM =="LG12" & POS <= 25*(10^6)) & !(CHROM =="LG19" & POS %in% seq(5500000,7500000)))
		
		v7.stats=array[["v7Stats"]]
		v7.stats$LOCUS=paste0(v7.stats$CHROM,":",v7.stats$POS)
		v7.stats$REFD=sapply(strsplit(v7.stats$DP4,split=","), function(x) sum(as.numeric(x[1]),as.numeric(x[2])))
		v7.stats$ALTD=sapply(strsplit(v7.stats$DP4,split=","), function(x) sum(as.numeric(x[3]),as.numeric(x[4])))
		v7.stats$vAD=v7.stats$ALTD/v7.stats$DP
		v7.stats.a=subset(v7.stats, LOCUS %in% v6.stats.a$LOCUS)
	    
	    
	    
		v6.stats.a$GROUP=""
		v7.stats.a$GROUP=""
		v6.stats.a$GROUP2=""
		v7.stats.a$GROUP2=""
		for(i in 1:nrow(v6.stats.a)){
		    locus=v6.stats.a[i,"LOCUS"]
	        v6.stats.a[i,"GROUP"]=v7[v7$LIFTLOCUS==locus,"CALLED"]
	        locus=v7.stats.a[i,"LOCUS"]
	        v7.stats.a[i,"GROUP"]=v7[v7$LIFTLOCUS==locus,"CALLED"]
	    }  
	    v6.stats.a[grepl("ver. 6",v6.stats.a$GROUP),"GROUP2"]="shared"
	    v7.stats.a[grepl("ver. 6",v7.stats.a$GROUP),"GROUP2"]="shared"
	    v6.stats.a[!(grepl("ver. 6",v6.stats.a$GROUP)) & grepl("h",v6.stats.a$GROUP),"GROUP2"]="ver. 7 haplot."
	    v7.stats.a[!(grepl("ver. 6",v7.stats.a$GROUP)) & grepl("h",v7.stats.a$GROUP),"GROUP2"]="ver. 7 haplot."
	    v6.stats.a[!(grepl("ver. 6",v6.stats.a$GROUP)) & !(grepl("h",v6.stats.a$GROUP)),"GROUP2"]="ver. 7"
	    v7.stats.a[!(grepl("ver. 6",v7.stats.a$GROUP)) & !(grepl("h",v7.stats.a$GROUP)),"GROUP2"]="ver. 7"
	  
		v6.stats.a$REF="ver. 6"
		v7.stats.a$REF="ver. 7"
		colnames(v7.stats.a)=paste0("v7",colnames(v7.stats.a))
		total=cbind(subset(v6.stats.a,select=c(QUAL,DP,MQ,LOCUS,REF,GROUP,GROUP2,vAD)),subset(v7.stats.a,select=paste0("v7",c("QUAL","DP","MQ","LOCUS","REF","GROUP","GROUP2","vAD"))))
		
		 ggplot(total, aes(x=DP,y=v7DP)) +
		  geom_point(aes(color=GROUP2),alpha=0.2,size=1,stroke=1)+xlim(0,70)+ylim(30,70) + 
		  xlab("Read depth in ver. 6") +ylab ("Read depth in ver. 7") +
		  theme(legend.title=element_blank(),legend.position="",legend.direction="horizontal",
		  legend.key.size = unit(0.06,"line"),
        legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),legend.key=element_blank(),
        legend.box.spacing=unit(0,"line"),
        legend.text=element_text(size=5),panel.background=element_rect(fill = NA), legend.background=element_blank(),
		  panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		  axis.ticks=element_line(size = 0.5),axis.ticks.length=unit(1,"mm"),
		  legend.spacing=unit(0,"mm"),
	     legend.margin=margin(-0.5, unit = "mm"),
		  plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm")) +
          scale_alpha(guide = "none") +
          scale_fill_manual(values = c("#D55E00","#009E73","#CC79A7")) +
        scale_colour_manual(values = c("#D55E00","#009E73","#CC79A7")) 
		ggsave("DP.png",width=75,height=75,units="mm",dpi=1000)
		  
				   
        ggplot(total, aes(x=v7vAD,color=GROUP2,fill=GROUP2)) +
        geom_density(alpha=0.1,size=1) + xlab("Allelic depth of variant allele in ver. 7") + ylab("Frequency") +
        theme(legend.title=element_blank(),legend.position ="none",legend.key.size = unit(0.06,"line"),
        legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),legend.key=element_blank(),
        legend.margin=margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),legend.box.spacing=unit(0,"line"),
        axis.ticks.x=element_line(size=0.5),axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),legend.text=element_text(size=1.5),panel.background=element_rect(fill = NA),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.ticks.length=unit(1,"mm"),plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")) +
        scale_alpha(guide = 'none') +
        scale_fill_manual(values = c("#D55E00","#009E73","#CC79A7")) +
        scale_colour_manual(values = c("#D55E00","#009E73","#CC79A7")) +
        scale_x_continuous(breaks=c(0,0.25,0.50,0.75,1),limits=c(0,1))
        ggsave("AD.png",height=75,width=75,units="mm",dpi=1000)		   
        ggsave("AD.pdf",dev="pdf",height=75,width=75,units="mm")		   

	    array.out=OUT
	    
	    """    
)


}


