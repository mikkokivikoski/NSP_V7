#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 16


object mapToLG19inversion {

val prepareReferenceSequences= INPUT("folder containing reference sequences")
val v6ref = INPUT("NSP_V6.fasta")
val v6refIndex = INPUT("NSP_V6.fasta.fai")
val v6agp = INPUT("NSP_V6.agp")
val v7ref = INPUT("V7.fasta")
val v7refIndex = INPUT("V7.fasta.fai")
val v7agp = INPUT("V7.agp")
val lepAnchorOrientation = INPUT("NSPv7genome/final/iall.la")
val individuals = Map("FINPYO33"->"FIN-PYO-33","FINPYO22"->"FIN-PYO-22","FINPYO20"->"FIN-PYO-20","FINPYO21"->"FIN-PYO-21","FINPYO0"->"FIN-PYO-0")

val WGSbwa = NamedMap[BinaryFolder]("bam")
val WGSsort = NamedMap[BinaryFolder]("sortedBams")
val WGSrmdup = NamedMap[BinaryFolder]("rmdup")


for((ind, indName) <- individuals){
	
	WGSbwa(ind) = BashEvaluate(
		    var1=prepareReferenceSequences,
		    var2=v6ref,
		    param1=indName,
		    param2="path/to/fastq/files",
		    command="""
                
		        cd @folder1@ 
				for i in v6 v7; do
	                 bwa mem -t4 -M -R "@RG\tID:@param1@\tSM:@param1@\tPL:illumina\tLB:@param1@\tPU:L0" \
					 @var1@/"$i".fa \
					 @param2@/@param1@"_1.fq.gz" \
					 @param2@/@param1@"_2.fq.gz" \
					 | samtools view -h -b -o "$i"_@param1@_aligned.bam -
			    done
		    """
		    ).folder1

		//STEP2: samtools sort
		WGSsort(ind) = BashEvaluate(
		    var1 = WGSbwa(ind),
		    param1 = indName,
		    command = """
		        cd @folder1@ 
		        for i in v6 v7; do
	              samtools sort -T @folder2@/"$i"_@param1@ -O bam -o "$i"_@param1@_sorted.bam @var1@/"$i"_@param1@_aligned.bam
		        done
		        
		    """
		).folder1     
	    //STEP3: samtools rmdup
	    WGSrmdup(ind) = BashEvaluate(
	        var1 = WGSsort(ind),
	        param1 = indName,
	        command="""
	            cd @folder1@
	            for i in v6 v7; do
	                samtools rmdup @var1@/"$i"_@param1@_sorted.bam "$i"_@param1@_dedup.bam
		            samtools index "$i"_@param1@_dedup.bam "$i"_@param1@_dedup.bai 
	            done
	            
	        """
	    ).folder1
}



val WGSvcf = NamedMap[BinaryFolder]("vcf")
val WGSvcfFiltered = NamedMap[BinaryFolder]("vcfFiltered")

for((ind, indName) <- individuals){
    
    WGSvcf(ind) = BashEvaluate(
		 var1=v6ref,
	     var2=v7ref,
	     var3=prepareReferenceSequences.folder1,
	     var4=WGSrmdup(ind),
	     param1=indName,
	     command="""
	     
	       cd @folder1@
    	   for i in v6 v7; do
			   bcftools mpileup --threads 3 -Ou -f @var3@/"$i".fa @var4@/"$i"_@param1@_dedup.bam | bcftools call -Ou -mv > "$i"_@param1@.vcf
           done
            
		 """
	).folder1

    WGSvcfFiltered(ind) = BashEvaluate(
    var1=prepareReferenceSequences.folder1,
    var2=WGSvcf(ind),
    param1=indName,
    param2="100000",
    command="""
      cd @folder1@
      for i in v6 v7; do
          bcftools filter -s LowQual -e '%QUAL<5' @var2@/"$i"_@param1@.vcf > "$i"_@param1@_filtered.vcf 
          bedtools makewindows -g @var1@/"$i".fa.fai -w @param2@ -i winnum |grep ^LG > windows.bed 
          bcftools view -f .,PASS -v snps -i 'GT~"1/1"' "$i"_@param1@_filtered.vcf | bedtools coverage -a windows.bed -b - > "$i"_@param1@_VarHomozygosity_coverage
          bcftools view -f .,PASS -v snps -i 'GT~"0/1"' "$i"_@param1@_filtered.vcf | bedtools coverage -a windows.bed -b - > "$i"_@param1@_VarHeterozygosity_coverage
           
      done
    """
    ).folder1
    
}

val vcfInFolder=Array2Folder(
  in = WGSvcfFiltered,
  fileMode = ".",
  keyCol = "Key",
  link = false
)

val plotHeterozygosity = REvaluate(
    var1=vcfInFolder.out,
    script="""
	    rm('optOut1')
		table.out=data.frame()
		fig.dir <- get.output(cf, 'optOut1')
		dir.create(fig.dir, recursive=TRUE)
		setwd(fig.dir)
		library(ggplot2)
	    
	    print(list.files(var1))
	    genotypes=c("aa","bb","bb","aa","ab")
  for(lg in paste0("LG",19)){      
        total=data.frame()
        j=1
		for (IND in c("FIN-PYO-20","FIN-PYO-21","FIN-PYO-22","FIN-PYO-33","FIN-PYO-0")){
			id=paste0(unlist(strsplit(IND,split="-")),collapse="")
			GT=genotypes[j]
			for(i in c("v6","v7")){
				HOMOZYGOSITY=read.table(paste0(var1,"/",i,"_",IND,"_VarHomozygosity_coverage"),header=F,stringsAsFactors=F)
				HETEROZYGOSITY=read.table(paste0(var1,"/",i,"_",IND,"_VarHeterozygosity_coverage"),header=F,stringsAsFactors=F)
				HOMOZYGOSITY=subset(HOMOZYGOSITY,V1==lg)
				HETEROZYGOSITY=subset(HETEROZYGOSITY,V1==lg)
				if(i=="v6"){HOMOZYGOSITY$REFERENCE="ver. 6"}
				if(i=="v7"){HOMOZYGOSITY$REFERENCE="ver. 7"}
				HOMOZYGOSITY$GENOTYPE=GT
				HOMOZYGOSITY$VARIANT="1/1"
				HOMOZYGOSITY$INDIVIDUAL=IND
				HOMOZYGOSITY$BIN=HOMOZYGOSITY$V4
				HOMOZYGOSITY$BINSIZE=HOMOZYGOSITY$V7
				HOMOZYGOSITY$NoSites=HOMOZYGOSITY$V8
				HOMOZYGOSITY$POS=HOMOZYGOSITY$BIN*HOMOZYGOSITY$BINSIZE/10^6
				if(i=="v6"){HETEROZYGOSITY$REFERENCE="ver. 6"}
				if(i=="v7"){HETEROZYGOSITY$REFERENCE="ver. 7"}
				HETEROZYGOSITY$GENOTYPE=GT
				HETEROZYGOSITY$VARIANT="0/1"
				HETEROZYGOSITY$INDIVIDUAL=IND
				HETEROZYGOSITY$BIN=HETEROZYGOSITY$V4
				HETEROZYGOSITY$BINSIZE=HETEROZYGOSITY$V7
				HETEROZYGOSITY$NoSites=HETEROZYGOSITY$V8
				HETEROZYGOSITY$POS=HETEROZYGOSITY$BIN*HETEROZYGOSITY$BINSIZE/10^6
	
				total=rbind(total,HETEROZYGOSITY,HOMOZYGOSITY)
		    }
		    j=j+1
		}
		
		if(lg=="LG19"){
		    tmp=subset(total,VARIANT %in% c("0/1","1/1") & REFERENCE %in% c("ver. 6","ver. 7"))
		    tmp$COLOR=paste0(tmp$REFERENCE,tmp$VARIANT)
		    ggplot(tmp, aes(POS, NoSites, colour = COLOR)) + geom_rect(aes(xmin=5.59,xmax=7.3,ymin=0,ymax=0.0065),inherit.aes=F,fill="grey80") +   
	 	    geom_line(aes(color=COLOR),key_glyph = "rect",size=0.5) + facet_grid(GENOTYPE ~ REFERENCE) + labs(y="No. variant sites / bp",x="Genomic position (Mbp)") +
		    scale_y_continuous(breaks=c(0,0.0025,0.005),limits=c(0,0.0065)) +
		    theme(panel.background=element_rect(fill = NA),legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),strip.background=element_blank(),
			strip.text=element_text(size=6),axis.text=element_text(size=6),axis.title=element_text(size=6),
			axis.ticks.length=unit(1.5,"mm"),axis.ticks=element_line(size=0.5),
			plot.margin = margin(t = 1, r = 1, b = 1, l = 1, unit = "mm")) +
            #scale_color_manual(values=c( "#D55E00","#E69F00","#56B4E9","#0072B2"))
            scale_color_manual(values=c( "tan2","tan4","steelblue2","steelblue4"))
			ggsave(paste0(lg,".png"),height=60,width=100,units="mm",dpi=1000)
		}
	}	
		  
    """
)

val getDepth = BashEvaluate(
	var1=WGSrmdup("FINPYO0"),
	command="""
	    cd @folder1@
		for ref in v6 v7; do
		    bam=FINPYO0/folder1/"$ref"_FIN-PYO-0_dedup.bam
			index="$ref".fa.fai
			for i in `cut -f1 $index`; do
		        samtools depth -a $bam -r "$i"| awk 'BEGIN{OFS="\t";c=0;s=0;w=100}{c+=1;s+=$3;if(c==w){print $1,$2,s/w;s=0;c=0}}' >> "$ref"_depth
			done
			sed -i '1i CHR\tBIN\tMEANDEPTH' "$ref"_depth 
		done
		cp v6_depth @out1@
		cp v7_depth @out2@
    """    
   )
   

val plotDepth = REvaluate(
    table1=getDepth.out1,
    table2=getDepth.out2,
    script="""
     rm('optOut1')
	 table.out=data.frame()
	 fig.dir <- get.output(cf, 'optOut1')
	 dir.create(fig.dir, recursive=TRUE)
	 setwd(fig.dir)
	 library(ggplot2)
	 
	 table1$REF="ver. 6"
	 table2$REF="ver. 7"
	 total=rbind(table1,table2)
	 
	 v6=subset(table1,grepl("LG",CHR))$MEANDEPTH
	 v7=subset(table2,grepl("LG",CHR))$MEANDEPTH
	 n <- max(length(v6), length(v7))
     length(v6) <- n                      
     length(v7) <- n
	 total2=data.frame(V6=v6,V7=v7)
	 ggplot(total2,aes(x=x))+ #geom_density(alpha = 0.1)+
	 geom_density( aes(x = V7, y = -..density..), fill="steelblue2",size=0.6) +
	 geom_density( aes(x = V6, y = ..density..), fill= "tan2",size=0.6)  + 
	 xlab("Read depth") + ylab("") + ylim(-0.04,0.04) + scale_x_continuous(breaks=c(0,50,100),limits=c(0,100)) + scale_y_continuous(breaks=NULL)+
	 theme(panel.background=element_rect(NA),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
	 axis.text=element_text(size=8),axis.title=element_text(size=8),plot.margin = unit(c(2,2,1,1), "mm"))
 	 ggsave("depthsLG.png",width=30,height=30,units="mm",dpi=1000)
 	 
    """
)   

   
}	
