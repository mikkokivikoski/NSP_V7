#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 16



object gacu {

  val gasAcuRef = INPUT("Gac-HiC_revised_genome_assembly.fa")
  val gasAcuRefIndex = INPUT("Gac-HiC_revised_genome_assembly.fa.fai")
  val gasAcuRefAgp = INPUT("gacuNEW/contigs.agp")

  val gasAcuRefLA10X = INPUT("result_gacuMappingVariantCalling/prepareReferenceGenomes/folder1/gasAcuLA10X.fasta")
  val gasAcuRefLA10XIndex = INPUT("result_gacuMappingVariantCalling/prepareReferenceGenomes/folder1/gasAcuLA10X.fasta.fai")
  val gasAcuRefLA10XAgp = INPUT("gacuNEW/3SP2.agp")

  val makeMask=INPUT("makeMappabilityMask script")
  val centromereFasta = INPUT("gacuCentromere.fa")
  
  
val fullHaplotypes = INPUT("gacuNEW/fullHaplotypes50.txt")

val gacuReferences = Map("gacuRefIndex" -> gasAcuRefIndex,"gacuRefLA10XIndex" -> gasAcuRefLA10XIndex,
                         "gacuRefAgp" -> gasAcuRefAgp,"gacuRefLA10XAgp" -> gasAcuRefLA10XAgp)

  
  val snpbility=NamedMap[TextFile]("snpbility")
  val repeat=NamedMap[TextFile]("repeat")
  val positiveMask=NamedMap[TextFile]("postiveMask")
  val refs = Map("gacu"-> gasAcuRef,"gasAcuLA10X"->gasAcuRefLA10X)




val summaryStatistics=BashEvaluate(
      array1=gacuReferences,
      var1=fullHaplotypes,
      var2=gasAcuRefAgp,
      command="""
        keys=( $( getarraykeys array1 ) )
        files=( $( getarrayfiles array1 ) )
        for i in $( seq 0 8 );do
         echo $( basename ${files[$i]} )
        done
		
		gacuInd=$( getarraykeyindex array1 gacuRefIndex )
		gacuLA10XInd=$( getarraykeyindex array1 gacuRefLA10XIndex )
		gacuAGPInd=$( getarraykeyindex array1 gacuRefAgp )
		gacuLA10XAGPInd=$( getarraykeyindex array1 gacuRefLA10XAgp )
		
		gaculength=`grep -v "chrUn" ${files[$gacuInd]}|grep -v "chrM" |awk '{sum+=$2}END{print sum}'`
		gaculengthUnGapped=`awk '$6 != "chrM" && $6 != "chrUn"' ${files[$gacuAGPInd]} |awk '{sum+=$3}END{print sum}'`
		gaculengthUn=`awk '$6 == "chrUn"' ${files[$gacuAGPInd]} |awk '{sum+=$3}END{print sum}'`
		
		gacuLA10Xlength=`grep -v "chr.Un" ${files[$gacuLA10XInd]}|grep -v "chrM" |awk '{sum+=$2}END{print sum}'` ##`awk '{sum+=$2}END{print sum}' ${files[$gacuLA10XInd]}`
		gacuLA10XlengthUnGapped=`awk '$5 == "W" {sum+=$8-$7+1}END{print sum}' ${files[$gacuLA10XAGPInd]}`
		gacucontigslg=`awk '$6 != "chrM" && $6 != "chrUn"' ${files[$gacuAGPInd]}|cut -f1 |sort|uniq |wc -l`
		gacucontigsUn=`awk '$6 == "chrUn"' ${files[$gacuAGPInd]}|cut -f1 |sort|uniq |wc -l`
		gacuLA10Xcontigslg=`awk '$5 =="W"'  ${files[$gacuLA10XAGPInd]} |cut -f6|sort|uniq |wc -l`
        gacuLA10XcontigslgUn=`cut -f1 ${files[$gacuLA10XInd]} |grep "chrUn"|wc -l`
        gacuLA10XcontigslgUnLength=`grep "chrUn" ${files[$gacuLA10XInd]}|awk '{sum+=$2}END{print sum}'`
        cut -f2 @var1@|sort|uniq > removedContigs
        awk '($3-$2<0)' gacuNEW/ichrb*.la |cut -f1 >> removedContigs
        removedContigs=`sort removedContigs|uniq|wc -l`
        cat removedContigs > gacuLA10XcontigsRemovedContigs
        cut -f1 gacuNEW/pruned2.txt >> gacuLA10XcontigsRemovedContigs
        gacuLA10XcontigsRemovedContigs=`sort gacuLA10XcontigsRemovedContigs|uniq|wc -l`
        gacuLA10XremovedContigsLength=`grep -wFf gacuLA10XcontigsRemovedContigs @var2@|awk '{sum+=$3} END {print sum}'`
		gacuLA10Xcontigschain=`awk 'BEGIN{OFS="\t";c=1;lg="CHR1"} {if($5 != "W" || $1 != lg) {c+=1;lg=$1} else {lg=$1}} END { print c }' ${files[$gacuLA10XAGPInd]}`
		awk '$6 != "chrM" && $6 != "chrUn"' ${files[$gacuAGPInd]} |cut -f1|sort|uniq > gaculgcontigs
        awk '$5 =="W"'  ${files[$gacuLA10XAGPInd]} |cut -f6|sort|uniq > gacuLA10Xlgcontigs
		gacucontigsnotingacuLA=`grep -vFf gacuLA10Xlgcontigs gaculgcontigs|wc -l`
		gacuLA10Xcontigsnotingacu=`grep -vFf gaculgcontigs gacuLA10Xlgcontigs|wc -l`
		echo -e Ungapped length of the linkage groups \(bp\)'\t'$gaculengthUnGapped'\t'$gacuLA10XlengthUnGapped >> @out1@
		echo -e Total length of the unassigned contigs \(bp\)'\t'$gaculengthUn'\t'$gacuLA10XcontigslgUnLength >> @out1@
		echo -e No. contigs in linkage groups'\t'$gacucontigslg'\t'$gacuLA10Xcontigslg >> @out1@
		echo -e No. contig chains in linkage groups'\t'NOT ASSIGNED'\t'$gacuLA10Xcontigschain >> @out1@
		echo -e Contigs not in linkage groups of the other assembly'\t'$gacucontigsnotingacuLA'\t'$gacuLA10Xcontigsnotingacu >> @out1@
		echo -e No. unassigned contigs'\t'$gacucontigsUn'\t'$gacuLA10XcontigslgUn >> @out1@
		echo -e No. removed contigs'\t'NOT APPLICABLE'\t'$gacuLA10XcontigsRemovedContigs >> @out1@
		echo -e Length of removed contigs'\t'NOT APPLICABLE'\t'$gacuLA10XremovedContigsLength >> @out1@	    
     """
     )
 val sankey=REvaluate(
     inArray=gacuReferences,
     asConnection=true,
     script="""
      table.out=data.frame()
	  OUT=list()
	  rm('optOut1')
      fig.dir <- get.output(cf, 'optOut1')
	  dir.create(fig.dir, recursive=TRUE)
	  setwd(fig.dir)
	  library(ggplot2)
	  library(ggforce)
	  
	  gacu.agp=read.table(array[["gacuRefAgp"]],header=F,stringsAsFactors=F,sep="\t") 
	  gacuNew.agp=read.table(array[["gacuRefLA10XAgp"]],header=F,stringsAsFactors=F,sep="\t")
	  unassigned=read.table("gacuNEW/propb0.la",header=F,sep="\t",stringsAsFactors=F)
	  pruned=read.table("gacuNEW/pruned2.txt",header=F,sep="\t",stringsAsFactors=F)
	  
	  print(head(gacu.agp))
	  print(head(gacuNew.agp))
	  print(head(unassigned))
	  print(head(pruned))
	  gacu.agp=subset(gacu.agp,V1!="chrM.c1")
	  newchr=c()
	  for(i in 1:nrow(gacu.agp)){
	      ctg=as.character(gacu.agp[i,1])
	      if(ctg %in% gacuNew.agp$V6){
	          new=as.character(gacuNew.agp[gacuNew.agp$V6==ctg,1])
		  } else if (ctg %in% unassigned$V1 | ctg %in% pruned$V1){
	          new="Unassigned"
	      } else{
	          new="Rem. haplot."
	      }
	      newchr=c(newchr,new)
      }
      print(length(newchr))
      
      gacu.agp$new=newchr
      gacu.agp$old=sapply(strsplit(gacu.agp$V1,split="\\."),function(x) unlist(x)[1])
      
      
      gacu.agp$group=paste0(gacu.agp$old,":",gacu.agp$new)
      print(head(gacu.agp))
      gacu.agp$length=(gacu.agp$V3-gacu.agp$V2)+1
      sum=aggregate(gacu.agp$length, by=list(gacu.agp$group),sum)
      print(head(sum))
      
      plot.table=data.frame(REF=rep(c("Gacu","Gacu_new"),each=nrow(sum)),ID=sum$Group.1,
      GROUP=c(sapply(strsplit(sum$Group.1,split=":"),function(x) x[1]),sapply(strsplit(sum$Group.1,split=":"),function(x) x[2])),
      LENGTH=sum$x, ID2=1:nrow(sum)) 
      
      plot.table$GROUP=factor(plot.table$GROUP,levels=c(paste0("chr",as.roman(1:21)),"chrUn","chrM",paste0("CHR",1:21),"Unassigned","Rem. haplot."))
      COLORGROUP=rep(as.character(plot.table[plot.table$REF=="Gacu","GROUP"]),times=2)
      COLORGROUP[which(grepl("chrUn:CHR",plot.table$ID))]="chrUn:CHR"
      COLORGROUP[which(grepl("chrUn:Un",plot.table$ID))]="chrUn:un"
      COLORGROUP[which(grepl("chrUn:Rem",plot.table$ID))]="chrUn:rem"
      plot.table$COLORGROUP=factor(COLORGROUP,levels=c(paste0("chr",as.roman(1:21)),"chrUn:CHR","chrUn:un","chrUn:rem"))
    
      sankey.palette1=c(rep(c("pink","pink3"),length.out=21),rep(c("steelblue4","steelblue2"),length.out=9),"purple","grey40","grey20")
 	  library(ggforce)
		ggplot(plot.table, aes(x=REF,id=ID2, split = GROUP, value = LENGTH)) +
		  geom_parallel_sets(aes(fill = COLORGROUP,colour=COLORGROUP), alpha = 0.9, axis.width = 0.1,size=0.1) +
		  geom_parallel_sets_axes(axis.width = 0.1,fill="grey50",colour="grey30",size=0.1) +
		  geom_parallel_sets_labels(size=1.5,colour = 'black',hjust=0,angle=0,position = position_nudge(x = 0.1), data=subset(plot.table,REF=="Gacu_new")) + 
		  geom_parallel_sets_labels(size=1.5,colour = 'black',hjust=1,angle=0,position = position_nudge(x = -0.1), data=subset(plot.table,REF=="Gacu")) + 
		  scale_fill_manual(values=sankey.palette1, aesthetics = "fill", breaks = waiver()) +
		  scale_colour_manual(values=sankey.palette1, aesthetics = "colour", breaks = waiver()) +
		  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_blank(),
		legend.position = "none", axis.text.y = element_blank(), axis.ticks=element_blank(),
		axis.text=element_text(size=8),axis.title=element_text(size=7,face="bold"),plot.margin = unit(c(2,2,2,2), "mm")) +
		labs( x = "Genome version")
		#ggsave("sankey1.png",units="mm",height=80,width=60)    
		ggsave("sankey1.png",units="mm",height=80,width=50,dpi=1000)    
      
     
     """)

val haplotypeContigs = REvaluate(
    var1 = fullHaplotypes,
    script="""
        d=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
        
        d$START=-10
        d$END=-10
        d[grepl("+",d$V10),"START"]=d[grepl("+",d$V10),"V11"]
        d[grepl("+",d$V10),"END"]=d[grepl("+",d$V10),"V12"]

        d[grepl("-",d$V10),"START"]=d[grepl("-",d$V10),"V7"]-d[grepl("-",d$V10),"V12"]
        d[grepl("-",d$V10),"END"]=d[grepl("-",d$V10),"V7"]-d[grepl("-",d$V10),"V11"]
        
        d=subset(d,select=c("V5","START","END"))
        colnames(d)=c("CONTIG","START","END")
        table.out=d
    """
    )


val searchHaplotypes=REvaluate(
var1=gasAcuRefLA10XAgp,
param1="gacuNEW/",
script="""
haplotypes=data.frame()
lgs=paste0("LG",1:21)
gacuNew.agp=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
for(i in 1:21){
    lg=paste0("CHR",i)
    agp=subset(gacuNew.agp,V1==lg & V5=="W")
    lep.anchor=read.table(paste0(param1,"ichrb",i,".la"),header=F,stringsAsFactors=F,sep="\t")
    lep.anchor=subset(lep.anchor, V1 %in% agp$V6)
    
    for(j in 1:nrow(agp)){
        ctg=as.character(agp[j,6])
        #tmp1=subset(lep.anchor, V1==ctg)
        tmp1=agp[j,]
        
        tmp=subset(lep.anchor, V1==ctg)
        if(nrow(tmp)==1){
             ctg.start=as.numeric(tmp[1,2])
             ctg.end=as.numeric(tmp[1,3])
             if(grepl("+",tmp[1,4])){ctg.orientation="+"}
             if(grepl("-",tmp[1,4])){ctg.orientation="-"}
             start=as.numeric(unlist(strsplit(as.character(tmp[1,"V7"]),split="/")))
             end=as.numeric(unlist(strsplit(as.character(tmp[1,"V8"]),split="/")))
           if(ctg.orientation=="+"){  
             if(ctg.start > max(start)){
                 first.trim.contig=tmp[1,9]
                 first.trim.pair=lep.anchor[lep.anchor$V1==first.trim.contig,]
                 if(nrow(first.trim.pair)==1 & first.trim.contig %in% agp$V6){
                         if(grepl("+",first.trim.pair[1,4])){pair.orientation="+"}
                         if(grepl("-",first.trim.pair[1,4])){pair.orientation="-"}
                     if(pair.orientation=="+"){
                       alt.end=as.numeric(unlist(strsplit(first.trim.pair[1,8],split="/")))
                             row=data.frame(CONTIG=first.trim.contig,START=min(alt.end),END=first.trim.pair[1,3])
                                         haplotypes=rbind(haplotypes,row)     
                                        } 
                                 if(pair.orientation=="-"){
                                     alt.end=as.numeric(unlist(strsplit(first.trim.pair[1,7],split="/")))
                             row=data.frame(CONTIG=first.trim.contig,START=first.trim.pair[1,2],END=max(alt.end))
                                         haplotypes=rbind(haplotypes,row)
                                     print(paste("PLUSSTARTMINUS",lg, ctg, ctg.orientation))
                                         
                                 }
                             }
                     }
                     
                     if(min(end) > ctg.end){
                         second.trim.contig=tmp[1,12]
                 second.trim.pair=lep.anchor[lep.anchor$V1==second.trim.contig,]
                 if(nrow(second.trim.pair)==1 & second.trim.contig %in% agp$V6){
                         if(grepl("+",second.trim.pair[1,4])){pair.orientation="+"}
                         if(grepl("-",second.trim.pair[1,4])){pair.orientation="-"}
                     if(pair.orientation=="+"){
                                 alt.end=max(as.numeric(unlist(strsplit(second.trim.pair[1,7],split="/"))))
                                 row=data.frame(CONTIG=second.trim.contig,START=second.trim.pair[1,2],END=alt.end)
                                         haplotypes=rbind(haplotypes,row)
                                         
                                 }
                     if(pair.orientation=="-"){
                                 alt.end=min(as.numeric(unlist(strsplit(second.trim.pair[1,8],split="/"))))
                         row=data.frame(CONTIG=second.trim.contig,START=alt.end,END=second.trim.pair[1,3])
                         haplotypes=rbind(haplotypes,row)
                         print(paste("PLUSENDMINUS",lg, ctg, ctg.orientation))

                                           print(row)
                                 }
                             }
                     
                     }
                 }
         
         if(ctg.orientation=="-"){  
             if(ctg.start > max(start)){
                 first.trim.contig=tmp[1,12]
                 first.trim.pair=lep.anchor[lep.anchor$V1==first.trim.contig,]
                 if(nrow(first.trim.pair)==1 & first.trim.contig %in% agp$V6){
                         if(grepl("+",first.trim.pair[1,4])){pair.orientation="+"}
                         if(grepl("-",first.trim.pair[1,4])){pair.orientation="-"}
                     if(pair.orientation=="+"){
                            alt.end=as.numeric(unlist(strsplit(first.trim.pair[1,7],split="/")))
                            row=data.frame(CONTIG=first.trim.contig,START=first.trim.pair[1,2],END=max(alt.end))
                                        haplotypes=rbind(haplotypes,row)
                                        
                                 }
                     if(pair.orientation=="-"){
                             alt.end=as.numeric(unlist(strsplit(first.trim.pair[1,8],split="/")))
                             row=data.frame(CONTIG=first.trim.contig,START=min(alt.end),END=first.trim.pair[1,3])
                                         haplotypes=rbind(haplotypes,row)
                                         
                                  }
                             }
                     }
                     
                     if(ctg.end<min(end)){
                         second.trim.contig=tmp[1,9]
                 second.trim.pair=lep.anchor[lep.anchor$V1==second.trim.contig,]
                 if(nrow(second.trim.pair)==1 & second.trim.contig %in% agp$V6){
                         if(grepl("+",second.trim.pair[1,4])){pair.orientation="+"}
                         if(grepl("-",second.trim.pair[1,4])){pair.orientation="-"}
                     if(pair.orientation=="+"){
                         alt.end=as.numeric(unlist(strsplit(second.trim.pair[1,8],split="/")))
                                 row=data.frame(CONTIG=second.trim.contig,START=min(alt.end),END=second.trim.pair[1,3])

                                 }
                     if(pair.orientation=="-"){
                                 alt.end=as.numeric(unlist(strsplit(second.trim.pair[1,7],split="/")))
                                 row=data.frame(CONTIG=second.trim.contig,START=second.trim.pair[1,2],END=max(alt.end))
                                         
                                 }
                             }
                     
                     }
                 }
                 
            } 
        
        }
}

haplotypes=subset(haplotypes, START>0)
haplotypes$START=haplotypes$START-1 #BED format
table.out=subset(haplotypes, START>=0)
"""   
)


 for((ref, refName) <- refs){
  snpbility(ref)=BashEvaluate(
      var1=refName,
      param1=ref,
      array1=refs,
      command="""
	        #files=( $( getarrayfiles array1 ) )    
	        #ind=$( getarraykeyindex array1 @param1@ )
	        #ref=${files[$ind]} 
			ref=@var1@
			cd @folder1@
			seqbility-20091110/splitfa $ref 35 |split -l 20000000
			for i in `ls x??`;do
			    echo $i
			    bwa aln -R 1000000 -O 3 -E 3 $ref $i > "$i".sai  
			done
			for i in `ls x??`;do
			    bwa samse $ref "$i".sai $i > "$i".sam
			done
			cat x??.sam| seqbility-20091110/gen_raw_mask.pl > "$ref"_rawMask_35.fa
			seqbility-20091110/gen_mask -l 35 -r 0.5 "$ref"_rawMask_35.fa > "$ref"_mask_35_50.fa			
			cp "$ref"_mask_35_50.fa @out1@
			rm x*
     """
  ).out1
}  
 for((ref, refName) <- refs){
   repeat(ref) = BashEvaluate(
     var1=refName,
     param1=ref,
     array1=refs,
     command="""
     ref=@var1@
     ls "$ref".fai
     echo $ref
     RepeatMasker -xsmall -gff -dir @folder1@ -pa 10 -species "Gasterosteus aculeatus" $ref
     
     awk '!/#/{OFS="\t";print $1,$4,$5}' @folder1@/@param1@.fa.out.gff | sort -k1,1 -k2,2n | bedtools merge -d 100|bedtools slop -i - -g "$ref".fai -b 100|bedtools merge > @out1@
     
     """
   ).out1
}
 for((ref, refName) <- refs){
  positiveMask(ref)=BashEvaluate(
    var1=snpbility(ref),
    var2=repeat(ref),
    var3=makeMask,
    param1=ref,
    command="""
      cd @folder1@
      cp @var1@ tmp_mask_35_50.fa
      @var3@
      zcat tmp_mask_35_50.fa*bed.gz |sort -k1,1 -k2,2n > mask.bed
      rm tmp*
      bedtools subtract -a mask.bed -b @var2@ > @param1@_mask.bed 
      cp @param1@_mask.bed @out1@
    """  
    ).out1
  }

  val masksInFolder=Array2Folder(
  in = positiveMask
  )
  
}
