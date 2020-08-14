#!/usr/bin/env anduril
import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 16

object Ppun {


//REFERENCE GENOME INFORMATION

   val referenceGenomes = INPUT("folder containing reference sequences")
   
   val v6ref = INPUT("NSP_V6.fasta")
   val v6refIndex = INPUT("NSP_V6.fasta.fai")
   val v6agp = INPUT("NSP_V6.agp")
   val v7ref = INPUT("V7.fasta")
   val v7refIndex = INPUT("V7.fasta.fai")
   val v7agp = INPUT("V7.agp")
   
   val gasAcuRef = INPUT("Gac-HiC_revised_genome_assembly.fa")
   val gasAcuRefIndex = INPUT("Gac-HiC_revised_genome_assembly.fa.fai")
   
   //Linkage maps
   val v6linkageMap = INPUT("NSP_V6_map169.txt")
   val v7linkageMap = INPUT("NSP_V7_map169.txt")
   
//Lep-Anchor files   
   val v7haplotypes = INPUT("cov/haplotype_regions.txt")
   val v7y = INPUT("cov/Y_regions.txt")
   val v7x = INPUT("cov/X_regions.txt")
   val v7autosomal = INPUT("cov/autosome_regions.txt")
   val lepAnchorOrientation = INPUT("final/iall.la")
   val fullHaplotypes = INPUT("final/fullHaplotypes50.txt")
   val laOut = INPUT("final")
//Other   
   val Ppuncentromere =INPUT("centromere.fa")
   val makeMask=INPUT("makeMappabilityMaskScript")
   val repeatLibrary = INPUT("Ninespine_RepeatLibrary.fasta")
   val liftover = INPUT("liftover_fix.awk")

   val individuals = Map("FINPYO33"->"FIN-PYO-33","FINPYO22"->"FIN-PYO-22","FINPYO20"->"FIN-PYO-20","FINPYO21"->"FIN-PYO-21","FINPYO0"->"FIN-PYO-0")
   val fastq = INPUT("folder containing fastq files")

val summaryStatistics=BashEvaluate(
     var1=v6refIndex,
     var2=v7refIndex,
     var3=v6agp,
     var4=v7agp,
     var5=lepAnchorOrientation,
     command="""
                v6length=`awk '{sum+=$2}END{print sum}' @var1@`
                v6lengthlg=`awk '{ if( $1 ~ "LG" ) {sum+=$2} }END{print sum}' @var1@`
                v6lengthlg12=`awk '{ if( $1 == "LG12" ) {sum+=$2} }END{print sum}' @var1@`
                v7length=`awk '{sum+=$2}END{print sum}' @var2@`
                v7lengthlg=`awk '{ if( $1 ~ "LG" ) {sum+=$2} }END{print sum}' @var2@`
                v7lengthlg12=`awk '{ if( $1 == "LG12" ) {sum+=$2} }END{print sum}' @var2@`
                v6contigslg=`cut -f6 @var3@ |grep "pilon"|sort|uniq |wc -l`
                v6contigslg12=`grep ^LG12 @var3@ |cut -f6|grep "pilon"|sort|uniq |wc -l`
                v6contigsunassigned=`grep -v "LG" @var1@ |cut -f1|sort|uniq|wc -l`
                v7contigslg=`cut -f6 @var4@ |grep "pilon"|sort|uniq |wc -l`
                v7contigschain=`awk 'BEGIN{OFS="\t";c=1;lg="LG1"} {if($5 != "W" || $1 != lg) {c+=1;lg=$1} else {lg=$1}} END { print c }' @var4@`
                v7contigslg12=`grep ^LG12 @var4@ |cut -f6|grep "pilon"|sort|uniq |wc -l`
                v7contigsunassigned=`grep -v "LG" @var2@ |cut -f1|sort|uniq|wc -l`
                cut -f6 @var3@ |grep "pilon"|sort|uniq > v6lgcontigs
                cut -f6 @var4@ |grep "pilon"|sort|uniq > v7lgcontigs
                grep "LG12" @var3@|cut -f6 |grep "pilon"|sort|uniq > v6lg12contigs
                grep "LG12" @var4@|cut -f6 |grep "pilon"|sort|uniq > v7lg12contigs
                v6contigsnotinv7=`grep -vFf v7lgcontigs v6lgcontigs|wc -l`
                v7contigsnotinv6=`grep -vFf v6lgcontigs v7lgcontigs|wc -l`
                v6lg12contigsnotinv7=`grep -vFf v7lg12contigs v6lg12contigs|wc -l`
                v7lg12contigsnotinv6=`grep -vFf v6lg12contigs v7lg12contigs|wc -l`
                v7knownorientation=`grep -v ^# @var5@ |awk '$4=="+" || $4=="-" { ++count } END { print count }'` #OF KNOWN ORIENTATION
                v7knownorientationl=`grep -v ^# @var5@ |awk '$4=="+" || $4=="-" { count+=$3-$2+1 } END { print count }'` #OF KNOWN ORIENTATION
                v7knownorientation2=`grep -v ^# @var5@ |awk '$4=="++" || $4=="+++" || $4=="--" || $4=="---" { ++count } END { print count }'` #ORIENTATION KNOWN RESPECT TO NEIGHBOR
                v7knownorientation3=`grep -v ^# @var5@ |awk '$4=="?" { ++count } END { print count }'` #ORIENTATION NO KNOWN RESPECT TO NEIGHBOR
                echo -e '\t'ver. 6'\t'ver. 7 > changetable
                echo -e Total length of the assembly \(bp\)'\t'$v6length'\t'$v7length >> @out1@
                echo -e Reference genome length in linkage groups'\t'$v6lengthlg'\t'$v7lengthlg >> @out1@
                echo -e Reference genome length in LG12'\t'$v6lengthlg12'\t'$v7lengthlg12 >> @out1@
                echo -e No. contigs in linkage groups'\t'$v6contigslg'\t'$v7contigslg >> @out1@
                echo -e No. contig chains in linkage groups'\t'NOT ASSIGNED'\t'$v7contigschain >> @out1@
                echo -e No. contigs in LG12'\t'$v6contigslg12'\t'$v7contigslg12 >> @out1@
                echo -e Contigs not assigned in linkage groups'\t'$v6contigsunassigned'\t'$v7contigsunassigned >> @out1@
                echo -e Contigs not in linkage groups of the other assembly'\t'$v6contigsnotinv7'\t'$v7contigsnotinv6 >> @out1@
                echo -e Contigs not in linkage LG12 of the other assembly'\t'$v6lg12contigsnotinv7'\t'$v7lg12contigsnotinv6 >> @out1@
            echo -e Contigs with known orientation'\t'NOT ASSSIGNED'\t'$v7knownorientation >> @out1@
            echo -e Length of contigs with known orientation'\t'NOT ASSSIGNED'\t'$v7knownorientationl >> @out1@     
     """
)
summaryStatistics._execute="once"

/// CONTIG CLASSIFICATION
val blastCentromere = BashEvaluate(
       var1 = Ppuncentromere,
       var2 = v7ref,
       var3 = v7refIndex,
       command = """
       cd @folder1@
       grep "quiver" @var3@ |cut -f1 > unAssignedContigs 
       samtools faidx @var2@ --region-file unAssignedContigs > unAssignedContigs.fa
       samtools faidx unAssignedContigs.fa
       for i in `cat unAssignedContigs`; do
               samtools faidx unAssignedContigs.fa "$i" > tmpref.fa
               makeblastdb -in tmpref.fa -dbtype nucl
               blastn -task blastn -db tmpref.fa -query @var1@ -outfmt 6 -out tmpOut
               cat tmpOut >> @out1@
       done 
       #filter hits
       echo -e QUERY "\t" SUBJECT "\t" PIDENT "\t"LENGTH "\t" MISMATCH "\t" GAPOPEN "\t" QSTART "\t" QEND "\t" SSTART "\t" SEND "\t" EVALUE "\t" BITSCORE > @out2@ 
       awk '$11 < 0.00001' @out1@ >> @out2@
       rm tmpref.fa*
       rm unAssignedContigs.fa*
       """
      )
     blastCentromere._execute="once" 
   val identifyRepeatContigs = BashEvaluate(
        var1 = repeatMaskerOutput,
        var2 = v7ref,
        var3 = v7refIndex,
        command = """
        cd @folder1@
        grep "quiver" @var3@ |cut -f1 > unAssignedContigs 
        samtools faidx @var2@ --region-file unAssignedContigs > unAssignedContigs.fa
        samtools faidx unAssignedContigs.fa
        for i in `cat unAssignedContigs`; do
           grep "$i" @var1@ |cut -f1,4,5 > tmpRepeats.bed 
           echo -e `grep "$i" @var3@|cut -f1`"\t"0"\t"`grep "$i" @var3@|cut -f2` > tmpContig.bed
           bedtools coverage -a tmpContig.bed -b tmpRepeats.bed >> @out1@
        done
        """
    )
   identifyRepeatContigs._execute="once"
   val identifyZeroCoverageContigs = BashEvaluate(
       param1="cov/split",
       command="""
       for FILE in @param1@/*718*.pb; do
          ctg=$( basename "$FILE" .pb)
          awk -v ctg="$ctg" '{s+=$4}END{print (ctg"\t"s)}' @param1@/"$ctg".female @param1@/"$ctg".male @param1@/"$ctg".pb >> tmp
       done
       awk '$2 == 0' tmp > @out1@    
       rm tmp
       """
   )
   identifyZeroCoverageContigs._execute="once"
    val references = Map("v6agp" -> v6agp, "v6refIndex" -> v6refIndex,"v7reference" -> v7ref, 
   "v7agp" -> v7agp, "v7refIndex" -> v7refIndex, "v7haplotypes" -> v7haplotypes, "v7y" -> v7y,"v7x" -> v7x,
   "v7autosomal" -> v7autosomal)
   val sankeyDiagram = REvaluate(
       var1 = blastCentromere.out2,
       var2 = identifyRepeatContigs.out1,
       var3 = identifyZeroCoverageContigs.out1,
       inArray = references,
       asConnection=true,
       script = """
              table.out = data.frame()
                  rm('optOut1')
                  fig.dir <- get.output(cf, 'optOut1')
                  dir.create(fig.dir, recursive=TRUE)
                  setwd(fig.dir)
                  library(ggforce)
                  
                  v6agp=read.table(array[["v6agp"]],header=F,stringsAsFactors=F)
          v6agp=v6agp[grepl(pattern="quiver",x=v6agp$V6),]
          v6agp$CONTIG = sapply(strsplit(v6agp$V6,split="\\|"), function(x) x[1])
          v6agp$GROUP=v6agp$V1
          v6agp$LENGTH=as.numeric(v6agp$V3)-as.numeric(v6agp$V2)+1

          v6index=read.table(array[["v6refIndex"]],header=F,stringsAsFactors=F)



v6unassigned=subset(v6index,grepl("quiver",V1))
          v6unassigned$CONTIG=sapply(strsplit(v6unassigned$V1,split="\\|"), function(x) x[1])
          v6unassigned$GROUP="Unassigned"
          v6unassigned$LENGTH=v6unassigned$V2

                  v6=data.frame(rbind(v6agp[,c("CONTIG","GROUP","LENGTH")],v6unassigned[,c("CONTIG","GROUP","LENGTH")]))
                  v6$REF="ver. 6"
                  colnames(v6)=paste0(colnames(v6),"v6")
                  
                  v7agp=read.table(array[["v7agp"]],header=F,stringsAsFactors=F)
                  v7agp=v7agp[grepl(pattern="quiver",x=v7agp$V6),]
                  v7agp$CONTIG = sapply(strsplit(v7agp$V6,split="\\|"), function(x) x[1])
                  v7agp$GROUP=v7agp$V1
                  v7agp$LENGTH=as.numeric(v7agp$V3)-as.numeric(v7agp$V2)+1
                  v7unassigned=read.table(array[["v7refIndex"]],header=F,stringsAsFactors=F)
                  v7unassigned$CONTIG=sapply(strsplit(v7unassigned$V1,split="\\|"), function(x) x[1])
                  v7unassigned$GROUP="Unassigned"
                  v7unassigned$LENGTH=v7unassigned$V2
                  v7haplotypes=read.table(array[["v7haplotypes"]],header=F,stringsAsFactors=F)
                  v7haplotypes$CONTIG=v7haplotypes$V1
                  v7haplotypes$GROUP="haplotype"
                  v7ychromosome=read.table(array[["v7y"]],header=F,stringsAsFactors=F)
                  v7ychromosome$CONTIG=v7ychromosome$V1
                  v7ychromosome$GROUP="Y-chromosome"
                  v7xchromosome=read.table(array[["v7x"]],header=F,stringsAsFactors=F)
                  v7xchromosome$CONTIG=v7xchromosome$V1
                  v7xchromosome$GROUP="X-chromosome"
          v7autosomal=read.table(array[["v7autosomal"]],header=F,stringsAsFactors=F)
                  v7autosomal$CONTIG=v7autosomal$V1
                  v7autosomal$GROUP="autosomal"
          
          v7centromeric = read.table(var1,header=T,stringsAsFactors=F)          
          v7centromeric$SUBJECT=sapply(strsplit(v7centromeric$SUBJECT,split="\\|"), function(x) x[1])   
                  v7repeatContigs=read.table(var2,header=F,stringsAsFactors=F)
                  v7repeatContigs=subset(v7repeatContigs, V7>0.2)
                  v7repeatContigs$CONTIG=sapply(strsplit(v7repeatContigs$V1,split="\\|"), function(x) x[1])
                  v7zeroCoverage = read.table(var3,header=F,stringsAsFactors=F)
                  
                  v7unassigned[which(v7unassigned$CONTIG %in% v7zeroCoverage$V1),"GROUP"]="Zero cov."
                  v7unassigned[which(v7unassigned$CONTIG %in% v7repeatContigs$CONTIG),"GROUP"]="Repeat"
                  v7unassigned[which(v7unassigned$CONTIG %in% v7haplotypes$CONTIG),"GROUP"]="Put. haplot."
                  v7unassigned[which(v7unassigned$CONTIG %in% v7ychromosome$CONTIG),"GROUP"]="Put. Y"
                  v7unassigned[which(v7unassigned$CONTIG %in% v7xchromosome$CONTIG),"GROUP"]="Put. X"
                  v7unassigned[which(v7unassigned$CONTIG %in% v7autosomal$CONTIG),"GROUP"]="Put. autos."
                  v7unassigned[which(v7unassigned$CONTIG %in% v7centromeric$SUBJECT),"GROUP"]="Put. centr."
                  
                  
                  v7=data.frame(rbind(v7agp[,c("CONTIG","GROUP","LENGTH")],v7unassigned[,c("CONTIG","GROUP","LENGTH")]))
                  v7$REF="ver. 7"
                  colnames(v7)=paste0(colnames(v7),"v7")
                  out=data.frame()
                 
                 for (i in unique(v6$CONTIG)){
                    source=v6[v6$CONTIGv6==i,c("CONTIGv6","GROUPv6","LENGTHv6","REFv6")]    
                    destination = v7[v7$CONTIGv7==i,c("CONTIGv7","GROUPv7","LENGTHv7","REFv7")]
                    
                    if(nrow(destination)==0){
                        dest=data.frame(CONTIGv7=source$CONTIGv6,GROUPv7="Rem. haplot.",LENGTHv7=source$LENGTHv6,REFv7="ver. 7")
                        tmp.out=cbind(source,dest) 
                    } else if(nrow(source)==1 & nrow(destination)==1 & source$LENGTHv6>destination$LENGTHv7){
                        src=data.frame(CONTIGv6=source$CONTIGv6,GROUPv6=source$GROUPv6,LENGTHv6=c(destination$LENGTHv7,source$LENGTHv6-destination$LENGTHv7),REFv6="ver. 6")
                        dest=data.frame(CONTIGv7=destination$CONTIGv7,GROUPv7=c(destination$GROUPv7,"Rem. overl."),LENGTHv7=c(destination$LENGTHv7,source$LENGTHv6-destination$LENGTHv7),REFv7="ver. 7")
                        tmp.out=cbind(src,dest)
                    } else if(nrow(source)==nrow(destination)){
                        tmp.out=cbind(source,destination) 
                    } else if(nrow(source)==1 & nrow(destination)>1){
                        src=data.frame(CONTIGv6=destination$CONTIGv7,GROUPv6=source$GROUPv6,LENGTHv6=(source$LENGTHv6/nrow(destination)),REFv6=source$REFv6)
                        tmp.out=cbind(src,destination)
                    } else if(nrow(source)>1 & nrow(destination)==1){
                        dest=data.frame(CONTIGv7=source$CONTIGv6,GROUPv7=destination$GROUPv7,LENGTHv7=(destination$LENGTHv7/nrow(source)),REFv7=destination$REFv7)
                        tmp.out=cbind(source,dest)
                    }
                    out=rbind(out,tmp.out)
                
                }
                destination=v7[grepl("canuA",v7$CONTIGv7),c("CONTIGv7","GROUPv7","LENGTHv7","REFv7")]
                src=data.frame(CONTIGv6=destination$CONTIGv7,GROUPv6="LG19 inv.",LENGTHv6=destination$LENGTHv7,REFv6="ver. 6")
                tmp.out=cbind(src,destination)
                out=rbind(out,tmp.out)
  
                out$CLASS=paste0(out$GROUPv6,":",out$GROUPv7)
                
                total=aggregate(out$LENGTHv6, by=list(out$CLASS), sum)
                
                colnames(total)=c("CHANGE","LENGTH")
                total$V6=sapply(strsplit(total$CHANGE,split=":"), function(x) x[1])
                total$V7=sapply(strsplit(total$CHANGE,split=":"), function(x) x[2])
                
                total$V6=factor(total$V6,levels=c(paste0("LG",1:21),"LG19 inv.","Unassigned"))
                total$V7=factor(total$V7,levels=c(paste0("LG",1:21),"Put. Y","Put. X","Put. autos.","Put. haplot.","Put. centr.","Repeat","Zero cov.","Unassigned","Rem. overl.","Rem. haplot."))
                
                
                total$GRANDGROUP="LG"
                total[!(grepl(pattern="LG",x=total$V7)),"GRANDGROUP"]=total[!(grepl(pattern="LG",x=total$V7)),as.character("V7")]
                
                ## Total in ggforce format
                
                tmpv6=data.frame(LENGTH=total$LENGTH,GROUP=total$V6,REF="ver. 6",ID=1:nrow(total),BEFORE=total$V6,AFTER=total$V7)
                tmpv7=data.frame(LENGTH=total$LENGTH,GROUP=total$V7,REF="ver. 7",ID=1:nrow(total),BEFORE=total$V6,AFTER=total$V7)
                total.ggforce=rbind(tmpv6,tmpv7)
                total.ggforce$GROUP=factor(total.ggforce$GROUP, c(paste0("LG",1:21),"Put. Y","Put. X","Put. autos.","Put. haplot.","Put. centr.","Repeat","Zero cov.","LG19 inv.","Unassigned","Rem. overl.","Rem. haplot."))
                total.ggforce$COLORGROUP=""
                
                total.ggforce[which(grepl("LG",total.ggforce$BEFORE)),"COLORGROUP"]=as.character(total.ggforce[which(grepl("LG",total.ggforce$BEFORE)),]$BEFORE)
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & grepl("LG",total.ggforce$AFTER),"COLORGROUP"]="unassignedToLG" 
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Put. Y","COLORGROUP"]="unassignedToY" 
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Put. X","COLORGROUP"]="unassignedToX"
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Put. haplot.","COLORGROUP"]="unassignedTohaplotype"
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Put. centr.","COLORGROUP"]="unassignedToCentromeric"
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Put. autos.","COLORGROUP"]="unassignedToAutosomal"
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Repeat","COLORGROUP"]="unassignedToRepeat"
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Unassigned","COLORGROUP"]="unassignedTounassinged" 
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Rem. haplot.","COLORGROUP"]="unassignedToremoved"
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Rem. overl.","COLORGROUP"]="unassignedTotrimmed" 
                total.ggforce[total.ggforce$BEFORE=="Unassigned" & total.ggforce$AFTER == "Zero cov.","COLORGROUP"]="unassignedToZerocov" 
                total.ggforce$COLORGROUP=factor(total.ggforce$COLORGROUP,levels=c(paste0("LG",1:21),"unassignedToLG","unassignedToY","unassignedToX","unassignedToAutosomal","unassignedTohaplotype","unassignedToCentromeric","unassignedToRepeat","unassignedToZerocov","unassignedTounassinged","LG19 inv.","unassignedTotrimmed","unassignedToremoved"))
        
                  
                sankey.palette1=c(rep(c("pink","pink3"),length.out=21),rep(c("steelblue4","steelblue2"),length.out=9),"purple","grey40","grey20")
                library(ggforce)
                ggplot(total.ggforce, aes(x=REF,id=ID, split = GROUP, value = LENGTH)) +
                  geom_parallel_sets(aes(fill = COLORGROUP,colour=COLORGROUP), alpha = 0.9, axis.width = 0.1,size=0.1) +
                  geom_parallel_sets_axes(axis.width = 0.1,fill="grey50",colour="grey30",size=0.1) +
                  geom_parallel_sets_labels(size=3,colour = 'black',hjust=0,angle=0,position = position_nudge(x = 0.1), data=subset(total.ggforce,REF=="ver. 7")) + 
                  geom_parallel_sets_labels(size=3,colour = 'black',hjust=1,angle=0,position = position_nudge(x = -0.1), data=subset(total.ggforce,REF=="ver. 6")) + 
                  scale_fill_manual(values=sankey.palette1, aesthetics = "fill", breaks = waiver()) +
                  scale_colour_manual(values=sankey.palette1, aesthetics = "colour", breaks = waiver()) +
                  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_blank(),
                legend.position = "none", axis.text.y = element_blank(), axis.ticks=element_blank(),
                axis.text=element_text(size=10),axis.title=element_text(size=7,face="bold"),plot.margin = unit(c(2,2,2,2), "mm")) +
                labs( x = "Genome version")
                ggsave("sankey1.png",units="mm",height=170,width=100,dpi=1000)    
                ggsave("sankey1.pdf",device="pdf",units="mm",height=170,width=100)    
                 
                table.out=out
       
       """
   )
//~ sankeyDiagram._execute="once"



val annotateCentromeres=BashEvaluate(
    var1=v6ref,
    var2=v7ref,
    var3=Ppuncentromere,
    command="""
      cd @folder1@   
          for i in $(seq 1 21 );do
                  samtools faidx @var1@ LG$i > tmprefV6.fa
                  samtools faidx @var2@ LG$i > tmprefV7.fa
                  makeblastdb -in tmprefV6.fa -dbtype nucl
                  makeblastdb -in tmprefV7.fa -dbtype nucl
                  
                  blastn -task blastn -db tmprefV6.fa -query @var3@ -outfmt 6 -out v6Out
                  blastn -task blastn -db tmprefV7.fa -query @var3@ -outfmt 6 -out v7Out
                  
                  cat v6Out >> v6Hits
                  cat v7Out >> v7Hits
      done
    """
    )
    annotateCentromeres._execute="once"
    val linkageMaps = Map("v6" -> v6linkageMap, "v7" -> v7linkageMap, "v6Index" -> v6refIndex, "v7Index" -> v7refIndex)
   val plotCentromeres=REvaluate(
	    var1 = annotateCentromeres.folder1,
	    inArray= linkageMaps,
	    asConnection=true,
	    script="""
	       rm('optOut1')
	       table.out=data.frame()
		   fig.dir <- get.output(cf, 'optOut1')
		   dir.create(fig.dir, recursive=TRUE)
		   setwd(fig.dir)
		   library(ggplot2)
		   set.seed=13

		   CENTROMERES.OUT=data.frame(CHR=paste0("LG",1:21),CENTROMERESTART=0,CENTROMEREEND=0,CENTROMERE=0)
		   blastHits.filtered=data.frame()
           blastHits.all=data.frame()
           marey.total=data.frame()
		   for(v in c("v6","v7")){
		     ref=read.table(array[[paste0(v,"Index")]],header=F,stringsAsFactors=F)
             blastHits=read.table(paste0(var1,"/",v,"Hits"),stringsAsFactors=F,header=F)
		     blastHits$CHRLENGTH=0
             blastHits$REF=v
             blastHits$MIDPOINT=rowMeans(blastHits[,c(9,10)])
             blastHits$loge=log10(blastHits$V11)
             blastHits$CHR=blastHits$V2

		     for(LG in paste0("LG",1:21)){
				 blastHits[blastHits$V2==LG,"CHRLENGTH"]=ref[ref$V1==LG,2]
				 if(LG %in% blastHits$CHR){
					 blast.hits=subset(blastHits,CHR==LG & loge < -5)
					 if(nrow(blast.hits)>0){
					  clustering=kmeans(blast.hits$MIDPOINT,centers=3,iter.max=100)
					  blast.hits$k=clustering$cluster
					  keep=which(prop.table(as.numeric(table(blast.hits$k))) > 0.1)
					  blast.hits=subset(blast.hits,k %in% keep)
					  blastHits.filtered=rbind(blastHits.filtered,blast.hits)  
		             }
				  }
                  if(v=="v7"){
                     CENTROMERES.OUT[CENTROMERES.OUT$CHR==LG,c(2,3,4)]=c(min(c(blast.hits$V9,blast.hits$V10)),max(c(blast.hits$V9,blast.hits$V10)),round(mean(c(min(c(blast.hits$V9,blast.hits$V10)),max(c(blast.hits$V9,blast.hits$V10))))))
                  } 
                   blastHits.all=rbind(blastHits.all,blastHits)
		      }
                   marey=read.table(array[[v]],header=T,skip=2,stringsAsFactors=F,sep="\t",fill=T)
                   marey=subset(marey, CONTIG!="CONTIG",select=c("CONTIG","POS","CHR","MALE_POS","FEMALE_POS"))
                   marey$CHR=marey$CONTIG
                   marey$REF=v

                   marey.total=rbind(marey.total,marey)
}
                   marey.total$POS=as.numeric(marey.total$POS)
                   marey.total$MALE_POS=as.numeric(marey.total$MALE_POS)
                   marey.total$FEMALE_POS=as.numeric(marey.total$FEMALE_POS)
                   marey.total$AVERAGEGENTICDISTANCE=rowMeans(marey.total[,c("MALE_POS","FEMALE_POS")])
                   marey.total$CHR=factor(marey.total$CHR,levels=paste0("LG",1:21))
                   blastHits.all$CHR=factor(blastHits.all$CHR,levels=paste0("LG",1:21))
                   blastHits.filtered$CHR=factor(blastHits.filtered$CHR,levels=paste0("LG",1:21))

           ## RE-NAME REFERENCE GENOMES FOR PLOTTING
           marey.total$REF=gsub(pattern="v6", replacement="ver. 6", x=marey.total$REF)
           marey.total$REF=gsub(pattern="v7", replacement="ver. 7", x=marey.total$REF)
           blastHits.all$REF=gsub(pattern="v6", replacement="ver. 6", x=blastHits.all$REF)
           blastHits.all$REF=gsub(pattern="v7", replacement="ver. 7", x=blastHits.all$REF)
           blastHits.filtered$REF=gsub(pattern="v6", replacement="ver. 6", x=blastHits.filtered$REF)
           blastHits.filtered$REF=gsub(pattern="v7", replacement="ver. 7", x=blastHits.filtered$REF)
           ##
           
		   p=ggplot(marey.total,aes(x=POS,y=AVERAGEGENTICDISTANCE,color=REF)) + geom_point(size=0.5) + facet_grid(CHR ~ REF) + xlab("Genetic position (bp)") + ylab("Genetic distance (cM)")+ scale_y_continuous(breaks=c(0,50,100),limits=c(0,125)) +
		     geom_vline(aes(xintercept = MIDPOINT),col=rgb(0,0,0,0.1), blastHits.all,size=0.25) +
		     scale_colour_manual(values=c("tan2","steelblue2")) +
		      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		      panel.background = element_rect(fill = "white"), strip.background.y=element_blank(),strip.background.x=element_blank(),legend.position = "none",
		      axis.ticks.y=element_line(size=0.5),axis.ticks.length=unit(1,"mm"),axis.text=element_text(size=5),axis.title=element_text(size=5),
		      strip.text.y=element_text(size=5))
		   ggsave("allCentromereHits.png", plot = p,device="png",width=75,height=180,units="mm",dpi=1000)
		   
		   p=ggplot(marey.total,aes(x=POS,y=AVERAGEGENTICDISTANCE,color=REF)) + geom_point(size=0.5) + facet_grid(CHR ~ REF) + xlab("Genetic position (bp)") + ylab("Genetic distance (cM)")+scale_y_continuous(breaks=c(0,50,100),limits=c(0,125)) + 
		      geom_vline(aes(xintercept = MIDPOINT),col=rgb(0,0,0,0.1),blastHits.filtered,show.legend=F,size=0.25) +
		      scale_colour_manual(values=c("tan2","steelblue2")) +
		      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
		      panel.background = element_rect(fill = "white"), strip.background.y=element_blank(),strip.background.x=element_blank(),legend.position = "none",
		      axis.ticks.y=element_line(size=0.5),axis.ticks.length=unit(1,"mm"),axis.text=element_text(size=5),axis.title=element_text(size=5),
		      strip.text.y=element_text(size=5))
		   ggsave("centromeresFiltered.png", plot = p,device="png",width=75,height=180,units="mm",dpi=1000)
		     
		   table.out=CENTROMERES.OUT
		   
	       """
	  )
plotCentromeres._execute="once"


val refArray = Map("gacuRef" -> gasAcuRef, "ppunv6Ref" -> v6ref, "ppunv7ref" -> v7ref)
val minimap2 = BashEvaluate(
    array1=refArray,
    command="""
         keys=( $( getarraykeys array1 ) )
     files=( $( getarrayfiles array1 ) )
     
     gacuInd=$( getarraykeyindex array1 gacuRef )
     v6Ind=$( getarraykeyindex array1 ppunv6Ref )
         v7Ind=$( getarraykeyindex array1 ppunv7Ref )
         
         cd @folder1@
         romans=(I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII XVIII XIX XX XXI XXII)
         threads=3
         v6=NSP_V6.fasta
         v7=V7.fasta
         gasacu=Gac-HiC_revised_genome_assembly.fa
         gasacuIndex=Gac-HiC_revised_genome_assembly.fa.fai
         for i in $( seq 1 21 ); do
            roman=${romans[$i-1]}
            samtools faidx $v6 LG$i > tmp-v6.fa 
            samtools faidx $v7 LG$i > tmp-v7.fa
            samtools faidx $gasacu chr$roman > tmp-gasacu.fa
            if [ $i -eq 12 ]; then
              samtools faidx $gasacu chrVII:1-14000000 chrXII > tmp-gasacu.fa
            fi  
            minimap2 -xasm10 tmp-gasacu.fa tmp-v6.fa > LG"$i"-v6.paf
            minimap2 -xasm10 tmp-gasacu.fa tmp-v7.fa > LG"$i"-v7.paf
           
        done
        rm tmp*
        """
)
minimap2._execute="once"
val processMinimap2 = BashEvaluate(
    var1=minimap2.folder1,
    var2=v6ref,
    var3=v7ref,
    var4=gasAcuRef,
    var5=gasAcuRefIndex,
    command="""
        cd @folder1@
        romans=(I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII XVIII XIX XX XXI XXII)
        threads=3
        v6=@var2@
        v7=@var3@
        gasacu=@var4@
        gasacuIndex=@var5@
        for i in $( seq 1 21 ); do
           roman=${romans[$i-1]}

            awk '$10>=5000' @var1@/LG"$i"-v6.paf > LG"$i"-v6Filtered.paf
            awk '$10>=5000' @var1@/LG"$i"-v7.paf > LG"$i"-v7Filtered.paf
            
            echo -e chr"$roman"\\t0\\t`grep -w ^chr"$roman" $gasacuIndex |cut -f2`\\tLG"$i"-v6\\t50\\t+ > v6ref.bed
            echo -e chr"$roman"\\t0\\t`grep -w ^chr"$roman" $gasacuIndex |cut -f2`\\tLG"$i"-v6\\t50\\t- >> v6ref.bed
            echo -e chr"$roman"\\t0\\t`grep -w ^chr"$roman" $gasacuIndex |cut -f2`\\tLG"$i"-v7\\t50\\t+ > v7ref.bed
            echo -e chr"$roman"\\t0\\t`grep -w ^chr"$roman" $gasacuIndex |cut -f2`\\tLG"$i"-v7\\t50\\t- >> v7ref.bed
            if [ $i -eq 12 ]; then
                   echo -e chrVII\\t0\\t14000000\\tLG"$i"-v6\\t50\\t+ > v6ref.bed
                   echo -e chrVII\\t0\\t14000000\\tLG"$i"-v6\\t50\\t- >> v6ref.bed
                   echo -e chrVII\\t0\\t14000000\\tLG"$i"-v7\\t50\\t+ > v7ref.bed
                   echo -e chrVII\\t0\\t14000000\\tLG"$i"-v7\\t50\\t- >> v7ref.bed
                   echo -e chrXII\\t0\\t`grep -w ^chrXII $gasacuIndex |cut -f2`\\tLG"$i"-v6\\t50\\t+ >> v6ref.bed
               echo -e chrXII\\t0\\t`grep -w ^chrXII $gasacuIndex |cut -f2`\\tLG"$i"-v6\\t50\\t- >> v6ref.bed
               echo -e chrXII\\t0\\t`grep -w ^chrXII $gasacuIndex |cut -f2`\\tLG"$i"-v7\\t50\\t+ >> v7ref.bed
               echo -e chrXII\\t0\\t`grep -w ^chrXII $gasacuIndex |cut -f2`\\tLG"$i"-v7\\t50\\t- >> v7ref.bed
           
           fi
            awk '{ print( $6,$8,$9,$1,50,$5 ) }' OFS='\t' LG"$i"-v6Filtered.paf > tmp-v6alns.bed
            awk '{ print( $6,$8,$9,$1,50,$5 ) }' OFS='\t' LG"$i"-v7Filtered.paf > tmp-v7alns.bed
            bedtools coverage -a v6ref.bed -b tmp-v6alns.bed -s > LG"$i"_v6_syntenyCoverage.bed
            bedtools coverage -a v7ref.bed -b tmp-v7alns.bed -s > LG"$i"_v7_syntenyCoverage.bed

            cat LG"$i"_v6_syntenyCoverage.bed >> @out1@
            cat LG"$i"_v7_syntenyCoverage.bed >> @out1@ 
           
        done
        """
)
processMinimap2._execute="once"
val plotMinimap2Synteny=REvaluate(
    var1=processMinimap2.out1,
    var2=processMinimap2.folder1,
    script="""
     rm('optOut1')
         table.out=data.frame()
         fig.dir <- get.output(cf, 'optOut1')
         dir.create(fig.dir, recursive=TRUE)
         setwd(fig.dir)
         library(ggplot2)
     
      dotplot=function(DATA,QUERY,FIRSTTARGET,SECONDTARGET,FWDCOLOR,REVCOLOR,LWD,IMAGEWIDTH,IMAGEHEIGHT,TWOTARGETS,PLOTNAME){    
             png(paste0(PLOTNAME,".png"),width=IMAGEWIDTH,height=IMAGEHEIGHT,units="mm",res=3000)
                 par(mar = c(3.1, 3.1, 2.1, 2.1)*0.15)
             YLIM=c(0,max(DATA$V2))
                 first.target.length=max(DATA[DATA$V6==FIRSTTARGET,7])
             if(!(TWOTARGETS)){
                XLIM=c(0,first.target.length)
                use.data=subset(DATA,V6 == FIRSTTARGET)
             } else{
                second.target.length=max(DATA[DATA$V6==SECONDTARGET,7])
                XLIM=c(0,first.target.length+second.target.length)
                DATA[DATA$V6==SECONDTARGET,"V8"]=DATA[DATA$V6==SECONDTARGET,"V8"]+first.target.length
                DATA[DATA$V6==SECONDTARGET,"V9"]=DATA[DATA$V6==SECONDTARGET,"V9"]+first.target.length
                use.data=subset(DATA,V6 %in% c(FIRSTTARGET,SECONDTARGET))
             }
             plot('n',xlim=XLIM,ylim=YLIM,ylab=NA,xlab=NA,xaxt='n',cex.lab=1,cex.axis=0.1,yaxt='n',lwd=0.1,bty="n")#cex.lab=3,cex.axis=2
             mtext(QUERY,side=2,line=0,outer=FALSE,adj=0.5,cex=0.2)
             for(i in 1:nrow(use.data)){
               strand=use.data[i,5]
               if(strand=="-"){
                       x0=use.data[i,9]
                       x1=use.data[i,8]
                       y0=use.data[i,3]
                       y1=use.data[i,4]
                       COLOR=REVCOLOR
                       segments(x0=x0,y0=y0,x1=x1,y1=y1,col=COLOR,lwd=LWD)
               } else{
                 x0=use.data[i,8]
                 x1=use.data[i,9]
                 y0=use.data[i,3]
                 y1=use.data[i,4]
                 COLOR=FWDCOLOR
                 segments(x0=x0,y0=y0,x1=x1,y1=y1,col=COLOR,lwd=LWD)
               }
            }
            mtext(unlist(strsplit(FIRSTTARGET,split=":"))[1], side = 1, line = -0.5, outer = FALSE, at = 0,adj=0, cex=0.2)
            if(TWOTARGETS){
                abline(v=first.target.length,lwd=LWD)
                mtext(SECONDTARGET, side = 1, line = -0.5, outer = FALSE, at = first.target.length+5*10^5,adj=0,cex=0.2)
            }
           dev.off()
        }
        ##EXAMPLE-RUN:dotplot(d2,"LG12","chrVII:1-14000000","chrXII","red","blue",4,300,300,TRUE,"LG12-v7")
     
     
     
     d=read.table(var1,header=F,stringsAsFactors=F)
     colnames(d)=c("CHR","START","END","QUERY","OMIT","STRAND","HITS","HITLENGTH","REFLENGTH","COVERAGE")
     d$GENOME=sapply(strsplit(d$QUERY,split="-"), function(x) x[2])
     d$PPUNCHR=sapply(strsplit(d$QUERY,split="-"), function(x) x[1])
     
     strand.change.data=c()
     for(i in 1:21){
         d=read.table(paste0(var2,"/LG",i,"-v6Filtered.paf"),header=F,stringsAsFactors=F,fill=T)
         d=d[order(d$V3),]
         strand=d[1,5]
         strand.changes=0
         for(r in 1:nrow(d)){
            tmp.strand=d[r,5]
            if(tmp.strand != strand){
                        strand.changes=strand.changes+1
                        strand=tmp.strand
                    }
             }
             strand.change.data=c(strand.change.data,strand.changes)
             
             d=read.table(paste0(var2,"/LG",i,"-v7Filtered.paf"),header=F,stringsAsFactors=F,fill=T)
         d=d[order(d$V3),]
         strand=d[1,5]
         strand.changes=0
         for(r in 1:nrow(d)){
            tmp.strand=d[r,5]
            if(tmp.strand != strand){
                        strand.changes=strand.changes+1
                        strand=tmp.strand
                    }
             }
             strand.change.data=c(strand.change.data,strand.changes)
     }
     plot.strand.change=data.frame(PPUNCHR=rep(paste0("LG",1:21),each=2),NO.CHANGES=strand.change.data,GENOME=rep(c("ver. 6","ver. 7")))
     plot.strand.change$LG=factor(plot.strand.change$PPUNCHR,levels=paste0("LG",1:21))
     p=ggplot(plot.strand.change,aes(x=LG,y=NO.CHANGES))+geom_bar(stat="identity",position=position_dodge(),aes(fill=GENOME)) + 
     ylab("No. orientation changes") +xlab("Linkage group") +  
     scale_fill_manual(values=c("tan2","steelblue2")) +
     theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_blank(), axis.ticks=element_blank(),
                legend.title = element_text(size = 3), legend.text = element_text(size = 3),
                legend.key.size = unit(0.5,"line"),
                axis.text=element_text(size=3),axis.title=element_text(size=7))
     ggsave("orientationChanges.png",p,width=100,height=50,units="mm")
     
     
     v6=read.table(paste0(var2,"/LG12-v6Filtered.paf"),header=F,stringsAsFactors=F,fill=T)
     v7=read.table(paste0(var2,"/LG12-v7Filtered.paf"),header=F,stringsAsFactors=F,fill=T)
     dotplot(v6,"LG12","chrVII:1-14000000","chrXII","red","blue",0.15,20,20,TRUE,"LG12-v6")
     dotplot(v7,"LG12","chrVII:1-14000000","chrXII","red","blue",0.15,20,20,TRUE,"LG12-v7")
       
    """
)
plotMinimap2Synteny._execute="once"

//COMPARE LG12 SEQUENCNG DEPTH
   val sequencingdepth=BashEvaluate(
		var1=referenceGenomes,
		var2=bam,
		command="""
		cd @folder1@
		for i in 20 21 22 33 0; do
			samtools depth -a @var2@/rmdup_FINPYO"$i"/folder1/v6_FIN-PYO-"$i"_dedup.bam -r LG12 > v6_FIN-PYO-"$i".depth
			samtools depth -a @var2@/rmdup_FINPYO"$i"/folder1/v7_FIN-PYO-"$i"_dedup.bam -r LG12 > v7_FIN-PYO-"$i".depth
			samtools depth -a @var2@/rmdup_FINPYO"$i"/folder1/v6_FIN-PYO-"$i"_dedup.bam -r LG1 |awk 'BEGIN {sum=0;c=0} {sum+=$3;c+=1} END {print sum/c}' > v6_FIN-PYO-"$i".autosomal.depth
			samtools depth -a @var2@/rmdup_FINPYO"$i"/folder1/v7_FIN-PYO-"$i"_dedup.bam -r LG1 |awk 'BEGIN {sum=0;c=0} {sum+=$3;c+=1} END {print sum/c}' > v7_FIN-PYO-"$i".autosomal.depth
		done
		"""
		)
	sequencingdepth._execute="once"
		
	val getDepth=REvaluate(
  var1 = sequencingdepth.folder1,
  script="""
	 rm('optOut1')
			 table.out=data.frame()
			 fig.dir <- get.output(cf, 'optOut1')
			 dir.create(fig.dir, recursive=TRUE)
			 setwd(fig.dir)
			 library(ggplot2)
		 total=data.frame()
		 bin.size=1000
		 for(i in c(0,20,21,22,33)){
		for(r in c("v6","v7")){
							ind=paste0("FIN-PYO-",i)
				tmp=read.table(paste0(var1,"/",r,"_",ind,".depth"),header=F,stringsAsFactors=F,sep="\t")
				norm=read.table(paste0(var1,"/",r,"_",ind,".autosomal.depth"),header=F,stringsAsFactors=F,sep="\t")   
				tmp$BIN=ceiling(tmp$V2/bin.size)
				tmp2=aggregate(tmp$V3,by=list(tmp$BIN),mean)
				colnames(tmp2)=c("BIN","DEPTH")
				tmp2$NORM=tmp2$DEPTH/as.numeric(norm[1,1])
				tmp2$REF=r
				tmp2$IND=ind
				tmp2$SEX="FEMALE"
				if(ind=="FIN-PYO-0"){tmp2$SEX="MALE"}
				total=rbind(total,tmp2)
					}
 }
     table.out=total 
    """    
) 
getDepth._execute="once"
    val plotDepth=REvaluate(
      table1 = getDepth.table,
      table2 = genotypeStats,
      script="""
         rm('optOut1')
		 table.out=data.frame()
		 fig.dir <- get.output(cf, 'optOut1')
		 dir.create(fig.dir, recursive=TRUE)
		 setwd(fig.dir)
		 library(ggplot2)
		 table1$REF=gsub(pattern="v6",replacement="ver. 6",x=table1$REF)
		 table1$REF=gsub(pattern="v7",replacement="ver. 7",x=table1$REF)
         table1$SEX="FEMALE"
         table1[table1$IND=="FIN-PYO-0","SEX"]="MALE"
         table1$BIN=1000*table1$BIN
          ggplot(table1, aes(y=NORM,x=BIN)) + geom_point(aes(color=SEX),alpha=0.1,size=0.05,stroke=0.05) +  
         facet_grid(IND~REF) + scale_color_manual(values=c("#009E73", "#CC79A7")) + ylab("Normalised read depth")+
         scale_y_continuous(breaks=c(0,1,2),limits=c(0,2)) + xlab("Genomic position") +
         theme(legend.position="none",panel.background=element_rect(fill="white"),panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),strip.background=element_rect(fill="white"),
         axis.title=element_text(size=2), axis.text=element_text(size=2),strip.text=element_text(size=1.8),
         axis.ticks=element_line(size=0.15),axis.ticks.length=unit(0.5,"mm"),
         plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"))
         ggsave("lg12depth.png",width=40,height=40,units="mm",dpi=1000)         
         
         ggplot(subset(table1,IND %in% paste0("FIN-PYO-",c(0,20))), aes(y=NORM,x=BIN)) + geom_point(aes(color=REF),alpha=0.1,size=0.1,stroke=0.1) +  
         facet_grid(IND~REF) + scale_color_manual(values=c("tan2", "steelblue2")) + ylab("Normalised read depth") +
         scale_y_continuous(breaks=c(0,1,2),limits=c(0,2)) + xlab("Genomic position") +
         theme(legend.position="none",panel.background=element_rect(fill="white"),panel.grid.major=element_blank(),
         panel.grid.minor=element_blank(),strip.background=element_rect(fill="white"),
         axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title.x=element_blank(),
         axis.title=element_text(size=4), axis.text=element_text(size=4),strip.text=element_text(size=4),
         axis.ticks=element_line(size=0.5),axis.ticks.length=unit(0.5,"mm"),
         plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"))
         ggsave("lg12depthb.png",width=100,height=25,units="mm",dpi=1000)               
          
    """
    )
 plotDepth._execute="once"
 

 //EXECUTE JOINT CALLING AND MOSAICISM ANALYSIS.

val v6haplotypecaller = NamedMap[BinaryFolder]("v6haplotypecaller")
val v7haplotypecaller = NamedMap[BinaryFolder]("v7haplotypecaller")

for((ind, indName) <- individuals){
	    v6haplotypecaller(ind) = BashEvaluate(
	     var3=v6ref,
	     var4=v7ref,
	     param1=indName,
	     param2=ind,
	     command="""
	       cd @folder1@
    	   	gatk4 HaplotypeCaller \
    	   	    --intervals LG1 \
    	   	    --intervals LG2 \
    	   	    --intervals LG3 \
    	   	    --intervals LG4 \
    	   	    --intervals LG5 \
    	   	    --intervals LG6 \
    	   	    --intervals LG7 \
    	   	    --intervals LG8 \
    	   	    --intervals LG9 \
    	   	    --intervals LG10 \
    	   	    --intervals LG11 \
    	   	    --intervals LG12 \
    	   	    --intervals LG13 \
    	   	    --intervals LG14 \
    	   	    --intervals LG15 \
    	   	    --intervals LG16 \
    	   	    --intervals LG17 \
    	   	    --intervals LG18 \
    	   	    --intervals LG19 \
    	   	    --intervals LG20 \
    	   	    --intervals LG21 \
	            -R @var3@ \
	            -I <path/to/ver6/ind.bam> \
	            -O v6_@param2@.g.vcf.gz \
	            -ERC GVCF \
	            --native-pair-hmm-threads 1 \
	            -ploidy 2
		 """
	    ).folder1
}

for((ind, indName) <- individuals){
	    v7haplotypecaller(ind) = BashEvaluate(
	     var3=v6ref,
	     var4=v7ref,
	     param1=indName,
	     param2=ind,
	     command="""
	       cd @folder1@
    	   	gatk4 HaplotypeCaller \
	            --intervals LG1 \
    	   	    --intervals LG2 \
    	   	    --intervals LG3 \
    	   	    --intervals LG4 \
    	   	    --intervals LG5 \
    	   	    --intervals LG6 \
    	   	    --intervals LG7 \
    	   	    --intervals LG8 \
    	   	    --intervals LG9 \
    	   	    --intervals LG10 \
    	   	    --intervals LG11 \
    	   	    --intervals LG12 \
    	   	    --intervals LG13 \
    	   	    --intervals LG14 \
    	   	    --intervals LG15 \
    	   	    --intervals LG16 \
    	   	    --intervals LG17 \
    	   	    --intervals LG18 \
    	   	    --intervals LG19 \
    	   	    --intervals LG20 \
    	   	    --intervals LG21 \
    	   	    -R @var4@ \
	            -I <path/to/ver7/ind.bam> \
	            -O v7_@param2@.g.vcf.gz \
	            -ERC GVCF \
	            --native-pair-hmm-threads 1 \
	            -ploidy 2
		 """
	    ).folder1
}



val v6combinegvcf=BashEvaluate(
    var1=v6haplotypecaller("FINPYO20"),
    var2=v6haplotypecaller("FINPYO21"),
    var3=v6haplotypecaller("FINPYO22"),
    var4=v6haplotypecaller("FINPYO33"),
    var5=v6haplotypecaller("FINPYO0"),
    var6=v6ref,
    command="""
     gatk4 CombineGVCFs \
       -R @var6@ \
	   --variant @var1@/v6_FINPYO20.g.vcf.gz \
	   --variant @var2@/v6_FINPYO21.g.vcf.gz \
	   --variant @var3@/v6_FINPYO22.g.vcf.gz \
	   --variant @var4@/v6_FINPYO33.g.vcf.gz \
	   --variant @var5@/v6_FINPYO0.g.vcf.gz \
	   -O @folder1@/v6.gvcf.gz
	
    """    
)
val v7combinegvcf=BashEvaluate(
    var1=v7haplotypecaller("FINPYO20"),
    var2=v7haplotypecaller("FINPYO21"),
    var3=v7haplotypecaller("FINPYO22"),
    var4=v7haplotypecaller("FINPYO33"),
    var5=v7haplotypecaller("FINPYO0"),
    var6=v7ref,
    command="""
     gatk4 CombineGVCFs \
	   -R @var6@ \
	   --variant @var1@/v7_FINPYO20.g.vcf.gz \
	   --variant @var2@/v7_FINPYO21.g.vcf.gz \
	   --variant @var3@/v7_FINPYO22.g.vcf.gz \
	   --variant @var4@/v7_FINPYO33.g.vcf.gz \
	   --variant @var5@/v7_FINPYO0.g.vcf.gz \
	   -O @folder1@/v7.gvcf.gz
	
    """    
)

val genotypegvcf = BashEvaluate(
   var1=v6combinegvcf.folder1,
   var2=v6ref,
   var3=v7combinegvcf.folder1,
   var4=v7ref,
	command="""
	gatk4 --java-options "-Xmx10G" GenotypeGVCFs \
	    -R @var2@ \
	    -V @var1@/v6.gvcf.gz \
	    -O @folder1@/v6.vcf.gz
	gatk4 --java-options "-Xmx10G" GenotypeGVCFs \
	    -R @var4@ \
	    -V @var3@/v7.gvcf.gz \
	    -O @folder1@/v7.vcf.gz
	"""
	
)

val getStats = BashEvaluate(
    var1 = genotypegvcf.folder1,
    var2 = v6ref,
    var3 = v7ref,
    param2="10000",
    command="""
    bcftools view -M2 -v snps @var1@/v6.vcf.gz |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |grep ^LG > @out1@    
    bcftools view -M2 -v snps @var1@/v7.vcf.gz |bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n' |grep ^LG > @out2@    
    
    """
)
getStats._execute="once"
val lg12HMM = REvaluate(
    var1=getStats.out1,
    var2=getStats.out2,
    table1=plotCentromeres.table,
    script="""
      rm('optOut1')
	 table.out=data.frame()
	 fig.dir <- get.output(cf, 'optOut1')
	 dir.create(fig.dir, recursive=TRUE)
	 setwd(fig.dir)
	 library(ggplot2)
	 library(ggforce)
	 library(HMM)
	 v6=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
	 v7=read.table(var2,header=F,stringsAsFactors=F,sep="\t")
	 colnames(v6)=c("CHROM","POS","REF","ALT",paste0("FIN-PYO-",c(0,20,21,22,33)))
	 colnames(v7)=c("CHROM","POS","REF","ALT",paste0("FIN-PYO-",c(0,20,21,22,33)))
	 BIN.SIZE=50000
	 v6$BIN=floor(v6$POS/BIN.SIZE)
	 v7$BIN=floor(v7$POS/BIN.SIZE)
	 stat2=function(x){
        s2=round(100*(x[2]/sum(x[2],x[4],1)))
        return(s2)
     }
         
         total2=data.frame()
         out=data.frame()
    	 tmp=subset(v7,CHROM=="LG12")
         for(ind in paste0("FIN-PYO-",c(0,20,21,22,33))){
	         tmp2=data.frame(table(tmp$BIN,tmp[,ind]))
	         tmp3=data.frame(CHROM="LG12",BIN=as.numeric(tmp2$Var1),GT=tmp2$Var2,COUNT=as.numeric(tmp2$Freq),IND=ind,REF="ver. 7")
	         tmp2stat2=as.numeric(apply(table(tmp[,ind],tmp$BIN),2,function(x) stat2(x)))
	         tmp2sum=as.numeric(apply(table(tmp[,ind],tmp$BIN),2,sum))
	         total2=rbind(total2,data.frame(CHROM="LG12",BIN=unique(as.numeric(tmp2$Var1)),STAT2=tmp2stat2,NLOCI=tmp2sum,IND=ind,REF="ver. 7"))
	     }
	    lg12total=data.frame()
		centromere=as.numeric(table1[table1$CHR=="LG12","CENTROMERESTART"])
		for(ind in paste0("FIN-PYO-",c(0,20,21,22,33))){
			 tmp=subset(total2,IND==ind & REF=="ver. 7" & CHROM=="LG12" & BIN < ceiling(centromere/BIN.SIZE)) #choose only the sex chromosome part before centromere
			 obs=tmp$STAT
			 SYMBOLS=seq(1,40)
		     xySYMBOLS=seq(-20,20)
		     xySYMBOLS=seq(0,100)
			 obs=tmp$STAT2
		     x.probs=c(rep(0.001,21),1:20)
		     x.probs=x.probs/sum(x.probs)
		     y.probs=rev(x.probs)
		     xy.mosaic.probs=as.numeric(prop.table(table(factor(round(rnorm(mean=0,sd=2,n=10000)),levels=-20:20))))+0.000001
		     xy.mosaic.probs=xy.mosaic.probs/sum(xy.mosaic.probs)
		     
		     x=round(100*rnorm(1000000,mean=1,sd=0.0625))
		     x.probs=as.numeric(prop.table(table(factor(x[x>=0 & x<=100],levels=0:100))))
		     y=round(100*rnorm(1000000,mean=0,sd=0.0625))
		     y.probs=as.numeric(prop.table(table(factor(y[y>=0 & y<=100],levels=0:100))))
		     m=round(100*rnorm(1000000,mean=0.5,sd=0.25))
		     xy.mosaic.probs=as.numeric(prop.table(table(factor(m[m>=0 & m<=100],levels=0:100))))

		     modellist=list()
		     modellist[[1]]=list(c(0.33,0.34,0.33),matrix(c(0.34,0.33,0.33,0.33,0.34,0.33,0.33,0.33,0.34),3,byrow=T),matrix(c(y.probs,xy.mosaic.probs,x.probs),3,byrow=T))
		     modellist[[2]]=list(c(0.33,0.34,0.33),matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),3,byrow=T),matrix(c(y.probs,xy.mosaic.probs,x.probs),3,byrow=T))
		     x=round(100*rnorm(1000000,mean=1,sd=0.25))
		     x.probs=as.numeric(prop.table(table(factor(x[x>=0 & x<=100],levels=0:100))))
		     y=round(100*rnorm(1000000,mean=0,sd=0.25))
		     y.probs=as.numeric(prop.table(table(factor(y[y>=0 & y<=100],levels=0:100))))
		     m=round(100*rnorm(1000000,mean=0.5,sd=0.25))
		     xy.mosaic.probs=as.numeric(prop.table(table(factor(m[m>=0 & m<=100],levels=0:100))))
		     modellist[[3]]=list(c(0.33,0.34,0.33),matrix(c(0.34,0.33,0.33,0.33,0.34,0.33,0.33,0.33,0.34),3,byrow=T),matrix(c(y.probs,xy.mosaic.probs,x.probs),3,byrow=T))
		     modellist[[4]]=list(c(0.33,0.34,0.33),matrix(c(0.8,0.1,0.1,0.1,0.8,0.1,0.1,0.1,0.8),3,byrow=T),matrix(c(y.probs,xy.mosaic.probs,x.probs),3,byrow=T))
		     
		     for(i in 1:length(modellist)){
			     xyhmm = initHMM(States=c("Y","MOSAIC","X"), Symbols=as.character(xySYMBOLS), startProbs=modellist[[i]][[1]], transProbs=modellist[[i]][[2]], emissionProbs=modellist[[i]][[3]])
	             bw = baumWelch(hmm=xyhmm, observation=as.character(obs), maxIterations=1000, delta=1E-9, pseudoCount=1)
			     vit = viterbi(hmm=bw$hmm, observation=as.character(obs))
			     values=rep(100,length(obs))
			     values[which(vit=="Y")]=-1
			     values[which(vit=="X")]=1
			     values[which(vit=="MOSAIC")]=0
			     lg12total=rbind(lg12total,data.frame(CHROM="LG12",REF="ver. 7",IND=ind,BIN=seq(1,length(obs)),
			     STATE=vit,STAT=obs,NLOCI=tmp$NLOCI,VALUE=values,ITERATION=i))
   			     out=rbind(out,data.frame(IND=ind,ITER=i,X=as.numeric(prop.table(table(vit))["X"])))

	        }
	     }
	 
	 lg12total$GENOMICPOSITION=lg12total$BIN*BIN.SIZE
	 lg12total=subset(lg12total,IND!="FIN-PYO-0")

	 lg12total$STATE=factor(lg12total$STATE,levels=c("MOSAIC","Y","X"))
     
	 ggplot(lg12total,aes(x=GENOMICPOSITION)) + geom_rect(aes(ymin=0,ymax=1,xmin=GENOMICPOSITION,xmax=GENOMICPOSITION+BIN.SIZE,fill=STATE)) + 
	 scale_fill_manual(values=rev(c("steelblue4","steelblue1","#CC79A7"))) + facet_grid(CHROM+REF+ITERATION~IND) + xlab("Genomic position") + ylab("Iteration")+
	 scale_x_continuous(breaks=c(0,15*10^6)) +
	 theme(axis.line.y=element_blank(),axis.text.y=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
	 panel.grid.minor = element_blank(), axis.ticks.y=element_blank(),axis.text=element_text(size=5),
	 strip.background.x = element_blank(),strip.text.y = element_blank(),strip.text.x = element_text(size=5),
	 legend.title = element_text(size = 5),legend.text = element_text(size = 5),legend.key.size = unit(0.1,"line")) 
	 ggsave("LG12HMM1.png",width=150,height=75,units="mm",dpi=1000)


	 ggplot(subset(lg12total),aes(x=1,fill=STATE)) + geom_bar(position="fill") + 
	 scale_fill_manual(values=rev(c("steelblue4","steelblue1","#CC79A7"))) + facet_grid(CHROM+REF+ITERATION~IND) + xlab("Iteration") + ylab("Proportion") +coord_flip()+
	 theme(axis.line.y=element_blank(),axis.text.y=element_blank(), panel.border = element_blank(), panel.grid.major = element_blank(),
	 panel.grid.minor = element_blank(),axis.text=element_text(size=5),
	 strip.background.x = element_blank(),strip.text.y = element_blank(),strip.text.x = element_text(size=5),
	 axis.ticks.y=element_blank())
	 ggsave("LG12HMM2.png",width=150,height=75,units="mm",dpi=1000)
     
     table.out=out

    """
)
lg12HMM._execute="once"
val plotMosaicism = REvaluate(
    var1=getStats.out1,
    var2=getStats.out2,
    script="""
     rm('optOut1')
	 table.out=data.frame()
	 fig.dir <- get.output(cf, 'optOut1')
	 dir.create(fig.dir, recursive=TRUE)
	 setwd(fig.dir)
	 library(ggplot2)
	 library(ggforce)
	 library(HMM)
	 OUT=list()
	 v6=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
	 v7=read.table(var2,header=F,stringsAsFactors=F,sep="\t")
	 colnames(v6)=c("CHROM","POS","REF","ALT",paste0("FIN-PYO-",c(0,20,21,22,33)))
	 colnames(v7)=c("CHROM","POS","REF","ALT",paste0("FIN-PYO-",c(0,20,21,22,33)))
	 BIN.SIZE=50000
	 v6$BIN=floor(v6$POS/BIN.SIZE)
	 v7$BIN=floor(v7$POS/BIN.SIZE)
	 total=data.frame()
     total2=data.frame()
     stat=function(x){
        s3=ceiling(-10*log10(((x[2]-x[4])^2+1)/(sum(x[2],x[4])^2+1)))
        if(s3==0){s3=1} # CONVERT 0 to 1 
        if(s3>40){s3=40} # CONVERT values above 40 to 40
        return(s3)
     }
     
     for(lg in paste0("LG",1:21)){
         tmp=subset(v6,CHROM==lg)
         for(ind in paste0("FIN-PYO-",c(0,20,21,22,33))){
	         tmp2=data.frame(table(tmp$BIN,tmp[,ind]))
	         tmp3=data.frame(CHROM=lg,BIN=as.numeric(tmp2$Var1),GT=tmp2$Var2,COUNT=as.numeric(tmp2$Freq),IND=ind,REF="ver. 6")
	         total=rbind(total,tmp3)
	         tmp2stat=as.numeric(apply(table(tmp[,ind],tmp$BIN),2,function(x) stat(x)))
	         tmp2sum=as.numeric(apply(table(tmp[,ind],tmp$BIN),2,sum))
	         total2=rbind(total2,data.frame(CHROM=lg,BIN=unique(as.numeric(tmp2$Var1)),STAT=tmp2stat,NLOCI=tmp2sum,IND=ind,REF="ver. 6"))
	     }
	 }
	 
     for(lg in paste0("LG",1:21)){
         tmp=subset(v7,CHROM==lg)
         for(ind in paste0("FIN-PYO-",c(0,20,21,22,33))){
	         tmp2=data.frame(table(tmp$BIN,tmp[,ind]))
	         tmp3=data.frame(CHROM=lg,BIN=as.numeric(tmp2$Var1),GT=tmp2$Var2,COUNT=as.numeric(tmp2$Freq),IND=ind,REF="ver. 7")
	         total=rbind(total,tmp3)
	         tmp2stat=as.numeric(apply(table(tmp[,ind],tmp$BIN),2,function(x) stat(x)))
	         tmp2sum=as.numeric(apply(table(tmp[,ind],tmp$BIN),2,sum))
	         total2=rbind(total2,data.frame(CHROM=lg,BIN=unique(as.numeric(tmp2$Var1)),STAT=tmp2stat,NLOCI=tmp2sum,IND=ind,REF="ver. 7"))
	     }
	 }
	 
	 table.out=total2

	 lg12tmp=subset(total,CHROM=="LG12")
	ggplot(lg12tmp,aes(x=BIN*BIN.SIZE,y=COUNT)) + geom_line(aes(color=GT),size=0.05) + facet_grid(IND~REF) +
	 xlab("GENOMIC POSITION") + ylab("Number of sites") +
	 theme(strip.text.x = element_text(size = 2),
	 strip.text.y = element_text(size = 1),
	 axis.text=element_text(size=1),
	 axis.title=element_text(size=2),
	 axis.ticks = element_line(size = 0.1),
	 axis.ticks.length=unit(0.5,"mm"),
	 panel.background = element_rect(fill = "white"),
	 strip.background = element_rect(fill = "white"),
	 legend.text=element_text(size=2),
	 legend.title=element_text(size=2),
	 legend.key=element_blank(),
	 legend.key.size=unit(0.1,"line"),
	 legend.key.height=unit(0.1,"line"),
	 legend.position="bottom",
	 legend.background=element_blank(),
	 legend.box.background=element_blank(),
	 legend.box.spacing=unit(0,"line"),
	 legend.spacing=unit(0,"mm"),
	 legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
	 legend.margin=margin(0, unit = "mm"),
	 plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit = "mm")) +
	 guides(color=guide_legend(title="GENOTYPE"))
	 ggsave("LG12segregatingSites.png",width=100,height=30,units="mm",dpi=1000)
	  
	ggplot(subset(lg12tmp,IND %in% paste0("FIN-PYO-",c(0,20))),aes(x=BIN*BIN.SIZE,y=COUNT)) + geom_line(aes(color=GT),size=0.1) + facet_grid(IND~REF) +
	 scale_x_continuous(breaks=c(0,10,20,30)*10^6) +
	 xlab("Genomic position") + ylab("Number of sites") +
	 theme(strip.text = element_text(size = 4),
	 axis.text=element_text(size=4),
	 axis.title=element_text(size=4),
	 axis.ticks = element_line(size = 0.25),
	 axis.ticks.length=unit(0.5,"mm"),
	 panel.background = element_rect(fill = "white"),
	 strip.background = element_rect(fill = "white"),
	 legend.text=element_text(size=2),
	 legend.title=element_text(size=2),
	 legend.key=element_blank(),
	 legend.key.size=unit(0.1,"line"),
	 legend.key.height=unit(0.1,"line"),
	 legend.background=element_blank(),
	 legend.box.background=element_blank(),
	 legend.box.spacing=unit(0,"line"),
	 legend.spacing=unit(0,"mm"),
	 legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
	 legend.margin=margin(0, unit = "mm"),
	 plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit = "mm")) 
	 ggsave("LG12segregatingSites2.png",width=100,height=30,units="mm",dpi=1000)
	  
	 OUT[["TOTAL"]]=total
	 
	
	 total3=data.frame()
for(lg in paste0("LG",c(1,4,8,12,13))){	 
	for(ind in paste0("FIN-PYO-",c(0,20,21,22,33))){
	 for(ref in c("ver. 6", "ver. 7")){
		 if(ref=="ver. 6"){
		     tmp=subset(total2,IND==ind & REF==ref & CHROM==lg)
		 }
		 if(ref=="ver. 7"){
		     tmp=subset(total2,IND==ind & REF==ref & CHROM==lg)
		 }
		 obs=tmp$STAT
		 SYMBOLS=seq(1,40)
		 
		 mosaic.probs=(1/exp(SYMBOLS))/sum(1/exp(SYMBOLS))
	     non.mosaic.probs=exp(SYMBOLS)/sum(exp(SYMBOLS))
		 non.mosaic.probs=(1/2^(SYMBOLS))/sum(1/2^(SYMBOLS))
	     mosaic.probs=2^(SYMBOLS)/sum(2^(SYMBOLS))
	     non.mosaic.probs=(1/SYMBOLS)/sum(1/SYMBOLS)
	     mosaic.probs=rev(non.mosaic.probs)
	     mosaic.probs=SYMBOLS/sum(SYMBOLS)
	     non.mosaic.probs=rev(mosaic.probs)
	     
	     
	     hmm = initHMM(States=c("MOSAIC","NORMAL"), Symbols=SYMBOLS, startProbs=c(0.01,0.99), transProbs=matrix(c(0.90,0.10,0.005,0.995),2,byrow=T), emissionProbs=matrix(c(mosaic.probs,non.mosaic.probs),2,byrow=T))
	     bw = baumWelch(hmm=hmm, observation=obs, maxIterations=1000000, delta=1E-9, pseudoCount=1)
	     vit = viterbi(hmm=bw$hmm, observation=obs)
	     values=rep(0,length(obs))
	     values[which(vit=="MOSAIC")]=1
	     total3=rbind(total3,data.frame(CHROM=lg,REF=ref,IND=ind,BIN=seq(1,length(obs)),STATE=vit,STAT=obs,NLOCI=tmp$NLOCI,VALUE=values))
	 
	 
	 }
 }
}
	    total3$GENOMICPOSITION=total3$BIN*BIN.SIZE
        total3$YMIN=0.5
        total3$YMAX=1.5
        total3[total3$REF=="ver. 6","YMIN"]=2.5
        total3[total3$REF=="ver. 6","YMAX"]=3.5
        total3$COLOR=paste(total3$REF,total3$STATE)
        total3$CHROM=factor(total3$CHROM,levels=paste0("LG",c(1,4,8,12,13)))

	 
        ggplot(subset(total3,CHROM!="LG12" & IND!="FIN-PYO-0"),aes(x=GENOMICPOSITION)) + geom_rect(aes(ymin=YMIN,ymax=YMAX,xmin=GENOMICPOSITION,xmax=GENOMICPOSITION+BIN.SIZE,fill=COLOR)) + 
        facet_grid(CHROM~IND) + xlab("Genomic position (Mbp)")+scale_x_continuous(breaks=c(0,10,20,30)*10^6,label=c("0","10","20","30")) + scale_fill_manual(values=c("#D55E00","#E69F00","#0072B2","#56B4E9"))+ #c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
	    theme(strip.text.y = element_text(size = 4),strip.text.x = element_text(size = 4),panel.background=element_rect(fill = NA),
        panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y=element_blank(),axis.text.x=element_text(size=4),
        axis.title=element_text(size=4),legend.title=element_blank(),legend.position ="none",legend.key.size = unit(0.1,"line"),
        legend.box.margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),legend.key=element_blank(),
        legend.margin=margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),legend.box.spacing=unit(0,"line"),legend.text=element_text(size=2), panel.spacing = unit(0.2, "lines"),
        axis.ticks.y=element_blank(),axis.ticks.x=element_line(size=0.3),axis.ticks.length=unit(1.5,"mm"),strip.background=element_blank(),plot.margin=margin(t = 1, r = 1, b = 1, l = 1, unit = "mm"))
	    ggsave("mosaicism.png",width=150,height=50,units="mm",dpi=1000)
	
OUT[["TOTAL3"]]=total3

array.out=OUT
    """
) 
  plotMosaicism._execute="once"     
  
/// GENERATE MASKS FOR VARIANT FILTERING 
  
  val snpbility=NamedMap[TextFile]("snpbility")
  val repeat=NamedMap[TextFile]("repeat")
  val positiveMask=NamedMap[TextFile]("postiveMask")
  val refs = Map("v6"->"v6","v7"->"v7")

 for((ref, refName) <- refs){
  snpbility(ref)=BashEvaluate(
      var1=referenceGenomes,
      param1=ref,
      command="""
	        ref=@var1@/@param1@.fa 
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
     var1=referenceGenomes,
     var2=repeatLibrary,
     param1=ref,
     command="""
     ref=@var1@/@param1@.fa
     ls "$ref".fai
     echo $ref
     RepeatMasker -xsmall -gff -dir @folder1@ -pa 10 -lib @var2@ $ref
     awk '!/#/{OFS="\t";print $1,$4,$5}' @folder1@/@param1@.fa.out.gff | sort -k1,1 -k2,2n | bedtools merge -d 100|bedtools slop -i - -g "$ref".fai -b 100|bedtools merge > @out1@
     #rm @folder1@/*
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
  
 //Collect information of regions which haplotype copy has  
  

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
haplotypeContigs._execute="once"
val searchHaplotypes=REvaluate(
var1=v7agp,
var2=laOut,
script="""
haplotypes=data.frame()
lgs=paste0("LG",1:21)
v7.agp=read.table(var1,header=F,stringsAsFactors=F,sep="\t")
for(i in 1:21){
    lg=paste0("LG",i)
    agp=subset(v7.agp,V1==lg & grepl("quiver",V6))
    lep.anchor=read.table(paste0(var2,"/ichr",i,".la"),header=F,stringsAsFactors=F,sep="\t")
    if(i==19){lep.anchor=read.table(paste0(var2,"/iichr",i,".la"),header=F,stringsAsFactors=F,sep="\t")}
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
searchHaplotypes._execute="once"
}    
