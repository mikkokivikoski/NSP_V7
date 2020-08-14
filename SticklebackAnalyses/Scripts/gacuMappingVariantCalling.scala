#!/usr/bin/env anduril

import anduril.builtin._
import anduril.tools._
import anduril.anima._
import org.anduril.runtime._
//$OPT --threads 16


object threespineReference {


val gasAcuRef = INPUT("Gac-HiC_revised_genome_assembly.fa")
val gasAcuRefLA10X = INPUT("gacuNEW/3SP2.fasta")
val gasAcuRefIndex = INPUT("Gac-HiC_revised_genome_assembly.fa.fai")
val gasAcuRefLA10XIndex = INPUT("gacuNEW/3SP2.fasta.fai")
val gasAcuRefAgp = INPUT("gacuNEW/contigs.agp")
val gasAcuRefLA10XAgp = INPUT("gacuNEW/3SP2.agp")

val unAssingnedContigs = INPUT("3SP2 unassigned contigs")

val contiglist= INPUT("gacuNEW/Gac-HiC_contigs.fa")

val gacuReferences = Map("gasAcuRef" -> gasAcuRef, "gasAcuRefLA10X" -> gasAcuRefLA10X,
                         "gacuRefIndex" -> gasAcuRefIndex,"gacuRefLA10XIndex" -> gasAcuRefLA10XIndex,
                         "gacuRefAgp" -> gasAcuRefAgp,"gacuRefLA10XAgp" -> gasAcuRefLA10XAgp)

val prepareReferenceGenomes = BashEvaluate(
    var1=gasAcuRefLA10X,
    var2=contiglist,
    var3=unAssingnedContigs,
    var5=gasAcuRef,
    command="""
    cd @folder1@
    for i in I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI XVII XVIII XIX XX XXI;do
        samtools faidx @var5@ chr"$i" >> gasAcu.fasta
    done
    samtools faidx gasAcu.fasta
    bwa index gasAcu.fasta
    samtools dict gasAcu.fasta > gasAcu.dict
    
    for i in `seq 1 21`;do
        samtools faidx @var1@ CHR"$i" >> gasAcuLA10X.fasta
    done
    samtools faidx @var2@ -r @var3@ >> gasAcuLA10X.fasta
    samtools faidx gasAcuLA10X.fasta
    bwa index gasAcuLA10X.fasta
    samtools dict gasAcuLA10X.fasta > gasAcuLA10X.dict
    
    echo DONE
    """

)
prepareReferenceGenomes._execute="once"

val WGSbwa = BashEvaluate(
		    var1=gasAcuRef,
		    command="""		        
		        cd @folder1@ 
				for i in SRR5626529; do
	                 bwa mem -t8 -M -R "@RG\tID:$i\tSM:$i\tPL:illumina\tLB:$i\tPU:L0" \
					 @var1@ \
					 "$i"_1.fastq.gz \
					 "$i"_2.fastq.gz \
					 | samtools view -h -b -o Gacu_"$i"_aligned.bam -
                 
                     samtools sort -T @folder2@/Gacu_"$i" -O bam -o Gacu_"$i"_sorted.bam Gacu_"$i"_aligned.bam
		             samtools rmdup Gacu_"$i"_sorted.bam Gacu_"$i"_dedup.bam
		             samtools index Gacu_"$i"_dedup.bam Gacu_"$i"_dedup.bai   
			    done
		    """
		    )
   WGSbwa._execute="once"    
val WGSbwaNew = BashEvaluate(
		    var1=gasAcuRefLA10X,
		    var2=prepareReferenceGenomes.folder1,
		    command="""		        
		        cd @folder1@ 
				for i in SRR5626529; do
	                 bwa mem -t8 -M -R "@RG\tID:$i\tSM:$i\tPL:illumina\tLB:$i\tPU:L0" \
					 @var2@/gasAcuLA10X.fasta \
					 "$i"_1.fastq.gz \
					 "$i"_2.fastq.gz \
					 | samtools view -h -b -o GacuNew_"$i"_aligned.bam -
                      
                    samtools sort -T @folder2@/GacuNew_"$i" -O bam -o GacuNew_"$i"_sorted.bam GacuNew_"$i"_aligned.bam
		            samtools rmdup GacuNew_"$i"_sorted.bam GacuNew_"$i"_dedup.bam
		            samtools index GacuNew_"$i"_dedup.bam GacuNew_"$i"_dedup.bai  
                      
			    done
		    """
		    )
   WGSbwaNew._execute="once"    


val individuals = Map("SRR5626529"->"SRR5626529")


val gacuMpileup = BashEvaluate(
    var1=WGSbwa.folder1,
    var2=gasAcuRef,
    command="""
      cd @folder1@
	  bcftools mpileup --threads 8 -Ou -f @var2@ @var1@/Gacu_SRR5626529_dedup.bam | bcftools call -Oz -mv > Gacu_SRR5626529.vcf.gz
	  bcftools index Gacu_SRR5626529.vcf.gz
    """
)
val gacuNewMpileup = BashEvaluate(
    var1=WGSbwaNew.folder1,
    var2=prepareReferenceGenomes.folder1,
    command="""
      cd @folder1@
	  bcftools mpileup --threads 8 -Ou -f @var2@/gasAcuLA10X.fasta @var1@/GacuNew_SRR5626529_dedup.bam | bcftools call -Oz -mv > GacuNew_SRR5626529.vcf.gz
	  bcftools index GacuNew_SRR5626529.vcf.gz
    """
)

}
