#!/bin/bash
##########################################################
#
#  Lep-Anchor wrapper
#   
#    This file is part of Lep-Anchor.
#
#    Lep-Anchor is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Lep-Anchor is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Lep-Anchor.  If not, see <http://www.gnu.org/licenses/>.
#
#    Copyright (C) 2020 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki
#
#  usage: lepanchor_wrapper.sh -f ref.fasta -n num_chr -c chain_file -p paf_file -m map_file1 -m map_file2, ... 
#
#  output: LA_REF.fa.gz anchored reference genome
#          LA_REF.agp   agp file describing the genome
#          marey*.png   Marey maps for visual verification
#          chr*.agp     agp files for each chromosome
#          scaffolds_chr*.agp agp files for each chromosome in scaffolds (each block of linked contigs as a scaffold)
##########################################################

#if [ "$#" -ne 3 ]; then
#    echo "At least three input parameters must be provided"
#    exit 1
#fi


#parse parameters

function print_usage()
{
echo "##########################################################"
echo "#"
echo "#  Lep-Anchor wrapper"
echo "#    This file is part of Lep-Anchor."
echo "#"
echo "#    Lep-Anchor is free software: you can redistribute it and/or modify"
echo "#    it under the terms of the GNU General Public License as published by"
echo "#    the Free Software Foundation, either version 3 of the License, or"
echo "#    (at your option) any later version.echo"
echo "#"
echo "#    Lep-Anchor is distributed in the hope that it will be useful,"
echo "#    but WITHOUT ANY WARRANTY; without even the implied warranty of"
echo "#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
echo "#    GNU General Public License for more details."
echo "#"
echo "#    You should have received a copy of the GNU General Public License"
echo "#    along with Lep-Anchor.  If not, see <http://www.gnu.org/licenses/>."
echo "#"
echo "#    Copyright (C) 2020 Pasi Rastas, pasi.rastas@gmail.com, University of Helsinki"
echo "#"
echo "#  usage: lepanchor_wrapper.sh -t threads -f ref.fasta -n num_chr -c chain_file -p paf_file -m map_file1 -m map_file2, ... "
echo "#"
echo "#  download Lep-Anchor by"
echo "#  wget https://sourceforge.net/projects/lep-anchor/files/binary%2Bcode.zip/download -O la.zip;unzip la.zip"
echo "#"
echo "#"
echo "##########################################################"
}

while getopts ":n:c:p:m:f:t:" OPTION; do
	case ${OPTION} in
	t)
	THREADS=$OPTARG;;
	n)
	CHR=$OPTARG;;
	c)
	CHAIN=$OPTARG;;
	p)
	PAF=$OPTARG;;
	m)
	MAP="$MAP $OPTARG";;
	f)
	REF=$OPTARG;;
	*)
        echo "Incorrect options provided"
        exit 1;;
    esac
done

if [[ $MAP =~ ^$ ]];then
	print_usage
	echo "Please provide at least one map"
	exit 1
fi

if [[ $REF =~ ^$ ]];then
    if [[ ! -e "contigs.length" ]];then
	print_usage
	echo "Please provide either reference fasta or contigs.length file"
	exit 1
    fi
fi

if [[ $CHR =~ ^[0-9]+$ ]];then
	echo "Number of chromosomes = $CHR"
else
	print_usage
	echo "Please provide parameter n (number of chromosomes)"
	exit 1
fi

if [[ $CHAIN =~ ^$ ]];then
	echo "No chain provided"
	CHAIN="/dev/null"
fi

if [[ $PAF =~ ^$ ]];then
	echo "No paf provided" 
	PAF="/dev/null"
fi

#number of threads
if [[ $THREADS =~ ^$ ]];then
	THREADS=8
fi

#get Lep-Anchor
#wget https://sourceforge.net/projects/lep-anchor/files/binary%2Bcode.zip/download -O la.zip
#unzip la.zip

#where Lep-Anchor binaries are located
LABIN=bin/

#java runtime located here
JAVA=java

#parallel 
if ! command -v "parallel"
then
	echo "command parallel not found, using only one thread"
	PARALLEL=bash
else
	PARALLEL="parallel --jobs $THREADS"
fi

if [[ ! $REF =~ ^$ ]];then
	echo "calculating contigs.length file..."
	echo "gunzip -fc $REF|awk -f contigLength  >contigs.length"|bash
fi


echo "finding full haplotypes..."
echo "gunzip -fc $CHAIN|awk -f findFullHaplotypes >fullHaplotypes50.txt"|bash
wc -l fullHaplotypes50.txt


echo "running liftoverHaplotypes for all input maps..."
for i in $MAP
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN LiftoverHaplotypes map=$i haplotypes=fullHaplotypes50.txt chain=- >$i.liftover"
done|$PARALLEL

#make input for CleanMap
for i in $MAP
do
cat $i.liftover
done|sort -V >map_all_sorted.liftover

#CleanMap
echo "running CleanMap..."
$JAVA -cp $LABIN CleanMap map=map_all_sorted.liftover >map_all.clean

#Map2Bed
echo "running Map2Bed..."
$JAVA -cp $LABIN Map2Bed map=map_all.clean minQuality=1 contigLength=contigs.length >map.bed

#find contigs not put into chromosomes
cut -f 1 contigs.length|grep -v -w -F -f <(cut -f 2 fullHaplotypes50.txt; cut -f 1 map.bed) >not_used.txt

grep -w -F -f not_used.txt contigs.length|awk '{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=21;++i) print s i}' >chr0.bed
cat map.bed chr0.bed >map_extra.bed


#PlaceAndOrientContigs
echo "running PlaceAndOrientContigs (first iteration)..."
for i in $(seq $CHR)
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN PlaceAndOrientContigs bed=map_extra.bed chromosome=$i map=$MAP chain=- paf=$PAF keepEmptyIntervals=1 >chr$i.la 2>chr$i.la.err"
done|$PARALLEL


#propagate
echo "running propagate..."

awk -f propagate chr*.la >tmp1.la
awk -f propagate tmp1.la >tmp2.la
i=2

while ! cmp -s "tmp$i.la" "tmp$(( $i-1 )).la" ;do
	awk -f propagate tmp$i.la >tmp$[$i+1].la
	i=$[$i+1]
done

#create prop*.la
awk '/^[^#]/{++d[$1 "\t" $7+0 "\t" $8+0]; data[++line]=$0}END{for (i = 1; i <= line; ++i) {$0=data[i];if (d[$1 "\t" $7+0 "\t" $8+0] == 1) fn="prop"$5".la"; else if ($5==1) fn="prop0.la"; else fn=""; if (fn != "") print $0>fn}}' tmp$i.la

#create a new bed by combining prop[1-9]*.la and map.bed
awk '(NR==FNR){print;c[$1]}(NR!=FNR && !($1 in c)){print $1 "\t" $7+0 "\t" $8+0"\t?\t"$5}' map.bed prop[1-9]*.la >map_prop.bed

#PlaceAndOrientContigs
echo "running PlaceAndOrientContigs (second iteration)..."
for i in $(seq $CHR)
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN PlaceAndOrientContigs bed=map_prop.bed chromosome=$i map=$MAP chain=- paf=$PAF keepEmptyIntervals=1 >ichr$i.la 2>ichr$i.la.err"
done|$PARALLEL

#pruning contig blocks without map support
for i in $(seq $CHR)
do
        awk -f prune ichr$i.la >ichr${i}_pruned.la
done 2>pruned.la

#remove overlap(s)
awk -f removeOverlaps.awk map_prop.bed ichr*_pruned.la >iall.la

#construct agp files
for i in $(seq $CHR)
do
awk -vn=$i '($5==n)' iall.la|awk -vprefix="LG" -vlg=$i -f makeagp_full2.awk - >chr$i.agp
awk -vn=$i '($5==n)' iall.la|awk -vprefix="LG" -vlg=$i -f makeagp2.awk - >scaffolds_chr$i.agp
done

#find contigs not used
cut -f 1 contigs.length|grep -v -w -F -f <(cut -f 2 fullHaplotypes50.txt;awk '($5!="U"){print $6}' chr*.agp) >not_used_final.txt

grep -F -w -f not_used_final.txt contigs.length|awk '{print $1,1,$2,1,"W",$1,1,$2,"+"}' >not_used.agp

cat chr*.agp not_used.agp >REF_LA.agp
#one could use scaffolds_chr*.agp as well instead of chr*.agp

#make final fasta
if [[ ! $REF =~ ^$ ]];then
	echo "constructing final fasta (REF_LA.fa.gz)..."
	echo "gunzip -fc $REF|awk -f makefasta - REF_LA.agp|gzip >REF_LA.fa.gz"|bash
fi

#construct Marey map
echo "constructing Marey maps... (marey*.png)"
j=1
for m in $MAP
do
	for c in $(seq $CHR)
	do
	awk -vn=$c '($3==n)' $m.liftover|awk -f liftover chr$c.agp -|awk -vm=$j '(/LG/ && NR>=4){if (NF>4) s=0.5; else s=1;print $1"\t"$2"\t"$3"\t"m"\t"s*($4+$5)}'
	done
	j=$[$j + 1]
done|gzip >marey.data.gz

Rscript plot_marey.R

#TODO: lepanchor_wrapper_step2.sh

