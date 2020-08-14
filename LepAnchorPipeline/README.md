# Lep-Anchor wrapper

## Summary
This is a wrapper for the Lep-Anchor pipeline. The minimal wrapper command requires: (1) path to the original reference genome (in fasta format), (2) expected number of chromosomes, and (3) path to at least one linkage map. If paths to liftover chains and long read alignments are provided, the pipeline will also identify and remove haplotypes. More details about the analysis steps and the Lep-Anchor modules can be found at https://sourceforge.net/p/lep-anchor/wiki/Home/

## Installation

Instructions for the installation of required and optional  software: 

- Lep-Anchor: https://sourceforge.net/p/lep-anchor/wiki/Home/
- Lep-MAP3: https://sourceforge.net/p/lep-map3/wiki/LM3%20Home
- minimap2: https://github.com/lh3/minimap2
- Haplomerger2: https://github.com/mapleforest/HaploMerger2/

## Usage
`lepanchor_wrapper.sh -f ref.fasta -n num_chr -c chain_file -p paf_file -m map_file1 -m map_file2 ... -m map_fileN`

## Input files

- `ref.fasta`: reference genome to be anchored (required, can be gz compressed)
- `num_chr`: expected number of chromosomes (required)
- `chain_file`: liftover chains generated with Haplomerger2. These are essential for identifying haplotype contigs (optional, can be gz compressed)
- `paf_file`: long read alignment file in paf format. Can be obtained with minimap2 (optional, compression supported by `<(gunzip -c file.paf.gz)`)
- `map_file1`: Linkage map built with Lep-MAP3 (required, a text file).
- `map_file2`: Additional linkage maps for more accurate anchoring can be provided (optional, text file(s)).

## Output files
- `LA_REF.fa.gz`: anchored reference genome  
- `LA_REF.agp`: agp file describing the genome anchoring  
- `marey*.png`: Marey maps for visual verification, one for each chromosome  
- `scaffolds_chr*.agp`: agp file for each chromosome in scaffolds (each block of linked contigs as a scaffold)


## Analysis steps

The pipeline code is documented in the wrapper script file. The analysis steps are the following: 

1. Identification of contigs that are haplotypes over full length. Executed if a chain file is provided.
2. Liftover of markers from the haplotypes to primary contigs
3. CleanMap: assigns contigs into linkage groups
4. Identification of contigs that cannot be placed in linkage groups
5. Placement and orientation of contigs: finds the correct order and orientation of the contigs that are assigned in the same linkage group. For the first iteration, all unplaced contigs are provided for each linkage group
6. Propagation: find contigs without markers that can be put uniquely to a linkage group based on chain alignments or paf links
7. Second iteration of placement and orientation of contigs: includes only contigs with markers or contigs propagated to exactly one chromosome
8. Pruning out of chains of linked contigs with no linkage map support
9. Removal of overlaps: fixes chimeric contigs so that no region is included multiple times
10. Construction of an agp file for the new assembly and for alternative scaffold assembly
11. Output of contigs that are not haplotypes and not included in the anchored assembly. Add these to the final agp.
12. Generation of fasta sequence of the new reference
13. Construction of Marey maps (R code)

## Replication of stickleback analyses

The analyses for improving the stickleback genomes can be replicated using the wrapper script and the ready-made data files (maps, alignments, chains).

**The data files are available at**  
- [CABVRH02.fasta.gz](https://www.ebi.ac.uk/ena/browser/view/GCA_902500615.2)  (nine-spine raw assembly)
- [9sp_\*_map.txt.gz](../SticklebackAnalyses/Data) (nine-spine maps)  
- [all.chain.9sp.gz](http://wasabi2.biocenter.helsinki.fi/data/nspv7/all.chain.9sp.gz) (nine-spine liftover chains)  
- [alnPB.paf.gz](http://wasabi2.biocenter.helsinki.fi/data/nspv7/alnPB.paf.gz) (nine-spine PacBio alignments)  
<!-- - [3sp_map.txt.gz](../SticklebackAnalyses/Data) (three-spine maps)  -->
<!-- - [all.chain.3sp.gz](http://wasabi2.biocenter.helsinki.fi/data/nspv7/all.chain.3sp.gz) (three-spine chains)  -->

**Alignment of PacBio reads (all.fastq.gz) with minimap2**  
`./minimap2 -x map-pb -t 16 Pungitius_RawAssembly.fa.gz all.fastq.gz > alnPB.paf`

**Addition of the LG19 inversion haplotypes (pilon\_allele\*.fasta) to the raw assembly**  
`cat Pungitius\_RawAssembly.fa.gz <(cat pilon_allele*.fasta|gzip) > Pungitius_RawAssemblyAB.fa.gz`

**Repeat masking of three-spine genome (winMasker is in the HaploMerger directory)**  
`../winMasker/windowmasker -mk_counts -in Gac-HiC_contigs.fa -out counts.txt`

`../winMasker/windowmasker -ustat counts.txt -in Gac-HiC_contigs.fa -out Gac_HiC_contigs_sm.fa -outfmt fasta`

`gzip Gac_HiC_contigs_sm.fa`

**Running of HaploMerger**  
File `haplomerger.tar.gz` contains a clean HaploMerger project. This was run separately for (1) nine-spine, (2) nine-spine+inversion haplotypes and (3) three-spine. 

To do that, change `genome_name` in `run_all.batch` to
1. `Pungitius_RawAssembly` and add `Pungitius_RawAssembly.fa.gz` to the project folder
2. `Pungitius_RawAssemblyAB` and add `Pungitius_RawAssemblyAB.fa.gz` to the project folder
3. `Gac_HiC_contigs_sm` and add `Gac_HiC_contigs_sm.fa.gz` to the project folder

**Running of Lep-Anchor wrapper for nine-spine**  
`./lepanchor_wrapper.sh -m 9sp_Hel_map.txt -m 9sp_O13_map.txt -m 9sp_O23_map.txt -n 21 -c  "all.chain.gz|awk -f addpilon.awk" -f NineSpine.fasta.gz`

**Manual postprocessing steps**  
Some non-standard postprocessing steps were performed manually: 
- versions with inversion haplotypes A and B
- fixing of chromosome orientation (chromosome orientation is arbitrary and manual curation is needed for comparisons between assembly versions)
- fixing of some chimeric contigs
- three-spine data

