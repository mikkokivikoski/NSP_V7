### SticklebackAnalyses/Scripts/

\*.scala: Anduril2 command files for downstream analyses of the nine-spined and three-spined stickleback genome assemblies.

For replication, modify the paths for the input files and the programs. Pipelines should be run in the following order:

**Nine-spined stickleback**
1. mapping.scala
2. ppunAnalyses.scala
3. variantAnalysis.scala

**Three-spined stickleback**
1. gacuMappingVariantCalling.scala
2. gacuAnalyses.scala
3. gacuVariantAnalysis.scala

lg19InversionAssembly.txt: description of the reassembly of LG19 inversion haplotypes.