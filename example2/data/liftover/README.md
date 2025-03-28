### Lifting archaic genomes over from hg19 to hg38

This directory contains a Snakefile which lifts GVCF files for four archaic humans over from genome build hg19 to hg38. GVCFs are VCF files that include  (genotyped) monomorphic sites.

To use this pipeline, you will need GVCF files in hg19 and the `bcftools` `liftover` plugin. We downloaded GVCFs from the `Vindija33.19/`, `Denisova/`, and `Altai/` directories from http://ftp.eva.mpg.de/neandertal/Vindija/VCF/ and the `VCF/` directory from http://ftp.eva.mpg.de/neandertal/Chagyrskaya/ into subdirectories of `./data/raw_genomes_hg19/` named `Vindija/`, `Denisova/`, `Altai/`, and `Chagyrskaya` respectively. 
The Snakefile should download other all required materials upon execution.
Documentation of `bcftools` and links to download it can be found at https://samtools.github.io/bcftools/. 
It can also be installed from several Linux package managers, but note that the `liftover` plugin requires version $\geq$ 1.20. 
The `liftover` plugin is documented at https://github.com/freeseek/score?tab=readme-ov-file#liftover-vcfs, and binaries can be downloaded from https://software.broadinstitute.org/software/score/score_1.20-20240927.zip. 
Once downloaded, extract the binaries to `~/bin/bcftools_plugins/` or another directory of your choice and invoke `export BCFTOOLS_PLUGINS=~/bin/bcftools_plugins/`.
