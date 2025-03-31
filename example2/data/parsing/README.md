### Estimating the SFS

This Snakefile prepares VCF files for SFS parsing by subsetting them and annotating sites with ancestral states, then estimates the SFS. You will need `bcftools` to run the pipeline- you can obtain it here https://samtools.github.io/bcftools/.

To run this script, first obtain lifted-over GVCFs (hg38) for archaic genome sequences using [../liftover/Snakefile](../liftover/Snakefile). Also download the resequenced 1000 genomes sequences available at https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/ and set the variable `vcfs_1kG` in [Snakefile](Snakefile) to the path to the directory that contains them. All other required materials should either be already present in the repository or downloaded automatically by [Snakefile](Snakefile).

The masks in [data/masks/exons/](data/masks/exons/) were obtained from Ensembl. The promoter annotations that we used are not published, but the pipeline we make available avoids masking this category of sites- thus the SFS and `L` it parses will have higher magnitudes than the ones used in the examples.