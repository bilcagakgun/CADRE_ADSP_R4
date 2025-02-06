#!/bin/bash
# step2_sbatch_merge_jf.sh
# OVERLAP 1000G/HGDP REF

# File 1 is JF's original pruned vcf files
# File 2 is the overlap of JF's pruned vcf files WITH gnomad 1000G/HGDP ref. # of SNPs drops slightly

# Now, need the hgdp1000g SNPs from JF's overlap, while including all the sample data from both
# File 1:
filepath_jf="/ld_pruned"
# File 2:
filepath_hgdp1000g_jfoverlap="/jf_1000g_overlap"

for i in {1..22}

do

# Merge

bcftools merge $filepath_jf/filtered.biallelic.genotypes.chr${i}.ld_pruned.ALL.vcf.bgz $filepath_hgdp1000g_jfoverlap/gnomad_jf_v4.0.sites.chr${i}.vcf.bgz -Oz -o merged_file_chr${i}.vcf.bgz

done

##mkdir -p ./JF_1000HGDP_merged
##mv merged_file_chr*.vcf.bgz ./JF_1000HGDP_merged
