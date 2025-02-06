#!/bin/bash
# step 1 FIND OVERLAP of QC/filtered biallelic genotypes
# step1_pull_jf_1000g_snp_overlap.sh


filepath_jf="/ld_pruned"
filepath_hgdp1000g="/phased_haplotypes_v2"

for i in {1..22}

do

bcftools view -R $filepath_jf/filtered.biallelic.genotypes.chr${i}.ld_pruned.ALL.vcf.bgz -Oz -o gnomad_jf_v4.0.sites.chr${i}.vcf.bgz $filepath_hgdp1000g/hgdp1kgp_chr${i}.filtered.SNV_INDEL.phased.shapeit5.bcf

bcftools index -t gnomad_jf_sites.chr${i}.vcf.bgz

#done

mkdir -p jf_1000g_overlap
mv gnomad*bgz ./jf_1000g_overlap