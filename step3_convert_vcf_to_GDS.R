# step3_convert_vcf_to_GDS.R
# GDS file is needed as input into the PC calculation script

# install and load SeqArray for conversion
#if (!requireNamespace('BiocManager', quietly = TRUE))
#    install.packages('BiocManager')
#BiocManager::install('SeqArray')

library(SeqArray)

filepath="/home/user/"

# If there are more than one files in vcf.fn, seqVCF2GDS will merge all VCF files together if they contain the same samples. It is useful to merge multiple VCF files if data are divided by chromosomes.

# enter ALL vcf_files into an object to be loaded into the seqVCF2GDS() conversion function

vcf_files <- paste0(filepath, "/merged_file_chr", 1:22, ".vcf.bgz")
output_filename <- paste0(filepath, "gds_files/36K_1000Ghgdp_JF_chr1-22.gds")

seqVCF2GDS(vcf.fn = vcf_files, out.fn = output_filename)