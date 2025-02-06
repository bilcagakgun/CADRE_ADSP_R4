## step4_run_PCA_36K.R
# RUN PCA analysis

#install.packages(c("data.table", "tidyr", "tidyverse", "rlist", "MASS"))
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install(c("GENESIS", "GWASdata", "gdsfmt", "SNPRelate", "GWASTools", "SeqArray", "SeqVarTools"))

library(data.table)
library(tidyr)
#library(tidyverse)
#library(rlist)
library(GENESIS)
#library(GWASdata)
library(gdsfmt)
library(SNPRelate)
library(MASS)
library(GWASTools)
library(SeqArray)
library(SeqVarTools)

# add parallel jobs spec
library(BiocParallel)          
register(MulticoreParam(workers = 2)) #should be same as in job header sbatch flag

showfile.gds(closeall=TRUE)

getwd()
#setwd("/s3buckets/ancestry_specific")
#######################################
# READ IN ancestry_specific lists

#dem = read.delim("ADSP_36K_samples_and_REF_IDS.txt")
dem = read.delim("./sep_lists/GRP1_over55_sampleIDS.txt")
ref = read.delim("./HGDP_1KG_IDs.high_quality.3922.txt")
###############################################
#Genotype data, Calculating PCs and GRM matrix

# files_list = list.files(recursive = F)
# Read GDS files from step 3
# open seqArray file
#seqClose("/s3buckets/gds_files/36K_1000Ghgdp_JF_chr1-22.gds")
gds = seqOpen("./36K_1000Ghgdp_JF_chr1-22.gds")

# filter for GRP1 samples
seqSetFilter(gds, sample.id=dem$SampleID)
seqExport(gds, "GRP1_36K_1000Ghgdp_JF_chr1-22_highqualityref.gds")
seqClose(gds)

gds = seqOpen("GRP1_36K_1000Ghgdp_JF_chr1-22_highqualityref.gds")

#----------------------------------
# KING relatedness

ibd_KING = snpgdsIBDKING(gds, sample.id=dem$SampleID, snp.id=NULL,
                         autosome.only = TRUE, maf=NaN, type="KING-robust", family.id=NULL, missing.rate = NaN)
saveRDS(ibd_KING, file = "GRP1_36K_hgdp1000g_IBDKING_snpgds.rds")
ibd_KING = readRDS("GRP1_36K_hgdp1000g_IBDKING_snpgds.rds")

data_ibdKING <- snpgdsIBDSelection(ibd_KING)
saveRDS(data_ibdKING, file = "GRP1_36K_hgdp1000g_IBDKING_snpgds_data.rds")
data_ibdKING = readRDS("GRP1_36K_hgdp1000g_IBDKING_snpgds_data.rds")

#Extract matrix of kinship coefficients
ibd_KING_matrix = ibd_KING$kinship
saveRDS(ibd_KING_matrix, file = "GRP1_36K_hgdp1000g_IBDKING_matrix.rds")
ibd_KING_matrix = readRDS("GRP1_36K_hgdp1000g_IBDKING_matrix.rds")

dimnames(ibd_KING_matrix) = list(ibd_KING$sample.id, ibd_KING$sample.id)

# close gds file
#seqClose(gds)
#showfile.gds(closeall=TRUE)
#---------------------------------
#PCAir calculation and GRM calculation


# trim down the gds file to the same high quality ref variants and then also reorder them to match 1:1 with sample-ids
### Annotated SeqVarData
annot <- dem
annot <- annot[match(seqGetData(gds,"sample.id"), annot$SampleID), ]
annot <- as.data.frame(annot)
colnames(annot)[1] <- c("sample.id")
metadata <- data.frame(labelDescription=c("sample.id"), row.names=names(annot))
annot <- AnnotatedDataFrame(annot, metadata)

all.equal(annot$sample.id, seqGetData(gds, "sample.id"))

seqData <- SeqVarData(gds, sampleData=annot)

###

# run pcair
mypcair <- pcair(seqData, 
                 kinobj = ibd_KING_matrix, kin.thresh=2^(-11/2),
                 divobj = ibd_KING_matrix, div.thresh=-2^(-11/2),
                 snp.include = NULL,
		 unrel.set = ref$IDs,  
                 sample.include = dem$SampleID)
pcair_sample_sets <- pcairPartition(kinobj = ibd_KING_matrix, divobj = ibd_KING_matrix)

# iterate results into blocks
iterator <- SeqVarBlockIterator(seqData, variantBlock=20000, verbose=FALSE)

#--------------------------------------
#PCRelate calculation => PCs and GRM
# pcrelate on seqData iterator

pcrelate_pcairKING <- pcrelate(iterator, pcs = mypcair$vectors[,1:3],
                               scale = "overall", training.set = mypcair$unrels,
                               sample.include = dem$SampleID)

grm_matrix <- pcrelateToMatrix(pcrelate_pcairKING)

mypcair <- pcair(iterator,
                 kinobj = grm_matrix, kin.thresh=2^(-11/2),
                 divobj = ibd_KING_matrix, div.thresh=-2^(-11/2),
                 snp.include = NULL,
                 sample.include = dem$SampleID)


pcs = mypcair$vectors
pcs = data.frame(pcs)
colnames(pcs) = paste0("pc", c(1:ncol(pcs)))
pcs$SampleID = row.names(pcs)

saveRDS(grm_matrix, file = "GRP1_36K_hgdp1000g.grm_matrix.rds")
saveRDS(pcs, file = "GRP1_36K_hgdp1000g.pcs.rds")
saveRDS(mypcair, file  = "GRP1_36K_hgdp1000g.mypcair.rds")