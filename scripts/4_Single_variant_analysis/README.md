---

## Description
### Input files
- WGS VCF and index files
- Plink files for the same WGS dataset (bed, bim fam)
- Covariates (this file should include diagnosis, age, sex, apoe4 dosages, and PCA results)

---

### Execution

1. Single variant testing [scripts/4_Single_variant_analysis/1_SingleVariantTesting.sh](1_SingleVariantTesting.sh)
	- The "saige_pheno_file.tsv" file generated in the PCA step can be used as a covar file.

---

## Dependencies
### Tools
- PLINK 1.9
- SAIGE

