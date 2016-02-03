# ABCC4 Analysis

Pipeline used to analyse data in ABCC4 analysis.

## Extract ABCCC4/SLC22A6-8 genes from VCF.

Extract only the regions from the gene of interest.

out.vcf is the final VCF file from the overall QC procedure.

```
    bgzip out.vcf
    tabix -p out.vcf
    bcftools -S beds/genes.bed
```

## Generate common variant P-values
