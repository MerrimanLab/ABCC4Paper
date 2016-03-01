# Quality control processing of SNP variants for Case/control.

Automates a approach for QC from the following paper.

Purcell, Shaun M., et al. "A polygenic burden of rare disruptive mutations in schizophrenia." Nature 506.7487 (2014): 185-190.

Similar approach will also work for indels, but currently only for biallelic ones. 

It is also worth noting that these steps can be done on a desktop computer. 

# Software Requirements

- PLINK/SEQ (v0.10) (http://atgu.mgh.harvard.edu/plinkseq/index.shtml)
    - with HG19 database
- siteQCresequencing (master) (https://github.com/smilefreak/siteQCresequencing)
- site_qc_pipeline (master) (https://github.com/smilefreak/site_qc_pipeline)
- indivresequencing_QC (master) (https://github.com/smilefreak/indivresequencing_QC)
- plinkseq\_utilies (master) (https://github.com/smilefreak/plinkseq_utilities)
- plink (v1.07) (http://pngu.mgh.harvard.edu/~purcell/plink/)
- tabix and bgzip (v1.2.1) (http://www.htslib.org/doc/tabix.html)


### Perform individual Level QC

#### Upfront filtering based on MAC and biallelic SNPs.

```
    cat ../gatk/vqsr_snp_names.vcf | bcftools filter -e "MAC == 0" | bcftools view -m2 -M2 | filter.vcf
```

####Generate input file for SHINY app

```

snpqc prep_indivqc -p  phenotypes/phenotype_pseq.txt -o phenotypes/phenotypes.txt -c HYPER -r ~/rm_G6174_G5913_AT0721/hg19 filter.vcf

```

Indiv contains all the metrics for Quality control processing.

####Use indivresequencing_QC SHINY App.

load ```indiv.txt``` into app.

Samples were filtered by visualising and discovering clear outliers for informative metrics. Some are listed below. 
- NHET (number of heterozygotes)
- NALT (number of alternate alleles)
- NSING (number of singleton variants in a sample)
- NVAR (number of variants)

Download samples remaining and save in a file ```samples_remaining.txt``` 

Extract samples from vcf file. 

```
    bgzip -c ../gatk/vqsr_snp_names.vcf  > all_samples.vcf.gz
    tabix -p all_samples.vcf.gz
    bcftools view -S samples_remaining all_samples.vcf.gz > remove_bad_samples.vcf
```

For this dataset, 3 samples G6174,G5913,AT0721 were removed as a result of this step.

### Perform Site Level QC

#### Generate input file for SHINY application

```
    snpqc prep_site_qc -i remove_bad_samples.vcf -p phenotypes/phenotype_pseq.txt -c HYPER  -f "geno=DP:ge:10" -o site_qc.txt 

```


```geno=DP:ge:10``` means than any variant with depth less than 10 or genotype quality less than 10 will be set to missing for analyses. 

This is just a mask from PLINK/Seq, see (https://atgu.mgh.harvard.edu/plinkseq/masks.shtml) for more information.

#### siteQCresequencing application

Load ```site_qc.txt``` into application.

Samples were filtered using a number of metrics, specifically for this dataset the options were.

- GQM=30
- SAMPLE_REF_AB=0.1 
- SAMPLE_HOM_AB=0.9
- SAMPLE_HET_AB=0.2
- SAMPLE_HET_AB=0.8
- PROP_DEVIANT_HET_AB_20_80=0.1
- PROP_DEVIANT_HET_AB_30_70=0.3

Sites were downloaded and saved in ```sites.txt```


#### Post-process VCF files to generate final VCFs.

Create final VCF file.

```
    snpqc postqc -s sites.txt  -v out.vcf -p plink -o out 
```

