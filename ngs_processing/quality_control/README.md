# Quality control processing of SNP variants for Case/control.

Automates a approach for QC from the following paper.

Purcell, Shaun M., et al. "A polygenic burden of rare disruptive mutations in schizophrenia." Nature 506.7487 (2014): 185-190.

Similar approach will also work for indels, but currently only for biallelic ones. 

# Software Requirements

- PLINK/SEQ (v0.10) (http://atgu.mgh.harvard.edu/plinkseq/index.shtml)
    - with HG19 database
- siteQCresequencing (master) (https://github.com/smilefreak/siteQCresequencing)
- site_qc_pipeline (master) (https://github.com/smilefreak/site_qc_pipeline)
- indivresequencing_QC (master) (https://github.com/smilefreak/indivresequencing_QC)
- plinkseq\_utilies (master) (https://github.com/smilefreak/plinkseq_utilities)

### Convert SNPmax phenotypes to pseq format. 

Takes tab-delimited output similar to that from SNPmax and writes it out in a format for adding to a
PLINK/seq project

```
    snpmax_to_pseq -i phenotypes/phenotypes.txt -o phenotypes/phenotype_pseq.txt
```

### Perform individual Level QC

### Perform Site Level QC







 
