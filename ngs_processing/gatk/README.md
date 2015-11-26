## Processing of the NGS data from wustl for the resequencing experiment.

Sequencing data was provided in SAM format, this first step was to convert all
the SAMs back into fastq files. 

http://gatkforums.broadinstitute.org/discussion/2908/howto-revert-a-bam-file-to-fastq-format

## Software Requirements

- GATK (v3.3) (https://www.broadinstitute.org/gatk/download/)
- samtools (v1.1-2) (http://www.htslib.org/)
- bcftools (v1.1-6) (https://samtools.github.io/bcftools/bcftools.html)
- bwa (v0.7.12-r1039) (http://bio-bwa.sourceforge.net/)
- Picard tools (v1.41) (http://broadinstitute.github.io/picard/) 

## Resources 

Broad institute bundle contains reference genome and auxilary files used in the GATK best-practices

Best practices - https://www.broadinstitute.org/gatk/guide/best-practices 

Bundle data - ftp://ftp.broadinstitute.org/bundle/2.8/b37/

- ```references/hapmap.vcf``` (HapMap genotypes and sites VCFs)
- ```references/omni.vcf``` (OMNI 2.5 genotypes for 1kg samples)
- ```references/mills.vcf``` (Mills gold standard set of indels)
- ```references/dbsnp.vcf``` (dbsnp v138.vcf)
- ```references/1000G_phase1.indels.b37.vcf``` (phase 1 indels for the 1000 genomes)
- ```references/human_g1k_v37_decoy.fasta``` GRCh37 decoy reference genome.

## Pipeline

### Bash Variables

Need to set the following variables to your environment.

```bash
    export PICARD=<path to picard.jar>
    export GATK=<path to GATK.jar>
    export REFERENCE=references/human_g1k_v37_decoy.fasta
```

### BAM2Fq
Sequencing data was provided in SAM format, this first step was to convert all
the SAMs back into fastq files. 

http://gatkforums.broadinstitute.org/discussion/2908/howto-revert-a-bam-file-to-fastq-format

Following assumes all the bams are located in your current working directory.

Shuffles bam and then turns them back into a fastq file.

```bash
    parallel "samtools bamshuf -uOn 128 {} tmp > {/.}.shuffled.bam" ::: *.bam
    parallel "samtools bam2fq -a {} > {/.}.fq" :::  *.shuffled.bam
    parallel "gzip {}" ::: *.fq
```

### Alignment

All alignments go into the bwa\_mem folder.

Align
```bash
     mkdir -p bwa_mem
     ./scripts/alignments.sh *.fq
```
Sort
```bash
    ./scripts/sort_sam.sh  bwa_mem/*sam
```
MarkDuplicates
```bash
    parallel "./scripts/mark_duplicates.sh" ::: bwa_mem/*bam
```
Places BAM files with duplicates marked in the ```rmdup``` directory. 

### Picard Metrics

Generate a number of metrics using picardmetrics. 

```
mkdir -p metrics
cd metrics
parallel "picardmetrics run -o picard_metrics {}" ::: ../bwa_mem/*bam
```

### Indel Realignment  

https://www.broadinstitute.org/gatk/guide/article?id=2800

First step generate realignment intervals.
```
    parallel  "java -Xmx16g -jar ${GATK} -T  RealignerTargetCreator -R ${REFERENCE} -I {} -o {/.}.intervals \
               --known references/1000G_phase1.indels.b37.vcf --known references/mills.vcf \
               -L references/resequencing_regions.bed " ::: rmdup/*bam
```
Perform the realignment.
```
    parallel  "java -Xmx16g -jar ${GATK} -T  IndelRealigner -R ${REFERENCE} -I {} --targetIntervals {/.}.intervals \ 
               -o realign/{/.}.realigned_indels.bam -known references/mills.vcf --known references/1000G_phase1.indels.b37.vcf \ 
                -L references/resequencing_regions.bed" ::: rmdup/*bam 
```

### Base recalibration.

https://www.broadinstitute.org/gatk/guide/article?id=2801

Analyse covariation
```
    mkdir -p recal
    cd recal
    parallel  "java -Xmx32g -jar ${GATK} -T BaseRecalibrator -R ${REFERENCE} -I {}  -o {/.}.table  -log {/.}.base_recal.log \ 
               -L references/resequencing_regions.bed" ::: ../realign/*.bam
```

Second pass

```
    parallel  "java jar ${GATK} -T BaseRecalibrator -R ${REFERENCE} -I {} -L  -o {/.}.post.table -BQSR {/.}.table \ 
                --knownSites references/dbsnp.vcf -knownSites references/mills.vcf -knownSites references/1000G_phase1.indels.b37.vcf \
                -log {/.}.second_base_recal.log -L references/resequencing_regions.bed" ::: ../realign/*.bam 
```

Plot generation

```
parallel  "java -Xmx32g -jar ${GATK} -T AnalyzeCovariates -R ${REFERENCE}  -plots {/.}.pdf \ 
           -before {/.}.table -after {/.}.post.table -log {/.}.plots " :::  ../realign/*.bam
```

Print reads    

```
parallel  "java -Xmx32g -jar ${GATK} -T PrintReads -R ${REFERENCE} -I {}  \
           -BQSR {/.}.table -o {/.}.recal_reads.bam" ::: ../realign/*.bam
cd .. 
```

### HaplotypeCaller 

https://www.broadinstitute.org/gatk/guide/article?id=3893

This analysis was actually performed on the NeSI cluster using (but I have added illustrative commands below)

Create Genotype GVCFs 

```
mkdir -p gvcfs
parallel  "java -Xmx32g -jar ${GATK} -T HaplotypeCaller -R ${REFERENCE}  --dbsnp references/dbsnp.vcf -I {} -o gvcfs/{/.}.gvcf \
              --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 \ 
           -L references/resequencing_regions.bed" ::: recal/*.recal_reads.bam
```

Run the HaplotypeCaller

```
    java -Xmx32g -jar ${GATK} -T GenotypeGVCFs  -R ${REFERENCE}  --dbsnp references/dbsnp.vcf \ 
    --max_alternate_alleles 20   \
    `for gvcf in *.gvcf
    do
        echo "--variant ${gvcf}"
    done` \
    -o raw.vcf -L -L references/resequencing_regions.bed 
```

### Variant Filtering

- VQSR (https://www.broadinstitute.org/gatk/guide/article?id=2805)
- Hard filters (https://www.broadinstitute.org/gatk/guide/article?id=2806)


```
    ./scripts/vqsr.sh raw.vcf 
```


### Reheader Variant Files.

Use bcftools to replace the headers in the variant files, because they were nonsense.

```
    bcftools reheader --samples sample_replace.txt recalibrated_snps_raw.vcf > \
    vqsr_snp_names.vcf
    bcftools reheader --samples sample_replace.txt filtered_indels.vcf > \
    filtered_indels_names.vcf
```

### Final output files

Running this analysis outputs two final VCF files, for further processing.

- ```filtered_indels.vcf``` - filtered indels. (not used in any further analyses)
- ```recalibrated_snps_raw.vcf``` - recalibrated SNPs.

