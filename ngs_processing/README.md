## Processing of NGS data from wustl for the resequencing experiment.

Sequencing data was provided in SAM format, this first step was to convert all
the SAMs back into fastq files. 

http://gatkforums.broadinstitute.org/discussion/2908/howto-revert-a-bam-file-to-fastq-format

## Software Requirements

- GATK (v3.3)
- samtools (v1.1-2)
- bcftools (v1.1-6)
- bwa (v0.7.12-r1039)
- Picard tools (v1.41) 


## Resources 

Broad institute bundle contains reference genome and auxilary files used in the GATK best-practices

Best practices - https://www.broadinstitute.org/gatk/guide/best-practices 

- ftp://ftp.broadinstitute.org/bundle/2.8/b37/

## Pipeline

### Bash Variables

Need to set the following variables to your environment.

```bash
    export PICARD=<path to picard.jar>
    export GATK=<path to GATK.jar>
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

###  alignment

Now we follow the GATK best-practices

All alignments go into the bwa\_mem folder.

```bash
     mkdir -p bwa_mem
     ./scripts/alignments.sh *.fq
```








