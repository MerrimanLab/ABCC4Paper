#!/usr/bin/bash
#
# 
VCF_INPUT=$1
REFERENCE=references/human_g1k_v37_decoy.fasta
java -jar $GATK \
    -T VariantRecalibrator \
    -R $REFERENCE \ 
    -input $VCF_INPUT    \
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 references/hapmap.vcf \
    -resource:omni,known=false,training=true,truth=true,prior=12.0 references/omni.vcf \
    -resource:1000G,known=false,training=true,truth=false,prior=10.0 references/1000G.vcf \
    -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 references/dbsnp.vcf \
    -an QD \
    -an FS \
    -an MQRankSum \
    -an ReadPosRankSum \
    -mode SNP \
    -nt 6 \
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -rscriptFile recalibrate_SNP_plots.R

java -jar $GATK \ 
    -T ApplyRecalibration \
    -R $REFERENCE \ 
    -input output_raw_vcf.vcf \
    -mode SNP \
    --ts_filter_level 99.0 \
    -recalFile recalibrate_SNP.recal \
    -tranchesFile recalibrate_SNP.tranches \
    -o recalibrated_snps_raw.vcf \
    -nt 6

java -jar $GATK \
    -T SelectVariants \
    -R $REFERENCE \
    -V recalibrated_snps_raw.vcf \
    -selectType INDEL \
    -o raw_indels.vcf


### Filter USing fisher strand and ReadPosRank Sum and Quality by depth. ### 
java -jar $GATK \ 
    -T VariantFiltration \
    -R $REFERENCE 
    -V raw_indels.vcf \
    --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
    --filterName "my_indel_filter" \
    -o filtered_indels.vcf
