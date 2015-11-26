#!/usr/bin/bash
#
# Script for running base recallibration.

parallel -j 8 "java -Xmx16g -jar ../GenomeAnalysisTK.jar -T  RealignerTargetCreator -R ../references/GRCh37-lite.fa -I {} -o {/.}.intervals --known ../references/ALL.wg
s.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.indels.sites.vcf.gz --known ../references/ALL.wgs.low_coverage_vqsr.20101123.indels.sites.vcf.gz" ::: 
am
