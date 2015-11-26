#!/bin/bash
METRICS_FILE="${1%.bam}.metrics.txt"
PICARD_MARK_OUTPUT=$(basename "${1%.bam}.bam")
java -jar -Xmx4g /Volumes/BiochemXsan/staff_users/jamesboocock/picard/MarkDuplicates.jar INPUT=$1 OUTPUT=rmdup/$PICARD_MARK_OUTPUT METRICS_FILE=$METRICS_FILE 
