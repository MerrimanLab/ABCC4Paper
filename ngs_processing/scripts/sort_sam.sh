#!/Volumes/BiochemXsan/staff_users/jamesboocock/bin/parallel --shebang-wrap bash
BAM_OUTPUT="${1%.sam}.bam"
java -jar $PICARD -Xmx4g SortSam INPUT=$1 OUTPUT=$BAM_OUTPUT SORT_ORDER=coordinate
