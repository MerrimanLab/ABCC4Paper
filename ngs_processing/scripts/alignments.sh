#!/Volumes/BiochemXsan/staff_users/jamesboocock/bin/parallel --shebang-wrap bash
#
# Script for running bwa mem.

BAM_OUTPUT=$(basename "${1%.fq}.bam")
RG=$(samtools view -H ${BAM_OUTPUT} | grep @RG -m1)
RG=${RG//	/\\t}
echo "$RG"
SAM_OUTPUT=$(basename "${1%.fq}.sam")
bwa mem  -p -R "${RG}" -M ../decoy/human_g1k_v37_decoy.fasta ${1} > bwa_mem/$SAM_OUTPUT
