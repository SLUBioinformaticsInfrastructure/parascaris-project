#!/bin/bash -l
module load bioinfo-tools trinity trimmomatic blast

in=$(realpath ../data/IVM_P_univalens)
long=$(realpath ../data/ps_329/ccsreads/ps_329_pool1.fasta)
out=$(realpath ../data/trinity)
proj=snic2021-22-60
mail=nicolas.delhomme@umu.se

[[ ! -d $out ]] && mkdir -p $out

sbatch -A $proj --mail-user=$mail -t 4-00:00:00 \
-o $out/trinity.out -e $out/trinity.err \
$(realpath ../UPSCb-common/pipeline/runTrinity.sh) \
-l $long -t -s RF -m 120G $out \
$(find $in -name "*R1_001.fastq.gz" | sort | paste -s --delimiters=,) \
$(find $in -name "*R2_001.fastq.gz" | sort | paste -s --delimiters=,)
