#!/bin/bash -l
module load bioinfo-tools  blast

in=$(realpath ../data/trinity/Trinity.fasta)
out=$(realpath ../reference/indices/blast)
proj=snic2021-22-60
mail=nicolas.delhomme@umu.se

[[ ! -d $out ]] && mkdir -p $out

sbatch -A $proj --mail-user=$mail \
-o $out/makeblastdb.out -e $out/makeblastdb.err \
$(realpath ../UPSCb-common/pipeline/runBlastPlusMakeblastdb.sh) -t trinity $in $out
