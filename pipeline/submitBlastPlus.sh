#!/bin/bash -l
module load bioinfo-tools blast

fasta=$(realpath ../reference/fasta/parascaris_univalens.PRJNA386823.WBPS15.Btubulin_transcripts.fa)
inx=$(realpath ../reference/indices/blast/Trinity.fasta)
out=$(realpath ../data)/blast
proj=snic2021-22-60
mail=nicolas.delhomme@umu.se

[[ ! -d $out ]] && mkdir -p $out

sbatch -A $proj --mail-user=$mail --mail-type=ALL \
-o $out/blast.out -e $out/blast.err -p core -n 8 -t 1:00:00 \
$(realpath ../UPSCb-common/pipeline/runBlastPlus.sh) -p 8 blastn $fasta $inx $out
