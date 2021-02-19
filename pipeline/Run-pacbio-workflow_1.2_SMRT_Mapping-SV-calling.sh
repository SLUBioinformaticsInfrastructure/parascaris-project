#!/bin/bash

#####################
# Script to map ccs reads to the reference genes and SV calling using SMRT tool
#####################

#Install SMRT tool and set memory 
#ulimit -n 65535
#ulimit -c unlimited

################

SMRT="/home/gala0002/softwares/SMRT/smrtlink/install/smrtlink-release_9.0.0.92188/bundles/smrttools/smrtcmds/bin/"
rawdata_sub="/home/gala0002/proj/proj_Parascaris/ps_329/rawdata/ps_329_pool1/r54259_20201218_083912/2_B01/"
rawdata_ccs="/home/gala0002/proj/proj_Parascaris/ps_329/ccsreads/ps_329_pool1/"
dir_sub_ccs="/home/gala0002/proj/proj_Parascaris/DATA_1.0_SUD-CCS_reads/"
dir_ccs_map="/home/gala0002/proj/proj_Parascaris/DATA_2.0_CCS_map/"
dir_pbsv="/home/gala0002/proj/proj_Parascaris/DATA_3.0_pbsv/"

ref="/home/gala0002/proj/proj_Parascaris/Ref_parascaris_univalens/parascaris_crop-seven-genes.fa"
ref2="/home/gala0002/proj/proj_Parascaris/Ref_parascaris_univalens/parascaris_crop-seven-genes.mmi"
ref3="/home/gala0002/proj/proj_Parascaris/Ref_parascaris_univalens/parascaris_crop-seven-genes_refsub.bam"

#### Map with SMRT pbmm2
# A. Generate index file for reference and reuse it to align reads
#${SMRT}pbmm2 index ${ref} 
#${SMRT}pbmm2 align ${ref2} ${rawdata_sub}m54259_201218_140834.subreads.bam ${ref3}

 #B. Align reads and sort on-the-fly, with 4 alignment and 2 sort threads
#${SMRT}pbmm2 align ${ref} ${rawdata_sub}m54259_201218_140834.subreads.bam ${ref3} --sort -j 4 -J 55

#### Map with SMRT 
#pbmm2 align [options] <ref.fa|xml|mmi> <in.bam|xml|fa|fq> [out.aligned.bam|xml]
#ConsensusReadSet > ConsensusAlignmentSet
#pbmm2 align hg38.referenceset.xml movie.consensusreadset.xml hg38.movie.consensusalignmentset.xml --preset CCS

#### one sample
#${SMRT}pbmm2 align ${ref} ${rawdata_ccs}demultiplex.bc1001--bc1093.hifi_reads.fastq ${dir_ccs_map}bc1001--bc1093.bam --preset CCS --sort --median-filter --sample bc1001--bc1093 -j 55

#####################
#### Mapping all independent samples 

for f in `ls /home/gala0002/proj/proj_Parascaris/ps_329/ccsreads/ps_329_pool1/demultiplex.*.fastq`; do
DIR=$(dirname $f)"/"
insample=$(basename $f)
sampleID=$(basename $f | cut -d. -f2)

echo "DIR_file:" ${f}
echo "Sample_filepath:" ${insample}
echo "Sample_ID:" ${sampleID}

${SMRT}pbmm2 align -j 55 --preset CCS --sort --median-filter --sample ${sampleID} ${ref} ${f} ${dir_ccs_map}${sampleID}".bam"

samtools index ${dir_ccs_map}${sampleID}".bam"

${SMRT}pbindex ${dir_ccs_map}${sampleID}".bam"

# ################
# SV calling on independent sample

${SMRT}pbsv discover ${dir_ccs_map}${sampleID}".bam" ${dir_pbsv}${sampleID}".svsig.gz"

done

################
# merging SVs from all samples 

cd ${dir_pbsv}

${SMRT}pbsv call ${ref} bc1001--bc1093.svsig.gz bc1001--bc1094.svsig.gz bc1001--bc1095.svsig.gz bc1001--bc1096.svsig.gz bc1002--bc1093.svsig.gz bc1002--bc1094.svsig.gz bc1002--bc1095.svsig.gz bc1002--bc1096.svsig.gz bc1003--bc1093.svsig.gz bc1003--bc1094.svsig.gz bc1003--bc1095.svsig.gz bc1003--bc1096.svsig.gz bc1004--bc1093.svsig.gz bc1004--bc1094.svsig.gz bc1004--bc1095.svsig.gz bc1004--bc1096.svsig.gz bc1005--bc1093.svsig.gz bc1005--bc1094.svsig.gz bc1005--bc1095.svsig.gz bc1005--bc1096.svsig.gz bc1006--bc1093.svsig.gz bc1006--bc1094.svsig.gz bc1006--bc1095.svsig.gz bc1006--bc1096.svsig.gz bc1007--bc1093.svsig.gz bc1007--bc1094.svsig.gz bc1007--bc1095.svsig.gz bc1007--bc1096.svsig.gz Sample_all_out.vcf

#####################


echo "done script..."

