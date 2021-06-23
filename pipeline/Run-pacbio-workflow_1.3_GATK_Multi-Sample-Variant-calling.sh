#!/bin/bash

#####################
# Script for Variant (SNP and Indels) calling using GATK workflow
#####################

######################
# GATK Variant Calling workflows
#https://www.pacb.com/wp-content/uploads/Rowell-CSHLBioData-2018-Comprehensive-Variant-Detection-in-a-Human-Genome-with-PacBio-High-Fidelity-Reads.pdf
#https://www.pacb.com/wp-content/uploads/Application-Brief-Variant-detection-using-whole-genome-sequencing-with-HiFi-reads-Best-Practices.pdf
#https://gatkforums.broadinstitute.org/gatk/discussion/13060/resequencing-using-pacbio-using-gatk-haplotypecaller-to-call-variants

######################
#### create index for the reference using picard
#
java -jar /bioinfo/picard/picard.jar CreateSequenceDictionary \
R=parascaris_crop-seven-genes.fa \
O=parascaris_crop-seven-genes.dict

#### Setting paths

dir_aln="/home/gala0002/proj/proj_Parascaris/DATA_2.0_CCS_map/"
dir_gatk="/home/gala0002/proj/proj_Parascaris/DATA_4.0_GATK/"
mkdir -p $dir_gatk

ref="/home/gala0002/proj/proj_Parascaris/Ref_parascaris_univalens/parascaris_crop-seven-genes.fa"
#GATK="/bioinfo/GATK/gatk-4.1.6.0/gatk"
GATK="/usr/lib/jvm/java-8-openjdk-amd64/bin/java -Xmx200g -jar /bioinfo/GATK/gatk-4.1.6.0/gatk-package-4.1.6.0-local.jar"

######################
# Processing independent samples

for f in `ls /home/gala0002/proj/proj_Parascaris/DATA_2.0_CCS_map/bc100*.bam`; do

BASEFL_CP=$(basename $f | cut -d. -f1 )
sampleID=${BASEFL_CP}

echo "Sample_ID:" ${sampleID}

#conda env create -f gatkcondaenv.yml

#Removing (marking) duplicates with GATK4
######################
cd ${dir_gatk}

${GATK} MarkDuplicatesSpark \
-I ${dir_aln}${sampleID}.bam \
-O ${dir_gatk}${sampleID}.sorted.dedup.bam \
-M ${dir_gatk}${sampleID}.sorted.dedup.metrics.txt

#From BAM files to population variants
######################

####
${GATK} HaplotypeCaller \
-R ${ref} \
-I ${sampleID}.sorted.dedup.bam \
--native-pair-hmm-threads 55 \
--minimum-mapping-quality 60 \
-ERC GVCF \
--pcr-indel-model AGGRESSIVE \
--read-filter MappingQualityReadFilter \
--read-filter NotSecondaryAlignmentReadFilter \
--read-filter NotSupplementaryAlignmentReadFilter \
-O ${sampleID}.g.vcf

done

######################
##### All samples ####
#Followed by combining g.vcf files
######################

cd ${dir_gatk}

${GATK} CombineGVCFs \
-R ${ref} \
--variant bc1001--bc1093.g.vcf --variant bc1001--bc1094.g.vcf --variant bc1001--bc1095.g.vcf --variant bc1001--bc1096.g.vcf --variant bc1002--bc1093.g.vcf --variant bc1002--bc1094.g.vcf --variant bc1002--bc1095.g.vcf --variant bc1002--bc1096.g.vcf --variant bc1003--bc1093.g.vcf --variant bc1003--bc1094.g.vcf --variant bc1003--bc1095.g.vcf --variant bc1003--bc1096.g.vcf --variant bc1004--bc1093.g.vcf --variant bc1004--bc1094.g.vcf --variant bc1004--bc1095.g.vcf --variant bc1004--bc1096.g.vcf --variant bc1005--bc1093.g.vcf --variant bc1005--bc1094.g.vcf --variant bc1005--bc1095.g.vcf --variant bc1005--bc1096.g.vcf --variant bc1006--bc1093.g.vcf --variant bc1006--bc1094.g.vcf --variant bc1006--bc1095.g.vcf --variant bc1006--bc1096.g.vcf --variant bc1007--bc1093.g.vcf --variant bc1007--bc1094.g.vcf --variant bc1007--bc1095.g.vcf --variant bc1007--bc1096.g.vcf \
-O allsample.g.vcf

#joint variant calling with GenotypeGVCFs
######################

${GATK} GenotypeGVCFs \
-R ${ref} \
-V allsample.g.vcf \
-stand-call-conf 5 \
-O allsample.vcf

######################

#Useful tool: VariantFiltration – hard filtering on various criteria Example:

${GATK} VariantFiltration \
-R ${ref} \
--filter-expression "MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)" \
--filter-expression "FS >= 10.0" \
--filter-expression "AN >= 4" \
--filter-expression "DP > 100 || DP < 4" \
--filter-name HARD_TO_VALIDATE \
--filter-name SNPSBFilter \
--filter-name SNPNalleleFilter \
--filter-name SNPDPFilter \
-cluster 3 \
-window 10 \
-V allsample.vcf \
-O allsample.filtered.vcf

######################

echo "done script..."





