# parascaris-project
Eva Tyden and Frida Martin horse's worm parasite - parascaris project

## Setup

* UPPMAX

Compute | Storage
---|---
SNIC 2021/22-60| SNIC 2020/16-13

* CLI
```{bash}
ln -s /proj/snic2020-16-13 data
ln -s /proj/snic2020-16-13/reference .
```

## Analysis

### _de-novo_ transcriptome assembly using the PacBio long reads

Data from Martis et al. was used for a _de-novo_ transcriptome assembly using Trinity version 2.11 (Grabherr, Haas, Yassour et al.,2011 - DOI: 10.1038/nbt.1883) including the PacBio beta-tubulin long-reads as extra support. Details about the parameters are available from the companion code repository.

Beta-tubulin transcripts identification and curation:

The transcriptome obtained at the previous stage was turned into a BLAST database using blast+ version 2.6.0+ (Camacho et al., 2009, - DOI:10.1186/1471-2105-10-421). Beta-tubulin mRNA sequences for the _Parascaris univalens_ beta-tubulin genes were retrieved from the WormBase ParaSite website (https://parasite.wormbase.org/index.html, Howe et al., 2017 - DOI: 10.1016/j.molbiopara.2016.11.005), release 15, Oct 2020 and aligned against the transcriptome database to identify beta-tubulin transcripts. The results were processed wrapping up code published in the UPSCb-common repository (DOI: 10.5281/zenodo.4001774), as described in this manuscript companion code repository. The protein sequences of the transcripts identified as putative beta-tubulin were obtained using Transdecoder v5.5.0 (https://github.com/TransDecoder/TransDecoder) and aligned using muscle on phylogeny.fr. The resulting alignments were manually curated to refine the existing protein structure retrieved from the ParaSite Website, prior to their use in the next step.

### Structural Variant (SV) Calling and Variant (SNP and Indel) calling

* SV-calling
```{bash}
./Run-pacbio-workflow_1.2_SMRT_Mapping-SV-calling.sh
```

* Variant-calling
```{bash}
./Run-pacbio-workflow_1.3_GATK_Multi-Sample-Variant-calling.sh
```

### Phylogenetic analysis

Orthologues protein sequences from the _Parascaris univalens_ beta-tubulin genes were retrieved from the WormBase ParaSite website (https://parasite.wormbase.org/index.html, Howe et al., 2017 - DOI: 10.1016/j.molbiopara.2016.11.005), release 15, Oct 2020. 
The orthologues sequences were kept for the species of interest provided they had at least a reciprocal percentage identity over 60% (the "Target %id” and "Query %id” reported on the website.). 

The sequences were further processed using the phylogeny.fr website using default settings apart from gblocks were the most lenient option was selected.

Finally, the R package ggtree (version 2.4.1, Yu et al., 2017 - DOI:  10.1111/2041-210X.12628) was used to visualise the reconstructed tree.
