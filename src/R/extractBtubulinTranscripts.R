# The sequence of all transcripts where retrieved from Wormbase ParaSite

library(Biostrings)
ref <- readDNAStringSet("parascaris_univalens.PRJNA386823.WBPS15.mRNA_transcripts.fa.gz")
IDs <- c("PgR003_g161","PgE153_g002","PgR007_g022","PgB10_g062","PgR045_g070","PgB04_g135","PgB04_g136")
writeXStringSet(ref[unlist(sapply(IDs,grep,names(ref)))],file="parascaris_univalens.PRJNA386823.WBPS15.Btubulin_transcripts.fa")

