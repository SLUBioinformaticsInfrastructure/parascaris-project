#' ---
#' title: "B-tubulin de-novo assembly blast"
#' author: "Nicolas Delhomme"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#'    code_folding: hide
#' ---
#' # Setup
#' * Library
suppressPackageStartupMessages({
  library(Biostrings)
  library(here)
  library(igraph)
  library(pander)
})

#' * Helper files
suppressPackageStartupMessages({
  source(here("UPSCb-common/src/R/blastUtilities.R"))
})

#' # Data
blt <- readBlast(here("data/blast/Trinity.fasta_parascaris_univalens.PRJNA386823.WBPS15.Btubulin_transcripts.blt"),
                 format=BM8ext)

df <- blt$df[order(blt$df$query.cum.cov,decreasing=TRUE),]

#' # Analysis
#' Find the minimum coverage of any of our queries (the tubulin genes)
min <- min(df[match(unique(df$query.id),df$query.id),"query.cum.cov"])

#' Create a graph from these
graf <- graph.edgelist(as.matrix(df[df$query.cum.cov > min, c("query.id","subject.id")]))

#' There are six clusters, very promising
cl <- clusters(graf)
cl$no

#' Let's take a visual look
barplot(cl$csize,main="number of transcripts per cluster")
plot.igraph(graf)

#' Looks definitely promising, let's zoom on the individual clusters.
#' 
#' Apart from cluster 1 grouping two beta-tubulin genes together, all the other ones group a single gene isoform.
dev.null <- sapply(1:cl$no,function(n,g){
  message(paste("Looking at cluster",n))
  vertices <- sort(names(cl$membership[cl$membership==n]))
  pander(vertices)
  plot.igraph(induced_subgraph(g,vertices))
},graf)

#' # Export 
#' We export the sequence so that we can subject them to a multiple sequence alignment
seq <- c(readDNAStringSet(here("data/trinity/Trinity.fasta.gz")),
         readDNAStringSet(here("reference/fasta/parascaris_univalens.PRJNA386823.WBPS15.Btubulin_transcripts.fa.gz")))

names(seq) <- sub(" .*","",names(seq))

dir.create(here("data/analysis/muscle"),showWarnings=FALSE,recursive=TRUE)

dev.null <- sapply(1:cl$no,function(n,m,s){
  nam <- sort(names(m[m==n]))
  sub("_t.*","",nam[grep("Pg.*",nam)[1]])
  writeXStringSet(s[nam],file=here("data/analysis/muscle",sub("_t.*",".fa",nam[grep("Pg.*",nam)[1]])))
},cl$membership,seq)

#' # Session Info
#' ```{r session info, echo=TRUE}
#' sessionInfo()
#' ```
