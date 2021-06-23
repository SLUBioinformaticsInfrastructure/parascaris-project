suppressPackageStartupMessages({
  library(Biostrings)
  library(here)
  library(tidyverse)
})

species <- scan(here("doc/species.txt"),what="character",sep="\n")

#' Some IDs in the 
tx2gene <- data.frame(TX=c("WBGene00228922Bma-mec-7",
                           "WBGene00003171mec-7",
                           "WBGene00238644Ovo-mec-7",
                           "WBGene00224994Bm4733",
                           "WBGene00247048Ovo-ben-1",
                           "WBGene00233027Bma-tbb-4",
                           "WBGene00006538tbb-4",
                           "WBGene00238643Ovo-tbb-4",
                           "WBGene00248365OVOC11556",
                           "Tsp_04895"),
                      GENEID=c("Bm8661",
                               "ZK154.3",
                               "OVOC1835",
                               "Bm4733",
                               "OVOC10239",
                               "Bm12766",
                               "B0272.1",
                               "OVOC1834",
                               "OVOC11556",
                               "EFV50889"
                      ))

a.galli.seq <- readAAStringSet(here("data/Martis_et_al_2017/selected_genes_TBB_pep.fa"))
names(a.galli.seq) <- sub("TRINITY","Agalli",sub("::.*","",names(a.galli.seq)))

p.univalens <- readAAStringSet(here("data/analysis/transdecoder/all.fna.transdecoder.pep"))
names(p.univalens) <- sub(" .*","",names(p.univalens))

# filter transcripts (this was done after an initial tree was constructed)
sel <- grep("TRINITY_DN28139_c0.*",names(p.univalens))
p.univalens <- p.univalens[ - sel[!grepl("TRINITY_DN28139_c0_g1_i45.p1",names(p.univalens)[sel])] ]

sel <- grep("TRINITY_DN1498_c0.*",names(p.univalens))
p.univalens <- p.univalens[ - sel[!grepl("TRINITY_DN1498_c0_g1_i11.p1",names(p.univalens)[sel])] ]

sel <- grep("PgR003_g161*",names(p.univalens))
p.univalens <- p.univalens[ - sel[!grepl("PgR003_g161_t01.p1",names(p.univalens)[sel])] ]

sel <- grep("PgB10_g062*",names(p.univalens))
p.univalens <- p.univalens[ - sel[!grepl("PgB10_g062_t01.p1",names(p.univalens)[sel])] ]

sel <- grep("TRINITY_DN11305_c0.*",names(p.univalens))
p.univalens <- p.univalens[ - sel[!grepl("TRINITY_DN11305_c0_g1_i10.p1",names(p.univalens)[sel])] ]

names(p.univalens) <- sub("\\.p1$","",names(p.univalens))

# outgroup
d.melanogaster <- readAAStringSet(here("data/ncbi/Dmelanogaster-b-tubulin.faa"))
names(d.melanogaster) <- sub(" .*","",names(d.melanogaster))


orth <- Reduce(bind_rows,lapply(dir(here("data/Wormbase-Parasite"),pattern="*.csv",full.names=TRUE),function(f){
  read_csv(f,
         col_names=c("Species","ID","TargetPercID","QueryPercID"),
         col_types=cols(.default=col_character())) %>% 
  filter(TargetPercID > 60 & QueryPercID > 60) %>% 
  mutate(Species=sub(" \\(.*","",Species)) %>% 
  filter(Species %in% species)
})) %>% distinct(ID,.keep_all=TRUE)

aa <- Reduce("c",lapply(dir(here("data/Wormbase-Parasite"),pattern="*faa.gz",full.names=TRUE),readAAStringSet))
aa <- aa[!duplicated(names(aa))]
names(aa) <- sub("-mRNA-1$|\\.1$|\\.t1$","",names(aa))

sel <- sapply(tx2gene$TX,grep,orth$ID)
orth$ID[sel] <- tx2gene$GENEID
lst <- sapply(sub("Tubu$|Tu$|TBB $|Tub$|Tubuli$|tbb-$|t$|Tubul$|No de?s?$|No$","",substr(orth$ID,1,14)),grep,names(aa))

# if that fires, need to extend the tx2gene
stopifnot(typeof(lst)=="integer")

orth$ID <- names(lst)
orth$Sequence <- as.character(aa[lst])
orth <- bind_rows(orth,tibble(Species="Ascaridia galli",ID=names(a.galli.seq),TargetPercID=NA,
                      QueryPercID=NA,Sequence=as.character(a.galli.seq)),
                  tibble(Species="Parascaris univalens",ID=names(p.univalens),TargetPercID=NA,
                         QueryPercID=NA,Sequence=as.character(p.univalens)),
                  tibble(Species="Drosophila melanogaster",ID=names(d.melanogaster),TargetPercID=NA,
                         QueryPercID=NA,Sequence=as.character(d.melanogaster))
                  )

writeXStringSet(c(aa[lst],a.galli.seq,
                  p.univalens, d.melanogaster),filepath=here("data/phylogeny.fr/orthologues.faa"))

write_delim(orth,file=here("data/phylogeny.fr/orthologues.tsv"))


