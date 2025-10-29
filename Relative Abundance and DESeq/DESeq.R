#Loading Packages
library('tidyverse')
library('phyloseq')
library('DESeq2')
library('ape')
library('vegan')

#Loading Files
#meta file
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

#otu file
otuFP <- "Mexico Dataset QIIME2 files/mexico-feature-table.txt"
otu <- read.delim(file=otuFP,sep = "\t", skip = 1)

#taxonomy table
taxFP <- "Mexico Dataset QIIME2 files/mexico-taxonomy.tsv"
tax <- read.delim(file=taxFP,sep = "\t")

#phylogenetic tree
phyFP <- "Mexico Dataset QIIME2 files/mexico-tree.nwk"
phy <- read.tree(phyFP)

#Loading Phyloseq object
load("Rarefaction Code/cm_final.RData")

#loading Rareified data
load("Rarefaction Code/cm_rare.RData")


#conversion to DESeq object
phyloseq_object_plus1 <- transform_sample_counts(cm_final, function(x) x+1)

deseq_object <- phyloseq_to_deseq2(phyloseq_object_plus1, ~sex)

  
DESeq_output <- DESeq(deseq_object)



results(DESeq_output, tidy=TRUE)






