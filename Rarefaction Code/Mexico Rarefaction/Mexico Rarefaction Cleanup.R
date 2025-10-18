#Loading Packages
library('tidyverse')
library('phyloseq')
library('ape')
library('vegan')

#Loading Files
#meta file
metaFP <- "././Mexico Dataset QIIME2 files/mexico_metadata.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

#otu file
otuFP <- "././Mexico Dataset QIIME2 files/mexico-feature-table.txt"
otu <- read.delim(file=otuFP,sep = "\t", skip = 1)

#taxonomy table
taxFP <- "././Mexico Dataset QIIME2 files/mexico-taxonomy.tsv"
tax <- read.delim(file=taxFP,sep = "\t")

#phylogenetic tree
phyFP <- "././Mexico Dataset QIIME2 files/mexico-tree.nwk"
phy <- read.tree(phyFP)