#Loading Packages
library('tidyverse')
library('phyloseq')
library('ape')
library('vegan')

#Loading Files
#meta file
metaFP <- "././China Dataset QIIME2 files/china_metadata.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

#otu file
otuFP <- "././China Dataset QIIME2 files/china-feature-table.txt"
otu <- read.delim(file=otuFP,sep = "\t", skip = 1)

#taxonomy table
taxFP <- "././China Dataset QIIME2 files/china-taxonomy.tsv"
tax <- read.delim(file=taxFP,sep = "\t")

#phylogenetic tree
phyFP <- "././China Dataset QIIME2 files/china-tree.nwk"
phy <- read.tree(phyFP)