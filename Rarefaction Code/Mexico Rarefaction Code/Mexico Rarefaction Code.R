#Loading Packages
library('tidyverse')
library('phyloseq')
library('ape')
library('vegan')

#Loading Files
#meta file
metaFP <- "../../Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

#otu file
otuFP <- "../../Mexico Dataset QIIME2 files/mexico-feature-table.txt"
otu <- read.delim(file=otuFP,sep = "\t", skip = 1)

#taxonomy table
taxFP <- "../../Mexico Dataset QIIME2 files/mexico-taxonomy.tsv"
tax <- read.delim(file=taxFP,sep = "\t")

#phylogenetic tree
phyFP <- "../../Mexico Dataset QIIME2 files/mexico-tree.nwk"
phy <- read.tree(phyFP)

#phyloseq data prep
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$`X.OTU.ID`
OTU <- otu_table(otu_mat,taxa_are_rows = TRUE)

meta_df <- as.data.frame(meta[,-1])
rownames(meta_df) <- meta$`sample.id`
META <- sample_data(meta_df)

tax_mat <- tax |>
  select(-Confidence) |>
  separate(col = Taxon, sep=';',
           into = c("Domain","Phylum",'Class','Order',
                    'Family','Genus','Species')) |>
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$`Feature.ID`
TAX <- tax_table(tax_mat)

#-----------------------------------------------------
#making phyloseq combination + verification + filtering / rarification
cm <- phyloseq(OTU,META,TAX,phy)
#-------
cm_filt <-subset_taxa(cm, Domain == "d__Bacteria" & Class!="c__Chloroplast" & Class!="f__Mitochondira")
cm_filt_nolow <- filter_taxa(cm_filt, function(x) sum(x)>5, prune = TRUE)
cm_filt_nolow_samps <- prune_samples(sample_sums(cm_filt_nolow)>100,cm_filt_nolow)
cm_final <- subset_samples(cm_filt_nolow_samps,!is.na(sex))

#-------
rarecurve(t(as.data.frame(otu_table(cm_final))),cex=0.1)
rare_cm <- rarefy_even_depth(cm_final,rngseed = 11,sample.size = 3000)

#-------
save(cm_final,file="cm_final.RData")
save(cm_rare,file="cm_rare.RData")
#-----------------------------------------------------