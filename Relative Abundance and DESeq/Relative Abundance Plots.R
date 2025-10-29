#Loading Packages
library('tidyverse')
library('phyloseq')
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

##-------------------------------------------------
## Adapted from https://www.youtube.com/watch?v=siIoupAnILk

## Adapting Taxon data
tax_df <- tax |>
  select(-Confidence) |>
  separate(col = Taxon, sep=';',
           into = c("Domain","Phylum",'Class','Order',
                    'Family','Genus','Species')) |>
  as.data.frame()

## Adapting Meta data
 meta_adapted <- meta |>
   rename(sample_id = sample.id)

## Building main data
dat <- otu |>
  pivot_longer(-X.OTU.ID,names_to = "sample_id", values_to = "count") |>
  rename(Feature.ID = X.OTU.ID) |>
  left_join(tax_df,by = "Feature.ID") |>
  left_join(meta_adapted,by = "sample_id") |>
  filter(!is.na(sex)) |>
  filter(!is.na(Phylum))

## Constructing the plot

rawabundanceplot <- ggplot(data = dat, aes(x=sample_id,y=count)) +
  facet_grid(~sex,space = "free_x", scale = "free_x") +
  geom_bar(aes(fill = Phylum),stat = "identity",position="fill")

rawabundanceplot

##----------------
#Loading Phyloseq object
load("Rarefaction Code/cm_final.RData")

#loading Rareified data
load("Rarefaction Code/cm_rare.RData")


cm_RA <- transform_sample_counts(rare_cm, function(x) x/sum(x))

cm_phylum <- tax_glom(cm_RA,taxrank = "Phylum", NArm = FALSE)

phyloseq_ra_plot <- plot_bar(cm_phylum,fill = "Phylum") + 
  facet_wrap(~sex, scales = "free_x")
phyloseq_ra_plot

#-------------
ggsave(filename ="Relative Abundance and DESeq/rawabundanceplot.png",rawabundanceplot)
ggsave(filename ="Relative Abundance and DESeq/phyloseq_ra_plot.png",phyloseq_ra_plot)
