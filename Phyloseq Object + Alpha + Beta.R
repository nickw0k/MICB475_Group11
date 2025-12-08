library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(dplyr)
library(picante)

#### Load Data ####
covid_mex_metadata_fp <- "COVID_Mexico_Final/covid_mex_metadata.tsv"
covid_mex_metadata <- read_delim(covid_mex_metadata_fp, delim="\t")

covid_mex_otu_fp <- "COVID_Mexico_Final/covid_mexico_final_export/table_export/feature-table.txt"
covid_mex_otu <- read_delim(covid_mex_otu_fp, delim="\t", skip=1)

covid_mex_tax_fp <- "COVID_Mexico_Final/covid_mexico_final_export/taxonomy_export/taxonomy.tsv"
covid_mex_taxonomy <- read_delim(covid_mex_tax_fp, delim="\t")

covid_mex_phylotree_fp <- "COVID_Mexico_Final/covid_mexico_final_export/rooted_tree_export/tree.nwk"
covid_mex_phylotree <- read.tree(covid_mex_phylotree_fp)

#### Format OTU table ####
# OTU tables should be a matrix with rownames and colnames as OTUs and SampleIDs respectively
# save everything except first column (OTU ID) into a matrix
covid_mex_otu_matrix <- as.matrix(covid_mex_otu[, -1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(covid_mex_otu_matrix) <- covid_mex_otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(covid_mex_otu_matrix, taxa_are_rows = TRUE)
class(OTU)

#### Format sample metadata ####
# Save everything except sample-id as new data frame
covid_mex_metadata_df <- as.data.frame(covid_mex_metadata[, -1])
# Make sample-id the rownames
rownames(covid_mex_metadata_df) <- covid_mex_metadata$'sample-id'
# Changing "Group" column name to "Severity" and standardizing "Group" column into 4 categories
# which are: Asymptomatic, Hospitalized Negative, Hospitalized Positive, and Deceased
covid_mex_metadata_edited <- covid_mex_metadata_df %>%
  mutate(Group = case_when(
    grepl("asymptomatic", Group, ignore.case = TRUE) ~ "Asymptomatic",
    grepl("ambulatory negative", Group, ignore.case = TRUE) ~ "Mild negative",
    grepl("ambulatory positive", Group, ignore.case = TRUE) ~ "Mild",
    grepl("hospitalized positive", Group, ignore.case = TRUE) &
      !grepl("deceased", Group, ignore.case = TRUE) ~ "Severe",
    grepl("deceased", Group, ignore.case = TRUE) ~ "Deceased",
    TRUE ~ "Unknown"   # fallback for unexpected entries
  )) %>%
  rename(severity = Group)  # rename "Group" to "severity"
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(covid_mex_metadata_edited)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
covid_mex_taxonomy_matrix <- covid_mex_taxonomy %>% select(-Confidence)%>%
  separate(col=Taxon, sep=";"
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs
covid_mex_taxonomy_matrix <- covid_mex_taxonomy_matrix[, -1]
# Make sample-ids the rownames
rownames(covid_mex_taxonomy_matrix) <- covid_mex_taxonomy$`Feature ID`
# Make taxa table
TAX <- tax_table(covid_mex_taxonomy_matrix)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
covid_mexico_phyloseq <- phyloseq(OTU, SAMP, TAX, covid_mex_phylotree)

#### Looking at phyloseq object ####
# View components of phyloseq object with the following commands
otu_table(covid_mexico_phyloseq)
sample_data(covid_mexico_phyloseq)
tax_table(covid_mexico_phyloseq)
phy_tree(covid_mexico_phyloseq)

# View sample names in our data set
sample_names(covid_mexico_phyloseq)

# Number of samples in our data set
nsamples(covid_mexico_phyloseq)

# Sum of reads in each sample
sample_sums(covid_mexico_phyloseq)

# 3 samples with the most reads from our data set
samps_most_reads <- names(sort(sample_sums(covid_mexico_phyloseq), decreasing = TRUE)[1:3])
get_taxa(covid_mexico_phyloseq, samps_most_reads)

# Names of taxa from our data set
taxa_names(covid_mexico_phyloseq)

# Number of taxa we have from our data set
ntaxa(covid_mexico_phyloseq)

# 3 most abundant taxa from our data set
taxa_most_abundant <- names(sort(taxa_sums(covid_mexico_phyloseq), decreasing = TRUE)[1:3])
get_sample(covid_mexico_phyloseq, taxa_most_abundant)

# Remove any non-bacterial sequences, if any
covid_mexico_phyloseq_filt <- subset_taxa(
  covid_mexico_phyloseq,
  !(Family %in% c("Mitochondria") |
      Order %in% c("Chloroplast")))

# Remove ASVs that have less than 5 counts total
covid_mexico_phyloseq_filt_nolow <- filter_taxa(covid_mexico_phyloseq_filt,
                                                function(x) sum(x)>5,
                                                prune = TRUE)

# Remove samples with less than 100 reads
covid_mexico_phyloseq_filt_nolow_samps <- prune_samples(sample_sums(covid_mexico_phyloseq_filt_nolow)>100,
                                                        covid_mexico_phyloseq_filt_nolow)

# Remove samples where Sex is NA
covid_mexico_phyloseq_no_na <- subset_samples(covid_mexico_phyloseq_filt_nolow_samps, !is.na(sex))

# Filter out "Mild Negative" group in Severity column
covid_mexico_phyloseq_final <- subset_samples(covid_mexico_phyloseq_no_na, severity != "Mild Negative")

# Rarefy samples
# rngseed sets a random number. To be able to reproduce this exact analysis each time, need to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(covid_mexico_phyloseq_final))), cex=0.1)
covid_mexico_phyloseq_rare <- rarefy_even_depth(covid_mexico_phyloseq_final, rngseed = 11, sample.size = 83384)

##### Saving #####
save(covid_mexico_phyloseq_final, file="covid_mexico_phyloseq_final.RData")
save(covid_mexico_phyloseq_rare, file="covid_mexico_phyloseq_rare.RData")

#### Alpha diversity ###
gg_richness_sex <- plot_richness(covid_mexico_phyloseq_rare, x = "sex", measures = c("Shannon")) +
  xlab("Sex") +
  geom_boxplot()
gg_richness_sex

gg_richness_severity <- plot_richness(covid_mexico_phyloseq_rare, x = "severity", measures = c("Shannon")) +
  xlab("Severity") +
  geom_boxplot()
gg_richness_severity

ggsave(filename = "plot_richness_sex.png", gg_richness_sex, height=4, width=6)
ggsave(filename = "plot_richness_severity.png", gg_richness_severity, height=4, width=6)

estimate_richness(covid_mexico_phyloseq_rare)

# Phylogenetic diversity
# Calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(t(otu_table(covid_mexico_phyloseq_rare)), phy_tree(covid_mexico_phyloseq_rare),
                 include.root=F)

# Add PD to metadata table
sample_data(covid_mexico_phyloseq_rare)$PD <- phylo_dist$PD

# Plot Sex and Severity metadata category against the PD
plot_pd_sex <- ggplot(sample_data(covid_mexico_phyloseq_rare), aes(sex, PD)) +
  geom_boxplot() +
  xlab("Sex") +
  ylab("Phylogenetic Diversity")

plot_pd_severity <- ggplot(sample_data(covid_mexico_phyloseq_rare), aes(severity, PD)) +
  geom_boxplot() +
  xlab("Severity") +
  ylab("Phylogenetic Diversity")

# View plots
plot_pd_sex
plot_pd_severity

ggsave("plot_pd_sex.png", plot_pd_sex, height = 4, width = 5)
ggsave("plot_pd_severity.png", plot_pd_severity, height = 4, width = 5)

#### Beta diversity #####
wunifrac_dm <- distance(covid_mexico_phyloseq_rare, method="wunifrac")

pcoa_bc <- ordinate(covid_mexico_phyloseq_rare, method="PCoA", distance=wunifrac_dm)

plot_ordination(covid_mexico_phyloseq_rare, pcoa_bc, color = "severity", shape="sex")

gg_pcoa <- plot_ordination(covid_mexico_phyloseq_rare, pcoa_bc, color = "severity", shape="sex") +
  labs(pch="Severity", col = "Sex")
gg_pcoa

ggsave("plot_pcoa.png", gg_pcoa, height=4, width=5)

#### Taxonomy bar plots ####
# Plot bar plot of taxonomy
plot_bar(covid_mexico_phyloseq_rare, fill = "Phylum")

# Convert to relative abundance
covid_mexico_phyloseq_RA <- transform_sample_counts(covid_mexico_phyloseq_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
covid_mexico_phylum <- tax_glom(covid_mexico_phyloseq_RA, taxrank = "Phylum", NArm=FALSE)

plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~sex, scales = "free_x")

plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~severity, scales = "free_x")

gg_taxa_sex <- plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~sex, scales = "free_x")
gg_taxa_sex

gg_taxa_severity <- plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~severity, scales = "free_x")
gg_taxa_severity

ggsave("plot_taxonomy_sex.png"
       , gg_taxa_sex
       , height=8, width =12)

ggsave("plot_taxonomy_severity.png"
       , gg_taxa_severity
       , height=8, width =12)