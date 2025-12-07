# --- 1. Load Libraries ---
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(microbiome)
library(indicspecies)
library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)

setwd("/Users/lunakurokawa/Desktop/Group_Project_11")

# --- 2. Load Data ---
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(metaFP, sep = "\t")

otuFP <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
otu <- read.delim(otuFP, sep = "\t", skip = 1)

taxFP <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
tax <- read.delim(taxFP, sep = "\t")

phyFP <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"
phy <- read.tree(phyFP)

# --- 3. Phyloseq Prep ---
# OTU table
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$`X.OTU.ID`
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

# Sample metadata
meta_df <- meta %>% column_to_rownames("sample.id")
META <- sample_data(meta_df)

# Taxonomy table
tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge") %>%
  as.data.frame(stringsAsFactors = FALSE)

tax_mat2 <- as.matrix(tax_mat[,-1])
rownames(tax_mat2) <- tax$Feature.ID

# Construct phyloseq object
TAX <- tax_table(tax_mat2)
ps <- phyloseq(OTU, TAX, META, phy)

# --- 4. Define sex-severity groups ---
meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
meta$Severity <- NA_character_
meta$Severity[grepl("asymptomatic", tolower(meta$Group))] <- "Control"
meta$Severity[grepl("ambulatory positive", tolower(meta$Group))] <- "Mild"
meta$Severity[grepl("hospitalized positive", tolower(meta$Group))] <- "Severe"
meta$Severity[grepl("deceased hospitalized", tolower(meta$Group))] <- "Deceased"

meta$Group_combined <- paste(meta$sex, meta$Severity, sep = "_")
sample_data(ps) <- sample_data(meta)

# --- 5. Collapse to Genus Level ---
ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = TRUE)

# --- 6. Prepare data for ISA ---
otu_mat <- as.data.frame(t(otu_table(ps_genus)))
meta_df <- data.frame(sample_data(ps_genus), stringsAsFactors = FALSE)

# Keep only samples with non-NA Group_combined
meta_clean <- meta_df[!is.na(meta_df$Group_combined), ]
otu_clean <- otu_mat[rownames(meta_clean), ]
otu_clean[is.na(otu_clean)] <- 0

# --- 7. Run Indicator Species Analysis ---
isa_mpt <- multipatt(
  otu_clean,
  meta_clean$Group_combined,
  func = "IndVal.g",
  duleg = TRUE,
  control = how(nperm = 999)
)

summary(isa_mpt, indvalcomp = TRUE, alpha = 0.05)

# --- 8. Extract Significant Results (p <= 0.05) ---
sig_table <- data.frame(isa_mpt$sign) %>%
  rownames_to_column("OTU") %>%
  filter(p.value <= 0.05)

group_names <- colnames(isa_mpt$comb)
sig_table$cluster <- group_names[sig_table$index]

significant_indicators <- sig_table %>%
  select(OTU, stat, p.value, cluster)

# --- 9. Annotate with Genus Names ---
tax_df <- as.data.frame(tax_table(ps_genus)) %>%
  rownames_to_column("OTU") %>%
  mutate(Label = Genus) %>% 
  filter(!is.na(Label) & Label != "" & Label != "g__")  # remove unclassified/missing

top_sps_named <- significant_indicators %>%
  left_join(tax_df %>% select(OTU, Label), by = "OTU") %>%
  filter(!grepl("NA", cluster))


# --- 10. Plot ---
custom_colors <- c(
  "male_Control" = "cyan",
  "female_Control" = "mediumpurple",
  "male_Mild" = "steelblue",
  "female_Mild" = "#d62728",
  "male_Severe" = "mediumturquoise",
  "female_Severe" = "#e377c2",
  "male_Deceased" = "springgreen",
  "female_Deceased" = "orchid1"
)

isa_plot <- ggplot(top_sps_named, aes(x = reorder(Label, stat), y = stat, fill = cluster)) +
  geom_col() +
  coord_flip() +
  scale_y_continuous(limits = c(0, 1.0)) +
  scale_fill_manual(values = custom_colors) +  
  labs(
    x = "Genus",
    y = "Indicator Value",
    fill = "Sex-Severity Group",
    title = "Significant Indicator Genera"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.minor.y = element_blank()
  )

print(isa_plot)
ggsave("indicator_species_genus_plot_custom_colors.png",
       plot = isa_plot, width = 10, height = 15, dpi = 300)
