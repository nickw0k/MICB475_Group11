# --- 1. Load Necessary Libraries ---
# Ensure all these libraries are installed
library(phyloseq)
library(ggplot2)
library(tidyverse)
library(microbiome)
library(indicspecies)
library(DESeq2) # Not used in this script, but in your original
library(ape)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble) # For rownames_to_column

# --- 2. Loading Files ---
# (This is the code you provided)
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

#otu file
otuFP <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
otu <- read.delim(file=otuFP,sep = "\t", skip = 1)

#taxonomy table
taxFP <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
tax <- read.delim(file=taxFP,sep = "\t")

#phylogenetic tree
phyFP <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"
phy <- read.tree(phyFP)

# --- 3. Phyloseq Data Prep ---
# (This is the code you provided)
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


#making phyloseq combination 
ps <- phyloseq(OTU,META,TAX,phy)

# --- 4. Prepare Metadata for Analysis ---
# (This is the code you provided)

# Separate phyloseq object into sex-severity groups

meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
meta$Severity <- NA_character_
meta$Severity[grepl("asymptomatic", tolower(meta$Group))] <- "Control"
meta$Severity[grepl("ambulatory negative", tolower(meta$Group))] <- "Mild-"
meta$Severity[grepl("ambulatory positive",   tolower(meta$Group))] <- "Mild+"
meta$Severity[grepl("hospitalized positive", tolower(meta$Group))] <- "Severe"
meta$Severity[grepl("Deceased hospitalized",   tolower(meta$Group))] <- "Fatal"

sample_data(ps) <- sample_data(meta)

# Create the combined group variable for ISA
meta <- data.frame(sample_data(ps))
meta$Group_combined <- paste(meta$sex, meta$Severity, sep = "_")
sample_data(ps) <- sample_data(meta)


# --- 5. Run Indicator Species Analysis (ISA) ---
# (This is the code you provided)

# Prepare OTU table (samples as rows, taxa as columns)
otu_mat <- as.data.frame(t(otu_table(ps)))

# Prepare metadata
meta_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

# Clean data: remove samples with NA groups
meta_clean <- meta_df[!is.na(meta_df$Group_combined), ]

# Ensure OTU table and metadata have the same samples
otu_clean <- otu_mat[rownames(meta_clean), ]

# Replace any potential NAs in the OTU table with 0
otu_clean[is.na(otu_clean)] <- 0

# Run the multipatt analysis
isa_mpt <- multipatt(otu_clean, meta_clean$Group_combined,
                     func = "IndVal.g", duleg = TRUE,
                     control = how(nperm = 999))

# Optional: Print a summary to the console to check
summary(isa_mpt, indvalcomp = TRUE, alpha = 0.05)


# --- 6. Programmatically Extract Significant Results ---
# (This is the code I provided earlier)
# This step REPLACES all the manual data.frame() blocks

# Get the table of significant results (p.value <= 0.05)
sig_table <- data.frame(isa_mpt$sign) %>%
  rownames_to_column(var = "OTU") %>%
  filter(p.value <= 0.05)

# Get the names of the groups from the analysis
# 'isa_mpt$comb' stores the cluster combinations
group_names <- colnames(isa_mpt$comb)

# Map the significant 'index' (a number) to the actual 'cluster' name
sig_table$cluster <- group_names[sig_table$index]

# Create the final, clean table of significant indicators
significant_indicators <- sig_table %>%
  select(OTU, stat, p.value, cluster)

# You can print the first few rows to see the result
print("--- Top Significant Indicators (Programmatically Extracted) ---")
print(head(significant_indicators))


# --- 7. Annotate Results with Taxonomy (UPDATED) ---

# Create a clean taxonomy data frame from the phyloseq object
tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column(var = "OTU") %>%
  # NEW: Clean up the taxonomy strings before coalescing
  mutate(
    Species = na_if(Species, "s__"), 
    Genus = na_if(Genus, "g__"),   
    Family = na_if(Family, "f__"), 
    Order = na_if(Order, "o__"),   
    Class = na_if(Class, "c__"),  
    Phylum = na_if(Phylum, "p__"), 
    Domain = na_if(Domain, "d__")  
  ) %>%
 
  
  mutate(
    Label = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain),
    Label = ifelse(is.na(Label), "Unclassified", Label)
  )


top_sps_named <- significant_indicators %>%
  left_join(tax_df %>% select(OTU, Label), by = "OTU")


# --- 8. Create the Final Plot ---

isa_plot <- ggplot(top_sps_named, aes(x = reorder(Label, stat), y = stat, fill = cluster)) +
  geom_col() +
  coord_flip() + 
  labs(
    x = "Taxonomy",
    y = "Indicator Value (stat)",
    fill = "Sex-Severity Group",
    title = "Significant Indicator Species by Group"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5) 
  )

# Display the plot
print(isa_plot)

# Save the plot to a file
ggsave("indicator_species_plot_all.png", plot = isa_plot, width = 10, height = 15, dpi = 300)

print("--- Plot saved as indicator_species_plot.png ---")
