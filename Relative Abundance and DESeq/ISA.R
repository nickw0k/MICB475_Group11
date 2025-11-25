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

# --- 2. Loading Files ---

metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(file=metaFP,sep = "\t")

otuFP <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
otu <- read.delim(file=otuFP,sep = "\t", skip = 1)

taxFP <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
tax <- read.delim(file=taxFP,sep = "\t")

phyFP <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"
phy <- read.tree(phyFP)

# --- 3. Phyloseq Data Prep ---

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

ps <- phyloseq(OTU,META,TAX,phy)

# --- 4. Separate phyloseq object into sex-severity groups ---

meta <- data.frame(sample_data(ps), stringsAsFactors = FALSE)
meta$Severity <- NA_character_
meta$Severity[grepl("asymptomatic", tolower(meta$Group))] <- "Control"
meta$Severity[grepl("ambulatory positive",   tolower(meta$Group))] <- "Mild"
meta$Severity[grepl("hospitalized positive", tolower(meta$Group))] <- "Severe"
meta$Severity[grepl("Deceased hospitalized",   tolower(meta$Group))] <- "Fatal"

sample_data(ps) <- sample_data(meta)

meta <- data.frame(sample_data(ps))
meta$Group_combined <- paste(meta$sex, meta$Severity, sep = "_")
sample_data(ps) <- sample_data(meta)


# --- 5. Run Indicator Species Analysis (ISA) ---

otu_mat <- as.data.frame(t(otu_table(ps)))

meta_df <- data.frame(sample_data(ps), stringsAsFactors = FALSE)

meta_clean <- meta_df[!is.na(meta_df$Group_combined), ]

otu_clean <- otu_mat[rownames(meta_clean), ]

otu_clean[is.na(otu_clean)] <- 0

isa_mpt <- multipatt(otu_clean, meta_clean$Group_combined,
                     func = "IndVal.g", duleg = TRUE,
                     control = how(nperm = 999))

summary(isa_mpt, indvalcomp = TRUE, alpha = 0.05)


# --- 6. Programmatically Extract Significant Results (p.value <= 0.05) ---

sig_table <- data.frame(isa_mpt$sign) %>%
  rownames_to_column(var = "OTU") %>%
  filter(p.value <= 0.05)

group_names <- colnames(isa_mpt$comb)

sig_table$cluster <- group_names[sig_table$index]

significant_indicators <- sig_table %>%
  select(OTU, stat, p.value, cluster)

print(head(significant_indicators))


# --- 7. Annotate Results with Taxonomy ---


tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column(var = "OTU") %>%
  
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
  left_join(tax_df %>% select(OTU, Label), by = "OTU") %>%

  filter(!grepl("NA", cluster))


# --- 8. Create the Final Plot ---
isa_plot <- ggplot(top_sps_named, aes(x = reorder(Label, stat), y = stat, fill = cluster)) +
  geom_col() +
  coord_flip() + 

  scale_y_continuous(limits = c(0, 1.0)) + # Map to the original y-axis (stat)
  labs(
    x = "Taxonomy (Taxa Sorted by Indicator Value)",
    y = "Indicator Value (stat)",
    fill = "Sex-Severity Group",
    title = "Significant Indicator Species (IndVal Max = 1.0)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom", 
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid.minor.y = element_blank() # Clean up the vertical grid lines
  )

print(isa_plot)
ggsave("indicator_species_plot.png", plot = isa_plot, width = 10, height = 15, dpi = 300)





## For top 20 Significant Taxa 

#top_sps_named <- significant_indicators %>%
# left_join(tax_df %>% select(OTU, Label), by = "OTU") %>%

# arrange(desc(stat)) %>%  
# slice_head(n = 20)     









# Table of Taxon + Indicator Value

#final_isa_table <- sig_table %>%
 # left_join(tax_df %>% select(OTU, Label), by = "OTU") %>%
  
 # filter(!grepl("NA", cluster)) %>%

 # select(Taxon_Label = Label, 
         Indicated_Group = cluster, 
         Indicator_Value_Stat = stat, 
         P_Value = p.value) %>%
  
 # filter(Taxon_Label != "Unclassified") %>%
 
 # arrange(desc(Indicator_Value_Stat)) %>%
 
 # mutate(
  #  Indicator_Value_Stat = round(Indicator_Value_Stat, 3),
  #  P_Value = format.pval(P_Value, digits = 4)
  )

#cat("\n\n#################################################################\n")
#cat("--- Significant Indicator Species Table (p-value <= 0.05) ---\n")
#cat("#################################################################\n\n")

#print(as.data.frame(final_isa_table), max.print = 99999) 
