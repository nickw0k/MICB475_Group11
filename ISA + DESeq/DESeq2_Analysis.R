library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ape)
library(forcats)
library(ashr)

# --- 1. Load Data ---
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
otuFP  <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
taxFP  <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
phyFP  <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"

meta <- read.delim(metaFP, sep="\t")
otu  <- read.delim(otuFP, sep="\t", skip=1)
tax  <- read.delim(taxFP, sep="\t")
phy  <- read.tree(phyFP)

# --- 2. Build Phyloseq ---
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$X.OTU.ID
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

meta_df <- meta[,-1]
rownames(meta_df) <- meta$sample.id
META <- sample_data(meta_df)

tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into=c("Domain","Phylum","Class","Order",
                         "Family","Genus","Species"),
           sep=";", fill="right") %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$Feature.ID
TAX <- tax_table(tax_mat)

ps <- phyloseq(OTU, META, TAX, phy)

# --- 3. Clean Metadata & Grouping ---
meta$sex[meta$sex %in% c("", " ", "NA")] <- NA

meta$Severity <- NA
meta$Severity[grepl("asymptomatic",           meta$Group, ignore.case=TRUE)] <- "Control"
meta$Severity[grepl("ambulatory positive",    meta$Group, ignore.case=TRUE)] <- "Mild"
meta$Severity[grepl("hospitalized positive",  meta$Group, ignore.case=TRUE)] <- "Severe"
meta$Severity[grepl("Deceased hospitalized",  meta$Group, ignore.case=TRUE)] <- "Fatal"

meta$Group_combined <- paste(meta$sex, meta$Severity, sep="_")
meta$Group_combined <- factor(
  meta$Group_combined,
  levels = c("male_Control","male_Mild","male_Severe","male_Fatal",
             "female_Control","female_Mild","female_Severe","female_Fatal")
)

rownames(meta) <- meta$sample.id  
meta_aligned <- meta[match(sample_names(ps), rownames(meta)), ]
sample_data(ps) <- sample_data(meta_aligned)
ps_clean <- subset_samples(ps, !is.na(Group_combined))

# --- 4. Run DESeq2 ---
dds <- phyloseq_to_deseq2(ps_clean, ~ Group_combined)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, fitType = "parametric")

# --- 5. Generate Results with SHRINKAGE ---
groups <- levels(dds$Group_combined)
comparisons <- combn(groups, 2, simplify = FALSE)

deseq_list <- list()

for (comp in comparisons) {
  cname <- paste(comp[1], "vs", comp[2], sep = "_")
  
  # lfcShrink fixes the "infinite" values statistically
  res <- tryCatch(
    lfcShrink(dds, contrast = c("Group_combined", comp[1], comp[2]), type = "ashr"),
    error = function(e) { message("Error in ", cname); return(NULL) }
  )
  
  if (!is.null(res)) {
    res_df <- as.data.frame(res) %>%
      rownames_to_column("OTU") %>%
      mutate(contrast = cname)
    deseq_list[[cname]] <- res_df
  }
}
deseq_all <- bind_rows(deseq_list)

# --- 6. Annotate and Fix LABELS (The "Stacking" Fix) ---
tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column(var = "OTU") %>%
  mutate(
    # Clean up names
    across(c(Species, Genus, Family, Order, Class, Phylum, Domain), ~na_if(., "s__")),
    across(c(Species, Genus, Family, Order, Class, Phylum, Domain), ~na_if(., "g__")), 
    # (Repeat for other levels if needed or use regex replacement)
    
    Label = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain),
    Label = ifelse(is.na(Label), "Unclassified", Label),
    
    # CRITICAL FIX: Make Label Unique so bars don't stack
    # We take the first 5 characters of the OTU ID to keep it short
    Label_Unique = paste0(Label, " (", substr(OTU, 1, 5), ")")
  )

deseq_annot <- deseq_all %>%
  left_join(tax_df %>% select(OTU, Label, Label_Unique), by = "OTU")

# --- 7. Filter Significant Hits ---
sig_deseq <- deseq_annot %>%
  filter(!is.na(padj) & padj <= 0.05) %>%
  filter(!is.na(log2FoldChange), !is.na(Label_Unique))

# --- 8. Loop to Save Vertical Plots ---

contrasts <- unique(sig_deseq$contrast)

for (c in contrasts) {
  
 
  df <- sig_deseq %>% filter(contrast == c)
  if(nrow(df) == 0) next
  

  df_top <- df %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = 20)
  
  
  df_top <- df_top %>% mutate(Label_Unique = fct_reorder(as.factor(Label_Unique), -log2FoldChange))
  
  
  p <- ggplot(df_top, aes(x = Label_Unique, y = log2FoldChange, fill = log2FoldChange)) +
    geom_col() +
   
    
    theme_minimal() +
    labs(
      title = paste0("Top 20 Significant Taxa: ", c),
      x = "Taxon",
      y = "log2 Fold Change"
    ) +
    scale_fill_gradient2(low = "blue", mid = "grey", high = "red", midpoint = 0) +
    coord_cartesian(ylim = c(-10, 10)) +
    theme(
     
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10)
    )
  
  filename <- paste0("DESeq2_Vertical_", c, ".png")
  ggsave(filename = filename, plot = p, width = 12, height = 8, dpi = 300)
  
}
