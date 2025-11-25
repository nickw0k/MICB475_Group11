library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ape)
library(forcats)

# --- 2. Load files ---
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(metaFP, sep="\t")

otuFP <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
otu <- read.delim(otuFP, sep="\t", skip=1)

taxFP <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
tax <- read.delim(taxFP, sep="\t")

phyFP <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"
phy <- read.tree(phyFP)

# --- 3. Build phyloseq object ---
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

# --- 4. Create severity + Group_combined ---

meta$sex[meta$sex %in% c("", " ", "NA")] <- NA

# --- 3. Create Severity column ---
meta$Severity <- NA
meta$Severity[grepl("asymptomatic",           meta$Group, ignore.case=TRUE)] <- "Control"
meta$Severity[grepl("ambulatory positive",    meta$Group, ignore.case=TRUE)] <- "Mild"
meta$Severity[grepl("hospitalized positive",  meta$Group, ignore.case=TRUE)] <- "Severe"
meta$Severity[grepl("Deceased hospitalized",  meta$Group, ignore.case=TRUE)] <- "Fatal"

# --- 4. Create Group_combined ---
meta$Group_combined <- paste(meta$sex, meta$Severity, sep="_")

# --- 5. Make it a factor in the desired order ---
meta$Group_combined <- factor(
  meta$Group_combined,
  levels = c(
    "male_Control","male_Mild","male_Severe","male_Fatal",
    "female_Control","female_Mild","female_Severe","female_Fatal"
  )
)

# --- 6. Align rownames with phyloseq sample names ---
# IMPORTANT: phyloseq sample_names must match rownames exactly
rownames(meta) <- meta$sample.id  # if not already
meta_aligned <- meta[match(sample_names(ps), rownames(meta)), ]

# Check for mismatches
if(any(is.na(meta_aligned$sample.id))){
  stop("Some phyloseq samples have no matching metadata! Check sample IDs.")
}

# --- 7. Assign cleaned metadata back to phyloseq object ---
sample_data(ps) <- sample_data(meta_aligned)

ps_clean <- subset_samples(ps, !is.na(Group_combined))

# --- 2. Double-check sample names alignment ---
#sample_names(ps_clean)  # only samples with valid Group_combined remain

# --- 3. Create DESeq2 object ---
dds <- phyloseq_to_deseq2(ps_clean, ~ Group_combined)

# --- 4. Estimate size factors (poscounts is good for sparse data) ---
dds <- estimateSizeFactors(dds, type = "poscounts")

# --- 5. Run DESeq ---
dds <- DESeq(dds, fitType = "parametric")

# --- 3. Create all pairwise contrasts between levels of Group_combined ---
groups <- levels(dds$Group_combined)
comparisons <- combn(groups, 2, simplify = FALSE)

deseq_list <- list()
for (comp in comparisons) {
  cname <- paste(comp[1], "vs", comp[2], sep = "_")
  res <- tryCatch(
    results(dds, contrast = c("Group_combined", comp[1], comp[2])),
    error = function(e) {
      message("results() error for ", cname, ": ", e$message)
      return(NULL)
    }
  )
  if (!is.null(res)) {
    res_df <- as.data.frame(res) %>%
      rownames_to_column("OTU") %>%
      mutate(contrast = cname)
    deseq_list[[cname]] <- res_df
  }
}
deseq_all <- bind_rows(deseq_list)

# --- 4. Annotate with taxonomy (Label column) ---
tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column(var = "OTU") %>%
  mutate(
    Species = na_if(Species, "s__"),
    Genus   = na_if(Genus,   "g__"),
    Family  = na_if(Family,  "f__"),
    Order   = na_if(Order,   "o__"),
    Class   = na_if(Class,   "c__"),
    Phylum  = na_if(Phylum, "p__"),
    Domain  = na_if(Domain, "d__")
  ) %>%
  mutate(
    Label = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain),
    Label = ifelse(is.na(Label), "Unclassified", Label)
  )

deseq_annot <- deseq_all %>%
  left_join(tax_df %>% select(OTU, Label), by = "OTU")

# --- 5. Clean extreme / infinite LFCs (these produce wild axes like 90000) ---
# Replace Inf / -Inf with NA and optionally filter extreme values
deseq_annot <- deseq_annot %>%
  mutate(
    log2FoldChange = as.numeric(log2FoldChange),
    log2FoldChange = ifelse(is.infinite(log2FoldChange), NA, log2FoldChange),
    is_extreme = ifelse(!is.na(log2FoldChange) & abs(log2FoldChange) > 100, TRUE, FALSE)
  )

# OPTIONAL: If you want to exclude absurd extreme values from plotting entirely uncomment:
# deseq_annot <- deseq_annot %>% filter(!is_extreme)

# NOTE: extreme values often mean problems with low counts or MLE instability.
# Recommended: use lfcShrink (apeglm) for more stable LFC estimates (see note below).

# --- 6. Keep only significant results (padj <= 0.05) ---
sig_deseq <- deseq_annot %>%
  filter(!is.na(padj) & padj <= 0.05)

# If no significant hits, inform and stop plot
if (nrow(sig_deseq) == 0) {
  message("No significant taxa found at padj <= 0.05 for any contrast.")
} else {
  sig_deseq <- sig_deseq %>%
    filter(!is.na(log2FoldChange), !is.na(Label)) %>%
    group_by(contrast) %>%
    mutate(Label = fct_reorder(as.factor(Label), abs(log2FoldChange))) %>%
    ungroup()
}

# --- 7. Plot: faceted colored bars, cap y-axis visually to ±10 with coord_cartesian ---
deseq_plot <- ggplot(sig_deseq, aes(x = Label, y = log2FoldChange, fill = contrast)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~ contrast, scales = "free_y", ncol = 4) +
  theme_minimal() +
  labs(
    title = "Significant Taxa (DESeq2, padj ≤ 0.05) — log2 Fold Change",
    x = "Taxon",
    y = "log2FoldChange"
  ) +
  # Use coord_cartesian to cap display without dropping points from data
  coord_cartesian(ylim = c(-10, 10)) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(hjust = 0.5)
  )

# print(deseq_plot)  #all plots in one 
# ggsave("DESeq2_pairwise_faceted_log2FC.png", plot = deseq_plot, width = 20, height = 20, dpi = 300)



# saves each plot individually
sig_deseq <- sig_deseq %>%
  filter(!is.na(log2FoldChange), !is.na(Label))

contrasts <- unique(sig_deseq$contrast)

for (c in contrasts) {
  df <- sig_deseq %>% filter(contrast == c)
  
  df <- df %>% mutate(Label = fct_reorder(as.factor(Label), abs(log2FoldChange)))
  
  p <- ggplot(df, aes(x = Label, y = log2FoldChange, fill = log2FoldChange)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(
      title = paste0("Significant Taxa — ", c),
      x = "Taxon",
      y = "log2 Fold Change"
    ) +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
  
 
  ggsave(filename = paste0("DESeq2_", c, ".png"),
         plot = p, width = 8, height = 6, dpi = 300)
}




# ---  Summary of significant counts per contrast  ---

#sig_deseq <- as.data.frame(sig_deseq)
#sig_deseq$contrast <- as.character(sig_deseq$contrast)
# sig_summary <- sig_deseq %>% 
#dplyr::count(contrast)

# print(sig_summary)

# ---  Summary of extreme values for troubleshooting ---


#deseq_annot <- as.data.frame(deseq_annot)

#extreme_report <- deseq_annot %>% 
# dplyr::filter(is_extreme) %>% 
#  dplyr::select(OTU, contrast, log2FoldChange, padj, Label)

# if (nrow(extreme_report) > 0) {
#   message("Warning: some contrasts had extreme |log2FoldChange| > 100. See 'extreme_report' for details.")
#  print(head(extreme_report, 20))
#  }
