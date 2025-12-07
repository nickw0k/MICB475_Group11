# --- 1. Load Libraries ---
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ape)
library(forcats)
library(dplyr)

# --- 2. Load Data ---
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
otuFP  <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
taxFP  <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
phyFP  <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"

meta <- read.delim(metaFP, sep = "\t", stringsAsFactors = FALSE)
otu  <- read.delim(otuFP, sep = "\t", skip = 1, stringsAsFactors = FALSE)
tax  <- read.delim(taxFP, sep = "\t", stringsAsFactors = FALSE)
phy  <- read.tree(phyFP)

# --- 3. Phyloseq Data Prep ---
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$X.OTU.ID
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge") %>%
  mutate(across(everything(), ~ str_replace_all(.x, "^[a-z]__", ""))) %>%
  as.data.frame(stringsAsFactors = FALSE)

tax_mat2 <- as.matrix(tax_mat[, -1])
rownames(tax_mat2) <- tax$Feature.ID
TAX <- tax_table(tax_mat2)

rownames(meta) <- meta$sample.id
META <- sample_data(meta)

ps <- phyloseq(OTU, TAX, META, phy)

# --- 4. Define sex-severity groups ---
sample_data(ps)$sex[sample_data(ps)$sex %in% c("", " ", "NA")] <- NA

sample_data(ps)$Severity <- NA
sample_data(ps)$Severity[grepl("asymptomatic", sample_data(ps)$Group, ignore.case = TRUE)] <- "Control"
sample_data(ps)$Severity[grepl("ambulatory positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Mild"
sample_data(ps)$Severity[grepl("hospitalized positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Severe"
sample_data(ps)$Severity[grepl("Deceased hospitalized", sample_data(ps)$Group, ignore.case = TRUE)] <- "Deceased"

sample_data(ps)$Group_combined <- paste(sample_data(ps)$sex, sample_data(ps)$Severity, sep = "_")

desired_levels <- c("male_Control","male_Mild","male_Severe","male_Deceased",
                    "female_Control","female_Mild","female_Severe","female_Deceased")
present_levels <- intersect(desired_levels, unique(sample_data(ps)$Group_combined))
sample_data(ps)$Group_combined <- factor(sample_data(ps)$Group_combined, levels = present_levels)

ps <- subset_samples(ps, !is.na(Group_combined))

# --- 5. Filter sparse taxa ---
min_prevalence <- 3
min_count      <- 5
ps_filt <- filter_taxa(ps, function(x) sum(x >= min_count) >= min_prevalence & sum(x) > 10, TRUE)

# --- 6. Collapse to genus level ---
ps_genus <- tax_glom(ps_filt, taxrank = "Genus", NArm = TRUE)

tax_df_genus <- as.data.frame(tax_table(ps_genus)) %>%
  rownames_to_column("OTU") %>%
  mutate(Label = ifelse(!is.na(Genus) & Genus != "", Genus, NA)) %>%
  filter(!is.na(Label)) %>%
  select(OTU, Label)

# Keep only OTUs with valid genus
ps_genus <- prune_taxa(tax_df_genus$OTU, ps_genus)

tax_table(ps_genus) <- tax_table(as.matrix(tax_df_genus %>% column_to_rownames("OTU")))

# --- 7. DESeq2 analysis ---
dds <- phyloseq_to_deseq2(ps_genus, ~ Group_combined)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, fitType = "parametric")

# --- 8. Pairwise contrasts + LFC shrinkage ---
groups <- levels(droplevels(sample_data(ps_genus)$Group_combined))
comparisons <- combn(groups, 2, simplify = FALSE)
results_list <- list()

for (comp in comparisons) {
  a <- comp[1]; b <- comp[2]
  cname <- paste0(a, "_vs_", b)
  
  res_raw <- tryCatch(results(dds, contrast = c("Group_combined", a, b)), error = function(e) NULL)
  if (is.null(res_raw)) next
  
  res_raw_df <- as.data.frame(res_raw) %>% rownames_to_column("OTU")
  
  res_shrunk <- tryCatch(lfcShrink(dds, contrast = c("Group_combined", a, b), type = "normal"), error = function(e) NULL)
  
  if (!is.null(res_shrunk)) {
    res_shrunk_df <- as.data.frame(res_shrunk) %>% rownames_to_column("OTU") %>% select(OTU, log2FoldChange)
    res_combined <- res_raw_df %>%
      left_join(res_shrunk_df, by = "OTU", suffix = c("_raw", "")) %>%
      mutate(contrast = cname) %>%
      select(OTU, log2FoldChange, everything())
  } else {
    res_combined <- res_raw_df %>% mutate(log2FoldChange = log2FoldChange) %>% mutate(contrast = cname)
  }
  
  results_list[[cname]] <- res_combined
}

deseq_all <- bind_rows(results_list)

# --- 9. Annotate genera ---
deseq_all <- left_join(deseq_all, tax_df_genus, by = "OTU")

# --- 10. Filter significant genera ---
if(!"padj" %in% colnames(deseq_all)) {
  sig_deseq <- deseq_all %>% filter(!is.na(pvalue) & pvalue <= 0.05)
} else {
  sig_deseq <- deseq_all %>% filter(!is.na(padj) & padj <= 0.05)
}

sig_deseq <- sig_deseq %>% filter(!is.na(log2FoldChange), !is.na(Label)) %>% mutate(Label = as.character(Label))
sig_deseq_plot <- sig_deseq
contrasts <- unique(sig_deseq_plot$contrast)

# --- 11a. Log2Fold Bar Plots ---
for (c in contrasts) {
  df <- sig_deseq_plot %>% filter(contrast == c)
  df <- df %>% mutate(Label = fct_reorder(Label, log2FoldChange),
                      direction = ifelse(log2FoldChange > 0, "Up", "Down"))
  
  groups <- unlist(strsplit(c, "_vs_"))
  group1 <- groups[1]; group2 <- groups[2]
  
  df <- df %>% mutate(direction_label = case_when(
    direction == "Up" ~ paste("Higher in", group1),
    direction == "Down" ~ paste("Higher in", group2)
  ))
  
  color_mapping <- setNames(c("cornflowerblue", "plum2"),
                            c(paste("Higher in", group1), paste("Higher in", group2)))
  
  p <- ggplot(df, aes(x = Label, y = log2FoldChange, fill = direction_label)) +
    geom_col(show.legend = TRUE) +
    coord_flip() +
    theme_minimal() +
    labs(title = paste0("Significant Genera — ", c),
         x = "Genus",
         y = "log2 Fold Change",
         fill = "Abundance Difference") +
    scale_fill_manual(values = color_mapping) +
    theme(axis.text.y = element_text(size = 9),
          axis.text.x = element_text(size = 9),
          plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = paste0("DESeq2_Genus_", c, ".png"), plot = p, width = 10, height = 6, dpi = 300)
}

# --- 11b. Volcano Plots ---
for (c in contrasts) {
  df <- sig_deseq_plot %>% filter(contrast == c)
  df_volcano <- df %>% filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    mutate(neglog10_padj = -log10(padj))
  
  groups <- unlist(strsplit(c, "_vs_"))
  group1 <- groups[1]; group2 <- groups[2]
  
  df_volcano <- df_volcano %>% mutate(significance = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ paste("Higher in", group1),
    padj < 0.05 & log2FoldChange < 0 ~ paste("Higher in", group2),
    TRUE ~ "Not Significant"
  ))
  
  color_mapping <- setNames(c("#1f78b4", "#e31a1c", "grey70"),
                            c(paste("Higher in", group1),
                              paste("Higher in", group2),
                              "Not Significant"))
  
  volcano_plot <- ggplot(df_volcano, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = significance), size = 2, alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    scale_color_manual(values = color_mapping) +
    theme_minimal() +
    labs(title = paste0("Volcano Plot — ", c),
         x = "log2 Fold Change",
         y = "-log10(adj p-value)",
         color = "Significance") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(filename = paste0("Volcano_Genus_", c, ".png"),
         plot = volcano_plot, width = 8, height = 8, dpi = 300)
}
