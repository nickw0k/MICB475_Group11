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
otu_mat <- as.matrix(otu[ , -1])
rownames(otu_mat) <- otu$X.OTU.ID
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge") %>%
  mutate(across(everything(), ~ str_replace_all(.x, "^[a-z]__", ""))) %>%
  as.data.frame(stringsAsFactors = FALSE)

tax_mat2 <- as.matrix(tax_mat[ , -1])
rownames(tax_mat2) <- tax$Feature.ID
TAX <- tax_table(tax_mat2)

rownames(meta) <- meta$sample.id
META <- sample_data(meta)

ps <- phyloseq(OTU, TAX, META, phy)


# --- 4. Separate phyloseq object into sex-severity groups ---

sample_data(ps)$sex[sample_data(ps)$sex %in% c("", " ", "NA")] <- NA

sample_data(ps)$Severity <- NA
sample_data(ps)$Severity[grepl("asymptomatic", sample_data(ps)$Group, ignore.case = TRUE)] <- "Control"
sample_data(ps)$Severity[grepl("ambulatory positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Mild"
sample_data(ps)$Severity[grepl("hospitalized positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Severe"
sample_data(ps)$Severity[grepl("Deceased hospitalized", sample_data(ps)$Group, ignore.case = TRUE)] <- "Fatal"

sample_data(ps)$Group_combined <- paste(sample_data(ps)$sex, sample_data(ps)$Severity, sep = "_")

desired_levels <- c("male_Control","male_Mild","male_Severe","male_Fatal",
                    "female_Control","female_Mild","female_Severe","female_Fatal")
present_levels <- intersect(desired_levels, unique(sample_data(ps)$Group_combined))
sample_data(ps)$Group_combined <- factor(sample_data(ps)$Group_combined, levels = present_levels)

ps <- subset_samples(ps, !is.na(Group_combined))


# --- 5. Filter Sparse Taxa --- 

min_prevalence <- 3
min_count      <- 5
ps_filt <- filter_taxa(ps, function(x) sum(x >= min_count) >= min_prevalence & sum(x) > 10, TRUE)


# --- 6. DESeq Analysis ---

dds <- phyloseq_to_deseq2(ps_filt, ~ Group_combined)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, fitType = "parametric")


# --- 7. Pairwise Contrasts + LFC Shrinkage ---

groups <- levels(droplevels(sample_data(ps_filt)$Group_combined))
comparisons <- combn(groups, 2, simplify = FALSE)

results_list <- list()

for (comp in comparisons) {
  a <- comp[1]; b <- comp[2]
  cname <- paste0(a, "_vs_", b)
  
  res_raw <- tryCatch(results(dds, contrast = c("Group_combined", a, b)), error = function(e) NULL)
  if(is.null(res_raw)) next
  
  res_raw_df <- as.data.frame(res_raw) %>% rownames_to_column("OTU")
  
  res_shrunk <- tryCatch(lfcShrink(dds, contrast = c("Group_combined", a, b), type = "normal"), error = function(e) NULL)
  
  if(!is.null(res_shrunk)) {
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


# --- 8. Annotate Taxonomy ---

tax_df <- as.data.frame(tax_table(ps_filt)) %>%
  rownames_to_column("OTU") %>%
  mutate(across(everything(), ~ na_if(.x, ""))) %>%

  mutate(
    Label = case_when(
      !is.na(Species) & Species != "" ~ paste0("s__", Species),
      !is.na(Genus)   & Genus   != "" ~ paste0("g__", Genus),
      !is.na(Family)  & Family  != "" ~ paste0("f__", Family),
      !is.na(Order)   & Order   != "" ~ paste0("o__", Order),
      !is.na(Class)   & Class   != "" ~ paste0("c__", Class),
      !is.na(Phylum)  & Phylum  != "" ~ paste0("p__", Phylum),
      !is.na(Domain)  & Domain  != "" ~ paste0("d__", Domain),
      TRUE                            ~ OTU
    )
  ) %>%
  
  select(OTU, Label)

deseq_all <- left_join(deseq_all, tax_df, by = "OTU")


# --- 9. Filter Significant OTUs ---

if(!"padj" %in% colnames(deseq_all)) {
  sig_deseq <- deseq_all %>% filter(!is.na(pvalue) & pvalue <= 0.05)
} else {
  sig_deseq <- deseq_all %>% filter(!is.na(padj) & padj <= 0.05)
}

sig_deseq <- sig_deseq %>% filter(!is.na(log2FoldChange), !is.na(Label)) %>% mutate(Label = as.character(Label))

sig_deseq_plot <- sig_deseq %>%
  filter(!is.na(padj) & padj <= 0.05) %>%
  filter(!is.na(Label), !is.na(log2FoldChange)) %>%
  mutate(Label = as.character(Label)) 

contrasts <- unique(sig_deseq_plot$contrast)

for (c in contrasts) {
  df <- sig_deseq_plot %>% filter(contrast == c)
  
  df <- df %>% mutate(Label = fct_reorder(Label, log2FoldChange))
  

  df <- df %>% mutate(direction = ifelse(log2FoldChange > 0, "Up", "Down"))
  
  
  # ---10. Plot ---
  p <- ggplot(df, aes(x = Label, y = log2FoldChange, fill = direction)) +
    geom_col(show.legend = TRUE) +
    coord_flip() +  
    theme_minimal() +
    labs(
      title = paste0("Significant OTUs — ", c),
      x = "OTU",
      y = "log2 Fold Change"
    ) +
    scale_fill_manual(values = c("Up" = "plum2", "Down" = "cornflowerblue")) +
    theme(
      axis.text.y = element_text(size = 6),  # smaller font if many OTUs
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(filename = paste0("DESeq2_", c, "_individual_OTUs.png"),
         plot = p, width = 10, height = 6, dpi = 300)}

