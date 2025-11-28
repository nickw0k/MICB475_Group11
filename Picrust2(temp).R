### Paths (relative to your Git Project root)
abundance_file <- "Mexico Dataset QIIME2 files/path_abun_unstrat.tsv.gz"
metadata_file  <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Create a list of all the packages you need to install 
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2") 

# Use the above list to install all the packages using a for loop 
for (pkg in pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg) } 

# when asked if you want to update all, some or none please type "n" for none 

# After installing all of its above dependencies, install ggpicrust2 (no need if already installed)
install.packages("ggpicrust2")

library(readr)
library(ggpicrust2)
library(tibble)
library(dplyr)
library(ggplot2)
library(ggprism)
library(patchwork)
library(DESeq2)
library(ggh4x)
library(stringr)
library(cowplot)
library(matrixStats)


# 1. Read PICRUSt2 unstratified pathway abundance table
abundance_data <- read_tsv(
  abundance_file,
  col_names = TRUE,
  trim_ws   = TRUE
) |> as.data.frame()

# 2. Read metadata
metadata <- read_tsv(metadata_file)

### 3. Build severity + severity_sex groups and drop MildNegative ####
metadata <- metadata |>
  dplyr::mutate(
    # 4 severity levels from the long Group text
    severity = dplyr::case_when(
      stringr::str_detect(Group, regex("asymptomatic",          ignore_case = TRUE)) ~ "Asymp",
      stringr::str_detect(Group, regex("ambulatory negative",   ignore_case = TRUE)) ~ "Mild_neg",
      stringr::str_detect(Group, regex("ambulatory positive",   ignore_case = TRUE)) ~ "Mild+",
      stringr::str_detect(Group, regex("Deceased",              ignore_case = TRUE)) ~ "Deceased",
      stringr::str_detect(Group, regex("hospitalized positive", ignore_case = TRUE)) ~ "Severe+",
      TRUE ~ NA_character_
    ),
    sex_short = dplyr::case_when(
      sex == "female" ~ "F",
      sex == "male"   ~ "M",
      TRUE            ~ NA_character_
    ),
    severity_sex = paste0(severity, "_", sex_short)
  ) |>
  # drop Mild negatives and any rows without clear severity/sex
  dplyr::filter(severity != "Mild_neg") |>
  dplyr::filter(!is.na(severity_sex))

# lock factor order (also forces Deceased F/M to be kept if present)
metadata$severity_sex <- factor(
  metadata$severity_sex,
  levels = c("Asymp_F","Asymp_M",
             "Mild+_F","Mild+_M",
             "Severe+_F","Severe+_M",
             "Deceased_F","Deceased_M")
)

table(metadata$severity_sex)

#### 4. Get sample IDs that remain after filtering ####
sample_names <- metadata$`sample-id`

#### 5. Keep only 'pathway' + those samples, in the same order ####
keep_cols <- c("pathway", sample_names)
abundance_data_filtered <- abundance_data[, colnames(abundance_data) %in% keep_cols]
abundance_data_filtered <- abundance_data_filtered[, keep_cols]

#### 6. Drop samples that are all zero (no counts) ####
non_zero <- colSums(abundance_data_filtered[, -1] != 0) > 0
abundance_data_filtered <- abundance_data_filtered[, c(TRUE, non_zero)]

#### 7. Align metadata rows and order with filtered abundance table ####
abun_samples <- colnames(abundance_data_filtered)[-1]
metadata <- metadata[metadata$`sample-id` %in% abun_samples, ]
metadata <- metadata[match(abun_samples, metadata$`sample-id`), ]
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$`sample-id`
metadata$sample_name <- metadata$`sample-id`

#### 8. Differential abundance testing at pathway level (DESeq2 on severity_sex) ####
abundance_daa_results_df <- pathway_daa(
  abundance  = abundance_data_filtered %>% tibble::column_to_rownames("pathway"),
  metadata   = metadata,
  group      = "severity_sex",
  daa_method = "DESeq2"
)

#### 9. Add MetaCyc pathway descriptions ####
metacyc_daa_annotated_results_df <- pathway_annotation(
  pathway        = "MetaCyc",
  daa_results_df = abundance_daa_results_df,
  ko_to_kegg     = FALSE
)

#### 10. Keep only significant pathways (p < 0.05) ####
feature_with_p_0.05 <- abundance_daa_results_df %>%
  dplyr::filter(p_values < 0.05)

#### 11. Replace feature IDs with descriptions in the results table ####
feature_desc <- dplyr::inner_join(
  feature_with_p_0.05,
  metacyc_daa_annotated_results_df,
  by = "feature"
)
feature_desc$feature <- feature_desc$description
feature_desc <- feature_desc[, 1:7]        # keep same columns as original results
colnames(feature_desc) <- colnames(feature_with_p_0.05)

#### 12. Subset abundance table to significant pathways only ####
abundance <- abundance_data_filtered %>%
  dplyr::filter(pathway %in% feature_with_p_0.05$feature)

colnames(abundance)[1] <- "feature"

#### 13. Join abundance with MetaCyc descriptions ####
abundance_desc <- dplyr::inner_join(
  abundance,
  metacyc_daa_annotated_results_df,
  by = "feature"
)
abundance_desc$feature <- abundance_desc$description

# keep only feature + all sample columns (drop annotation-only columns)
n_samples <- ncol(abundance_data_filtered) - 1  # minus the "pathway" column
abundance_desc <- abundance_desc[, 1:(1 + n_samples)]

#### 14. Collapse abundances to one column per severity×sex group ####

# drop any rows that don't have a pathway description
abundance_desc <- abundance_desc |> dplyr::filter(!is.na(feature))

# pathways × samples matrix
abun_mat <- abundance_desc %>%
  tibble::column_to_rownames("feature")

# ensure that columns are in the same order as metadata rows
abun_mat <- abun_mat[, metadata$`sample-id`]

group_vec <- metadata$severity_sex
names(group_vec) <- metadata$`sample-id`
group_levels <- levels(group_vec)

# mean abundance per group (pathways × groups)
abun_group_mat <- sapply(group_levels, function(g) {
  idx <- group_vec == g
  rowMeans(abun_mat[, idx, drop = FALSE])
})

abundance_group <- as.data.frame(abun_group_mat)
abundance_group$feature <- rownames(abundance_group)
abundance_group <- abundance_group %>% dplyr::relocate(feature)

# one-row metadata for the 8 groups
metadata_group <- data.frame(
  `sample-id`  = group_levels,
  severity_sex = factor(group_levels, levels = group_levels)
)
rownames(metadata_group) <- metadata_group$`sample-id`

#### 15. Heatmap of significant pathways (group-level) + legend + export ####

#### 15a. Matrix for heatmap (unchanged) ####
abun_mat <- abundance_desc |>
  dplyr::filter(!is.na(feature)) |>
  tibble::column_to_rownames("feature")

#### 15b. Draw heatmap grouped by severity_sex ####
heat_severity <- pathway_heatmap(
  abundance = abun_mat,
  metadata  = metadata,
  group     = "severity_sex"
)

# Remove ALL x-axis labels/ticks and keep right side padded for Z-score bar
heat_severity <- heat_severity +
  theme(
    axis.text.x  = element_blank(),   # <- removes “Asymp_F … Deceased_M”
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 10),
    plot.margin  = margin(5.5, 40, 5.5, 5.5)  # top, right, bottom, left
  ) +
  coord_cartesian(clip = "off")

#### 15c. Build a clean legend for severity_sex ####
dummy_leg <- ggplot(metadata, aes(x = severity_sex, fill = severity_sex)) +
  geom_bar() +
  theme_minimal(base_size = 10) +
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    panel.grid  = element_blank(),
    legend.position      = "bottom",
    legend.title         = element_text(size = 10),
    legend.text          = element_text(size = 9),
    legend.key.width     = unit(0.8, "cm"),
    legend.box.margin    = margin(t = 0, r = 0, b = 0, l = 0)
  )

leg_severity <- cowplot::get_legend(dummy_leg)

#### 15d. Combine heatmap + legend (legend pulled UP closer to plot) ####
final_heatmap <- cowplot::plot_grid(
  heat_severity,
  leg_severity,
  ncol        = 1,
  rel_heights = c(4, 0.7),   # increase 2nd value if you want a bit more legend space
  align       = "v"
)

ggsave(
  filename = "Fig1_pathway_heatmap_severity_sex.png",
  plot     = final_heatmap,
  width    = 14,
  height   = 7,   # extra height so legend is fully in frame
  dpi      = 300
)

# 16. PCA on all pathways

pca_severity <- pathway_pca(
  abundance = abundance_data_filtered %>% column_to_rownames("pathway"),
  metadata  = metadata,
  group     = "severity_sex"
)

pca_severity <- pca_severity +
  geom_point(size = 2.5, alpha = 0.6) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    legend.title = element_text(size = 12),
    legend.text  = element_text(size = 10)
  )

ggsave(
  "Fig2_pathway_PCA_severity_sex.png",
  plot   = pca_severity,
  width  = 8, height = 6, dpi = 300
)

### 17. DESeq2: all pairwise contrasts on severity_sex ###

# Function: run DESeq2 for all pairwise combinations of a metadata factor
DEseq2_function <- function(abundance_table, metadata, col_of_interest) {
  
  # 17.1 Copy metadata and rename the grouping column to a stable name ----
  DESeq2_metadata <- as.data.frame(metadata)
  
  DESeq2_colnames <- colnames(DESeq2_metadata)
  DESeq2_colnames[DESeq2_colnames == col_of_interest] <- "Group_group_nonsense"
  colnames(DESeq2_metadata) <- DESeq2_colnames
  DESeq2_metadata$Group_group_nonsense <- factor(DESeq2_metadata$Group_group_nonsense)
  
  # 17.2 Build count matrix with pathways as rows, samples as columns ----
  DESeq2_abundance_mat <- abundance_table |>
    tibble::column_to_rownames("pathway")
  
  # ensure column order of counts == metadata row order
  DESeq2_abundance_mat <- DESeq2_abundance_mat[, rownames(DESeq2_metadata), drop = FALSE]
  
  # 17.3 Generate all pairwise combinations of groups (as CHARACTERS) ----
  groups_char <- unique(as.character(DESeq2_metadata$Group_group_nonsense))
  groups_char <- groups_char[!is.na(groups_char)]
  DESeq2_combinations <- utils::combn(groups_char, 2)  # 2 x N matrix
  
  DESeq2_results_list <- vector("list", ncol(DESeq2_combinations))
  
  message("Performing pairwise comparisons with DESeq2...")
  
  # 17.4 Loop over each pair of groups and run DESeq2 ----
  for (i in seq_len(ncol(DESeq2_combinations))) {
    g1 <- DESeq2_combinations[1, i]
    g2 <- DESeq2_combinations[2, i]
    
    message("  Contrast: ", g2, " vs ", g1)
    
    # subset samples in these two groups
    sub_idx <- DESeq2_metadata$Group_group_nonsense %in% c(g1, g2)
    DESeq2_metadata_sub  <- droplevels(DESeq2_metadata[sub_idx, , drop = FALSE])
    DESeq2_abundance_sub <- DESeq2_abundance_mat[, sub_idx, drop = FALSE]
    DESeq2_abundance_sub <- round(DESeq2_abundance_sub)
    
    # DESeq2 object
    DESeq2_object <- DESeq2::DESeqDataSetFromMatrix(
      countData = DESeq2_abundance_sub,
      colData   = DESeq2_metadata_sub,
      design    = ~ Group_group_nonsense
    )
    
    DESeq2_object <- BiocGenerics::estimateSizeFactors(DESeq2_object, type = "poscounts")
    DESeq2_object <- DESeq2::DESeq(DESeq2_object)
    
    # results: g2 vs g1
    res <- as.data.frame(
      DESeq2::results(
        DESeq2_object,
        contrast = c("Group_group_nonsense", g2, g1)
      )
    )
    
    res$feature  <- rownames(res)
    res$contrast <- paste0(g2, "_vs_", g1)
    
    DESeq2_results_list[[i]] <- res
  }
  
  # 17.5 Combine all contrasts into one data frame ----
  all_results <- dplyr::bind_rows(DESeq2_results_list)
  return(all_results)
}

### 17a. Run DESeq2 on your PICRUSt2 pathway table (severity_sex) ----

deseq_res <- DEseq2_function(
  abundance_table = abundance_data_filtered,  # has "pathway" + samples
  metadata        = metadata,
  col_of_interest = "severity_sex"
)

# Check what contrasts you have:
unique(deseq_res$contrast)

### 17b. Annotate with MetaCyc pathway descriptions ----

deseq_res_annot <- deseq_res %>%
  dplyr::inner_join(metacyc_daa_annotated_results_df, by = "feature") %>%
  dplyr::select(
    feature,
    description,
    contrast,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  )

### 17c. Filter to FDR-significant pathways (padj < 0.05) ----

sig_res <- deseq_res_annot %>%
  dplyr::filter(!is.na(padj), padj < 0.05)

# See how many sig pathways per contrast
table(sig_res$contrast)

### 17d. Make and save a log2FC barplot for *every* contrast with ≥1 sig pathway ----
all_contrasts <- unique(sig_res$contrast)

for (cname in all_contrasts) {
  plot_data <- sig_res %>%
    dplyr::filter(contrast == cname) %>%
    # order by significance, then effect size
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
    # keep only the top 20 pathways
    dplyr::slice_head(n = 20)
  
  if (nrow(plot_data) == 0) next
  
  p <- ggplot(plot_data,
              aes(y    = reorder(description, log2FoldChange),
                  x    = log2FoldChange,
                  fill = padj)) +
    geom_col() +
    theme_bw() +
    labs(
      title = paste("Differential pathways:", cname),
      x     = "log2 fold change",
      y     = "MetaCyc pathway"
    ) +
    scale_fill_gradient(
      name  = "FDR (padj)",
      low   = "red",
      high  = "grey80",
      trans = "reverse"
    ) +
    theme(
      axis.text.y  = element_text(size = 5),   # smaller pathway labels
      axis.text.x  = element_text(size = 8),
      axis.title   = element_text(size = 9),
      plot.title   = element_text(size = 10, face = "bold"),
      legend.title = element_text(size = 8),
      legend.text  = element_text(size = 7)
    )
  
  print(p)
  
  # Safe filename (remove weird symbols)
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", cname)
  
  ggsave(
    filename = paste0("Fig_DESeq2_pathways_", safe_name, ".png"),
    plot     = p,
    width    = 6,
    height   = 7,
    dpi      = 300
  )
}

### 18. Prevalence filter + TSS+CLR transform for exploratory plots (3F) ----

# counts matrix: all pathways (rows) × samples (cols)
counts_mat_all <- abundance_data_filtered %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# 18a. Prevalence filter: keep pathways present in ≥10% of samples
prevalence_vec <- rowMeans(counts_mat_all > 0)
keep_rows      <- prevalence_vec >= 0.10
counts_mat_pf  <- counts_mat_all[keep_rows, , drop = FALSE]

# 18b. Total-sum scaling (TSS)
col_totals   <- colSums(counts_mat_pf)
rel_abun     <- sweep(counts_mat_pf, 2, col_totals, "/")

# 18c. CLR transform (add small pseudocount to avoid log(0))
rel_abun[rel_abun == 0] <- 1e-6
log_rel  <- log(rel_abun)
clr_mat  <- sweep(log_rel, 2, colMeans(log_rel), "-")  # center per sample

### 18.d CLR-based heatmap on top 25 most variable pathways ----

# choose how many pathways to show
top_n <- 25

# compute variance for each pathway across samples
row_var <- matrixStats::rowVars(as.matrix(clr_mat))

# get indices for top N most variable pathways
keep_idx <- order(row_var, decreasing = TRUE)[1:top_n]

# subset CLR matrix
clr_mat_top <- clr_mat[keep_idx, , drop = FALSE]

# draw heatmap grouped by severity_sex
heat_clr <- pathway_heatmap(
  abundance = clr_mat_top,
  metadata  = metadata,
  group     = "severity_sex"
)

heat_clr <- heat_clr +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 5),  # smaller y-axis labels
    plot.margin  = margin(5.5, 40, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = "Fig1b_pathway_heatmap_CLR_severity_sex.png",
  plot     = heat_clr,
  width    = 10,
  height   = 6,
  dpi      = 300
)

# 18e. CLR-based heatmap using ggpicrust2
heat_clr <- pathway_heatmap(
  abundance = clr_mat,
  metadata  = metadata,
  group     = "severity_sex"
)

heat_clr <- heat_clr +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 8),
    plot.margin  = margin(5.5, 40, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = "Fig1c_pathway_heatmap_CLR_severity_sex.png",
  plot     = heat_clr,
  width    = 14,
  height   = 7,
  dpi      = 300
)


### 19. Sex effects within Mild+ and Severe+ COVID groups (3G) ----

run_sex_within_severity <- function(severity_label) {
  meta_sub <- metadata %>%
    dplyr::filter(severity == severity_label)
  
  message("Running DESeq2 for sex within severity = ", severity_label,
          " (n = ", nrow(meta_sub), " samples)")
  
  # counts for those samples only
  counts_sub <- abundance_data_filtered[, c("pathway", meta_sub$`sample-id`)]
  mat_sub    <- counts_sub %>%
    tibble::column_to_rownames("pathway") %>%
    as.matrix()
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = round(mat_sub),
    colData   = meta_sub,
    design    = ~ sex_short
  )
  
  dds <- BiocGenerics::estimateSizeFactors(dds, type = "poscounts")
  dds <- DESeq2::DESeq(dds)
  
  # contrast: Male vs Female
  res <- as.data.frame(
    DESeq2::results(dds, contrast = c("sex_short", "M", "F"))
  )
  
  res$feature  <- rownames(res)
  res$severity <- severity_label
  res
}

sex_mildplus  <- run_sex_within_severity("Mild+")
sex_severeplus <- run_sex_within_severity("Severe+")

sex_within_severity <- dplyr::bind_rows(sex_mildplus, sex_severeplus)

sex_within_severity_annot <- sex_within_severity %>%
  dplyr::inner_join(metacyc_daa_annotated_results_df, by = "feature") %>%
  dplyr::select(
    feature,
    description,
    severity,
    baseMean,
    log2FoldChange,
    lfcSE,
    stat,
    pvalue,
    padj
  )

readr::write_tsv(
  sex_within_severity_annot,
  "Table_DESeq2_sex_within_MildPlus_SeverePlus.tsv"
)


### 20. Volcano plots for key sex contrasts (3H) ----

volcano_targets <- c(
  "Asymp_M_vs_Asymp_F",
  "Mild+_M_vs_Mild+_F",
  "Severe+_M_vs_Severe+_F",
  "Deceased_M_vs_Deceased_F"
)

existing_contrasts <- intersect(volcano_targets,
                                unique(deseq_res_annot$contrast))

for (cname in existing_contrasts) {
  vdat <- deseq_res_annot %>%
    dplyr::filter(contrast == cname, !is.na(padj))
  
  if (nrow(vdat) == 0) next
  
  vdat <- vdat %>%
    dplyr::mutate(
      neg_log10_padj = -log10(padj),
      sig            = padj < 0.05 & abs(log2FoldChange) >= 0.5
    )
  
  p_volcano <- ggplot(vdat,
                      aes(x = log2FoldChange,
                          y = neg_log10_padj,
                          colour = sig)) +
    geom_point(alpha = 0.7, size = 1.5) +
    scale_colour_manual(
      values = c("FALSE" = "grey70", "TRUE" = "red"),
      name   = "Significant\n(padj < 0.05,\n|log2FC| ≥ 0.5)"
    ) +
    theme_bw() +
    labs(
      title = paste("Volcano plot:", cname),
      x     = "log2 fold change (M vs F)",
      y     = expression(-log[10](padj))
    )
  
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", cname)
  
  ggsave(
    filename = paste0("Fig_DESeq2_volcano_", safe_name, ".png"),
    plot     = p_volcano,
    width    = 7,
    height   = 6,
    dpi      = 300
  )
}


### 21. Top-10 up- and down-regulated pathway tables per contrast (3J) ----

all_contrasts_full <- unique(deseq_res_annot$contrast)

for (cname in all_contrasts_full) {
  cdat <- deseq_res_annot %>%
    dplyr::filter(contrast == cname, !is.na(padj))
  
  if (nrow(cdat) == 0) next
  
  top_up <- cdat %>%
    dplyr::arrange(dplyr::desc(log2FoldChange)) %>%
    dplyr::slice_head(n = 10)
  
  top_down <- cdat %>%
    dplyr::arrange(log2FoldChange) %>%
    dplyr::slice_head(n = 10)
  
  safe_name <- gsub("[^A-Za-z0-9_]+", "_", cname)
  
  readr::write_tsv(
    top_up,
    paste0("Table_DESeq2_top10_up_", safe_name, ".tsv")
  )
  
  readr::write_tsv(
    top_down,
    paste0("Table_DESeq2_top10_down_", safe_name, ".tsv")
  )
}


### 22. OPTIONAL: Pathway–taxa mapping hook for Aim 2 (3I) ----
# This will ONLY run if you have already created an object
# `aim2_taxa_pathway_map` with a column `feature` (MetaCyc ID) and
# one or more taxonomic columns (e.g. genus, species, KO, etc.).

if (exists("aim2_taxa_pathway_map")) {
  pathway_taxa_sig <- sig_res %>%
    dplyr::select(feature, description, contrast,
                  log2FoldChange, padj) %>%
    dplyr::left_join(aim2_taxa_pathway_map, by = "feature")
  
  readr::write_tsv(
    pathway_taxa_sig,
    "Table_Pathway_Taxa_mapping_sig_pathways.tsv"
  )
} else {
  message("Note: to generate the Aim 2 ↔ Aim 3 pathway–taxa mapping, ",
          "first load a data frame named `aim2_taxa_pathway_map` ",
          "with a `feature` column matching MetaCyc IDs, then re-run section 22.")
}




