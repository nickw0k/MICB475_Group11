### Paths (relative to your Git Project root)
abundance_file <- "Mexico Dataset QIIME2 files/path_abun_unstrat.tsv.gz"
metadata_file  <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Create a list of all the packages you need to install 
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2", "vegan", "ashr") 

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
library(vegan)
library(ashr)

#Read PICRUSt2 unstratified pathway abundance table
abundance_data <- read_tsv(
  abundance_file,
  col_names = TRUE,
  trim_ws   = TRUE
) |> as.data.frame()

#Read metadata
metadata <- read_tsv(metadata_file)

#Build severity + severity_sex groups and drop MildNegative #
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
    # >>> key change: propagate NA instead of making a "NA_*" string
    severity_sex = dplyr::if_else(
      is.na(severity) | is.na(sex_short),
      NA_character_,
      paste0(severity, "_", sex_short)
    )
  ) |>
  # drop Mild negatives and any rows without clear severity/sex
  dplyr::filter(severity != "Mild_neg") |>
  dplyr::filter(!is.na(severity_sex))

metadata$severity_sex <- factor(
  metadata$severity_sex,
  levels = c("Asymp_F","Asymp_M",
             "Mild+_F","Mild+_M",
             "Severe+_F","Severe+_M",
             "Deceased_F","Deceased_M")
)

table(metadata$severity_sex)

#Get sample IDs that remain after filtering #
sample_names <- metadata$`sample-id`

#Keep only 'pathway' + those samples, in the same order #
keep_cols <- c("pathway", sample_names)
abundance_data_filtered <- abundance_data[, colnames(abundance_data) %in% keep_cols]
abundance_data_filtered <- abundance_data_filtered[, keep_cols]

# Drop samples that are all zero (no counts) #
non_zero <- colSums(abundance_data_filtered[, -1] != 0) > 0
abundance_data_filtered <- abundance_data_filtered[, c(TRUE, non_zero)]

#Align metadata rows and order with filtered abundance table#
abun_samples <- colnames(abundance_data_filtered)[-1]
metadata <- metadata[metadata$`sample-id` %in% abun_samples, ]
metadata <- metadata[match(abun_samples, metadata$`sample-id`), ]
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$`sample-id`
metadata$sample_name <- metadata$`sample-id`

# Differential abundance testing at pathway level (DESeq2 on severity_sex) #
abundance_daa_results_df <- pathway_daa(
  abundance  = abundance_data_filtered %>% tibble::column_to_rownames("pathway"),
  metadata   = metadata,
  group      = "severity_sex",
  daa_method = "DESeq2"
)

# Add MetaCyc pathway descriptions #
metacyc_daa_annotated_results_df <- pathway_annotation(
  pathway        = "MetaCyc",
  daa_results_df = abundance_daa_results_df,
  ko_to_kegg     = FALSE
)

#PCA on all pathways

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
  "Figure2.png",
  plot   = pca_severity,
  width  = 8, height = 6, dpi = 300
)

#Sex-specific DESeq2 plots (F vs M within each severity) -------------
# Uses: abundance_data_filtered, metadata, metacyc_daa_annotated_results_df

#Aim 3 QC: group sample sizes  ----
# Make a dplyr-friendly copy of metadata 
metadata_df <- as.data.frame(metadata)

n_sev_sex <- metadata_df %>%
  dplyr::count(severity, sex_short, name = "n") %>%
  dplyr::arrange(severity, sex_short)

n_severity_sex <- metadata_df %>%
  dplyr::count(severity_sex, name = "n") %>%
  dplyr::arrange(severity_sex)

write.csv(n_sev_sex, "Table_n_bySeverity_bySex.csv", row.names = FALSE)
write.csv(n_severity_sex, "Table_n_bySeveritySex.csv", row.names = FALSE)

print(n_sev_sex)
print(n_severity_sex)

# Build DESeq2 object using severity_sex as the design 

# counts matrix: pathways x samples
count_mat_sex <- abundance_data_filtered %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# ensure sample order matches metadata
count_mat_sex <- count_mat_sex[, metadata$`sample-id`, drop = FALSE]
stopifnot(identical(colnames(count_mat_sex), metadata$`sample-id`))

# Prefilter low-abundance / ultra-rare pathways (global)
prev_sex <- rowSums(count_mat_sex > 0) / ncol(count_mat_sex)
tot_sex  <- rowSums(count_mat_sex)

keep_rows <- (prev_sex >= 0.10) & (tot_sex >= 10)
count_mat_sex <- count_mat_sex[keep_rows, , drop = FALSE]

dds_sex <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(count_mat_sex),
  colData   = metadata,
  design    = ~ severity_sex
)

# enforce group order
dds_sex$severity_sex <- factor(dds_sex$severity_sex, levels = levels(metadata$severity_sex))
stopifnot(all(levels(dds_sex$severity_sex) == levels(metadata$severity_sex)))

# size factors (poscounts handles sparse)
dds_sex <- BiocGenerics::estimateSizeFactors(dds_sex, type = "poscounts")
dds_sex <- DESeq2::DESeq(dds_sex)

# Sex contrasts within each severity: shrunken LFC + FDR-filtered plots 
if (!requireNamespace("ashr", quietly = TRUE)) install.packages("ashr")

alpha <- 0.05

# Plot/count “universe” filter for sex-within-severity contrasts (contrast-specific)
min_meanNorm_sev <- 10   # within-severity mean normalized count threshold
min_prev_sev     <- 0    # set to 0.10 if you want a within-severity prevalence filter too

top_n <- 20

# MetaCyc description map (for nicer y-axis labels)
# Expecting metacyc_daa_annotated_results_df has columns: feature, description
name_map <- metacyc_daa_annotated_results_df %>%
  dplyr::select(feature, description) %>%
  dplyr::distinct() %>%
  dplyr::rename(pathway = feature)

# Precompute normalized counts for within-severity mean/prevalence filters
norm_mat <- DESeq2::counts(dds_sex, normalized = TRUE)

# must exist in metadata/colData; used to select columns for each severity
sev_vec <- as.character(SummarizedExperiment::colData(dds_sex)$severity)

# helper: get DESeq2 results + shrink log2FC for an arbitrary contrast
get_res_shrunk <- function(dds, numerator, denominator, alpha = 0.05) {
  res <- DESeq2::results(dds, contrast = c("severity_sex", numerator, denominator), alpha = alpha)
  res_s <- DESeq2::lfcShrink(
    dds,
    contrast = c("severity_sex", numerator, denominator),
    res = res,
    type = "ashr"
  )
  res_tbl <- as.data.frame(res_s)
  res_tbl$pathway <- rownames(res_tbl)
  res_tbl
}

# helper: annotate + plot sig-only barplot (padj < alpha AND meanNorm_sev >= min_meanNorm_sev)
plot_sex_bar <- function(res_tbl, title_stub, out_png,
                         alpha = 0.05, top_n = 20,
                         min_meanNorm_sev = 10, min_prev_sev = 0,
                         name_map = NULL) {
  
  # annotate labels if map provided
  if (!is.null(name_map)) {
    res_tbl <- res_tbl %>%
      left_join(name_map, by = "pathway") %>%
      mutate(label = if_else(is.na(description) | description == "", pathway, description))
  } else {
    res_tbl$label <- res_tbl$pathway
  }
  
  res_f <- res_tbl %>%
    filter(!is.na(padj)) %>%
    filter(meanNorm_sev >= min_meanNorm_sev) %>%
    filter(prev_sev >= min_prev_sev) %>%
    filter(padj < alpha) %>%
    mutate(
      direction = ifelse(log2FoldChange > 0, "Enriched in females", "Enriched in males"),
      direction = factor(direction, levels = c("Enriched in females", "Enriched in males"))
    ) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n) %>%
    mutate(label = stringr::str_wrap(label, width = 55))
  
  if (nrow(res_f) == 0) return(NULL)
  
  p <- ggplot(res_f, aes(x = log2FoldChange, y = reorder(label, log2FoldChange), fill = direction)) +
    geom_col() +
    scale_fill_manual(
      values = c("Enriched in females" = "#E91E63", "Enriched in males" = "#1E88E5"),
      drop = FALSE
    ) +
    labs(
      title = paste0(title_stub, " (shrunken log2FC; FDR < ", alpha,
                     ", meanNorm(sev) ≥ ", min_meanNorm_sev,
                     ifelse(min_prev_sev > 0, paste0(", prev(sev) ≥ ", min_prev_sev), ""),
                     ")"),
      x = "shrunken log2 fold change (F vs M; >0 enriched in females)",
      y = "MetaCyc pathway",
      fill = ""
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(out_png, p, width = 10, height = 7, dpi = 300)
  p
}

# helper: annotate + plot sig-only DOT plot with 95% CI (requires lfcSE)
plot_sex_dot <- function(res_tbl, title_stub, out_png,
                         alpha = 0.05, top_n = 20,
                         min_meanNorm_sev = 10, min_prev_sev = 0,
                         name_map = NULL) {
  
  # annotate labels if map provided
  if (!is.null(name_map)) {
    res_tbl <- res_tbl %>%
      left_join(name_map, by = "pathway") %>%
      mutate(label = if_else(is.na(description) | description == "", pathway, description))
  } else {
    res_tbl$label <- res_tbl$pathway
  }
  
  df <- res_tbl %>%
    filter(!is.na(padj), padj < alpha) %>%
    filter(meanNorm_sev >= min_meanNorm_sev) %>%
    filter(prev_sev >= min_prev_sev) %>%
    mutate(
      direction = ifelse(log2FoldChange > 0, "Enriched in females", "Enriched in males"),
      direction = factor(direction, levels = c("Enriched in females", "Enriched in males")),
      label = stringr::str_wrap(label, width = 55)
    ) %>%
    arrange(desc(abs(log2FoldChange))) %>%
    slice_head(n = top_n)
  
  if (nrow(df) == 0) return(NULL)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = reorder(label, log2FoldChange), color = direction)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(
      aes(xmin = log2FoldChange - 1.96 * lfcSE,
          xmax = log2FoldChange + 1.96 * lfcSE),
      height = 0.2
    ) +
    geom_point(size = 2.2) +
    scale_color_manual(values = c("Enriched in females" = "#E91E63",
                                  "Enriched in males"   = "#1E88E5"),
                       drop = FALSE) +
    labs(
      title = paste0(title_stub, " (shrunken log2FC ±95% CI; FDR < ", alpha,
                     ", meanNorm(sev) ≥ ", min_meanNorm_sev,
                     ifelse(min_prev_sev > 0, paste0(", prev(sev) ≥ ", min_prev_sev), ""),
                     ")"),
      x = "shrunken log2 fold change (F vs M; >0 enriched in females)",
      y = "MetaCyc pathway",
      color = ""
    ) +
    theme_bw() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(out_png, p, width = 10, height = 7, dpi = 300)
  p
}

# Define your 4 within-severity sex contrasts (MUST match levels(metadata$severity_sex))
sex_contrasts <- tibble::tribble(
  ~label,     ~num,           ~den,
  "Asymp",    "Asymp_F",      "Asymp_M",
  "Mild+",    "Mild+_F",      "Mild+_M",
  "Severe+",  "Severe+_F",    "Severe+_M",
  "Deceased", "Deceased_F",   "Deceased_M"
)

# sanity check: contrasts exist
missing_groups <- setdiff(c(sex_contrasts$num, sex_contrasts$den), levels(dds_sex$severity_sex))
if (length(missing_groups) > 0) {
  stop("These severity_sex levels are missing from metadata/dds_sex: ", paste(missing_groups, collapse = ", "))
}

all_res <- list()

for (i in seq_len(nrow(sex_contrasts))) {
  lab <- sex_contrasts$label[i]
  num <- sex_contrasts$num[i]
  den <- sex_contrasts$den[i]
  
  res_tbl <- get_res_shrunk(dds_sex, num, den, alpha = alpha)
  
  # Add contrast-specific filter variables (within this severity only)
  cols_sev <- which(sev_vec == lab)
  if (length(cols_sev) == 0) {
    stop("No samples matched severity label '", lab, "' in colData(dds_sex)$severity. Check spelling/levels.")
  }
  
  nm <- norm_mat[res_tbl$pathway, cols_sev, drop = FALSE]
  res_tbl$meanNorm_sev <- rowMeans(nm)
  res_tbl$prev_sev     <- rowMeans(nm > 0)
  
  # Save full shrunken table (includes non-sig but stabilized log2FC)
  out_full <- paste0("DESeq2_", gsub("[^A-Za-z0-9_]+", "_", lab), "_F_vs_M_shrunken_FULL.csv")
  write.csv(res_tbl, file = out_full, row.names = FALSE)
  
  # Save sig-only table (MATCHES plot filters)
  res_sig <- res_tbl %>%
    filter(!is.na(padj), padj < alpha, meanNorm_sev >= min_meanNorm_sev, prev_sev >= min_prev_sev) %>%
    left_join(name_map, by = "pathway") %>%
    mutate(label = if_else(is.na(description) | description == "", pathway, description))
  
  out_sig <- paste0("DESeq2_", gsub("[^A-Za-z0-9_]+", "_", lab), "_F_vs_M_shrunken_FDR", alpha, ".csv")
  write.csv(res_sig, file = out_sig, row.names = FALSE)
  
  # Plot sig-only barplot (MATCHES the same filters used in res_sig)
  plot_sex_bar(
    res_tbl,
    title_stub        = paste0("Differential pathways: ", lab, " (F vs M)"),
    out_png           = paste0("Fig_DESeq2_", gsub("[^A-Za-z0-9_]+", "_", lab), "_F_vs_M_shrunken_FDR", alpha, ".png"),
    alpha             = alpha,
    top_n             = top_n,
    min_meanNorm_sev  = min_meanNorm_sev,
    min_prev_sev      = min_prev_sev,
    name_map          = name_map
  )
  
  # Plot sig-only dotplot (MATCHES the same filters used in res_sig)
  plot_sex_dot(
    res_tbl,
    title_stub        = paste0("Differential pathways: ", lab, " (F vs M)"),
    out_png           = paste0("Fig_DESeq2_", gsub("[^A-Za-z0-9_]+", "_", lab), "_F_vs_M_shrunken_DOT_FDR", alpha, ".png"),
    alpha             = alpha,
    top_n             = top_n,
    min_meanNorm_sev  = min_meanNorm_sev,
    min_prev_sev      = min_prev_sev,
    name_map          = name_map
  )
  
  all_res[[lab]] <- res_tbl
}

# Sig-count table derived from the SAME 17b shrunken results (recommended) ----
# Counts are among pathways that:
# (i) have padj (not NA) AND (ii) pass meanNorm_sev >= min_meanNorm_sev (and prev_sev >= min_prev_sev)
# so n_total/n_sig/prop_sig matches the plotted universe.

sex_sig_counts <- purrr::pmap_dfr(
  sex_contrasts,
  function(label, num, den) {
    res_tbl <- all_res[[label]]
    
    tested <- !is.na(res_tbl$padj) &
      (res_tbl$meanNorm_sev >= min_meanNorm_sev) &
      (res_tbl$prev_sev >= min_prev_sev)
    
    sig <- tested & (res_tbl$padj < alpha)
    
    tibble::tibble(
      contrast = paste0(num, "_vs_", den),
      severity = label,
      n_total  = sum(tested),
      n_sig    = sum(sig),
      prop_sig = ifelse(sum(tested) > 0, sum(sig) / sum(tested), NA_real_)
    )
  }
)

print(sex_sig_counts)
readr::write_tsv(sex_sig_counts, "Table_DESeq2_sex_sig_pathway_counts_FvM.tsv")

# Prevalence filter + TSS + CLR transformation (for severity-level plots) 
# counts matrix: pathways x samples
count_mat <- abundance_data_filtered %>%
  tibble::column_to_rownames("pathway") %>%
  as.matrix()

# ensure sample order matches metadata
count_mat <- count_mat[, metadata$`sample-id`, drop = FALSE]

# prevalence (fraction of samples where count > 0) and total counts
prev       <- rowSums(count_mat > 0) / ncol(count_mat)
tot_counts <- rowSums(count_mat)

# keep pathways present in ≥10% of samples and with total counts ≥10
keep_idx       <- (prev >= 0.10) & (tot_counts >= 10)
count_mat_filt <- count_mat[keep_idx, , drop = FALSE]

# total-sum scaling (TSS) → relative abundances
col_totals <- colSums(count_mat_filt)
rel_abun   <- sweep(count_mat_filt, 2, col_totals, "/")

# centred log-ratio (CLR) per sample
pseudo   <- 1e-6
log_rel  <- log(rel_abun + pseudo)
geo_mean <- colMeans(log_rel)                      # geometric mean per sample (on log scale)
clr_mat  <- sweep(log_rel, 2, geo_mean, "-")       # pathways x samples

# make sure severity is an ordered factor
metadata$severity <- factor(
  metadata$severity,
  levels = c("Asymp", "Mild+", "Severe+", "Deceased")
)

# keep CLR columns in metadata order
clr_mat <- clr_mat[, metadata$`sample-id`, drop = FALSE]

# PERMANOVA on CLR-transformed pathways (severity vs sex) #

# Safety check: columns of clr_mat must match metadata sample order
stopifnot(all(colnames(clr_mat) == metadata$`sample-id`))

# Make sure severity and sex_short are factors
metadata$severity  <- factor(
  metadata$severity,
  levels = c("Asymp", "Mild+", "Severe+", "Deceased")
)
metadata$sex_short <- factor(metadata$sex_short, levels = c("F", "M"))

# Distance matrix on CLR-transformed pathways
# (Euclidean on CLR is equivalent to Aitchison distance for compositional data)
clr_dist <- vegdist(t(clr_mat), method = "euclidean")

# PERMANOVA without interaction: "how much variance is explained by severity vs sex?"
permanova_main <- adonis2(
  clr_dist ~ severity + sex_short,
  data = metadata
)

# PERMANOVA with interaction term (optional, but nice to report)
permanova_int <- adonis2(
  clr_dist ~ severity * sex_short,
  data = metadata
)

# Save results to text files for the paper / supplement
capture.output(permanova_main,
               file = "PERMANOVA_CLR_severity_sex_main.txt")
capture.output(permanova_int,
               file = "PERMANOVA_CLR_severity_sex_interaction.txt")

# Also print to console when you run the script
print(permanova_main)
print(permanova_int)

# PERMANOVA term-level effects + dispersion diagnostics 

# Term-level variance explained (severity vs sex)
permanova_main_terms <- adonis2(clr_dist ~ severity + sex_short,
                                data = metadata, by = "margin")
permanova_int_terms  <- adonis2(clr_dist ~ severity * sex_short,
                                data = metadata, by = "margin")

print(permanova_main_terms)
print(permanova_int_terms)

capture.output(
  permanova_main_terms,
  file = "PERMANOVA_CLR_severity_sex_main_byMargin.txt"
)

capture.output(
  permanova_int_terms,
  file = "PERMANOVA_CLR_severity_sex_interaction_byMargin.txt"
)

# Check homogeneity of dispersion (PERMANOVA assumption)
bd_severity <- betadisper(clr_dist, metadata$severity)
bd_sex      <- betadisper(clr_dist, metadata$sex_short)

print(permutest(bd_severity, permutations = 999))
print(permutest(bd_sex,      permutations = 999))

capture.output(
  permutest(bd_severity, permutations = 999),
  file = "PERMDISP_CLR_bySeverity.txt"
)

capture.output(
  permutest(bd_sex, permutations = 999),
  file = "PERMDISP_CLR_bySex.txt"
)

# Visualize and save dispersion diagnostics 
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2)
})

# Extract distances to centroid for plotting
disp_severity <- data.frame(
  severity = metadata$severity,
  dist_to_centroid = bd_severity$distances
)

p_disp <- ggplot(disp_severity, aes(x = severity, y = dist_to_centroid)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  theme_bw() +
  labs(title = "PERMDISP: Distances to centroid by severity", x = "Severity", y = "Distance to centroid")

ggsave("Fig_PERMDISP_bySeverity.png", p_disp, width = 9, height = 5, dpi = 300)


# PCA of CLR-transformed pathways by severity #

# PCA on samples (transpose: samples x pathways)
clr_pca <- prcomp(t(clr_mat), center = TRUE, scale. = FALSE)

# variance explained
var_expl <- (clr_pca$sdev^2) / sum(clr_pca$sdev^2)
pc1_var <- round(100 * var_expl[1], 1)
pc2_var <- round(100 * var_expl[2], 1)

pca_df <- as.data.frame(clr_pca$x[, 1:2])
pca_df$`sample-id` <- rownames(pca_df)

pca_df <- dplyr::left_join(
  pca_df,
  metadata[, c("sample-id", "severity", "sex_short")],
  by = "sample-id"
)

FigA_pca_severity <- ggplot(
  pca_df,
  aes(
    x      = PC1,
    y      = PC2,
    colour = severity,
    shape  = sex_short
  )
) +
  # NEW: ellipses per severity group
  stat_ellipse(
    aes(group = severity),
    type       = "norm",
    level      = 0.68,      
    linewidth  = 0.7,
    alpha      = 0.6,
    show.legend = FALSE
  ) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_shape_manual(values = c(F = 16, M = 17), na.translate = FALSE) +
  theme_bw() +
  labs(
    title  = "PCA of CLR-transformed pathways",
    x      = paste0("PC1 (", pc1_var, "%)"),
    y      = paste0("PC2 (", pc2_var, "%)"),
    colour = "Severity",
    shape  = "Sex"
  ) +
  theme(
    legend.position = "right",
    axis.text       = element_text(size = 10),
    axis.title      = element_text(size = 11),
    plot.title      = element_text(size = 12, face = "bold")
  )

# Add inference annotation box (PERMANOVA + PERMDISP) 
# Use the by-margin model (already computed):
# permanova_main_terms <- adonis2(..., by="margin")
# and the PERMDISP permutest output:
# permutest(bd_severity, permutations=999)

permanova_terms <- permanova_main_terms

# Extract PERMDISP p-value for severity
permdisp_sev_p <- permutest(bd_severity, permutations = 999)$tab[1, "Pr(>F)"]

lab <- paste0(
  "PERMANOVA (CLR-Aitchison; by=margin)\n",
  "severity: R²=", sprintf("%.3f", permanova_terms["severity", "R2"]),
  ", p=", signif(permanova_terms["severity", "Pr(>F)"], 2), "\n",
  "sex: R²=", sprintf("%.3f", permanova_terms["sex_short", "R2"]),
  ", p=", signif(permanova_terms["sex_short", "Pr(>F)"], 2), "\n",
  "PERMDISP (severity): p=", signif(permdisp_sev_p, 2)
)

FigA_pca_severity <- FigA_pca_severity +
  labs(caption = lab) +
  theme(
    plot.caption = element_text(hjust = 1, size = 9),
    plot.caption.position = "plot",
    plot.margin = margin(10, 10, 18, 10)  # extra bottom space for caption
  )

ggsave(
  "Figure3.png",
  plot   = FigA_pca_severity,
  width  = 6.5,
  height = 5,
  dpi    = 300
)

# CLR heatmap of top severity-associated pathways 

sev_levels <- levels(metadata$severity)

# mean CLR per severity group
mean_by_sev <- sapply(
  sev_levels,
  function(s) {
    samp <- metadata$`sample-id`[metadata$severity == s]
    rowMeans(clr_mat[, samp, drop = FALSE])
  }
)

# difference between Deceased and Asymp to rank pathways
diff_D_vs_A <- mean_by_sev[, "Deceased"] - mean_by_sev[, "Asymp"]

topN   <- min(30, nrow(clr_mat))
top_idx <- order(abs(diff_D_vs_A), decreasing = TRUE)[1:topN]
clr_mat_top <- clr_mat[top_idx, , drop = FALSE]

# order samples by severity (within-severity keep original order)
sample_order     <- order(metadata$severity)
metadata_ordered <- metadata[sample_order, ]
clr_mat_top      <- clr_mat_top[, metadata_ordered$`sample-id`, drop = FALSE]

heat_B <- pathway_heatmap(
  abundance = clr_mat_top,
  metadata  = metadata_ordered,
  group     = "severity"
)

heat_B <- heat_B +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y  = element_text(size = 6),
    plot.margin  = margin(5.5, 40, 5.5, 5.5)
  ) +
  coord_cartesian(clip = "off")

ggsave(
  filename = "Figure4.png",
  plot     = heat_B,
  width    = 7,
  height   = 6,
  dpi      = 300
)
