##### Install packages #####
# Start by installing all necessary packages when asked if you want to install
# from source, please just type Yes in the terminal below

### Paths (relative to your Git Project root)
abundance_file <- "Mexico Dataset QIIME2 files/path_abun_unstrat.tsv.gz"
metadata_file  <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Create a list of all the packages you need to install 
pkgs <- c("ALDEx2", "SummarizedExperiment", "Biobase", "devtools", "ComplexHeatmap", "BiocGenerics", "metagenomeSeq", "Maaslin2", "edgeR", "lefser", "limma", "KEGGREST", "DESeq2") 

# Use the above list to install all the packages using a for loop 
for (pkg in pkgs) { if (!requireNamespace(pkg, quietly = TRUE)) BiocManager::install(pkg) } 

# when asked if you want to update all, some or none please type "n" for none 

# After installing all of its above dependencies, install ggpicrust2 
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
    sex = dplyr::if_else(sex %in% c("male", "female"), sex, NA_character_)
  ) |>
  # drop Mild negatives and any rows without clear severity/sex
  dplyr::filter(severity != "Mild_neg") |>
  dplyr::filter(!is.na(severity), !is.na(sex)) |>
  dplyr::mutate(
    severity_sex = paste0(severity, "_", substr(sex, 1, 1))
  )

# keep only samples that are in the abundance table and preserve order
abun_samples <- colnames(abundance_data_filtered)[-1]
metadata <- metadata[metadata$`sample-id` %in% abun_samples, ]
metadata <- metadata[match(abun_samples, metadata$`sample-id`), ]
rownames(metadata) <- metadata$`sample-id`

# lock factor order (also forces Deceased F/M to be kept if present)
metadata$severity_sex <- factor(
  metadata$severity_sex,
  levels = c("Asymp_F","Asymp_M",
             "Mild+_F","Mild+_M",
             "Severe+_F","Severe+_M",
             "D_F","D_M")
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

# If DESeq2_function.R lives in Git Project:
source("DESeq2_function.R")

res <- DEseq2_function(abundance_data_filtered, metadata, "severity_sex")
res$feature <- rownames(res)

res_desc <- dplyr::inner_join(
  res,
  metacyc_daa_annotated_results_df,
  by = "feature"
)

# Keep only needed columns
res_desc <- res_desc[, -c(8:13)]

# Filter significant and plot
sig_res <- res_desc %>% dplyr::filter(pvalue < 0.05)
sig_res <- sig_res[order(sig_res$log2FoldChange), ]

ggplot(sig_res, aes(
  y    = reorder(description, sort(as.numeric(log2FoldChange))),
  x    = log2FoldChange,
  fill = pvalue
)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(x = "Log2FoldChange", y = "Pathways")


