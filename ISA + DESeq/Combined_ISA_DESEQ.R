# --- 1. Load libraries ---
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tibble)
library(forcats)


# --- 2. Get OTUs significant in ISA and DESeq2 ---


isa_otus <- significant_indicators$OTU 

deseq_sig_otus <- deseq_all %>%
  filter(!is.na(padj) & padj <= 0.05) %>%
  pull(OTU) %>%
  unique()

both_otus <- intersect(isa_otus, deseq_sig_otus)

# Prune phyloseq to only these significant OTUs
ps_both <- prune_taxa(both_otus, ps)

# --- 3. Clean metadata and define Group_combined ---
meta_df <- data.frame(sample_data(ps_both))

# Remove any samples with missing sex or Severity
meta_df <- meta_df %>% filter(!is.na(sex) & !is.na(Severity))

# Build Group_combined
meta_df <- meta_df %>% mutate(Group_combined = paste(sex, Severity, sep = "_"))

# Prune samples in phyloseq to match cleaned metadata
ps_both <- prune_samples(rownames(meta_df), ps_both)
sample_data(ps_both) <- sample_data(meta_df)

# --- 4. Clean taxonomy (Genus level) ---
tax_df <- as.data.frame(tax_table(ps_both)) %>%
  rownames_to_column("OTU") %>%
  # Keep only valid genera (remove unclassified)
  filter(!is.na(Genus) & Genus != "" & Genus != "g__") %>%
  mutate(Label = gsub("^g__", "", Genus))

# Keep only OTUs with valid genera
ps_both <- prune_taxa(tax_df$OTU, ps_both)

# --- 5. Melt phyloseq and join genus labels ---
ps_df <- psmelt(ps_both) %>%
  left_join(tax_df %>% select(OTU, Label), by = "OTU") %>%
  filter(!is.na(Label) & !is.na(Group_combined))

# --- 5b. FILTER: Manually exclude known high-abundance commensals ---

# Define the bulk commensal genera to exclude based on literature/high baseline abundance.
commensal_genera_to_exclude <- c("Veillonella", "Prevotella_7")

ps_df_filtered <- ps_df %>%
  filter(!Label %in% commensal_genera_to_exclude)
# NOTE: If your plot still looks dominated by one genus, add that genus's name here.

# --- 6. Summarize ABSOLUTE abundance per sample ---
plot_data_abs <- ps_df_filtered %>% # **** Using the filtered object ****
  group_by(Sample, Group_combined, Label) %>%
  summarize(sample_abundance = sum(Abundance), .groups = "drop") %>%
  ungroup() %>%
  # Summarize the mean ABSOLUTE count per group and genus
  group_by(Group_combined, Label) %>%
  summarize(mean_absolute_abundance = mean(sample_abundance), .groups = "drop") %>%
  # Add a pseudocount of 1 to handle zero values for the Log transformation
  mutate(log10_absolute_abundance = log10(mean_absolute_abundance + 1))


# --- 7. Set x-axis order ---
group_levels <- c(
  "male_Control", "female_Control",
  "male_Mild", "female_Mild",
  "male_Severe", "female_Severe",
  "male_Deceased", "female_Deceased"
)
plot_data_abs$Group_combined <- factor(plot_data_abs$Group_combined, levels = group_levels)

# --- 8. Color palette ---
unique_genera <- sort(unique(plot_data_abs$Label))
palette <- grDevices::colorRampPalette(c( "lightgreen", "#0072B8", "#009E73", "turquoise", "#56B4E8", "#3366CC", "lightblue"))(length(unique_genera))




# --- 9. Plot Log10 Stacked Bar Chart  ---
log_stacked_plot <- ggplot(plot_data_abs, aes(
  x = Group_combined,
  # Use the log10-transformed mean absolute abundance for the stack height
  y = log10_absolute_abundance,
  fill = Label
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  labs(
    title = "Absolute Abundance of Significant, Non-Commensal Genera (Log10 Stack)",
    x = "Sex–Severity Group",
    # Y-axis label reflects the transformation
    y = "Mean Absolute Abundance (Log10(Counts + 1))",
    fill = "Genus"
  ) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
  )

print(log_stacked_plot)
ggsave("stacked_genus_log10_absolute_abundance_non_commensal.png", log_stacked_plot, width = 12, height = 9, dpi = 300)

