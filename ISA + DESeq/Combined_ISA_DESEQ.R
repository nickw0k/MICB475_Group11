# --- 1. Load libraries ---
library(dplyr)
library(ggplot2)
library(phyloseq)
library(tidyr)
library(tibble)
library(forcats)

# --- 2. Get OTUs significant in BOTH ISA and DESeq2 ---
isa_otus <- significant_indicators$OTU 
deseq_otus <- deseq_all$OTU 
both_otus <- intersect(isa_otus, deseq_otus)

# Prune phyloseq to only these OTUs
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

# --- 6. Summarize relative abundance per sample ---
plot_data <- ps_df %>%
  group_by(Sample, Group_combined, Label) %>%
  summarize(sample_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(rel_abundance = sample_abundance / sum(sample_abundance)) %>%
  ungroup() %>%
  group_by(Group_combined, Label) %>%
  summarize(mean_rel_abundance = mean(rel_abundance), .groups = "drop")

# --- 7. Set x-axis order ---
group_levels <- c(
  "male_Control", "female_Control",
  "male_Mild", "female_Mild",
  "male_Severe", "female_Severe",
  "male_Deceased", "female_Deceased"
)
plot_data$Group_combined <- factor(plot_data$Group_combined, levels = group_levels)

# --- 8. Color palette ---
unique_genera <- sort(unique(plot_data$Label))
palette <- grDevices::colorRampPalette(c("lightskyblue", "springgreen4", "yellow", "magenta"))(length(unique_genera))

# --- 9. Plot 100% stacked bar chart ---
stacked_plot <- ggplot(plot_data, aes(
  x = Group_combined,
  y = mean_rel_abundance,
  fill = Label
)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  labs(
    title = "Genus-Level ASVs Significant in Both ISA & DESeq2",
    x = "Sex–Severity Group",
    y = "Mean Relative Abundance",
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

print(stacked_plot)
ggsave("stacked_genus_relative_abundance_noNA.png", stacked_plot, width = 12, height = 9, dpi = 300)

