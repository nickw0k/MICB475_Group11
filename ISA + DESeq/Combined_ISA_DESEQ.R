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
ps_both <- prune_taxa(both_otus, ps)

# --- 3. Clean metadata and build Group_combined ---
meta_df <- data.frame(sample_data(ps_both))
meta_df$Group_combined <- paste(meta_df$sex, meta_df$Severity, sep = "_")
meta_df <- meta_df %>% filter(!is.na(Group_combined))

ps_both <- prune_samples(rownames(meta_df), ps_both)
sample_data(ps_both) <- sample_data(meta_df)

# --- 4. Extract taxonomy table (family-level only) ---
tax_df <- as.data.frame(tax_table(ps_both)) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Family = ifelse(Family %in% c("f__", "", NA), NA, Family),
    Label = Family
  ) %>%
  filter(!is.na(Label)) %>%
  mutate(Label = gsub("^f_", "", Label))

# --- 5. Melt phyloseq and join family labels ---
ps_df <- psmelt(ps_both)
ps_df <- left_join(ps_df, tax_df %>% select(OTU, Label), by = "OTU") %>%
  filter(!is.na(Label))

# --- 6. Summarize relative abundance per group ---
plot_data <- ps_df %>%
  group_by(Group_combined, Label, Sample) %>%
  summarize(sample_abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Sample) %>%
  mutate(rel_abundance = sample_abundance / sum(sample_abundance)) %>%  # relative abundance per sample
  ungroup() %>%
  group_by(Group_combined, Label) %>%
  summarize(mean_rel_abundance = mean(rel_abundance), .groups = "drop")  # average across samples

# --- 7. Order x-axis ---
group_levels <- c(
  "male_Control", "female_Control",
  "male_Mild", "female_Mild",
  "male_Severe", "female_Severe",
  "male_Fatal", "female_Fatal"
)
plot_data$Group_combined <- factor(plot_data$Group_combined, levels = group_levels)

# --- 8. Color palette ---
unique_families <- sort(unique(plot_data$Label))
palette <- grDevices::colorRampPalette(
  c( "#4CA", "#1B5", "lightblue", "#29B6C2", "#006070" )
)(length(unique_families))

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
    title = "Family-Level ASVs Significant in Both ISA & DESeq2",
    x = "Sex–Severity Group",
    y = "Mean Relative Abundance",
    fill = "Family"
  ) +
  theme(
    title = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 40, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(face = "bold", size = 11),
    axis.title.y = element_text(face = "bold", size = 11 ),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.box = "vertical"
  )

print(stacked_plot)
ggsave("stacked_family_relative_abundance.png", stacked_plot, width = 12, height = 9, dpi = 300)
