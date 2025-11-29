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

# --- 4. Extract taxonomy table (genus-level only) ---
tax_df <- as.data.frame(tax_table(ps_both)) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Genus = ifelse(Genus %in% c("g__", "", NA), NA, Genus),
    Label = Genus
  ) %>%
  filter(!is.na(Label)) %>%
  mutate(Label = gsub("^g_", "", Label))

# --- 5. Melt phyloseq and join genus labels ---
ps_df <- psmelt(ps_both)
ps_df <- left_join(ps_df, tax_df %>% select(OTU, Label), by = "OTU") %>%
  filter(!is.na(Label))

# --- 6. Summarize genus prevalence for plotting ---
plot_data <- ps_df %>%
  mutate(present = Abundance > 0) %>%
  group_by(Group_combined, Label, Sample) %>%
  summarize(sample_present = any(present), .groups = "drop") %>%
  group_by(Group_combined, Label) %>%
  summarize(count = sum(sample_present), .groups = "drop")

# --- 7. Order x-axis: M+F Control → Mild → Severe → Fatal ---
group_levels <- c(
  "male_Control", "female_Control",
  "male_Mild", "female_Mild",
  "male_Severe", "female_Severe",
  "male_Fatal", "female_Fatal"
)

plot_data$Group_combined <- factor(plot_data$Group_combined, levels = group_levels)

# --- 8. Generate soft green/blue color palette ---
unique_genera <- sort(unique(plot_data$Label))
palette <- grDevices::colorRampPalette(
  c("#1B5", "#4CA", "#29B6F2", "#006049")
)(length(unique_genera))

# --- 9. Plot stacked bar chart (genus-level only) ---
stacked_plot <- ggplot(plot_data, aes(
  x = Group_combined,
  y = count,
  fill = Label
)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette) +
  theme_minimal() +
  labs(
    title = "Genus-Level ASVs Significant in Both ISA & DESeq2",
    x = "Sex–Severity Group",
    y = "Number of Samples with Taxon Present",
    fill = "Genus"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.box = "vertical"
  )
p <- p + theme(
  axis.title.x = element_text(face = "bold"),
  axis.title.y = element_text(face = "bold")
)


print(stacked_plot)
ggsave("stacked_genus_combined.png", stacked_plot, width = 12, height = 9, dpi = 300)
