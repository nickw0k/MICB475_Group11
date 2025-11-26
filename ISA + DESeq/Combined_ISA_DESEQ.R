library(dplyr)
library(ggplot2)
library(phyloseq)
library(pheatmap)

# --- Load data from ISA and DESeq Analysis ---
isa_otus   <- significant_indicators$OTU
deseq_otus <- deseq_all$OTU
both_otus  <- intersect(isa_otus, deseq_otus)

# --- Prune phyloseq object to OTUs in both ISA and DESeq2 ---
ps_both <- prune_taxa(both_otus, ps)

# --- Clean metadata and create Group_combined ---
meta_df <- data.frame(sample_data(ps_both))
meta_df$Group_combined <- paste(meta_df$sex, meta_df$Severity, sep = "_")

# Remove samples with NA in Group_combined
meta_df <- meta_df %>% filter(!is.na(Group_combined))
ps_both <- prune_samples(rownames(meta_df), ps_both)
sample_data(ps_both) <- sample_data(meta_df)

# --- Taxonomy table with fallback labels ---
tax_df <- as.data.frame(tax_table(ps_both)) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Species = na_if(Species, "s__"),
    Genus   = na_if(Genus, "g__"),
    Family  = na_if(Family, "f__"),
    Order   = na_if(Order, "o__"),
    Class   = na_if(Class, "c__"),
    Phylum  = na_if(Phylum, "p__"),
    Domain  = na_if(Domain, "d__")
  ) %>%
  mutate(
    Label = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain),
    Label = ifelse(is.na(Label), "Unclassified", Label)
  )

# --- Melt phyloseq to long dataframe ---
ps_df <- psmelt(ps_both)
ps_df <- left_join(ps_df, tax_df %>% select(OTU, Label), by = "OTU")

# --- Summarize for plotting ---
plot_data <- ps_df %>%
  group_by(Group_combined, Label) %>%
  summarize(count = sum(Abundance > 0), .groups = "drop")

# --- Stacked bar plot ---
ISA_DESeq2_combined <- ggplot(plot_data, aes(x = Group_combined, y = count, fill = Label)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Taxa Significant in Both ISA and DESeq2",
    x = "Sex-Severity Group",
    y = "Number of Samples with Taxon Present",
    fill = "Taxon"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(ISA_DESeq2_combined)
ggsave("ISA_DESeq2_combined1.png", plot = ISA_DESeq2_combined, width = 15, height = 15, dpi = 300)

# --- Heatmap using pheatmap ---
heat_df <- plot_data %>%
  pivot_wider(names_from = Group_combined, values_from = count, values_fill = 0) %>%
  column_to_rownames("Label")

pheatmap(heat_df, 
         color = colorRampPalette(c("white", "blue"))(50),
         cluster_rows = TRUE, cluster_cols = FALSE,
         main = "OTU prevalence across groups",
         filename = "heatmap.png",
         width = 2000, height = 1500, res = 150)

# --- Bubble plot ---
gg_heat <- ggplot(plot_data, aes(x = Group_combined, y = Label, size = count)) +
  geom_point(color = "steelblue") +  # single color, legend removed
  guides(color = "none") +
  theme_minimal() +
  labs(
    title = "OTU Prevalence Across Groups",
    x = "Group",
    y = "OTU"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(gg_heat)
ggsave("bubble_plot.png", plot = gg_heat, width = 12, height = 15, dpi = 300)
