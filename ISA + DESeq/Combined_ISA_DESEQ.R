# --- 1. Load Library ---
library(dplyr)
library(ggplot2)
library(phyloseq)
library(pheatmap)
library(tidyr)
library(tibble)

# --- 2. Load ISA + DESeq Data ---
isa_otus <- significant_indicators$OTU 
deseq_otus <- deseq_all$OTU 
both_otus <- intersect(isa_otus, deseq_otus) 
ps_both <- prune_taxa(both_otus, ps)

# --- 3. Clean metadata and remove samples with NA groups ---
meta_df <- data.frame(sample_data(ps_both))
meta_df$Group_combined <- paste(meta_df$sex, meta_df$Severity, sep = "_")
meta_df <- meta_df %>% filter(!is.na(Group_combined))
ps_both <- prune_samples(rownames(meta_df), ps_both)
sample_data(ps_both) <- sample_data(meta_df)

# --- 4. Taxonomy table ---
tax_df <- as.data.frame(tax_table(ps_both)) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Species = ifelse(Species %in% c("s__", "", NA), NA, Species),
    Genus   = ifelse(Genus   %in% c("g__", "", NA), NA, Genus),
    Family  = ifelse(Family  %in% c("f__", "", NA), NA, Family),
    Order   = ifelse(Order   %in% c("o__", "", NA), NA, Order),
    Class   = ifelse(Class   %in% c("c__", "", NA), NA, Class),
    Phylum  = ifelse(Phylum  %in% c("p__", "", NA), NA, Phylum),
    Domain  = ifelse(Domain  %in% c("d__", "", NA), NA, Domain)
  ) %>%
  mutate(
    Label = case_when(
      !is.na(Species) ~ paste0("s_", Species),
      !is.na(Genus)   ~ paste0("g_", Genus),
      !is.na(Family)  ~ paste0("f_", Family),
      !is.na(Order)   ~ paste0("o_", Order),
      !is.na(Class)   ~ paste0("c_", Class),
      !is.na(Phylum)  ~ paste0("p_", Phylum),
      !is.na(Domain)  ~ paste0("d_", Domain),
      TRUE            ~ "Unassigned"
    )
  )

# --- 5. Melt phyloseq to long dataframe ---
ps_df <- psmelt(ps_both)
ps_df <- left_join(ps_df, tax_df %>% select(OTU, Label), by = "OTU")

# --- 6. Summarize prevalence for plotting ---
plot_data <- ps_df %>%
  group_by(Group_combined, Label) %>%
  summarize(count = sum(Abundance > 0), .groups = "drop")

# --- 7. Keep top 50 OTUs by total prevalence ---
top_otus <- plot_data %>%
  group_by(Label) %>%
  summarize(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 50) %>%
  pull(Label)

plot_data_sub <- plot_data %>% filter(Label %in% top_otus)

# --- 8. Bubble plot ---
gg_heat <- ggplot(plot_data_sub, aes(x = Group_combined, y = Label, size = count)) +
  geom_point(color = "steelblue") +  # single color
  guides(color = "none") +           # remove legend
  theme_minimal() +
  labs(
    title = "Top 50 OTU Prevalence Across Groups",
    x = "Group",
    y = "OTU"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(gg_heat)
ggsave("bubble_plot_top50.png", plot = gg_heat, width = 12, height = 10, dpi = 300)

# --- 9. Stacked bar plot ---
ISA_DESeq2_combined <- ggplot(plot_data_sub, aes(x = Group_combined, y = count, fill = Label)) +
  geom_bar(stat = "identity") +
  labs(
    title = "Top 50 OTUs Significant in Both ISA and DESeq2",
    x = "Sex-Severity Group",
    y = "Number of Samples with Taxon Present",
    fill = "Taxon"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(ISA_DESeq2_combined)
ggsave("stacked_bar_top50.png", plot = ISA_DESeq2_combined, width = 12, height = 10, dpi = 300)

