library(dplyr)
library(ggplot2)
library(phyloseq)
library(pheatmap)
library(tidyr)
library(tibble)
library(forcats)

# --- 1. Get OTUs significant in BOTH ISA and DESeq2 ---
isa_otus <- significant_indicators$OTU 
deseq_otus <- deseq_all$OTU 
both_otus <- intersect(isa_otus, deseq_otus)
ps_both <- prune_taxa(both_otus, ps)

# --- 2. Clean metadata and build Group_combined ---
meta_df <- data.frame(sample_data(ps_both))
meta_df$Group_combined <- paste(meta_df$sex, meta_df$Severity, sep = "_")
meta_df <- meta_df %>% filter(!is.na(Group_combined))

ps_both <- prune_samples(rownames(meta_df), ps_both)
sample_data(ps_both) <- sample_data(meta_df)

# --- 3. Taxonomy table ---
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

# --- 4. Melt phyloseq to long dataframe ---
ps_df <- psmelt(ps_both)
ps_df <- left_join(ps_df, tax_df %>% select(OTU, Label), by = "OTU")

# --- 5. Summarize prevalence for plotting ---
plot_data <- ps_df %>%
  mutate(present = Abundance > 0) %>%
  group_by(Group_combined, Label, Sample) %>%
  summarize(sample_present = any(present), .groups = "drop") %>%
  group_by(Group_combined, Label) %>%
  summarize(count = sum(sample_present), .groups = "drop")


# --- 5. Keep top 50 OTUs by total prevalence ---
top_otus <- plot_data %>%
  group_by(Label) %>%
  summarize(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count)) %>%
  slice_head(n = 50) %>%
  pull(Label)

plot_data_sub <- plot_data %>% filter(Label %in% top_otus)

# --- 6. Order groups so male/female are side-by-side ---
group_levels <- c(
  "male_Control","female_Control",
  "male_Mild","female_Mild",
  "male_Severe","female_Severe",
  "male_Fatal","female_Fatal"
)
plot_data_sub$Group_combined <- factor(plot_data_sub$Group_combined, levels = group_levels)

# --- 7. Add Sex for bubble-plot coloring ---
plot_data_sub$Sex <- ifelse(grepl("^male", plot_data_sub$Group_combined), "Male", "Female")

# --- 8. Bubble plot ---
bubble_plot <- ggplot(plot_data_sub, aes(
  x = Group_combined,
  y = Label,
  size = count,
  color = Sex
)) +
  geom_point() +
  scale_color_manual(values = c("Male" = "deepskyblue", "Female" = "plum")) +
  theme_minimal() +
  labs(
    title = "Top 50 OTU Prevalence (Male vs Female)",
    x = "Group",
    y = "OTU",
    color = "Sex"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(bubble_plot)
ggsave("bubble_plot_top50_sexcolor.png", bubble_plot, width = 12, height = 10, dpi = 300)

# --- 9. Stacked bar plot ---
stacked_plot <- ggplot(plot_data_sub, aes(
  x = Group_combined,
  y = count,
  fill = Label
)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(
    title = "Top 50 OTUs Significant in both ISA & DESeq2",
    x = "Sex-Severity Group",
    y = "Number of Samples with Taxon Present",
    fill = "Taxon"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(stacked_plot)
ggsave("stacked_top50.png", stacked_plot, width = 12, height = 10, dpi = 300)

