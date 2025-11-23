library(dplyr)
library(ggplot2)
library(phyloseq)

# --- 1. Identify taxa significant in both ISA and DESeq2 ---
isa_otus   <- significant_indicators$OTU
deseq_otus <- sig_deseq_annot$OTU
both_otus  <- intersect(isa_otus, deseq_otus)

ps_both <- prune_taxa(both_otus, ps)

meta_df <- data.frame(sample_data(ps_both))
meta_df$Group_combined <- paste(meta_df$sex, meta_df$Severity, sep="_")
sample_data(ps_both) <- sample_data(meta_df)

# --- 2. Taxonomy table ---
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


ps_df <- psmelt(ps_both)

ps_df <- left_join(ps_df, tax_df %>% select(OTU, Label), by="OTU")

plot_data <- ps_df %>%
  group_by(Group_combined, Label) %>%
  summarize(count = sum(Abundance > 0), .groups = "drop") 


# --- 3. Stacked bar plot ---
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

print (ISA_DESeq2_combined)
ggsave("ISA_DESeq2_combined.png", plot = ISA_DESeq2_combined, width = 15, height = 15, dpi = 300)
