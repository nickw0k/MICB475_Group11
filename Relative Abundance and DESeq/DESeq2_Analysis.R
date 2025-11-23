# --- 1. Load libraries ---
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(DESeq2)
library(ape)

# --- 2. Load files ---
metaFP <- "Mexico Dataset QIIME2 files/covid_mex_metadata.tsv"
meta <- read.delim(metaFP, sep="\t")

otuFP <- "Mexico_Dataset_QIIME2_files/mexico-feature-table.txt"
otu <- read.delim(otuFP, sep="\t", skip=1)

taxFP <- "Mexico_Dataset_QIIME2_files/mexico-taxonomy.tsv"
tax <- read.delim(taxFP, sep="\t")

phyFP <- "Mexico_Dataset_QIIME2_files/mexico-tree.nwk"
phy <- read.tree(phyFP)

# --- 3. Build phyloseq object ---
otu_mat <- as.matrix(otu[,-1])
rownames(otu_mat) <- otu$X.OTU.ID
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

meta_df <- meta[,-1]
rownames(meta_df) <- meta$sample.id
META <- sample_data(meta_df)

tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into=c("Domain","Phylum","Class","Order",
                         "Family","Genus","Species"),
           sep=";", fill="right") %>%
  as.matrix()
tax_mat <- tax_mat[,-1]
rownames(tax_mat) <- tax$Feature.ID
TAX <- tax_table(tax_mat)

ps <- phyloseq(OTU, META, TAX, phy)

# --- 4. Create severity + Group_combined ---

meta2 <- as.data.frame(sample_data(ps))

meta2$sex[meta2$sex %in% c("", " ", "NA")] <- NA


meta2$Severity <- NA
meta2$Severity[grepl("asymptomatic",        meta2$Group, ignore.case=TRUE)] <- "Control"
meta2$Severity[grepl("ambulatory positive", meta2$Group, ignore.case=TRUE)] <- "Mild"
meta2$Severity[grepl("hospitalized positive", meta2$Group, ignore.case=TRUE)] <- "Severe"
meta2$Severity[grepl("deceased hospitalized", meta2$Group, ignore.case=TRUE)] <- "Fatal"

keep <- !is.na(meta2$sex) & !is.na(meta2$Severity)
names(keep) <- rownames(meta2)
ps <- prune_samples(keep, ps)

meta2 <- as.data.frame(sample_data(ps))

meta2$Group_combined <- paste(meta2$sex, meta2$Severity, sep="_")
meta2$Group_combined <- factor(meta2$Group_combined)

sample_data(ps) <- sample_data(meta2)


# --- 5. Prepare taxonomy names ---
tax_df <- as.data.frame(tax_table(ps)) %>%
  rownames_to_column("OTU") %>%
  mutate(
    Species = na_if(Species, "s__"),
    Genus   = na_if(Genus,   "g__"),
    Family  = na_if(Family,  "f__"),
    Order   = na_if(Order,   "o__"),
    Class   = na_if(Class,   "c__"),
    Phylum  = na_if(Phylum, "p__"),
    Domain  = na_if(Domain, "d__")
  ) %>%
  mutate(
    Label = coalesce(Species, Genus, Family, Order, Class, Phylum, Domain),
    Label = ifelse(is.na(Label), "Unclassified", Label)
  )

# --- 6. DESeq2 Analysis + Pairwise Comparison ---
dds <- phyloseq_to_deseq2(ps, ~ Group_combined)
dds <- estimateSizeFactors(dds, type="poscounts")
dds <- DESeq(dds, fitType="parametric")

groups <- levels(dds$Group_combined)
comparisons <- combn(groups, 2, simplify=FALSE)

deseq_list <- list()
for (comp in comparisons) {
  res <- results(dds, contrast=c("Group_combined", comp[1], comp[2]))
  res <- as.data.frame(res) %>%
    rownames_to_column("OTU") %>%
    mutate(contrast = paste(comp[1], "vs", comp[2], sep="_"))
  deseq_list[[paste(comp[1], comp[2], sep="_")]] <- res
}

deseq_all <- bind_rows(deseq_list)

# --- 7. ALL significant taxa (padj ≤ 0.05) ---
sig_deseq <- deseq_all %>%
  filter(!is.na(padj) & padj <= 0.05)

sig_deseq_annot <- sig_deseq %>%
  left_join(tax_df %>% select(OTU, Label), by="OTU")

# --- 8. Plot all significant taxa ---
DESeq2_plot <- ggplot(sig_deseq_annot,
       aes(x=reorder(Label, abs(log2FoldChange)),
           y=log2FoldChange,
           fill=contrast)) +
  geom_col() +
  coord_flip() +
  theme_minimal() +
  labs(
    title="All Significant Taxa (DESeq2, padj ≤ 0.05)",
    x="Taxon",
    y="log2 Fold Change",
    fill="Comparison"
  )

ggsave("DESeq2.png", plot = DESeq2_plot, width = 15, height = 15, dpi = 300)

