library(phyloseq)
library(ape)
library(tidyverse)
library(vegan)
library(dplyr)
library(picante)
library(ggpubr)
library(rstatix)

#### Load Data ####
covid_mex_metadata_fp <- "COVID_Mexico_Final_Updated/covid_mex_metadata_edited.tsv"
covid_mex_metadata <- read_delim(covid_mex_metadata_fp, delim="\t")

covid_mex_otu_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/table_export/feature-table.txt"
covid_mex_otu <- read_delim(covid_mex_otu_fp, delim="\t", skip=1)

covid_mex_tax_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/taxonomy_export/taxonomy.tsv"
covid_mex_taxonomy <- read_delim(covid_mex_tax_fp, delim="\t")

covid_mex_phylotree_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/rooted_tree_export/tree.nwk"
covid_mex_phylotree <- read.tree(covid_mex_phylotree_fp)

#### Format OTU table ####
# OTU tables should be a matrix with rownames and colnames as OTUs and SampleIDs respectively
# save everything except first column (OTU ID) into a matrix
covid_mex_otu_matrix <- as.matrix(covid_mex_otu[, -1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(covid_mex_otu_matrix) <- covid_mex_otu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
OTU <- otu_table(covid_mex_otu_matrix, taxa_are_rows = TRUE)
class(OTU)

#### Format sample metadata ####
# Save everything except sample-id as new data frame
covid_mex_metadata_df <- as.data.frame(covid_mex_metadata[, -1])
# Make sample-id the rownames
rownames(covid_mex_metadata_df) <- covid_mex_metadata$'sample-id'
# Changing "Group" column name to "Severity" and standardizing "Group" column into 4 categories
# which are: Asymptomatic, Hospitalized Negative, Hospitalized Positive, and Deceased
covid_mex_metadata_edited <- covid_mex_metadata_df %>%
  mutate(Group = case_when(
    grepl("asymptomatic", Group, ignore.case = TRUE) ~ "Asymptomatic",
    grepl("ambulatory negative", Group, ignore.case = TRUE) ~ "Mild negative",
    grepl("ambulatory positive", Group, ignore.case = TRUE) ~ "Mild positive",
    grepl("hospitalized positive", Group, ignore.case = TRUE) &
      !grepl("deceased", Group, ignore.case = TRUE) ~ "Severe positive",
    grepl("deceased", Group, ignore.case = TRUE) ~ "Deceased",
    TRUE ~ NA_character_     # <-- make unmatched truly missing
  )) %>%
  rename(severity = Group)  # rename "Group" to "severity"
# Make phyloseq sample data with sample_data() function
SAMP <- sample_data(covid_mex_metadata_edited)
class(SAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
covid_mex_taxonomy_matrix <- covid_mex_taxonomy %>% select(-Confidence)%>%
  separate(col=Taxon, sep=";"
           , into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs
covid_mex_taxonomy_matrix <- covid_mex_taxonomy_matrix[, -1]
# Make sample-ids the rownames
rownames(covid_mex_taxonomy_matrix) <- covid_mex_taxonomy$`Feature ID`
# Make taxa table
TAX <- tax_table(covid_mex_taxonomy_matrix)
class(TAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
covid_mexico_phyloseq <- phyloseq(OTU, SAMP, TAX, covid_mex_phylotree)

#### Looking at phyloseq object ####
# View components of phyloseq object with the following commands
otu_table(covid_mexico_phyloseq)
sample_data(covid_mexico_phyloseq)
tax_table(covid_mexico_phyloseq)
phy_tree(covid_mexico_phyloseq)

# View sample names in our data set
sample_names(covid_mexico_phyloseq)

# Number of samples in our data set
nsamples(covid_mexico_phyloseq)

# Sum of reads in each sample
sample_sums(covid_mexico_phyloseq)

# 3 samples with the most reads from our data set
samps_most_reads <- names(sort(sample_sums(covid_mexico_phyloseq), decreasing = TRUE)[1:3])
get_taxa(covid_mexico_phyloseq, samps_most_reads)

# Names of taxa from our data set
taxa_names(covid_mexico_phyloseq)

# Number of taxa we have from our data set
ntaxa(covid_mexico_phyloseq)

# 3 most abundant taxa from our data set
taxa_most_abundant <- names(sort(taxa_sums(covid_mexico_phyloseq), decreasing = TRUE)[1:3])
get_sample(covid_mexico_phyloseq, taxa_most_abundant)

# Remove any non-bacterial sequences, if any
covid_mexico_phyloseq_filt <- subset_taxa(
  covid_mexico_phyloseq,
  !(Family %in% c("Mitochondria") |
      Order %in% c("Chloroplast")))

# Remove ASVs that have less than 5 counts total
covid_mexico_phyloseq_filt_nolow <- filter_taxa(covid_mexico_phyloseq_filt,
                                                function(x) sum(x)>5,
                                                prune = TRUE)

# Remove samples with less than 100 reads
covid_mexico_phyloseq_filt_nolow_samps <- prune_samples(sample_sums(covid_mexico_phyloseq_filt_nolow)>100,
                                                        covid_mexico_phyloseq_filt_nolow)

# Remove samples where sex is NA, severity is NA, and filter out "Mild negative" group in severity column
covid_mexico_phyloseq_final <- subset_samples(
  covid_mexico_phyloseq_filt_nolow_samps,
  !is.na(sex) & !is.na(severity) & severity != "Mild negative")

# Rarefy samples
# rngseed sets a random number. To be able to reproduce this exact analysis each time, need to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(covid_mexico_phyloseq_final))), cex=0.1)
covid_mexico_phyloseq_rare <- rarefy_even_depth(covid_mexico_phyloseq_final, rngseed = 11, sample.size = 48792)

#### Saving ####
save(covid_mexico_phyloseq_final, file="covid_mexico_phyloseq_final.RData")
save(covid_mexico_phyloseq_rare, file="covid_mexico_phyloseq_rare.RData")

#### Alpha diversity ###

## 1) Build the base phyloseq richness plot (your code)
base_alpha <- plot_richness(covid_mexico_phyloseq_rare, x = "sex", measures = "Shannon") +
  xlab("Sex") +
  geom_boxplot(outlier.shape = NA) +
  facet_wrap(~severity) +
  theme_minimal()

base_alpha$data <- subset(base_alpha$data, !is.na(severity) & severity != "NA" & !is.na(sex))

## 2) Compute Kruskal–Wallis p-values within each severity (Shannon ~ sex)
alpha_df <- estimate_richness(covid_mexico_phyloseq_rare, measures = "Shannon") %>%
  tibble::rownames_to_column("SampleID")

meta_df <- as.data.frame(sample_data(covid_mexico_phyloseq_rare)) %>%
  tibble::rownames_to_column("SampleID")

df <- left_join(alpha_df, meta_df, by = "SampleID") %>%
  mutate(
    severity = as.factor(severity),
    sex = as.factor(sex)
  ) %>%
  filter(!is.na(severity), severity != "NA", !is.na(sex))

kw_df <- base_alpha$data %>%
  dplyr::filter(!is.na(severity), severity != "NA", !is.na(sex)) %>%
  dplyr::group_by(severity) %>%
  dplyr::summarise(
    p = if (dplyr::n_distinct(sex) < 2) NA_real_
    else stats::kruskal.test(value ~ sex, data = dplyr::cur_data())$p.value,
    y = max(value, na.rm = TRUE) * 1.06,   # vertical position inside each facet
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    label = ifelse(is.na(p), "KW p = NA", paste0("KW p = ", signif(p, 3))),
    sex = "male"  # anchor text at the 'male' x position in every facet
  )

## 3) Add facet-wise labels + significance brackets

# Compute max y per severity for bracket placement
gg_richness_severity <- base_alpha +
  stat_compare_means(
    method = "wilcox.test",
    comparisons = list(c("female", "male")),
    label = "p.signif",
    hide.ns = FALSE,
    bracket.size = 0.6,
    tip.length = 0.05,
    vjust = 1
  ) +
  geom_text(
    data = kw_df,
    aes(x = sex, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.15,  # push to the right of 'male'
    vjust = 0,
    size = 3
  ) +
  coord_cartesian(clip = "off")

gg_richness_severity_final <- gg_richness_severity +
  theme_classic(base_size = 11) +
  theme(
    # Panel appearance
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.spacing = unit(1.2, "lines"),
    
    # Facet strip (severity labels)
    strip.background = element_blank(),
    strip.text = element_text(
      size = 11,
      face = "bold",
      margin = margin(b = 6)
    ),
    
    # Axis text & titles
    axis.title.x = element_text(size = 11, face = "bold", margin = margin(t = 8)),
    axis.title.y = element_text(size = 11, face = "bold", margin = margin(r = 8)),
    axis.text = element_text(size = 10),
    
    # Remove legend (not needed)
    legend.position = "none",
    
    # Clean look
    axis.line = element_line(color = "black"),
    plot.margin = margin(10, 12, 10, 10)
  )

gg_richness_severity_final

ggsave(filename = "plot_richness_severity.png", gg_richness_severity_final, height=4, width=6)

estimate_richness(covid_mexico_phyloseq_rare)

# Phylogenetic diversity
# Calculate Faith's phylogenetic diversity as PD
phylo_dist <- pd(
  t(otu_table(covid_mexico_phyloseq_rare)),
  phy_tree(covid_mexico_phyloseq_rare),
  include.root = FALSE
)

# Add PD to metadata table
sample_data(covid_mexico_phyloseq_rare)$PD <- phylo_dist$PD

# Convert sample_data to a data.frame for ggplot + stats
# Force sample_data to a plain data.frame first (avoids <sample_data> pipe issues)
df_pd <- data.frame(sample_data(covid_mexico_phyloseq_rare), check.names = FALSE) %>%
  tibble::rownames_to_column("SampleID") %>%
  dplyr::mutate(
    severity = as.factor(severity),
    sex = as.factor(sex)
  ) %>%
  dplyr::filter(!is.na(severity), severity != "NA", !is.na(sex), !is.na(PD))

# 1) Kruskal–Wallis per severity (PD ~ sex)
kw_brackets <- df_pd %>%
  group_by(severity) %>%
  kruskal_test(PD ~ sex) %>%
  add_significance("p") %>%
  # Define the groups for the brackets
  mutate(
    group1 = "female",
    group2 = "male",
    # Create the label for the text: "KW p = 0.XXX"
    kw_label = paste0("KW p = ", format.pval(p, digits = 3, eps = 0.001))
  )

# 2) Calculate position for brackets and text per facet
pos_df <- df_pd %>%
  group_by(severity) %>%
  summarise(
    # INCREASED spacing from max PD value:
    # Set y.position for the bracket/ns label
    y_ns_signif = max(PD, na.rm = TRUE) + 0.15 * IQR(PD, na.rm = TRUE),
    # Set y.position for the KW p-value text (slightly higher to prevent overlap)
    y_kw_text   = max(PD, na.rm = TRUE) + 0.22 * IQR(PD, na.rm = TRUE),
    .groups = "drop"
  )

# Join the positions to the brackets data frame
kw_brackets <- kw_brackets %>%
  left_join(pos_df, by = "severity")

# 3) Create a data frame for the Kruskal-Wallis p-value text label
# This will be positioned specifically for the KW p-value.
kw_text_df <- kw_brackets %>%
  mutate(
    # Use the separate higher position
    y.position = y_kw_text,
    # Pin the text position to the 'male' boxplot column (position x=2)
    xmin = "male",
    xmax = "male",
    # Use the combined label for the text
    label = kw_label
  )

# Generate the final plot
plot_pd_severity_final <- ggplot(df_pd, aes(x = sex, y = PD)) +
  geom_jitter(
    width = 0.15,          # Maximum horizontal jitter. Keep it narrow.
    shape = 21,            # Circle shape (21) is often preferred for visibility
    size = 2,              # Point size
    fill = "gray70",       # Fill color for the points
    color = "black",       # Border color for the points
    alpha = 0.7            # Slightly transparent
  ) +
  geom_boxplot(width = 0.45, outlier.shape = NA, fill = NA, color = "black") +
  facet_wrap(~severity, scales = "free_y") +
  
  # Add significance brackets (using p.signif like ns, *, **)
  stat_pvalue_manual(
    kw_brackets,
    label = "p.signif",
    xmin = "group1",
    xmax = "group2",
    y.position = "y_ns_signif", # Use the lower position for the bracket/ns
    tip.length = 0.02,
    bracket.size = 0.5
  ) +
  
  # Add Kruskal-Wallis p-value text (KW p = 0.XXX)
  # Using stat_pvalue_manual for precise placement
  stat_pvalue_manual(
    kw_text_df,
    label = "label",  # The column with the "KW p = ..." string
    xmin = "xmin",    # Pinned to 'male'
    xmax = "xmax",    # Pinned to 'male'
    y.position = "y.position", # Use the higher position for KW text
    tip.length = 0,
    bracket.size = 0,
    size = 3,
    hjust = 1.05 # Move the text slightly to the right of the 'male' position
  ) +
  
  # Tweak scales for better appearance and spacing
  # Increased the top expansion value to ensure enough room for both bracket and text
  scale_x_discrete(expand = expansion(mult = c(0.6, 0.6))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) + # <-- Increased top room to 0.3
  coord_cartesian(clip = "off") + # Changed to "off" to prevent clipping if needed
  
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    panel.spacing = unit(0.4, "lines"),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold", margin = margin(b = 6)),
    legend.position = "none"
  )

plot_pd_severity_final

ggsave("plot_pd_severity.png", plot_pd_severity_final, height = 4, width = 6)

#### Beta diversity #####
dm_unifrac <- UniFrac(covid_mexico_phyloseq_rare, weighted = TRUE)

samp_dat_wdiv <- data.frame(sample_data(covid_mexico_phyloseq_rare),
                            estimate_richness(covid_mexico_phyloseq_rare))

permanova <- adonis2(
  dm_unifrac ~ severity * sex,
  data = samp_dat_wdiv,
  permutations = 999,
  by = "terms")

permanova

## Extract PERMANOVA (severity) results safely
R2 <- unname(permanova["severity", "R2"])
p  <- unname(permanova["severity", "Pr(>F)"])

label <- sprintf(
  "PERMANOVA (Weighted UniFrac)\nSeverity: R² = %.3f, p = %.3g",
  R2, p
)

pcoa_wu <- ordinate(covid_mexico_phyloseq_rare, method="PCoA", distance=dm_unifrac)

base_beta <- plot_ordination(covid_mexico_phyloseq_rare, pcoa_wu, color="severity", shape="sex") +
  labs(shape = "Sex", colour = "Severity") +
  theme_minimal()

gg_pcoa <- base_beta +
  stat_ellipse(
    data = base_beta$data,
    aes(x = Axis.1, y = Axis.2, group = severity, color = severity),
    type = "norm",
    level = 0.68,
    linewidth = 0.7,
    alpha = 0.6,
    show.legend = FALSE,
    inherit.aes = FALSE
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = label,
    hjust = 1.05,
    vjust = 1.15,
    size = 3.5
  ) +
  coord_cartesian(clip = "off")

gg_pcoa_final <- gg_pcoa +
  # Make ellipses a subtle fill instead of thin outlines (cleaner)
  stat_ellipse(
    aes(group = severity, fill = severity),
    geom = "polygon",
    type = "norm",
    level = 0.68,
    inherit.aes = TRUE,
    linewidth = 0.0,
    alpha = 0.12,
    show.legend = FALSE
  ) +
  # Keep your outline ellipses (optional; comment out if you prefer fill only)
  stat_ellipse(
    aes(group = severity, color = severity),
    type = "norm",
    level = 0.68,
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  # Make points slightly larger + readable
  geom_point(size = 2.4, alpha = 0.9) +
  # Put PERMANOVA label in a neat corner "caption" style
  annotate(
    "label",
    x = Inf, y = Inf,
    label = label,   # <- your PERMANOVA label string
    hjust = 1.02, vjust = 1.15,
    size = 3.5,
    label.size = 0.25
  ) +
  # Clean axes/labels (keep your % variance labels if you already set them)
  labs(
    x = "PCoA1 [96.3%]",
    y = "PCoA2 [3.6%]",
    color = "Severity",
    fill  = "Severity",
    shape = "Sex"
  ) +
  # Make legends compact and non-distracting
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(order = 2, override.aes = list(size = 3, alpha = 1)),
    fill  = "none"
  ) +
  # Manuscript-style theme
  theme_classic(base_size = 11) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.title = element_text(face = "bold"),
    axis.text  = element_text(color = "black"),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.box = "vertical",
    legend.spacing.y = unit(3, "pt"),
    plot.margin = margin(10, 12, 10, 10),
    plot.title.position = "plot"
  ) +
  coord_cartesian(clip = "off")

gg_pcoa_final

ggsave("plot_pcoa.png", gg_pcoa_final, height=4, width=6)

#### Taxonomy bar plots ####
# Plot bar plot of taxonomy
plot_bar(covid_mexico_phyloseq_rare, fill = "Phylum")

# Convert to relative abundance
covid_mexico_phyloseq_RA <- transform_sample_counts(covid_mexico_phyloseq_rare, function(x) x/sum(x))

# To remove black bars, "glom" by phylum first
covid_mexico_phylum <- tax_glom(covid_mexico_phyloseq_RA, taxrank = "Phylum", NArm=FALSE)

plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~sex, scales = "free_x")

plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~severity, scales = "free_x")

gg_taxa_sex <- plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~sex, scales = "free_x")
gg_taxa_sex

gg_taxa_severity <- plot_bar(covid_mexico_phylum, fill="Phylum") + 
  facet_wrap(.~severity, scales = "free_x")
gg_taxa_severity

ggsave("plot_taxonomy_sex.png"
       , gg_taxa_sex
       , height=8, width =12)

ggsave("plot_taxonomy_severity.png"
       , gg_taxa_severity
       , height=8, width =12)




#### --- ISA/DESeq2 - Differential Abundance/Integrated Analysis ---####
#### --- ISA --- ####
# Load libraries 
#library(phyloseq)
#library(ape)
#library(tidyverse)
#library(vegan)
#library(dplyr)
library(picante)
library(ggpubr)
library(rstatix)
library(DESeq2)
library(indicspecies)
library(forcats)
library(tibble)
library(ggplot2)



# --- 1. Load phyloseq objects ---
load("covid_mexico_phyloseq_final.RData")  # unrarefied
load("covid_mexico_phyloseq_rare.RData")   # rarefied for ISA

# --- 2. Define sex-severity groups for both phyloseq objects ---
prepare_groups <- function(ps) {
  meta <- data.frame(sample_data(ps))
  if(!"Group_combined" %in% colnames(meta)) {
    meta$Severity <- NA_character_
    meta$Severity[grepl("asymptomatic", tolower(meta$severity))] <- "Control"
    meta$Severity[grepl("ambulatory positive", tolower(meta$severity))] <- "Mild"
    meta$Severity[grepl("hospitalized positive", tolower(meta$severity))] <- "Severe"
    meta$Severity[grepl("deceased", tolower(meta$severity))] <- "Deceased"
    meta$Group_combined <- paste(meta$sex, meta$Severity, sep="_")
    sample_data(ps) <- sample_data(meta)
  }
  return(ps)
}

covid_mexico_phyloseq_final <- prepare_groups(covid_mexico_phyloseq_final)
covid_mexico_phyloseq_rare <- prepare_groups(covid_mexico_phyloseq_rare)

# --- 3. Collapse to Genus level ---
ps_genus_rare <- tax_glom(covid_mexico_phyloseq_rare, taxrank="Genus", NArm=TRUE)
ps_genus_unrare <- tax_glom(covid_mexico_phyloseq_final, taxrank="Genus", NArm=TRUE)

# --- 4. Indicator Species Analysis (ISA) - rarefied ---
otu_mat <- as.data.frame(t(otu_table(ps_genus_rare)))
meta_df <- data.frame(sample_data(ps_genus_rare))
meta_clean <- meta_df[!is.na(meta_df$Group_combined), ]
otu_clean <- otu_mat[rownames(meta_clean), ]
otu_clean[is.na(otu_clean)] <- 0

isa_mpt <- multipatt(
  otu_clean,
  meta_clean$Group_combined,
  func = "IndVal.g",
  duleg = TRUE,
  control = how(nperm = 999)
)

sig_table <- data.frame(isa_mpt$sign) %>%
  rownames_to_column("OTU") %>%
  filter(p.value <= 0.05)

group_names <- colnames(isa_mpt$comb)
sig_table$cluster <- group_names[sig_table$index]

significant_indicators <- sig_table %>%
  select(OTU, stat, p.value, cluster)

# --- 5. Add Genus names ---
tax_df <- as.data.frame(tax_table(ps_genus_rare)) %>%
  rownames_to_column("OTU") %>%
  mutate(Label = Genus) %>%
  filter(!is.na(Label) & Label != "" & Label != "g__")

isa_plot_df <- significant_indicators %>%
  left_join(tax_df %>% select(OTU, Label), by="OTU") %>%
  arrange(cluster, desc(stat)) %>%
  group_by(cluster) %>% slice_max(stat, n=10) %>% ungroup()

# --- 6. Plot ---

isa_plot <- ggplot(isa_plot_df, aes(x=reorder(Label, stat), y=stat, fill=cluster)) +
  geom_col() +
  coord_flip() +
  labs(title="Top Indicator Genera by Sex-Severity Cluster (ISA)",
       x="Genus", y="Indicator Value", fill="Sex-Severity Group") +
  theme_minimal() +
  theme(axis.text=element_text(size=11),
        axis.title=element_text(size=12, face="bold"),
        legend.title=element_text(size=11),
        legend.text=element_text(size=10))

print(isa_plot)
ggsave("ISA_top_indicator_genera.png", isa_plot, width=10, height=6, dpi=300)


#### ---- Differential Abundance Analysis ---- ####

# --- 1. Load Libraries ---
# should be loaded in ISA code

# --- 2. Load Data ---
metadata_fp <- "COVID_Mexico_Final_Updated/covid_mex_metadata_edited.tsv"
metadata <- read_delim(covid_mex_metadata_fp, delim="\t")

otu_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/table_export/feature-table.txt"
otu <- read_delim(covid_mex_otu_fp, delim="\t", skip=1)

tax_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/taxonomy_export/taxonomy.tsv"
taxonomy <- read_delim(covid_mex_tax_fp, delim="\t")

phylotree_fp <- "COVID_Mexico_Final_Updated/covid_mexico_final_updated_export/rooted_tree_export/tree.nwk"
phylotree <- read.tree(covid_mex_phylotree_fp)

meta <- read.delim(metadata_fp, sep = "\t", stringsAsFactors = FALSE)
otu  <- read.delim(otu_fp, sep = "\t", skip = 1, stringsAsFactors = FALSE)
tax  <- read.delim(tax_fp, sep = "\t", stringsAsFactors = FALSE)
phy  <- read.tree(phylotree_fp)

# --- 3. Phyloseq Data Prep ---
otu_mat <- as.matrix(otu[, -1])
rownames(otu_mat) <- otu$X.OTU.ID
OTU <- otu_table(otu_mat, taxa_are_rows = TRUE)

tax_mat <- tax %>%
  select(-Confidence) %>%
  separate(Taxon, into = c("Domain","Phylum","Class","Order","Family","Genus","Species"),
           sep = ";", fill = "right", extra = "merge") %>%
  mutate(across(everything(), ~ str_replace_all(.x, "^[a-z]__", ""))) %>%
  as.data.frame(stringsAsFactors = FALSE)

tax_mat2 <- as.matrix(tax_mat[, -1])
rownames(tax_mat2) <- tax$Feature.ID
TAX <- tax_table(tax_mat2)

rownames(meta) <- meta$sample.id
META <- sample_data(meta)

ps <- phyloseq(OTU, TAX, META, phy)

# --- 4. Define sex-severity groups ---
sample_data(ps)$sex[sample_data(ps)$sex %in% c("", " ", "NA")] <- NA

sample_data(ps)$Severity <- NA
sample_data(ps)$Severity[grepl("asymptomatic", sample_data(ps)$Group, ignore.case = TRUE)] <- "Control"
sample_data(ps)$Severity[grepl("ambulatory positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Mild"
sample_data(ps)$Severity[grepl("hospitalized positive", sample_data(ps)$Group, ignore.case = TRUE)] <- "Severe"
sample_data(ps)$Severity[grepl("Deceased hospitalized", sample_data(ps)$Group, ignore.case = TRUE)] <- "Deceased"

sample_data(ps)$Group_combined <- paste(sample_data(ps)$sex, sample_data(ps)$Severity, sep = "_")

desired_levels <- c("male_Control","male_Mild","male_Severe","male_Deceased",
                    "female_Control","female_Mild","female_Severe","female_Deceased")
present_levels <- intersect(desired_levels, unique(sample_data(ps)$Group_combined))
sample_data(ps)$Group_combined <- factor(sample_data(ps)$Group_combined, levels = present_levels)

ps <- subset_samples(ps, !is.na(Group_combined))

# --- 5. Filter sparse taxa ---
min_prevalence <- 3
min_count      <- 5
ps_filt <- filter_taxa(ps, function(x) sum(x >= min_count) >= min_prevalence & sum(x) > 10, TRUE)

# --- 6. Collapse to genus level ---
ps_genus <- tax_glom(ps_filt, taxrank = "Genus", NArm = TRUE)

tax_df_genus <- as.data.frame(tax_table(ps_genus)) %>%
  rownames_to_column("OTU") %>%
  mutate(Label = ifelse(!is.na(Genus) & Genus != "", Genus, NA)) %>%
  filter(!is.na(Label)) %>%
  select(OTU, Label)

# Keep only OTUs with valid genus
ps_genus <- prune_taxa(tax_df_genus$OTU, ps_genus)

tax_table(ps_genus) <- tax_table(as.matrix(tax_df_genus %>% column_to_rownames("OTU")))

# --- 7. DESeq2 analysis ---
dds <- phyloseq_to_deseq2(ps_genus, ~ Group_combined)
dds <- estimateSizeFactors(dds, type = "poscounts")
dds <- DESeq(dds, fitType = "parametric")

# --- 8. Pairwise contrasts + LFC shrinkage ---
groups <- levels(droplevels(sample_data(ps_genus)$Group_combined))
comparisons <- combn(groups, 2, simplify = FALSE)
results_list <- list()

for (comp in comparisons) {
  a <- comp[1]; b <- comp[2]
  cname <- paste0(a, "_vs_", b)
  
  res_raw <- tryCatch(results(dds, contrast = c("Group_combined", a, b)), error = function(e) NULL)
  if (is.null(res_raw)) next
  
  res_raw_df <- as.data.frame(res_raw) %>% rownames_to_column("OTU")
  
  res_shrunk <- tryCatch(lfcShrink(dds, contrast = c("Group_combined", a, b), type = "normal"), error = function(e) NULL)
  
  if (!is.null(res_shrunk)) {
    res_shrunk_df <- as.data.frame(res_shrunk) %>% rownames_to_column("OTU") %>% select(OTU, log2FoldChange)
    res_combined <- res_raw_df %>%
      left_join(res_shrunk_df, by = "OTU", suffix = c("_raw", "")) %>%
      mutate(contrast = cname) %>%
      select(OTU, log2FoldChange, everything())
  } else {
    res_combined <- res_raw_df %>% mutate(log2FoldChange = log2FoldChange) %>% mutate(contrast = cname)
  }
  
  results_list[[cname]] <- res_combined
}

deseq_all <- bind_rows(results_list)

# --- 9. Annotate genera ---
deseq_all <- left_join(deseq_all, tax_df_genus, by = "OTU")

# --- 10. Filter significant genera  ---

if(!"padj" %in% colnames(deseq_all)) {
  sig_deseq <- deseq_all %>% filter(!is.na(pvalue) & pvalue <= 0.05)
} else {
  sig_deseq <- deseq_all %>% filter(!is.na(padj) & padj <= 0.05)
}

sig_deseq <- sig_deseq %>% filter(!is.na(log2FoldChange), !is.na(Label)) %>% mutate(Label = as.character(Label))
sig_deseq_plot <- sig_deseq


severity_levels <- c("Control", "Mild", "Severe", "Deceased")

desired_contrasts <- c(
  paste0("male_", severity_levels, "_vs_female_", severity_levels),
  paste0("female_", severity_levels, "_vs_male_", severity_levels)
)


filtered_contrasts <- intersect(desired_contrasts, unique(sig_deseq_plot$contrast))

sig_deseq_plot <- sig_deseq_plot %>% filter(contrast %in% filtered_contrasts)

contrasts <- unique(sig_deseq_plot$contrast)


# --- 11a. Log2Fold Bar Plots ---
if (length(contrasts) == 0) {
  print("Skipping Bar Plots: No significant genera found for the Male vs Female same-severity comparisons.")
} else {
  for (c in contrasts) {
    df <- sig_deseq_plot %>% filter(contrast == c)
    
    # Extract groups and sex for dynamic labeling
    groups <- unlist(strsplit(c, "_vs_"))
    group1 <- groups[1]; group2 <- groups[2]
    sex1 <- strsplit(group1, "_")[[1]][1] 
    sex2 <- strsplit(group2, "_")[[1]][1]
    
    df <- df %>% mutate(Label = fct_reorder(Label, log2FoldChange),
                        direction = ifelse(log2FoldChange > 0, "Up", "Down"))
    
    df <- df %>% mutate(direction_label = case_when(
      direction == "Up" ~ paste("Higher in", sex1),
      direction == "Down" ~ paste("Higher in", sex2)
    ))
    
    color_mapping <- setNames(c("cornflowerblue", "plum2"),
                              c(paste("Higher in", sex1), paste("Higher in", sex2)))
    
    p <- ggplot(df, aes(x = Label, y = log2FoldChange, fill = direction_label)) +
      geom_col(show.legend = TRUE) +
      coord_flip() +
      theme_minimal() +
      labs(title = paste0("Significant Genera — ", c),
           x = "Genus",
           y = "log2 Fold Change",
           fill = "Abundance Difference") +
      scale_fill_manual(values = color_mapping) +
      theme(axis.text.y = element_text(size = 9),
            axis.text.x = element_text(size = 9),
            plot.title = element_text(hjust = 0.5))
    
    ggsave(filename = paste0("DESeq2_Genus_", c, ".png"), plot = p, width = 10, height = 6, dpi = 300)
  }
}

   


#### ---- Integrated Absolute abundance plot ---- ####
# --- 1. Load library ---
# Loaded in ISA code

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

# --- 9. Plot Log10 Stacked Bar Chart (Updated Title) ---
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
    title = "Absolute Abundance of Significant (Log10 Stack)",
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
ggsave("stacked_genus_log10_absolute_abundance.png", log_stacked_plot, width = 12, height = 9, dpi = 300)


