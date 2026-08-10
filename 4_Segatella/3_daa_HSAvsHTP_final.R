############################################################################
## DIFFERENTIAL ABUNDANCE ANALYSIS — HSA vs HTP CONTROLS
## =======================================================
## Paired differential abundance analysis comparing High Segatella Abundance (HSA)
## samples with their closest preceding healthy (HTP) controls.
##
## Output figures:
##   Supplementary Figure S17
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(phyloseq)
library(dplyr)
library(tidyr)
library(microbiome)
library(ggplot2)

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")
prev_events_all <- openxlsx::read.xlsx("Data/HSA_meta.xlsx")

#agglomerate to genus
genus_agglom<-tax_glom(bbgut_ps, taxrank="Genus")

genus_agglomcomp <- transform_sample_counts(genus_agglom, function(OTU) OTU/sum(OTU))
genus_agglom_abundance_comp <- psmelt(genus_agglomcomp)

#hsa samples
hsa<-prev_events_all$Sample_sequencing_ID

###############################################################
# 1) Metadata table
###############################################################
sample_meta <- genus_agglom_abundance_comp2 %>%
  distinct(Sample_sequencing_ID, Child_ID, Day, sample_category)

################################################################
# 2) Define HSA samples and HTP controls to use later
#    HTP must NOT include any HSA samples(some are HTP)
###############################################################
hsa_meta <- sample_meta %>%
  filter(Sample_sequencing_ID %in% hsa) %>%
  select(Sample_sequencing_ID, Child_ID, Day)

htp_meta <- sample_meta %>%
  filter(sample_category == "HTP") %>%
  filter(!Sample_sequencing_ID %in% hsa) %>%
  select(Sample_sequencing_ID, Child_ID, Day)

###############################################################
# 3) Match each HSA sample to its closest preceding HTP
# within the same child to use as a control check
#  for the composition before the high Segatella event
###############################################################
matched <- hsa_meta %>%
  inner_join(htp_meta, by = "Child_ID", suffix = c("_hsa", "_htp")) %>%
  mutate(time_diff = Day_hsa - Day_htp) %>%
  filter(time_diff > 0) %>%   # keep only true "before spike" samples
  group_by(Sample_sequencing_ID_hsa) %>%
  slice_min(order_by = time_diff, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    pair_id = Sample_sequencing_ID_hsa,
    HSA = Sample_sequencing_ID_hsa,
    HTP = Sample_sequencing_ID_htp
  )

###############################################################
# 4) Subset abundance table to only matched samples
###############################################################
keep_ids <- c(matched$HSA, matched$HTP)

htp_hsa_paired_meta <- meta_ps %>%
  filter(Sample_sequencing_ID %in% keep_ids) %>%
  select(Sample_sequencing_ID, Day, sample_category)

###############################################################
# 5) Create analysis dataset with group labels (HSA, HTP)
###############################################################
melted_htp_hsa <- genus_agglom_abundance_comp %>%
  filter(Sample_sequencing_ID %in% keep_ids) %>%
  mutate(Group = ifelse(Sample_sequencing_ID %in% hsa,
                        "HSA",
                        "HTPcontrol"))


###############################################################
# 6) Modify so that that each comparison pair has an ID
###############################################################
melted_htp_hsa2 <- melted_htp_hsa %>%
  mutate(
    pair_id = ifelse(
      Sample_sequencing_ID %in% hsa,
      Sample_sequencing_ID,
      matched$pair_id[match(Sample_sequencing_ID, matched$HTP)]
    )
  )

###############################################################
# 7) Prepare abundance table for statistical testing
###############################################################
daa_df <- melted_htp_hsa2 %>%
  mutate(Abundance_adj = ifelse(Abundance <= 0.005, 0, Abundance)) %>%
  select(Abundance_adj, Group, Genus, pair_id)

###############################################################
# 8) Collapse duplicated taxa within each pair
###############################################################
daa_df_clean <- daa_df %>%
  group_by(pair_id, Genus, Group) %>%
  summarise(Abundance_adj = sum(Abundance_adj), .groups = "drop")

###############################################################
# 9) Turn to wide format (for Wilcoxon test)
###############################################################
wide_df <- daa_df_clean %>%
  pivot_wider(
    names_from = Group,
    values_from = Abundance_adj
  )

###############################################################
# 10) Remove NAs or pairs with 0 abundance
###############################################################
wide_df <- wide_df %>%
  filter(!is.na(HSA) & !is.na(HTPcontrol)) %>%
  filter(HSA > 0 | HTPcontrol > 0)

###############################################################
# 11) Paired Wilcoxon test per genus
###############################################################
wilcox_results <- wide_df %>%
  group_by(Genus) %>%
  summarise(
    p = wilcox.test(HSA, HTPcontrol, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p.adj = p.adjust(p, method = "BH"),
    sig = ifelse(p.adj < 0.05, "significant", "non-significant")
  )

#### plot 

sig_genera <- wilcox_results %>% filter(sig=="significant") %>% pull(Genus)
sig_genera

group_plot <- ggplot(daa_df_clean %>% 
                        filter(Genus %in% sig_genera), 
                      aes(Genus, Abundance_adj, fill = Group)) +
  geom_boxplot(lwd = 0.2, outlier.color = "darkgrey", outlier.size = 0.7, outlier.alpha = 0.3) +
  labs(x = "", y = "Relative Abundance", fill = "Group") +
  theme_light() +
  scale_y_continuous(trans='log10')+
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, face = "italic", size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "darkgrey"),
    axis.ticks = element_line(size = 0.5, color = "darkgrey"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 1)
  ) +
  scale_fill_manual(values = c("#CBBED9" ,"#6B96B6"))

group_plot
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Fig_S17_HSAvsHTP.pdf", sep=""),width=9,height=6)
print(group_plot)
dev.off()
