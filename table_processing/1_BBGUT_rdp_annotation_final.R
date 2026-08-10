############################################################################
## RDP ANNOTATION OF ASVs — TRAINING SET 19
## ==========================================
## Annotates all ASVs using the RDP classifier with training set 19.
## Handles unclassified taxa at Family, Order, Class, Genus levels.
##
## Input:
##   table_processing/BBGUT_phyloseq_preprocessed.RData
## Output:
##   bbgut_ps_rdpannotated.RData (save manually after running)
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dada2)
library(phyloseq)
library(dplyr)

# ---- LOAD DATA ----
## Load preprocessed phyloseq object from table_processing directory
load("table_processing/BBGUT_phyloseq_preprocessed.RData")

bbgut_ab_table_ASV <- as.data.frame(otu_table(bbgut_ps))

# ---- RE-ANNOTATE ASVs ----
ASV_seq <- taxa_names(bbgut_ps)

set.seed(1234)
ASV_seq_rdp_set19 <- assignTaxonomy(ASV_seq, "rdp_19_toGenus_trainset.fa.gz", multithread = TRUE)

bbgut_rdp_set19 <- merge_phyloseq(sample_data(bbgut_ps), tax_table(ASV_seq_rdp_set19), otu_table(bbgut_ps))

# ---- RENAME UNCLASSIFIED GENERA ----
mat <- as.data.frame(bbgut_rdp_set19@tax_table)

mat <- mat %>%
  mutate(
    Family = ifelse(is.na(Family), paste0("uncl_", Order), Family),
    Order  = ifelse(is.na(Order),  paste0("uncl_", Class), Order),
    Class  = ifelse(is.na(Class),  paste0("uncl_", Phylum), Class),
    Genus  = ifelse(is.na(Genus),  paste0("uncl_", Family), Genus)
  )

# ---- BUILD NEW PHYLOSEQ OBJECT ----
bbgut_rdp_set19_genus_names_edited <- phyloseq(
  otu_table(bbgut_rdp_set19, taxa_are_rows = FALSE),
  sample_data(bbgut_rdp_set19),
  tax_table(as.matrix(mat))
)

bbgut_rdp_set19_genus_names_edited

# # View(tax_table(bbgut_rdp_set19_genus_names_edited))

bbgut_ps <- bbgut_rdp_set19_genus_names_edited

# # View(tax_table(bbgut_ps))

## Save as environment: save(bbgut_ps, file = "table_processing/bbgut_ps_rdpannotated.RData")
