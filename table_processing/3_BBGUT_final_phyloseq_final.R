############################################################################
## FINAL BBGUT PHYLOSEQ OBJECT CONSTRUCTION
## ==========================================
## Builds the final phyloseq object from decontaminated & rarefied data.
## Filters out Archaea and applies abundance threshold (0.0001 relative).
##
## Input:
##   table_processing/bbgut_data_decontamrarefied.RData  (from script 2)
## Output:
##   bbgut_phyloseq_final.RData (save manually after running)
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(phyloseq)
library(dplyr)

# ---- LOAD DATA ----
## Load decontaminated & rarefied data from step 2
load("table_processing/bbgut_data_decontamrarefied.RData")

# ---- FILTER OUT NON-BACTERIA ----
bbgut_taxa_ASVtable_final2 <- bbgut_taxa_ASVtable_final %>% filter(Kingdom != "Archaea")

# ---- BUILD PHYLOSEQ OBJECT (INFANT + MATERNAL) ----
bbgut_ab <- t(bbgut_ab_final) %>% as.data.frame()

bbgut_meta <- bbgut_meta_final %>% filter(fecalsample_type %in% c("infant", "maternal")) %>% drop_na(., Type)

## Anonymize sample IDs
bbgut_meta$Sample_sequencing_ID <- bbgut_meta$Sample_sequencing_ID_new
bbgut_meta$Sample_sequencing_ID_new <- NULL

bbgut_meta$Child_ID <- bbgut_meta$Child_ID_new
bbgut_meta$Child_ID_new <- NULL

bbgut_otu <- otu_table(bbgut_ab, taxa_are_rows = FALSE)
bbgut_meta_ps <- sample_data(bbgut_meta)

bbgut_taxa <- tax_table(as.matrix(bbgut_taxa_ASVtable_final2))

bbgut_ps <- phyloseq(bbgut_otu, bbgut_meta_ps, bbgut_taxa)

# ---- FILTER BY RELATIVE ABUNDANCE THRESHOLD (0.0001 ACROSS ALL SAMPLES) ----
bbgut_ab_final2 <- as.data.frame(t(otu_table(bbgut_ps)))
bbgut_ab_final2$total_Counts <- rowSums(bbgut_ab_final2)
bbgut_ab_final2$relative_totalabundance <- bbgut_ab_final2$total_Counts / sum(bbgut_ab_final2$total_Counts)

bbgut_ab_final_filtered <- bbgut_ab_final2 %>% filter(relative_totalabundance > 0.0001)
bbgut_ab_final_filtered$total_Counts <- NULL
bbgut_ab_final_filtered$relative_totalabundance <- NULL

bbgut_ab_final_filtered <- t(bbgut_ab_final_filtered) %>% as.data.frame()

## Replace OTU table
bbgut_otu <- otu_table(bbgut_ab_final_filtered, taxa_are_rows = FALSE)
# # View(bbgut_otu)

## Filter taxa table
bbgut_taxa_ASVtable_final2_filtered <- bbgut_taxa_ASVtable_final2 %>% filter(rownames(.) %in% colnames(bbgut_otu))
bbgut_taxa <- tax_table(as.matrix(bbgut_taxa_ASVtable_final2_filtered))

## Final phyloseq object
bbgut_ps <- phyloseq(bbgut_otu, bbgut_meta_ps, bbgut_taxa)

## Save final object: save(bbgut_ps, file = "Data/bbgut_phyloseq_final.RData")
