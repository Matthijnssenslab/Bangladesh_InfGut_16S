############################################################################
## DECONTAMINATION & RAREFACTION
## ================================
## Identifies contaminant ASVs using prevalence method (negative controls),
## rarefies to equal sequencing depth, then removes contaminants & zero-ASVs.
##
## Input:
##   table_processing/bbgut_ps_rdpannotated.RData  (from script 1)
## Output:
##   bbgut_data_decontamrarefied.RData (save manually after running)
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(ggplot2)
library(dplyr)
library(phyloseq)
library(vegan)

# ---- LOAD DATA ----
## Load RDP-annotated phyloseq object from step 1
load("table_processing/bbgut_ps_rdpannotated.RData")

# ---- PREPARE PHYLOSEQ OBJECT ----
bbgut_otu <- otu_table(bbgut_rdp_set19_genus_names_edited, taxa_are_rows = TRUE)

bbgut_meta <- bbgut_meta_final
bbgut_meta <- sample_data(bbgut_meta)

bbgut_taxa <- tax_table(bbgut_rdp_set19_genus_names_edited)

ps <- phyloseq(bbgut_otu, bbgut_meta, bbgut_taxa)
ps

# ---- 1. DECONTAMINATION ----
## Identify contaminant sequences without removing them immediately

library(decontam)

head(sample_data(ps))
df <- as.data.frame(sample_data(ps))
df$LibrarySize <- sample_sums(ps)
df <- df[order(df$LibrarySize), ]
ggplot(data = df, aes(x = Index, y = LibrarySize, color = fecalsample_type)) + geom_point()

## Identify contaminant nodes
sample_data(ps)$is.neg <- sample_data(ps)$fecalsample_type == "control"
contamdf.prev <- isContaminant(ps, method = "prevalence", neg = "is.neg", batch = "plate.Name", threshold = 0.1)
table(contamdf.prev$contaminant)
head(which(contamdf.prev$contaminant))

## Prevalence plot: positive vs negative samples
ps.pa <- transform_sample_counts(ps, function(abund) 1 * (abund > 0))
ps.pa.neg <- prune_samples(sample_data(ps.pa)$fecalsample_type == "control", ps.pa)
ps.pa.pos <- prune_samples(sample_data(ps.pa)$fecalsample_type %in% c("infant", "maternal"), ps.pa)
df.pa <- data.frame(
  pa.pos = taxa_sums(ps.pa.pos), pa.neg = taxa_sums(ps.pa.neg),
  contaminant = contamdf.prev$contaminant
)
ggplot(data = df.pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

## Contaminant ASV list
vector_contam <- contamdf.prev %>% dplyr::filter(contamdf.prev$contaminant == "TRUE")
vector_contam
vector_contam_names1 <- tibble::rownames_to_column(vector_contam, "NODES")
vector_cont_names <- vector_contam_names1$NODES

## Inspect contaminant taxonomy
contamonly <- bbgut_taxa[rownames(bbgut_taxa) %in% vector_cont_names, ]
# # contamonly %>% as.data.frame() %>% View()

# ## Save contamination report: openxlsx::write.xlsx(contamonly, "table_processing/contaminant_ASVs_decontam.xlsx")

# ---- 2. FILTER TO STOOL SAMPLES ONLY ----
bbgut_ab_table_ASV <- as.data.frame(bbgut_otu)

bbgut_stool_ab_samples <- bbgut_meta_final %>% filter(Type == "Sample") %>% drop_na(., Type)

bbgutstool_ab_table <- bbgut_ab_table_ASV %>% select(bbgut_stool_ab_samples$Sample_sequencing_ID_new)

# ---- 3. RAREFACTION ----
## Check read count distribution
list_unrarefied <- colSums(bbgutstool_ab_table)
mapped_read_counts <- data.frame(tot_reads = list_unrarefied)
ggplot(mapped_read_counts, aes(x = tot_reads)) + geom_histogram(bins = 100) +
  scale_x_continuous(breaks = seq(0, 250000, 10000))

## Rarefaction curve (sample first 150 ASVs)
abudances2rarefy <- t(otu_table(bbgut_ab_table_ASV, taxa_are_rows = TRUE))
abudances2rarefy <- as.data.frame(abudances2rarefy)
colnames(abudances2rarefy) <- 1:length(colnames(abudances2rarefy))
abudances2rarefy <- as.matrix(abudances2rarefy)
raremax <- min(rowSums(abudances2rarefy[1:150, ]))
rarecurve(abudances2rarefy, step = 50, label = FALSE, col = "lightblue")
abline(v = 16000)

# # View(mapped_read_counts)

## Rarefy to 15086 reads (remove samples below this threshold)
bbgutstool_ab_rarefied <- bbgutstool_ab_table[, colSums(bbgutstool_ab_table) >= 15086]

set.seed(3313)
bbgutstool_ab_rarefied_sum <- t(rrarefy(t(bbgutstool_ab_rarefied), sample = 15086))

bbgutstool_ab_rarefied <- as.data.frame(bbgutstool_ab_rarefied_sum)

# ---- REMOVE CONTAMINANTS AFTER RAREFACTION ----
bbgut_ab_table_rarefied_decontam <- bbgutstool_ab_rarefied %>% filter(!rownames(.) %in% vector_cont_names)

## Remove ASVs that became 0 after rarefaction
bbgut_ab_table_rarefied_decontam <- bbgut_ab_table_rarefied_decontam %>% filter(rowSums(bbgut_ab_table_rarefied_decontam) > 0)

bbgut_ab_final <- bbgut_ab_table_rarefied_decontam

## Filter taxa table to matching ASVs
asv_names <- rownames(bbgut_ab_final)
bbgut_taxa_ASVtable <- tax_table(bbgut_rdp_set19_genus_names_edited) %>% as.data.frame()
bbgut_taxa_ASVtable2 <- bbgut_taxa_ASVtable %>% filter(rownames(.) %in% asv_names)
bbgut_taxa_ASVtable_final <- bbgut_taxa_ASVtable2

## Filter metadata to matching samples
bbgut_meta_final2 <- bbgut_meta_final %>% filter(rownames(.) %in% colnames(bbgut_ab_final))
bbgut_meta_final <- bbgut_meta_final2

## Save environment: save(list = c("bbgut_ab_final", "bbgut_taxa_ASVtable_final", "bbgut_meta_final"), file = "table_processing/bbgut_data_decontamrarefied.RData")
