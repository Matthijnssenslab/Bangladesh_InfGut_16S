############################################################################
## SEGATELLA ASV ABUNDANCE ACROSS HEALTH GROUPS
## ==============================================
## Compares individual Segatella ASV abundance between disease, healthy, and maternal samples.
## Identifies ASVs more prevalent in specific health groups.
##
## Output figures:
##   Supplementary Figure S16
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(phyloseq)
library(tidyverse)
library(microViz)
library(r2r)
library(broom)

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")

## Process data
bbgut_asv_table<-tax_table(bbgut_ps) %>% as.data.frame()

bbgut_ps_abundance <- transform_sample_counts(bbgut_ps, function(OTU) OTU/sum(OTU))
bbgut_abundance_melt<-ps_melt(bbgut_ps_abundance)

##############################
# OTU - Genus_NUMBER hashmap
##############################

#BBGUT
bbgut_unique_OTU_Genus_df <- bbgut_abundance_melt %>% 
  select(OTU,Genus) %>%
  unique() %>% 
  group_by(Genus) %>%
  mutate(OTU_short = paste0(Genus, "_", row_number())) 

bbgut_otu_genus_hashmap <- hashmap()
bbgut_otu_genus_hashmap[bbgut_unique_OTU_Genus_df$OTU] <- bbgut_unique_OTU_Genus_df$OTU_short

#change otu sequnce into segatella_1,2,3,... (short name)
bbgut_abundance_melt_otushort<-bbgut_abundance_melt %>% 
  as.data.frame() %>%
  mutate(OTU_short = unlist(bbgut_otu_genus_hashmap[OTU]))


########################################################
# Find if some ASVs are more prevalent in 
# ADHC disease or HTP samples
# *Supplementary Figure 14*
########################################################

seg_sampletype<-read.csv("Data/segatella_acrossclusters_healthstatus.csv")
seg_highabund<-openxlsx::read.xlsx("Data/HSA_meta.xlsx")

seg_sampletype<-seg_sampletype %>% filter(Presence==1) %>% select(-X)

# #find which ASVs are present in the disease, healthy, maternal samples
disease<-seg_sampletype %>% filter(Group=="disease") %>% pull(Sample_sequencing_ID)
mother<-seg_sampletype %>% filter(Group=="maternal") %>% pull(Sample_sequencing_ID)
health<-seg_sampletype %>% filter(Group=="healthy") %>% pull(Sample_sequencing_ID)

#group by health status only relevant samples that have Segatella
segatella_grouped <- bbgut_abundance_melt_otushort %>%
  filter(Genus == "Segatella") %>%
  #0.05% presence threshold for ASV
  filter(Abundance >= 0.0005) %>%
  filter(Sample_sequencing_ID %in%seg_sampletype$Sample_sequencing_ID) %>%
  mutate(Group = case_when(
    Sample_sequencing_ID %in% disease ~ "disease",
    Sample_sequencing_ID %in% health ~ "healthy",
    Sample_sequencing_ID %in% mother ~ "mother"
  )) %>%
  select(OTU,Abundance, Sample_sequencing_ID, Group, OTU_short, Child_ID) %>%
  filter(!is.na(Group))

#order
ordered_otus <- paste0("Segatella_", 1:14)

#Set OTU_short as a factor with the levels
segatella_grouped$OTU_short <- factor(segatella_grouped$OTU_short, levels = ordered_otus)


seg_forwilcox<- segatella_grouped %>%
  #keep HSA samples
  filter(Sample_sequencing_ID %in% seg_highabund$Sample_sequencing_ID) %>%
  filter(Group %in% c("disease", "healthy"))%>%
  group_by(OTU_short, Child_ID, Group) %>%
  dplyr::summarise(Abundance_median=median(Abundance))

#define the order
ordered_otus <- paste0("Segatella_", 1:14)

#set OTU_short as a factor with the levels
seg_forwilcox$OTU_short <- factor(seg_forwilcox$OTU_short, levels = ordered_otus)

seg_wide<-seg_forwilcox %>% pivot_wider(id_cols=c(OTU_short, Child_ID), names_from = Group, values_from = Abundance_median)

otus<-unique(seg_wide$OTU_short)

wil_results2 <- seg_wide %>%
  # Replace NA with 0 in disease and healthy columns
  mutate(
    disease = ifelse(is.na(disease), 0, disease),
    healthy = ifelse(is.na(healthy), 0, healthy)
  ) %>%
  group_by(OTU_short) %>% 
  #filter groups to keep only those with > 2 paired (non-NA) observations
  filter(sum(!is.na(disease) & !is.na(healthy)) > 2) %>% 
  summarise(
    wilcox = list(wilcox.test(disease, healthy, paired = TRUE))
  ) %>% 
  mutate(tidied = map(wilcox, tidy)) %>%
  unnest(tidied) %>% 
  select(OTU_short, statistic, p.value, method) %>%
  mutate(p.adj = p.adjust(p.value, method = "BH"))

#plot
fplot<-seg_forwilcox %>%
  filter(OTU_short %in% wil_results2$OTU_short) %>%
  ggplot(aes(x=factor(OTU_short), 
             y=Abundance_median, fill = as.factor(Group))) +
  geom_boxplot() +
  theme_light() +
  scale_fill_manual(values= c(disease="#CBBED9", healthy="#6B96B6"))+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))
fplot
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_figure_S16_segatella_asvs.pdf", sep=""),width=7,height=5)
print(fplot)
dev.off()
