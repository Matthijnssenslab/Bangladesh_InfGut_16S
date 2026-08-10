############################################################################
## MICROBIOTA COMPOSITION COMPARISON: BBGUT vs BABEL
## ===================================================
## Compares prevalence of common bacterial genera over time between the
## Bangladeshi (BBGUT) and Belgian (BABEL) infant cohorts.
##
## Output figures:
##   Figure 2a
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(tibble)
library(microbiome)
library(microViz)

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")
load("Data/babel_phyloseq_final_newtaxonomy.RData")


##########################
# BBGUT data processing
##########################
# Filter infant samples, keep year 1
bbgut_ps_y1<-bbgut_ps %>% 
  ps_filter(fecalsample_type== "infant") %>%
  ps_filter(Day<=365)

#define healthy samples
HTP_samples<-sample_data(bbgut_ps_y1)$Sample_sequencing_ID[sample_data(bbgut_ps_y1)$sample_category == "HTP"] %>%na.omit()

#filter to keep only healthy first year samples
bbgut_ps_healthy_y1<-bbgut_ps_y1 %>%
  ps_filter(Sample_sequencing_ID %in% HTP_samples)

samples_y1_bbgut<-sample_data(bbgut_ps_healthy_y1)$Sample_sequencing_ID

# transform
bbgut_ps_abundance <- transform_sample_counts(bbgut_ps_healthy_y1, 
                                              function(OTU) OTU/sum(OTU))
#define infant ids and remove - so that you can match a pattern later
Infants<-unique(sample_data(bbgut_ps_healthy_y1)$Child_ID)[1:20] %>% str_replace("-","")
Infants2<-unique(sample_data(bbgut_ps_healthy_y1)$Child_ID)[1:20]

# agglomerate by genus                
genus_agglom <- tax_glom(bbgut_ps_healthy_y1, taxrank="Genus")

#top15 most abundant BBGUT genera in year 1
# get total abundance per genus
htp_filtered_ps_melt <- psmelt(genus_agglom)

tot_abund<-htp_filtered_ps_melt %>% 
  group_by(Genus) %>%
  dplyr::summarize(total_abundance=sum(Abundance)) %>% 
  # sort descending
  arrange(desc(total_abundance)) %>% 
  # add new column with relative ab
  mutate(relative_ab=round(total_abundance/sum(total_abundance)*100,2))

tot_abund %>% # View( )
Genera_tokeep <- head(tot_abund$Genus,15)


##########################
# BABEL data processing
##########################
# filter ps to keep only 1 year samples BABEL
babel_inf_y1<- babel_ps %>% 
  ps_filter(InfantID !="S011") %>% #S011 not considered in this study
  ps_filter(X.days <= 365)

#define healthy samples
LDA_samples<-sample_data(babel_inf_y1)$Sample.ID2[sample_data(babel_inf_y1)$LDA == 1] %>%na.omit()

#filter to keep healthy samples
htp_filtered_babel_ps  <- ps_filter(babel_inf_y1, 
                                    Sample.ID2 %in% 
                                      LDA_samples)

# agglomerate by genus                
babel_gen_ps <- tax_glom(htp_filtered_babel_ps, taxrank="Genus") 

#top15 BABEL genera in year 1
htp_filtered_babel_ps_melt <- psmelt(babel_gen_ps)

tot_abund_bab<-htp_filtered_babel_ps_melt %>% 
  group_by(Genus) %>%
  dplyr::summarize(total_abundance=sum(Abundance)) %>% 
  # sort descending
  arrange(desc(total_abundance)) %>% 
  # add new column with relative ab
  mutate(relative_ab=round(total_abundance/sum(total_abundance)*100,2))

tot_abund_bab %>% # View( )
Genera_tokeep_bab <- head(tot_abund_bab$Genus,15)

#all common genera (top15 bbgut + top15 babel)
comp_genera<-purrr::reduce(list(Genera_tokeep,Genera_tokeep_bab),union)

###########################################
## Mean prevalence of most abundant genera
## BBGUT vs BABEL
## Publication Figure 2a 
###########################################

#Rel. abundance ps
genus_agglom_rel<- genus_agglom %>% 
  # calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU))

babel_gen_ps_rel<- babel_gen_ps %>% 
  # calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU))

###combine cohorts into 1
merged_gen_ps<-merge_phyloseq(genus_agglom_rel,babel_gen_ps_rel)

#merge infant IDs
sample_data(merged_gen_ps)$InfantID_merged <- as.character(coalesce(sample_data(merged_gen_ps)$InfantID, sample_data(merged_gen_ps)$Child_ID))

#merge age in days
sample_data(merged_gen_ps)$Age_merged <- as.character(coalesce(sample_data(merged_gen_ps)$X.days, sample_data(merged_gen_ps)$Day))

#merge age in months
get_num_month <- function(days){
  return(ifelse(days<=10,0, ((days%/%31)+1) ))
}

# modify month to be days divided by 31
new_months_merged<- as.numeric(sample_data(merged_gen_ps)$Age_merged) %>% get_num_month() %>% as.factor()

sample_data(merged_gen_ps)$Month_merged <- new_months_merged

#merge sample IDs
sample_data(merged_gen_ps)$Sample_IDs_merged<-rownames(sample_data(merged_gen_ps))

#make new column in data that assigns a cohort name per sample
sample_data(merged_gen_ps)$Cohort <- ifelse(startsWith(sample_data(merged_gen_ps)$InfantID_merged, "I"), "BBGUT", ifelse(startsWith(sample_data(merged_gen_ps)$InfantID_merged, "S"), "BABEL", NA))

#check data is correct
sample_data(merged_gen_ps) %>% 
  data.frame() %>% 
  filter(Cohort=='BBGUT') %>%
  pull(InfantID_merged) %>% #gets values from a column
  unique()

sample_data(merged_gen_ps) %>% 
  data.frame() %>% 
  filter(Cohort=='BABEL') %>%
  pull(InfantID_merged) %>%
  unique()

#keep genera to compare in the object and transform to relative abundances
 merged_gen_ps_toplot<-merged_gen_ps %>%
  # 3. select only top 15 genus
  tax_select(comp_genera, ranks_searched = "Genus", 
             strict_matches=TRUE,n_typos = 0) 
 
 merged_gen_ps_toplot_melt<-ps_melt(merged_gen_ps_toplot)
 
# Heatmap plot

grad_colors <- c("#E6F4FE","#8FB6D1","#6D8FA9", "#4A6980","#284258","#061B2F")

order<-c("Bifidobacterium" , 
"Staphylococcus", 
"Bacteroides",   
"Dorea",
"Megasphaera" ,
"Veillonella"  ,
"Enterococcus",
"Enterobacter" ,
"Escherichia/Shigella",
"Phocaeicola",
"Sutterella" ,
"Parabacteroides",
"Thomasclavelia",
"Clostridium_sensu_stricto",
"Haemophilus",
"Parolsenella",
"Ligilactobacillus", 
"Limosilactobacillus",
"Segatella",
"Lactobacillus",  
"Streptococcus"
) 
                          

#mean prevalence by month 
p<-merged_gen_ps_toplot_melt %>%
  mutate(Prevalence = ifelse(Abundance >= 0.005, 1, 0)) %>%
  group_by(Genus,Month_merged,Cohort) %>%
  dplyr::summarise(mean_prevalence=mean(Prevalence)) %>% 
  mutate(
    Month_merged = as.factor(Month_merged), 
    Genus = factor(Genus, levels = order), 
    Cohort = as.factor(Cohort), 
    mean_prevalence = as.numeric(mean_prevalence)# create new column with presence/absence
  ) %>% 
  ggplot(aes(x = Month_merged, y = Genus, fill = mean_prevalence)) +
  geom_tile(color=NA) +
  facet_grid(rows = vars(Genus), cols = vars(Cohort), scales = "free_y", space = "free_x") + 
  scale_fill_gradientn(colors = grad_colors)+
  labs(
    title = "Mean prevalence of top genera per month",
    x = "Month",
    y = NULL, # Remove duplicate y-axis labels
    fill = "Prevalence"
  ) +
  theme_minimal() +
  scale_x_discrete(expand = c(0, 0)) + # Remove extra space on the x-axis (horizontal white space)
  scale_y_discrete(expand = c(0, 0)) + 
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1),
    axis.text.y = element_text(face = "italic"),
    strip.text.y = element_blank(), 
    panel.spacing = unit(0,"cm"), 
    strip.text = element_text(size = 10), 
    panel.grid = element_blank() 
  )
p
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_2a_Meanprevalence_topgen_BABELvsBBGUT.pdf", sep=""),
    width=8,height=6)
print(p)
dev.off()
