############################################################################
## ALPHA DIVERSITY BY AGE BINS & SETBACKS COMPARISON
## ==================================================
## Analyzes alpha diversity across age bins and compares setbacks vs HTP.
##
## Output figures:
##   Supplementary Figure S6b  |  Supplementary Figure S12
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(tidyverse)
library(ggplot2)
library(microViz)
library(dplyr)
library(phyloseq)
library(vegan)
library(tidyr)
library(microbiome)

# ---- COLOR PALETTE ----
pastel_colors <- c("#6B96B6", "#8FC3D4", "#C4E3F0", "#F8BDAF", "#FFD9AC",
                   "#8EB9B7", "#AEC89E", "#CBBE76", "#F5DDB7", "#75AEBD",
                   "#A6CFCF", "#E78F8D", "#FFC1BE", "#A39F96", "#D3A9B9",
                   "#FED4E3", "#C694B3", "#E7C6D5", "#BAAFA5", "#E4CEC2",
                   "lightgrey")

# ---- HELPER FUNCTIONS ----
age_levels <- function(){
  lvls <- c("0-3mo", "3-6mo", "6-9mo", "9-12mo", "12-15mo",
            "15-18mo", "18-21mo", "21-24mo", "Mother")
  return(lvls)
}

add_agecat_agglom <- function(ps){
  start_days <- seq(from = 0, to = 730, by = 92)
  age_lvls <- age_levels()

  ps_agecat <- ps %>% ps_mutate(Day = as.integer(Day)) %>% ps_mutate(
    age_cat = case_when(
      Day >= start_days[1] & Day < start_days[2] ~ age_lvls[1],
      Day >= start_days[2] & Day < start_days[3] ~ age_lvls[2],
      Day >= start_days[3] & Day < start_days[4] ~ age_lvls[3],
      Day >= start_days[4] & Day < start_days[5] ~ age_lvls[4],
      Day >= start_days[5] & Day < start_days[6] ~ age_lvls[5],
      Day >= start_days[6] & Day < start_days[7] ~ age_lvls[6],
      Day >= start_days[7] & Day < start_days[8] ~ age_lvls[7],
      Day >= start_days[8] & Day < 730          ~ age_lvls[8],
      .default                                  = age_lvls[9]
    )) %>% ps_mutate(age_cat = factor(age_cat, levels = age_lvls))

  ps_genus <- tax_glom(ps_agecat, taxrank = "Genus")
  return(ps_genus)
}

calc_alpha_manual <- function(pseq, age_lvls){
  alpha_man <- microbiome::alpha(pseq, index = c("Observed", "Shannon")) %>%
    as.data.frame() %>% mutate(Sample_sequencing_ID = rownames(.))

  meta <- data.frame(sample_data(pseq)) %>%
    select(Sample_sequencing_ID, Day, age_cat, Child_ID, sample_category)

  diversity_manual <- left_join(meta, alpha_man, by = "Sample_sequencing_ID")
  diversity_manual$age_cat <- factor(diversity_manual$age_cat, levels = age_lvls)
  return(diversity_manual)
}

plot_richness_manual <- function(diversity_manual, title){
  diversity_manual <- pivot_longer(diversity_manual, cols = c(observed, diversity_shannon))
  plt <- ggplot(data = diversity_manual, aes(x = age_cat, y = value)) +
    geom_jitter(color = "darkgrey", size = 2, alpha = 0.15) +
    geom_boxplot(aes(fill = age_cat), outlier.shape = NA) +
    scale_fill_manual(values = pastel_colors) +
    ggtitle(title) + theme_light() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
  return(plt + facet_grid(name ~ ., scales = "free"))
}

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")

# data processing
bbgut_meta_final<-sample_data(bbgut_ps) %>% data.frame()
bbgut_ab<-otu_table(bbgut_ps) %>% data.frame()
bbgut_taxa_ASVtable_final2<-tax_table(bbgut_ps) %>% data.frame()

start_days <- seq(from=0,to=730,by=92)

baby_ids <- bbgut_meta_final %>%filter(fecalsample_type=='infant') %>% rownames(.)

mom_ids <- bbgut_meta_final %>% filter(fecalsample_type=='maternal') %>% rownames(.)
seq_ids <- c(baby_ids,mom_ids)

combined_meta <- bbgut_meta_final %>% filter(Sample_sequencing_ID %in% seq_ids)
combined_ab <- bbgut_ab %>% filter(rownames(.) %in% seq_ids)

ps_combined <- phyloseq(otu_table(combined_ab, taxa_are_rows=F), 
                        sample_data(combined_meta), 
                        tax_table(as.matrix(bbgut_taxa_ASVtable_final2))) %>% 
                        prune_taxa(taxa_sums(.)>0,.)

ps_agecat <- add_agecat_agglom(ps_combined)

diversity_all <- calc_alpha_manual(ps_agecat)

all_plot <- plot_richness_manual(diversity_all,"Alpha Diversity - all samples")
all_plot

ggsave(filename = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Supplementary_Figure_S6b_alpha_diversity_time_ALL.pdf"),
       plot = all_plot, height = 6, width = 7, units = "in")

####################################################################################
# Alpha diversity comparison Setbacks vs Rest of the samples
# in order to make a meaningful comparison, here i take for each 
# setback the closest HTP sample before that to do a pairwise comparison per child 
# and per setback event
# Publication *Supplementary Figure S12*
####################################################################################

#read setback samples
setbacks<-openxlsx::read.xlsx("Data/setbacks.xlsx")
setback_samples<-setbacks$Sample_sequencing_ID

#filter phyloseq object to keep only infant samples
bbgut_ps2<-bbgut_ps %>% ps_filter(Sample_sequencing_ID %in% c(baby_ids))

ps_agecat2 <- add_agecat_agglom(bbgut_ps2)

#calculate alpha diversity
diversity_all2 <- calc_alpha_manual(ps_agecat2)

#create new columns
diversity_all2 <- diversity_all2 %>%
  filter(Sample_sequencing_ID %in% c(baby_ids)) %>%
  mutate(setback=ifelse(Sample_sequencing_ID %in% 
                     setback_samples, "setback","all_samples"))


#make a function that picks only the sample before the setback 
pick_htp<-function(data,setback_id){
  setback_data <- data %>% filter(Sample_sequencing_ID == setback_id)
  child_data<-data %>% filter(Child_ID == setback_data$Child_ID)
  setback_day<-setback_data$Day
  htp_data<- child_data %>% filter(sample_category=="HTP") 
  closest_htp<-htp_data %>% filter(Day < setback_day) %>% 
    arrange(Day) %>% 
    pull(Sample_sequencing_ID) %>%
    tail(.,n=1)
  return(closest_htp)
}

#get matching HTP samples
comparison_div_data<-c()
for (s in setback_samples){
  htp_match<-pick_htp(diversity_all2,s)
  comparison_div_data<-c(comparison_div_data,htp_match)
} 

#create dataframe and deduplicate two setbacks that have the same HTP sample
com<-diversity_all2 %>% filter(Sample_sequencing_ID %in% c(comparison_div_data,setback_samples)) %>% filter(Sample_sequencing_ID != "I18_164")

pairs<-data.frame(setback_id=setback_samples,htp_id = comparison_div_data) %>% 
  mutate(event_id = setback_id) %>% pivot_longer(!event_id, names_to='type', values_to = "Sample_sequencing_ID") %>% filter(event_id != "I18_164")

#create final df 
comparison_div_pairs<-left_join(com,pairs,by="Sample_sequencing_ID")

#plot observed diversity
obs_set<-comparison_div_pairs %>% ggplot(aes(x=type, y=observed,fill=type))+  geom_boxplot() +
  scale_fill_manual(values = c("setback_id" = "#E78F8D", "htp_id" = "#8EB9B7")) +
  ylab("Observed diversity")+
  theme_light()
obs_set
ggsave(filename = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Supplementary_Figure_10_observed_diversity_setbacks_vs_HTP.pdf"),
       plot = obs_set, height = 5, width = 6, units = "in")

#plot shannon diversity
sh_set<-comparison_div_pairs %>% ggplot(aes(x=type, y=diversity_shannon,fill=type))+  geom_boxplot() +
  scale_fill_manual(values = c("setback_id" = "#E78F8D", "htp_id" = "#8EB9B7")) +
  ylab("Shannon diversity")+
  theme_light()
sh_set
ggsave(filename = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Supplementary_Figure_S12_shannon_diversity_setbacks_vs_HTP.pdf"),
       plot = sh_set, height = 5, width = 6, units = "in")


### turn wide for wilcoxon pairwise test each setback with each HTP sample before 
# the event per child

#observed
ob_comp<-comparison_div_pairs %>% select(observed,event_id, type) %>% pivot_wider(id_cols=event_id, names_from = type, values_from=observed)


wilcox_test_obs <- wilcox.test(ob_comp$htp_id, ob_comp$setback_id, paired = TRUE)
wilcox_test_obs
# Wilcoxon signed rank test with continuity correction
# 
# data:  ob_comp$htp_id and ob_comp$setback_id
# V = 428, p-value = 5.501e-06
# alternative hypothesis: true location shift is not equal to 0

#shannon
sh_comp<-comparison_div_pairs %>% select(diversity_shannon,event_id, type) %>% pivot_wider(id_cols=event_id, names_from = type, values_from=diversity_shannon)


wilcox_test_sh <- wilcox.test(sh_comp$htp_id, sh_comp$setback_id, paired = TRUE)
wilcox_test_sh
# Wilcoxon signed rank exact test
# 
# data:  sh_comp$htp_id and sh_comp$setback_id
# V = 404, p-value = 0.0001886
# alternative hypothesis: true location shift is not equal to 0