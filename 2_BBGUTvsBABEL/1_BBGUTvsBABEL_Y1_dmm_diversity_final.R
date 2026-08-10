############################################################################
## GUT MICROBIOTA DEVELOPMENT COMPARISON: BBGUT vs BABEL (YEAR 1)
## ================================================================
## DMM clustering for BABEL cohort, alpha diversity comparison, and beta
## diversity comparison between the Bangladeshi and Belgian infant cohorts.
##
## Output figures & tables:
##   Figure 2c  |  Supplementary Figure S9  |  Supplementary Tables S4, S5
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(tibble)
library(reshape2)
library(stringr)
library(phyloseq)
library(vegan)
library(tidyr)
library(synchrony)
library(rlist)
library(DirichletMultinomial)
library(parallel)
library(microbiome)
library(ggplot2)
library(microViz)

# ---- COLOR PALETTES ----
pastel_colors <- c("#6B96B6", "#8FC3D4", "#C4E3F0", "#F8BDAF", "#FFD9AC",
                   "#8EB9B7", "#AEC89E", "#CBBE76", "#F5DDB7", "#75AEBD",
                   "#A6CFCF", "#E78F8D", "#FFC1BE", "#A39F96", "#D3A9B9",
                   "#FED4E3", "#C694B3", "#E7C6D5", "#BAAFA5", "#E4CEC2",
                   "lightgrey")
pastel_colors_2 <- c("#BFE1ED", "#7FAFBF", "#778899")
pastel_colors_4 <- c("#BFE1ED", "#7FAFBF", "#778899", "#99A9BD")

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")
load("Data/babel_phyloseq_final_newtaxonomy.RData")



########################
# BBGUT data processing
########################

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

#load BBGUT data with cluster information from initial DMM (full dataset)
#also provided in Data folder
per_sample_stats_bbgut<-read.csv("Data/Per_sample_stats.csv")

clust_bbgut<-per_sample_stats_bbgut %>% select(X, Cluster)
clust_bbgut$Sample_sequencing_ID<-clust_bbgut$X
clust_bbgut$X<-NULL

#filter to keep y1 samples only
clust_bbgut2<-clust_bbgut %>% filter(Sample_sequencing_ID %in% samples_y1_bbgut)

#merge ps object with cluster info
meta_bbgut<-sample_data(bbgut_ps_healthy_y1) %>% data.frame() %>% 
  left_join(., clust_bbgut2, by="Sample_sequencing_ID") %>% tibble::column_to_rownames(.,"Sample_sequencing_ID")

meta_bbgut$Sample_sequencing_ID<-rownames(meta_bbgut)
sample_data(bbgut_ps_healthy_y1)<-meta_bbgut
meta_bbgut %>% group_by(Cluster) %>% tally()

# transform rel. abundance
bbgut_ps_abundance <- transform_sample_counts(bbgut_ps_healthy_y1, 
                                              function(OTU) OTU/sum(OTU))
#define infant ids and remove - so that you can match a pattern later
Infants<-unique(sample_data(bbgut_ps_healthy_y1)$Child_ID)[1:20] %>% str_replace("-","")
Infants2<-unique(sample_data(bbgut_ps_healthy_y1)$Child_ID)[1:20]

# agglomerate by genus                
genus_agglom <- tax_glom(bbgut_ps_healthy_y1, taxrank="Genus")

################################################
#### 1. DMM for BABEL dataset ####
## 1 year samples =<365 days, 7 infants
################################################

# set seed
set.seed(1)
seeds <- round(runif(20, 0, .Machine$integer.max))
# parallel dmn runs

run_dmn <- function(maxnumberofclusters,pseq,seed){
  set.seed(1)
  fit <- mclapply(1:maxnumberofclusters, dmn,
                  count = otu_table(pseq), verbose=TRUE,
                  seed = seed,
                  mc.cores = 4)
  result <- reshape2::melt(t(as.data.frame(sapply(fit, goodnessOfFit))))
  result$randomseed <- seed
  result$Var1 <- result$Var1 %>% str_replace("V","") %>% as.integer()
  colnames(result) <- c("numCluster","score",'value','seed')
  return(result)
}

## BABEL data processing 

# filter ps to keep only 1 year samples BABEL. (there are 3 samples collected after
# day 365, remove for comparison)
babel_inf_y1<- babel_ps %>% 
  ps_filter(InfantID !="S011") %>% #removed S011 as an outlier (see Beller et al. 2021)
  ps_filter(X.days <= 365)

#define healthy samples
LDA_samples<-sample_data(babel_inf_y1)$Sample.ID2[sample_data(babel_inf_y1)$LDA == 1] %>%na.omit()

# agglomerate by genus                
babel_gen_ps <- tax_glom(babel_inf_y1, taxrank="Genus")

# set seed
set.seed(1)
seeds <- round(runif(20, 0, .Machine$integer.max))


list_seedsBAB <- list()
set.seed(1)
for (seed in 1:20){
  print(paste('seed: ',seeds[seed]))
  #set.seed(1)
  fit <- run_dmn(10,babel_gen_ps,seeds[seed])
  list_seedsBAB[[seed]]  <-  fit
}
DMM_out_allBAB <- list.rbind(list_seedsBAB)

# check per seed where the minimum is reached
DMM_out_min_BICBAB <- DMM_out_allBAB %>% 
  filter(score=='BIC') %>% 
  group_by(seed) %>% 
  slice_min(value) %>%
  as.data.frame()

# get best seed and seed plots
best_seedBAB <- DMM_out_min_BICBAB %>% filter(value == min(value)) %>% pull(seed)

numClust_bestseedBAB<-DMM_out_min_BICBAB %>% 
  filter(value == min(value)) %>% 
  pull(numCluster) %>% 
  as.numeric()

plot_min_BIC_clustersBAB <- ggplot(DMM_out_min_BICBAB, aes(x= numCluster))+
  geom_histogram()
plot_min_BIC_clustersBAB
plot_best_seed_clustersBAB <- ggplot(DMM_out_allBAB %>% 
                                    filter(score=='BIC' & seed==best_seedBAB), 
                                  aes(x=numCluster, y=value)) +
  geom_point()+
  geom_path()+
  scale_x_continuous(breaks = c(1:10))+
  xlab("Number of Dirichlet Components") + 
  ylab("Model fit (BIC)")
plot_best_seed_clustersBAB
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "SPLM_BIC_best_clusternumber_BABEL_y1.pdf", sep=""),width=7,height=6)
print(plot_best_seed_clustersBAB)
dev.off()


#### Best fit stats ####
fit_optimalBAB <- dmn(otu_table(babel_gen_ps), numClust_bestseedBAB , 
                   verbose = F , seed = best_seedBAB) 
sample_data_babel_gen_agglom <- sample_data(babel_gen_ps) %>% data.frame()
sample_data_babel_gen_agglom$ClusterNum <- as.character(
  mixture(fit_optimalBAB, assign=TRUE))
sample_data_babel_gen_agglom$Cluster <- chartr("123","ABC",
                                         as.character(sample_data_babel_gen_agglom$ClusterNum))
sample_data(babel_gen_ps)<-sample_data_babel_gen_agglom


# statistics per cluster
per_clust_statsBAB <- mixturewt(fit_optimalBAB)
rownames(per_clust_statsBAB) <- c(LETTERS[1:numClust_bestseedBAB])
num_per_clustBAB <- table(sample_data_babel_gen_agglom$Cluster)
per_clust_statsBAB['A','cluster'] <- num_per_clustBAB['A']
per_clust_statsBAB['B','cluster'] <- num_per_clustBAB['B']
per_clust_statsBAB['C','cluster'] <- num_per_clustBAB['C']

#write.csv(per_clust_statsBAB, "Per_cluster_stats_BABEL_y1.csv")

# per sample statistics
per_sample_statsBAB <- as.data.frame(
  cbind(mixture(fit_optimalBAB), 
        mixture(fit_optimalBAB, assign=TRUE)))
colnames(per_sample_statsBAB) <- c(LETTERS[1:numClust_bestseedBAB],'selection')
per_sample_statsBAB$Cluster <- chartr("123","ABC",
                                   as.character(per_sample_statsBAB$selection))
per_sample_statsBAB <- per_sample_statsBAB %>% mutate(selection.prob=pmax(A,B,C))
per_sample_statsBAB %>% group_by(Cluster) %>% tally()

#write.csv(per_sample_statsBAB, "Per_sample_stats_BABEL_y1.csv")


#Kendall test for babel clusters
# kendalltable
meta_babel<-data.frame(sample_data(babel_gen_ps))
meta_day<- meta_babel %>% select(Sample.ID2,X.days,InfantID)
per_sample_statsBAB$Sample.ID2<-rownames(per_sample_statsBAB)
babel_clustersandage<-left_join(per_sample_statsBAB,meta_day, by="Sample.ID2")

babel_healthy<-babel_clustersandage %>% filter(Sample.ID2 %in% LDA_samples)

cluster_appearanceball <- babel_healthy %>%
  data.frame(.) %>%
  group_by(InfantID, Cluster) %>%
  dplyr::summarise(appearance_day = min(X.days)) %>% 
  pivot_wider(id_cols = InfantID,names_from = Cluster,values_from = appearance_day) %>%
  column_to_rownames("InfantID")
  
KENDALOUTCOMEall <- kendall.w(t(cluster_appearanceball), 
                           nrands = 10000, type = 1, quiet = T)
KENDALOUTCOME_tableall <- as.data.frame(matrix(nrow = 4, ncol=2))
KENDALOUTCOME_tableall[,1] <- c("Kendall's W (uncorrected for ties)",
                             "Kendall's W (corrected for ties)",
                             "Spearman's ranked correlation",
                             "Kendall's W p-value (one-tailed test [greater])")
KENDALOUTCOME_tableall[,2] <- c(KENDALOUTCOMEall$w.uncorrected, KENDALOUTCOMEall$w.corrected, 
                             KENDALOUTCOMEall$spearman.corr, KENDALOUTCOMEall$pval.rand)
colnames(KENDALOUTCOME_tableall) <- c('Kendall test, all infants','value')


############################################################################
## 2. Alpha diversity measures (observed richness, Shannon diversity) 
## of the each cohort per over time
## Publication *Figure 2c*
## Supplementary Table S5
#############################################################################

#add clusters to babel
babel_ps_healthy<-babel_inf_y1 %>%
  ps_filter(LDA==1)

clust_bab<-per_sample_statsBAB %>% 
  #tibble::rownames_to_column(.,"Sample.ID2") %>% 
  select(Sample.ID2, Cluster)
meta_bab<-sample_data(babel_ps_healthy) %>% data.frame() %>% 
  left_join(., clust_bab, by="Sample.ID2") %>% tibble::column_to_rownames(.,"Sample.ID2")
meta_bab$Sample.ID2<-rownames(meta_bab)
meta_bab$Day<-meta_bab$X.days
sample_data(babel_ps_healthy)<-meta_bab

#merge phyloseq objects
ps_combined<-merge_phyloseq(bbgut_ps_healthy_y1,babel_ps_healthy)
ps_combined_gen<-tax_glom(ps_combined, taxrank="Genus")


#merge sample IDs
sample_data(ps_combined_gen)$Sample_IDs_merged<-rownames(sample_data(ps_combined_gen))
# 
# #make new column in data that assigns a cohort name per sample
sample_data(ps_combined_gen)$Cohort <- ifelse(startsWith(sample_data(ps_combined_gen)$Sample_IDs_merged, "I"), "BBGUT", ifelse(startsWith(sample_data(ps_combined_gen)$Sample_IDs_merged, "S"), "BABEL", NA))
#merge age in days
sample_data(ps_combined_gen)$Age_merged <- as.character(coalesce(sample_data(ps_combined_gen)$X.days, sample_data(ps_combined_gen)$Day))

#differentiate the cohorts in each cluster per sample
ps_combined_cluster<-data.frame(sample_data(ps_combined_gen))
ps_combined_cluster$Cluster_Cohort <- paste(ps_combined_cluster$Cluster, ps_combined_cluster$Cohort, sep = "_")
ps_combined_cluster <-sample_data(ps_combined_cluster)

ps_combined_gen2<-phyloseq(ps_combined_cluster,otu_table(ps_combined_gen),tax_table(ps_combined_gen))


# Alpha diversity with time bins

### Add age category ###
age_levels <- function(){
  lvls <- c("0-3mo",
            "3-6mo",
            "6-9mo",
            "9-12mo")
  return(lvls)
}

add_agecat_agglom <- function(df){
  start_days <- seq(from=0,to=365,by=91)
  age_lvls <- age_levels() 
  
  ps_agecat <- df %>% ps_mutate(Day=as.integer(Day)) %>% ps_mutate(
    age_cat = case_when(
      Day >= start_days[1] & Day < start_days[2] ~ age_lvls[1],
      Day >= start_days[2] & Day < start_days[3] ~ age_lvls[2],
      Day >= start_days[3] & Day < start_days[4] ~ age_lvls[3],
      Day >= start_days[4] & Day < 365 ~ age_lvls[4],
      .default = age_lvls[4]
    )) %>% ps_mutate(age_cat=factor(age_cat,levels=age_lvls))
  
  #ps_genus <- tax_glom(ps_agecat,taxrank = "Genus")
  return(ps_agecat)
}

### Calculate alpha diversity for age categories ###
calc_alpha_manual_age<- function(pseq,age_lvls){
  alpha_man <- microbiome::alpha(pseq,index = c("Observed","Shannon")) %>% 
    as.data.frame() %>% mutate(Sample_IDs_merged=rownames(.))
  
  meta <- data.frame(sample_data(pseq)) %>% select(Sample_IDs_merged,Cohort,Day,age_cat)
  
  diversity_manual <- left_join(meta,alpha_man,by="Sample_IDs_merged")
  diversity_manual$age_cat <- factor(diversity_manual$age_cat, levels=age_lvls)
  
  return(diversity_manual)
}

### plotting function ###
plot_richness_manual <- function(diversity_manual,title){
  diversity_manual <- pivot_longer(diversity_manual,cols = c(observed,diversity_shannon))
  set.seed(123)
  plt <- ggplot(data=diversity_manual, aes(x=age_cat_cohort,y=value)) + 
    geom_jitter(color="darkgrey", size=2, alpha=0.15) +
    geom_boxplot(fill="#6B96B6",outlier.shape = NA) +  
    ggtitle(title)  + theme_light() +
    theme(axis.text.x = element_text(
      angle = 45, vjust = 1, hjust=1))
  return(plt + facet_grid(name ~ ., scales='free'))
}  

#run functions
ps_agecat <- add_agecat_agglom(ps_combined_gen2)

diversity_coh <- calc_alpha_manual_age(ps_agecat)
diversity_coh$age_cat_cohort<-paste(diversity_coh$age_cat,
                                    diversity_coh$Cohort,sep="_")

## Statistics
## Wilcoxon test for comparison
#shannon
wilcox_results_babvsb_shannon <- diversity_coh %>%
  group_by(age_cat) %>%
  rstatix::wilcox_test(diversity_shannon ~ age_cat_cohort) %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p.adj) %>%
  dplyr::mutate(sig=ifelse(p.adj<0.05,"significant","non-significant"))

write.csv(wilcox_results_babvsb_shannon, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "wilcox_babelvsbbgut_shannon.csv"))

#observed
wilcox_results_babvsb_observed <- diversity_coh %>%
  group_by(age_cat) %>%
  rstatix::wilcox_test(observed ~ age_cat_cohort) %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p.adj) %>%
  dplyr::mutate(sig=ifelse(p.adj<0.05,"significant","non-significant"))

write.csv(wilcox_results_babvsb_observed, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "wilcox_babelvsbbgut_observed.csv"))


all_plot <- plot_richness_manual(diversity_coh,"Alpha Diversity - cohorts")
all_plot
ggsave(filename = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Figure_2c_alpha_diversity_time_cohorts.pdf"),
       plot = all_plot, height = 6, width = 7, units = "in")


#Cluster proportion per time bin

cluster_infoagebins<-data.frame(sample_data(ps_agecat))

#prepare data for plotting
cl_summary <- cluster_infoagebins %>%
  group_by(Cohort, age_cat, Cluster_Cohort) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  group_by(Cohort, age_cat) %>%
  mutate(Proportion = round(Count / sum(Count),2))

bccolor_mapping <- c(
  "A_BABEL" = "#E4CEC2",
  "B_BABEL" = "#BAAFA5",
  "A_BBGUT" = "#BFE1ED",
  "B_BBGUT" = "#7FAFBF",
  "C_BBGUT" = "#778899"
)
#plot
bc<-ggplot(cl_summary, aes(x = Cohort, y = Proportion, fill = Cluster_Cohort)) +
  geom_bar(stat = "identity", position = "stack") +
  facet_wrap(~ age_cat,nrow=1) +
  labs(
    title = "Proportion of Clusters by Cohort and Age Category",
    x = "Age Category",
    y = "Proportion"
  ) +
  theme_minimal() +
  scale_fill_manual(values = bccolor_mapping)
bc
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_2c_Clusterprop_peragecat_cohort.pdf", sep=""),width=8,height=3)
print(bc)
dev.off()

###########################################################################
## Alpha diversity across months
## Supplementary Figure S9
###########################################################################
get_num_month <- function(days){
  return(ifelse(days<=10,0, ((days%/%31)+1) ))
}

# modify month to be days divided by 31
new_months_merged<- as.numeric(sample_data(ps_combined_gen2)$Age_merged) %>% get_num_month() %>% as.factor()

sample_data(ps_combined_gen2)$Month_merged <- new_months_merged

#function
calc_alpha_manual_month <- function(pseq){
  alpha_man <- microbiome::alpha(pseq,index = c("Observed","Shannon")) %>%
    as.data.frame() %>%
    mutate(Sample_IDs_merged=rownames(.))
  
  meta <- data.frame(sample_data(pseq)) %>%
    select(Sample_IDs_merged, Cohort, Day, Month_merged)
  
  diversity_manual <- left_join(meta, alpha_man, by="Sample_IDs_merged")
  
  diversity_manual$Month_merged <- factor(
    diversity_manual$Month_merged,levels = 0:12)
  
  return(diversity_manual)
}


### Plotting function ###
plot_richness_manual_month <- function(diversity_manual, title){
  
  diversity_manual <- pivot_longer(diversity_manual,
                                   cols = c(observed, diversity_shannon),
                                   names_to = "name")
  
  set.seed(123)
  
  plt <- ggplot(data=diversity_manual,aes(x=Month_merged, y=value)) +
    geom_jitter(aes(color=Cohort),size=2,alpha=0.25,width=0.15) +
    geom_boxplot(aes(fill=Cohort),alpha=0.8,outlier.shape = NA) +
    ggtitle(title) +
    scale_fill_manual(values=c("#8FC3D4" ,"#778899"))+
    scale_color_manual(values=c("#8FC3D4" ,"#778899"))+
    theme_light() +
    theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
  
  return(plt + facet_grid(name ~ ., scales='free'))
}


### Run ###
diversity_month <- calc_alpha_manual_month(ps_combined_gen2)

all_plot_month <- plot_richness_manual_month(
  diversity_month,"Alpha Diversity - cohorts by month")
all_plot_month


## Statistics ##
## Wilcoxon test for comparison between cohorts within each month ##

# Shannon
wilcox_results_month_shannon_month <- diversity_month %>%
  group_by(Month_merged) %>%
  rstatix::wilcox_test(diversity_shannon ~ Cohort) %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p.adj) %>%
  dplyr::mutate(
    sig = ifelse(p.adj < 0.05, "significant", "non-significant") )

#write.csv(wilcox_results_month_shannon_month,"wilcox_babelvsbbgut_shannon_month.csv")

# Observed richness
wilcox_results_month_observed_month <- diversity_month %>%
  group_by(Month_merged) %>%
  rstatix::wilcox_test(observed ~ Cohort) %>%
  dplyr::mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p.adj) %>%
  dplyr::mutate(
    sig = ifelse(p.adj < 0.05, "significant", "non-significant"))

#write.csv(wilcox_results_month_observed_month,"wilcox_babelvsbbgut_observed_month.csv")

ggsave(filename = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Supplementary_Figure_S9_alpha_diversity_timeallmonths_cohorts.pdf"),
       plot = all_plot_month, height = 6, width = 10, units = "in")


############################################################################
## 3. Beta diversity comparison
## Supplementary Table S5
#############################################################################

#merge age (in days) in the combined phyloseq object
# sample_data(ps_combined_gen2)$Age_merged <- as.character(coalesce(sample_data(ps_combined_gen2)$X.days, sample_data(ps_combined_gen2)$Day))
# 
# get_num_month <- function(days){
#   return(ifelse(days<=10,0, ((days%/%31)+1) ))
# }
# 
# # modify month to be days divided by 31
# new_months_merged<- as.numeric(sample_data(ps_combined_gen2)$Age_merged) %>% get_num_month() %>% as.factor()
# 
# sample_data(ps_combined_gen2)$Month_merged <- new_months_merged

#calculate
cohort_ps_abundance <- transform_sample_counts(ps_combined_gen2, function(OTU) OTU / sum(OTU))
distance2 <- "bray"


#Compare using same number of samples per cohort
#same #n of samples from 7 random infants, 1 per month

#extract sample metadata
meta <- data.frame(sample_data(cohort_ps_abundance))

#create Children_IDs_merged by splitting Sample_IDs_merged
meta <- meta %>%
  dplyr::mutate(Children_IDs_merged = str_extract(Sample_IDs_merged, "^[^_]+"))

#randomly select 7 children per cohort
set.seed(1)  # For reproducibility
selected_children <- meta %>%
  group_by(Cohort) %>%
  distinct(Children_IDs_merged) %>%
  slice_sample(n = 7) %>%
  pull(Children_IDs_merged)

#filter metadata to only include selected children
meta_filtered <- meta %>%
  filter(Children_IDs_merged %in% selected_children)

#for each child and each month (0–12), pick 1 random sample
meta_final <- meta_filtered %>%
  filter(Month_merged %in% 0:12) %>%
  group_by(Cohort, Children_IDs_merged, Month_merged) %>%
  slice_sample(n = 1) %>%
  ungroup()
meta_final %>% filter(Cohort=="BBGUT") %>% # View( )
#subset the phyloseq object using selected sample names
samples_to_keep <- meta_final$Sample_IDs_merged
samples_to_keep2<-meta_filtered$Sample_IDs_merged

#filter ps
cohort_ps_filtered <- prune_samples(samples_to_keep, cohort_ps_abundance)
cohort_ps_filtered2<- prune_samples(samples_to_keep2, cohort_ps_abundance)


#test beta diversity differences:

#Bray-Curtis distance
bray_dist2 <- phyloseq::distance(cohort_ps_filtered2, method = "bray")

#metadata
metadata2 <- as(sample_data(cohort_ps_filtered2), "data.frame")

#PERMANOVA
adonis_result2 <- adonis2(bray_dist2 ~ Cohort, data = metadata2)

print(adonis_result2)

# Permutation test for adonis under reduced model
# Terms added sequentially (first to last)
# Permutation: free
# Number of permutations: 999

# adonis2(formula = bray_dist2 ~ Cohort, data = metadata2)
#            Df  SumOfSqs     R2     F      Pr(>F)    
# Cohort     1    2.919   0.05867 15.082  0.001 ***
# Residual  242   46.837  0.94133                  
# Total     243   49.756  1.00000               
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

