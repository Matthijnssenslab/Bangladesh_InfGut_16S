############################################################################
## HEALTHY GUT MICROBIOTA DEVELOPMENT IN BANGLADESHI INFANTS
## ==========================================================
## This script produces:
##   - DMM clustering & maturation stages (GMMs)
##   - Alpha diversity analysis
##   - Composition analysis (top 10 genera per cluster & by month)
##   - COVID-19 lockdown period comparisons
##
## Output figures:
##   Figure 1b, 1c, 1e  |  Figure 4a, 4b
##   Supplementary Figures S6a, S7, S13
##   Supplementary Tables S2, S10
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder 
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(reshape2)
library(scales)
library(stringr)
library(phyloseq)
library(tidyr)
library(data.table)
library(synchrony)
library(rlist)
library(DirichletMultinomial)
library(parallel)
library(microbiome)
library(ggplot2)
library(glue)
library(vegan)
library(tidyverse)

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

############################################################################
## Data Processing ##
############################################################################
# transform bbgut phyloseq object
bbgut_ps_abundance <- transform_sample_counts(bbgut_ps, 
                                              function(OTU) OTU/sum(OTU))

#define infant ids and remove - so that you can match a pattern later
Infants<-unique(sample_data(bbgut_ps)$Child_ID)[1:20] %>% str_replace("-","")
Infants2<-unique(sample_data(bbgut_ps)$Child_ID)[1:20]

#define healthy samples
HTP_samples<-sample_data(bbgut_ps)$Sample_sequencing_ID[sample_data(bbgut_ps)$sample_category == "HTP"] %>%na.omit()

# agglomerate by genus                
genus_agglom <- tax_glom(bbgut_ps, taxrank="Genus")
genus_agglom


################
#### 1. DMM ####
################

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

list_seeds <- list()
set.seed(1)
for (seed in 1:20){
  print(paste('seed: ',seeds[seed]))
  fit <- run_dmn(10,genus_agglom,seeds[seed])
  list_seeds[[seed]]  <-  fit
}
DMM_out_all <- list.rbind(list_seeds)

# check per seed where the minimum is reached
DMM_out_min_BIC <- DMM_out_all %>% 
  filter(score=='BIC') %>% 
  group_by(seed) %>% 
  slice_min(value) %>%
  as.data.frame()

# get best seed and seed plots
most_freq_numClust <- (sort(table(DMM_out_min_BIC$numCluster),decreasing = TRUE) %>% names(.) %>%
                         as.numeric())[1]

best_seed <- DMM_out_min_BIC %>% filter(value == min(value)) %>% pull(seed)

plot_min_BIC_clusters <- ggplot(DMM_out_min_BIC, aes(x= numCluster))+
  geom_histogram()
plot_min_BIC_clusters
plot_best_seed_clusters <- ggplot(DMM_out_all %>% 
                                    filter(score=='BIC' & seed==best_seed), 
                                  aes(x=numCluster, y=value)) +
  geom_point()+
  geom_path()+
  scale_x_continuous(breaks = c(1:10))+
  xlab("Number of Dirichlet Components") + 
  ylab("Model fit (BIC)")
plot_best_seed_clusters
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "SPLM_BIC_best_clusternumber.pdf", sep=""),width=7,height=6)
print(plot_best_seed_clusters)
dev.off()


#### Best fit stats ####
fit_optimal <- dmn(otu_table(genus_agglom), most_freq_numClust , 
                   verbose = F , seed = best_seed) 
sample_data_gen_agglom <- sample_data(genus_agglom) %>% data.frame()
sample_data_gen_agglom$ClusterNum <- as.character(
  mixture(fit_optimal, assign=TRUE))
sample_data_gen_agglom$Cluster <- chartr("123","ABC",
                                         as.character(sample_data_gen_agglom$ClusterNum))
sample_data(genus_agglom)<-sample_data_gen_agglom



# statistics per cluster
per_clust_stats <- mixturewt(fit_optimal)
rownames(per_clust_stats) <- c(LETTERS[1:most_freq_numClust])
num_per_clust <- table(sample_data_gen_agglom$Cluster)
per_clust_stats['A','cluster'] <- num_per_clust['A']
per_clust_stats['B','cluster'] <- num_per_clust['B']
per_clust_stats['C','cluster'] <- num_per_clust['C']

write.csv(per_clust_stats, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Per_cluster_stats.csv"))

# per sample statistics
per_sample_stats <- as.data.frame(
  cbind(mixture(fit_optimal), 
        mixture(fit_optimal, assign=TRUE)))
colnames(per_sample_stats) <- c(LETTERS[1:most_freq_numClust],'selection')
per_sample_stats$Cluster <- chartr("123","ABC",
                                   as.character(per_sample_stats$selection))
per_sample_stats <- per_sample_stats %>% mutate(selection.prob=pmax(A,B,C))
per_sample_stats %>% group_by(Cluster) %>% tally()

write.csv(per_sample_stats, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Per_sample_stats.csv"))


############################
##### Cluster Changes ######
############################
clusters <- sample_data_gen_agglom %>% 
  select(Child_ID,Sample_sequencing_ID,Day,Cluster) %>%
  arrange(Child_ID,Day) 

##### Helper functions #####

#define maturation shifts:
# Year 1
# Change if cluster changes in the next sample
# Setback if cluster sample2<cluster sample1
# Significant change: if setback after 2 samples with the same cluster
# 
# Year 2
# Significant changes: all setbacks (due to the big gap between samples)

sliding_window_concat <- function(input_vector) {
  n <- length(input_vector)
  result <- character(n)
  for (i in 1:n) {
    if (i < 2) {
      result[i] <- NA
    } else if(i==2){
      result[i] <- paste(input_vector[(i - 1):i], collapse = "")
    }
    else {
      result[i] <- paste(input_vector[(i - 2):i], collapse = "")
    }
  }
  return(result)
}

isChange <- function(input_XXX){
  if (is.na(input_XXX)){
    return("No")
  }
  last_two <- unlist(str_split(input_XXX,"")) %>% tail(2)
  current <- last_two[2]
  previous <- last_two[1]
  change <- ifelse(current == previous,"No","Yes")
  return(change)
}

isStepback <- function(input_XXX){
  if (is.na(input_XXX)){
    return("No")
  }
  last_two <- unlist(str_split(input_XXX,"")) %>% tail(2)
  current <- last_two[2]
  previous <- last_two[1]
  # no change so no stepback
  if(current < previous){
    return ('Yes')
  } else {
    return ('No')
  }
}

isImportantChange <- function(input_XXX,day){
  if (is.na(input_XXX)){
    return("No")
  }
  split_clusters <- unlist(str_split(input_XXX,""))
  if (length(split_clusters) == 3){
    current <- split_clusters[3]
    previous_1 <- split_clusters[2]
    previous_2 <- split_clusters[1]
    if (current != previous_1 & 
        previous_1 == previous_2 &
        day < 366) {
      return ("Yes")
    }
    # year 2 criteria
    else if (current != previous_1 & 
             day >= 366) {
      return ("Yes")
    }
  }
  return ('No')
}

addChanges <- function(Child_ID_df) {
  Child_ID_df$cluster_hist <- sliding_window_concat(Child_ID_df$Cluster)
  Child_ID_df$change <- sapply(Child_ID_df$cluster_hist, isChange)
  Child_ID_df$stepback <- sapply(Child_ID_df$cluster_hist, isStepback)
  Child_ID_df$importance <- mapply(isImportantChange, 
                                   Child_ID_df$cluster_hist, 
                                   Child_ID_df$Day)
  return(Child_ID_df)
}


##### Add change cols to cluster ####
clusters_change <- list.rbind(
  lapply(split(clusters, clusters$Child_ID),addChanges)) 
rownames(clusters_change) <- clusters_change$Sample_sequencing_ID

write.csv(clusters_change, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "Maturation_changes_persample_clusters.csv"))

#################################################################################
### 2. Maturation Development timeline (MATURATION STAGES - healthy samples only) ##
### Publication *Figure 1e* ###
#################################################################################

test_healthynew<- clusters_change %>% 
  filter(Sample_sequencing_ID %in% HTP_samples) %>%
  drop_na(., Child_ID)

test_healthynew<-sample_data(test_healthynew)

#Plot 
cluster_stages_plot <- ggplot(test_healthynew, 
                              aes(x=Day, y=Cluster, color=Cluster, group=1)) + 
  geom_point(size = 4) +
  scale_color_manual(values= as.character(pastel_colors_2))  +
  #facet_wrap(~Child_ID, nrow = 2, scales = "free_y",) +  # Use nrow = 2 for two rows
  facet_wrap(~Child_ID, nrow=4,ncol=5)+
  geom_line(color = 'grey') +
  geom_vline(xintercept=365, linetype="dotted")+
  theme(axis.text=element_text(size=18))+
  theme_light()+
  ggtitle(glue("Healthy samples, clusters = {most_freq_numClust}, seed = {best_seed}")) +
  scale_y_discrete(limits = c("A", "B", "C"),
                   expand = c(0.3, 0)) +
  xlab('Days after birth')+
  ylab('GM Maturation Stage')
cluster_stages_plot
# Print the plot
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_1e_Maturation_stages_BBGUT_dmm.pdf", sep=""),width=6,height=8)
print(cluster_stages_plot)
dev.off()


# Statistics to confirm temporal ranking of the clusters is significant
cluster_appearance <- test_healthynew %>%
  data.frame(.) %>%
  group_by(Child_ID, Cluster) %>%
  dplyr::summarise(appearance_day = min(Day)) %>% 
  pivot_wider(id_cols = Child_ID,names_from = Cluster,values_from = appearance_day) %>%
  column_to_rownames("Child_ID")

KENDALOUTCOME <- kendall.w(t(cluster_appearance), 
                           nrands = 10000, type = 1, quiet = T)
KENDALOUTCOME_table <- as.data.frame(matrix(nrow = 4, ncol=2))
KENDALOUTCOME_table[,1] <- c("Kendall's W (uncorrected for ties)",
                             "Kendall's W (corrected for ties)",
                             "Spearman's ranked correlation",
                             "Kendall's W p-value (one-tailed test [greater])")
KENDALOUTCOME_table[,2] <- c(KENDALOUTCOME$w.uncorrected, KENDALOUTCOME$w.corrected, 
                             KENDALOUTCOME$spearman.corr, KENDALOUTCOME$pval.rand)
colnames(KENDALOUTCOME_table) <- c('Kendall test, all infants','value')


write.csv(KENDALOUTCOME_table, file = paste(workdir,format(Sys.time(), "%Y-%m-%d"),
                                            "maturationstages_kendalltest.csv", sep=""))

#####################################################
## Transition day from AtoB and BtoC
## Publication *Supplementary Figure S6a*
#####################################################

# Find the first day of transition from A to B and B to C
transition_data <- test_healthynew %>%
  data.frame(.) %>%
  filter(importance=="Yes" & stepback=="No") %>%
  group_by(Child_ID, Cluster) %>%
  dplyr::summarise(TransitionDay = min(Day))

# Create a new dataframe for plotting
plot_data <- transition_data %>%
  mutate(TransitionType = ifelse(Cluster == "B", "AtoB", "BtoC"))

# Create a boxplot
t_time<-ggplot(plot_data, aes(x = TransitionType, y = TransitionDay, fill = TransitionType)) +
  geom_boxplot() +
  labs(title = "Average Day of Transition from A to B and B to C",
       x = "Transition Type",
       y = "Day") +
  geom_point(aes(color=Child_ID), size=3)+
  scale_fill_manual(values = c("#F2E8DF","#AFC1D5")) +
  scale_color_manual(values=pastel_colors)+
  scale_y_continuous(breaks=pretty_breaks())+
  theme_light()
t_time
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_S6a_Cluster_transition_time_BBGUT.pdf", sep=""),
    width=7,height=7)
print(t_time)
dev.off()

############################################################################
## 3. Alpha diversity measures (observed richness and Shannon diversity) 
## of the samples within every stage
## Publication *Figure 1d*
#############################################################################
genus_agglom_infonly<-subset_samples(genus_agglom, fecalsample_type!="maternal")
Alpha_percluster <- plot_richness(genus_agglom_infonly, x="Cluster", color="Cluster" ,nrow =1, 
                                  measures = c("Observed", "Shannon"))+
  geom_boxplot()+
  scale_color_manual(values= pastel_colors_2)  +
  ggtitle("Alpha Diversity") + 
  theme_light()
Alpha_percluster


# plot mothers as separate cluster C_m

bbgut_meta_mothercluster<-data.frame(sample_data(genus_agglom))
bbgut_meta_mothercluster=bbgut_meta_mothercluster %>% mutate(Cluster = ifelse(fecalsample_type == "maternal", "C_m", Cluster))
bbgut_meta_mothercluster <-sample_data(bbgut_meta_mothercluster)
bbgut_ps3<-phyloseq(bbgut_meta_mothercluster,otu_table(bbgut_ps),tax_table(bbgut_ps))

#Agglomerate table at the Genus level
genus_agglom3  <- tax_glom(bbgut_ps3, taxrank="Genus")
genus_agglom3

#plot alpha diversity per cluster (including mothers as a separate cluster)
Alpha_percluster_M <- plot_richness(genus_agglom3, x="Cluster", color="Cluster" ,nrow =1, measures = c("Observed", "Shannon"))+
  geom_boxplot()+
  scale_color_manual(values= as.character(pastel_colors_4))  +
  ggtitle("Alpha Diversity") + 
  theme_light()
Alpha_percluster_M
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_1d_Alpha_diversity_clusters.pdf", sep=""),width=7,height=5)
print(Alpha_percluster_M)
dev.off()

#NOTE:
#Statistics for alpha diversity in *Alpha_gmm_lmm.R*
#calculate alpha diversity

calc_alpha_manual<- function(pseq,cluster_lvls){
  alpha_man <- microbiome::alpha(pseq,index = c("Observed","Shannon")) %>% 
    as.data.frame() %>% mutate(Sample_sequencing_ID=rownames(.))
  
  meta <- data.frame(sample_data(pseq)) %>% select(Sample_sequencing_ID,Day,Cluster)
  
  diversity_manual <- left_join(meta,alpha_man,by="Sample_sequencing_ID")
  diversity_manual$Cluster <- factor(diversity_manual$Cluster, levels=cluster_lvls)
  
  return(diversity_manual)
}

#infants only
alpha_infonly<-calc_alpha_manual(genus_agglom_infonly,c("A","B","C"))
#full dataset
alpha_fulldataset<-calc_alpha_manual(genus_agglom3,c("A","B","C","C_m"))


##########################################
## Sample time distribution (per cluster)
## Publication *Figure 1b*
##########################################

#see how many samples you have per cluster (c infant samples = -18 mothers)
clusters %>% group_by(Cluster) %>% tally()

#density plot, showing the probability of distribution of the sample being around 
#a specific age

t_distrib<-clusters %>%
  filter(!is.na(Child_ID)) %>%
  ggplot(aes(x=Day, color=Cluster, fill=Cluster))+
  #  geom_histogram(aes(alpha=0.5,position="identity"),bins = 104)+
  geom_density(alpha=0.5)+
  scale_fill_manual(values=pastel_colors_2)+
  scale_color_manual(values=pastel_colors_2)+
  geom_vline(xintercept=365,linetype="dotted")+
  scale_x_continuous(breaks=pretty_breaks())+
  labs(title = "Distribution of sample age by maturation stage",
       x = "Age (days)",
       y = "Density") +
  theme_light()
t_distrib
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_1b_Clusters_timedistribution.pdf", sep=""),width=7,height=5)
print(t_distrib)
dev.off()


############################################################
## Cluster composition (mean relative abundance per cluster)
## top 10 most abundant genera per cluster
## Publication *Figure 1c*
############################################################
genus_agglom3_abund<-transform_sample_counts(genus_agglom3, function(OTU) OTU/sum(OTU))
genus_agglom_abundance_comp <- psmelt(genus_agglom3_abund)
#calculate mean abundance per genus per cluster
mean_abund_perclust<-genus_agglom_abundance_comp %>% 
  group_by(Genus,Cluster) %>%
  dplyr::summarise(mean_abundance=mean(Abundance)) %>%
  pivot_wider(id_cols=Genus, names_from=Cluster, values_from=mean_abundance) # pivot to get mean
#abundance per genus per cluster in columns

#get top10 per cluster
top10_A<-mean_abund_perclust %>% arrange(desc(A)) %>% head(10) %>% pull(Genus)
top10_B<-mean_abund_perclust %>% arrange(desc(B)) %>% head(10) %>% pull(Genus)
top10_C<-mean_abund_perclust %>% arrange(desc(C)) %>% head(10) %>% pull(Genus)
top10_M<-mean_abund_perclust %>% arrange(desc(C_m)) %>% head(10) %>% pull(Genus)

#Take a look at % of total abundance per cluster for the top10 genera
mean_abund_perclust %>% arrange(desc(C_m)) %>% head(10) %>% pull(C_m) %>% sum()

#genera to keep, the union between all top10s
gen_tokeep_top10percl<-purrr::reduce(list(top10_A,top10_B,top10_C,top10_M),union)

#turn all the rest to others
#calculate sum of Other genera per cluster
other_sum <- mean_abund_perclust %>% 
  mutate(Genus = ifelse(!Genus %in% gen_tokeep_top10percl, "Other", Genus)) %>%
  filter(Genus=="Other") %>% as.data.frame() %>%
  select(-Genus) %>% colSums()

#subset genus abundance df with only genera to keep (from union)
top10_mean_abund_perclust <- mean_abund_perclust %>% 
  filter(Genus %in% gen_tokeep_top10percl) %>%
  column_to_rownames("Genus")

#add "Others" to subset df above
top10_mean_abund_perclust["Other",] <- other_sum

#change rownames to colnames again
top10_mean_abund_perclust<-top10_mean_abund_perclust %>% rownames_to_column("Genus")

#define specific colors per genus
genus_colors_cl <- c(
  "Bacteroides"      ="#BAAFA5",    
  "Bifidobacterium"  = "#6B96B6",    
  "Blautia"           = "#B9DCE4", 
  "Catenibacterium"   = "#E0F1F7", 
  "Dorea"             ="#E3E1DA",
  "Enterococcus"      = "#FFD9AC",  
  "Escherichia/Shigella" = "#8EB9B7",
  "Faecalibacterium"     = "#AEC89E",
  "Holdemanella"        = "#D1CF99", 
  "Enterobacter"          = "#C694B3",
  "Lactobacillus"        = "#F8BDAF",
  "Ligilactobacillus"    = "#F5DDB7",
  "Limosilactobacillus"  = "#75AEBD",
  "Megasphaera"          ="#A6CFCF",
  "Parolsenella"        ="#FFC1BE",
  "Segatella"           = "#E78F8D",
  "Staphylococcus"       = "#FED4E3",
  "Streptococcus"       = "#D3A9B9",
  "Veillonella"         =  "#E4CEC2",
  "Other"                = "lightgrey")

top10_mean_abund_perclust$Genus <- factor(top10_mean_abund_perclust$Genus, 
                              levels = c("Bacteroides", "Bifidobacterium", 
"Blautia", "Catenibacterium", "Dorea", "Enterococcus", "Escherichia/Shigella", 
"Faecalibacterium", "Holdemanella","Enterobacter", "Lactobacillus", 
"Ligilactobacillus", "Limosilactobacillus","Megasphaera", "Parolsenella", "Segatella", 
"Staphylococcus", "Streptococcus","Veillonella","Other"))

#plot
top10_comp<-top10_mean_abund_perclust %>%
  pivot_longer(cols=c(A,B,C,C_m),names_to="Cluster",values_to="mean_abundance") %>% #pivot longer for plotting
  ggplot(aes(x=Cluster, y=mean_abundance, fill=Genus)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_fill_manual(values= genus_colors_cl)+
  labs(title = "Rel. abundance of bacterial genera by cluster",
       x = "Cluster",
       y = "Mean abundance") +
  theme_light()
top10_comp                  
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_1c_Composition_clusters_top10genera.pdf", sep=""),width=6,height=7)
print(top10_comp)
dev.off()


############################################################
## 4. Relative abundance per month in healthy samples 
## top 10 most abundant genera across clusters
## Publication *Figure 4a*
############################################################
#modify months function (first 10 days of life = 0,1:24)
get_num_month <- function(days){
  return(ifelse(days<=10,0, ((days%/%31)+1) ))
}

# modify month
genus_agglom_abundance_comp_3<-genus_agglom_abundance_comp[, c('Abundance','Sample_sequencing_ID','Genus','Cluster','Day','fecalsample_type')]
genus_agglom_abundance_comp_month <- genus_agglom_abundance_comp_3 %>%
  filter(fecalsample_type=="infant") %>%
  mutate(Month=get_num_month(Day))


#make df with mean abundance of a genus per month
mean_abund_perclust_permonth<-genus_agglom_abundance_comp_month %>%
  filter(Sample_sequencing_ID %in% HTP_samples) %>%
  group_by(Genus,Month) %>%
  dplyr::summarise(mean_abundance=mean(Abundance)) %>%
  pivot_wider(id_cols=Genus, names_from=Month, values_from=mean_abundance)

#count samples per month
genus_agglom_abundance_comp_month %>% 
  filter(Sample_sequencing_ID %in% HTP_samples) %>%
  group_by(Month) %>% 
  dplyr::summarize(nsamplesperclust=n_distinct(Sample_sequencing_ID)) %>% # View( )


#calculate sum of Other genera per cluster
other_sum_m <- mean_abund_perclust_permonth %>% 
  mutate(Genus = ifelse(!Genus %in% gen_tokeep_top10percl, "Other", Genus)) %>%
  filter(Genus=="Other") %>% as.data.frame() %>%
  select(-Genus) %>% colSums()

#subset genus abundance df with only genera to keep (from union)
top10_mean_abund_per_m<- mean_abund_perclust_permonth %>% 
  filter(Genus %in% gen_tokeep_top10percl) %>%
  column_to_rownames("Genus")

#add "Others" to subset df above
top10_mean_abund_per_m["Other",] <- other_sum_m

top10_mean_abund_per_m<-top10_mean_abund_per_m %>% 
  rownames_to_column("Genus")

top10_mean_abund_per_m$Genus <- factor(top10_mean_abund_per_m$Genus, 
                         levels = c("Bacteroides", "Bifidobacterium", 
"Blautia", "Catenibacterium", "Dorea", "Enterococcus", "Escherichia/Shigella", 
"Faecalibacterium", "Holdemanella","Enterobacter", "Lactobacillus", 
"Ligilactobacillus", "Limosilactobacillus","Megasphaera", "Parolsenella", "Segatella", 
"Staphylococcus", "Streptococcus","Veillonella","Other"))


#plot all months
comp_months_healthy<- top10_mean_abund_per_m %>% 
  pivot_longer(cols=c(-Genus),names_to="Month",values_to="mean_abundance") %>% 
  mutate(Month = factor(Month, levels = 0:24)) %>%
  ggplot(aes(x=Month, y=mean_abundance, fill=Genus))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=genus_colors_cl)+
  labs(title = "Rel. abundance of bacterial genera by month",
       x = "Time (months)",
       y = "Mean Relative Abundance") +
  theme_light()
comp_months_healthy
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_4a_Composition_bymonth_HEALTHY.pdf", sep=""),width=12,height=6)
print(comp_months_healthy)
dev.off()

############################################################
### 4. Maternal sample composition 
### Relative abundance of the top 10 most common genera found 
### in maternal samples
### Publication *Supplementary Figure S7*
############################################################

genus_agglom_abundance_comp_M <- genus_agglom_abundance_comp_3 %>%
  filter(fecalsample_type=="maternal") 

#make df with mean abundance of a genus per month
abund_permother<-genus_agglom_abundance_comp_M %>%
  group_by(Genus, Sample_sequencing_ID) %>%
  #dplyr::summarise(mean_abundance=mean(Abundance))%>%
  pivot_wider(id_cols=Genus,names_from=Sample_sequencing_ID, values_from=Abundance)

# get most abundant genera in maternal sample group
top15_M<-genus_agglom_abundance_comp %>% 
  group_by(Genus,Cluster) %>%
  filter(Cluster=="C_m") %>%
  dplyr::summarise(t_abundance=sum(Abundance)) %>%
  #dplyr::summarise(mean_abundance=mean(Abundance)) %>%
  pivot_wider(id_cols=Genus, names_from=Cluster, values_from=t_abundance) %>%
  arrange(desc(C_m)) %>% 
  head(15) %>% 
  pull(Genus)


#calculate sum of Other genera per cluster
other_sum_M <- abund_permother %>% 
  mutate(Genus = ifelse(!Genus %in% top15_M, "Other", Genus)) %>%
  filter(Genus=="Other") %>% as.data.frame() %>%
  select(-Genus) %>% colSums()

#subset genus abundance df with only genera to keep (from union)
top15_abund_per_mother<- abund_permother %>% 
  filter(Genus %in% top15_M) %>%
  column_to_rownames("Genus")

#add "Others" to subset df above
top15_abund_per_mother["Other",] <- other_sum_M

top15_abund_per_mother<-top15_abund_per_mother %>% 
  rownames_to_column("Genus")

# Manually assign each color to a genus
genus_colors_m <- c(
  "Bifidobacterium"      = "#6B96B6",
  "Blautia"              = "#B9DCE4",
  "Catenibacterium"      = "#E0F1F7",
  "Dialister"            = "#FBD8CB",
  "Enterococcus"         = "#FFD9AC",
  "Escherichia/Shigella" = "#8EB9B7",
  "Faecalibacterium"     = "#AEC89E",
  "Holdemanella"         = "#D1CF99",
  "Ligilactobacillus"    = "#F5DDB7",
  "Hallella"             = "#D8C8E8", 
  "Leyella"              = "#C0E0E0", 
  "Segatella"            = "#E78F8D",
  "Prevotellamassilia"   = "#F1B2B1",
  "Ruminococcoides"      = "#C8BFB5",
  "Streptococcus"        = "#D3A9B9",
  "Other"                = "lightgrey"
)


top15_abund_per_mother$Genus <- factor(top15_abund_per_mother$Genus, 
                                       levels = c("Bifidobacterium","Blautia","Catenibacterium","Dialister","Enterococcus","Escherichia/Shigella",
"Faecalibacterium","Holdemanella","Ligilactobacillus","Hallella", "Leyella",
"Segatella","Prevotellamassilia","Ruminococcoides","Streptococcus","Other"))

#plot maternal abundances
m_comp<-top15_abund_per_mother %>% 
  pivot_longer(cols=c(-Genus),names_to="Sample_sequencing_ID",values_to="mean_abundance") %>%
  ggplot(aes(x=Sample_sequencing_ID, y=mean_abundance, fill=Genus))+
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values=genus_colors_m)+
  labs(title = "Rel. abundance of bacterial genera in maternal samples",
       x = "Maternal samples",
       y = "Relative Abundance") +
  theme_light()
m_comp
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_S7_Composition_maternal_samples_top15.pdf", sep=""),width=14,height=6)
print(m_comp)
dev.off()

#--- How many mothers have Bacteroides? ---#

gen_prevalence <- genus_agglom_abundance_comp_M %>%
  group_by(Genus, Sample_sequencing_ID) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  group_by(Genus) %>%
  dplyr::summarise(
    n_individuals_present = sum(Abundance >= 0.005),
    .groups = "drop"
  )


#--- Check for Prevotellaceae members in mothers ---#

target_genera <- c("Segatella", "Leyella", "Prevotellamassilia")

resultprev <- abund_permother %>%
  filter(Genus %in% target_genera) %>%
  pivot_longer(
    -Genus,
    names_to = "Sample_sequencing_ID",
    values_to = "Abundance"
  ) %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  mutate(Abundance_percentage=round(Abundance,2)*100)

#total abundance of the 3 genera per individual
total_abundanceprev <- resultprev %>%
  group_by(Sample_sequencing_ID) %>%
  summarise(
    total_abundance = sum(Abundance),
    .groups = "drop"
  )
total_abundanceprev<- total_abundanceprev %>% 
  dplyr::mutate(total_abundance_percent=round(total_abundance,2)*100)

#average of abundance across individuals
total_abundance_positive <- total_abundanceprev %>%
  filter(total_abundance>=0.005)
mean_total_abundance <- mean(total_abundance_positive$total_abundance)

#genera present per individual
prevgenera_present <- resultprev %>%
  filter(Abundance >= 0.005) %>%
  group_by(Sample_sequencing_ID) %>%
  dplyr::summarise(
    genera_present = paste(Genus, collapse = ", "),
    n_genera = n(),
    .groups = "drop"
  )

############################################################
### 5. Comparison of Rel. Abundances in healthy samples 
### BEFORE and AFTER the COVID-19 lockdown period
### Publication *Figure 4b*
############################################################
#covid samples
lockdown_samples<-read.csv("Data/lockdown_samples_y2.csv") %>% pull(x)
#sum(lockdown_samples %in% HTP_samples)
prelockdown_samples<-read.csv("Data/prelockdown_samples_y2.csv") %>% pull(x)
#sum(prelockdown_samples %in% HTP_samples)

#group in two categories
median_pergen_month<- genus_agglom_abundance_comp %>%
  filter(Day>365) %>%
  #keep 16 infants with data for both y1 and y2
  filter(Child_ID %in% c("I01","I02","I05","I06","I07","I08","I09","I10",
                  "I11","I12","I13","I14","I15","I17","I18","I20")) %>%
  select(Abundance,Sample_sequencing_ID,Child_ID, Month, Genus) %>%
  filter(Sample_sequencing_ID %in% HTP_samples) %>%
  filter(Genus %in% gen_tokeep_top10percl) %>% 
  dplyr::mutate(Month_grouped = case_when(
    Sample_sequencing_ID %in% prelockdown_samples ~ "months13-16",
    Sample_sequencing_ID %in% lockdown_samples ~ "months17-24",
    TRUE ~ NA_character_)) %>% 
  group_by(Month_grouped,Genus,Child_ID) %>%
  dplyr::summarise(median_abundance=median(Abundance)) %>%
  mutate(Month_grouped = factor(Month_grouped, levels = c("months13-16", "months17-24")))


#load 
library(broom)

#empty list to store results
results_list <- list()

#list of unique genera
unique_genera <- unique(median_pergen_month$Genus)

for (genus in unique_genera) {
  
  #Subset data for the current genus
  subset_data <- median_pergen_month %>% filter(Genus == genus)
  
  #Check if we have enough data to perform the test
  if (n_distinct(subset_data$Child_ID) > 1) {  # Need at least 2 pairs
    
    #Reshape data so we can compute paired differences
    paired_wide <- subset_data %>%
      select(Genus, Child_ID, median_abundance, Month_grouped) %>%
      pivot_wider(names_from = Month_grouped, values_from = median_abundance)
    
    #Perform paired Wilcoxon signed-rank test
    wilcox_test <- wilcox.test(paired_wide$`months13-16`, 
                               paired_wide$`months17-24`, paired = TRUE)
    test_results2 <- data.frame(
      Genus = genus,
      Wilcoxon_p_value = wilcox_test$p.value
    )
    
    results_list[[genus]] <- test_results2
  }
}

#Combine all results into a single dataframe
final_results_df <- bind_rows(results_list)

#Adjust p-values using Benjamini-Hochberg correction
final_results_df$Adjusted_p_value <- p.adjust(final_results_df$Wilcoxon_p_value, method = "BH")

#Add significance label for readablity
final_results_df$Significance <- ifelse(final_results_df$Adjusted_p_value < 0.05, "significant", "non-significant")
#write.csv(final_results_df,"2yagegroup_comparison_generaabundances.csv")

#Filter only significant results
significant_results <- final_results_df %>% filter(Adjusted_p_value < 0.05)

#Which genera show significant differences?
gen_forplot <- unique(significant_results$Genus)


#plot
genc<-median_pergen_month %>%
  filter(Genus %in% gen_forplot) %>%
  ggplot(aes(x = Month_grouped, y = median_abundance, fill = Month_grouped)) +
  geom_boxplot() +
  scale_fill_manual(values = c("months13-16" = "#F8BDAF", 
                               "months17-24" = "#E78F8D")) + 
  labs(
    x = "Month Grouped",
    y = "Median Abundance",
    fill = "Month Grouped",
  ) +
  facet_wrap(~ Genus, scales = "free_y",ncol=1,nrow=7) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold.italic")
  )
genc
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_4b_Genera_changes_2ndyear_agegroups_final.pdf", sep=""),width=3,height=10)
print(genc)
dev.off()


################################################################
## 5. Normalized 2-month bacterial bundance changes across periods
## E-E, E-L, L-L
## Publication *Supplementary Figure 13*
################################################################

# PARAMETERS
months_keep  <- 9:24  # month window
lag_months   <- 2     # equidistant gap
random_seed  <- 1     # for reproducibility


#Prepare the base dataframe with period labels
df <- genus_agglom_abundance_comp %>%
  filter(
    Sample_sequencing_ID %in% HTP_samples,
    #keep 16 infants with data for both y1 and y2
    Child_ID %in% c("I01","I02","I05","I06","I07","I08","I09","I10",
                    "I11","I12","I13","I14","I15","I17","I18","I20"),
    Genus %in% gen_tokeep_top10percl,
    Month %in% months_keep) %>%
  #define pre-lockdown and lockdown periods
  dplyr::mutate(
    Month_grouped = case_when(
      Sample_sequencing_ID %in% lockdown_samples ~ "months17-24",
      TRUE ~ "months13-16")) %>%
  group_by(Month, Child_ID, Genus) %>%
  dplyr::mutate(Abundance_median=median(Abundance)) %>%
  ungroup() %>%
  select(Child_ID, Genus, Month, Abundance_median, Month_grouped)
#openxlsx::write.xlsx(df,"median_abundance_lockdown_prelockdown.xlsx")

#Generate all valid 2‑month pairs
all_pairs <- inner_join(
  df, df,
  by     = c("Child_ID","Genus"),
  suffix = c("_1","_2")) %>% 
  filter(Month_2 == Month_1 + lag_months) %>%
  #define intervals
  dplyr::mutate(
    Period_Comparison = case_when(
      Month_grouped_1 == "months13-16" & Month_grouped_2 == "months13-16" ~ "E-E",
      Month_grouped_1 == "months13-16" & Month_grouped_2 == "months17-24"    ~ "E-L",
      Month_grouped_1 == "months17-24" & Month_grouped_2 == "months17-24"    ~ "L-L",
      TRUE ~ NA_character_)
  ) %>%
  filter(!is.na(Period_Comparison))
#openxlsx::write.xlsx(all_pairs,"differences_abundance_allpairs_2months.xlsx")

#Normalized difference
set.seed(random_seed)

nall_pairs <- all_pairs %>%
  group_by(Child_ID, Genus, Period_Comparison) %>%
  ungroup() %>%
  mutate(
    NormDiff = (Abundance_median_2 - Abundance_median_1) / ((Abundance_median_1 + Abundance_median_2)/2))


#load library

# Set factor levels to control order
nall_pairs$Period_Comparison <- factor(nall_pairs$Period_Comparison,
                                         levels = c("E-E", "E-L", "L-L"))

# Define custom colors
colors <- c("E-E" = "#F8BDAF",
            "E-L" = "#FFD9AC", 
            "L-L" = "#E78F8D")


#boxplot
boxch<-ggplot(nall_pairs, aes(x = Period_Comparison, y = NormDiff, fill = Period_Comparison)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Genus, scales = "free_y") +
  theme_minimal() +
  scale_fill_manual(values = colors) +
  labs(
    title = sprintf("Normalized %d-Month Abundance Changes", lag_months),
    x = NULL,
    y = "Δabundance / mean(abundances)"
  )
boxch
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_13_2month_abundance_changes.pdf", sep=""),width=10,height=6)
print(boxch)
dev.off()

#NOTE:
#Statistical differences between groups can be calculated 
#in *covidmonthpairs_comparison_lmm.R*
