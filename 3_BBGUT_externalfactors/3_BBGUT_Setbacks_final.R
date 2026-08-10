############################################################################
## MATURATION SETBACKS — BBGUT COHORT
## ====================================
## - Identifies and visualizes maturation setbacks (cluster changes)
## - Calculates maturation scores per sample
## - Tests setback association with disease events
## - Compares maturation scores between 2nd-year age groups
##
## Output figures & tables:
##   Figure 3c, 3d, 4c  |  Supplementary Figures S11a, S11b  |  Supplementary Table S8
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(reshape2)
library(stringr)
library(phyloseq)
library(tidyr)
library(data.table)
library(rlist)
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


#read in cluster changes file that was created for figure 1
#this file includes all the cluster information per cluster including the
# cluster changes (setbacks) for all infant samples
cluster_changes<-read.csv("Data/Maturation_changes_persample_clusters.csv")
cluster_changesinf <- cluster_changes %>%
  filter(str_starts(Child_ID, "I"))

#####################################################
## 1. Visualizing GMM changes (important setbacks)
## Part of Publication *Figure 3c*
#####################################################

#plot to visualize all changes (part of Publication *Figure 3c* example)
setback_plot <- ggplot(cluster_changesinf, aes(x=Day, y=Cluster, color=Cluster, group=1)) + 
  geom_point(data = cluster_changesinf[cluster_changesinf$setback == "Yes" & cluster_changesinf$importance== "Yes",], 
             aes(x=Day,y=Cluster, color=Cluster, group=Child_ID),size = 4)+# 
  geom_point(size = 1)+
  scale_color_manual(values= as.character(pastel_colors_2))  +
  facet_wrap(~Child_ID, nrow = 2)+
  geom_vline(xintercept=365, linetype="dotted")+
  theme_light() +
  scale_y_discrete(expand=c(0.3, 0))+
  xlab('Days after birth')
setback_plot
# pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
#           "allsetbacks_BBGUT_dmm.pdf", sep=""),width=16,height=6)
# print(setback_plot)
# dev.off()


#####################################################
## 2. Maturation score calculation
#####################################################

#Data processing
bbgut_ps_inf <- subset_samples(bbgut_ps, fecalsample_type == "infant") 


#define healthy samples
HTP_samples<-sample_data(bbgut_ps)$Sample_sequencing_ID[sample_data(bbgut_ps)$sample_category == "HTP"] %>%na.omit()


htp_filtered_bbgut_ps  <- subset_samples(bbgut_ps_inf, 
                                         Sample_sequencing_ID %in% 
                                           HTP_samples)

htp_filtered_bbgut_ps_melt <- psmelt(htp_filtered_bbgut_ps)

infants<-htp_filtered_bbgut_ps_melt$Child_ID %>% unique()


## Steps for calculating maturation score 

# TOP 15 GENERA BBGUT
tot_ab_htp_filtered <- htp_filtered_bbgut_ps_melt %>% 
  group_by(Genus) %>% 
  # get total abundance per genus
  dplyr::summarize(total_abundance=sum(Abundance)) %>% 
  # sort descending
  arrange(desc(total_abundance)) %>% 
  # add new column with relative ab
  mutate(relative_ab=round(total_abundance/sum(total_abundance)*100,4))

tot_ab_htp_filtered %>% # View( )
Genera_tokeep <- head(tot_ab_htp_filtered$Genus,15)


# per sample relative ab of top 15 genera

# phyloseq object for top 15 relative ab genera
bbgut_ps_relative_top15 <- htp_filtered_bbgut_ps %>% 
  # 1. agglomerate healthy_ps into Genus
  tax_glom(taxrank="Genus") %>% 
  # 2. calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU)) %>%
  # 3. select only top 15 genus
  tax_select(Genera_tokeep, ranks_searched = "Genus", 
             strict_matches=TRUE,n_typos = 0)

bbgut_ps_relative_melt<-psmelt(bbgut_ps_relative_top15)

#set threshold
ab_detect_threshold <- 0.005 
#Table showing how much rel. abundance per genus is in each sample
PresenceTable <- bbgut_ps_relative_melt %>% 
  select(Genus,Sample_sequencing_ID,Child_ID,Day,Abundance) %>% 
  filter(Abundance >= ab_detect_threshold) %>%
  arrange(Genus,Child_ID,Day) #sort by genus, child and day

First_appearance_day <- PresenceTable %>% 
  pivot_wider(id_cols = "Child_ID",
              names_from = "Genus",
              values_from = "Day",
              values_fn = min) %>%
  tibble::column_to_rownames("Child_ID") 

#check in how many children you see each of the top 15 genera. If <25% remove genus from top15.
First_appearance_day %>% 
  mutate(across(colnames(.),~ ifelse(.x>0,1,0))) %>% 
  mutate(across(colnames(.),~ replace_na(.x,0))) %>%  colSums(.)

First_appearance_rank <- First_appearance_day %>%
  t() %>% as.data.frame() %>% 
  mutate(across(colnames(.),
                ~ frank(.x, ties.method = "dense"))) %>% #turn day value into rank,ties.method="dense" you break the ties if ranks are 1,1,1,4 then 4 is written as 2. 
  t() %>% as.data.frame()

First_appearance_rank_melt <- First_appearance_rank %>% 
  tibble::rownames_to_column() %>% 
  reshape2::melt() %>% 
  setnames(new=c("Child_ID","Genus","Rank"))

First_appearance_day_melt <- First_appearance_day %>% 
  tibble::rownames_to_column() %>% 
  reshape2::melt() %>% 
  setnames(new=c("Child_ID","Genus","Day"))

First_appearance_MELT <- First_appearance_day_melt %>% 
  left_join(First_appearance_rank_melt,by=c("Child_ID",'Genus'))

median_first_rank <- First_appearance_rank_melt %>% 
  group_by(Genus) %>%
  dplyr::summarize(median_rank=median(Rank, na.rm=T)) %>% #summarize one group into 1 row
  arrange(median_rank)
#write.csv(median_first_rank,"median_ranks_top15_y1and2_bbgut.csv")

#use melted df with all samples
bbgut_ps_relative_melted_allsamples<- bbgut_ps_inf %>% 
  # 1. agglomerate into Genus
  tax_glom(taxrank="Genus") %>% 
  # 2. calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU)) %>%
  # 3. select only top 15 genus
  tax_select(Genera_tokeep, ranks_searched = "Genus", 
             strict_matches=TRUE,n_typos = 0) %>%
  #melt
  ps_melt()

PresenceTable_allgenera <- bbgut_ps_relative_melted_allsamples %>% 
  select(Genus,Sample_sequencing_ID,Child_ID,Day,Abundance) %>% 
  arrange(Genus,Child_ID,Day) #sort by genus, child and day

#make a new column that defines whether the top 15 genera are present in each sample
PresenceTable_allgenera$Genus_presence<- ifelse(PresenceTable_allgenera$Abundance>=0.005,1,0)

#make new column with ranks per genus
library(r2r)
rankings_pergen<-hashmap()
gen<-as.character(median_first_rank$Genus)
ran<-median_first_rank$median_rank

rankings_pergen[gen]<-ran

PresenceTable_allgenera$Genus_rank<- unlist(rankings_pergen[PresenceTable_allgenera$Genus])


#make new column for maturation score per sample
maturation_scores_persample<-PresenceTable_allgenera %>%
  group_by(Sample_sequencing_ID) %>%
  dplyr::summarise(average_weighted_maturationscore=
                     weighted.mean(Genus_rank, w=Genus_presence))

#merge with cluster change df
clusters_change_maturationscore<-left_join(cluster_changesinf,maturation_scores_persample,by="Sample_sequencing_ID")

#calculate maturation score line on df without the setbacks
important_setbacks<- clusters_change_maturationscore %>% filter(setback =="Yes" & importance=="Yes")

#filter out important setbacks for plotting line
clusters_change_maturationscore_ns<-clusters_change_maturationscore %>% filter(!clusters_change_maturationscore$Sample_sequencing_ID %in% important_setbacks$Sample_sequencing_ID)

############################################
## Find setbacks with low maturation score
## Publication *Supplementary Figure S11a, b*
## Part of Publication *Figure 3c*
############################################

#define samples that don't have an important setback according to the first criterium
no_important_sb_IDs <- clusters_change_maturationscore_ns$Sample_sequencing_ID

#empty df
allinfs_mat_score_stepbcks<-data.frame()

#loop over all infants to do a linear model predicting the maturation score per sample/infant and create new columns indicating if a setback has a maturation score lower than the lowest predicted value (=below the confidence interval)
for (i in infants){
  #filter per infant id (i) and make new column with first criterium for important setback
  child_df<-clusters_change_maturationscore %>%
    filter(Child_ID== i) %>%
    mutate(important_setback=ifelse(Sample_sequencing_ID %in% 
           important_setbacks$Sample_sequencing_ID,"Yes","No")) 
  #keep only samples without an important setback for the linear model
  no_imp_df <- child_df %>% 
    filter(Sample_sequencing_ID %in% no_important_sb_IDs)
  #linear model
  m_matscore<-lm(data=no_imp_df,
                 average_weighted_maturationscore ~ Day + I(Day^2))
  #predict maturation score
  pred_mat <- predict(m_matscore,child_df,interval='confidence') %>%
              as.data.frame()
  pred_mat <- cbind(pred_mat,child_df) %>% data.frame()
  #create new column for samples with lower maturation score
  pred_mat2<-pred_mat %>% mutate(important_mat_setback=
                                  ifelse(important_setback =="Yes" & 
                                         average_weighted_maturationscore<lwr, 
                                         "Yes","No"))
  allinfs_mat_score_stepbcks<-rbind(allinfs_mat_score_stepbcks,pred_mat2)
}


#plot for only samples with a lower maturation score and setback
mscore<-ggplot(allinfs_mat_score_stepbcks,aes(x=Day)) +
  geom_point(aes(y=average_weighted_maturationscore,
                      color=Cluster,
                      size=important_mat_setback))+
  scale_size_manual(values = c(1,4),guide="none")+
  geom_line(aes(y=fit),color="darkgrey") +
  facet_wrap(~Child_ID, nrow=2)+
  geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.3)+
  theme_light() +
  scale_color_manual(values= as.character(pastel_colors_2))  +
  ggtitle(paste('Weigted Averages per sample, presence as >', 
                ab_detect_threshold,   sep=""))+
  xlab("Days after birth")
mscore 
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_S11b_Maturation_score_BBGUT_dmm.pdf", sep=""),width=16,height=6)
print(mscore)
dev.off()


# Repeat setback plot with only the samples that have a lower maturation score
setback_plot2 <- ggplot(allinfs_mat_score_stepbcks, aes(x=Day, y=Cluster, color=Cluster, group=1)) + 
  geom_point(data = allinfs_mat_score_stepbcks[allinfs_mat_score_stepbcks$important_mat_setback=="Yes",], 
             aes(x=Day,y=Cluster, color=Cluster, group=Child_ID),size = 4)+# 
  geom_point(size = 1)+
  scale_color_manual(values= as.character(pastel_colors_2))  +
  facet_wrap(~Child_ID, nrow = 2)+
  geom_vline(xintercept=365, linetype="dotted")+
  theme_light() +
  scale_y_discrete(expand=c(0.3, 0))+
  xlab('Days after birth')
setback_plot2
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_S11a_Setbackslowmatscore_BBGUT_dmm.pdf", sep=""),
    width=16,height=6)
print(setback_plot2)
dev.off()


#plot to visualize all maturation score changes 
#(part of Publication *Figure 3c* example)
mscore2<-ggplot(allinfs_mat_score_stepbcks,aes(x=Day)) +
  geom_point(aes(y=average_weighted_maturationscore,
                 color=Cluster,
                 size=important_setback))+
  scale_size_manual(values = c(1,4),guide="none")+
  geom_line(aes(y=fit),color="grey34") +
  facet_wrap(~Child_ID, nrow=2)+
  geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.3)+
  theme_light() +
  scale_color_manual(values= as.character(pastel_colors_2))  +
  ggtitle(paste('Weigted Averages per sample, presence as >', 
                ab_detect_threshold,   sep=""))+
  xlab("Days after birth")
mscore2
# pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
#           "Allmatscore_BBGUT_dmm.pdf", sep=""),width=16,height=6)
# print(mscore2)
# dev.off()

####################################################
### Comparison of maturation scores
### 2nd year age groups (months 13-16 vs months 17-24)
### Publication *Figure 4c*
####################################################
lockdown_samples<-read.csv("Data/lockdown_samples_y2.csv") %>% pull(x)
prelockdown_samples<-read.csv("Data/prelockdown_samples_y2.csv") %>% pull(x)

get_num_month <- function(days){
  return(ifelse(days<=10,0, ((days%/%31)+1) ))}

mat_score2<-allinfs_mat_score_stepbcks %>% 
  filter(Sample_sequencing_ID %in% HTP_samples) %>%
  mutate(Month=get_num_month(Day))

### maturation score per month category
mat4 <- mat_score2 %>%
  filter(Month %in% c(13,14,15,16,17,18,19,20,21,22,23,24)) %>%
  mutate(Month_grouped = case_when(
    Sample_sequencing_ID %in% prelockdown_samples ~ "months13-16",
    Sample_sequencing_ID %in% lockdown_samples ~ "months17-24",
    TRUE ~ as.character(Month))) %>%
  dplyr::mutate(Month_grouped=factor(Month_grouped,
                                     levels=c("months13-16","months17-24"))) %>%
  group_by(Month_grouped, Child_ID) %>%
  dplyr::summarise(mean_scorebymonth = mean(average_weighted_maturationscore, na.rm = TRUE))

pastel_colors_ids <- c("I01" = "#6B96B6", "I02" = "#8FC3D4", "I03" = "#C4E3F0", "I04" = "#F8BDAF", 
                       "I05" = "#FFD9AC", "I06" = "#8EB9B7", "I07" = "#AEC89E", "I08" = "#CBBE76", 
                       "I09" = "#F5DDB7", "I10" = "#75AEBD", "I11" = "#A6CFCF", "I12" = "#E78F8D", 
                       "I13" = "#FFC1BE", "I14" = "#A39F96", "I15" = "#D3A9B9", "I16" = "#FED4E3", 
                       "I17" = "#C694B3", "I18" = "#E7C6D5", "I19" = "#BAAFA5", "I20" = "#E4CEC2")

#plot
mscorebycat<-mat4 %>% ggplot(aes(x = as.factor(Month_grouped), y = mean_scorebymonth, group = Child_ID, color = as.factor(Child_ID))) +
  geom_point(aes(color = as.factor(Child_ID)), size = 5) +
  geom_line(aes(group = Child_ID), alpha = 0.8) +
  geom_boxplot(aes(group = Month_grouped), fill="lightgrey",alpha = 0.6, outlier.shape = NA) +
  scale_color_manual(values=pastel_colors_ids)+
  theme_light() +
  labs(
    x = "Month",
    y = "Maturation Score",
    colour = "Child ID"
  )
mscorebycat
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_4c_Maturation_score_agegroupsBBGUT.pdf", sep=""),width=10,height=5)
print(mscorebycat)
dev.off()


### stats ###
# Reshape data so we can compute paired differences
paired_wide <- mat4 %>%
  select(Child_ID, mean_scorebymonth, Month_grouped) %>%
  pivot_wider(names_from = Month_grouped, values_from = mean_scorebymonth)

# Perform paired Wilcoxon signed-rank test
wilcox_test <- wilcox.test(paired_wide$`months13-16`, paired_wide$`months17-24`,paired = TRUE)
print(wilcox_test)


#####################################################
## Setback association with disease events
## Publication *Figure 3d*
#####################################################

#Explore Setback metadata in Data/setbacks.xlsx
#Figure is made according to setback.xlsx information

#explore metadata
# meta<- data.frame(sample_data(bbgut_ps_inf))
# meta %>%
#   select(Child_ID, disease_event, disease_event_type) %>%
#   # keep only unique combinations
#   distinct() %>%
#   # count disease event types
#   count("disease_event_type")

#plot barplot
event_order<-c("ADHC_disease","HTP","AB+GI_events","AB_events",
                 "GI_events","Other")
event_colors<-c(
  "HTP" = "#6B96B6",
  "AB_events" = "#CBBED9",
  "GI_events"      = "#DED6E5", 
  "AB+GI_events"          =  "#B7A5C6", 
  "Other"  = "#F2EBF5"
)


df_plot<-openxlsx::read.xlsx("Data/plot_setbacks.xlsx")

change_bar<-df_plot %>% 
mutate(CATEGORY=factor(CATEGORY,levels=event_order)) %>%
 ggplot(aes(x=GROUP, y=as.numeric(PERCENTAGE), fill=CATEGORY)) +
 geom_bar(stat="identity", position="stack") +
 ylab("Setback correspondence (%) to each group")+
 xlab("") +
 scale_fill_manual(values=event_colors)+
scale_y_continuous(breaks = seq(0, 13, by = 1)) + 
theme_light()
change_bar
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
           "Figure_3d_change&events_barplot_percent.pdf", sep=""),width=6,height=6)
print(change_bar)
dev.off() 


# CONTIGENCY TEST (Supplementary Table S8)
#determine if there is a significant association between two categorical variables.In this case, we are interested in whether the different categories (HTP, ADHC disease) are associated with the outcomes ("YES" and "NO").
# You want to compare the counts that you observe with counts you would expect if there was no association

#prepare data HTP vs ADHC DISEASE
#based on the results
table_stp<-data.frame(
    category = c("HTP", "ADHC"),
    yes = c(15, 16),
    no = c(471, 111))

row.names(table_stp)<-NULL
table_stp2<-table_stp %>% tibble::column_to_rownames( "category")


#Chi-square test: <0.005 There is enough evidence to suggest that the type of event is associated with the outcome
chisq_stp<-chisq.test(table_stp2,correct=TRUE)
chisq_stp


# Look into residuals:Residuals larger than ±2 (absolute value) suggest significant deviations. For example, if a particular event has a residual of 3, it means that the observed count for that event is significantly higher than expected under the null hypothesis.

residuals <- chisq_stp$residuals
print(residuals)

#standardizes res
std_residuals <- chisq_stp$stdres
print(std_residuals)
