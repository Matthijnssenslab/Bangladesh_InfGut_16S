############################################################################
## SEGATELLA ANALYSIS — BBGUT COHORT
## ===================================
## - Segatella prevalence across GMM stages and health groups
## - High Segatella abundance (HSA) threshold calculation
## - Monthly proportion of HSA in ad-hoc disease and HTP samples
## - HSA metadata analysis
## - Heatmap of HSA spikes per child over time
##
## Output figures & tables:
##   Figure 5a, 5b, 5c  |  Supplementary Figures S5a, S5b, S15a, S15b
##   Supplementary Table S11
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(phyloseq)
library(vegan)
library(tidyr)
library(microbiome)
library(microViz)
library(tibble)
library(ggplot2)
library(reshape2)
library(stringr)
library(data.table)

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")

meta_ps<-data.frame(sample_data(bbgut_ps))
   
sample_data(bbgut_ps)<-sample_data(meta_ps)    

#define healthy samples
HTP_samples<-sample_data(bbgut_ps)$Sample_sequencing_ID[sample_data(bbgut_ps)$sample_category == "HTP"] %>%na.omit()

maternal_samples <- sample_data(bbgut_ps)$Sample_sequencing_ID[sample_data(bbgut_ps)$fecalsample_type == "maternal"]

#Define ADHC disease samples
adhc_disease<-sample_data(bbgut_ps)$Sample_sequencing_ID[sample_data(bbgut_ps)$sample_type %in% c("D","N")] %>%na.omit()


#add cluster information per sample to ps object

#read per sample stats file generated for Figure 1 that includes GMM info
per_sample_stats<-read.csv("Data/Per_sample_stats.csv")

clust_bbgut<-per_sample_stats %>% select(X, Cluster)
clust_bbgut$Sample_sequencing_ID<-clust_bbgut$X
clust_bbgut$X<-NULL

#merge ps object with cluster info
meta_bbgut<-sample_data(bbgut_ps) %>% data.frame() %>% 
  left_join(., clust_bbgut, by="Sample_sequencing_ID") %>% tibble::column_to_rownames(.,"Sample_sequencing_ID")

meta_bbgut$Sample_sequencing_ID<-rownames(meta_bbgut)
sample_data(bbgut_ps)<-meta_bbgut

#agglomerate to genus
genus_agglom<-tax_glom(bbgut_ps, taxrank="Genus")

genus_agglomcomp <- transform_sample_counts(genus_agglom, function(OTU) OTU/sum(OTU))
genus_agglom_abundance_comp <- psmelt(genus_agglomcomp)
#replace maternal samples falling in cluster C with C_m in Cluster column
genus_agglom_abundance_comp2 <- genus_agglom_abundance_comp %>%
  mutate(Cluster = ifelse(fecalsample_type == "maternal", "C_m", Cluster))

genus_agglom_abundance_comp<- genus_agglom_abundance_comp2[, c('Abundance','Sample_sequencing_ID','Genus','Cluster','Child_ID')]


##########################################################
# 1. Visualize Segatella prevalence across GMM stages 
# in all, healthy and ad hoc disease associated samples
# Publication *Supplementary Figure S15a*
##########################################################  

sega <- genus_agglom_abundance_comp %>% 
  filter(Genus == "Segatella") %>%
  filter(Sample_sequencing_ID %in% c(HTP_samples, adhc_disease, maternal_samples)) %>%
  mutate(
      Group = case_when(
        Sample_sequencing_ID %in% adhc_disease     ~ "disease",
        Sample_sequencing_ID %in% HTP_samples      ~ "healthy",
        TRUE                                       ~ "maternal"
      ),
      Presence = ifelse(Abundance >= 0.005, 1, 0)
    )
#write.csv(sega,"segatella_acrossclusters_healthstatus.csv")

# Calculate prevalence
prevalence_df <- sega%>%
    group_by(Cluster, Group) %>%
    dplyr::summarise(prevalence = mean(Presence == 1), .groups = "drop")
  
#make new df with prevalence per group for all samples
prevalence_all <-genus_agglom_abundance_comp %>% 
    filter(!Sample_sequencing_ID %in% maternal_samples) %>%
    filter(Genus == "Segatella") %>%
    dplyr::mutate(Presence = ifelse(Abundance >= 0.005, 1, 0)) %>%
    group_by(Cluster) %>%
    dplyr::summarise(prevalence=mean(Presence==1),.groups = "drop") %>%
    mutate(Group="all samples")
  
#merge
prevalence_merged <- bind_rows(prevalence_df, prevalence_all)

#plot
sp <- ggplot(prevalence_merged, aes(x = Cluster, y = prevalence, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  labs(
    y = "Segatella Prevalence",
    x = "GMM stage"
  ) +
  scale_fill_manual(values = c(
    "all samples" = "#99A9BD",
    "disease"     = "#CBBED9",
    "healthy"     = "#6B96B6",
    "maternal"    = "#7FAFBF"
  )) +
  theme_light()
sp
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_S15a_prevelanceGMMs.pdf", sep=""),width=14,height=5)
print(sp)
dev.off()

##########################################################
# PcoA to visualize Segatella differences by GMM stage
# in all, healthy and ad hoc disease associated samples
# Publication *Supplementary Figure S15b*
##########################################################

SHAPES <- c('infant' = 16,'maternal' = 24)

bbgut_relabund<- transform_sample_counts(bbgut_ps, function(OTU) OTU/sum(OTU))
bbgut_relabund_agglom<-tax_glom(bbgut_relabund, taxrank="Genus")                       
allinfant_pcoa_mi <- ordinate(physeq=bbgut_relabund_agglom,method="PCoA",distance = 'bray')

#use melted df to filter genus_agglomcomp
melt_bbgut<-ps_melt(bbgut_relabund_agglom)
Abund <- melt_bbgut[,c("Sample_sequencing_ID","Abundance","Genus")]
Abund_prev <- Abund %>% filter(Genus=="Segatella")

#check order first
#add column in sample data of phyloseq object
sample_data(bbgut_relabund_agglom)$Abund_prev<- Abund_prev$Abundance

bbgut_ps_healthy <-bbgut_relabund_agglom %>% 
  ps_filter(Sample_sequencing_ID %in% c(HTP_samples,maternal_samples))

bbgut_ps_AD <-bbgut_relabund_agglom %>% ps_filter(Sample_sequencing_ID %in% c(adhc_disease,maternal_samples))


### plot PcoA ###

#ALL SAMPLES
p_prev_all<-  plot_ordination(physeq=bbgut_relabund_agglom, 
                          ordination = allinfant_pcoa_mi, color="Abund_prev", shape = "fecalsample_type")+
  scale_color_gradient(low = "#6B96B6",
                      high = "#FF5252",
                      guide = "colorbar")+
  theme_light()+
  scale_shape_manual(values = SHAPES)+
  scale_size_manual(values=2)+
  labs( color ="Prevotella Abundance")+
  coord_fixed()  ## need aspect ratio of 1!
p_prev_all
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_S15b_ALLsamples.pdf", sep=""),width=6,height=5)
print(p_prev_all)
dev.off()

#HEALTHY SAMPLES
p_prev_healthy<-  plot_ordination(physeq=bbgut_ps_healthy, 
                              ordination = allinfant_pcoa_mi, color="Abund_prev", shape = "fecalsample_type")+
  scale_color_gradient(low = "#6B96B6",
                       high = "#FF5252",
                       guide = "colorbar")+
  theme_light()+
  scale_shape_manual(values = SHAPES)+
  scale_size_manual(values=2)+
  labs( color ="Prevotella Abundance")+
  coord_fixed()  ## need aspect ratio of 1!
p_prev_healthy
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_S15b_HTPsamples.pdf", sep=""),width=6,height=5)
print(p_prev_healthy)
dev.off()

#AD HOC SAMPLES
p_prev_ad<-  plot_ordination(physeq=bbgut_ps_AD, 
                                  ordination = allinfant_pcoa_mi, color="Abund_prev", shape = "fecalsample_type")+
  scale_color_gradient(low = "#6B96B6",
                       high = "#FF5252",
                       guide = "colorbar")+
  theme_light()+
  scale_shape_manual(values = SHAPES)+
  scale_size_manual(values=2)+
  labs( color ="Prevotella Abundance")+
  coord_fixed()  ## need aspect ratio of 1!
p_prev_ad
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_S15b_ADHOCsamples.pdf", sep=""),width=6,height=5)
print(p_prev_ad)
dev.off()

#Data processing
get_num_month <- function(days){
  return(ifelse(days<=10,0, ((days%/%31)+1) ))
}

# modify month to be days divided by 31
genus_agglom_abundance_comp_month <- genus_agglom_abundance_comp2 %>%
  filter(fecalsample_type=="infant") %>%
  mutate(Month=get_num_month(Day))

healthy_and_adhocdisease<-c(HTP_samples,adhc_disease)


###########################################################################
# 2. Setting threshold for what is a sample with high Segatella abundance
# Publication *Supplementary Figure S5a,b*
###########################################################################

#Setting threshold
#Empirical cummulative density function (ECDF plot) that will show 
#for each Segatella abundance threshold, what's the percentage of corresponding samples. 
#Here, as a high Segatella abundance, all samples above 25% (above the mean) were chosen
#that makes up around 100-75=25% of the samples

#ecdf plot
ec<-genus_agglom_abundance_comp2 %>% 
  filter(fecalsample_type=="infant") %>%
  filter(Genus =="Segatella") %>% 
  select(Sample_sequencing_ID,Abundance,Day) %>% 
  filter(Abundance >= 0.005) %>%  
  ggplot(aes(x=Abundance)) +
  stat_ecdf() + scale_x_continuous(breaks=seq(0,1,0.05))+
  scale_y_continuous(breaks=seq(0,1,0.05)) + 
  geom_vline(xintercept = 0.25,linetype='dashed')+
  theme_light()
ec
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Suppl_figure_S5a_ec.pdf", sep=""),width=7,height=5)
print(ec)

#histogram
h<-genus_agglom_abundance_comp2 %>% 
  filter(fecalsample_type=="infant") %>%
  filter(Genus =="Segatella") %>% 
  filter(!Sample_sequencing_ID %in% HTP_samples) %>%
  select(Sample_sequencing_ID,Abundance,Day) %>% 
  filter(Abundance >= 0.005)%>%
  ggplot(aes(x=Abundance)) +
  geom_histogram(bins=50) + 
  #geom_density() +
  theme_light()+
  scale_x_continuous(breaks=seq(0,1,0.05))+
  geom_vline(xintercept = 0.25,linetype='dashed')
h
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Suppl_figure_S5b_hist.pdf", sep=""),width=7,height=5)
print(h)

#########################################################################
# 3. Monthly proportion of HSA in adhoc and HTP samples
# Publication *Figure 5a*
#########################################################################

# in every month
month_order<-c(1:12)

prev_highsamplespermonth<-genus_agglom_abundance_comp_month %>% 
  filter(fecalsample_type=="infant") %>%
  filter(Genus=="Segatella") %>% 
  filter(Abundance>=0.25) %>%
  filter(Sample_sequencing_ID %in% c(healthy_and_adhocdisease)) %>%
  mutate(Month=ifelse(Month %in% c(0,1),1,Month)) %>% #group months 0,1 into 1
  mutate(Month=factor(Month,levels=month_order)) %>%
  filter(Month %in% c(1:12)) %>%
  select(Sample_sequencing_ID,Month,Genus,Abundance) %>% 
  mutate(healthy=ifelse(Sample_sequencing_ID %in% adhc_disease,'adhoc','healthy')) %>% 
  select(Month,Sample_sequencing_ID,healthy) %>% 
  group_by(Month,healthy) %>%
  dplyr::summarise(count=n())#count how many healthy/ad hoc samples per month

#normalize the number of adhc/healthy samples per month by the total number of 
#samples in those months
samples_permonth<-genus_agglom_abundance_comp_month %>% 
  filter(fecalsample_type=="infant") %>%
  filter(Sample_sequencing_ID %in% c(healthy_and_adhocdisease)) %>% 
  mutate(Month=ifelse(Month %in% c(0,1),1,Month)) %>% #group months 0,1 into 1
  mutate(Month=factor(Month,levels=month_order)) %>%
  filter(Month %in% c(1:12)) %>%
  select(Sample_sequencing_ID,Month,Genus,Abundance) %>%
  mutate(healthy=ifelse(Sample_sequencing_ID %in% adhc_disease,'adhoc','healthy')) %>%
  select(Month,Sample_sequencing_ID,healthy) %>%
  distinct() %>% 
  group_by(Month,healthy)  %>% tally()

#normalize
plot_highprevsamples_normal<-left_join(prev_highsamplespermonth,samples_permonth,by=c("Month","healthy")) %>%
  mutate(normalized=count/n)

sum(plot_highprevsamples_normal$n)


#Define the range of months and the target health category
health_cat_target <- c('adhoc','healthy')

#Create a complete df for the given range of months
complete_df1 <- expand.grid(
  Month = as.factor(month_order),
  healthy = health_cat_target
)

# Join with the original dataset to find missing months
df_complete2 <- complete_df1 %>%
  left_join(plot_highprevsamples_normal, by = c("Month", "healthy")) %>% 
  mutate(
    normalized = ifelse(is.na(normalized), 0, normalized),
  ) %>%
  arrange(Month)

#plot
prev_normal<-ggplot(df_complete2,aes(x=as.factor(Month), y= normalized, fill=healthy)) +
  geom_col(width = 0.8, position = position_dodge(0.7),alpha=0.8)+
  scale_fill_manual(values=c("#CBBED9","#6B96B6"))+
  ylab("Proportion HSA samples (%)") +
  xlab("Time (Months)")+
  theme_light()
prev_normal
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_5a_hsa_monthlyproportions.pdf", sep=""),width=8,height=5)
print(prev_normal)
dev.off()


# Select HSA samples considering appearance of Segatella

#find first day of appearance of Segatella for every child
First_appearance_day_prev<- genus_agglom_abundance_comp_month %>% 
  filter(Genus =="Segatella") %>% 
  select(Genus,Sample_sequencing_ID,Child_ID,Day,Abundance) %>% 
  filter(Abundance >= 0.005) %>% 
  arrange(Genus,Child_ID,Day) %>%  #sort by genus, child and day
  pivot_wider(id_cols = "Child_ID",
              names_from = "Genus",
              values_from = "Day",
              values_fn = min) 
names(First_appearance_day_prev) <- c("Child_ID", "First_day")

#HSA samples
#filter all samples based on whether prevotella has appeared or not
#allsamples
prevotella_total_allsamples<-genus_agglom_abundance_comp_month %>% 
  filter(fecalsample_type=="infant") %>%
  select(Genus, Sample_sequencing_ID, Child_ID, Cluster,Day,Month, 
         Abundance, disease_event,sample_category,sample_type) %>%
  filter(Genus=="Segatella")


#merge and filter for all samples
prev_all_total<-merge(prevotella_total_allsamples,First_appearance_day_prev, by="Child_ID")
#filter samples > first day of appearance
prev_all_total<-prev_all_total %>%
  filter(Day >= First_day)

# further filter those to keep the most abundant prevotella samples > 25%
prev_all_total_mostabund<-prev_all_total %>% 
  mutate(across("Abundance", round, 2)) %>%
  filter(Abundance>=0.25)

#EXTRA FILTER: Keep only one representative sample per high prevotella abundance
#event. Keep sample with highest abundance per event
prev_all_total_mostabund_eventr<-prev_all_total_mostabund %>%
  arrange(Child_ID, disease_event, desc(Abundance)) %>% 
  mutate(disease_event=ifelse(is.na(disease_event),Sample_sequencing_ID,disease_event)) %>%
  distinct(Child_ID, disease_event, .keep_all = TRUE)


#write.csv(prev_all_total_mostabund_eventr, "HSA_meta.csv")


##############################################
# 4. HSA across life events
# Publication *Figure 5c*
##############################################

#Explore HSA metadata in Data/HSA_meta.xlsx
#Figure is made according to HSA_meta information

#plot barplot

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

df_plot<-openxlsx::read.xlsx("Data/plot_barsegatella25.xlsx")

change_bar<-df_plot %>% 
  mutate(CATEGORY=factor(CATEGORY,levels=event_order)) %>%
  ggplot(aes(x=GROUP, y=as.numeric(PERCENTAGE), fill=CATEGORY)) +
  geom_bar(stat="identity", position="stack") +
  ylab("HSA correspondence (%) to each group")+
  xlab("") +
  scale_fill_manual(values=event_colors)+
  scale_y_continuous(breaks = seq(0, 11, by = 1)) +  
  theme_light()
change_bar
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_5c_hsa_groupcorr.pdf", sep=""),width=4,height=6)
print(change_bar)
dev.off() 


###########################################################
# 5. Heatmap for samples with HSA spikes per child over time
# Publication *Figure 5b*
###########################################################

library(ggnewscale)

#1. turn segatella abundance table into presence absence table
#2. Keep Child_ID, Sample_sequencing_ID, Day,Month,Presence

prev_melted<-genus_agglom_abundance_comp_month %>% 
  filter(fecalsample_type=="infant") %>%
  filter(Genus %in% "Segatella") %>%
  #mutate Abundance instead of filtering in order to keep the 0s for plotting the   differences between months
  dplyr::mutate(Abundance=ifelse(Abundance>=0.25,Abundance,0)) %>%
  filter(Sample_sequencing_ID %in% c(healthy_and_adhocdisease)) %>%
  dplyr::mutate(Presence=ifelse(Abundance>0,1,0)) %>%
  mutate(health_cat=ifelse(Sample_sequencing_ID %in% 
                             adhc_disease,'adhoc','healthy')) %>% 
  select(Child_ID, Sample_sequencing_ID, Day,Month,Presence,health_cat)

#make heatmap that shows presence of prevotella in each health_cat over time per child
heatmap_month<-prev_melted %>% 
  mutate(Month=ifelse(Month %in% c(0,1),1,Month)) %>% #group months 0,1 into 1
  group_by(health_cat,Child_ID,Month) %>%
  dplyr::summarise(Presence_monthly=sum(Presence)) %>% 
  dplyr::mutate(Presence_monthly=ifelse(Presence_monthly>0,1,0)) 


#remake dataset for heatmap  
h<-heatmap_month %>% 
  mutate(Month_health_stat = paste(Month, health_cat, sep = "_"),
         Presence_monthly = as.factor(Presence_monthly),
         fill_color = ifelse(grepl("adhoc", Month_health_stat) & Presence_monthly == "1", "#CBBED9", ifelse(Presence_monthly == "1", "#6B96B6", "white")))


#define the range of months and the target health category
months <- 1:12
days<-1:365
health_cat_target <- c('adhoc','healthy')

#create a complete df for the given range of months
complete_df <- expand.grid(
  Month = months,
  health_cat = health_cat_target,
  Child_ID = unique(h$Child_ID)
)


#Join with the original dataset to find missing months
df_complete <- complete_df %>%
  left_join(h, by = c("Month", "health_cat", "Child_ID")) %>% 
  mutate(
    Presence_monthly = ifelse(is.na(Presence_monthly), 0, Presence_monthly),
    Month_health_stat = ifelse(is.na(Month_health_stat), paste0(Month, "_", health_cat), Month_health_stat),
    fill_color = ifelse(is.na(fill_color), "white", fill_color)
  ) %>%
  arrange(Child_ID, Month)

#sort month categories
h2<-str_sort(x=unique(df_complete$Month_health_stat), numeric = TRUE)
#order infants by grouping together the ones for which high prevotella appears first in disease samples
hpd_order<-c("I03","I05","I10","I13","I01","I17",
             "I15","I09","I06","I11","I04","I18","I20",
             "I19","I02","I07","I16","I08","I12","I14")

#plot
hpres<-  df_complete %>%
  mutate(Month_health_stat=factor(Month_health_stat,levels=h2)) %>%
  mutate(Child_ID=factor(Child_ID,levels=hpd_order)) %>%
  ggplot(aes(x = as.factor(Month_health_stat), y = as.factor(Child_ID))) +
  geom_tile(aes(fill = fill_color), color = "white") +
  geom_hline(yintercept = seq(1.5, 30, 1),color="lightgrey") +
  geom_vline(xintercept = seq(2.5, 50, 2),color="darkgrey") +
  scale_fill_identity() +  # Use the exact color values provided in fill_color
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust = 1),
        panel.grid.major = element_blank())+
  xlab("Month/health category") +
  ylab ("Infant ID")
hpres  
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_5b_hsa_infantpresence.pdf", sep=""),width=14,height=10)
print(hpres)
dev.off()  


#-----Contingency test------#
#Supplementary Table S11
#prepare data HTP vs ADHC DISEASE
#based on HSA_meta.xlsx

ptable_stp<-data.frame(
  category = c("HTP", "ADHC"),
  yes = c(19, 14),
  no = c(467, 113))
row.names(ptable_stp)<-NULL
ptable_stp2<-ptable_stp %>% tibble::column_to_rownames( "category")


# 1. Chi-square test: <0.005 There is enough evidence to suggest that the type of event is associated with the outcome
pchisq_stp<-chisq.test(ptable_stp2,correct=TRUE)
pchisq_stp
# data:  ptable_stp2
# X-squared = 8.6568, df = 1, p-value = 0.003258

presiduals <- pchisq_stp$residuals
print(presiduals)
# yes         no
# HTP  -1.400419  0.3340420
# ADHC  2.739519 -0.6534575

std_residuals <- pchisq_stp$stdres
print(std_residuals)
# yes        no
# HTP  -3.163026  3.163026
# ADHC  3.163026 -3.163026

