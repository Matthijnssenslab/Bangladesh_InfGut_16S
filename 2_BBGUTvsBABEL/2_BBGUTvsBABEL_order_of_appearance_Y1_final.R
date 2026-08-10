############################################################################
## BBGUT vs BABEL: ORDER OF GENUS APPEARANCE (YEAR 1)
## ====================================================
## - Top 15 most abundant genera per cohort during year 1
## - Order of appearance of most abundant genera
## - Rank correlations between cohorts
## - Bacterial prevalence correlations
## - Differential abundance analysis
##
## Output figures:
##   Figure 2b, 2d  |  Supplementary Figures S10a, S10b, S8a, S8b
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(phyloseq)
library(tidyr)
library(microbiome)
library(ggplot2)
library(r2r)
library(glue)
library(data.table)
library(scales)
library(microViz)

# ---- COLOR PALETTE ----
pastel_colors <- c("#6B96B6", "#8FC3D4", "#C4E3F0", "#F8BDAF", "#FFD9AC",
                   "#8EB9B7", "#AEC89E", "#CBBE76", "#F5DDB7", "#75AEBD",
                   "#A6CFCF", "#E78F8D", "#FFC1BE", "#A39F96", "#D3A9B9",
                   "#FED4E3", "#C694B3", "#E7C6D5", "#BAAFA5", "#E4CEC2",
                   "lightgrey")

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")
load("Data/babel_phyloseq_final_newtaxonomy.RData")

#######################
#BBGUT data processing
#######################
#Keep 1st year healthy samples only
bbgut_ps_y1<-bbgut_ps %>% ps_filter(fecalsample_type== "infant") %>%
  ps_filter(Day<=365)

HTP_samples<-sample_data(bbgut_ps_y1)$Sample_sequencing_ID[sample_data(bbgut_ps_y1)$sample_category == "HTP"] %>%na.omit()

htp_filtered_bbgut_ps  <- ps_filter(bbgut_ps_y1, 
                                         Sample_sequencing_ID %in% 
                                           HTP_samples)

htp_filtered_bbgut_ps_melt <- psmelt(htp_filtered_bbgut_ps)

#######################
#BABEL data processing
#######################
#Keep 1st year healthy samples only
babel_ps_y1<-babel_ps %>% 
  ps_filter(X.days<=365) %>%
  ps_filter(InfantID !="S011")

LDA_Samples <-sample_data(babel_ps_y1)$Sample.ID[sample_data(babel_ps_y1)$LDA == 1] %>%na.omit()
htp_filtered_babel_ps  <- ps_filter(babel_ps_y1, 
                                         Sample.ID %in% 
                                           LDA_Samples)

htp_filtered_babel_ps_melt <- psmelt(htp_filtered_babel_ps)

#######################
# 1.TOP 15 GENERA BBGUT
#######################

tot_ab_htp_filtered <- htp_filtered_bbgut_ps_melt %>% 
  group_by(Genus) %>% 
  # get total abundance per genus
  dplyr::summarize(total_abundance=sum(Abundance)) %>% 
  # sort descending
  arrange(desc(total_abundance)) %>% 
  # add new column with relative ab
  mutate(relative_ab=round(total_abundance/sum(total_abundance)*100,2))

tot_ab_htp_filtered %>% # View( )
Genera_tokeep <- head(tot_ab_htp_filtered$Genus,15)


#calculate % of read abundance that the top 15 genera make overall
top15_perc <- head(tot_ab_htp_filtered$relative_ab,15) %>% sum() %>% round(2)
print(glue("Top 15 genera: {top15_perc}% of total read abundance")) # 95.9%

######################
# 1.TOP 15 GENERA BABEL
######################

tot_ab_htp_filtered_babel <- htp_filtered_babel_ps_melt %>% 
  group_by(Genus) %>% 
  # get total abundance per genus
  dplyr::summarize(total_abundance=sum(Abundance)) %>% 
  # sort descending
  arrange(desc(total_abundance)) %>% 
  # add new column with relative ab
  mutate(relative_ab=round(total_abundance/sum(total_abundance)*100,2))

tot_ab_htp_filtered_babel %>% # View( )
Genera_tokeep_babel <- head(tot_ab_htp_filtered_babel$Genus,15)

#calculate % of abundance that the top 15 genera make overall
top15_perc_bab <- head(tot_ab_htp_filtered_babel$relative_ab,15) %>% sum() %>% round(2)
print(glue("Top 15 genera babel: {top15_perc_bab}% of total abundance")) #95.3%


##################################################################################
### UNION TOP15 BBGUT TOP15 BABEL (calculate top15 per cohort first) #####
##################################################################################
#union takes what doesn't overlap in the combined vectors
gen_tokeep_top15bbgutbabel<-purrr::reduce(list(Genera_tokeep,Genera_tokeep_babel),union)

################################################
# BBGUT per sample relative ab of top 15 genera
################################################

# phyloseq object for top 15 relative ab genera
bbgut_ps_relative_top15 <- htp_filtered_bbgut_ps %>% 
  # 1. agglomerate healthy_ps into Genus
  tax_glom(taxrank="Genus") %>% 
  # 2. calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU)) %>%
  # 3. select only top 15 genus
  tax_select(gen_tokeep_top15bbgutbabel, ranks_searched = "Genus", 
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


#manually remove the genera that are common between bbgut and babel because they are present in < 25% of the infants (5/20) - ok if present in at least 6 infants in BBGUT
First_appearance_day$Parabacteroides<-NA
First_appearance_day$Sutterella<-NA
First_appearance_day$Thomasclavelia<-NA

First_appearance_rank <- First_appearance_day %>%
  t() %>% as.data.frame() %>% 
  mutate(across(colnames(.),
                ~ frank(.x, ties.method='dense',na.last = "keep"))) %>% #turn day value into rank,ties.method="dense" you break the ties if ranks are 1,1,1,4 then 4 is written as 2. 
  t() %>% as.data.frame()

First_appearance_rank_melt <- First_appearance_rank %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  setnames(new=c("Child_ID","Genus","Rank"))

First_appearance_day_melt <- First_appearance_day %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  setnames(new=c("Child_ID","Genus","Day"))

First_appearance_MELT <- First_appearance_day_melt %>% 
  left_join(First_appearance_rank_melt,by=c("Child_ID",'Genus'))


#write.csv(First_appearance_MELT,"BBGUT_order_of_appearance.csv")

median_first_rank <- First_appearance_rank_melt %>% 
  group_by(Genus) %>%
  dplyr::summarize(median_rank=median(Rank, na.rm=T)) %>% #summarize one group into 1 row
  mutate(ranked_medians=frank(median_rank,ties.method='dense',na.last = "keep")) %>%
  arrange(ranked_medians) %>%
  select(-median_rank)

median_first_day <- First_appearance_day_melt %>% 
  group_by(Genus) %>%
  dplyr::summarize(median_day=median(Day, na.rm=T)) %>% 
  arrange(median_day)

##################################################
#### BABEL per sample relative ab of top 15 genera
##################################################

# phyloseq object for top 15 relative ab genera
babel_ps_relative_un <- htp_filtered_babel_ps %>% 
  # 1. agglomerate healthy_ps into Genus
  tax_glom(taxrank="Genus") %>% 
  # 2. calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU)) %>% 
  # 3. select only union genera
  tax_select(gen_tokeep_top15bbgutbabel, ranks_searched = "Genus", 
             strict_matches=TRUE,n_typos = 0)

babel_ps_relative_melt<-psmelt(babel_ps_relative_un)

#take a look at which genera are missing in the BABEL cohort
setdiff(gen_tokeep_top15bbgutbabel, babel_ps_relative_melt %>% pull(Genus) %>% unique())
#set threshold
ab_detect_threshold_babel <- 0.005 
#Table showing how much rel. abundance per genus is in each sample
PresenceTable_babel <- babel_ps_relative_melt %>% 
  select(Genus,Sample.ID,InfantID,X.days,Abundance) %>% 
  filter(Abundance >= ab_detect_threshold_babel) %>%
  arrange(Genus,InfantID,X.days) #sort by genus, child and day

First_appearance_day_babel <- PresenceTable_babel %>% 
  pivot_wider(id_cols = "InfantID",
              names_from = "Genus",
              values_from = "X.days",
              values_fn = min) %>%
  tibble::column_to_rownames("InfantID")

#check in how many children you see each of the top 15 genera. If <25% remove genus from top15.
First_appearance_day_babel %>% 
  mutate(across(colnames(.),~ ifelse(.x>0,1,0))) %>% 
  mutate(across(colnames(.),~ replace_na(.x,0))) %>%  colSums(.)

#manually add missing genera from union "g_Prevotella" "g_Ligilactobacillus"   "g_Limosilactobacillus" "g_Parolsenella" 
First_appearance_day_babel$Segatella<-NA
First_appearance_day_babel$Ligilactobacillus<-NA
First_appearance_day_babel$Limosilactobacillus<-NA
First_appearance_day_babel$Parolsenella<-NA

First_appearance_rank_babel <- First_appearance_day_babel %>%
  t() %>% as.data.frame() %>% 
  mutate(across(colnames(.),
                ~ frank(.x, ties.method = "dense",na.last="keep"))) %>% #turn day value into rank,ties.method="dense" you break the ties if ranks are 1,1,1,4 then 4 is written as 2. 
  t() %>% as.data.frame()

First_appearance_rank_babel_melt <- First_appearance_rank_babel %>% 
  tibble::rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  setnames(new=c("InfantID","Genus","Rank"))

First_appearance_day_babel_melt <- First_appearance_day_babel %>% 
  as.data.frame() %>%
  tibble::rownames_to_column() %>% 
  pivot_longer(-rowname) %>% 
  setnames(new=c("InfantID","Genus","Day"))

First_appearance_babel_MELT <- First_appearance_day_babel_melt %>% 
  left_join(First_appearance_rank_babel_melt,by=c("InfantID",'Genus'))

#write.csv(First_appearance_babel_MELT,"BABEL_order_of_appearance.csv")

median_first_rank_babel <- First_appearance_rank_babel_melt %>% 
  group_by(Genus) %>% 
  dplyr::summarize(median_rank=median(Rank, na.rm=T)) %>% #summarize one group into 1 row
  mutate(ranked_medians=frank(median_rank,ties.method='dense',na.last = "keep")) %>%
  arrange(ranked_medians) %>%
  select(-median_rank)
  
median_first_day_babel <- First_appearance_day_babel_melt %>% 
  group_by(Genus) %>%
  dplyr::summarize(median_day=median(Day, na.rm=T)) %>% 
  arrange(median_day)


####################################################
## 2. Order of Appearance of most abundant genera (y1)
## Publication *Supplementary Figure S10a,b*
####################################################

# #color by phylum so create new column in dataframe that assigns phylum per genus
phylumname_bbgut<-hashmap()
phylumname_babel<-hashmap()

gen_name_bbgut<-as.character(as.data.frame(tax_table(bbgut_ps))$Genus)
gen_name_babel<-as.character(as.data.frame(tax_table(babel_ps))$Genus)

phy_bbgut<-as.data.frame(tax_table(bbgut_ps))$Phylum
phy_babel<-as.data.frame(tax_table(babel_ps))$Phylum

phylumname_bbgut[gen_name_bbgut]<-phy_bbgut
phylumname_babel[gen_name_babel]<-phy_babel


First_appearance_MELT$Phylum<-unlist(phylumname_bbgut[as.character(First_appearance_MELT$Genus)])
First_appearance_babel_MELT$Phylum<-unlist(phylumname_bbgut[as.character(First_appearance_babel_MELT$Genus)])

#BBGUT 
ranksday_boxplot_bbgutnoday<-First_appearance_MELT %>%
  gather(key='variable',value='measurement',Rank) %>%
  mutate(variable=factor("Rank")) %>%
  mutate(Genus = factor(Genus,levels=median_first_rank$Genus)) %>%
  filter(Genus %in% Genera_tokeep) %>%
  ggplot(aes(x=measurement,y=Genus)) +
  geom_boxplot(aes(fill=Phylum),na.rm = TRUE,outlier.color="darkgrey") +
  scale_fill_manual(values=pastel_colors)+
  theme_light()+
  scale_x_continuous(breaks= pretty_breaks())
ranksday_boxplot_bbgutnoday
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_S10b_Order_of_appearance_BBGUT_y1.pdf", sep=""),width=7,height=7)
print(ranksday_boxplot_bbgutnoday)
dev.off()


# #BABEL
ranksday_boxplot_babel<-First_appearance_babel_MELT %>%
  gather(key='variable',value='measurement',Rank) %>%
  mutate(variable=factor(variable,levels=c("Rank","Day"))) %>%
  mutate(Genus = factor(Genus,levels=median_first_rank_babel$Genus)) %>%
  filter(Genus %in% Genera_tokeep_babel) %>%
  ggplot(aes(x=measurement,y=Genus)) +
  facet_wrap(~variable,  ncol=2, scales='free_x') +
  geom_boxplot(aes(fill=Phylum),na.rm = TRUE, outlier.color="darkgrey") +
  scale_fill_manual(values=pastel_colors)+
  theme_light()+
  scale_x_continuous(breaks= pretty_breaks())
ranksday_boxplot_babel
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figura_S10a_Order_of_appearance_BABEL_y1.pdf", sep=""),width=7,height=7)
print(ranksday_boxplot_babel)
dev.off()


#--- STATISTICS for consistancy in rankings ---#


### BABEL ###

First_appearance_rank_babel_kendall<-First_appearance_rank_babel %>%
  #keep only genera to consider for the test
  select(-c(Parolsenella,Segatella,Limosilactobacillus,
            Ligilactobacillus))
  
KENDALOUTCOMEbabel <-kendall.w(t(First_appearance_rank_babel_kendall), nrands = 10000,
                          type = 1, quiet = T)
KENDALOUTCOME_tablebabel <- as.data.frame(matrix(nrow = 4, ncol=2))
KENDALOUTCOME_tablebabel[,1] <- c("Kendall's W (uncorrected for ties)",
                             "Kendall's W (corrected for ties)",
                             "Spearman's ranked correlation",
                             "Kendall's W p-value (one-tailed test [greater])")
KENDALOUTCOME_tablebabel[,2] <- c(KENDALOUTCOMEbabel$w.uncorrected, KENDALOUTCOMEbabel$w.corrected,
                             KENDALOUTCOMEbabel$spearman.corr, KENDALOUTCOMEbabel$pval.rand)
colnames(KENDALOUTCOME_tablebabel) <- c('Kendall test, all infants','value')
write.csv(KENDALOUTCOME_tablebabel, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "kendall_orderofappearance_BABEL_y1.csv"))


### BBGUT ### 
First_appearance_rankkendall<- First_appearance_rank %>%
  #keep only genera to consider for the test
  select(-c(Parabacteroides,Sutterella,Thomasclavelia))

KENDALOUTCOMEbbgut <-kendall.w(t(First_appearance_rankkendall), nrands = 10000,
                          type = 1, quiet = T)
KENDALOUTCOME_tablebbgut <- as.data.frame(matrix(nrow = 4, ncol=2))
KENDALOUTCOME_tablebbgut[,1] <- c("Kendall's W (uncorrected for ties)",
                             "Kendall's W (corrected for ties)",
                             "Spearman's ranked correlation",
                             "Kendall's W p-value (one-tailed test [greater])")
KENDALOUTCOME_tablebbgut[,2] <- c(KENDALOUTCOMEbbgut$w.uncorrected, KENDALOUTCOMEbbgut$w.corrected, 
                             KENDALOUTCOMEbbgut$spearman.corr, KENDALOUTCOMEbbgut$pval.rand)
colnames(KENDALOUTCOME_tablebbgut) <- c('Kendall test, all infants','value')
write.csv(KENDALOUTCOME_tablebbgut, file = paste0(workdir, format(Sys.time(), "%Y-%m-%d"), "kendall_orderofappearance_BBGUT_y1.csv"))

##################################################################
## 3. Rank correlations between the two cohorts
## Publication *Figure 2d*
##################################################################

#find common genera between the two cohorts but 
#gen_tokeep_top15bbgutbabel

#merge median rank datasets
median_rank_babel<-First_appearance_babel_MELT %>% 
  group_by(Genus) %>% 
  dplyr::summarize(median_rank_BABEL=median(Rank, na.rm=T))%>%
  mutate(ranked_medians_BABEL=frank(median_rank_BABEL,ties.method='dense',na.last = "keep")) %>%
  mutate(cohort="BABEL") %>%
  select(-median_rank_BABEL)

#all infants babel
median_rank_babel_allinf<-First_appearance_babel_MELT %>% 
  group_by(Genus,InfantID) %>% 
  dplyr::summarize(median_rank=median(Rank, na.rm=T))%>%
  mutate(ranked_medians=frank(median_rank,ties.method='dense',na.last = "keep")) %>%
  mutate(cohort="BABEL") %>%
  dplyr::rename(Child_ID=InfantID)

median_rank_bbgut<-First_appearance_MELT %>% 
  group_by(Genus) %>% 
  dplyr::summarize(median_rank_BBGUT=median(Rank, na.rm=T)) %>%
  mutate(ranked_medians_BBGUT=frank(median_rank_BBGUT,ties.method='dense',na.last = "keep")) %>%
  mutate(cohort="BBGUT") %>%
  select(-median_rank_BBGUT)

#all infants bbgut
median_rank_bbgutallinfs<-First_appearance_MELT %>% 
  group_by(Genus, Child_ID) %>% 
  dplyr::summarize(median_rank=median(Rank, na.rm=T)) %>%
  mutate(ranked_medians=frank(median_rank,ties.method='dense',na.last = "keep")) %>%
  mutate(cohort="BBGUT")
  

#merge
allranks<-left_join(median_rank_bbgut, median_rank_babel, by="Genus")

#fix NAs by turning missing genera ranks into highest rank for plotting
median_maxbbgut<-max(allranks$ranked_medians_BBGUT,na.rm = TRUE)
median_maxbabel<-max(allranks$ranked_medians_BABEL,na.rm = TRUE)

missing_bbgut<- c("Thomasclavelia","Parabacteroides","Sutterella")
missing_babel<-c("Segatella","Ligilactobacillus",
                 "Limosilactobacillus", "Parolsenella")

allranks_fix<- allranks %>%
  mutate(ranked_medians_BBGUT=ifelse(is.na(ranked_medians_BBGUT),median_maxbbgut+5,ranked_medians_BBGUT)) %>%
 mutate(ranked_medians_BABEL=ifelse(is.na(ranked_medians_BABEL),median_maxbabel+5,ranked_medians_BABEL))    



### Rank correlations plot


genus_colors <- c(
  "Bacteroides"                = "#BAAFA5",   
  "Bifidobacterium"            = "#6B96B6",   
  "Clostridium_sensu_stricto"  = "#8FC3D4",
  "Enterococcus"               = "#FFD9AC", 
  "Thomasclavelia"             = "#99A9BD",
  "Escherichia/Shigella"       = "#8EB9B7",
  "Enterobacter"               = "#C694B3",
  "Faecalibacterium"           = "#AEC89E",
  "Haemophilus"                = "#fff399",
  "Klebsiella"                 = "#E3E1DA",
  "Lactobacillus"              = "#F8BDAF",
  "Ligilactobacillus"          = "#F5DDB7",
  "Limosilactobacillus"        = "#75AEBD",
  "Megasphaera"                = "#A6CFCF",
  "Parolsenella"               = "#FFC1BE",
  "Phocaeicola"                = "#C4E3F0",
  "Segatella"                  = "#E78F8D",
  "Staphylococcus"             = "#FED4E3",
  "Streptococcus"              = "#D3A9B9",
  "Sutterella"                 = "#CBBE76",
  "Veillonella"                = "#E4CEC2",
  "Parabacteroides"            = "#778899",
  "Dorea"                      = "#A39F96")

# Create a new column for shape based on conditions. Different shape is given for genera present in only one of the two cohorts.
allranks_fix$shape <- ifelse(
  allranks_fix$Genus %in% missing_babel & allranks_fix$ranked_medians_BABEL == 16, "Present in BBGUT only",
  ifelse(allranks_fix$Genus %in% missing_bbgut & allranks_fix$ranked_medians_BBGUT == 11, "Present in BABEL only", "Present in both cohorts")
)

# Map shapes to specific values
shape_values <- c("Present in both cohorts" = 21, "Present in BBGUT only" = 22, "Present in BABEL only" = 23)

## predict the regression for all the data although it's calculated only on the common + FORCE REGRESSION TO 0,0 for plotting

allranksfiltered<-allranks_fix %>% filter(!Genus %in% c(missing_babel,missing_bbgut))

allranks_missingonly<-allranks_fix %>% filter(Genus %in% c(missing_babel,missing_bbgut))

# Fit linear model based on allranks_filtered
model <- lm(ranked_medians_BBGUT ~ 0+ranked_medians_BABEL, data = allranksfiltered)

# Predict over a range slightly beyond the data range
new_data <- data.frame(
  ranked_medians_BABEL = seq(min(allranksfiltered$ranked_medians_BABEL) -0.1, max(allranksfiltered$ranked_medians_BABEL) + 5, length.out = 100)
)
new_data$ranked_medians_BBGUT <- predict(model, newdata = new_data)

#additional plotting fixes for overlapping genera on the plot
#fix missing data so that they don't overlap
fix_m <- allranks_missingonly %>%
  group_by(cohort.x, ranked_medians_BBGUT) %>%  # Group by cohort and rank
  mutate(ranks_plot_bbgut = ranked_medians_BBGUT + 0.06 * (row_number() - 1)) %>%
  ungroup() %>%
  group_by(cohort.y, ranked_medians_BABEL) %>%
  mutate(ranks_plot_babel = ranked_medians_BABEL + 0.06 * (row_number() - 1)) %>%
  ungroup()

fix_common <- allranksfiltered %>%
  mutate(
    ranks_plot_bbgut = ifelse(Genus %in% c("Staphylococcus", "Haemophilus"),
                              ranked_medians_BBGUT + 0.01 * (row_number() - 1),
                              ranked_medians_BBGUT),  # Keep others unchanged
    ranks_plot_babel = ifelse(Genus %in% c("Staphylococcus", "Haemophilus"),
                              ranked_medians_BABEL + 0.01 * (row_number() - 1),
                              ranked_medians_BABEL)   # Keep others unchanged
  )

#plot
cor3 <- ggplot(allranksfiltered, aes(x = ranked_medians_BABEL, y = ranked_medians_BBGUT)) +
  geom_point(data=fix_common,aes(x = ranks_plot_babel, y = ranks_plot_bbgut,fill = Genus, shape = shape, color = Genus), size = 7, alpha = 0.8, stroke = 1.2) +
  geom_smooth(method = "lm", se = TRUE, formula = y ~ 0 + x, 
              linetype = "solid", color = "black", fill = "grey82") +
  geom_line(data = new_data, aes(x = ranked_medians_BABEL, y = ranked_medians_BBGUT), 
            color = "black", linetype = "dotted") +
  geom_point(data = fix_m, 
             aes(x = ranks_plot_babel, y = ranks_plot_bbgut, 
                 fill = Genus, shape = shape, color = Genus), 
             size = 7, alpha = 0.8, stroke = 1.5) +
  labs(x = "Median Rank (BABEL)", y = "Median Rank (BBGUT)", 
       fill = "Genus", color = "Genus", shape = "Genus") +
  scale_fill_manual(values = genus_colors) +
  scale_color_manual(values = genus_colors) +
  scale_shape_manual(values = shape_values) +  # Adjust the shape values as needed
  scale_y_continuous(breaks = 1:13) +
  scale_x_continuous(breaks = 1:16) +
  theme_light()
cor3

pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_2d_Order_of_appearance_BBGUTvsBABEL_y1FIXED_predictedforced.pdf", 
          sep=""),width=12,height=8)
print(cor3)
dev.off()



#----------- KENDALL test between the bbbgut/babel ranks -------------#
#pivot df to get it in the right format
allranks_wide<-allranks %>% select(Genus,ranked_medians_BABEL, ranked_medians_BBGUT) %>%
  filter(!Genus %in% c("Thomasclavelia","Ligilactobacillus",
                       "Limosilactobacillus","Parolsenella",
                       "Segatella",
                       "Parabacteroides","Sutterella")) %>% # remove genera present in only one cohort for the kendall test
  pivot_longer(cols=-Genus, values_to = "Rank") %>% 
  pivot_wider(names_from = Genus,
              values_from = Rank) %>%
  tibble::column_to_rownames("name")


KENDALOUTCOMEcommon <-kendall.w(t(allranks_wide), nrands = 10000,
                               type = 1, quiet = T)
KENDALOUTCOME_tablecommon <- as.data.frame(matrix(nrow = 4, ncol=2))
KENDALOUTCOME_tablecommon[,1] <- c("Kendall's W (uncorrected for ties)",
                                  "Kendall's W (corrected for ties)",
                                  "Spearman's ranked correlation",
                                  "Kendall's W p-value (one-tailed test [greater])")
KENDALOUTCOME_tablecommon[,2] <- c(KENDALOUTCOMEcommon$w.uncorrected, KENDALOUTCOMEcommon$w.corrected, 
                                  KENDALOUTCOMEcommon$spearman.corr, KENDALOUTCOMEcommon$pval.rand)
colnames(KENDALOUTCOME_tablecommon) <- c('Kendall test, all cohorts','value')
#write.csv(KENDALOUTCOME_tablecommon, "Kendall_test_bbgutvsbabel_correlationgeneraappearance_onlycommongenera.csv")


########################################################################
### COMBINED DATASET ###
#1. merge phyloseq objects
#2. look into abundance patterns of most abundant genera
## Publication *Supplementary Figure S8b*
########################################################################
merged_ASV_ps<-merge_phyloseq(bbgut_ps_y1,babel_ps_y1)

#merge infant IDs
sample_data(merged_ASV_ps)$InfantID_merged <- as.character(coalesce(sample_data(merged_ASV_ps)$InfantID, sample_data(merged_ASV_ps)$Child_ID))

#merge age in days
sample_data(merged_ASV_ps)$Age_merged <- as.character(coalesce(sample_data(merged_ASV_ps)$X.days, sample_data(merged_ASV_ps)$Day))

#merge sample IDs
sample_data(merged_ASV_ps)$Sample_IDs_merged<-rownames(sample_data(merged_ASV_ps))

#make new column in data that assigns a cohort name per sample
sample_data(merged_ASV_ps)$Cohort <- ifelse(startsWith(sample_data(merged_ASV_ps)$InfantID_merged, "I"), "BBGUT", ifelse(startsWith(sample_data(merged_ASV_ps)$InfantID_merged, "S"), "BABEL", NA))

#check data is correct
sample_data(merged_ASV_ps) %>% 
  data.frame() %>% 
  filter(Cohort=='BBGUT') %>%
  pull(InfantID_merged) %>% #gets values from a column
  unique()

sample_data(merged_ASV_ps) %>% 
  data.frame() %>% 
  filter(Cohort=='BABEL') %>%
  pull(InfantID_merged) %>%
  unique()


#combine healthy sample IDs
healthy_all<-c(LDA_Samples, HTP_samples) %>% str_replace("-","_")


# phyloseq object common relative ab genera to plot
merged_asv_ps_commontoplot <- merged_ASV_ps %>% 
  #0. filter only healthy samples to keep
  subset_samples(Sample_IDs_merged %in% healthy_all) %>% #keep healthy samples only
  # 1. agglomerate healthy_ps into Genus
  tax_glom(taxrank="Genus") %>% 
  # 2. calculate relative ab
  transform_sample_counts(.,function(OTU) OTU/sum(OTU)) %>%
  # 3. select only top 15 genus
  tax_select(gen_tokeep_top15bbgutbabel, ranks_searched = "Genus", 
             strict_matches=TRUE,n_typos = 0)


# 0 is "birth" and for this I will consider the first 10 days of life based on the paper Vatanen et al. 2022
get_num_month <- function(days){
  return(ifelse(days<=10,0, ((days%/%31)+1) ))
}

# modify month to be days divided by 31
new_months_merged<- as.numeric(sample_data(merged_asv_ps_commontoplot)$Age_merged) %>% get_num_month() %>% as.factor()

sample_data(merged_asv_ps_commontoplot)$Month_merged <- new_months_merged

# group by month
months_wanted <- c(0:12)

#melt
merged_gen_commontoplot_MELT<-ps_melt(merged_asv_ps_commontoplot)


## plot rel. abundance patterns for joint genera BBGUT+BABEL during year 1
abund_21<-merged_gen_commontoplot_MELT %>%
  filter(Abundance>=0.005) %>%
  filter(Genus %in% gen_tokeep_top15bbgutbabel) %>%
  group_by(Month_merged,Cohort,Genus) %>%
  filter(Month_merged %in% months_wanted) %>%
  dplyr::summarize(median_abundance=median(Abundance)) %>%
  mutate(Month_merged=as.numeric(as.character(Month_merged)))  %>%
  
  filter(!(Genus== "Lactobacillus" & Cohort=="BABEL")) %>%
  filter(!(Genus== "Ligilactobacillus" & Cohort=="BABEL")) %>%
  filter(!(Genus== "Limosilactobacillus" & Cohort=="BABEL")) %>%
  filter(!(Genus== "Segatella" & Cohort=="BABEL")) %>%
  filter(!(Genus== "Parolsenella" & Cohort=="BABEL")) %>% 
  filter(!(Genus== "Sutterella" & Cohort=="BBGUT")) %>%
  filter(!(Genus== "Parabacteroides" & Cohort=="BBGUT")) %>%
  filter(!(Genus== "Thomasclavelia" & Cohort=="BBGUT")) %>%

  ggplot(aes(x=Month_merged,y=median_abundance, color=Cohort))+
  geom_line(size=1, alpha=0.5)+
  geom_point(size=2.5)+
  theme_light()+
  facet_wrap(~Genus, nrow=4, ncol=6, scales='free') +
  scale_color_manual(values=c("#8FC3D4" ,"#778899"))+
  scale_x_continuous(breaks=0:12)+
  ylab("Median Rel.Abundance")+
  xlab("Time(Months)")
abund_21
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_S8b_RelAbundance_BBGUTvsBABEL_21gen_y1.pdf", sep=""),
    width=18,height=10)
print(abund_21)
dev.off()
# 

# CALCULATE PERCENTAGE OF READS PER COHORT FOR most abundant GENERA
#calculate % of abundance that the top 15 genera make overall FOR COMMON GENERA = 8
top_perc_bab <- tot_ab_htp_filtered_babel %>% filter(Genus %in% gen_tokeep_top15bbgutbabel) %>% select(relative_ab) %>% head() %>% sum() %>% round(2)
print(glue("Common genera babel: {top_perc_bab}% of total abundance"))

top_perc <- tot_ab_htp_filtered %>% filter(Genus %in% gen_tokeep_top15bbgutbabel) %>% select(relative_ab) %>% head() %>% sum() %>% round(2)
print(glue("Common genera bbgut: {top_perc}% of total abundance"))


#########################################
# 4. Genera Prevalence correlations
# Publication *Figure 2b*
##########################################
#Spearman correlations between the genera bbgut vs babel
#Get correlations together with pvalues

### Correlations in prevalence of specific genera ###
#dataframe with correlations per genus (mean prevalence by month)
df_cor_cohorts_prev<-merged_gen_commontoplot_MELT %>% 
  mutate(Prevalence = ifelse(Abundance >= 0.005, 1, 0)) %>%
  filter(Genus %in% gen_tokeep_top15bbgutbabel) %>%
  group_by(Month_merged,Cohort,Genus) %>%
  filter(Month_merged %in% months_wanted) %>%
  dplyr::summarise(mean_prevalence=mean(Prevalence)) %>% 
  mutate(Month_merged=as.numeric(as.character(Month_merged))) %>%
  pivot_wider(id_cols=c(Month_merged, Genus),names_from = Cohort, 
              values_from = mean_prevalence) %>%
  #Remove genera for which there are not enough points to compute correlation
  filter(!Genus %in% c(missing_bbgut,missing_babel)) %>% 
  #replace NAs with 0 
  mutate(BABEL=replace_na(BABEL,0),
         BBGUT=replace_na(BBGUT,0)) %>%
  group_by(Genus) %>%
  dplyr::summarise(spearman_prev=cor(BABEL, BBGUT, method="spearman"),
                   spearman_p_prev = cor.test(BABEL, BBGUT, 
                                              method = "spearman",exact=FALSE)$p.value) %>%
  mutate(spearman_p_prev_adj = p.adjust(spearman_p_prev, method = "BH")) %>% 
  mutate(spearman_p_prev_adj=round(spearman_p_prev_adj,2),
         spearman_prev=round(spearman_prev,2)) %>%
  select(Genus, spearman_prev, spearman_p_prev_adj) %>% 
  arrange(spearman_p_prev_adj)


#adjust data to include missing genera from the two cohorts by adding a random non significant pvalue
df_missing_prev<- data.frame(Genus=c(missing_babel,missing_bbgut),
                             spearman_prev=0,
                             spearman_p_prev_adj=1) 
#merge with df
df_cor_cohorts_merged_prev<-rbind(df_cor_cohorts_prev,df_missing_prev)

#add column specifying significant
df_cor_cohorts2_prev<- df_cor_cohorts_merged_prev %>%
  mutate(significance=ifelse(spearman_p_prev_adj<0.05, "significant", "non-significant"))

#levels for genera (plotting)
genera_corr_sorted_prev<-df_cor_cohorts2_prev %>% arrange(spearman_prev) %>% pull(Genus)


#plot correlations
corgen_prev<-df_cor_cohorts2_prev %>%
  pivot_longer(cols = -c(Genus,significance)) %>%
  filter(name=="spearman_prev") %>%
  #filter(!is.na(value) & name=="spearman") %>%
  mutate(Genus=factor(Genus,levels=genera_corr_sorted_prev)) %>%
  ggplot(aes(x=value, y = Genus, fill=significance))+
  geom_col(position = "dodge")+
  geom_vline(xintercept =c(-0.7,0.7),linetype="dotted",color="#6B96B6",size=1.2)+
  xlim(-1,1)+
  theme_light()+
  scale_fill_manual(values=c("lightgrey","#8EB9B7"))
corgen_prev
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_2b_genera_prevalencecorrelations_y1_fixed.pdf", sep=""),width=8,height=6)
print(corgen_prev)
dev.off()

##########################################
### 5. Differential abundance analysis ###
### Publication *Supplementary Figure S8a*
##########################################
daa_df<-merged_gen_commontoplot_MELT %>%
  mutate(Abundance_adj = ifelse(Abundance <= 0.005, 0, Abundance)) %>%
  select(Abundance_adj, Cohort, Genus)


## Wilcoxon test for comparison
wilcox_results <- daa_df %>%
  group_by(Genus) %>%
  rstatix::wilcox_test(Abundance_adj ~ Cohort) %>% 
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p.adj) %>%
  dplyr::mutate(sig=ifelse(p.adj<0.05,"significant","non-significant"))

#write.csv(wilcox_results,"daa_wilcoxonresults.csv")

# Subset genera that have a significant difference
sig_genera <- wilcox_results %>% filter(p.adj < 0.05) %>% pull(Genus)
sig_genera


ordered_names <- wilcox_results$Genus
daa_df$Genus<- factor(daa_df$Genus, levels = ordered_names)

# Plot the significantly different genera
cohort_plot <- ggplot(daa_df %>% 
                        filter(!Genus %in% c("Ligilactobacillus","Limosilactobacillus",
                                             "Segatella","Parolsenella",
                                             "Sutterella","Parabacteroides", "Thomasclavelia")) %>%
                        filter(Genus %in% sig_genera), 
                      aes(Genus, Abundance_adj, fill = Cohort)) +
  geom_boxplot(lwd = 0.2, outlier.color = "darkgrey", outlier.size = 0.7, outlier.alpha = 0.3) +
  labs(x = "", y = "Relative Abundance", fill = "Cohort") +
  theme_light() +
  scale_y_continuous(trans='log10')+
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, face = "italic", size = 10),
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "darkgrey"),
    axis.ticks = element_line(size = 0.5, color = "darkgrey"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title.y = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    plot.caption = element_text(size = 10, hjust = 1)
  ) +
  scale_fill_manual(values = c("#8FC3D4" ,"#778899"))

cohort_plot
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_S8a_DAA_significantgenera.pdf", sep=""),
    width=12,height=8)
print(cohort_plot)
dev.off()







