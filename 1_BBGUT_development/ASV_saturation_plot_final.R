############################################################################
## ASV SATURATION PLOT WITH INCREASING SAMPLE SIZE
## =================================================
## Shows how ASV richness increases with the number of infants sampled,
## comparing BBGUT alone and BBGUT vs BABEL cohorts.
##
## Output figures:
##   Supplementary Figure S1a  |  Supplementary Figure S1b
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(tidyverse)
library(dplyr)
library(phyloseq)
library(tidyr)
library(microbiome)
library(ggplot2)
library(microViz)

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")
load("Data/babel_phyloseq_final_newtaxonomy.RData")

# ---- COLORS ----
cohort_cols <- c("#8FC3D4", "#778899")

### melt phyloseq object ###
#bbgut
bbgut_ps_healthy<-bbgut_ps %>%
  ps_filter(fecalsample_type=="infant") %>%
  ps_filter(sample_category=="HTP")

bbgut_ps_healthy_y1<-bbgut_ps %>%
  ps_filter(fecalsample_type=="infant") %>%
  ps_filter(sample_category=="HTP") %>%
  ps_filter(Day <= 365)
  
bbgut_ps_melt<- ps_melt(bbgut_ps_healthy)
bbgut_ps_melty1<- ps_melt(bbgut_ps_healthy_y1)
  

infants<-bbgut_ps_melt$Child_ID %>% unique()

#babel
babel_ps_healthy_y1<- babel_ps %>%
  ps_filter(LDA==1) %>%
  ps_filter(InfantID !="S011") %>% #remove S011 (not considered in this study)
  ps_filter(X.days <= 365)

babel_ps_melt<- ps_melt(babel_ps_healthy_y1)
babel_ps_melt$Child_ID<-babel_ps_melt$InfantID
infants_babel<-babel_ps_melt$Child_ID %>% unique()


#make function that gets the unique number of asvs per random number of infants
get_unique_asv<- function(df,inflist,infnum){
  unique_asvs<-df %>% 
    filter(Child_ID %in% sample(inflist,infnum) & Abundance > 0) %>% 
    pull(OTU) %>% 
    unique() %>% length()
  return(unique_asvs)
}

##################################################################  
## BBGUT & BABEL comparison year 1 & HTP samples #############
## Publication *Supplementary Figure 1b*
##################################################################  
#create empty variable
plot_asvsatur<-c()

#loop through 20 bbgut infants using 7 iterations (=babel max number of infants)
set.seed(10)
for (n in 1:20) {
  for (i in 1:7){
    num_inf<-n
    repl<-i
    num_asv_bbgut<-get_unique_asv(bbgut_ps_melty1,infants,n)
    bbgut_row<-c("numinf"=num_inf,"numASV"=num_asv_bbgut,"cohort"="BBGUT","rep"=repl)
    plot_asvsatur<-rbind(plot_asvsatur,bbgut_row)
  }
}

#loop through 7 babel infants
set.seed(10)
for (n in 1:7) {
  for (i in 1:7){
    num_inf<-n
    repl<-i
    num_asv_babel<-get_unique_asv(babel_ps_melt,infants_babel,n)
    babel_row<-c("numinf"=num_inf,"numASV"=num_asv_babel,"cohort"="BABEL","rep"=repl)
    plot_asvsatur<-rbind(plot_asvsatur,babel_row)
  }
}

#make df
plot_asvsatur_df<-as.data.frame(plot_asvsatur)
#determine levels for plotting
plot_asvsatur_df$numinf<-factor(plot_asvsatur_df$numinf,levels=(1:20))

#boxplot
asv<-plot_asvsatur_df %>% ggplot(aes(x=numinf, y=as.numeric(numASV), fill=cohort)) +
  geom_boxplot(outlier.color="darkgrey")+
  scale_fill_manual(values= cohort_cols)+
  theme_light()+
  ylab("ASV richness")+
  xlab("Number of Infants")
asv
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_1b_ASV_saturationplot_HEALTHY_Y1.pdf", sep=""),
    width=10,height=5)
print(asv)
dev.off()

################################################ 
## BBGUT HTP samples 0-2 years only        
## Publication *Supplementary Figure 1a*
################################################
#create empty variable
plot_asvsaturbbgut<-c()

#loop through 20 bbgut infants using 20 iterations (=bbgut max number of infants)
set.seed(10)
for (n in 1:20) {
  for (i in 1:20){
    num_inf2<-n
    repl2<-i
    num_asv_bbgut2<-get_unique_asv(bbgut_ps_melt,infants,n)
    bbgut_row2<-c("numinf"=num_inf2,"numASV"=num_asv_bbgut2,"cohort"="BBGUT","rep"=repl2)
    plot_asvsaturbbgut<-rbind(plot_asvsaturbbgut,bbgut_row2)
  }
}

#make df
plot_asvsaturbbgut_df<-as.data.frame(plot_asvsaturbbgut)
#determine levels for plotting
plot_asvsaturbbgut_df$numinf<-factor(plot_asvsaturbbgut_df$numinf,levels=(1:20))

#boxplot
asvB<-plot_asvsaturbbgut_df %>% ggplot(aes(x=numinf, y=as.numeric(numASV), fill=cohort)) +
  geom_boxplot(outlier.color="darkgrey")+
  scale_fill_manual(values= "#778899")+
  theme_light()+
  ylab("ASV richness")+
  xlab("Number of Infants")
asvB
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Supplementary_Figure_1a_ASV_saturationplot_BBGUT_HEALTHY.pdf", sep=""),
    width=10,height=5)
print(asvB)
dev.off()
 