############################################################################
## PCoA & dbRDA ANALYSIS — BBGUT COHORT
## ======================================
## - Genus-level PCoA for 3 GMM stages
## - Distance-based Redundancy Analysis on metadata covariates
## - PCoA of infant and maternal samples together
##
## Output figures & tables:
##   Figure 1a  |  Figure 3a  |  Supplementary Table S6
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
library(microViz)
library(glue)
library(tibble)
library(ggplot2)
library(microViz)

# ---- COLOR PALETTES ----
pastel_colors_2 <- c("#BFE1ED", "#7FAFBF", "#778899")
pastel_colors_cat <- c("#D3A9B9", "#E4CEC2", "#A6CFCF", "#75AEBD")

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")

########################################
#Data processing and analysis
########################################
bbgut_ps_inf <- subset_samples(bbgut_ps, fecalsample_type == "infant") 
bbgut_ps_abundance_inf <- transform_sample_counts(bbgut_ps_inf, 
                                                       function(OTU) OTU/sum(OTU))

# match cluster to sample sequencing ID
dmm_cluster_stats <- read.csv("Data/Per_sample_stats.csv") %>% mutate(
  Sample_sequencing_ID=X) %>% select(Sample_sequencing_ID,Cluster) %>% 
  filter(Sample_sequencing_ID %in% 
           sample_data(bbgut_ps_abundance_inf)$Sample_sequencing_ID)

dmm_cluster_stats<-dmm_cluster_stats[match(sample_data(bbgut_ps_abundance_inf)$Sample_sequencing_ID,
                        dmm_cluster_stats$Sample_sequencing_ID),]

sample_data(bbgut_ps_abundance_inf)$Cluster <- dmm_cluster_stats$Cluster

# set up variables we want to include in the dbRDA
EFFECTS_simple_16S <- c('Child_ID', 'Day', 'Siblings','Delivery_mode', 
                        'Delivery_place','Gestation','hyg_level','Antibiotic_use',
                        'Solid_food', 'GI_disease', "Diarrhea","immu_day",
                        "V_Pentavalent","V_Polio_oral","V_Pneumonococcal",
                        "Fever","Other_med_use","other_meds_totalcount","A03",
                        "A07","A11","A12","N02","R01","R03","R05","R06",
                        "Cephalosporin_t","Macrolide","Penicillin","Nitroimidazole")

effects_categories <- c('Infant',  'Infant', 'Infant','Infant',
                        'Infant', 'Infant', 'Infant', 'Antibiotics',
                        'Diet',  'GI_disease','Diarrhea','Vaccination',
                        'Vaccination','Vaccination','Vaccination','Fever',
                        'Other_medication','Other_medication','Other_medication',
                        'Other_medication','Other_medication','Other_medication',
                        'Other_medication','Other_medication','Other_medication',
                        'Other_medication','Other_medication','Antibiotics',
                        'Antibiotics','Antibiotics','Antibiotics')


EFFECTS_Categories_simple_16S <- data.frame(Effect=EFFECTS_simple_16S, 
                                            Category=effects_categories)

RANDOMEFFECT <- c('Child_ID')
FIXEDEFFECTS_simple_16S <- setdiff(EFFECTS_simple_16S,RANDOMEFFECT)
VARIABLES_simple_16S <- c(FIXEDEFFECTS_simple_16S,RANDOMEFFECT)


# (1)PCoA analysis
# (1.1) Genus agglomerate phyloseq object without the mothers
bbgut_inf_agglom <- tax_glom(bbgut_ps_abundance_inf, 
                                  taxrank = 'Genus')

bbgut_inf_agglom_OTU <- otu_table(bbgut_inf_agglom)

# (1.2) PCoA ordination with Bray-Curtis distance
distance = "bray"
set.seed(123)
allinfant_pcoa <- ordinate(physeq=bbgut_inf_agglom,
                           method="PCoA",
                           distance = distance)

# (2) Set up dbRDA analysis

all_effects <- c(FIXEDEFFECTS_simple_16S,"Child_ID")
covariates_metadata <- sample_data(bbgut_inf_agglom) %>% 
  data.frame() %>% 
  select(matches(all_effects)) %>% 
  # Replace NA with No in immu_day
  mutate(immu_day= ifelse(is.na(immu_day), "No", immu_day))


# (2.1) Univariable dbRDA on every metadata variable
# function for univariable dbrda on every metadata var
# otu_genus: Genus agglomerated OTU table
# meta: dataframe with all metadata variables
# d: distance metric
univar_dbrda <- function(otu_genus,meta,d) {
  set.seed(4)
  all_vars <- c() #empty vector where you store results
  for (var in colnames(meta)) { # for each variable
    print(paste("capscale: otu_table ~ ",var))
    capsc <- capscale(otu_genus ~ meta[,var],
                      distance = d, na.action=na.omit)
    an <- anova.cca(capsc) #anova of the fit
    pval <- an["Pr(>F)"][[1]][[1]] #get p-value
    Fa <- an["F"][[1]][[1]] #get F score
    r2 <- RsquareAdj(capsc)[[1]] #get R2 (effect size)
    r2adj <- RsquareAdj(capsc)[[2]] #get adjusted R2
    all_vars <- rbind(all_vars,cbind(Fa,r2,r2adj,pval))
  }
  FDR = p.adjust(all_vars[,"pval"],method="BH")
  all_vars = cbind(all_vars,FDR)
  colnames(all_vars) <- c("F","r2","r2adj","p-value","FDR") #generate table
  row.names(all_vars) <- colnames(meta)
  #sort by r2
  all_vars <- as.data.frame(all_vars) %>% 
    arrange(desc(all_vars))
  return(all_vars)
}

# match order of samples in otu_table to metadata
bbgut_inf_agglom_OTU <- bbgut_inf_agglom_OTU[
                              match(rownames(covariates_metadata), 
                             rownames(bbgut_inf_agglom_OTU)),]

# Run univariable analysis
univar_capsc <- univar_dbrda(bbgut_inf_agglom_OTU,
                             covariates_metadata, distance)

#significant variables
fdr = 0.05
sig.vars <- univar_capsc %>% 
  filter(FDR < fdr & `p-value` < 0.05) %>%
  rownames()

# (2.2) Run OrdiR2step with the significant variables

# Helper function multivariable stepwise forward selected dbRDA
# otu_genus: genums agglomerated OTU table
# meta: metadata of significant variables
# dist: distance metric
multivar_dbrda <- function(otu_genus,meta){
  set.seed(9)
  # Base model
  print('  starting base model')
  m0 <- capscale(otu_genus ~ 1, data=meta ,distance='bray')
  # model with all significant vars
  print('  starting full model')
  m1 <- capscale(otu_genus ~ ., data=meta, distance='bray')
  m1_r2adj <- RsquareAdj(m1)$adj.r.squared
  # stepwise forward var selection
  print("  starting stepwise")
  step.res <- ordiR2step(m0, scope=formula(m1), data=meta,
                         direction="forward", Pin = 0.05,
                         R2scope = m1_r2adj, pstep = 1000, 
                         trace = F)
  print("done step.res")
  print(glue("correcting for multiple testing (n={ncol(meta)} variables)"))
  step.res$anova$adj.pval <- p.adjust (step.res$anova$`Pr(>F)`, 
                                       method = 'BH', n = ncol (meta))
  print(step.res$anova)
  if (length(step.res$anova) > 1){
    ordiR2step.tab <- as.data.frame(step.res$anova)
    rownames(ordiR2step.tab) <-gsub("\\+ ","",rownames(ordiR2step.tab))
    ordiR2step.tab <- independent_stepwise_r2(ordiR2step.tab)
    ordiR2step.tab$Effect <- rownames(ordiR2step.tab)
  } else {
    print("COULD NOT RUN STEPWISE dbRDA")
    ordiR2step.tab <- NULL
  }
  return(ordiR2step.tab)
}

# Helper function calculate independent contributions in stepwise dbRDA
# anova_table: anova table output from ordiR2step
independent_stepwise_r2 <- function(anova_table){
  anova_table <- anova_table[order(anova_table$R2.adj),]
  lagged_R2 <- dplyr::lag(anova_table$R2.adj) %>% replace_na(0)
  anova_table$R2.adj_independent <- anova_table$R2.adj - lagged_R2
  return(anova_table)
}

# Get metadata for significant variables
sig_var_meta <- get_variable(bbgut_inf_agglom,union(sig.vars,"Child_ID"))
sig_var_meta <- na.exclude(sig_var_meta) 
sig_var_meta <- sig_var_meta[match(rownames(bbgut_inf_agglom_OTU),
                                   rownames(sig_var_meta)),]

# Run multivariable analysis
multivar_capsc <- multivar_dbrda(bbgut_inf_agglom_OTU,
                                 sig_var_meta)
# merge results
univar_sig <- univar_capsc %>% 
  filter(FDR < fdr & `p-value` < 0.05) %>% 
  select(r2adj,`p-value`,FDR)
colnames(univar_sig) <- paste(colnames(univar_sig),'.univar',sep='')
univar_sig <- tibble::rownames_to_column(univar_sig, "Effect") 
combined_results <- dplyr::left_join(multivar_capsc,univar_sig,by="Effect") %>%
                    filter(adj.pval < 0.05)

#pivot result for plot
combined_results_pivot <- combined_results %>% select(Effect,R2.adj,r2adj.univar) %>%
  arrange(R2.adj) %>%
  pivot_longer(-Effect,names_to="var") %>% 
  mutate(var=ifelse(var=="r2adj.univar",
                         "Capscale","dbRDA")) %>%
  # add effect category
  left_join(x=.,y=EFFECTS_Categories_simple_16S,by='Effect')

effect_lvls<-combined_results_pivot$Effect %>% unique()
combined_results_pivot$var <- factor(combined_results_pivot$var,levels=c('dbRDA','Capscale'))
combined_results_pivot$Effect <- factor(combined_results_pivot$Effect,
                                        levels=effect_lvls)


## PLOT covariate effect sizes plot
p_cov<-ggplot(combined_results_pivot,aes(x=Effect,y=value))+
  geom_bar(position = "dodge", stat = "identity", 
      aes(fill = Category, group=var, alpha=var)) +
  scale_alpha_manual(values=c(0.55,1)) +
  scale_fill_manual(values=pastel_colors_cat)+  
  theme_light()+  labs(x = "Covariate", y="Adjusted R2")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_cov
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Pcoa_covariate_effect_sizes.pdf", sep=""),width=10,height=3)
print(p_cov)
dev.off()

## plotting processing PCoA Analysis with covariates

# ARROWS FOR PCOA COVARIATES

padj_sig_var <- combined_results_pivot$Effect %>% unique() %>% 
                as.character()

padj_sig_var_meta <- sig_var_meta[,padj_sig_var]
# make ordering match with ordination
padj_sig_var_meta <- padj_sig_var_meta[match(rownames(allinfant_pcoa$vectors),
                        rownames(padj_sig_var_meta)),]
# replace yes/no with 1/0
padj_sig_var_meta[padj_sig_var_meta=="Yes"] <- 1
padj_sig_var_meta[padj_sig_var_meta=="No"] <- 0

padj_sig_var_meta$Day <- padj_sig_var_meta$Day %>% as.numeric()
padj_sig_var_meta$Solid_food <- padj_sig_var_meta$Solid_food %>% as.numeric()
padj_sig_var_meta$Diarrhea <- padj_sig_var_meta$Diarrhea %>% as.numeric()
padj_sig_var_meta$Penicillin <- padj_sig_var_meta$Penicillin %>% as.numeric()
padj_sig_var_meta$Antibiotic_use <- padj_sig_var_meta$Antibiotic_use %>% as.numeric()

envfit_result <- envfit(allinfant_pcoa$vectors,
                    env = padj_sig_var_meta, perm = 10000, 
                    na.rm=T)

# continuous vars -> vector
en_coord_cont <- as.data.frame(scores(envfit_result, "vectors")) * ordiArrowMul(envfit_result)
en_coord_cont$Effect <- rownames(en_coord_cont) 

##################################################################
# PCoA plot with arrows
# Publication *Figure 3a*
##################################################################

#plot: adjust axis limits and arrow size for better visualization
p2 <- plot_ordination(physeq = bbgut_inf_agglom,
                     ordination = allinfant_pcoa,
                     color = "Cluster") +
  scale_color_manual(values = pastel_colors_2) +
  theme_light() +
  labs(color = "Cluster") +
  geom_segment(data = en_coord_cont,
               aes(x = 0, xend = Axis.1 * 0.15,   # Scale down xend
                   y = 0, yend = Axis.2 * 0.15),  # Scale down yend
               arrow = arrow(length = unit(0.08, 'cm')),
               colour = 'gray35') +
  geom_text(data = en_coord_cont,
            aes(x = Axis.1 * 0.17,  # Adjust text positions accordingly
                y = Axis.2 * 0.17, 
                label = Effect),
            color = "gray35", size = 4) +
  coord_fixed(ratio = 0.85) +
  scale_y_reverse() +
  xlim(-0.5, 0.7)
p2
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_3a_pcoa_dbRDA_correctedsize.pdf", sep=""),width=10,height=6)
print(p2)
dev.off()

##################################################################
# (3) PcoA with to visualize infant and maternal samples together
# Publication *Figure 1a*
##################################################################

#include bbgut_ps_abundance(has all samples)
bbgut_ps_abundance<-transform_sample_counts((bbgut_ps), 
                        function(OTU) OTU/sum(OTU))
bbgut_ps_Genus_infmat<-tax_glom(bbgut_ps_abundance, taxrank="Genus")
distance = "bray"

allinfant_pcoa_mi <- ordinate(physeq=bbgut_ps_Genus_infmat,method="PCoA",distance = distance)

# match cluster to sample sequencing ID for all samples
dmm_cluster_stats_allsamples <- read.csv("Data/Per_sample_stats.csv") %>% mutate(
  Sample_sequencing_ID=X) %>% select(Sample_sequencing_ID,Cluster) %>%
  filter(Sample_sequencing_ID %in%
            sample_data(bbgut_ps_Genus_infmat)$Sample_sequencing_ID)

dmm_cluster_stats_allsamples<-dmm_cluster_stats_allsamples[match(sample_data(bbgut_ps_Genus_infmat)$Sample_sequencing_ID,                                                     dmm_cluster_stats_allsamples$Sample_sequencing_ID),]

sample_data(bbgut_ps_Genus_infmat)$Cluster <- dmm_cluster_stats_allsamples$Cluster


SHAPES <- c('infant' = 16,'maternal' = 24)

p_mi <- plot_ordination(physeq=bbgut_ps_Genus_infmat, 
                        ordination = allinfant_pcoa_mi, color="Cluster", shape = "fecalsample_type")+
  ggtitle(paste("PCoA ,",distance,", All Samples",  sep=""))+
  scale_color_manual(values= as.character(pastel_colors_2))+
  theme_light()+
  scale_shape_manual(values = SHAPES)+
  geom_point(size=2)+
  labs( color ="Cluster")+
  coord_fixed()  ## need aspect ratio of 1!
p_mi
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_1a_pcoa_mother_infants.pdf", sep=""),width=7,height=7)
print(p_mi)
dev.off()
