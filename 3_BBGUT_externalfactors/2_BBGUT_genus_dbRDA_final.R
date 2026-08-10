############################################################################
## GENUS-LEVEL dbRDA — COVARIATES EXPLAINING BACTERIAL VARIATION
## ===============================================================
## Identifies which metadata covariates explain variation in the most
## abundant bacterial genera using dbRDA analysis.
##
## Output figures & tables:
##   Figure 3b  |  Supplementary Table S7
############################################################################

# ---- WORKING DIRECTORY ----
## Set this to the root project folder (the folder containing `Data/`)
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
pastel_colors <- c("#6B96B6", "#8FC3D4", "#C4E3F0", "#F8BDAF", "#FFD9AC",
                   "#8EB9B7", "#AEC89E", "#CBBE76", "#F5DDB7", "#75AEBD",
                   "#A6CFCF", "#E78F8D", "#FFC1BE", "#A39F96", "#D3A9B9",
                   "#FED4E3", "#C694B3", "#E7C6D5", "#BAAFA5", "#E4CEC2",
                   "lightgrey")
pastel_colors_cat <- c("#75AEBD", "#E4CEC2", "#A6CFCF")

# ---- LOAD DATA ----
load("Data/bbgut_phyloseq_final.RData")


#metadata to test
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

FIXEDEFFECTS_simple_16S_categorical <- setdiff(FIXEDEFFECTS_simple_16S,"Day")
FIXEDEFFECTS_simple_16S_ctu <- "Day"


#Data processing

bbgut_ps_inf <- subset_samples(bbgut_ps, fecalsample_type == "infant") 


#define healthy samples
HTP_samples<-sample_data(bbgut_ps)$Sample_sequencing_ID[sample_data(bbgut_ps)$sample_category == "HTP"] %>%na.omit()

bbgut_ps_abundance_filtered <- transform_sample_counts(bbgut_ps_inf, 
                                    function(OTU) OTU/sum(OTU))

bbgut_ps_Genus <- tax_glom(bbgut_ps_abundance_filtered, taxrank="Genus")


#top 15 most abundant genera in bbgut cohort (years 1-2)
top15_genera <- bbgut_ps_Genus %>% 
  psmelt() %>% 
  filter(Sample_sequencing_ID %in% HTP_samples) %>%
  group_by(Genus) %>%
  dplyr::summarize(total_abundance=sum(Abundance)) %>%
  # add new column with relative ab
  mutate(relative_ab=round(total_abundance/sum(total_abundance)*100,2)) %>%
  # sort descending
  arrange(desc(total_abundance)) %>% 
  head(15) %>%
  pull(Genus)

######################
## dbRDA analysis ##
######################

#helper functions
# function for univariable dbrda on every metadata var
univar_dbrda <- function(otu_genus,meta,d) {
  all_vars <- c() #empty vector where you store results
  for (var in colnames(meta)) { # for each variable
    print(paste(" ",var))
    set.seed(4)
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

multivar_dbrda <- function(otu_genus,meta){
  # model with Child_ID only
  print('  starting base model')
  mod0 <- capscale(otu_genus ~ Child_ID, data=meta ,distance='euclidean')
  # model with all vars
  print('  starting full model')
  mod1 <- capscale(otu_genus ~ ., data=meta, distance='euclidean')
  set.seed(2)
  # stepwise forward var selection
  print("  starting stepwise")
  step.res <- ordiR2step(mod0, scope=formula(mod1), data=meta,
                         direction="forward", Pin = 0.05,
                         R2scope = FALSE, pstep = 1000, trace = F)
  # BH-adjust p-values from forward selection
  step.res$anova$adj.pval <- p.adjust(
    step.res$anova$`Pr(>F)`,
    method = "BH",
    n = ncol(meta))
  print("done step.res")
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

independent_stepwise_r2 <- function(anova_table){
  anova_table <- anova_table[order(anova_table$R2.adj),]
  lagged_R2 <- dplyr::lag(anova_table$R2.adj) %>% replace_na(0)
  anova_table$R2.adj_independent <- anova_table$R2.adj - lagged_R2
  return(anova_table)
}


distance_metric <- 'euclidean'
fdr = 0.05
vars_fixed <- c(FIXEDEFFECTS_simple_16S,"Child_ID")
bbgut_genus_fixed <- sample_data(bbgut_ps_Genus) %>% 
  data.frame() %>% 
  select(matches(vars_fixed)) %>% 
  # Replace NA with No in immu_day
  mutate(immu_day= ifelse(is.na(immu_day), "No", immu_day))


all_genus_results <- c() # save all genus results here
# for each genus
for (genus in top15_genera){
  print("########################")
  print(genus)
  print("########################")
  ps_subset_genus <- subset_taxa(bbgut_ps_Genus, Genus==genus) # subset ps
  otu_genus <- otu_table(ps_subset_genus) #otu table for genus
  # match order to metadata
  otu_genus <- otu_genus[match(rownames(bbgut_genus_fixed), 
                  rownames(otu_genus)),] 
  # run a univariate analysis using euclidian distances
  print('RUNNING UNIVARIABLE ANALYSIS')
  univar_capsc <- univar_dbrda(otu_genus,bbgut_genus_fixed,distance_metric)
  
  #significant variables
  sig.vars <- univar_capsc %>% 
    filter(FDR < fdr & `p-value` < 0.05) %>%
    rownames()
  
  # if number of significant variables > 0 run OrdiR2step per genus
  if (length(sig.vars) > 0 ){
    sig_var_meta <- get_variable(bbgut_ps_Genus,union(sig.vars,"Child_ID"))
    sig_var_meta <- na.exclude(sig_var_meta) 
    sig_var_meta <- sig_var_meta[match(rownames(otu_genus),
                                rownames(sig_var_meta)),]
    
    print("RUNNING ORDIR2STEP")
    multivar_capsc <- multivar_dbrda(otu_genus,sig_var_meta)
    # if multivar_capsc results are not NULL
    if (!identical(multivar_capsc,NULL)){
      # merge results
      univar_sig <- univar_capsc %>% 
        filter(FDR < fdr & `p-value` < 0.05) %>% 
        select(r2adj,`p-value`,FDR)
      colnames(univar_sig) <- paste(colnames(univar_sig),'.univar',sep='')
      univar_sig <- tibble::rownames_to_column(univar_sig, "Effect") 
      combined_results <- dplyr::left_join(multivar_capsc,univar_sig,by="Effect")%>%
        filter(adj.pval < 0.05)
      combined_results$Genus <- genus
      # 
      all_genus_results <- rbind(all_genus_results,combined_results)
    }
  }
  print("")
}

#### select covariate with highest r2 value per genus ####
main_var_gen<-all_genus_results %>% 
  group_by(Genus) %>% 
  filter(r2adj.univar==max(r2adj.univar)) %>% 
  ungroup() %>%
  arrange(desc(r2adj.univar))

#write.csv(main_var_gen,"dbrda_results_genus.csv")

### plot ###
# Fist convert 'Genus' column to a factor with the desired order
main_var_gen$Genus <- factor(main_var_gen$Genus, levels = unique(main_var_gen$Genus))
main_var_gen<-main_var_gen %>% mutate(Genus = forcats::fct_rev(fct_inorder(Genus)))

cov_gen<-ggplot(main_var_gen, aes(x=Genus, y=r2adj.univar, fill=Effect)) +    
  geom_bar(stat = "identity") +
  scale_fill_manual(values=pastel_colors_cat)+
  theme_light()+
  coord_flip()
cov_gen
pdf(paste(workdir,format(Sys.time(), "%Y-%m-%d"),
          "Figure_3b_Genera_covariates.pdf", sep=""),width=7,height=7)
print(cov_gen)
dev.off()