############################################################################
## NORMALIZED 2-MONTH BACTERIAL ABUNDANCE CHANGES ACROSS PERIODS
## ===============================================================
## Linear mixed-effects model comparing E-E, E-L, L-L period transitions.
## Provides statistics for Supplementary Figure S13.
##
## Output:
##   Supplementary Table S9
############################################################################

rm(list = ls())
set.seed(2)

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(dplyr)
library(lme4)
library(lsmeans)

# ---- LOAD DATA ----
data <- openxlsx::read.xlsx("Data/differences_abundance_allpairs_2months.xlsx")
#normalize abundance differences
data<-data %>%
  mutate(
    NormDiff = (Abundance_median_2 - Abundance_median_1) / ((Abundance_median_1 + Abundance_median_2)/2))
#calculate extra age parameter
data$age <- (data$Month_2+data$Month_1)/2

genera <- unique(data$Genus)

outcomelist <- vector("list", length = length(genera))
outcomedf <- data.frame()

for (i in 1:length(genera)){
  genus <- genera[i]
  testdata <- data[which(data$Genus == genus),]
  testdata <- testdata[which(!is.na(testdata$NormDiff)), ]
  rownames(testdata) <- 1:nrow(testdata)
  mod <- lmer(
    data = testdata,
    formula = "NormDiff ~ Period_Comparison + age + (1|Child_ID)",
    REML = FALSE
  )
  modsum <- summary(mod)
  
  ################################
  ### Predictions/Pairwise tests # getting pvalue
  ################################
  
  
  
  #observed
  means <- lsmeans(mod, "Period_Comparison")
  #means #should look like the plot
  contrasts <- contrast(
    means, 
    list(
      #specify contrasts according to your formula
      `L-L vs E-L`     = c(-1, 1, 0), #-post-post, +pre-post, 0pre-pre --> pre-post - post-post (always second one - first one)
      `L-L vs E-E`     = c(-1, 0, 1),
      `E-L vs E-E`     = c(0, -1, 1)
    )
  )
  #contrasts
  
  sublist <- list(
    inferred_means = summary(means), #means of groups if you want to check simplified contrasts
    contrasts  = summary(contrasts)
  )
  outcomelist[[i]] <- sublist
  names(outcomelist)[i] <- genus
  
  consum <- summary(contrasts)
  tempdf <- data.frame(
    genus = genus,
    comp = consum$contrast,
    diff = consum$estimate,
    pval = consum$p.value
  )
  outcomedf <- rbind(outcomedf, tempdf)
}
  
outcomedf$sig <- ifelse(outcomedf$pval<0.05, "Significant", "")
#openxlsx::write.xlsx(outcomedf, "outcome_mixedmodel_2monthchanges.xlsx")

library(purrr)
out_all <- map_dfr(names(outcomelist), function(g) {
  outcomelist[[g]][["contrasts"]] %>%
    mutate(genus = g)
})
out_all$sig <- ifelse(out_all$p.value<0.05, "Significant", "")
#openxlsx::write.xlsx(out_all, "outcome_mixedmodel_2monthchanges_allinfo.xlsx")

