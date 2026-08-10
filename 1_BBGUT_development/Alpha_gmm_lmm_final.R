############################################################################
## LINEAR MIXED EFFECTS MODEL — ALPHA DIVERSITY ACROSS GMMs
## =========================================================
## Tests whether alpha diversity (observed richness & Shannon) differs
## significantly between maturation clusters (A, B, C, C_m).
##
## Output:
##   Supplementary Table S4
############################################################################

rm(list = ls())
set.seed(2)

# ---- WORKING DIRECTORY ----
## Set this to the root project folder
workdir <- ""
setwd(workdir)

# ---- LOAD PACKAGES ----
library(ggplot2)
library(lme4)
library(lsmeans)

# ---- LOAD DATA ----
load("Data/alpha_diversity_bbgut.RData")

data <- alpha_fulldataset
rm(alpha_fulldataset)
gc()

data$ID<-unlist(lapply(data$Sample_sequencing_ID,function(x) strsplit(x,"_")[[1]][1]))


################################
### Make the model
################################

data$Cluster <- factor(data$Cluster, levels = c("A", "B", "C","C_m"))


model_formula <- "observed ~ Cluster + (1 | ID)"
model_formula_shannon<-"diversity_shannon ~ Cluster + (1 | ID)"

#observed
model <- lmer(
  data = data,
  formula = model_formula,
  REML = FALSE #restricted maximum likelihood
)

modsum <- summary(model)
model_coefficients <- modsum[["coefficients"]]
modsum

#shannon
model_shannon <- lmer(
  data = data,
  formula = model_formula_shannon,
  REML = FALSE
)

modsum_shannon <- summary(model_shannon)
model_coefficients_shannon <- modsum_shannon[["coefficients"]]
modsum_shannon

################################
### Some assumption checking 
# how much do your values differ from the mean?
################################

resplot_pdata <- data.frame(
  diversity = data$diversity_shannon,
  fitted = fitted(model_shannon),
  residual = modsum_shannon$residuals
)

#heteroschedasticity check = you want to see that all your residuals have same variances
ggplot() +
  geom_point(
    data = resplot_pdata,
    aes(x = fitted, y = residual)
  ) + 
  geom_hline(yintercept = 0, col = "grey") +
  ggtitle("Residuals vs Fitted") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = .5)
  )

#normality check histogram
ggplot() +
  geom_histogram(
    data = resplot_pdata,
    aes(x = residual),
    bins = 10
  ) +
  ggtitle("Residual histogram") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(hjust = .5)
  )


################################
### Significance testing 
################################

#observed reduced formula
reduced_formula <- "observed ~ (1 | ID)"

reduced_model <- lmer(
  data = data,
  formula = reduced_formula,
  REML = FALSE
)

test <- anova(model, reduced_model)

test


#shannon reduced formula
reduced_formula_shannon <- "diversity_shannon ~ (1 | ID)"

reduced_model_shannon <- lmer(
  data = data,
  formula = reduced_formula_shannon,
  REML = FALSE
)

test_shannon <- anova(model_shannon, reduced_model_shannon)

test_shannon

################################
### Predictions/Pairwise tests # getting pvalue
################################

#observed
means <- lsmeans(model, "Cluster")
contrasts <- contrast(
  means, 
  list(
    AvB = c(1,-1,0,0),
    AvC = c(1, 0, -1,0),
    AvC_m = c(1,0,0,-1),
    BvC = c(0, 1, -1,0),
    BvC_m = c(0, 1,0, -1),
    CvC_m = c(0, 0,1, -1)
  )
)
contrasts

#shannon
means_shannon <- lsmeans(model_shannon, "Cluster")
contrasts_shannon <- contrast(
  means_shannon, 
  list(
    AvB = c(1,-1,0,0),
    AvC = c(1, 0, -1,0),
    BvC = c(0, 1, -1,0),
    AvC_m = c(1,0,0,-1),
    BvC_m = c(0, 1,0, -1),
    CvC_m = c(0, 0,1, -1)
  )
)
contrasts_shannon



