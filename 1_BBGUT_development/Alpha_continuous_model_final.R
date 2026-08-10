############################################################################
## LINEAR MIXED EFFECTS MODEL — ALPHA DIVERSITY OVER TIME
## =======================================================
## Tests whether alpha diversity (observed richness & Shannon) changes
## significantly with infant age, and compares last trimester infants vs mothers.
##
## Output:
##   Supplementary Table S3
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

# ---- LOAD DATA ----
data <- read.csv("Data/alpha_infants_monthlyAPRIL25.csv")

rm(alpha_infonly)
gc()

data$Child_ID<-unlist(lapply(data$Sample_sequencing_ID,function(x) strsplit(x,"_")[[1]][1]))

#filter data for keeping infants only
data<-data[which(!grepl("IM", data$Child_ID)),]
data$Month <- ifelse(data$Day <= 10,0,data$Month)


################################
### Make the model
################################


model_formula <- "observed ~ Month + (1 | Child_ID)"
model_formula_shannon<-"diversity_shannon ~ Month + (1 | Child_ID)"

#observed
model <- lmer(
  data = data,
  formula = model_formula,
  REML = FALSE #restricted maximum likelihood
)

modsum <- summary(model)
model_coefficients <- modsum[["coefficients"]]
modsum
# Fixed effects:
#             Estimate   Std. Error t value
# (Intercept)  12.2358     0.9731   12.57
# Month         1.6416     0.0494   33.23 
#if Month increases by 1 the diversity increases by 1.6416 

#shannon
model_shannon <- lmer(
  data = data,
  formula = model_formula_shannon,
  REML = FALSE #restricted maximum likelihood
)

modsum_shannon <- summary(model_shannon)
model_coefficients_shannon <- modsum_shannon[["coefficients"]]
modsum_shannon

################################
### Some assumption checking 
################################

resplot_pdata <- data.frame(
  diversity = data$diversity_shannon,
  fitted = fitted(model_shannon),
  residual = modsum_shannon$residuals
)

#heteroschedasticity check
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
reduced_formula <- "observed ~ (1 | Child_ID)"

reduced_model <- lmer(
  data = data,
  formula = reduced_formula,
  REML = FALSE
)

test <- anova(model, reduced_model)

test #pvalue is significant p< 2.2e-16 *** the slope is significant 
#so diversity increases with time (Month)



#shannon reduced formula
reduced_formula_shannon <- "diversity_shannon ~ (1 | Child_ID)"

reduced_model_shannon <- lmer(
  data = data,
  formula = reduced_formula_shannon,
  REML = FALSE
)

test_shannon <- anova(model_shannon, reduced_model_shannon)

test_shannon #same with shannon diversity, p< 2.2e-16 ***


#########################################################################
### Comparing diversity differences between the last 3mester and mothers
##########################################################################

#prepare data
data_full <- read.csv("Data/alpha_infants_monthlyAPRIL25.csv")

#Create Group column based on pattern
data_full$Group <- ifelse(grepl("^IM\\d+", data_full$Sample_sequencing_ID), "mother",
                          ifelse(grepl("^I\\d+", data_full$Sample_sequencing_ID), "infant", NA))

#Create ID column by extracting the number after IM or I and before the first underscore
data_full$ID <- sub("^IM(\\d+)_.*", "\\1", data_full$Sample_sequencing_ID)
data_full$ID <- sub("^I(\\d+)_.*", "\\1", data_full$ID)

#keep infant data from months 21-24 and mother data
data_full <- subset(data_full, 
                    (Group == "mother") | 
                      (Group == "infant" & Month %in% c(21, 22, 23, 24)))

################################
### Make the model
################################


model_formula_mi <- "observed ~ Group + (1 | ID)"
model_formula_shannon_mi<-"diversity_shannon ~ Group + (1 | ID)"

#observed
model_mi <- lmer(
  data = data_full,
  formula = model_formula_mi,
  REML = FALSE #restricted maximum likelihood
)

modsum_mi <- summary(model_mi)
model_coefficients_mi <- modsum_mi[["coefficients"]]
modsum_mi

#shannon
model_shannon_mi <- lmer(
  data = data_full,
  formula = model_formula_shannon_mi,
  REML = FALSE #restricted maximum likelihood
)

modsum_shannon_mi <- summary(model_shannon_mi)
model_coefficients_shannon_mi <- modsum_shannon_mi[["coefficients"]]
modsum_shannon_mi

################################
### Some assumption checking 
################################

resplot_pdata <- data.frame(
  diversity = data_full$observed,
  fitted = fitted(model_mi),
  residual = modsum_mi$residuals
)

#heteroschedasticity check
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
reduced_formula_mi <- "observed ~ (1 | ID)"

reduced_model_mi <- lmer(
  data = data_full,
  formula = reduced_formula_mi,
  REML = FALSE
)

test <- anova(model_mi, reduced_model_mi)

test #pvalue is significant 1.241e-06 *** there is a difference in diversity 
#between the two groups.


#shannon reduced formula
reduced_formula_shannon_mi <- "diversity_shannon ~ (1 | ID)"

reduced_model_shannon_mi <- lmer(
  data = data_full,
  formula = reduced_formula_shannon_mi,
  REML = FALSE
)

test_shannon <- anova(model_shannon_mi, reduced_model_shannon_mi)

test_shannon #with shannon diversity, p=0.09876 so it's not 
#significantly different between the groups

