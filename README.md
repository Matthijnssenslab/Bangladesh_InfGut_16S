# Gut microbiota development and dynamics in a Bangladeshi infant cohort during early life.
Code to perform analyses on the 16S rRNA data of the Bangladeshi infant
gut microbiota project

For questions please contact Maria Ioanna Papadaki (`papadaki.mg@gmail.com`)

This repository contains the analysis scripts for our study on infant
gut microbiota development in a Bangladeshi infant population (BBGUT
cohort) and how it compares to a Belgian infant cohort (BABEL cohort)
during early life. The STORMS checklist for this study is also provided.
Each script is organized to guide you through the
downstream analysis and visualization of results. All main figures and
supplementary data can be reproduced with the following scripts:

## 1_BBGUT_DEVELOPMENT
#### <ins>`1_BBGUT_healthygutprofile_AUG25.R`</ins>

**Function:**
Exploration of the healthy gut microbiota development in Bangladeshi infants:
1. DMM
2. GMMs
3. Alpha diversity
4. Composition
5. Covid-19

**Data generated:**
* Figure 1a,1b,1c,1e
* Figure 4a,b
* Supplementary Figure 5a
* Supplementary Figure 6
* Supplementary Figure 12
* Supplementary Table 1
* Supplementary Table 9

#### <ins>`2_BBGUT_alpha_div_timebinsAUG25.R`</ins>

**Function:**
Alpha diversity and Supplementary figures (BBGUT cohort)

**Data generated:**
* Supplementary Figure 5b
* Supplementary Figure 10

#### <ins>`ASV_saturation_plotAUG25.R`</ins>

**Function:**
ASV saturation plot with increasing sample size (# infants)

**Data generated:**
* Supplementary Figure 1a,b

#### <ins>`Alpha_gmm_lmm.R`</ins>

**Function:**
Alpha diversity across GMMs (observed, shannon) LMM model

**Data generated:**
* Supplementary Table 2


#### <ins>`Alpha_continuous_model.R`</ins>

**Function:**
Alpha diversity over time (observed, shannon) LMM model

**Data generated:**
* Supplementary Table 2

#### <ins>`covidmonthpairs_comparison_lmm.R`</ins>

**Function:**
Normalized 2-month bacterial abundance changes across periods (P-P, P-L, L-L)

**Data generated:**
* Supplementary Table 8


## 2_BBGUTVSBABEL
#### <ins>`1_BBGUTvsBABEL_Y1_dmm_diversity_AUG25.R`</ins>

**Function:** 
BBGUT – BABEL cohort comparison
1. DMM (for BABEL, for BBGUT use DMM results from full dataset)
2. Alpha diversity (comparison)
3. Beta diversity (comparison)

**Data generated:**
* Figure 2c
* Supplementary Table 3

#### <ins>`1_BBGUTvsBABEL_compdifferences_AUG25.R`</ins>

**Function:** 
BBGUT – BABEL comparison
1. Prevalence changes of common bacterial genera over time
2. Alpha diversity (comparison)

**Data generated:**
* Figure 2a
* Supplementary Table 4


#### <ins>`2_BBGUTvsBABEL_order_of_appearance_Y1_AUG25.R`</ins>

**Function:** 
BBGUT – BABEL comparison
1. Top15 most abundant genera per cohort during year 1
2. Order of Appearance of most abundant genera (y1)
3. Rank (order of appearance) correlations between the two cohorts (y1)
4. Bacterial prevalence correlations between the two cohorts
5. Differential abundances


**Data generated:**
* Figure 2b, d
* Supplementary Figure 7a, b
* Supplementary Figure 8a, b


## 3_BBGUT_EXTERNALFACTORS
#### <ins>`1_BBGUT_PCoA_dbRDA_AUG25.R`</ins>

**Function:** 
PCoA, dbRDA analysis in BBGUT cohort
1. genus level PCoA analysis for the 3 GMM stages
2. Distance-based Redundancy Analysis on for metadata covariates
3. PcoA infant and maternal samples


**Data generated:**
* Figure 1a
* Figure 3a
* Supplementary Table 5


#### <ins>`2_BBGUT_genus_dbRDA_AUG25.R`</ins>

**Function:** 
Finding covariates that explain the variation of the most abundant bacterial genera
dbRDA analysis in BBGUT cohort


**Data generated:**
* Figure 3b
* Supplementary Table 6


#### <ins>`3_BBGUT_Setbacks_AUG25.R`</ins>

**Function:** 
Maturation setbacks in the BBGUT cohort
1. Setbacks in BBGUT cohorts (all changes)
2. Calculate maturation score
3. Setback association to disease events


**Data generated:**
* Figure 3c
* Figure 3d
* Figure 4c
* Supplementary Figure 9a, b
* Supplementary Table 7


## 4_SEGATELLA
#### <ins>`1_BBGUT_Segatella_AUG25.R`</ins>

**Function:** 
Exploring Segatella in the BBGUT cohort
1. Segatella prevalence across GMM stages and health groups
2. HSA threshold calculation
3. Monthly proportion of HSA in adhoc disease and HTP samples
4. HSA metadata 
5. Heatmap for samples with HSA spikes per child over time


**Data generated:**
* Figure 5 a, b, c
* Supplementary Figure 4a, b
* Supplementary Figure 13a, b
* Supplementary Table 10


#### <ins>`2_BBGUT_Segatella_Supplementaryfigure14.R`</ins>

**Function:** 
Segatella ASV abundance across health groups


**Data generated:**
* Supplementary Figure 14



## DATA

Folder containing the data for running all scripts of the main analysis



## DATA_PREPROCESSING

Scripts used to create the final ASV table using the pre-processed raw data


#### <ins>`1_BBGUT_rdp_annotation.R`</ins>

**Function:** 
OTU Annotation with RDP taxonomy (rdp set 19)


#### <ins>`2_BBGUT_decontamination.R`</ins>

**Function:** 
Decontamination using decontam R package

#### <ins>`3_BBGUT_final_phyloseq.R`</ins>

**Function:** 
Create final phyloseq object


