# **Global and local trajectories of infant gut microbiota development in Bangladesh across the COVID-19 pandemic.**
Code to perform analyses on the 16S rRNA data of the Bangladeshi infant gut microbiota project.

**For questions please contact:** Maria Ioanna Papadaki (papadaki.mg@gmail.com)

This repository contains the analysis scripts for our study on infant gut microbiota development in a Bangladeshi infant population (**BBGUT cohort**) and how it compares to a Belgian infant cohort (**BABEL cohort**) during early life. Each script is organized to guide you through the downstream analysis and visualization of results. All main figures and supplementary data can be reproduced with the following scripts:

---

## **1_BBGUT_development** — Healthy Gut Microbiota Development (BBGUT)

| Script | Function | DATA generated |
|--------|----------|----------------|
| *1_BBGUT_healthygutprofile_final.R* | The healthy gut microbiota development in Bangladeshi infants:<br>1. DMM<br>2. GMMs<br>3. Alpha diversity<br>4. Composition<br>5. Covid-19 | Figure 1b, 1c, 1e<br>Figure 4a, b<br>Supplementary Figure S6a<br>Supplementary Figure S7<br>Supplementary Figure S13<br>Supplementary Table S2<br>Supplementary Table S10 |
| *2_BBGUT_alpha_div_timebins_final.R* | Alpha diversity Supplementary figures (BBGUT cohort) | Supplementary Figure S6b<br>Supplementary Figure S12 |
| *ASV_saturation_plot_final.R* | ASV saturation plot with increasing sample size (# infants) | Supplementary Figure S1a, b |
| *Alpha_gmm_lmm_final.R* | Alpha diversity across GMMs (observed, shannon) LMM model | Supplementary Table S3 |
| *Alpha_continuous_model_final.R* | Alpha diversity over time (observed, shannon) LMM model | Supplementary Table S3 |
| *covidmonthpairs_comparison_lmm_final.R* | Normalized 2-month bacterial abundance changes across periods (P-P, P-L, L-L) | Supplementary Table S9 |

---

## **2_BBGUTvsBABEL** — BBGUT vs. BABEL Cohort Comparison

| Script | Function | DATA generated |
|--------|----------|----------------|
| *1_BBGUTvsBABEL_Y1_dmm_diversity_final.R* | BBGUT – BABEL cohort comparison:<br>3. Alpha diversity (comparison)<br>4. Alpha diversity by month (comparison)<br>5. Beta diversity (comparison) | Figure 2c<br>Supplementary Figure S9<br>Supplementary Table S4, S5 |
| *1_BBGUTvsBABEL_compdifferences_final.R* | BBGUT – BABEL comparison | Figure 2a |
| *2_BBGUTvsBABEL_order_of_appearance_Y1_final.R* | 1. Top 15 most abundant genera per cohort during year 1<br>2. Order of Appearance of most abundant genera (y1)<br>3. Rank (order of appearance) correlations between the two cohorts (y1)<br>4. Bacterial prevalence correlations between the two cohorts<br>5. Differential abundances | Figure 2b, d<br>Supplementary Figure S10a, b<br>Supplementary Figure S8a, b |

---

## **3_BBGUT_externalfactors** — External Factors & Metadata Covariates

| Script | Function | DATA generated |
|--------|----------|----------------|
| *1_BBGUT_PCoA_dbRDA_final.R* | PCoA, dbRDA analysis in BBGUT cohort:<br>2. Distance-based Redundancy Analysis for metadata covariates<br>3. PCoA infant and maternal samples | Figure 1a<br>Figure 3a<br>Supplementary Table S6 |
| *2_BBGUT_genus_dbRDA_final.R* | Covariates that explain the variation of the most abundant bacterial genera — dbRDA analysis in BBGUT cohort | Figure 3b<br>Supplementary Table S7 |
| *3_BBGUT_Setbacks_final.R* | Maturation setbacks in the BBGUT cohort:<br>2. Calculate maturation score<br>3. Setback association to disease events | Figure 3c<br>Figure 3d<br>Figure 4c<br>Supplementary Figure S11a, b<br>Supplementary Table S8 |

---

## **4_Segatella** — Segatella Analysis in the BBGUT Cohort

| Script | Function | DATA generated |
|--------|----------|----------------|
| *1_BBGUT_Segatella_final.R* | Exploring Segatella in the BBGUT cohort:<br>1. Segatella prevalence across GMM stages and health groups<br>2. HSA threshold calculation<br>3. Monthly proportion of HSA in adhoc disease and HTP samples<br>4. HSA metadata<br>5. Heatmap for samples with HSA spikes per child over time | Figure 5a, b, c<br>Supplementary Figure S5a, b<br>Supplementary Figure S15a, b<br>Supplementary Table S11 |
| *2_BBGUT_Segatella_SupplementaryfigureS16_final.R* | Segatella ASV abundance across health groups | Supplementary Figure S16 |
| *3_daa_HSAvsHTP_final.R* | Bacterial genera that differ significantly between the HSA and HTP control groups | Supplementary Figure S17 |

---

## **Data**

No scripts — data preprocessing on processed reads/OTU table to get the final data.

### **table_processing**

| Script | Function |
|--------|----------|
| *1_BBGUT_rdp_annotation_final.R* | OTU Annotation with RDP taxonomy (RDP set 19) |
| *2_BBGUT_decontamination_final.R* | Decontamination using decontam R package |
| *3_BBGUT_final_phyloseq_final.R* | Create final phyloseq object |

---

## **QC**

Scripts and data to replicate the extended quality control analysis across sequencing runs.

| Script | DATA generated |
|--------|----------------|
| *01_QC_and_FigS3_final2.R*<br>*02_reproduce_Response1_final2.R* | Supplementary Figure S3 |
