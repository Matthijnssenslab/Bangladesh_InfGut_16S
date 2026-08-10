################################################################################
# Reproduce Response 1 — Reviewer #2 revision analyses
##
## Inputs
##   data/bbgut_phyloseq_final.RData            (object: bbgut_ps)
##   outputs/TableS_D_control_library_QC.csv    (per-run control summary, from QC script)
##   outputs/TableS_I_positive_control_libraries.csv  (per-library positive controls)
##
## Output: outputs/reproduce_Response1_results.txt
################################################################################

suppressMessages({
  library(phyloseq)
  library(dplyr)
  library(vegan)
  library(splines)
})

## ---- configuration -----------------------------------------------------------
## Primary path (original location); falls back to project root if unavailable.
a <- commandArgs(trailingOnly = FALSE)
f <- grep("^--file=", a, value = TRUE)
if (length(f)) {
  script_dir <- dirname(sub("^--file=", "", f[1]))
} else {
  script_dir <- getwd()
}

hardcoded_path <- ""
if (dir.exists(hardcoded_path)) {
  workdir <- hardcoded_path
} else {
  workdir <- dirname(script_dir)   # project root (parent of scripts/)
}

outdir <- file.path(workdir, "outputs")
setwd(workdir)

NPERM <- 999         # permutations for PERMANOVA
DET   <- 0.005       # 0.5% detection threshold used in the manuscript
LOCK  <- as.Date("2020-03-26")   # Bangladesh general holiday / shutdown onset

sink(file.path(outdir, "reproduce_Response1_results.txt"), split = TRUE)

################################################################################
## 0. DATA LOADING 
################################################################################
e <- new.env()
load("data/bbgut_phyloseq_final.RData", envir = e)
ps <- e$bbgut_ps

## Remove plant and host contaminant taxa
PAT <- "chloroplast|mitochondri|streptophyt|cyanobacteri"
tt  <- as.data.frame(as(tax_table(ps), "matrix"))
isp <- apply(tt, 1, function(x) any(grepl(PAT, x, ignore.case = TRUE)))
ps  <- prune_taxa(!isp, ps)
ps  <- prune_taxa(taxa_sums(ps) > 0, ps)

## The phyloseq object is already rarefied — use it directly.
## Compute relative abundance matrices at genus and family rank.
relmat <- function(p, rank) {
  g <- tax_glom(p, taxrank = rank)
  g <- transform_sample_counts(g, function(x) x / sum(x))
  M <- as(otu_table(g), "matrix")
  if (!taxa_are_rows(g)) M <- t(M)
  rownames(M) <- as.data.frame(as(tax_table(g), "matrix"))[[rank]]
  M
}

G  <- relmat(ps, "Genus")
F_ <- relmat(ps, "Family")

## Extract and clean sample metadata
md <- data.frame(sample_data(ps)@.Data, stringsAsFactors = FALSE)
names(md) <- sample_variables(ps)
rownames(md) <- sample_names(ps)
md$id    <- rownames(md)
md$run   <- sub("p[0-9]+[.].*$", "", md$ID)
md$plate <- regmatches(md$ID, regexpr("p[0-9]+", md$ID))
md$Day   <- as.numeric(md$Day)
md$Month <- as.numeric(as.character(md$Month))

## Define sample subsets
mat_ids <- md$id[md$fecalsample_type == "maternal"]
inf_ids <- md$id[md$fecalsample_type == "infant" & !is.na(md$Day)]
htp_ids <- md$id[md$fecalsample_type == "infant" & md$sample_category == "HTP" &
                  !is.na(md$Day) & !is.na(md$Month)]

rownames(md) <- md$id

cat("################################################################\n")
cat("REPRODUCTION OF RESPONSE 1\n")
cat("  samples           :", length(inf_ids), "infant (", length(htp_ids),
    "healthy-timepoint ),", length(mat_ids), "maternal\n")
cat("################################################################\n")


################################################################################
## [ii] MATERNAL BASELINE
################################################################################
cat("\n==============================================================\n")
cat("[ii] THE MATERNAL BASELINE\n")
cat("==============================================================\n")

bm <- 100 * G["Bacteroides", mat_ids]
bact_mean <- mean(bm)
bact_med  <- median(bm)
bact_det  <- sum(bm >= 100 * DET)
bact_n    <- length(bm)

cat(sprintf("  Bacteroides in mothers: mean %.2f%%, median %.2f%%, detected (>=%.1f%%) in %d of %d\n",
            bact_mean, bact_med, 100 * DET, bact_det, bact_n))

pv <- 100 * F_["Prevotellaceae", mat_ids]
cat(sprintf("  Prevotellaceae in mothers: mean %.1f%%, median %.1f%%\n", mean(pv), median(pv)))

sg <- 100 * G["Segatella", mat_ids]
sg_det <- sum(sg >= 100 * DET)
sg_n   <- length(sg)
cat(sprintf("  Segatella in mothers: detected (>=%.1f%%) in %d of %d, mean %.1f%%\n",
            100 * DET, sg_det, sg_n, mean(sg)))

## Compute Bray-Curtis distance matrix (needed for section [v])
D <- as.matrix(vegdist(t(G), "bray"))

################################################################################
## [v] SEQUENCING DESIGN AND THE RUN x AGE PARTITION
################################################################################
cat("\n==============================================================\n")
cat("[v] SEQUencing DESIGN, RANDOMISATION AND THE RUN x AGE PARTITION\n")
cat("==============================================================\n")

md$study_year <- ifelse(is.na(md$Day), "maternal", ifelse(md$Day <= 365, "year1", "year2"))

cat("\nRun x sample type:\n")
print(table(md$run, md$fecalsample_type))

cat("\nRun x study year:\n")
print(table(md$run, md$study_year))

infants_per_run <- tapply(
  md$Child_ID[md$fecalsample_type == "infant"],
  md$run[md$fecalsample_type == "infant"],
  n_distinct
)

cat(sprintf("  %d runs, %d plates; infants per run: %s\n",
            n_distinct(md$run), n_distinct(md$plate),
            paste(range(infants_per_run), collapse = "-")))

mat_on_r28 <- sum(md$run == "DI22R28" & md$fecalsample_type == "maternal")
inf_on_r28 <- sum(md$run == "DI22R28" & md$fecalsample_type == "infant")

cat(sprintf("  Maternal on DI22R28: %d (of %d total), infant on DI22R28: %d\n",
            mat_on_r28, length(mat_ids), inf_on_r28))

n06 <- sum(md$run == "DI22R06" & md$study_year == "year2")
ny2 <- sum(md$study_year == "year2")

cat(sprintf("  DI22R06 carries %d of %d year-2 samples (%.0f%%)\n", n06, ny2, 100 * n06 / ny2))

## ---- infant Bray-Curtis: run vs age ------------------------------------------
mdi <- md[inf_ids, ]
mdi$run <- factor(mdi$run)
Di  <- as.dist(D[inf_ids, inf_ids])

aov_age <- aov(Day ~ run, data = mdi)
ss <- summary(aov_age)[[1]]
eta2 <- ss[["Sum Sq"]][1] / sum(ss[["Sum Sq"]])

cat(sprintf("  Run accounts for %.0f%% of variance in infant age (eta2 = %.3f, F = %.1f, p = %.3g)\n",
            100 * eta2, eta2, ss[["F value"]][1], ss[["Pr(>F)"]][1]))

set.seed(42)
a_run <- adonis2(Di ~ run,       data = mdi, permutations = NPERM, by = "terms")

set.seed(42)
a_seq <- adonis2(Di ~ Day + run, data = mdi, permutations = NPERM, by = "terms")

set.seed(42)
a_sp  <- adonis2(Di ~ ns(Day, df = 6) + run, data = mdi,
                 permutations = NPERM, by = "terms")

cat(sprintf("  PERMANOVA — run alone: R2 = %.4f (p = %.3f)\n", a_run$R2[1], a_run$`Pr(>F)`[1]))
cat(sprintf("    age term:            R2 = %.4f (p = %.3f)\n", a_seq$R2[1], a_seq$`Pr(>F)`[1]))
cat(sprintf("    run | age:           R2 = %.4f (p = %.3f)\n", a_seq$R2[2], a_seq$`Pr(>F)`[2]))
cat(sprintf("  Equal df — age (6-df spline): R2 = %.4f; run | age: R2 = %.4f (p = %.3f)\n",
            a_sp$R2[1], a_sp$R2[2], a_sp$`Pr(>F)`[2]))

## ---- within-run age significance ---------------------------------------------
cat("\n  Within-run age effect:\n")
wr <- lapply(levels(mdi$run), function(r) {
  ii <- mdi$id[mdi$run == r]
  set.seed(42)
  aa <- adonis2(as.dist(D[ii, ii]) ~ Day, data = mdi[ii, ], permutations = NPERM, by = "terms")
  data.frame(run = r, n = length(ii), R2_age = round(aa$R2[1], 4), p = aa$`Pr(>F)`[1])
}) %>% bind_rows() %>% as.data.frame()

print(wr)

n_sig <- sum(wr$p < 0.05)
cat(sprintf("  Significant in %d of %d runs; R2 range %.3f-%.3f\n",
            n_sig, nrow(wr), min(wr$R2_age), max(wr$R2_age)))

## ---- residual run effect among age-mixed runs --------------------------------
i1 <- mdi$id[mdi$Day <= 365]
sub <- mdi[i1, ]
sub$run <- droplevels(sub$run)

set.seed(42)
ay1 <- adonis2(as.dist(D[i1, i1]) ~ Day + run, data = sub, permutations = NPERM, by = "terms")

cat(sprintf("  Six age-mixed runs (DI22R06 excluded): n = %d, runs = %d\n",
            nrow(sub), nlevels(sub$run)))
cat(sprintf("    run | age R2 = %.4f (p = %.3f)\n", ay1$R2[2], ay1$`Pr(>F)`[2]))

## ---- second-year contrast internal to DI22R06 --------------------------------
y2 <- md %>% filter(id %in% htp_ids, Month >= 13) %>%
  mutate(win = ifelse(Month <= 16, "m13_16", "m17_24"))

cat("\n  Months 13-16 vs 17-24 by run:\n")
print(table(y2$run, y2$win))

on_r06 <- sum(y2$run == "DI22R06")
cat(sprintf("  %d of %d (%.0f%%) contrast samples are on DI22R06\n",
            on_r06, nrow(y2), 100 * mean(y2$run == "DI22R06")))

################################################################################
## [vi] CONTROLS AND THE AGE-ADJUSTED WITHIN-RUN TEST
################################################################################
cat("\n==============================================================\n")
cat("[vi] CONTROL LIBRARIES AND THE WITHIN-RUN Bacteroides TEST\n")
cat("==============================================================\n")

qc <- read.csv(file.path(outdir, "TableS_D_control_library_QC.csv"), stringsAsFactors = FALSE)
pc <- read.csv(file.path(outdir, "TableS_I_positive_control_libraries.csv"),
               stringsAsFactors = FALSE)

cat(sprintf("  Positive-control Bacteroides: mean %.1f%%, range %.1f-%.1f%%\n",
            mean(qc$PC_Bacteroides_pct), min(qc$PC_Bacteroides_pct), max(qc$PC_Bacteroides_pct)))
cat(sprintf("    DI22R28 : %.1f%%\n", qc$PC_Bacteroides_pct[qc$MiSeq_run == "DI22R28"]))
cat(sprintf("    Plate p18 library: %.1f%%\n", pc$Bacteroides_pct[pc$plate == "p18"]))

cat(sprintf("  Runella off-target (per-run): %.3f-%.3f%%\n",
            min(qc$Runella_offtarget_pct), max(qc$Runella_offtarget_pct)))
cat(sprintf("    implied on-target: %.2f-%.2f%%\n",
            100 - max(qc$Runella_offtarget_pct), 100 - min(qc$Runella_offtarget_pct)))

if ("Mock_offtarget_pct" %in% names(qc)) {
  cat(sprintf("  Mock off-target (per-run): %.3f-%.3f%%\n",
              min(qc$Mock_offtarget_pct, na.rm = TRUE), max(qc$Mock_offtarget_pct, na.rm = TRUE)))
  cat(sprintf("    implied on-target: %.2f-%.2f%%\n",
              100 - max(qc$Mock_offtarget_pct, na.rm = TRUE), 100 - min(qc$Mock_offtarget_pct, na.rm = TRUE)))
  cat(sprintf("    Mock libraries per run: %d\n", median(qc$Mock_n, na.rm = TRUE)))
}

qc$nc_ratio <- 100 * qc$NC_median_reads / qc$study_median_reads
cat(sprintf("  Negative control depth: %.2f%% (range %.2f-%.2f%%) of study-sample median\n",
            mean(qc$nc_ratio, na.rm = TRUE), min(qc$nc_ratio, na.rm = TRUE), max(qc$nc_ratio, na.rm = TRUE)))
cat(sprintf("  Runs with no negative control: %s\n",
            paste(qc$MiSeq_run[is.na(qc$NC_n)], collapse = ", ")))

################################################################################
## session info
################################################################################
cat("\n==============================================================\n")
cat("Session information\n")
cat("==============================================================\n")
print(sessionInfo())

sink()
cat("\nWrote reproduce_Response1_results.txt to", outdir, "\n")
