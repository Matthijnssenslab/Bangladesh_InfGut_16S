################################################################################
## Figure S3 — Expanded quality control across the seven MiSeq runs
##
## ENTRY POINT: data/FigS_3_data.RData  (the only input)
##
## This script starts from a pre-built phyloseq object and reproduces the QC
## figure, tables, and every number quoted in the figure legend. It needs
## nothing else: no raw run directories, no taxonomy database, no network.
##
## Provenance of the object — raw ASV tables -> DADA2 + RDP 19 -> plant/host/
##  non-prokaryote removal -> genus agglomeration
##
## INPUT   data/FigS_3_data.RData       phyloseq object `qc_ps`
##                                  417 genera x 1148 libraries, integer counts,
##                                  NOT rarefied (rarefaction happens below)
##
## OUTPUT  outputs/FigS_3_expanded_QC.pdf  
##         outputs/FigS_3_expanded_QC.tiff 
##         outputs/TableS_D_control_library_QC.csv        per run
##         outputs/TableS_I_positive_control_libraries.csv per positive-control library
##         (stdout)                 every value quoted in the figure legend
##
##
## USAGE   Rscript scripts/01_QC_and_FigS3.R [/path/to/FigS_3_data.RData]
##         (with no argument the object is looked for in data/)
##
## R >= 4.1. Packages: phyloseq, vegan, ggplot2, patchwork, dplyr.
################################################################################

## ---- 0. configuration --------------------------------------------------------
SEED         <- 20260716   # rarefaction / sampling seed; fixes reproducibility
PC_DEPTH     <- NULL       # NULL = rarefy positive controls to their minimum
STUDY_DEPTH  <- 11070      # depth for the biological-variation comparison
N_STUDY      <- 200        # study samples sampled for that comparison
FOCAL_GENUS  <- "Bacteroides"
MATERNAL_RUN <- "DI22R28"  # the run carrying all 18 maternal samples
MATERNAL_PLATE <- "p18"    # the plate within it carrying the maternal samples
RESTORE_MOCK <- TRUE      # TRUE re-admits the mock to the read-depth panel

## ---- 1. dependencies ---------------------------------------------------------
need <- c("phyloseq", "vegan", "ggplot2", "patchwork", "dplyr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) {
  stop("missing package(s): ", paste(miss, collapse = ", "),
       "\n  install.packages(c(", paste0('"', setdiff(miss, "phyloseq"), '"', collapse = ", "), "))",
       if ("phyloseq" %in% miss) "\n  BiocManager::install('phyloseq')" else "")
}

suppressMessages({
  library(phyloseq)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

## ---- 2. locate and load the object -------------------------------------------
script_dir <- function() {
  a <- commandArgs(FALSE)
  f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(normalizePath(dirname(sub("^--file=", "", f[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- rstudioapi::getActiveDocumentContext()$path
    if (nzchar(p)) return(dirname(p))
  }
  getwd()
}

args   <- commandArgs(trailingOnly = TRUE)
here   <- script_dir()
DATA   <- if (length(args) >= 1) args[1] else file.path(here, "..", "data", "FigS_3_data.RData")
OUTDIR <- file.path(here, "..", "outputs")

if (!file.exists(DATA)) {
  stop("cannot find the phyloseq object: ", DATA,
       "\n  pass its path as the first argument, e.g.",
       "\n  Rscript scripts/01_QC_and_FigS3.R /path/to/FigS_3_data.RData")
}

e <- new.env()
loaded <- load(DATA, envir = e)

if (!"qc_ps" %in% loaded) {
  stop(DATA, " does not contain an object named 'qc_ps' (found: ",
       paste(loaded, collapse = ", "), ")")
}

qc_ps <- e$qc_ps

stopifnot(inherits(qc_ps, "phyloseq"))
cat("loaded", basename(DATA), ":", ntaxa(qc_ps), "genera x",
    nsamples(qc_ps), "libraries\n")

## ---- 3. validate ------------------------------------------------------------
M   <- as(otu_table(qc_ps), "matrix")
if (!taxa_are_rows(qc_ps)) M <- t(M)

smd <- data.frame(sample_data(qc_ps)@.Data, stringsAsFactors = FALSE)
names(smd) <- sample_variables(qc_ps)
rownames(smd) <- sample_names(qc_ps)

req <- c("run", "plate", "library_class", "depth", "Source")
if (!all(req %in% names(smd))) {
  stop("sample_data is missing: ", paste(setdiff(req, names(smd)), collapse = ", "))
}

if (!FOCAL_GENUS %in% rownames(M)) {
  stop("genus '", FOCAL_GENUS, "' is not present in the object")
}

if (any(M %% 1 != 0)) {
  warning("counts are not integers — was this object already normalised? ",
          "rarefaction assumes raw counts.")
}

bad <- grep("chloroplast|mitochondri|streptophyt|cyanobacteri",
            rownames(M), ignore.case = TRUE, value = TRUE)
if (length(bad)) {
  stop("plant/host taxa are still present: ", paste(bad, collapse = ", "),
       "\n  these should have been removed before normalisation")
}

cat("validated: no plant/host taxa; '", FOCAL_GENUS, "' present; counts are integers\n", sep = "")

RUNS <- sort(unique(smd$run))
cat("runs:", paste(RUNS, collapse = ", "), "\n")
cat("library classes:\n")
print(table(smd$library_class))

## ---- 4. helpers --------------------------------------------------------------
rare_rel <- function(sids, depth = NULL, seed = SEED) {
  sub <- M[, sids, drop = FALSE]
  sub <- sub[rowSums(sub) > 0, , drop = FALSE]
  if (is.null(depth)) depth <- min(colSums(sub))
  if (any(colSums(sub) < depth)) {
    stop("depth ", depth, " exceeds the smallest library (", min(colSums(sub)), ")")
  }
  set.seed(seed)
  r <- t(vegan::rrarefy(t(sub), sample = depth))
  attr(r, "depth") <- depth
  sweep(r, 2, colSums(r), "/")
}

getrow <- function(mat, g) {
  if (g %in% rownames(mat)) mat[g, ] else setNames(rep(0, ncol(mat)), colnames(mat))
}

short <- function(x) sub("^DI22", "", x)

## ---- 5. palette and theme ----------------------------------------------------
PAL     <- c(blue = "#2a78d6", amber = "#eda100", pink = "#e87ba4")
INK     <- c(primary = "#1a1a19", secondary = "#4a4a48", muted = "#767472")
GRIDCOL <- "#e6e4e1"
SURFACE <- "#fcfcfb"

theme_pub <- function(base = 7) {
  theme_minimal(base_size = base, base_family = "Helvetica") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = GRIDCOL, linewidth = 0.25),
      panel.background = element_rect(fill = SURFACE, colour = NA),
      plot.background  = element_rect(fill = SURFACE, colour = NA),
      axis.line   = element_line(colour = INK["muted"], linewidth = 0.3),
      axis.ticks  = element_line(colour = INK["muted"], linewidth = 0.3),
      axis.text   = element_text(colour = INK["secondary"], size = base - 0.5),
      axis.title  = element_text(colour = INK["primary"], size = base),
      plot.title    = element_text(colour = INK["primary"], size = base + 1, face = "bold"),
      plot.subtitle = element_text(colour = INK["secondary"], size = base - 0.5),
      plot.tag      = element_text(colour = INK["primary"], size = base + 3, face = "bold"),
      legend.text  = element_text(colour = INK["secondary"], size = base - 0.5),
      legend.key.size = unit(3, "mm"),
      strip.background = element_blank()
    )
}

## ---- 6. panel a — Bacteroides in the identical positive-control DNA ----------
pc  <- smd %>% filter(library_class == "Positive control")
if (!nrow(pc)) stop("no libraries with library_class == 'Positive control'")

pcR <- rare_rel(rownames(pc), depth = PC_DEPTH)
PCD <- attr(pcR, "depth")

pcb <- data.frame(
  run = pc$run,
  plate = pc$plate,
  raw_reads = pc$depth,
  Bacteroides = 100 * getrow(pcR, FOCAL_GENUS),
  row.names = NULL
) %>% mutate(
  grp = ifelse(run == MATERNAL_RUN,
               paste0(MATERNAL_RUN, " (carries all mothers)"),
               "Other runs"),
  run_short = short(run)
)

gmean <- mean(pcb$Bacteroides)
ref   <- pcb %>% group_by(run) %>% summarise(m = mean(Bacteroides), .groups = "drop")

oth   <- ref$m[ref$run != MATERNAL_RUN]
z28   <- (ref$m[ref$run == MATERNAL_RUN] - mean(oth)) / sd(oth)

pA <- ggplot(pcb, aes(run_short, Bacteroides, colour = grp)) +
  geom_hline(yintercept = gmean, linetype = "dashed",
             colour = INK["muted"], linewidth = 0.25) +
  geom_point(size = 1.6, stroke = 0.4) +
  annotate("text", x = 0.55, y = gmean - 0.9, hjust = 0, size = 1.9,
           colour = INK["muted"], label = sprintf("cross-run mean %.1f%%", gmean)) +
  scale_colour_manual(
    values = setNames(c(PAL[["pink"]], PAL[["blue"]]),
                      c(paste0(MATERNAL_RUN, " (carries all mothers)"), "Other runs")),
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(floor(min(pcb$Bacteroides)) - 1.5, ceiling(max(pcb$Bacteroides)) + 1.5)
  ) +
  labs(
    x = "MiSeq run",
    y = bquote(italic(.(FOCAL_GENUS)) ~ "(% of reads)"),
    subtitle = "Identical faecal DNA, sequenced on all seven runs"
  ) +
  theme_pub() +
  theme(panel.grid.major.x = element_blank(),
        legend.position = "top", legend.margin = margin(b = -4))

## ---- 7. panel b — technical vs biological dissimilarity ----------------------
d_pc <- as.matrix(vegdist(t(pcR), "bray"))

bio  <- rownames(smd)[smd$library_class == "Study sample" &
                      smd$Source == "FEC" & smd$depth >= STUDY_DEPTH]
if (length(bio) < N_STUDY) {
  stop("only ", length(bio), " study samples at depth >= ", STUDY_DEPTH)
}

set.seed(1)
bio <- sample(bio, N_STUDY)
d_bio <- as.matrix(vegdist(t(rare_rel(bio, depth = STUDY_DEPTH)), "bray"))

bc <- rbind(
  data.frame(bc = d_pc[upper.tri(d_pc)],   src = "Technical\n(control vs control)"),
  data.frame(bc = d_bio[upper.tri(d_bio)], src = "Biological\n(sample vs sample)")
)
bc$src <- factor(bc$src, levels = c("Technical\n(control vs control)",
                                     "Biological\n(sample vs sample)"))

pB <- ggplot(bc, aes(src, bc, fill = src)) +
  geom_violin(colour = NA, width = 0.75, alpha = 0.85) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = SURFACE,
               colour = INK["primary"], linewidth = 0.25) +
  scale_fill_manual(values = c(PAL[["blue"]], PAL[["amber"]]), guide = "none") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = NULL, y = "Bray-Curtis dissimilarity",
       subtitle = "Technical vs biological variation") +
  theme_pub() +
  theme(panel.grid.major.x = element_blank())

## ---- 8. panel c — read depth by library class (mock excluded) ----------------
classes <- c("Study sample", "Positive control", "Runella slithyformis",
             "Negative control")
if (RESTORE_MOCK) classes <- append(classes, "Mock community", after = 2)

labs_c <- lapply(classes, function(x) {
  if (grepl("Runella", x)) bquote(italic(.(x))) else bquote(.(x))
})

dep <- smd %>% filter(library_class %in% classes) %>%
  mutate(library_class = factor(library_class, levels = rev(classes)))

pC <- ggplot(dep, aes(depth + 1, library_class)) +
  geom_jitter(
    position = position_jitter(width = 0, height = 0.18, seed = 123),
    size = 0.35,
    alpha = 0.4,
    colour = PAL[["blue"]])+
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
               width = 0.5, linewidth = 0.45, colour = INK["primary"],
               orientation = "y") +
  scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5),
                labels = c("100", "1,000", "10,000", "100,000")) +
  scale_y_discrete(labels = rev(labs_c)) +
  labs(x = "Reads per library (log scale)", y = NULL,
       subtitle = "Read depth by library class") +
  theme_pub() +
  theme(panel.grid.major.y = element_blank())

## ---- 9. assemble and export --------------------------------------------------
fig <- ((pA | pB) + plot_layout(widths = c(1.5, 1))) / pC +
  plot_layout(heights = c(1, 0.72)) +
  plot_annotation(
    tag_levels = "a",
    title = "Expanded quality control across the seven MiSeq runs",
    theme = theme(
      plot.title = element_text(size = 8, face = "bold", colour = INK["primary"]),
      plot.background = element_rect(fill = SURFACE, colour = NA)
    )
  )

ggsave(file.path(OUTDIR, "FigS_3_expanded_QC.pdf"), fig, width = 174, height = 118,
       units = "mm", device = "pdf", useDingbats = FALSE)
ggsave(file.path(OUTDIR, "FigS_3_expanded_QC.tiff"), fig, width = 174, height = 118,
       units = "mm", dpi = 300, compression = "lzw")
cat("\nwrote FigS_3_expanded_QC.pdf / .tiff\n")

## ---- 10. tables --------------------------------------------------------------
rs <- smd %>% filter(library_class == "Runella slithyformis")
rs_off <- NULL

if (nrow(rs)) {
  rsR  <- rare_rel(rownames(rs))
  RSD <- attr(rsR, "depth")
  rrow <- grep("Runella", rownames(rsR), value = TRUE, ignore.case = TRUE)
  rs_off <- data.frame(run = rs$run,
                        off = 100 * (1 - colSums(rsR[rrow, , drop = FALSE])))
}

## ---- mock community off-target ------------------------------------------------
MOCK_TARGETS <- c("Limosilactobacillus", "Salmonella", "Bacillus",
                  "Escherichia/Shigella", "Staphylococcus", "Pseudomonas",
                  "Listeria", "Enterococcus")

mock_ctrl <- smd %>% filter(library_class == "Mock community")
mock_off  <- NULL

if (nrow(mock_ctrl) >= 2) {
  mockR  <- rare_rel(rownames(mock_ctrl))
  MCD    <- attr(mockR, "depth")
  present <- MOCK_TARGETS[MOCK_TARGETS %in% rownames(mockR)]
  mock_off <- data.frame(run = mock_ctrl$run,
                          off = 100 * (1 - colSums(mockR[present, , drop = FALSE])))
}

neg <- smd %>% filter(library_class == "Negative control")

tabD <- data.frame(MiSeq_run = RUNS) %>%
  left_join(
    pcb %>% group_by(MiSeq_run = run) %>%
      summarise(PC_n = n(), PC_Bacteroides_pct = round(mean(Bacteroides), 2), .groups = "drop"),
    by = "MiSeq_run"
  ) %>%
  left_join(
    if (!is.null(rs_off)) {
      rs_off %>% group_by(MiSeq_run = run) %>%
        summarise(Runella_n = n(), Runella_offtarget_pct = round(mean(off), 3), .groups = "drop")
    } else {
      data.frame(MiSeq_run = character(), Runella_n = integer(),
                 Runella_offtarget_pct = numeric())
    },
    by = "MiSeq_run"
  ) %>%
  left_join(
    if (!is.null(mock_off)) {
      mock_off %>% group_by(MiSeq_run = run) %>%
        summarise(Mock_n = n(), Mock_offtarget_pct = round(mean(off), 3), .groups = "drop")
    } else {
      data.frame(MiSeq_run = character(), Mock_n = integer(),
                 Mock_offtarget_pct = numeric())
    },
    by = "MiSeq_run"
  ) %>%
  left_join(
    neg %>% group_by(MiSeq_run = run) %>%
      summarise(NC_n = n(), NC_median_reads = median(depth), .groups = "drop"),
    by = "MiSeq_run"
  ) %>%
  left_join(
    smd %>% filter(library_class == "Study sample") %>%
      group_by(MiSeq_run = run) %>%
      summarise(study_n = n(), study_median_reads = median(depth), .groups = "drop"),
    by = "MiSeq_run"
  )

write.csv(tabD, file.path(OUTDIR, "TableS_D_control_library_QC.csv"), row.names = FALSE)
cat("wrote TableS_D_control_library_QC.csv (", nrow(tabD), "rows)\n")

tabI <- pcb %>%
  transmute(MiSeq_run = run, plate, raw_reads,
            rarefied_to = PCD, Bacteroides_pct = round(Bacteroides, 3)) %>%
  arrange(MiSeq_run, plate)

write.csv(tabI, file.path(OUTDIR, "TableS_I_positive_control_libraries.csv"), row.names = FALSE)
cat("wrote TableS_I_positive_control_libraries.csv (", nrow(tabI), "rows)\n")

## ---- 11. values quoted in the legend -----------------------------------------
cat("\n================ values for the figure legend ================\n")

cat("panel a — positive control\n")
cat("  libraries        :", nrow(pcb), "( rarefied to", PCD, "reads )\n")
cat("  cross-run mean   :", sprintf("%.2f %%", gmean), "\n")
cat("  run means, range :", sprintf("%.2f-%.2f %%", min(ref$m), max(ref$m)), "\n")
cat("  library range    :", sprintf("%.1f-%.1f %%", min(pcb$Bacteroides), max(pcb$Bacteroides)), "\n")
cat("  ", MATERNAL_RUN, " mean       : ", sprintf("%.2f %%", ref$m[ref$run == MATERNAL_RUN]),
    "   z vs other runs = ", sprintf("%.2f", z28), "\n", sep = "")
cat("  plate", MATERNAL_PLATE, "library :",
    sprintf("%.2f %%", pcb$Bacteroides[pcb$plate == MATERNAL_PLATE]), "\n")

cat("\n  leave-one-out z for every run (each scored against the other six):\n")
zz <- vapply(ref$run, function(r) {
  o <- ref$m[ref$run != r]
  (ref$m[ref$run == r] - mean(o)) / sd(o)
}, numeric(1))
print(round(sort(zz), 2))

cat("\npanel b — dissimilarity\n")
cat("  technical  (control vs control), median BC:",
    sprintf("%.3f", median(d_pc[upper.tri(d_pc)])), "\n")
cat("  biological (sample vs sample),   median BC:",
    sprintf("%.3f", median(d_bio[upper.tri(d_bio)])), "\n")
cat("  ratio: technical is ~", sprintf("%.1f", median(d_bio[upper.tri(d_bio)]) /
                                      median(d_pc[upper.tri(d_pc)])), "x smaller\n", sep = "")

cat("\n  PERMANOVA — does run explain positive-control composition?\n")
set.seed(42)
print(adonis2(as.dist(d_pc) ~ run, data = pc, permutations = 999))

cat("\npanel c — read depth\n")
print(smd %>% filter(library_class %in% classes) %>%
        group_by(library_class) %>%
        summarise(n = n(), median_reads = median(depth), min = min(depth),
                  max = max(depth), .groups = "drop") %>% as.data.frame())

if (nrow(neg)) {
  cat("  negative controls are ",
      sprintf("%.2f", 100 * median(neg$depth) /
                    median(smd$depth[smd$library_class == "Study sample"])),
      "% of median study-sample depth\n", sep = "")
}

if (!is.null(rs_off)) {
  cat("\n  Runella off-target per run (%; rarefied to ", RSD, " reads):\n", sep = "")
  print(rs_off %>% group_by(run) %>%
          summarise(mean_offtarget_pct = round(mean(off), 3), .groups = "drop") %>% as.data.frame())
}

if (!is.null(mock_off)) {
  cat("\n  Mock off-target per run (%; rarefied to ", MCD, " reads):\n", sep = "")
  print(mock_off %>% group_by(run) %>%
          summarise(mean_offtarget_pct = round(mean(off), 3), .groups = "drop") %>% as.data.frame())
}

## ---- 12. mock community reproducibility ------------------------------------
cat("\n================ MOCK COMMUNITY REPRODUCIBILITY ================\n")
mock <- smd %>% filter(library_class == "Mock community")

if (nrow(mock) < 2) {
  cat("SKIPPED: fewer than 2 mock libraries found (", nrow(mock), ")\n\n")
} else {
  cat("  Mock libraries :", nrow(mock), "\n")
  cat("  runs           :", paste(sort(unique(mock$run)), collapse = ", "), "\n")

  mockR <- rare_rel(rownames(mock))
  MCD   <- attr(mockR, "depth")
  cat("  rarefied to    :", MCD, "reads\n\n")

  ## ---- pairwise Pearson correlation ----------------------------------------
  cat("  Pairwise Pearson correlation between mock profiles:\n")
  mock_cor <- cor(mockR, method = "pearson")
  diag(mock_cor) <- NA
  mock_cor_values <- mock_cor[upper.tri(mock_cor)]

  cat(sprintf("    median : %.4f\n", median(mock_cor_values, na.rm = TRUE)))
  cat(sprintf("    min    : %.4f\n", min(mock_cor_values,   na.rm = TRUE)))
  cat(sprintf("    max    : %.4f\n\n", max(mock_cor_values, na.rm = TRUE)))

  ## ---- correlation to consensus profile ------------------------------------
  cat("  Correlation of each mock library to the consensus profile:\n")
  mock_consensus <- rowMeans(mockR)
  mock_consistency <- data.frame(
    sample = colnames(mockR),
    correlation_to_consensus = apply(
      mockR, 2, function(x) cor(x, mock_consensus, method = "pearson")))

  print(mock_consistency, row.names = FALSE)
  cat(sprintf("    lowest : %.4f\n\n", min(mock_consistency$correlation_to_consensus)))

  ## ---- Bray-Curtis dissimilarity -------------------------------------------
  cat("  Bray-Curtis dissimilarity between mock libraries:\n")
  mock_bc <- vegdist(t(mockR), method = "bray")

  cat(sprintf("    median : %.4f\n", median(mock_bc)))
  cat(sprintf("    min    : %.4f\n", min(mock_bc)))
  cat(sprintf("    max    : %.4f\n\n", max(mock_bc)))

  ## ---- RDA ordination (visualisation only) ---------------------------------
  ord <- rda(t(mockR))
  plot(ord, display = "sites", type = "n")
  points(ord, display = "sites", pch = 19, cex = 1.5)
  text(ord,   display = "sites", cex = 0.8)
}

## ---- session info ----------------------------------------------------------
cat("\n================ session ================\n")
cat(R.version.string, "\n")
for (p in need) cat(sprintf("  %-10s %s\n", p, as.character(packageVersion(p))))
cat("seed:", SEED, "\n")
cat("\ndone.\n")

