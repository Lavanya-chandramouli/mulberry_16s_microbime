# =============================================================================
# 00_config.R
# Shared paths, parameters and helper functions for the Paper 1 pipeline.
# Mulberry-associated bacterial microbiome: 16S rRNA amplicon reanalysis of
# SRP411931 and PRJNA822613.
#
# Every other script starts with source("00_config.R").
# Edit the PATHS block below before the first run; nothing else should need
# changing unless you deliberately change an analysis parameter.
# =============================================================================

# ---- PATHS ------------------------------------------------------------------
dir_project  <- normalizePath(".", mustWork = TRUE)
dir_trimmed  <- file.path(dir_project, "data", "trimmed")   # Cutadapt output, one sub-folder per dataset
dir_results  <- file.path(dir_project, "results")
dir_rds      <- file.path(dir_results, "rds")
dir_tables   <- file.path(dir_results, "tables")
dir_figs     <- file.path(dir_results, "figures")
dir_faprotax <- file.path(dir_results, "faprotax")
dir_picrust  <- file.path(dir_results, "picrust2")

# Sample sheet: one row per SRA run. Required columns:
#   sample_id  (must match the FASTQ file prefix, e.g. SRR22345678)
#   dataset    (SRP411931 or PRJNA822613)
# Optional: compartment (rhizosphere / endosphere / unknown), any other metadata
file_metadata <- file.path(dir_project, "data", "metadata.csv")

# SILVA v138.1 DADA2-formatted training set (genus level)
# https://zenodo.org/records/4587955  ->  silva_nr99_v138.1_train_set.fa.gz
file_silva <- file.path(dir_project, "ref", "silva_nr99_v138.1_train_set.fa.gz")

# ---- DATASETS ---------------------------------------------------------------
datasets <- c("SRP411931", "PRJNA822613")
# Trimmed read file naming (written by 01_cutadapt.sh)
fwd_suffix <- "_R1.trim.fastq.gz"
rev_suffix <- "_R2.trim.fastq.gz"

# ---- ANALYSIS PARAMETERS (as reported in the Methods) -------------------------
par <- list(
  # DADA2 (Section 2.3)
  truncLen     = c(270, 250),
  maxEE        = c(2, 5),
  truncQ       = 2,
  maxN         = 0,
  rm.phix      = TRUE,
  minOverlap   = 10,
  maxMismatch  = 1,
  chimera_method = "consensus",
  threads      = TRUE,          # TRUE = all cores; or an integer, e.g. 8
  # Taxonomy (Section 2.4)
  minBoot      = 80,
  # Read depth / rarefaction (Section 2.5)
  min_reads    = 10000,
  rare_seed    = 123,
  rare_replace = TRUE,          # phyloseq default; state it in the Methods
  # Beta diversity
  n_perm       = 999,
  # Core microbiome (Section 2.5)
  core_detection   = 0.0001,    # 0.01 % relative abundance
  core_prevalence  = 0.80,
  core_sensitivity = c(0.75, 0.85),
  # Figures
  heatmap_top_genera = 30,
  bubble_top_core    = 25,
  faprotax_top_fun   = 20,
  heatmap_pseudo_pct = 0.001,   # pseudocount (%) added before log10 in genus heatmap
  # Network (Section 2.8)
  net_top_genera = 30,
  net_r_cutoff   = 0.6,
  net_exclude_unclassified = TRUE,
  # Figure export
  dpi = 300
)

# Candidate genera carried into the integrated (chord / Sankey) figures.
# Names are SILVA v138.1 genus labels; display names on the right.
candidate_genera <- c(
  "Pseudomonas"         = "Pseudomonas",
  "Bacillus"            = "Bacillus",
  "Bdellovibrio"        = "Bdellovibrio",
  "Noviherbaspirillum"  = "Noviherbaspirillum",
  "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" = "Rhizobium group"
)

# FAPROTAX functions used in Table 3 / chord / Sankey, with display category
faprotax_focus <- c(
  nitrogen_fixation          = "Nitrogen cycling",
  nitrification              = "Nitrogen cycling",
  nitrate_reduction          = "Nitrogen cycling",
  denitrification            = "Nitrogen cycling",
  ureolysis                  = "Nitrogen cycling",
  chitinolysis               = "Antagonism-associated",
  predatory_or_exoparasitic  = "Antagonism-associated"
)

# PICRUSt2 KOs used in Table 4.  The "label" column is what the manuscript
# currently calls each KO.  09_picrust2_post.R adds the official KEGG
# description next to it so mislabels are visible (see README: K02014, K01580).
picrust_focus <- data.frame(
  KO       = c("K02014", "K01915", "K01077", "K00466", "K01580"),
  label    = c("Siderophore biosynthesis (manuscript label)",
               "Glutamine synthetase (N assimilation)",
               "Alkaline phosphatase (P mobilisation)",
               "Tryptophan 2-monooxygenase (IAA, IAM pathway)",
               "IAA biosynthesis pathway 2 (manuscript label)"),
  category = c("Iron acquisition", "Nitrogen cycling", "Phosphorus mobilisation",
               "Phytohormone", "Phytohormone"),
  stringsAsFactors = FALSE
)

# ---- SETUP ------------------------------------------------------------------
for (d in c(dir_results, dir_rds, dir_tables, dir_figs, dir_faprotax, dir_picrust))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(phyloseq)
  library(ggplot2)
})

# ---- HELPERS ----------------------------------------------------------------

# Save a ggplot or base-graphics figure as PDF + 300-dpi PNG
save_fig <- function(name, plot = NULL, width = 7, height = 5, draw = NULL) {
  pdf_f <- file.path(dir_figs, paste0(name, ".pdf"))
  png_f <- file.path(dir_figs, paste0(name, ".png"))
  if (!is.null(plot)) {
    ggsave(pdf_f, plot, width = width, height = height, bg = "white")
    ggsave(png_f, plot, width = width, height = height, dpi = par$dpi, bg = "white")
  } else {
    pdf(pdf_f, width = width, height = height); draw(); dev.off()
    png(png_f, width = width, height = height, units = "in", res = par$dpi); draw(); dev.off()
  }
  message("  figure: ", name, " (.pdf, .png)")
}

# Write one supplementary table as CSV and keep a copy for the combined workbook
write_supp <- function(df, name) {
  write.csv(df, file.path(dir_tables, paste0(name, ".csv")), row.names = FALSE)
  saveRDS(df, file.path(dir_rds, paste0("supp_", name, ".rds")))
  message("  table:  ", name, " (", nrow(df), " rows)")
}

# Label each taxon by its most precise assigned rank, e.g.
# "Pseudomonas" or "Unclassified Comamonadaceae"
lowest_label <- function(tax, level = "Genus") {
  tax <- as.data.frame(tax, stringsAsFactors = FALSE)
  ranks <- colnames(tax)[seq_len(match(level, colnames(tax)))]
  apply(tax[, ranks, drop = FALSE], 1, function(r) {
    r <- r[!is.na(r) & r != ""]
    if (length(r) == 0) return("Unclassified")
    if (names(r)[length(r)] == level) r[length(r)] else paste("Unclassified", r[length(r)])
  })
}

# Aggregate a phyloseq object to a rank, using lowest_label() names
glom_named <- function(ps, level = "Genus") {
  tax <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
  tax$label <- lowest_label(tax, level)
  otu <- as(otu_table(ps), "matrix")
  if (!taxa_are_rows(ps)) otu <- t(otu)
  agg <- rowsum(otu, group = tax$label[match(rownames(otu), rownames(tax))])
  agg
}

rel_abund <- function(m) sweep(m, 2, colSums(m), "/")   # columns = samples

is_unclassified <- function(x) grepl("^Unclassified", x)

log_step <- function(...) message("\n== ", ..., " ==")
