# =============================================================================
# 06_core_microbiome.R
# Candidate core microbiome (Methods 2.5): microbiome::core_members on
# relative abundances, detection 0.0001 (0.01 %), prevalence 0.80, with
# sensitivity at 0.75 and 0.85 and per-dataset genus-level prevalence.
#
# Output: Table S3 (ASV-level core, genus-level per-dataset check, sensitivity),
#         core bubble plot (top 25), Supplementary_File_S1_core_ASVs.fasta
#
# NOTE: microbiome::core_members keeps taxa with prevalence STRICTLY GREATER
# than the threshold and abundance STRICTLY GREATER than detection.
# With n = 21, "> 0.80" means present in >= 17 samples.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(Biostrings))

ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
meta    <- data.frame(sample_data(ps_rare), check.names = FALSE)
ps_rel  <- transform_sample_counts(ps_rare, function(x) x / sum(x))

otu <- as(otu_table(ps_rel), "matrix"); if (!taxa_are_rows(ps_rel)) otu <- t(otu)
prev_base <- function(m, det = par$core_detection) rowMeans(m > det)

core_ids <- function(prev_thr) {
  if (requireNamespace("microbiome", quietly = TRUE)) {
    microbiome::core_members(ps_rel, detection = par$core_detection, prevalence = prev_thr)
  } else {
    message("  microbiome package not installed - using the identical base-R rule")
    names(which(prev_base(otu) > prev_thr))
  }
}

log_step("Core microbiome: detection ", par$core_detection, ", prevalence ", par$core_prevalence)
core <- core_ids(par$core_prevalence)
check <- names(which(prev_base(otu) > par$core_prevalence))
stopifnot(setequal(core, check))   # independent cross-check of the package result
message("  candidate core ASVs: ", length(core), "  (manuscript: 80)")

tax <- as.data.frame(tax_table(ps_rare), stringsAsFactors = FALSE)
tax$Genus_label <- lowest_label(tax, "Genus")
per_ds_prev <- sapply(datasets, function(d) prev_base(otu[, meta$dataset == d, drop = FALSE]))

core_tab <- data.frame(
  ASV_ID = core,
  tax[core, c("Phylum", "Class", "Order", "Family", "Genus")],
  Genus_label = tax[core, "Genus_label"],
  prevalence_pooled_pct = round(100 * prev_base(otu)[core], 1),
  round(100 * per_ds_prev[core, , drop = FALSE], 1),
  mean_rel_abund_pct = round(100 * rowMeans(otu[core, , drop = FALSE]), 4),
  sequence = as.character(refseq(ps_rare)[core]),
  check.names = FALSE, row.names = NULL)
names(core_tab)[names(core_tab) %in% datasets] <- paste0("prevalence_", datasets, "_pct")
core_tab <- core_tab[order(-core_tab$mean_rel_abund_pct), ]

# Genus-level prevalence, pooled and per dataset (the "core in both datasets" check)
gm  <- rel_abund(glom_named(ps_rare, "Genus"))
gen_tab <- data.frame(
  Genus = rownames(gm),
  prevalence_pooled_pct = round(100 * prev_base(gm), 1),
  sapply(datasets, function(d) round(100 * prev_base(gm[, meta$dataset == d, drop = FALSE]), 1)),
  mean_rel_abund_pct = round(100 * rowMeans(gm), 4),
  n_core_ASVs = as.vector(table(factor(core_tab$Genus_label, levels = rownames(gm)))),
  check.names = FALSE, row.names = NULL)
names(gen_tab)[names(gen_tab) %in% datasets] <- paste0("prevalence_", datasets, "_pct")
gen_tab <- gen_tab[gen_tab$n_core_ASVs > 0, ]
gen_tab <- gen_tab[order(-gen_tab$mean_rel_abund_pct), ]
print(head(gen_tab, 15))

# Sensitivity analysis (Supplementary Note 2)
sens <- data.frame(prevalence_threshold = c(par$core_sensitivity[1], par$core_prevalence,
                                            par$core_sensitivity[2]))
sens$n_core_ASVs <- sapply(sens$prevalence_threshold, function(p) length(core_ids(p)))
sens$candidate_genera_retained <- sapply(sens$prevalence_threshold, function(p) {
  g <- unique(tax[core_ids(p), "Genus"])
  paste(intersect(names(candidate_genera), g), collapse = "; ")
})
print(sens)

write_supp(core_tab, "Table_S3a_core_ASVs")
write_supp(gen_tab,  "Table_S3b_core_genera_per_dataset")
write_supp(sens,     "Table_S3c_core_sensitivity")

# Supplementary File: representative sequences of the core ASVs (needed by Paper 2)
fa <- refseq(ps_rare)[core_tab$ASV_ID]
names(fa) <- paste(core_tab$ASV_ID, core_tab$Genus_label)
writeXStringSet(fa, file.path(dir_tables, "Supplementary_File_S1_core_ASVs.fasta"))
message("  FASTA: Supplementary_File_S1_core_ASVs.fasta (", length(fa), " sequences)")

# ---- Bubble plot: top 25 core ASVs --------------------------------------------
top <- head(core_tab, par$bubble_top_core)
top$label <- paste0(sub("Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium", "Rhizobium group",
                        top$Genus_label), " (", top$ASV_ID, ")")
top$label <- factor(top$label, levels = rev(top$label))
p_bub <- ggplot(top, aes(x = prevalence_pooled_pct, y = label)) +
  geom_point(aes(size = mean_rel_abund_pct, fill = Phylum), shape = 21, colour = "grey20", alpha = 0.85) +
  scale_size_area(max_size = 12, name = "Mean relative\nabundance (%)") +
  scale_x_continuous(limits = c(80, 100.5), breaks = seq(80, 100, 5)) +
  labs(x = "Prevalence across samples (%)", y = NULL) +
  theme_bw(base_size = 10) +
  guides(fill = guide_legend(override.aes = list(size = 4)))
save_fig("Fig_core_bubble_top25", p_bub, width = 7.5, height = 7)

saveRDS(core, file.path(dir_rds, "core_asvs.rds"))
message("\nDone: 06_core_microbiome.R")
