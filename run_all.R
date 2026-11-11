# =============================================================================
# run_all.R  - runs the R part of the pipeline in order.
# Step 01 (shell) must be run first; FAPROTAX and PICRUSt2 (commands printed
# by step 07) must be run before steps 08-10.
#
#   Rscript run_all.R            # steps 02-07
#   Rscript run_all.R post       # steps 08-12, after FAPROTAX + PICRUSt2
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
pre  <- c("02_dada2.R", "03_phyloseq_filter.R", "04_alpha_beta_diversity.R",
          "05_composition_heatmap.R", "06_core_microbiome.R", "07_export_for_faprotax_picrust2.R")
post <- c("08_faprotax_post.R", "09_picrust2_post.R", "10_chord_sankey.R",
          "11_cooccurrence_network.R", "12_export_supplementary.R")
steps <- if (length(args) && args[1] == "post") post else pre
for (s in steps) {
  message("\n########## ", s, " ##########")
  source(s, echo = FALSE, local = new.env())
}
