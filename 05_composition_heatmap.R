# =============================================================================
# 05_composition_heatmap.R
# Phylum and genus relative abundance (Tables S6, S7) and the clustered heatmap
# of the top 30 genera (log10 relative abundance, %).
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(pheatmap))

ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
meta <- data.frame(sample_data(ps_rare), check.names = FALSE)

abund_table <- function(level) {
  m  <- glom_named(ps_rare, level)
  ra <- 100 * rel_abund(m)
  out <- data.frame(taxon = rownames(ra),
                    mean_rel_abund_pct = rowMeans(ra),
                    sd_rel_abund_pct   = apply(ra, 1, sd),
                    prevalence_pct     = 100 * rowMeans(m > 0),
                    round(ra, 4), check.names = FALSE, row.names = NULL)
  names(out)[1] <- level
  out[order(-out$mean_rel_abund_pct), ]
}

log_step("Phylum / genus relative abundance")
phylum <- abund_table("Phylum")
genus  <- abund_table("Genus")
print(head(phylum[, 1:3], 8))
print(head(genus[, 1:3], 15))
write_supp(phylum, "Table_S6_phylum_relative_abundance")
write_supp(genus,  "Table_S7_genus_relative_abundance")

# ---- Heatmap: top 30 genera ---------------------------------------------------
log_step("Heatmap: top ", par$heatmap_top_genera, " genera")
g <- genus[!is_unclassified(genus$Genus), ]
top <- head(g, par$heatmap_top_genera)
mat <- as.matrix(top[, sample_names(ps_rare)]); rownames(mat) <- top$Genus
rownames(mat)[rownames(mat) == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium group"
mat_log <- log10(mat + par$heatmap_pseudo_pct)

ann <- data.frame(Dataset = meta[colnames(mat), "dataset"], row.names = colnames(mat))
ann_col <- list(Dataset = setNames(c("#0072B2", "#D55E00"), datasets))

save_fig("Fig_heatmap_top30_genera", width = 9, height = 8, draw = function() {
  pheatmap(mat_log, annotation_col = ann, annotation_colors = ann_col,
           clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",
           clustering_method = "complete",
           color = colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100),
           fontsize_row = 8, fontsize_col = 8, border_color = NA,
           main = sprintf("log10 relative abundance (%%) + %g", par$heatmap_pseudo_pct),
           silent = FALSE)
})
message("\nDone: 05_composition_heatmap.R")
