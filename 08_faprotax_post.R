# =============================================================================
# 08_faprotax_post.R
# Reads FAPROTAX v1.2.12 output: Table S4 (all predicted functions), the
# assignment rate, per-genus contributions to each function (used for Table 3
# "Dominant contributing taxa" and for the chord / Sankey figures), and the
# heatmap of the top 20 functions.
# All values are taxonomically inferred predictions, not measured activities.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(pheatmap))

ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
meta <- data.frame(sample_data(ps_rare), check.names = FALSE)
otu <- as(otu_table(ps_rare), "matrix"); if (!taxa_are_rows(ps_rare)) otu <- t(otu)
tax <- as.data.frame(tax_table(ps_rare), stringsAsFactors = FALSE)

# ---- Function x sample table ---------------------------------------------------
f_out <- file.path(dir_faprotax, "faprotax_functions.tsv")
fun <- read.delim(f_out, check.names = FALSE, comment.char = "", stringsAsFactors = FALSE)
names(fun)[1] <- "function_group"; fun$function_group <- sub("^#", "", fun$function_group)
fm <- as.matrix(fun[, -1]); rownames(fm) <- fun$function_group
fm <- fm[rowSums(fm) > 0, , drop = FALSE]

# ---- Which ASVs were assigned to which function (sub-tables directory) -------
sub_dir <- file.path(dir_faprotax, "sub_tables")
sub_files <- list.files(sub_dir, pattern = "\\.(tsv|txt)$", full.names = TRUE)
read_ids <- function(f) {
  l <- readLines(f); l <- l[!grepl("^#", l) & nzchar(l)]
  ids <- sub("\t.*$", "", l); ids[ids %in% rownames(otu)]
}
asv2fun <- do.call(rbind, lapply(sub_files, function(f)
  data.frame(function_group = sub("\\.(tsv|txt)$", "", basename(f)), ASV_ID = read_ids(f))))
asv2fun <- asv2fun[asv2fun$function_group %in% rownames(fm), ]

n_assigned <- length(unique(asv2fun$ASV_ID))
message(sprintf("FAPROTAX: %d of %d ASVs assigned (%.1f %%); %d functions with non-zero abundance  (manuscript: 941 of 2,752, 34.2 %%; 41 functions)",
                n_assigned, nrow(otu), 100 * n_assigned / nrow(otu), nrow(fm)))

# Per-genus contribution = reads of that genus' ASVs assigned to the function
asv2fun$Genus <- lowest_label(tax[asv2fun$ASV_ID, ], "Genus")
asv2fun$reads <- rowSums(otu)[asv2fun$ASV_ID]
contrib <- aggregate(reads ~ function_group + Genus, data = asv2fun, FUN = sum)
contrib$share_of_function <- contrib$reads / ave(contrib$reads, contrib$function_group, FUN = sum)
contrib <- contrib[order(contrib$function_group, -contrib$reads), ]
top3 <- tapply(seq_len(nrow(contrib)), contrib$function_group, function(i) {
  i <- head(i, 3)
  paste(sprintf("%s (%.0f %%)", contrib$Genus[i], 100 * contrib$share_of_function[i]), collapse = "; ")
})

s4 <- data.frame(function_group = rownames(fm),
                 total_predicted_abundance = rowSums(fm),
                 mean_per_sample = round(rowMeans(fm), 1),
                 n_ASVs_assigned = as.vector(table(factor(asv2fun$function_group, levels = rownames(fm)))),
                 top_contributing_genera = top3[rownames(fm)],
                 fm, check.names = FALSE, row.names = NULL)
s4 <- s4[order(-s4$total_predicted_abundance), ]
print(head(s4[, 1:5], 12))
write_supp(s4, "Table_S4_FAPROTAX_functions")
write_supp(contrib, "Table_S4b_FAPROTAX_genus_contributions")
saveRDS(contrib, file.path(dir_rds, "faprotax_genus_contrib.rds"))

# ---- Heatmap: top 20 functions -------------------------------------------------
top <- head(s4$function_group, par$faprotax_top_fun)
mat <- log10(fm[top, , drop = FALSE] + 1)
rownames(mat) <- gsub("_", " ", rownames(mat))
ann <- data.frame(Dataset = meta[colnames(mat), "dataset"], row.names = colnames(mat))
save_fig("Fig_FAPROTAX_heatmap_top20", width = 9, height = 7, draw = function() {
  pheatmap(mat, annotation_col = ann,
           annotation_colors = list(Dataset = setNames(c("#0072B2", "#D55E00"), datasets)),
           color = colorRampPalette(c("#fff5eb", "#fd8d3c", "#7f2704"))(100),
           clustering_method = "complete", border_color = NA, fontsize_row = 9, fontsize_col = 8,
           main = "FAPROTAX predicted functions: log10(read count + 1)")
})
message("\nDone: 08_faprotax_post.R")
