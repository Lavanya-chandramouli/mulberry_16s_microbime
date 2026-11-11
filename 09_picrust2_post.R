# =============================================================================
# 09_picrust2_post.R
# Reads PICRUSt2 v2.6.3 output: NSTI summary, predicted abundance of the
# focal KOs (Table 4 / Table S5), their official KEGG descriptions, and the
# genera contributing to each (stratified output).
# All values are computationally inferred gene-copy abundances, not measurements.
# =============================================================================
source("00_config.R")

dir_out <- file.path(dir_picrust, "out")
find1 <- function(pattern) {
  f <- list.files(dir_out, recursive = TRUE, full.names = TRUE)
  f <- f[grepl(pattern, f)]
  if (!length(f)) stop("Not found in ", dir_out, ": ", pattern)
  f[1]
}
ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
tax <- as.data.frame(tax_table(ps_rare), stringsAsFactors = FALSE)
depth <- unique(sample_sums(ps_rare))

# ---- NSTI ----------------------------------------------------------------------
nsti <- read.delim(find1("marker_predicted_and_nsti\\.tsv(\\.gz)?$"), check.names = FALSE)
names(nsti)[1] <- "ASV_ID"
nsti_col <- grep("NSTI", names(nsti), value = TRUE)[1]
nsti <- data.frame(ASV_ID = nsti$ASV_ID, NSTI = nsti[[nsti_col]])
nsti$Genus <- lowest_label(tax[nsti$ASV_ID, ], "Genus")
ns <- summary(nsti$NSTI)
message(sprintf("NSTI: n = %d, range %.4f-%.4f, mean %.3f, median %.3f; %d (%.1f %%) < 0.1; %d > 2  (manuscript: 0.0001-0.1319, mean 0.030, 99.4 %% < 0.1)",
                nrow(nsti), ns["Min."], ns["Max."], ns["Mean"], ns["Median"],
                sum(nsti$NSTI < 0.1), 100 * mean(nsti$NSTI < 0.1), sum(nsti$NSTI > 2)))
message("ASVs in input but not in PICRUSt2 output: ", length(setdiff(taxa_names(ps_rare), nsti$ASV_ID)))

# ---- KO table (with KEGG descriptions) ---------------------------------------
f_ko <- tryCatch(find1("KO_metagenome_out/pred_metagenome_unstrat_descrip\\.tsv(\\.gz)?$"),
                 error = function(e) find1("KO_metagenome_out/pred_metagenome_unstrat\\.tsv(\\.gz)?$"))
ko <- read.delim(f_ko, check.names = FALSE, stringsAsFactors = FALSE)
names(ko)[1] <- "KO"
has_desc <- "description" %in% names(ko)
smp <- intersect(names(ko), sample_names(ps_rare))
message("KOs predicted (non-zero): ", sum(rowSums(ko[, smp]) > 0), "  (manuscript: 5,072)")

foc <- merge(picrust_focus, ko, by = "KO", all.x = TRUE)
foc$KEGG_description <- if (has_desc) foc$description else "run add_descriptions.py (see 07)"
foc$total_predicted_abundance <- rowSums(foc[, smp, drop = FALSE])
foc$mean_per_sample <- round(rowMeans(foc[, smp, drop = FALSE]), 1)

# ---- Sanity checks: these caught real problems in the draft -----------------
dup <- foc$KO[duplicated(foc$total_predicted_abundance) | duplicated(foc$total_predicted_abundance, fromLast = TRUE)]
if (length(dup)) warning("Identical totals for KOs: ", paste(dup, collapse = ", "),
                         " - check they are not copies of the same column.")
eq_depth <- foc$KO[foc$total_predicted_abundance %in% c(depth, depth * length(smp))]
if (length(eq_depth)) warning("Total equals the rarefaction depth for: ", paste(eq_depth, collapse = ", "),
                              " - likely an artefact, not a prediction.")
if (has_desc) {
  message("\nCheck the manuscript label against KEGG for each focal KO:")
  op <- options(width = 200); print(foc[, c("KO", "label", "KEGG_description")], right = FALSE); options(op)
}

# ---- Contributing genera (stratified output) ---------------------------------
contrib_f <- find1("KO_metagenome_out/pred_metagenome_contrib\\.tsv(\\.gz)?$")
ct <- read.delim(contrib_f, stringsAsFactors = FALSE, check.names = FALSE)
names(ct)[names(ct) == "function"] <- "KO"   # column name in PICRUSt2 >= 2.5
ct <- ct[ct$KO %in% picrust_focus$KO, ]
ct$Genus <- lowest_label(tax[ct$taxon, ], "Genus")
gc <- aggregate(taxon_function_abun ~ KO + Genus, data = ct, FUN = sum)
gc$share_of_function <- gc$taxon_function_abun / ave(gc$taxon_function_abun, gc$KO, FUN = sum)
gc <- gc[order(gc$KO, -gc$taxon_function_abun), ]
top3 <- tapply(seq_len(nrow(gc)), gc$KO, function(i) {
  i <- head(i, 3)
  paste(sprintf("%s (%.0f %%)", gc$Genus[i], 100 * gc$share_of_function[i]), collapse = "; ")
})
foc$top_contributing_genera <- top3[foc$KO]

s5 <- foc[, c("KO", "label", "KEGG_description", "category", "total_predicted_abundance",
              "mean_per_sample", "top_contributing_genera", smp)]
write_supp(s5, "Table_S5a_PICRUSt2_focal_KOs")
write_supp(gc, "Table_S5b_PICRUSt2_genus_contributions")
write_supp(nsti[order(nsti$NSTI), ], "Table_S5c_PICRUSt2_NSTI_per_ASV")
saveRDS(gc, file.path(dir_rds, "picrust_genus_contrib.rds"))
message("\nDone: 09_picrust2_post.R")
