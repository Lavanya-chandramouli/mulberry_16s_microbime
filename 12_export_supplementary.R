# =============================================================================
# 12_export_supplementary.R
# Combines Tables S2-S8 into one Excel workbook (one sheet per table) and
# records package versions for the Methods / data-availability statement.
# Each table also exists as a CSV in results/tables/.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(openxlsx))

order <- c(
  "Table_S2_alpha_diversity"                = "S2 Alpha diversity",
  "Table_S3a_core_ASVs"                     = "S3a Core ASVs",
  "Table_S3b_core_genera_per_dataset"       = "S3b Core genera by dataset",
  "Table_S3c_core_sensitivity"              = "S3c Core sensitivity",
  "Table_S4_FAPROTAX_functions"             = "S4 FAPROTAX functions",
  "Table_S4b_FAPROTAX_genus_contributions"  = "S4b FAPROTAX contributors",
  "Table_S5a_PICRUSt2_focal_KOs"            = "S5a PICRUSt2 focal KOs",
  "Table_S5b_PICRUSt2_genus_contributions"  = "S5b PICRUSt2 contributors",
  "Table_S5c_PICRUSt2_NSTI_per_ASV"         = "S5c NSTI per ASV",
  "Table_S6_phylum_relative_abundance"      = "S6 Phylum abundance",
  "Table_S7_genus_relative_abundance"       = "S7 Genus abundance",
  "Table_S8a_network_edges"                 = "S8a Network edges",
  "Table_S8b_network_nodes"                 = "S8b Network nodes")

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7", border = "bottom")
for (nm in names(order)) {
  f <- file.path(dir_rds, paste0("supp_", nm, ".rds"))
  if (!file.exists(f)) { message("  missing (script not run yet): ", nm); next }
  df <- readRDS(f)
  sh <- substr(order[[nm]], 1, 31)
  addWorksheet(wb, sh)
  writeData(wb, sh, df, headerStyle = hdr)
  freezePane(wb, sh, firstRow = TRUE)
  setColWidths(wb, sh, cols = 1:min(ncol(df), 12), widths = "auto")
}
out <- file.path(dir_tables, "Supplementary_Tables_S2-S8.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)
message("Workbook: ", out)

writeLines(capture.output(sessionInfo()), file.path(dir_results, "sessionInfo.txt"))
message("Package versions: results/sessionInfo.txt")
