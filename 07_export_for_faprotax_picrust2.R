# =============================================================================
# 07_export_for_faprotax_picrust2.R
# Writes the input files for the two external (non-R) tools and prints the
# exact commands to run them.  Run the commands, then continue with 08 and 09.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(Biostrings))

ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
otu <- as(otu_table(ps_rare), "matrix"); if (!taxa_are_rows(ps_rare)) otu <- t(otu)
tax <- as.data.frame(tax_table(ps_rare), stringsAsFactors = FALSE)

# ---- FAPROTAX input: ASV x sample counts + semicolon-delimited taxonomy -------
ranks <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
tax_str <- apply(tax[rownames(otu), ranks], 1, function(r) paste(r[!is.na(r)], collapse = ";"))
f_in <- file.path(dir_faprotax, "faprotax_input.tsv")
con <- file(f_in, "w")
writeLines(paste(c("#ASV_ID", colnames(otu), "taxonomy"), collapse = "\t"), con)
write.table(data.frame(rownames(otu), otu, tax_str, check.names = FALSE),
            con, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
close(con)

# ---- PICRUSt2 input: FASTA + count table ---------------------------------------
writeXStringSet(refseq(ps_rare)[rownames(otu)], file.path(dir_picrust, "asv_seqs.fna"))
f_counts <- file.path(dir_picrust, "asv_counts.tsv")
con <- file(f_counts, "w")
writeLines(paste(c("ASV_ID", colnames(otu)), collapse = "\t"), con)
write.table(otu, con, sep = "\t", quote = FALSE, col.names = FALSE)
close(con)

message("Inputs written: ", nrow(otu), " ASVs x ", ncol(otu), " samples\n")
cat('
## ---- Run FAPROTAX v1.2.12 (Python 3) ------------------------------------
cd results/faprotax
python3 /path/to/FAPROTAX_1.2.12/collapse_table.py \\
  -i faprotax_input.tsv -o faprotax_functions.tsv \\
  -g /path/to/FAPROTAX_1.2.12/FAPROTAX.txt \\
  -c "#" -d "taxonomy" --omit_columns 0 --column_names_are_in last_comment_line \\
  -r faprotax_report.txt -n columns_after_collapsing \\
  --out_sub_tables_dir sub_tables -v
cd ../..

## ---- Run PICRUSt2 v2.6.3 (conda env) ------------------------------------
conda activate picrust2
picrust2_pipeline.py -s results/picrust2/asv_seqs.fna -i results/picrust2/asv_counts.tsv \\
  -o results/picrust2/out -p 8 --stratified --max_nsti 2
add_descriptions.py -m KO \\
  -i results/picrust2/out/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz \\
  -o results/picrust2/out/KO_metagenome_out/pred_metagenome_unstrat_descrip.tsv.gz
')
message("\nDone: 07_export_for_faprotax_picrust2.R")
