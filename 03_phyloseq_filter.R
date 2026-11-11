# =============================================================================
# 03_phyloseq_filter.R
# Build the phyloseq object, remove non-bacterial / chloroplast / mitochondrial
# ASVs (Methods 2.4), apply the >= 10,000-read filter and rarefy (Methods 2.5).
# Prints every count reported in Table S0a / Figure 1 so they can be checked.
#
# Output: results/rds/ps_bact.rds (52 samples), ps_depth.rds (21 samples),
#         ps_rare.rds (21 samples, rarefied), core FASTA written in 06.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(Biostrings))

seqtab <- readRDS(file.path(dir_rds, "seqtab_nochim.rds"))
taxa   <- readRDS(file.path(dir_rds, "taxa.rds"))
meta   <- read.csv(file_metadata, stringsAsFactors = FALSE)
stopifnot(all(c("sample_id", "dataset") %in% names(meta)))
meta <- meta[meta$sample_id %in% rownames(seqtab), , drop = FALSE]
rownames(meta) <- meta$sample_id
missing_meta <- setdiff(rownames(seqtab), meta$sample_id)
if (length(missing_meta)) stop("Samples missing from metadata.csv: ", paste(missing_meta, collapse = ", "))

# Short ASV IDs; full sequences kept in refseq()
seqs <- colnames(seqtab)
asv_ids <- sprintf("ASV_%04d", seq_along(seqs))
colnames(seqtab) <- asv_ids
rownames(taxa)   <- asv_ids[match(rownames(taxa), seqs)]
dna <- DNAStringSet(seqs); names(dna) <- asv_ids

ps <- phyloseq(otu_table(seqtab, taxa_are_rows = FALSE),
               tax_table(taxa), sample_data(meta[rownames(seqtab), , drop = FALSE]),
               refseq(dna))

counts <- list()
counts$after_dada2 <- c(samples = nsamples(ps), asvs = ntaxa(ps))

# ---- Taxonomic filtering ----------------------------------------------------
log_step("Taxonomic filtering")
tx <- as.data.frame(tax_table(ps), stringsAsFactors = FALSE)
drop <- is.na(tx$Kingdom) | tx$Kingdom != "Bacteria" |
        (!is.na(tx$Order)  & tx$Order  == "Chloroplast") |
        (!is.na(tx$Family) & tx$Family == "Mitochondria")
message(sprintf("  removed: %d non-bacterial/unassigned kingdom, %d chloroplast, %d mitochondria",
                sum(is.na(tx$Kingdom) | tx$Kingdom != "Bacteria"),
                sum(!is.na(tx$Order)  & tx$Order  == "Chloroplast"),
                sum(!is.na(tx$Family) & tx$Family == "Mitochondria")))
ps_bact <- prune_taxa(taxa_names(ps)[!drop], ps)
counts$after_taxfilter <- c(samples = nsamples(ps_bact), asvs = ntaxa(ps_bact))
message(sprintf("  %d ASVs removed; %d bacterial ASVs in %d samples  (manuscript: 254 removed; 2,752 ASVs; 52 samples)",
                sum(drop), ntaxa(ps_bact), nsamples(ps_bact)))

# ---- Read-depth filtering -----------------------------------------------------
log_step("Read-depth filter (>= ", par$min_reads, " reads)")
depth <- sample_sums(ps_bact)
ps_depth <- prune_samples(depth >= par$min_reads, ps_bact)
n_zero <- sum(taxa_sums(ps_depth) == 0)
message(sprintf("  %d samples excluded; %d retained (%s)  (manuscript: 31 excluded; 21 retained)",
                sum(depth < par$min_reads), nsamples(ps_depth),
                paste(names(table(sample_data(ps_depth)$dataset)),
                      table(sample_data(ps_depth)$dataset), sep = ": ", collapse = ", ")))
message("  ASVs with zero counts in the retained samples: ", n_zero,
        "  (if > 0, the ASV count used downstream is lower than 2,752 - update the Methods)")
ps_depth <- prune_taxa(taxa_sums(ps_depth) > 0, ps_depth)
counts$after_depth <- c(samples = nsamples(ps_depth), asvs = ntaxa(ps_depth))

# ---- Rarefaction ------------------------------------------------------------
log_step("Rarefaction")
depth_rare <- min(sample_sums(ps_depth))
ps_rare <- rarefy_even_depth(ps_depth, sample.size = depth_rare, rngseed = par$rare_seed,
                             replace = par$rare_replace, trimOTUs = TRUE, verbose = TRUE)
counts$after_rarefaction <- c(samples = nsamples(ps_rare), asvs = ntaxa(ps_rare))
otu_r <- as(otu_table(ps_rare), "matrix"); if (!taxa_are_rows(ps_rare)) otu_r <- t(otu_r)
message(sprintf("  depth = %d reads/sample (manuscript: 10,862); %d ASVs; singletons = %d, doubletons = %d (manuscript: 214, 250)",
                depth_rare, ntaxa(ps_rare), sum(rowSums(otu_r) == 1), sum(rowSums(otu_r) == 2)))

flow <- data.frame(step = names(counts),
                   samples = sapply(counts, `[`, "samples"),
                   asvs = sapply(counts, `[`, "asvs"), row.names = NULL)
print(flow)
write.csv(flow, file.path(dir_tables, "Table_S0a_sample_flow_summary.csv"), row.names = FALSE)

saveRDS(ps,       file.path(dir_rds, "ps_all.rds"))
saveRDS(ps_bact,  file.path(dir_rds, "ps_bact.rds"))
saveRDS(ps_depth, file.path(dir_rds, "ps_depth.rds"))
saveRDS(ps_rare,  file.path(dir_rds, "ps_rare.rds"))
message("\nDone: 03_phyloseq_filter.R")
