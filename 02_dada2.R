# =============================================================================
# 02_dada2.R
# ASV inference with DADA2 v1.30 (Methods 2.3) and taxonomy with SILVA v138.1
# (Methods 2.4). Each dataset is denoised separately (separate error models,
# as the runs differ), then the sequence tables are merged before chimera
# removal and taxonomy.
#
# Output: results/rds/seqtab_nochim.rds, results/rds/taxa.rds,
#         results/tables/Table_S0a_read_tracking_per_sample.csv
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(dada2))
set.seed(par$rare_seed)

getN <- function(x) sum(getUniques(x))
seqtabs <- list(); tracks <- list()

for (ds in datasets) {
  log_step("DADA2: ", ds)
  path <- file.path(dir_trimmed, ds)
  fnFs <- sort(list.files(path, pattern = paste0(fwd_suffix, "$"), full.names = TRUE))
  fnRs <- sort(list.files(path, pattern = paste0(rev_suffix, "$"), full.names = TRUE))
  stopifnot(length(fnFs) > 0, length(fnFs) == length(fnRs))
  samples <- sub(fwd_suffix, "", basename(fnFs), fixed = TRUE)
  stopifnot(identical(samples, sub(rev_suffix, "", basename(fnRs), fixed = TRUE)))

  filtFs <- file.path(path, "filtered", paste0(samples, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(samples, "_R_filt.fastq.gz"))
  names(filtFs) <- names(filtRs) <- samples

  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,
                       truncLen = par$truncLen, maxEE = par$maxEE,
                       truncQ = par$truncQ, maxN = par$maxN,
                       rm.phix = par$rm.phix, compress = TRUE,
                       multithread = par$threads)
  rownames(out) <- samples

  # Samples with no reads surviving filtering cannot be denoised.
  # (If one dataset was sequenced 2x250, every read is shorter than
  #  truncLen = 270 and the whole sample drops out here; check this table.)
  keep <- out[, "reads.out"] > 0
  if (any(!keep)) message("  ", sum(!keep), " sample(s) lost all reads in filterAndTrim: ",
                          paste(samples[!keep], collapse = ", "))
  filtFs <- filtFs[keep]; filtRs <- filtRs[keep]

  errF <- learnErrors(filtFs, multithread = par$threads)
  errR <- learnErrors(filtRs, multithread = par$threads)
  save_fig(paste0("QC_errors_", ds, "_F"), plotErrors(errF, nominalQ = TRUE), 8, 6)

  derepFs <- derepFastq(filtFs); derepRs <- derepFastq(filtRs)
  dadaFs  <- dada(derepFs, err = errF, multithread = par$threads)
  dadaRs  <- dada(derepRs, err = errR, multithread = par$threads)
  mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs,
                        minOverlap = par$minOverlap, maxMismatch = par$maxMismatch)
  if (inherits(dadaFs, "dada")) {     # only one sample: DADA2 returns bare objects
    dadaFs  <- setNames(list(dadaFs),  names(filtFs))
    dadaRs  <- setNames(list(dadaRs),  names(filtFs))
    mergers <- setNames(list(mergers), names(filtFs))
  }

  st <- makeSequenceTable(mergers)
  seqtabs[[ds]] <- st

  tr <- data.frame(sample_id = samples, dataset = ds,
                   input = out[, "reads.in"], filtered = out[, "reads.out"],
                   denoisedF = NA, denoisedR = NA, merged = NA)
  ok <- tr$sample_id %in% names(dadaFs)
  tr$denoisedF[ok] <- sapply(dadaFs, getN)[tr$sample_id[ok]]
  tr$denoisedR[ok] <- sapply(dadaRs, getN)[tr$sample_id[ok]]
  tr$merged[ok]    <- sapply(mergers, getN)[tr$sample_id[ok]]
  tracks[[ds]] <- tr
}

log_step("Merge datasets, remove chimeras")
seqtab <- mergeSequenceTables(tables = seqtabs)
message("  merged table: ", nrow(seqtab), " samples x ", ncol(seqtab), " ASVs")
message("  amplicon length distribution:"); print(table(nchar(getSequences(seqtab))))

seqtab_nochim <- removeBimeraDenovo(seqtab, method = par$chimera_method,
                                    multithread = par$threads, verbose = TRUE)
message(sprintf("  after chimera removal: %d ASVs; %.1f %% of merged reads retained",
                ncol(seqtab_nochim), 100 * sum(seqtab_nochim) / sum(seqtab)))

track <- do.call(rbind, tracks)
track$nonchim <- NA
idx <- match(rownames(seqtab_nochim), track$sample_id)
track$nonchim[idx] <- rowSums(seqtab_nochim)
track$retained_after_dada2 <- !is.na(track$nonchim) & track$nonchim > 0
rownames(track) <- NULL
write.csv(track, file.path(dir_tables, "Table_S0a_read_tracking_per_sample.csv"), row.names = FALSE)
message("  samples retained after DADA2: ", sum(track$retained_after_dada2),
        " of ", nrow(track), "  (manuscript: 52 of 70)")

log_step("Taxonomy: SILVA v138.1, minBoot = ", par$minBoot)
taxa <- assignTaxonomy(seqtab_nochim, file_silva, minBoot = par$minBoot,
                       multithread = par$threads)
# Species-level assignment deliberately not attempted (Methods 2.4)

saveRDS(seqtab_nochim, file.path(dir_rds, "seqtab_nochim.rds"))
saveRDS(taxa,          file.path(dir_rds, "taxa.rds"))
saveRDS(track,         file.path(dir_rds, "track.rds"))
message("\nDone: 02_dada2.R")
