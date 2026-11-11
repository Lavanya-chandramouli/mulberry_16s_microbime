# Paper 1 analysis scripts: mulberry-associated bacterial microbiome

These scripts reanalyse two public 16S rRNA (V3–V4) amplicon datasets, SRP411931 and PRJNA822613, in R 4.3. The parameters are the ones given in the manuscript Methods (Sections 2.2–2.9).

> **Before you deposit these on GitHub or Zenodo:** these scripts were rebuilt from the parameters in the Methods, not copied from the code that produced the manuscript's numbers. Run them on the real data and check that the key counts match (see the checklist below). If something doesn't match, use the original code for the deposit, or update the manuscript to the numbers these scripts produce. The deposited code has to be the code behind the published results.

## Run order

| Step | Script | Does | Manuscript output |
|---|---|---|---|
| 0 | `00_config.R` | Paths, parameters and helpers. **Edit the PATHS block first.** | — |
| 1 | `01_qc_cutadapt.sh` | FastQC, MultiQC, Cutadapt primer removal (shell) | Sec. 2.2, 99.7 % retention |
| 2 | `02_dada2.R` | Filtering, error models, denoising, merging (per dataset); merge; chimeras; SILVA taxonomy | Table S0a, 3,006 ASVs, 52 samples |
| 3 | `03_phyloseq_filter.R` | Removes non-bacterial, chloroplast and mitochondrial ASVs; ≥ 10,000-read filter; rarefaction (seed 123) | 2,752 ASVs, 21 samples, depth 10,862 |
| 4 | `04_alpha_beta_diversity.R` | Alpha indices; Bray–Curtis PCoA; adonis2; betadisper | **Table S2**, Table 2, violin plot, PCoA |
| 5 | `05_composition_heatmap.R` | Phylum and genus relative abundance; top-30-genera heatmap | **Tables S6, S7**, heatmap |
| 6 | `06_core_microbiome.R` | Core ASVs (0.01 % detection, 80 % prevalence); per-dataset check; 75 %/85 % sensitivity; core FASTA | **Table S3a–c**, bubble plot, Supplementary File S1 |
| 7 | `07_export_for_faprotax_picrust2.R` | Writes tool inputs and **prints the FAPROTAX and PICRUSt2 commands** | — |
| — | *run the printed commands* | FAPROTAX v1.2.12, PICRUSt2 v2.6.3 | — |
| 8 | `08_faprotax_post.R` | Function table, assignment rate, contributing genera, top-20 heatmap | **Table S4**, FAPROTAX heatmap |
| 9 | `09_picrust2_post.R` | NSTI, focal KOs with KEGG descriptions, contributing genera | **Table S5a–c** |
| 10 | `10_chord_sankey.R` | Links between candidate genera and predicted functions | Chord diagram, Sankey |
| 11 | `11_cooccurrence_network.R` | Spearman network of the top 30 genera, \|r\| ≥ 0.6, igraph | **Table S8a–b**, network |
| 12 | `12_export_supplementary.R` | Combines Tables S2–S8 into one `.xlsx`; writes `sessionInfo.txt` | `Supplementary_Tables_S2-S8.xlsx` |

Run the steps in two batches:

```bash
bash 01_qc_cutadapt.sh
Rscript run_all.R          # steps 02–07
# run the FAPROTAX + PICRUSt2 commands that step 07 prints
Rscript run_all.R post     # steps 08–12
```

Outputs go to `results/`: `tables/` holds each table as a CSV plus the combined workbook, `figures/` holds each figure as a PDF plus a 300-dpi PNG, and `rds/` holds intermediate objects.

## Inputs you provide

- `data/raw/<DATASET>/<RUN>_1.fastq.gz`, `_2.fastq.gz`: SRA downloads (`fasterq-dump --split-files`).
- `data/metadata.csv` with the columns `sample_id` (the SRA run accession) and `dataset` (`SRP411931` or `PRJNA822613`). This becomes Table S1.
- `ref/silva_nr99_v138.1_train_set.fa.gz`: the DADA2-formatted SILVA 138.1 training set (Zenodo record 4587955).
- FAPROTAX 1.2.12 (`collapse_table.py`, `FAPROTAX.txt`) and a PICRUSt2 2.6.3 conda environment.

## R packages

From Bioconductor: dada2 1.30, phyloseq 1.46, Biostrings, and microbiome 1.24. From CRAN: vegan 2.6, ggplot2 3.4, pheatmap 1.0.12, circlize 0.4.x, ggalluvial 0.12, igraph 1.5+, openxlsx.

If `microbiome` isn't installed, `06_core_microbiome.R` uses the identical base-R rule (relative abundance > detection in > 80 % of samples). When the package *is* installed, the script also runs that rule and stops if the two disagree.

## Checklist: compare with the manuscript

Each script prints its result next to the manuscript value (lines marked `(manuscript: …)`). Check these in particular:

| Check | Where it prints | Why it matters |
|---|---|---|
| Samples after DADA2 = 52 | step 02 | Section 2.5 says "53". Use whichever the script reports. |
| ASVs after the depth filter and rarefaction | step 03 | The manuscript says 2,752 at every stage. ASVs found only in the 31 excluded samples drop out, so the count used downstream is probably lower. Report the real number. |
| Rarefaction `replace` | `00_config.R` | phyloseq samples with replacement by default. State which you used. |
| 80 core ASVs; *Pseudomonas*, *Bacillus*, *Bdellovibrio* at 100 % in both datasets | step 06 | Core definition |
| 941/2,752 FAPROTAX-assigned ASVs; 41 functions | step 08 | Table 3 |
| KEGG description vs. manuscript label | step 09 | K02014 is an iron-complex outer-membrane receptor (uptake, not biosynthesis). K01580 is glutamate decarboxylase (not IAA). The script warns if a KO total equals the rarefaction depth or if two KOs have identical totals. |
| Edges with \|r\| ≥ 0.6; nodes/edges after removing isolated nodes | step 11 | The text says both "76 correlations" and "29 edges". Removing isolated nodes doesn't remove edges, so reconcile these. Table S8 lists raw and BH-adjusted p-values. |

## Testing done

All scripts ran end-to-end on the R 4.3 package versions above:

- `02_dada2.R` was tested on the example reads bundled with DADA2.
- Steps 03–12 were tested on a synthetic 70-sample dataset.
- FAPROTAX and PICRUSt2 couldn't be installed in the test environment, so steps 08–10 were tested on mock outputs written in those tools' documented formats.

If a later PICRUSt2 release renames its output files, adjust the file patterns at the top of `09_picrust2_post.R`. Check `collapse_table.py -h` if FAPROTAX rejects an option.
