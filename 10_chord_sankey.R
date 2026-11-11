# =============================================================================
# 10_chord_sankey.R
# Integrated taxon -> predicted function figures for the candidate genera.
#
# Link weight = the genus' share of each function's predicted total (0-1),
# computed separately within FAPROTAX (read counts) and PICRUSt2 (gene-copy
# abundance).  Shares are used because the two tools' raw values are on
# different scales; say so in the figure legends.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages({ library(circlize); library(ggalluvial) })

fx <- readRDS(file.path(dir_rds, "faprotax_genus_contrib.rds"))
pc <- readRDS(file.path(dir_rds, "picrust_genus_contrib.rds"))

fx <- fx[fx$function_group %in% names(faprotax_focus) & fx$Genus %in% names(candidate_genera), ]
pc <- pc[pc$KO %in% picrust_focus$KO & pc$Genus %in% names(candidate_genera), ]

links <- rbind(
  data.frame(genus = candidate_genera[fx$Genus],
             fun = paste0(gsub("_", " ", fx$function_group), " [F]"),
             category = faprotax_focus[fx$function_group],
             weight = fx$share_of_function, source = "FAPROTAX"),
  data.frame(genus = candidate_genera[pc$Genus],
             fun = paste0(pc$KO, " [P]"),
             category = picrust_focus$category[match(pc$KO, picrust_focus$KO)],
             weight = pc$share_of_function, source = "PICRUSt2"))
links <- links[links$weight >= 0.01, ]          # drop links < 1 % share
rownames(links) <- NULL
write.csv(links, file.path(dir_tables, "chord_sankey_links.csv"), row.names = FALSE)
message("Links: ", nrow(links), " (", paste(names(table(links$source)), table(links$source),
                                            sep = " = ", collapse = ", "), ")")

# ---- Chord diagram -----------------------------------------------------------
genus_cols <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00"),
                       unname(candidate_genera))
fun_cols <- setNames(rep("grey70", length(unique(links$fun))), unique(links$fun))
save_fig("Fig_chord_taxon_function", width = 8, height = 8, draw = function() {
  circos.clear()
  circos.par(gap.after = c(rep(2, length(unique(links$genus)) - 1), 10,
                           rep(2, length(unique(links$fun)) - 1), 10))
  chordDiagram(links[, c("genus", "fun", "weight")],
               order = c(intersect(names(genus_cols), links$genus), sort(unique(links$fun))),
               grid.col = c(genus_cols, fun_cols), transparency = 0.35,
               annotationTrack = "grid", preAllocateTracks = list(track.height = 0.25))
  circos.trackPlotRegion(track.index = 1, bg.border = NA, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.7)
  })
  title(sub = "[F] FAPROTAX, [P] PICRUSt2; ribbon width = share of predicted function", cex.sub = 0.8)
  circos.clear()
})

# ---- Sankey: genus -> function -> category -----------------------------------
links$genus <- factor(links$genus, levels = intersect(unname(candidate_genera), links$genus))
p_sk <- ggplot(links, aes(axis1 = genus, axis2 = fun, axis3 = category, y = weight)) +
  geom_alluvium(aes(fill = genus), width = 1/6, alpha = 0.7) +
  geom_stratum(width = 1/6, fill = "grey95", colour = "grey40") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6) +
  scale_x_discrete(limits = c("Genus", "Predicted function", "Category"), expand = c(0.08, 0.08)) +
  scale_fill_manual(values = genus_cols) +
  labs(y = "Summed share of predicted function", fill = "Genus",
       caption = "[F] FAPROTAX, [P] PICRUSt2. All links are computational inferences.") +
  theme_minimal(base_size = 10) + theme(panel.grid = element_blank())
save_fig("Fig_sankey_taxon_function", p_sk, width = 10, height = 7)
message("\nDone: 10_chord_sankey.R")
