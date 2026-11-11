# =============================================================================
# 11_cooccurrence_network.R
# Genus-level Spearman co-occurrence network (Methods 2.8): top 30 genera by
# mean relative abundance, |r| >= 0.6, undirected weighted igraph network,
# isolated nodes removed, Fruchterman-Reingold layout.
#
# The Methods report no FDR correction; Table S8 includes raw and BH-adjusted
# p-values so reviewers can see both.  Treat the network as exploratory.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(igraph))

ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
gm <- 100 * rel_abund(glom_named(ps_rare, "Genus"))
if (par$net_exclude_unclassified) gm <- gm[!is_unclassified(rownames(gm)), , drop = FALSE]
top <- names(sort(rowMeans(gm), decreasing = TRUE))[seq_len(min(par$net_top_genera, nrow(gm)))]
x <- t(gm[top, , drop = FALSE])                     # samples x genera
colnames(x)[colnames(x) == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium"] <- "Rhizobium group"

log_step("Spearman correlations: ", ncol(x), " genera x ", nrow(x), " samples")
pairs <- t(combn(colnames(x), 2))
res <- do.call(rbind, lapply(seq_len(nrow(pairs)), function(i) {
  ct <- suppressWarnings(cor.test(x[, pairs[i, 1]], x[, pairs[i, 2]], method = "spearman", exact = FALSE))
  data.frame(from = pairs[i, 1], to = pairs[i, 2], r = unname(ct$estimate), p = ct$p.value)
}))
res$p_BH <- p.adjust(res$p, method = "BH")
edges <- res[!is.na(res$r) & abs(res$r) >= par$net_r_cutoff, ]
edges$sign <- ifelse(edges$r > 0, "positive", "negative")
message(sprintf("  pairs tested: %d; |r| >= %.1f: %d (positive %d, negative %d); of these p < 0.05: %d; BH q < 0.05: %d",
                nrow(res), par$net_r_cutoff, nrow(edges), sum(edges$r > 0), sum(edges$r < 0),
                sum(edges$p < 0.05), sum(edges$p_BH < 0.05)))
message("  (manuscript reports 21 nodes and 29 edges, and separately '76 correlations' - ",
        "use these counts to reconcile the two numbers in the text)")

g <- graph_from_data_frame(edges[, c("from", "to", "r", "sign")], directed = FALSE,
                           vertices = data.frame(name = colnames(x)))
E(g)$weight <- abs(E(g)$r)
g <- delete_vertices(g, which(degree(g) == 0))
message(sprintf("  network after removing isolated nodes: %d nodes, %d edges", vcount(g), ecount(g)))

nodes <- data.frame(genus = V(g)$name, degree = degree(g),
                    betweenness = round(betweenness(g, weights = NA), 2),
                    closeness = round(closeness(g, weights = NA), 4),
                    mean_rel_abund_pct = round(colMeans(x)[V(g)$name], 4), row.names = NULL)
nodes <- nodes[order(-nodes$degree), ]
print(head(nodes, 8))

edges_out <- edges[edges$from %in% V(g)$name & edges$to %in% V(g)$name, ]
edges_out[, c("r", "p", "p_BH")] <- lapply(edges_out[, c("r", "p", "p_BH")], signif, 4)
write_supp(edges_out, "Table_S8a_network_edges")
write_supp(nodes,     "Table_S8b_network_nodes")

if (ecount(g) == 0) {
  saveRDS(g, file.path(dir_rds, "network.rds"))
  message("No correlations pass |r| >= ", par$net_r_cutoff, "; network figure skipped.")
  quit(save = "no")
}
set.seed(par$rare_seed)
lay <- layout_with_fr(g, weights = E(g)$weight)
save_fig("Fig_cooccurrence_network", width = 8, height = 8, draw = function() {
  par(mar = c(1, 1, 1, 1))
  plot(g, layout = lay,
       vertex.size = 6 + 3 * degree(g), vertex.color = "#9ecae1", vertex.frame.color = "grey30",
       vertex.label.cex = 0.75, vertex.label.color = "black", vertex.label.family = "sans",
       edge.width = 1 + 4 * (E(g)$weight - par$net_r_cutoff) / (1 - par$net_r_cutoff),
       edge.color = ifelse(E(g)$sign == "positive", "#0072B2", "#D55E00"))
  legend("bottomleft", c("Positive", "Negative"), col = c("#0072B2", "#D55E00"),
         lwd = 2, bty = "n", cex = 0.85, title = sprintf("Spearman |r| >= %.1f", par$net_r_cutoff))
})
saveRDS(g, file.path(dir_rds, "network.rds"))
message("\nDone: 11_cooccurrence_network.R")
