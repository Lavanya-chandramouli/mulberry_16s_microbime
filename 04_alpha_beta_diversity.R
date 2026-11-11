# =============================================================================
# 04_alpha_beta_diversity.R
# Alpha diversity (Table 2, Table S2, Figure: violin plots) and beta diversity
# (Bray-Curtis PCoA, PERMANOVA with adonis2, betadisper)  - Methods 2.5.
# =============================================================================
source("00_config.R")
suppressPackageStartupMessages(library(vegan))

ps_rare <- readRDS(file.path(dir_rds, "ps_rare.rds"))
meta <- data.frame(sample_data(ps_rare), check.names = FALSE)

# ---- Alpha diversity ----------------------------------------------------------
log_step("Alpha diversity")
alpha <- estimate_richness(ps_rare, measures = c("Observed", "Chao1", "Shannon", "Simpson"))
alpha <- data.frame(sample_id = sample_names(ps_rare), dataset = meta$dataset,
                    alpha[, c("Observed", "Chao1", "se.chao1", "Shannon", "Simpson")],
                    row.names = NULL)
write_supp(alpha, "Table_S2_alpha_diversity")

# Table 2: pooled descriptive statistics (no inferential test - see Methods)
summ <- function(x) c(mean = mean(x), sd = sd(x), min = min(x), max = max(x))
tab2 <- t(sapply(alpha[, c("Observed", "Chao1", "Shannon", "Simpson")], summ))
tab2 <- data.frame(index = rownames(tab2), round(tab2, 2), row.names = NULL)
print(tab2)
write.csv(tab2, file.path(dir_tables, "Table2_alpha_summary.csv"), row.names = FALSE)

alpha_long <- rbind(
  data.frame(sample_id = alpha$sample_id, index = "Observed ASVs",   value = alpha$Observed),
  data.frame(sample_id = alpha$sample_id, index = "Chao1",           value = alpha$Chao1),
  data.frame(sample_id = alpha$sample_id, index = "Shannon",         value = alpha$Shannon))
alpha_long$index <- factor(alpha_long$index, levels = c("Observed ASVs", "Chao1", "Shannon"))

p_alpha <- ggplot(alpha_long, aes(x = "", y = value)) +
  geom_violin(fill = "#9ecae1", colour = "grey30", width = 0.8, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.06, size = 1.4, alpha = 0.7) +
  facet_wrap(~ index, scales = "free_y") +
  labs(x = NULL, y = "Value") +
  theme_bw(base_size = 11) +
  theme(axis.ticks.x = element_blank(), strip.background = element_rect(fill = "grey92"))
save_fig("Fig_alpha_violin", p_alpha, width = 8, height = 4)

# ---- Beta diversity ------------------------------------------------------------
log_step("Beta diversity: Bray-Curtis PCoA")
ord  <- ordinate(ps_rare, method = "PCoA", distance = "bray")
dist <- phyloseq::distance(ps_rare, method = "bray")
eig  <- ord$values$Relative_eig
message(sprintf("  PC1 = %.1f %%, PC2 = %.1f %%, total = %.1f %%",
                100 * eig[1], 100 * eig[2], 100 * sum(eig[1:2])))

set.seed(par$rare_seed)
perm <- adonis2(dist ~ dataset, data = meta, permutations = par$n_perm)
print(perm)
set.seed(par$rare_seed)
bd   <- betadisper(dist, meta$dataset)
bd_t <- permutest(bd, permutations = par$n_perm)
print(bd_t)

stats <- data.frame(
  test = c("PERMANOVA (adonis2)", "betadisper (permutest)"),
  R2   = c(perm$R2[1], NA),
  F    = c(perm$F[1], bd_t$tab$F[1]),
  p    = c(perm$`Pr(>F)`[1], bd_t$tab$`Pr(>F)`[1]),
  permutations = par$n_perm)
write.csv(stats, file.path(dir_tables, "beta_diversity_tests.csv"), row.names = FALSE)

pc <- data.frame(ord$vectors[, 1:2], sample_id = rownames(ord$vectors))
pc$dataset <- meta[pc$sample_id, "dataset"]
p_pcoa <- ggplot(pc, aes(Axis.1, Axis.2, colour = dataset, shape = dataset)) +
  geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
  stat_ellipse(type = "t", linetype = 2, show.legend = FALSE) +
  geom_point(size = 3) +
  scale_colour_manual(values = c("#0072B2", "#D55E00")) +
  labs(x = sprintf("PC1 (%.1f %%)", 100 * eig[1]), y = sprintf("PC2 (%.1f %%)", 100 * eig[2]),
       colour = "Dataset", shape = "Dataset",
       caption = sprintf("PERMANOVA R2 = %.3f, F = %.3f, p = %.3f; betadisper p = %.3f",
                         perm$R2[1], perm$F[1], perm$`Pr(>F)`[1], bd_t$tab$`Pr(>F)`[1])) +
  theme_bw(base_size = 11)
save_fig("Fig_PCoA_bray", p_pcoa, width = 6.5, height = 5)

saveRDS(list(alpha = alpha, ord = ord, perm = perm, betadisper = bd_t),
        file.path(dir_rds, "diversity_results.rds"))
message("\nDone: 04_alpha_beta_diversity.R")
