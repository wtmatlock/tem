############################################################
# FIGURE 2: TWO SEPARATE PDFS (SIMPLE & CLEAN)
############################################################

library(ape)
library(phytools)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtree)
library(bayestestR)
library(RColorBrewer)

setwd("~/desktop/final_project")

cat("\n===== FIGURE 2 AS TWO SEPARATE PDFS =====\n\n")

# Load model and data
chains <- readRDS("output/mic_chains_final.rds")
model <- chains$chain.1
df.model <- readRDS("output/04_df_model_FINAL.rds")
phylo <- read.tree("data/ml_tree.txt")
df <- read_delim("data/supplementary_data_1.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)

# Prepare tree
phylo$tip.label <- as.numeric(phylo$tip.label)
phylo$node.label <- NULL
phylo <- midpoint.root(phylo)
phylo.u <- chronos(phylo, lambda=1, model="correlated")
phylo.plot <- keep.tip(phylo.u, as.character(df.model$phylo))
phylo.plot$tip.label <- as.character(phylo.plot$tip.label)

# Prepare metadata
isolates <- unique(df.model$isolate.assembly)
df.plot.phylo <- df %>%
  filter(contig.type=="chromosome") %>%
  filter(isolate.assembly %in% isolates) %>%
  select(isolate.assembly, isolate.id, ezclermont.phylogroup, mlst.st)

colnames(df.plot.phylo) <- c("tip.label", "id", "phylogroup", "mlst.st")
df.plot.phylo$tip.label <- as.character(df.plot.phylo$tip.label)

# Extract phylogenetic effects
phylo.effects <- model$Sol[, grep("phylo", colnames(model$Sol))]
phylo.effects <- as.data.frame(phylo.effects)

phylo.effects.long <- pivot_longer(
  phylo.effects, 
  cols = everything(), 
  names_to = "tip", 
  values_to = "value"
)

phylo.effects.long$tip <- gsub("phylo.", "", phylo.effects.long$tip)

# Calculate summary statistics
phylo.effects.summary <- phylo.effects.long %>%
  group_by(tip) %>%
  summarise(
    mean = mean(value),
    median = median(value),
    sd = sd(value),
    .groups = "drop"
  )

hdi_results <- phylo.effects.long %>%
  group_by(tip) %>%
  summarise(
    hdi = list(ci(value, method = "HDI")),
    .groups = "drop"
  ) %>%
  unnest_wider(hdi)

phylo.effects.summary <- phylo.effects.summary %>%
  left_join(hdi_results %>% select(tip, CI_low, CI_high), by = "tip")

cat("Data loaded and effects calculated.\n\n")

# ===================================================================
# PDF 1: PHYLOGENETIC TREE
# ===================================================================

cat("Creating PDF 1: Phylogenetic tree...\n")

ggsave("output/Figure2a_tree.pdf", 
       {
         p <- ggtree(phylo.plot) +
           geom_treescale(x=0, y=1, width=0.01) +
           theme_tree2()
         
         p <- p %<+% df.plot.phylo + 
           geom_tippoint(aes(fill = phylogroup, colour=NULL), 
                         shape = 21, size = 1.5, stroke = 0.3) +
           scale_fill_brewer(palette = "Set3", name="Phylogroup") +
           theme(legend.position = c(0.05, 0.95),
                 legend.background = element_rect(fill="white", color="black", size=0.5))
         p
       },
       width = 14, height = 20, dpi = 300)

cat("✓ Saved: Figure2a_tree.pdf (14 x 20 inches)\n\n")

# ===================================================================
# PDF 2: FOREST PLOT
# ===================================================================

cat("Creating PDF 2: Forest plot of phylogenetic effects...\n")

# Sort by mean effect
phylo.effects.summary <- phylo.effects.summary %>%
  arrange(mean) %>%
  mutate(tip = factor(tip, levels = tip))

p.forest <- ggplot(phylo.effects.summary, aes(x = mean, y = tip)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                 height = 0.2, color = "gray30", size = 0.4) +
  geom_point(aes(fill = mean), shape = 21, colour = "black", 
             size = 2, stroke = 0.3) +
  scale_fill_gradient2(
    midpoint = 0, 
    high = "#c51b7d",
    mid = "#f7f7f7",
    low = "#4d9221",
    name = "Posterior mean"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "Posterior mean (95% HDI)", y = NULL, 
       title = "Phylogenetic effects on co-amoxiclav MIC") +
  xlim(-2.5, 5)

ggsave("output/Figure2b_forest.pdf", p.forest, width = 10, height = 16, dpi = 300)

cat("✓ Saved: Figure2b_forest.pdf (10 x 16 inches)\n\n")

# ===================================================================
# SUMMARY
# ===================================================================

cat("===== FIGURE 2 COMPLETE =====\n\n")
cat("Figure 2a: Phylogenetic tree colored by phylogroup (14 x 20 inches)\n")
cat("Figure 2b: Forest plot of phylogenetic effects (10 x 16 inches)\n\n")

cat("Summary:\n")
print(phylo.effects.summary %>% head(10))

write.csv(phylo.effects.summary, "output/phylo_effects_summary.csv", row.names=FALSE)
cat("\n✓ Saved summary statistics to: phylo_effects_summary.csv\n\n")