library(readr)
library(dplyr)
library(tidyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(ape)
library(ggtree)
library(tibble)
library(gridExtra)
library(phytools)
library(RColorBrewer)

### wrangle data and create main dataframe ###

metadata <- read_delim("path/to/supplementary_table_1.tsv", 
                                    delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)

igr_pa <- read_csv("path/to/piggy_out/IGR_presence_absence.csv")

gene_pa <- read_csv("path/to/panaroo_output/gene_presence_absence_roary.csv")

igr_df <- igr_pa %>%
  as.data.frame() %>%                     
  select(-c(2,3,5:14)) %>%
  pivot_longer(
    cols = -c(Gene, `No. isolates`), 
    names_to = "isolate",
    values_to = "value"
  ) %>%
  drop_na() %>%
  rename("cluster_id" = Gene,
         "cluster_n" = `No. isolates`,
         "orientation" = value) %>%
  mutate(
    orientation = str_extract(orientation, "(DT|DP|CO_F|CO_R)"),
    type = "IGR"
  ) %>%
  select(isolate, cluster_id, cluster_n, type, orientation)

gene_df <- gene_pa %>%
  as.data.frame() %>%                     
  select(-c(2:14)) %>%
  pivot_longer(
    cols = -c(Gene), 
    names_to = "isolate",
    values_to = "value"
  ) %>%
  mutate(value = if_else(str_detect(value, ";"), NA, value)) %>% # filter out fragmented genes
  drop_na() %>%
  group_by(Gene) %>%
  mutate(cluster_n = n(),
         orientation = NA,
         type = "gene") %>%
  ungroup() %>%
  rename("cluster_id" = Gene) %>%
  select(isolate, cluster_id, cluster_n, type, orientation)

combined_df <- rbind(igr_df, gene_df)

isolate_df <- metadata %>%
  filter(contig.type=="chromosome") %>%
  select(isolate.id, isolate.assembly, ezclermont.phylogroup, mlst.st)

combined_df <- merge(combined_df, isolate_df, by.x="isolate", by.y="isolate.assembly", all.x=TRUE)

### orientation analysis ###

orientation_summary <- combined_df %>%
  filter(type=="IGR" & cluster_n==377) %>%
  group_by(cluster_id) %>%
  filter(n_distinct(orientation) > 1) %>%
  group_by(cluster_id, orientation) %>%
  summarise(n=n())

ggplot() +
  geom_bar()

### cluster frequencies ###

unique_gene_cluster <- combined_df %>%
  filter(type=="gene") %>%
  select(cluster_id) %>%
  unique() %>%
  nrow()

gene_cluster_frequencies <- combined_df %>%
  filter(type=="gene") %>%
  select(cluster_id, cluster_n) %>%
  unique() %>%
  group_by(cluster_n) %>%
  summarise(n = n()) %>%
  mutate(perc = (n/unique_gene_cluster)*100)

unique_igr_cluster <- combined_df %>%
  filter(type=="IGR") %>%
  select(cluster_id) %>%
  unique() %>%
  nrow()

igr_cluster_frequencies <- combined_df %>%
  filter(type=="IGR") %>%
  select(cluster_id, cluster_n) %>%
  unique() %>%
  group_by(cluster_n) %>%
  summarise(n = n()) %>%
  mutate(perc = (n/unique_igr_cluster)*100)

### IGR cluster conservation ###

igr_divergences <- read_csv("path/to/piggy_out/cluster_IGR_divergences.csv")
  
divergence_df <- igr_divergences %>%
  mutate(cluster_id = gsub("_aligned", "", Gene)) %>%
  select(cluster_id, Nuc_identity, Length_identity) 

divergence_df <- merge(combined_df, divergence_df, by="cluster_id", all.y=TRUE) %>%
  select(cluster_id, cluster_n, Nuc_identity, Length_identity) %>%
  unique()

mod_div <- lm(Nuc_identity ~ cluster_n, data = divergence_df)
summary(mod_div)
  
### accumulation curve analysis ###

# function to compute accumulation curves, filtering clusters by min_freq

evaluate_accumulation_curves <- function(df, min_freq = 1, n_reps = 100, phylo_min = 20) {
  
  # identify phylogroups with at least phylo_min isolates of this feature type
  good_phylo <- df %>%
    distinct(isolate.id, ezclermont.phylogroup) %>%
    count(ezclermont.phylogroup, name = "n_iso") %>%
    filter(n_iso >= phylo_min) %>%
    pull(ezclermont.phylogroup)
  
  # loop over phylogroups
  map_dfr(good_phylo, function(pg) {
    sub_df <- df %>% filter(ezclermont.phylogroup == pg)
    isos <- sub_df %>% distinct(isolate.id) %>% pull(isolate.id)
    
    # compute cluster frequency for this feature globally
    freq_tbl <- df %>%
      distinct(isolate.id, cluster_id) %>%
      count(cluster_id, name = "freq")
    valid <- freq_tbl %>% filter(freq >= min_freq) %>% pull(cluster_id)
    
    # presence/absence of valid clusters in this phylogroup
    pa <- sub_df %>%
      filter(cluster_id %in% valid) %>%
      distinct(isolate.id, cluster_id)
    
    # random sampling replicates
    map_dfr(seq_len(n_reps), function(rep) {
      ord <- sample(isos)
      seen <- character(0)
      tibble(
        replicate    = rep,
        step         = seq_along(ord),
        phylogroup   = pg,
        cum_clusters = map_int(step, function(i) {
          seen <<- union(seen, pa$cluster_id[pa$isolate.id == ord[i]])
          length(seen)
        })
      )
    }) %>% mutate(min_freq = min_freq)
  })
}

# split input by feature type
genes_df <- combined_df %>% filter(type == "gene")
igr_df   <- combined_df %>% filter(type == "IGR")

# generate curves for all clusters and clusters in >=2 isolates
gc_all   <- evaluate_accumulation_curves(genes_df, min_freq = 1)
gc_2plus <- evaluate_accumulation_curves(genes_df, min_freq = 2)
ig_all   <- evaluate_accumulation_curves(igr_df,   min_freq = 1)
ig_2plus <- evaluate_accumulation_curves(igr_df,   min_freq = 2)

# annotate feature and threshold for plotting
all_curves <- bind_rows(
  gc_all   %>% mutate(feature_type = "gene", threshold = "1"),
  gc_2plus   %>% mutate(feature_type = "gene", threshold = "2"),
  ig_all   %>% mutate(feature_type = "IGR",  threshold = "1"),
  ig_2plus %>% mutate(feature_type = "IGR",  threshold = "2")
)

mean_curves <- all_curves %>%
  group_by(phylogroup, threshold, feature_type, step) %>%
  summarise(mean_cum = mean(cum_clusters), .groups = "drop")

# compute crossover points: first step where mean IGR >= mean gene
gcross <- mean_curves %>%
  pivot_wider(names_from = feature_type, values_from = mean_cum) %>%
  group_by(phylogroup, threshold) %>%
  arrange(step) %>%
  summarise(
    crossover_step = if (any(IGR >= gene, na.rm = TRUE)) {
      min(step[which(IGR >= gene)])
    } else {
      NA_integer_
    },
    .groups = "drop"
  )

# fit power law models
exponent_df <- mean_curves %>%
  filter(step >= 2) %>%
  group_by(phylogroup, threshold, feature_type) %>%
  summarise(
    alpha = coef(lm(log(mean_cum) ~ log(step)))[2],
    kappa = exp(coef(lm(log(mean_cum) ~ log(step)))[1]),
    .groups = "drop"
  )

phylogroup_n <- combined_df %>%
  select(isolate.id, ezclermont.phylogroup) %>%
  unique() %>%
  group_by(ezclermont.phylogroup) %>%
  summarise(n=n())

exponent_df <- merge(exponent_df, phylogroup_n, all.x=TRUE, by.x="phylogroup", by.y="ezclermont.phylogroup")

mod_alpha <- lm(alpha ~ n + feature_type, data = exponent_df)
summary(mod_alpha)

#write.csv(exponent_df, '~/Desktop/heaps.csv')

# Plot
ggplot(all_curves, aes(x = step, y = cum_clusters, color = feature_type)) +
  geom_line(aes(group = interaction(replicate, feature_type, threshold, phylogroup)), alpha = 0.1) +
  stat_summary(aes(group = interaction(feature_type, threshold, phylogroup)), fun = mean, geom = "line", size = 1.1) +
  geom_vline(data = gcross, aes(xintercept = crossover_step), linetype = "dashed", color = "black") +
  scale_y_continuous(breaks = seq(0, max(all_curves$cum_clusters, na.rm = TRUE), by = 3000)) +
  facet_grid(phylogroup ~ threshold, scales = "free") +
  scale_color_manual(
    values = c("gene" = "#E66100", "IGR" = "#5D3A9B"),
    labels = c("gene" = "Gene", "IGR" = "IGR")
  ) +
  labs(
    x = "Number of isolates sampled",
    y = "Cumulative clusters",
    color = "Cluster\ntype"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 20, color = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 20, color = 'black'),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    strip.background = element_rect(fill = "black"),
    axis.ticks = element_line(linewidth = 1, color = 'black'),
    axis.text.x = element_text(color = 'black'),
    axis.text.y = element_text(color = 'black'),
    legend.position = "right"
  )


### plot phylogeny ###

igr_heatmap <- igr_pa %>%
  as.data.frame() %>%                     
  select(-c(1:14)) %>%
  t() %>%
  as.data.frame() %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, 1))) %>%
  select(order(-colSums(.))) %>%
  mutate(isolate = rownames(.))

ordered_clusters <- head(colnames(igr_heatmap), -1)

phylo <- read.tree("path/to/GTR_F_I_R4.treefile") 
phylo <- midpoint.root(phylo)
ordered_isolates <- phylo$tip.label 

igr_heatmap_mat <- igr_heatmap %>%
  arrange(match(isolate, ordered_isolates)) %>%
  select(-isolate) %>%
  select(1:10000) %>% # most common 10000 clusters
  as.matrix()

igr_heatmap_long <- as.data.frame(igr_heatmap_mat) %>%
  rownames_to_column("isolate") %>%
  pivot_longer(
    cols      = -isolate,
    names_to  = "variable",
    values_to = "value"
  ) %>%
  mutate(isolate = factor(isolate, levels = ordered_isolates),
         variable = factor(variable, levels = ordered_clusters))

tip_df <- combined_df %>%
  distinct(isolate, ezclermont.phylogroup) %>%
  unique()

color_palette <- brewer.pal(n = 9, name = "Set3")
tip_df$ezclermont.phylogroup <- factor(tip_df$ezclermont.phylogroup, 
                                       levels = sort(unique(tip_df$ezclermont.phylogroup)))
color_mapping <- setNames(color_palette, levels(tip_df$ezclermont.phylogroup))

tree_plot <- ggtree(phylo, layout = "rectangular") %<+% tip_df +
  geom_tiplab(size = 0, align = TRUE, linesize = 0.3) +
  geom_tippoint(aes(color = ezclermont.phylogroup), 
                size = 1,
                show.legend = FALSE) +
  scale_color_manual(values = color_mapping) +
  theme_tree2()

heatmap_plot <- ggplot(igr_heatmap_long, aes(x = variable, y = isolate, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "black") +
  scale_y_discrete(limits = rev(ordered_isolates)) + 
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  
    axis.text.y = element_blank(),  
    axis.ticks = element_blank(),  
    axis.title = element_blank(),   
    legend.position = "none",       
    panel.grid = element_blank(),   
    panel.border = element_blank(), 
    plot.margin = margin(l=10, r=10, t=10, b=10) 
  )


pdf("Figure4.pdf", width = 10, height = 8)  # Adjust width and height as needed
grid.arrange(tree_plot, heatmap_plot, ncol = 2, widths = c(1, 3))  # Adjust widths as needed
dev.off() 
