library(ape)
library(readr)
library(MCMCglmm)
library(phytools)
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggtext)
library(tidyr)
library(parameters)
library(coda)
library(stringr)

### read in data

# Main metadata file
df <- read_delim("path/to/supplementary_table_1.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
# blaTEM-1 and linked promoter annotations
tem1.report <- read_delim("path/to/supplementary_table_2.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
phylo <- read.tree("path/to/GTR_F_I_R4.treefile")

### prepare phylogeny

phylo$tip.label <- as.numeric(phylo$tip.label)
phylo$node.label <- NULL

phylo <- midpoint.root(phylo)

phylo.u <- chronos(phylo, lambda=1, model="correlated")

is.rooted(phylo.u)
is.ultrametric(phylo.u)

inv.phylo <- inverseA(phylo.u,nodes="TIPS",scale=TRUE)

### prepare metadata

df.model <- df %>% 
  mutate(tem1.contig.copy.number = contig.copy.number * tem1.replicon) %>%
  group_by(isolate.id) %>%
  mutate(tem1.isolate.copy.number = sum(tem1.contig.copy.number)) %>%
  ungroup() %>%
  mutate(tem1.isolate.copy.number.scaled = (tem1.isolate.copy.number - mean(tem1.isolate.copy.number)) / sd(tem1.isolate.copy.number),
         tem1.isolate.scaled = ifelse(tem1.isolate>1, TRUE, FALSE)) %>%
  mutate(tem1.isolate.scaled = as.factor(tem1.isolate.scaled)) %>%
  filter(contig.type=="chromosome") %>%
  ungroup() %>%
  select(isolate.id, isolate.assembly, coamox.mic, tem1.isolate.scaled, tem1.isolate.copy.number.scaled, ampc.promoter.snv) %>%
  mutate(coamox.mic = factor(coamox.mic, ordered = TRUE, levels = c("<=2.2", "4.2", "8.2",
                                                                    "16.2", "32.2", ">32.2"))) 

upper_percentile_value <- quantile(df.model$tem1.isolate.copy.number.scaled, 0.95)
df.model$tem1.isolate.copy.number.scaled[df.model$tem1.isolate.copy.number.scaled > upper_percentile_value] <- upper_percentile_value

tem1.model <- tem1.report %>%
  group_by(isolate.assembly) %>%
  filter(!is.na(promoter.start)) %>%
  mutate(n.linked = n()) %>%
  filter(n.linked == tem1.isolate) %>%
  group_by(isolate.assembly) %>%
  mutate(promoter.snv.types = n_distinct(promoter.snv)) %>%
  filter(promoter.snv.types==1) %>%
  group_by(promoter.snv) %>%
  mutate(promoter.snv.n = n()) %>%
  filter(promoter.snv.n >= 10) %>%
  mutate(promoter.snv = as.factor(toupper(promoter.snv))) %>%
  select(isolate.assembly, promoter.snv) %>%
  distinct()

df.model <- merge(df.model, tem1.model, by="isolate.assembly", all.x=TRUE) %>%
  filter(!is.na(promoter.snv)) 

df.model <- df.model %>%
  group_by(ampc.promoter.snv) %>%
  filter(n()>=10)

summary(as.factor(df.model$promoter.snv))
df.model$promoter.snv = factor(df.model$promoter.snv, ordered = FALSE, 
                               levels = c("CGGCGG", "CGGCGA", "TGGCGA",
                                          "TGGCGG"))

summary(as.factor(df.model$ampc.promoter.snv))
df.model$ampc.promoter.snv = factor(df.model$ampc.promoter.snv, ordered = FALSE, 
                               levels = c("GGCTCCTAGGG", "AGCTTCTAGGG", "AGCTCCTAGGG",
                                          "GATTCCTAGGG"))

df.model <- as.data.frame(df.model)
df.model$phylo <- as.factor(df.model$isolate.assembly)

# run model

prior <- list(G=list(G1=list(V=1,nu=0.02, alpha.mu = 0, alpha.V = 1e+3)),
                R = list(fix=1, V=1, nu=0.02))

set.seed(1)
chain.1 <- MCMCglmm(coamox.mic ~ tem1.isolate.copy.number.scaled + tem1.isolate.scaled + 
                      ampc.promoter.snv + promoter.snv,
                    random=~phylo,
                    family="ordinal",
                    ginverse=list(phylo=inv.phylo$Ainv),
                    prior=prior,
                    data=df.model,
                    nitt=10000000,
                    burnin=1000000,
                    thin=100,
                    DIC=FALSE,
                    trunc=TRUE,
                    pr=TRUE)

set.seed(2)
chain.2 <- MCMCglmm(coamox.mic ~ tem1.isolate.copy.number.scaled + tem1.isolate.scaled + 
                      ampc.promoter.snv + promoter.snv,
                    random=~phylo,
                    family="ordinal",
                    ginverse=list(phylo=inv.phylo$Ainv),
                    prior=prior,
                    data=df.model,
                    nitt=10000000,
                    burnin=1000000,
                    thin=100,
                    DIC=FALSE,
                    trunc=TRUE,
                    pr=TRUE)

summary(chain.1)
autocorr.diag(chain.1$VCV) 
plot(chain.1)

summary(chain.2)
autocorr.diag(chain.2$VCV) 
plot(chain.2)

mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
gelman.diag(mclist)

model <- chain.1

model_parameters(model, centrality="mean", ci=0.95, ci_method="hdi",
                 component = "all")

# phylo

library(ggtree)
library(ggnewscale)
library(bayestestR)
library(RColorBrewer)

phylo.plot <- keep.tip(phylo, as.character(df.model$phylo))
phylo.plot$tip.label <- as.character(phylo.plot$tip.label)

isolates <- unique(df.model$isolate.assembly)
df.plot.phylo <- df %>%
  filter(contig.type=="chromosome") %>%
  filter(isolate.assembly %in% isolates) %>%
  select(isolate.assembly, isolate.id, ezclermont.phylogroup, mlst.st)
colnames(df.plot.phylo) <- c("tip.label", "id", "phylogroup", "mlst.st")
df.plot.phylo$tip.label <- as.character(df.plot.phylo$tip.label)

p <- ggtree(phylo.plot) +
  geom_treescale(x=0, y=180, width=0.01)

p <- p %<+% df.plot.phylo + 
  geom_tiplab(align=TRUE,aes(label=id), offset=0.002, size=1) +
  geom_tippoint(aes(fill = phylogroup, colour=NULL), shape = 21, size = 1, stroke = 0.5) +
  scale_fill_brewer(palette = "Set3", name="Phylogroup")

p <- p + theme(legend.position = c(0.1,0.85)) + ggplot2::xlim(0, 0.03)


p

get_taxa_name(p)

phylo.effects <- model$Sol[, grep("phylo", colnames(model$Sol))]
phylo.effects <- as.data.frame(phylo.effects)
phylo.effects.long <- pivot_longer(phylo.effects, cols = everything(), names_to = "tip", values_to = "value")
phylo.effects.long$tip <- gsub("phylo.", "", phylo.effects.long$tip)
phylo.effects.long$tip <- factor(phylo.effects.long$tip, levels = get_taxa_name(p))

write.csv(phylo.effects.long, 'mic-full-effects_rev.csv')
# used in modelExpression.R
# "mic.phylo.effects <- read_csv("path/to/mic-effects_rev.csv")"
# for plotting Figure 2

phylo.effects.long <- phylo.effects.long %>%
  group_by(tip) %>%
  mutate(value = value) %>%
  mutate(hpd = ci(value, method = "HDI")) %>%
  mutate(mean = mean(value)) %>%
  select(tip, mean, hpd) %>%
  distinct()
phylo.effects.long <- cbind(phylo.effects.long, phylo.effects.long$hpd)

write.csv(phylo.effects.long, 'mic-effects_rev.csv')

p.2 <- ggplot(phylo.effects.long, aes(x = tip, y = mean)) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0) +
  geom_point(aes(fill = mean), shape = 21, colour = "black", size = 1, stroke = 0.5) +
  scale_fill_gradient2(midpoint = 0, high = "#c51b7d", mid = "#f7f7f7", low = "#4d9221") +
  theme_minimal() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),    
        axis.ticks.y = element_blank()) +
  scale_x_discrete(limits=rev) +
  labs(fill = "Posterior mean") +
  ylab(label="Posterior 95% HDI")

p.2 <- p.2 + theme(legend.position = c(-1.25,0.89))

library(cowplot)
plot_grid(p, NULL, p.2, align='hv', 
          rel_widths = c(1, 0, 0.6), nrow=1,
          labels=c("a", "b"))

###
### additional analysis for tip effect in terms of MIC category ###
###   

# latent cutpoint for 4/2 is 1.222
# latent cutpoint for 

tips <- as.character(df.model$phylo)

cutpoint_samples <- as.data.frame(model$CP)
phylo_samples <- as.data.frame(model$Sol[,-(1:9)])
colnames(phylo_samples) <- gsub("phylo.", "", colnames(phylo_samples))
phylo_samples <- phylo_samples %>% select(any_of(tips))
data <- cbind(cutpoint_samples, phylo_samples)

cp2 <- "cutpoint.traitcoamox.mic.2" 
cp3 <- "cutpoint.traitcoamox.mic.3" 

midpoint <- (data[[cp3]]+data[[cp2]])/2

p_up <- sapply(tips, function(col) mean(data[[col]] + midpoint > data[[cp3]]))

results <- data.frame(
  tip     = tips,
  p_up    = p_up,
  stringsAsFactors = FALSE
)

sum(results$p_up > 0.5)

plot(density(results$p_up))

df_to_merge <- df %>%
  filter(contig.type == "chromosome") %>%
  select(isolate.assembly, ezclermont.phylogroup, mlst.st)

results <- merge(results, df_to_merge, by.x="tip", by.y="isolate.assembly", all.x=TRUE)

isolates <- results %>%
  filter(p_up > 0.5) %>%
  pull(tip)

phylogroup_up <- df %>%
  filter(contig.type == "chromosome") %>%
  filter(isolate.assembly %in% isolates) %>%
  group_by(ezclermont.phylogroup) %>%
  summarise(n=n())

phylogroup_up

mlst_up <- df %>%
  filter(contig.type == "chromosome") %>%
  filter(isolate.assembly %in% isolates) %>%
  group_by(mlst.st) %>%
  summarise(n=n())

mlst_up

plot_df <- results %>%
  mutate(perc_50 = ifelse(p_up > 0.5, TRUE, FALSE),
         perc_60 = ifelse(p_up > 0.6, TRUE, FALSE),
         perc_70 = ifelse(p_up > 0.7, TRUE, FALSE))

to_merge <- df %>%
  filter(contig.type=="chromosome") %>%
  filter(isolate.assembly %in% df.model$isolate.assembly) %>%
  select(isolate.assembly, mlst.st)

plot_df <- merge(plot_df, to_merge, by.x='tip', by.y='isolate.assembly', all.x=TRUE)

plot_df_long <- plot_df %>% 
  pivot_longer(
    cols      = starts_with("perc_"),
    names_to  = "perc",
    values_to = "flag"
  ) %>% 
  filter(flag)

plot_mlsts <- plot_df_long %>% 
  filter(mlst.st != "-") %>% 
  count(mlst.st) %>% 
  filter(n >= 10) %>% 
  pull(mlst.st)

plot_df_long <- plot_df_long %>% 
  mutate(
    mlst_group = if_else(mlst.st %in% plot_mlsts,
                         as.character(mlst.st),
                         "Other")
  )

df_count <- plot_df_long %>% 
  count(perc, mlst_group)

perc_labels <- c(
  perc_50 = "50",
  perc_60 = "60",
  perc_70 = "70"
)

my_cols <- c(
  "Other" = "grey70",
  "127"   = "#FE6100",
  "69"    = "#FFB000",
  "12"  = "#DC267F"
)

ggplot(df_count, aes(x = perc, y = n, fill = mlst_group)) +
  geom_col() +
  scale_x_discrete(
    labels = perc_labels
  ) +
  scale_fill_manual(
    values = my_cols,
    na.value = "black",      # in case there are any unexpected levels
    guide    = guide_legend(title = "MLST")
  ) +
  labs(
    x = "Percent threshold",
    y = "Number of isolates"
  ) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=20, color='black'),
        panel.border = element_rect(color = "black", fill = NA, size = 2),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white"), 
        axis.ticks = element_line(size = 1, color='black'),
        axis.text.x = element_text(color='black'),
        axis.text.y = element_text(color='black'),
        axis.title.x = element_markdown(),
        axis.title.y = element_markdown())
