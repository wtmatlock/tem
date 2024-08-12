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
library(readxl)

### read in data

df <- read_delim("~/Desktop/ecoli548/metadata/metadata.tsv", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)
tem1.report <- read_delim("~/Desktop/ecoli548/metadata/tem1_report.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)
phylo <- read.tree("/Users/willmatlock/Desktop/GTR_F_I_R4.treefile")
exp <- read_excel("~/Desktop/ecoli548/delta_CtTEM_Ct16S.xlsx")


### prepare phylogeny

phylo$tip.label <- as.numeric(phylo$tip.label)
phylo$node.label <- NULL

phylo <- midpoint.root(phylo)

phylo.u <- chronos(phylo, lambda=1, model="correlated")

is.rooted(phylo.u)
is.ultrametric(phylo.u)

inv.phylo <- inverseA(phylo.u,nodes="TIPS",scale=TRUE)

# prepare data

df.model <- df %>%
  filter(tem1.isolate==1 & tem1.replicon==1) %>%
  select(accession, isolate.assembly, contig.assembly, contig.copy.number, contig.type)

tem1.report.less <- tem1.report[c("isolate.assembly", "contig.assembly", "promoter.snv", "tem1.snv")]

df.model <- merge(df.model, tem1.report.less, 
                    by=c("isolate.assembly", "contig.assembly"),
                    all.x=TRUE)

df.model$promoter.type <- ifelse(df.model$promoter.snv=="tggcgg", "TG",
                                   ifelse(df.model$promoter.snv=="tggcga", "TA",
                                   ifelse(df.model$promoter.snv=="cggcgg", "CG",
                                   ifelse(df.model$promoter.snv=="cggcga", "CA",
                                   NA))))

df.model <- df.model %>% filter(!is.na(promoter.type))

exp$Isolate <- sapply(strsplit(as.character(exp$Isolate), " "), function(x) x[1])
exp <- exp[c("Isolate", "delta Ct (control)")]
colnames(exp) <- c("accession", "exp")

df.model <- merge(df.model, exp, by="accession", all.y=TRUE)
df.model <- df.model %>% filter(!is.na(exp)) %>% filter(!is.na(isolate.assembly)) %>% distinct()
df.model <- df.model %>% group_by(accession) %>% mutate(run = row_number())

df.model$contig.copy.number <- as.numeric(df.model$contig.copy.number )

mean.exp <- mean(df.model$exp)
sd.exp <- sd(df.model$exp)
mean.cp <- mean(df.model$contig.copy.number)
sd.cp <- sd(df.model$contig.copy.number)

df.model <- df.model %>%
  mutate(exp.scaled = (exp - mean.exp) / sd.exp,
         contig.copy.number.scaled = (contig.copy.number - mean.cp)/sd.cp)

up.exp <- quantile(df.model$exp.scaled, 0.95)
df.model$exp.scaled[df.model$exp.scaled > up.exp] <- up.exp

new.cols <- data.frame(t(data.frame(strsplit(df.model$promoter.type, ""))))
df.model$pos1 <- new.cols$X1
df.model$pos2 <- new.cols$X2
df.model$pos1.bool <- ifelse(df.model$pos1=="T", TRUE, FALSE)
df.model$pos2.bool <- ifelse(df.model$pos2=="A", TRUE, FALSE)

df.model$phylo <- as.character(df.model$isolate.assembly)
length(unique(df.model$phylo))
phylo.trim <- keep.tip(phylo.u, df.model$phylo)
inv.phylo.trim<-inverseA(phylo.trim,nodes="TIPS",scale=TRUE)

df.model$phylo <- as.character(df.model$phylo)

df.model$isolate.assembly <- as.factor(df.model$isolate.assembly)
df.model$tem1.snv <- as.factor(df.model$tem1.snv)

df.model <- as.data.frame(df.model)

sample <- df.model %>% 
  select(isolate.assembly, promoter.type, contig.type) %>% 
  distinct()
table(sample$promoter.type, sample$contig.type)

pg <- df %>%
  filter(contig.type=="chromosome") %>%
  filter(isolate.assembly %in% sample$isolate.assembly)
summary(as.factor(pg$ezclermont.phylogroup))

replicates <- df.model %>% group_by(isolate.assembly) %>%
  summarise(n())

summary(as.factor(replicates$`n()`))

# run chains

prior.1 <- list(G=list(G1=list(V=1,nu=0.02, alpha.mu = 0, alpha.V = 1e+3),
                       G2 = list(V=1, nu=0.02)),
                R = list(V=1, nu=0.02))

set.seed(1)
chain.1 <- MCMCglmm(exp.scaled ~ pos1.bool*pos2.bool + contig.copy.number.scaled,
                    random= ~ phylo + isolate.assembly,
                    family="gaussian",
                    ginverse=list(phylo=inv.phylo.trim$Ainv),
                    data=df.model,
                    prior=prior.1,
                    nitt=10000000,
                    burnin=1000000,
                    thin=100,
                    DIC=FALSE,
                    pr=TRUE)

#save.image("~/Desktop/exp.RData")

set.seed(2)
chain.2 <- MCMCglmm(exp.scaled ~ pos1.bool*pos2.bool + contig.copy.number.scaled,
                    random= ~ phylo + isolate.assembly,
                    family="gaussian",
                    ginverse=list(phylo=inv.phylo.trim$Ainv),
                    data=df.model,
                    prior=prior.1,
                    nitt=10000000,
                    burnin=1000000,
                    thin=100,
                    DIC=FALSE,
                    pr=TRUE)

#save.image("~/Desktop/exp.RData")

summary(chain.1)
plot(chain.1)

summary(chain.2)
plot(chain.2)

mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
gelman.diag(mclist)

model <- chain.1

# make fig 2

library(ggtree)
library(ggnewscale)
library(bayestestR)
library(RColorBrewer)

phylo.plot <- keep.tip(phylo, df.model$phylo)
phylo.plot$tip.label <- as.character(phylo.plot$tip.label)

isolates <- unique(df.model$isolate.assembly)
df.plot.phylo <- df %>%
  filter(isolate.assembly %in% isolates) %>%
  group_by(isolate.assembly) %>%
  fill(ezclermont.phylogroup, .direction = "downup") %>%
  fill(mlst.st, .direction = "downup") %>%
  filter(tem1.replicon==1) %>%
  select(isolate.assembly, isolate.id, ezclermont.phylogroup, mlst.st, contig.type)
colnames(df.plot.phylo) <- c("tip.label", "id", "phylogroup", "mlst.st", "contig.type")
df.plot.phylo$tip.label <- as.character(df.plot.phylo$tip.label)

df.plot.phylo$ST <- ifelse(df.plot.phylo$mlst.st=="12", "12",
                           ifelse(df.plot.phylo$mlst.st=="372", "372",
                           "Other"))
df.plot.phylo$Replicon <- ifelse(df.plot.phylo$contig.type=="chromosome", "Chromosome",
                           ifelse(df.plot.phylo$contig.type=="plasmid", "Plasmid", NA))

p <- ggtree(phylo.plot) +
  geom_treescale(x=0, y=40, width=0.01)

p <- p %<+% df.plot.phylo + 
  geom_tiplab(align=TRUE,aes(label=id), offset=0.002) +
  geom_tippoint(aes(fill = phylogroup, colour=ST, shape=Replicon), size = 3, stroke = 1) +
  scale_fill_brewer(palette = "Set3", name="Phylogroup") +
  scale_colour_manual(values=c("#1b9e77", "#d95f02", "black")) +
  scale_shape_manual(values = c("Chromosome"=21, "Plasmid"=22)) +
  theme(legend.position = NA)

p <- p + theme(legend.position = c(0.1,0.9)) + ggplot2::xlim(0, 0.03)


p



get_taxa_name(p)

phylo.effects <- model$Sol[, grep("phylo", colnames(model$Sol))]
phylo.effects <- as.data.frame(phylo.effects)
phylo.effects.long <- pivot_longer(phylo.effects, cols = everything(), names_to = "tip", values_to = "value")
phylo.effects.long$tip <- gsub("phylo.", "", phylo.effects.long$tip)
phylo.effects.long$tip <- factor(phylo.effects.long$tip, levels = get_taxa_name(p))

#write.csv(phylo.effects.long, '~/Desktop/exp-full-effects.csv')

phylo.effects.long <- phylo.effects.long %>%
  group_by(tip) %>%
  mutate(value = value) %>%
  mutate(hpd = ci(value, method = "HDI")) %>%
  mutate(mean = mean(value)) %>%
  select(tip, mean, hpd) %>%
  distinct()
phylo.effects.long <- cbind(phylo.effects.long, phylo.effects.long$hpd)

#write.csv(phylo.effects.long, '~/Desktop/exp-effects.csv')

p.2 <- ggplot(phylo.effects.long, aes(x = tip, y = -mean)) +
  geom_errorbar(aes(ymin=-CI_low, ymax=-CI_high), width=0) +
  geom_point(aes(fill = mean), shape = 21, colour = "black", size = 3, stroke = 0.8) +
  scale_fill_gradient2(midpoint = 0, low = "#b2182b", mid = "#f7f7f7", high = "#2166ac") +
   theme_minimal() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),    
        axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits=c(-1, 1)) +
  scale_x_discrete(limits=rev) +
  labs(fill = "Posterior mean") +
  ylab(label="")

library(ggridges)

p.2.alt <- ggplot(phylo.effects.long, aes(x = tip, y = -mean)) +
  geom_density_ridges(aes(fill = mean)) +
  scale_fill_gradient2(midpoint = 0, low = "#b2182b", mid = "#f7f7f7", high = "#2166ac") +
  theme_minimal() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),    
        axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits=c(-1, 1)) +
  scale_x_discrete(limits=rev) +
  labs(fill = "Posterior mean") +
  ylab(label="")

mic.phylo.effects <- read_csv("Desktop/tem1_2024/models/mic-effects.csv")

mic.phylo.effects <- mic.phylo.effects %>%
  filter(tip %in% phylo.effects.long$tip) %>%
  select(tip, mean, CI_low, CI_high)
mic.phylo.effects <- rbind(mic.phylo.effects, c("490", NA, NA, NA))
mic.phylo.effects <- rbind(mic.phylo.effects, c("238", NA, NA, NA))
mic.phylo.effects$tip <- factor(mic.phylo.effects$tip, levels = get_taxa_name(p))
mic.phylo.effects$mean <- as.numeric(mic.phylo.effects$mean)
mic.phylo.effects$CI_low <- as.numeric(mic.phylo.effects$CI_low)
mic.phylo.effects$CI_high <- as.numeric(mic.phylo.effects$CI_high)

p.3 <- ggplot(mic.phylo.effects, aes(x = tip, y = mean)) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0) +
  geom_point(aes(fill = mean), shape = 21, colour = "black", size = 3, stroke = 0.8) +
  scale_fill_gradient2(midpoint = 0, high = "#c51b7d", mid = "#f7f7f7", low = "#4d9221") +
  theme_minimal() +
  coord_flip() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.y = element_blank(),  
        axis.text.y = element_blank(),    
        axis.ticks.y = element_blank()) +
  scale_y_continuous(expand = c(0, 0), limits=c(-4, 4)) +
  scale_x_discrete(limits=rev) +
  labs(fill = "Posterior mean") +
  ylab(label="")

p.2 <- p.2 + theme(legend.position = c(-1,0.9))
p.3 <- p.3 + theme(legend.position = c(-1,0.9))


library(cowplot)
plot_grid(p, NULL, p.2, NULL ,p.3, align='hv', 
          rel_widths = c(1, -0.18, 0.6, -0.15, 0.6), nrow=1,
          labels=c("a", "", "b", "", "c"))