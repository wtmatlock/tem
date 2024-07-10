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
                    DIC=TRUE,
                    pr=TRUE)

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
                    DIC=TRUE,
                    pr=TRUE)

#save.image("~/Desktop/models/exp.RData")

summary(chain.1)
plot(chain.1)

summary(chain.2)
plot(chain.2)

mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
gelman.diag(mclist)

model <- chain.1

# R2 table

mFixed <- mean(model$Sol[,1]) * model$X[, 1] +
  mean(model$Sol[,2]) * model$X[, 2] +
  mean(model$Sol[,3]) * model$X[, 3] +
  mean(model$Sol[,4]) * model$X[, 4] +
  mean(model$Sol[,5]) * model$X[, 5]
  
mVarF<- var(mFixed)

mean_VCV <- apply(model$VCV, 2, mean)

fixed_r2 <- mVarF / (mVarF + sum(mean_VCV))
random_r2 <- mean_VCV[1] / (sum(mean_VCV) + mVarF)
conditional_r2 <- (mean_VCV[1] + mVarF) / (sum(mean_VCV) + mVarF)

R2_table_model <- matrix(c(fixed_r2,random_r2,conditional_r2), ncol=1)
R2_table_model <- format(R2_table_model, digits=4)
rownames(R2_table_model) <- c("Fixed effects", "Random effects","Total")
colnames(R2_table_model) <- c("R-squared value")

R2_table_model

# plot bootstrap 

r <- seq(from=min(df.model$contig.copy.number.scaled), to=max(df.model$contig.copy.number.scaled), length.out=100)
pred <- data.frame()

sol.samples <- model$Sol
vcv.samples <- model$VCV
cp.samples <- model$CP

latent.sd <- sqrt(1 + sum(colMeans(vcv.samples)))

ci <- model_parameters(model, centrality="mean", ci_method="hdi", ci=0.95, bootstrap=TRUE)

ci.intercept <- c(ci$CI_low[1], ci$Mean[1], ci$CI_high[1])
ci.pos1 <- c(ci$CI_low[2], ci$Mean[2], ci$CI_high[2])
ci.pos2 <- c(ci$CI_low[3], ci$Mean[3], ci$CI_high[3])
ci.contig.copy.number.scaled <- c(ci$CI_low[4], ci$Mean[4], ci$CI_high[4])
ci.pos1.pos2 <- c(ci$CI_low[5], ci$Mean[5], ci$CI_high[5])

n_runs <- 1000

# Function to generate predictions accounting for random effects
generate_predictions <- function(pos1.value, pos2.value, n_runs) {
  pred <- data.frame(contig.copy.number.scaled = r)
  pred_list <- vector("list", n_runs)
  
  for (run in 1:n_runs) {
    run_pred <- data.frame()
    for (i in r) {
      linear.predictor <- ci$Mean[1] + ci$Mean[4] * i + pos1.value + pos2.value
      # Incorporate random effects variability
      prediction <- linear.predictor + rnorm(1, 0, latent.sd)
      run_pred <- rbind(run_pred, c(i, prediction))
    }
    colnames(run_pred) <- c("contig.copy.number.scaled", paste("predicted.expression", run, sep = "_"))
    pred_list[[run]] <- run_pred
  }
  
  # Combine all runs into a single data frame
  pred_combined <- Reduce(function(x, y) merge(x, y, by="contig.copy.number.scaled"), pred_list)
  
  # Calculate the mean prediction across all runs
  #pred_combined$predicted.expression <- rowMeans(pred_combined[, -1])
  
  # Select only the relevant columns
  #pred_combined <- pred_combined[, c("contig.copy.number.scaled", "predicted.expression")]
  
  return(pred_combined)
}

pred_CG <- generate_predictions(0, 0, n_runs)
pred_CG$promoter.type <- "CG"

pred_TG <- generate_predictions(ci$Mean[2], 0, n_runs)
pred_TG$promoter.type <- "TG"

pred_CA <- generate_predictions(0, ci$Mean[3], n_runs)
pred_CA$promoter.type <- "CA"

pred <- rbind(pred_CG, pred_TG, pred_CA)

pred.melt <- melt(pred, id=c("contig.copy.number.scaled", "promoter.type"))

pred.melt$contig.copy.number.scaled <- as.numeric(pred.melt$contig.copy.number.scaled)

pred.melt$promoter.type <- factor(pred.melt$promoter.type, levels=c("CG", "CA", "TG"),
                             labels=c("CG (wildtype)", "CA", "TG"))

ggplot(data=pred.melt, aes(x=contig.copy.number.scaled, y=value, colour=promoter.type)) + 
  stat_summary(fun.data=mean_cl_normal, geom="linerange", linewidth=1, alpha=0.5) +
  stat_summary(fun.y=mean, geom="line", linewidth=1) +
  scale_colour_brewer(palette = "Set1") + 
  labs(x = "Contig copy number (normalised)",
       y = "Predicted *bla*<sub>TEM-1</sub> ΔCt (normalised and *P*<sub>95</sub> truncated)",
       colour = "Promoter SNV") +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        text = element_text(size=20, color='black'),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 2),
        strip.background = element_rect(fill="black"),
        strip.text = element_text(color="white"), 
        axis.ticks = element_line(size = 1, color='black'),
        axis.text.x = element_text(color='black'),
        axis.text.y = element_text(color='black'),
        axis.title.y = element_markdown())

# explore phylo effects

library(ggtree)
library(ggnewscale)
library(bayestestR)
library(RColorBrewer)

phylo.plot <- keep.tip(phylo, df.model$phylo)
phylo.plot$tip.label <- as.character(phylo.plot$tip.label)

isolates <- unique(df.model$isolate.assembly)
df.plot.phylo <- df %>%
  filter(contig.type=="chromosome") %>%
  filter(isolate.assembly %in% isolates) %>%
  select(isolate.assembly, isolate.id, ezclermont.phylogroup, mlst.st)
colnames(df.plot.phylo) <- c("tip.label", "id", "phylogroup", "mlst.st")
df.plot.phylo$tip.label <- as.character(df.plot.phylo$tip.label)

df.plot.phylo$ST <- ifelse(df.plot.phylo$mlst.st=="12", "12",
                           ifelse(df.plot.phylo$mlst.st=="372", "372",
                           "Other"))

p <- ggtree(phylo.plot) +
  geom_treescale(x=0, y=40, width=0.01)

p <- p %<+% df.plot.phylo + 
  geom_tiplab(align=TRUE,aes(label=id), offset=0.002) +
  geom_tippoint(aes(fill = phylogroup, colour=ST), shape = 21, size = 3, stroke = 1.2) +
  scale_fill_brewer(palette = "Set3", name="Phylogroup") +
  scale_colour_manual(values=c("#1b9e77", "#d95f02", "black"))

p <- p + theme(legend.position = c(0.1,0.8)) + ggplot2::xlim(0, 0.03)


p

get_taxa_name(p)

phylo.effects <- model$Sol[, grep("phylo", colnames(model$Sol))]
phylo.effects <- as.data.frame(phylo.effects)
phylo.effects.long <- pivot_longer(phylo.effects, cols = everything(), names_to = "tip", values_to = "value")
phylo.effects.long$tip <- gsub("phylo.", "", phylo.effects.long$tip)
phylo.effects.long$tip <- factor(phylo.effects.long$tip, levels = get_taxa_name(p))
phylo.effects.long <- phylo.effects.long %>%
  group_by(tip) %>%
  mutate(value = value) %>%
  mutate(hpd = ci(value, method = "HDI")) %>%
  mutate(mean = mean(value)) %>%
  select(tip, mean, hpd) %>%
  distinct()
phylo.effects.long <- cbind(phylo.effects.long, phylo.effects.long$hpd)

p.2 <- ggplot(phylo.effects.long, aes(x = tip, y = mean)) +
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=0) +
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
  ylab(label="Posterior 95% HDI")

p.2 <- p.2 + theme(legend.position = c(-1.25,0.9))

library(cowplot)
plot_grid(p, NULL, p.2, align='hv', 
          rel_widths = c(1, 0, 0.6), nrow=1,
          labels=c("a", "b"))


####


phylo.effects <- model$Sol[, grep("phylo", colnames(model$Sol))]
phylo.effects <- as.data.frame(phylo.effects)
phylo.effects.long <- pivot_longer(phylo.effects, cols = everything(), names_to = "tip", values_to = "value")

anova_result <- aov(value ~ tip, data = phylo.effects.long)
summary(anova_result)

tukey_result <- TukeyHSD(anova_result)
print(tukey_result)
tukey.df <- as.data.frame(tukey_result$tip)
tukey.df$sig <- ifelse(tukey.df$`p adj` <0.001, TRUE, FALSE)
summary(tukey.df$sig)


