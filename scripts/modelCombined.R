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

df <- read_delim("path/to/supplementary_table_1.tsv", 
                 delim = "\t", escape_double = FALSE,
                 trim_ws = TRUE)
tem1.report <- read_delim("path/to/supplementary_table_2.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

phylo <- read.tree("path/to/GTR_F_I_R4.treefile")

exp <- read_excel("path/to/supplementary_table_3.xlsx")

### prepare MIC data

df.mic <- df %>% 
  mutate(tem1.contig.copy.number = contig.copy.number * tem1.replicon) %>%
  group_by(isolate.id) %>%
  mutate(tem1.isolate.copy.number = sum(tem1.contig.copy.number)) %>%
  ungroup() %>%
  mutate(tem1.isolate.copy.number.scaled = (tem1.isolate.copy.number - mean(tem1.isolate.copy.number)) / sd(tem1.isolate.copy.number),
         tem1.isolate.scaled = ifelse(tem1.isolate>1, TRUE, FALSE)) %>%
  mutate(tem1.isolate.scaled = as.factor(tem1.isolate.scaled)) %>%
  filter(contig.type=="chromosome") %>%
  ungroup() %>%
  select(isolate.id, isolate.assembly, coamox.mic, tem1.isolate.scaled, tem1.isolate.copy.number.scaled) %>%
  mutate(coamox.mic = factor(coamox.mic, ordered = TRUE, levels = c("<=2.2", "4.2", "8.2",
                                                                    "16.2", "32.2", ">32.2"))) 

upper_percentile_value <- quantile(df.mic$tem1.isolate.copy.number.scaled, 0.95)
df.mic$tem1.isolate.copy.number.scaled[df.mic$tem1.isolate.copy.number.scaled > upper_percentile_value] <- upper_percentile_value

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

df.mic <- merge(df.mic, tem1.model, by="isolate.assembly", all.x=TRUE) %>%
  filter(!is.na(promoter.snv)) 

summary(as.factor(df.mic$promoter.snv))
df.mic$promoter.snv = factor(df.mic$promoter.snv, ordered = FALSE, 
                               levels = c("CGGCGG", "CGGCGA", "TGGCGA",
                                          "TGGCGG"))

df.mic$trait <- "mic"
df.mic$family <- "gaussian"

df.mic <- as.data.frame(df.mic)
df.mic$phylo <- as.factor(df.mic$isolate.assembly)

### prepare expression data

df.exp <- df %>%
  filter(tem1.isolate==1 & tem1.replicon==1) %>%
  select(accession, isolate.assembly, contig.assembly)

exp$Isolate <- sapply(strsplit(as.character(exp$Isolate), " "), function(x) x[1])
exp <- exp[c("Isolate", "delta Ct (control)")]
colnames(exp) <- c("accession", "exp")

df.exp <- merge(df.exp, exp, by="accession", all.y=TRUE)

df.exp <- df.exp %>% filter(!is.na(exp)) %>% filter(!is.na(isolate.assembly)) %>% distinct()
df.exp <- df.exp %>% group_by(accession) %>% mutate(run = row_number())

mean.exp <- mean(df.exp$exp)
sd.exp <- sd(df.exp$exp)

df.exp <- df.exp %>%
  mutate(exp.scaled = (exp - mean.exp) / sd.exp)

up.exp <- quantile(df.exp$exp.scaled, 0.95)
df.exp$exp.scaled[df.exp$exp.scaled > up.exp] <- up.exp

df.exp$trait <- "exp"
df.exp$family <- "gaussian"

### wrangle

df.model <- merge(df.mic, df.exp, all=TRUE, by=c("isolate.assembly", "trait", "family"))

df.model <- df.model %>%
  group_by(isolate.assembly) %>%
  filter(!(!("mic" %in% trait) & "exp" %in% trait)) %>%
  fill(isolate.id, .direction = "downup") %>%
  fill(tem1.isolate.scaled, .direction = "downup") %>%
  fill(tem1.isolate.copy.number.scaled, .direction = "downup") %>%
  fill(promoter.snv, .direction = "downup") %>%
  fill(phylo, .direction = "downup") %>%
  ungroup() %>%
  select(isolate.id, phylo, trait, family, tem1.isolate.scaled, tem1.isolate.copy.number.scaled,
         promoter.snv, coamox.mic, exp.scaled, run)

df.model$y <- ifelse(!is.na(df.model$coamox.mic), as.factor(df.model$coamox.mic), as.numeric(df.model$exp.scaled))

df.model$trait <- as.factor(df.model$trait)

### prepare phylogeny

df.model$phylo <- as.character(df.model$phylo)
phylo$tip.label <- as.character(phylo$tip.label)
phylo$node.label <- NULL

phylo <- keep.tip(phylo, df.model$phylo)

phylo <- midpoint.root(phylo)

phylo.u <- chronos(phylo, lambda=1, model="correlated")

is.rooted(phylo.u)
is.ultrametric(phylo.u)

inv.phylo <- inverseA(phylo.u,nodes="TIPS",scale=TRUE)

### run model

df.model <- as.data.frame(df.model)

prior <- list(G=list(G1=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1e+3),
                     G2=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1e+3),
                     G3=list(V=diag(1),nu=1, alpha.mu = 0, alpha.V = 1e+3)),
              R = list(V=diag(2), nu=0.002),#, fix=2),
              S = list(mu=0, V=1e+3))

set.seed(1)
chain.1 <- MCMCglmm(y ~ -1 + trait:(1 + tem1.isolate.scaled + tem1.isolate.copy.number.scaled + promoter.snv),
                  random= ~ phylo + isolate.id + us(at.level(trait, "mic")):phylo,
                  rcov= ~ idh(trait):units,
                  family=NULL,
                  theta_scale = list(factor="trait", level="mic", random=1:2),
                  ginverse=list(phylo=inv.phylo$Ainv),
                  data=df.model,
                  prior=prior,
                  nitt=10000000,
                  burnin=1000000,
                  thin=100,
                  DIC=FALSE,
                  pr=TRUE)

set.seed(2)
chain.2 <- MCMCglmm(y ~ -1 + trait:(1 + tem1.isolate.scaled + tem1.isolate.copy.number.scaled + promoter.snv),
                    random= ~ phylo + isolate.id + us(at.level(trait, "mic")):phylo,
                    rcov= ~ idh(trait):units,
                    family=NULL,
                    theta_scale = list(factor="trait", level="mic", random=1:2),
                    ginverse=list(phylo=inv.phylo$Ainv),
                    data=df.model,
                    prior=prior,
                    nitt=10000000,
                    burnin=1000000,
                    thin=100,
                    DIC=FALSE,
                    pr=TRUE)

summary(chain.1)
plot(chain.1)

summary(chain.2)
plot(chain.2)

mclist <- mcmc.list(chain.1$Sol, chain.2$Sol)
gelman.diag(mclist)

