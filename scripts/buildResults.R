###
### Import libraries and data
###

library(readr)
library(vcfR)
library(dplyr)
library(tidyr)
library(reshape2)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(ordinal)
library(lme4)
library(clevr)
library(cowplot)
library(ggtext)

setwd('/Users/willmatlock/Desktop/tem1_2024')

# Main metadata file
df <- read_delim("metadata/metadata.tsv", 
                 delim = "\t", escape_double = FALSE, 
                 trim_ws = TRUE)
# Efflux annotations
efflux <- read_delim("metadata/amrfinder_report.tsv", 
                     delim = "\t", escape_double = FALSE, 
                     trim_ws = TRUE)
# blaTEM-1 and linked promoter annotations
tem1.report <- read_delim("metadata/tem1_report_final.csv", 
                          delim = "\t", escape_double = FALSE, 
                          trim_ws = TRUE)

###
### A curated sample of E. coli isolates with hybrid assemblies and co-amoxiclav MICs
###

# Sample size
length(unique(as.factor(df$isolate.id)))
# Sampling date range
summary(as.Date(df[df$contig.type=="chromosome",]$date, format = "%d/%m/%Y"))
# Complete assemblies
summary(as.factor(df[df$contig.type=="chromosome",]$is.complete))
# Plasmid stats
df %>% 
  filter(contig.type != "chromosome") %>%
  group_by(isolate.id) %>%
  dplyr::summarise(n.contigs=n()) %>%
  dplyr::summarise(median.plasmids = median(n.contigs),
                   lq.plasmids = quantile(n.contigs, 0.25),
                   uq.plasmids = quantile(n.contigs, 0.75))

# blaTEM-1 isolate copies
summary(df[df$tem1.isolate>0,]$tem1.isolate)
# blaTEM-1-positive contigs
nrow(df[df$tem1.replicon>0,])
nrow(df[df$tem1.replicon>1,])
# blaTEM-1 contig types
summary(as.factor(df[df$tem1.replicon>0,]$contig.type))
# blaTEM-1 SNV profiles
summary(as.factor(tem1.report$tem1.snv))
# blaTEM-1 repeated SNV profiles
tem1.report %>%
  group_by(isolate.assembly, contig.assembly) %>%
  dplyr::summarise(n.tem1 = n(),
            n.tem1.snv = n_distinct(tem1.snv)) %>%
  filter(n.tem1 > 1)

# Sequence length vs blaTEM-1 carriage
sum(df[df$contig.type=="plasmid",]$contig.length)
sum(df[df$contig.type=="chromosome",]$contig.length)
sum(df[df$contig.type=="plasmid",]$tem1.replicon)
sum(df[df$contig.type=="chromosome",]$tem1.replicon)

# Phylogroup stats
summary(as.factor(df[df$contig.type=="chromosome",]$ezclermont.phylogroup))

# Promoter stats
nrow(tem1.report[!is.na(tem1.report$promoter.snv),])
summary(as.factor(tem1.report$promoter.snv))


# MIC stats
summary(as.factor(df[df$contig.type=="chromosome",]$coamox.mic))

###
### Plot Figure 1
###

add_na_to_list <- function(lst) {
  if (length(lst) == 0 || !grepl("^chromosome-", lst[1])) {
    c(NA, lst)
  } else {
    lst
  }
}

df.fig1 <- df %>%
  filter(tem1.replicon>0) %>%
  select(isolate.id, contig.type, tem1.replicon) %>%
  group_by(isolate.id, contig.type) %>%
  mutate(contig.copies = paste(contig.type, tem1.replicon, sep="-")) %>%
  arrange(desc(tem1.replicon), contig.type) %>%
  group_by(isolate.id) %>%
  mutate(distribution = list(contig.copies)) %>%
  distinct(isolate.id, distribution) %>%
  mutate(distribution = lapply(distribution, add_na_to_list)) %>%
  tidyr::unnest_wider(distribution, names_sep = "_") %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate_all(~str_replace_all(., "chromosome-|plasmid-", "")) %>%
  rename(chromosome = distribution_1,  plasmid.1 = distribution_2,
           plasmid.2 = distribution_3, plasmid.3 = distribution_4)

df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  dplyr::summarise(count = n()) %>%
  arrange(desc(count))

df.fig1 <- merge(df.fig1, df[c("isolate.id", "ezclermont.phylogroup", "coamox.mic")], by = "isolate.id", all.x= TRUE)
df.fig1 <- unique(df.fig1[!is.na(df.fig1$ezclermont.phylogroup),])
    
df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df.fig1.summary <- df.fig1.summary[,-c(5)]
df.fig1.summary$arrangement.id <- nrow(df.fig1.summary):1
df.fig1.summary.melt <- melt(df.fig1.summary, id.vars = "arrangement.id")

# Plot TEM-1 arrangements for Figure 1a

p.fig1.a <- ggplot(df.fig1.summary.melt, aes(x = variable, y = arrangement.id, fill = value)) +
  geom_tile(color = "black",  width = 0.9, height = 0.9, size=0.5) +
  scale_fill_manual(values = c("white", rep("black", 4))) +
  geom_text(aes(label = ifelse(value != 0, value, "")), color = "white") +
  theme_minimal() +
  theme(
    legend.position = "none",  
    panel.grid = element_blank(), 
    axis.title.y = element_blank(),  
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1.5, color = "black")) +
  scale_x_discrete(labels = c("chromosome" = "Chromosome", 
                              "plasmid.1" = "Plasmid 1", "plasmid.2" = "Plasmid 2", "plasmid.3" = "Plasmid 3")) +  # Rename x-axis labels
  coord_fixed(ratio = 1, ylim = c(0.5, 12.5)) 

df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df.fig1.summary$arrangement.id <- nrow(df.fig1.summary):1
df.fig1.summary <- df.fig1.summary[c("count", "arrangement.id")]
df.fig1.summary$count <- as.numeric(df.fig1.summary$count)

# Plot TEM-1 arrangement counts for Figure 1b

p.fig1.b <- ggplot(df.fig1.summary, aes(x=1, y = arrangement.id, fill = count)) +
  geom_tile(color = "black",  width = 0.9, height = 0.9, size=0.5) +
  geom_text(aes(label = ifelse(count != 0, count, "")), color = "Black") +
  theme_minimal() +
  scale_fill_distiller(palette = "OrRd", direction=1) +
  theme(
    legend.position = "none",  
    panel.grid = element_blank(), 
    axis.title.y = element_blank(),  
    axis.title.x = element_blank(), 
    axis.text.y = element_blank(),  
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank()) +
  coord_fixed(ratio = 1, ylim = c(0.5, 12.5))

df.fig1.summary <- df.fig1 %>%
  group_by(chromosome, plasmid.1, plasmid.2, plasmid.3) %>%
  summarise(count = n()) %>%
  arrange(desc(count))

df.fig1.summary <- df.fig1.summary[,-c(5)]
df.fig1.summary$arrangement.id <- nrow(df.fig1.summary):1

df.fig1[, c("chromosome", "plasmid.1", "plasmid.2", "plasmid.3")] <- lapply(df.fig1[, c("chromosome", "plasmid.1", "plasmid.2", "plasmid.3")], as.character)
df.fig1 <- merge(df.fig1, df.fig1.summary, by=c("chromosome", "plasmid.1", "plasmid.2", "plasmid.3"), all.x=TRUE)

# Plot TEM-1 arrangement group phylogroups for Figure 1c

color_palette <- brewer.pal(n = 9, name = "Set3")

p.fig1.c <- ggplot(df.fig1, aes(x = reorder(interaction(chromosome, plasmid.1, plasmid.2, plasmid.3), arrangement.id), fill = ezclermont.phylogroup)) +
  geom_bar(position = "fill", show.legend = FALSE) +
  geom_bar(position = "fill", color = "black", size = 0.5) +
  scale_fill_manual(values = color_palette) +
  labs(y = "Combination", x = "", fill = "Phylogroup") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_blank(),  # Remove background color
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1.5, color = "black")
  ) +
  coord_flip() +
  theme(legend.position = "none")

# Plot TEM-1 arrangement group co-amoxicillin MIC for Figure 1c

df.fig1$coamox.mic <- factor(df.fig1$coamox.mic, levels = c("<=2.2", "4.2", "8.2", ">8.2",
                                                            "16.2", "32.2", ">32.2"))

color_palette <- brewer.pal(n = 7, name = "PuBu") 

p.fig1.d <- ggplot(df.fig1, aes(x = reorder(interaction(chromosome, plasmid.1, plasmid.2, plasmid.3), arrangement.id), fill = coamox.mic)) +
  geom_bar(position = "fill", show.legend = FALSE) +
  geom_bar(position = "fill", color = "black", size = 0.5) +
  scale_fill_manual(values = color_palette) +
  labs(y = "Combination", x = "", fill = "Co-amoxiclav MIC") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),  # Remove background grid
    panel.background = element_blank(),  # Remove background color
    axis.text.y = element_blank(),
    axis.title.x = element_blank(), 
    axis.text.x = element_text(angle = 45, hjust=1, vjust=1.5, color = "black")
  ) +
  coord_flip() +
  theme(legend.position = "none")

###
### blaTEM-1 associated with conjugative plasmids.
###

df.plasmid <- df[df$contig.type=="plasmid" & df$contig.topology=="circular",]

df.plasmid$predicted.mobility <- factor(df.plasmid$predicted.mobility, levels = c("non-mobilizable",
                                                                                  "mobilizable", "conjugative"))
df.plasmid$carries.tem1 <- ifelse(df.plasmid$tem1.replicon>0, 1, 0)
summary(as.factor(df.plasmid$carries.tem1))

tem1.counts <- c(200, 286)
totals <- c(117 + 386 + 200, 20 + 27 + 286)
prop.test.result <- prop.test(tem1.counts, totals)
print(prop.test.result)

df.plasmid.2 <- df[df$contig.type=="plasmid",]
summary(df.plasmid.2$contig.copy.number)
nrow(df.plasmid.2[df.plasmid.2$contig.copy.number<1,])

df.plasmid.3 <- df.plasmid[df.plasmid$contig.copy.number>=1,]
summary(df.plasmid.3[df.plasmid.3$contig.length>10000,]$contig.copy.number)
summary(df.plasmid.3[df.plasmid.3$contig.length<=10000,]$contig.copy.number)

p.fig.s1 <- ggplot(df.plasmid.3, aes(x = contig.copy.number, y = contig.length)) +
  geom_density_2d(h=NULL, colour='#ef8a62', size=1) +
  geom_point(alpha=0.5, size=2, stroke=0) +
  labs(x = "Plasmid copy number (log<sub>10</sub>)",
       y = "Plasmid length (bp, log<sub>10</sub>)") +
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
        axis.title.y = element_markdown()) +
  scale_y_continuous(
    trans = "log10",
    breaks = 10^seq(4,5),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_x_continuous(
    trans = "log10",
    breaks = c(0.1, 1, 10),
    labels = scales::trans_format("log10", scales::math_format(10^.x)))

###
### E. coli phylogeny structures ampC diversity.
###

df.ampc <- df %>%
  filter(contig.type=="chromosome")

homogeneity(df.ampc$ampc.snv, df.ampc$ezclermont.phylogroup)
completeness(df.ampc$ampc.snv, df.ampc$ezclermont.phylogroup)

sum(df.ampc$mlst.st=="-")

df.ampc <- df.ampc %>%
  filter(mlst.st != "-")

homogeneity(df.ampc$ampc.snv, df.ampc$mlst.st)
completeness(df.ampc$ampc.snv, df.ampc$mlst.st)

# promoter info
summary(as.factor(df[df$contig.type=="chromosome",]$ampc.promoter.snv))

promoter.mutations <- df %>%
  filter(contig.type == "chromosome") %>%
  select(isolate.id, ampc.promoter.snv) %>%
  mutate(position = row_number()) %>%
  separate(ampc.promoter.snv, into = paste0("position_", 1:12), sep = "") %>%
  select(isolate.id, paste0("position_", 2:12)) %>%
  mutate_all(as.factor)
summary(promoter.mutations)

homogeneity(df.ampc$ampc.snv, df.ampc$ampc.promoter.snv)
completeness(df.ampc$ampc.snv, df.ampc$ampc.promoter.snv)

###
### Efflux results
###

df.efflux <- efflux %>%
  filter(Class=="EFFLUX") %>%
  mutate(cov = as.numeric(`% Coverage of reference sequence`),
         id = as.numeric(`% Identity to reference sequence`)) %>%
  filter(cov >= 95 & id >=95) %>%
  select(`1`, `Contig id`, `Gene symbol`, `Sequence name`) %>%
  filter(`1` %in% unique(df$isolate.assembly))
colnames(df.efflux) <- c("isolate.assembly", "contig.assembly", "gene", "description")

df.efflux <- merge(df.efflux, df, by=c("isolate.assembly", "contig.assembly"), all.x=TRUE)

df.efflux %>% 
  group_by(isolate.assembly) %>%
  summarise(n.efflux = n()) %>%
  summary()
  
df.efflux.acrf <- df.efflux %>%
  filter(gene=="acrF") %>% 
  group_by(isolate.assembly) %>%
  mutate(n.efflux = n())

summary(as.factor(df.efflux.acrf$ezclermont.phylogroup))/nrow(df.efflux.acrf)* 100
  
  
  
  
  
  
  
  
