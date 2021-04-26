# code to generate figure 2B
# necessary data is available upon request
#
# author Sascha Schäuble

library("tidyverse")
library("magrittr")
library("GOplot")
library("ggpubr")
library("pheatmap")

PROJ_PATH <- "."
DAT_PATH <- paste0(PROJ_PATH, "dat/")
RES_PATH <- paste0(PROJ_PATH, "res/")

FN1 <- "clusterProfiler_Supplementary_File_S6.csv" 
FN2 <- "sel_degs.tsv"

# ================================================================ #



#### Data wrangling ################################################
dat.enrich <- read_delim(paste0(DAT_PATH, FN1), delim = ",")
# put to expected structure
dat.enrich.v2 <- dat.enrich %>% select(term_id, term_name, intersections, adjusted_p_value)
dat.enrich.v2 %<>% mutate(Category = "BP") 
dat.enrich.v2 %<>% select(Category, term_id, term_name, intersections, adjusted_p_value)
colnames(dat.enrich.v2) <- c("Category", "ID", "Term", "Genes", "adj_pval")

summary(dat.enrich.v2)

dat.genes <- read_tsv(paste0(DAT_PATH, FN2))
dat.genes.v2 <- dat.genes %>% select(Symbol, `fold change`, min_pvalue)
colnames(dat.genes.v2) <- c("ID", "logFC", "adj.P.Val") # of note with MA data we already should have log fold changes, since cycles get subtracted

dat.clustPrfler <- read_tsv(paste0(DAT_PATH, FN1))
dat.clustPrfler.v2 <- dat.clustPrfler %>% select(ID, Description, geneID, p.adjust)
dat.clustPrfler.v2 %<>% mutate(Category = "BP") 
dat.clustPrfler.v2 %<>% select(Category, ID, Description, geneID, p.adjust)
colnames(dat.clustPrfler.v2) <- c("Category", "ID", "Term", "Genes", "adj_pval")
summary(dat.clustPrfler.v2)
dat.clustPrfler.v2$Genes %<>% str_replace_all(pattern = "/", replacement = ",")
# ================================================================ #



#### Go plots ######################################################
dat.circ <- circle_dat(terms = dat.clustPrfler.v2[1:200,], genes = dat.genes.v2)
dim(dat.circ)
# the following performs poorly, but reduces redundant GO categories
dat.circ.red <- reduce_overlap(dat.circ, overlap = 0.33)
dim(dat.circ.red)
# GOBar
GOBar(dat.circ.red)

# bubble
p <- GOBubble(dat.circ.red, labels = 5, ID = F, table.legend = F)
p %>%
  ggexport(filename = paste0(RES_PATH, "bubbleChart.pdf"), width = 8, height = 5)
