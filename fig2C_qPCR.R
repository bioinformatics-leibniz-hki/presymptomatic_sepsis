# code to generale Figure 2C
# necessary data is available upon request

# author: Albert Garcia Lopez

pacman::p_load(
  devtools,
  tidyverse,
  BiocManager,
  tibble,
  gdata,
  ggfortify,
  ggpubr,
  cowplot,
  ggsci,
  enrichR,
  gplots,
  RColorBrewer,
  data.table,
  openxlsx,
  readxl,
  readr,
  viridis,
  extrafont,
  org.Hs.eg.db
  )

# Select the common selected genes ----------------------------------------

genes <- read.xlsx("./dat/common-genes_and_qPCRvalues_MA.xlsx", sheet = 1)

# Load qPCR data ----------------------------------------------------------

qpcr <- read.xlsx("./dat/common-genes_and_qPCRvalues_MA.xlsx", sheet = 2)
is.na(qpcr) %>% sum()

# Select only those that Microarray information. qpcr$MA == "Y"

qpcr$MA %>% unique
qpcr <- qpcr %>% filter(MA == "Y")
qpcr$MA %>% unique

# Select only relevant columns

qpcr <- qpcr[, c(1, 2, 4, 7:ncol(qpcr))]
qpcr <- qpcr %>% dplyr::rename(SampleID = sample.ID, Condition = S_C)

# Remove "I" since it is SIRS! --> NOT PRESENT HERE

grep("SIRS", qpcr$Class)
unique(qpcr$Condition)
qpcr$Condition <- as.factor(qpcr$Condition)
levels(qpcr$Condition)
levels(qpcr$Condition) <- c("Control", "Infection")
levels(qpcr$Condition)
qpcr <- as.data.frame(qpcr)

# Select DEGs ----------------------------------

colnames(qpcr)[which(colnames(qpcr) %in% genes$name.of.column)] %>% length
colnames(qpcr)[which(colnames(qpcr) %in% genes$Symbol)] %>% length

qpcr <- qpcr %>% dplyr::select(SampleID, Condition, genes$name.of.column)
qpcr <- na.omit(qpcr)


grep("\\.DCq", colnames(qpcr))
colnames(qpcr) <- gsub("\\.DCq", "", colnames(qpcr))

table(qpcr$Condition) # Ratio 1:1 between Control (116) and Infection (124) patients

# Check if they are normally distributed
norm_qpcr <- apply(qpcr[,3:ncol(qpcr)], 2, shapiro.test) # Apply a Shapiro-Wilk test to each column of the data set.
norm_qpcr <- sapply(norm_qpcr, `[`, "p.value") # Just get the p-value resulting from the Shapiro-Wilk test
norm_qpcr <- unlist(norm_qpcr) # Turn list values into a numeric vector of p-values.
adj_qpcr <- p.adjust(norm_qpcr, method = "fdr", n = length(norm_qpcr)) # adjust for multiple comparisons by FDR correction.
which(adj_qpcr <= 0.05) # 6/8 values are smaller than 0.05. Use Wilcoxon-Test.

nrow(qpcr)

qpcr_final <- qpcr %>%
  tidyr::gather(Genes, Value, 3:ncol(.))

# Create barplot ----------------------------------------------------------

# To tune the stars showing p-value significance from 4 to 3 stars.

symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), 
                    symbols = c("***", "**", "*", "x", "ns"))

# Plot now the barplot

barplot <-
  ggbarplot(
    qpcr_final,
    x = "Genes",
    y = "Value",
    add = "mean_se",
    color = "black",
    position = position_dodge(0.75),
    fill = "Condition",
    palette = c("#6996AA", "#C7533B")
  ) +
  stat_compare_means(
    aes(group = Condition),
    label = "p.signif",
    method = "wilcox.test",
    symnum.args = symnum.args,
    size = 10 # To make significant labels bigger
  ) +
  geom_abline(
    slope = 0,
    intercept = 0,
    col = "black",
    lty = 2
  ) +
  labs(y = expression(Delta*"Cq")) +
  theme(
    text = element_text(size = 20),
    axis.title.x = element_blank(),
    legend.title = element_blank(),
    axis.title.y = element_text(face = "bold", size = 20, vjust = 1.8),
    legend.text = element_text(size = 20),
    legend.position = "bottom",
    axis.text.x = element_text(face = "italic")
  )

barplot

ggsave(filename = "./res/Fig_2C_barplot.pdf",
       plot = barplot,
       scale = 1,
       width = 50,
       height = 20,
       units = "cm")
ggsave(filename = "./res/Fig_2C_barplot.svg",
       plot = barplot,
       device = "svg",
       scale = 1,
       width = 50,
       height = 20,
       units = "cm")

# To get the table of p-values / FDR-adjusted p-values

signf_table <- compare_means(Value ~ Condition, 
                             data = qpcr_final, 
                             group.by = "Genes", 
                             ref.group = "Control", 
                             symnum.args = symnum.args,
                             method = "wilcox.test", 
                             p.adjust.method = "fdr")
write.csv(signf_table, "./res/p-values_significance_table.csv")

# Enrichment analysis of selected genes -------------------------------------------------------

# In clusterProfiler, groupGO is designed for gene classification based on GO distribution at a specific level.

# CC

cc_enr <- enrichGO(
  gene     = unique(qpcr_final$Genes),
  OrgDb    = org.Hs.eg.db,
  keyType  = "SYMBOL",
  ont      = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
cc_table <- data.frame(cc_enr@result)
cc_table <- cc_table[grep("/", cc_table$geneID),]

# MF

mf_enr <- enrichGO(
  gene     = unique(qpcr_final$Genes),
  OrgDb    = org.Hs.eg.db,
  keyType  = "SYMBOL",
  ont      = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
mf_table <- data.frame(mf_enr@result)
mf_table <- mf_table[grep("/", mf_table$geneID),]

# BP

bp_enr <- enrichGO(
  gene     = unique(qpcr_final$Genes),
  OrgDb    = org.Hs.eg.db,
  keyType  = "SYMBOL",
  ont      = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
bp_table <- data.frame(bp_enr@result)
bp_table <- bp_table[grep("/", bp_table$geneID),]

# KEGG

# It is needed to convert SYMBOLS to ENTREZ ID first

entrez_degs <- mapIds(org.Hs.eg.db, unique(qpcr_final$Genes), 'ENTREZID', 'SYMBOL')
entrez_degs

kegg_enr <- enrichKEGG(
  gene = entrez_degs,
  organism = "hsa",
  keyType = "kegg",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
kegg_table <- as.data.frame(kegg_enr@result)
kegg_table <- kegg_table[grep("/", kegg_table$geneID),]

# To select enriched genes in ENTREZ format and convert it to SYMBOL
kegg_split <- kegg_table$geneID[grep("/", kegg_table$geneID)]
kegg_split <- strsplit(kegg_split, split = "/") %>% unlist()
kegg_trans <- bitr(kegg_split, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
kegg_table$SYMBOL <- paste(kegg_trans$SYMBOL, sep = "/", collapse = "/")

# Export all results into the same Excel file
library(openxlsx)
tables_enr <- createWorkbook()
addWorksheet(tables_enr, "BP")
addWorksheet(tables_enr, "CC")
addWorksheet(tables_enr, "MF")
addWorksheet(tables_enr, "KEGG")
writeData(tables_enr, 1, bp_table)
writeData(tables_enr, 2, cc_table)
writeData(tables_enr, 3, mf_table)
writeData(tables_enr, 4, kegg_table)
saveWorkbook(tables_enr, file = "./res/enrichment_tables_clusterprofiler.xlsx", overwrite = TRUE)
