# code to generate figure 2A 
# necessary data is available upon request

# author: Albert Garcia Lopez


pacman::p_load(
  tidyverse,
  data.table,
  gdata,
  ggpubr,
  ggsci,
  gplots,
  pheatmap,
  RColorBrewer,
  readxl,
  readr,
  openxlsx,
  viridis,
  Cairo
)

# Load logFC data info and select 1500 DEGs -------------------------------

FC_logFC_values <- read_excel("dat/Groups_of_Degs_Illumina_dataset.xlsx")
FC_logFC_values <- FC_logFC_values[, c(2:10)]
colnames(FC_logFC_values) <- FC_logFC_values[1,]
FC_logFC_values <- FC_logFC_values[-1,]

colnames(FC_logFC_values)[4:ncol(FC_logFC_values)]
colnames(FC_logFC_values)[4:ncol(FC_logFC_values)] <- c("Fold3", "pval3", "Fold2", "pval2", "Fold1", "pval1")
colnames(FC_logFC_values)[4:ncol(FC_logFC_values)]
FC_logFC_values[, 4:ncol(FC_logFC_values)] <- apply(FC_logFC_values[,4:ncol(FC_logFC_values)], 2, as.numeric)
stopifnot(apply(FC_logFC_values[,4:ncol(FC_logFC_values)], 2, class) == "numeric")

# Calculate the mean logFC
# We use the Mean Log FC to select the genes

FC_logFC_values <- FC_logFC_values %>%
  dplyr::mutate(MeanFC = rowMeans(subset(FC_logFC_values, select = c(Fold1, Fold2, Fold3)), na.rm = T)) 

# Select only the genes with a p-value <= 0.05

FC_logFC_values$min_pvalue <- apply(FC_logFC_values %>% select(pval3, pval2, pval1), 1, min)
length(which(FC_logFC_values$min_pvalue <= 0.05)) # 8692 genes with a p-value below <= 0.05
FC_logFC_values <- FC_logFC_values %>% filter(min_pvalue <= 0.05)

# Select only the relevant columns
degs <- FC_logFC_values
length(which(duplicated(degs$Symbol))) # 1384 Symbols have different illuID for the same name.

test <- degs[order(degs$Symbol),] # QC

# Remove duplicated instances but keeping the highest absolute values of each pair.
# 1) Sort in the order putting the less desired items last within the 'illuIDs' groups.
degs <- degs[order((degs$Symbol), -abs(degs$MeanFC)),]
# 2) Remove items afterthe first within 'illuIDs' groups
degs <- degs[!duplicated(degs$Symbol),]
unique(duplicated(degs))

# Select 1500 Symbols with highest absolute Sum_FC value

length(degs$Symbol) # 7308 Symbols in total
write.csv(degs, "./res/all_list_of_degs_Fig2A.csv")
degs <- degs %>% dplyr::mutate(Abs = abs(MeanFC)) %>% arrange(desc(Abs))
sel_degs <- degs[1:1500,]
write.csv(sel_degs, "./res/sel_degs.csv")

# How many of the 1500 are upregulated/downregulated. Export too.

upreg_sel_degs <- sel_degs %>% filter(MeanFC > 0)
upreg_sel_degs$Symbol %>% length # 863 upregulated
write.csv(upreg_sel_degs, "./res/upreg_sel_degs.csv")

down_sel_degs <- sel_degs %>% filter(MeanFC < 0)
down_sel_degs$Symbol %>% length # 637 downregulated
write.csv(down_sel_degs, "./res/down_sel_degs.csv")

# How many DEGs (whole group) above 1.2 and 1.3 Abs Mean FC

degs_1.2_or_more <- degs %>% filter(Abs >= 1.2) # 2594 DEGs
degs_1.2_or_more$Symbol %>% length # 2594 DEGs
write.csv(degs_1.2_or_more, "./res/degs_1.2_or_more.csv")

degs_1.3_or_more <- degs %>% filter(Abs >= 1.3) # 2594 DEGs
degs_1.3_or_more$Symbol %>% length # 2594 DEGs
write.csv(degs_1.3_or_more, "./res/degs_1.3_or_more.csv")


# Load Illumina Expression Data -------------------------------------------

illumina <- read.table("./dat/Expression_Illu.txt", sep = "")
illumina_2 <- illumina %>% rownames_to_column("IL-ID")
illu <- left_join(sel_degs, illumina_2, by = "IL-ID")
apply(illu[,13:ncol(illu)], 2, class) %>% unique # All added columns are already numeric

illu_filtered <- illu %>%
  select(Symbol, 13:ncol(.)) %>%
  column_to_rownames("Symbol") %>%
  t(.) %>%
  as.data.frame(.)

is.na(illu_filtered) %>% sum
rownames(illu_filtered) %>% duplicated %>% unique

# Scale expression data
# Z-score has to be column wise, which means that genes are in columns.

apply(illu_filtered, 2, class) %>% unique # All columns (genes) in the expression dataset are numeric.
expr_scaled <- scale(illu_filtered, center = TRUE, scale = TRUE) # Genes have to be in columns, samples have to be in rows.
apply(expr_scaled, 2, mean) %>% mean %>% round(., 5) # Has to be close to 0
apply(expr_scaled, 2, sd) %>% mean # Has to be 1

# Remove NAs (if any)
is.na(expr_scaled) %>% sum # No NAs found

# Represent 90% CI only (NOT APPLIED!)
quantile(expr_scaled, probs = c(0.05, 0.95))
expr_scaled[expr_scaled < -1.622206] <- -1.622206
expr_scaled[expr_scaled > 1.634455] <- 1.634455

# Load metadata -----------------------------------------------------------

# pheno192_388_ANNO_corrP1P2.txt: meta data for expression values
# Special interest -> column kzID:
# S1 - Infection Day-1
# S2 - Infection Day-2
# S3 - Infection Day-3
# C1 - Control Day-1
# C2 - Control Day-2
# C3 - Control Day-3

# Important: only MetaYes = 1

metadata <- read_delim("./dat/pheno_anno.txt", "\t",
                       escape_double = FALSE,
                       trim_ws = TRUE)

# Select those with MetaYes = 1

length(which(metadata$MetaYes == 1))
metadata$MetaYes %>% unique
metadata <- metadata %>% filter(MetaYes == 1)
metadata$MetaYes %>% unique

# Select those that are either S1, S2, S3, C1, C2, C3. Change column name to Timepoint and rename values.

timepoint_selection <- c('S1', 'S2', 'S3', 'C1', 'C2', 'C3')
length(grep('S1|S2|S3|C1|C2|C3', metadata$kzID))
metadata <- metadata %>% filter(kzID %in% timepoint_selection) %>% rename(Timepoint = kzID)
metadata$Timepoint <- gsub("S1", "Day -1", metadata$Timepoint)
metadata$Timepoint <- gsub("S2", "Day -2", metadata$Timepoint)
metadata$Timepoint <- gsub("S3", "Day -3", metadata$Timepoint)
metadata$Timepoint <- gsub("C1", "Day -1", metadata$Timepoint)
metadata$Timepoint <- gsub("C2", "Day -2", metadata$Timepoint)
metadata$Timepoint <- gsub("C3", "Day -3", metadata$Timepoint)
metadata$Timepoint <- as.factor(metadata$Timepoint)

# Condition column

metadata <- metadata %>% rename(Condition = Typ)
metadata$Condition <- as.factor(metadata$Condition)
levels(metadata$Condition)
levels(metadata$Condition) <- c("Control", "Infection")
levels(metadata$Condition)

# Select only relevant columns

metadata_selection <- metadata %>%
  dplyr::select(sample_name, Condition, Timepoint) %>%
  rename(ExtBarCode = sample_name)

# Relevel Condition and Timepoint

levels(metadata_selection$Condition)
metadata_selection$Condition <- fct_relevel(metadata_selection$Condition, "Infection")
levels(metadata_selection$Condition)

levels(metadata_selection$Timepoint)
metadata_selection$Timepoint <- fct_relevel(metadata_selection$Timepoint, "Day -3", "Day -2", "Day -1")
levels(metadata_selection$Timepoint)

# Arrange by Condition and Timepoint
metadata_selection <- metadata_selection %>%
  arrange(Condition, Timepoint) %>%
  as.data.frame() %>%
  column_to_rownames("ExtBarCode")


# Color palette -----------------------------------------------------------

colors_heatmap <- list("Condition" = c("Infection" = "#7876B1FF", "Control" = "#7876B199"),
                       "Timepoint" = c("Day -3" = "#E18727FF", "Day -2" = "#0072B5FF", "Day -1" = "#BC3C29FF"))

color_gradient <- viridis(1000)


# Heatmap -----------------------------------------------------------------

expr_heatmap <- t(expr_scaled) # Genes are always represented in rows in a heatmap (samples in columns)!
dim(expr_heatmap) # 1500 genes in 580 samples (30/30 patients?)

samples_order <- match(rownames(metadata_selection), colnames(expr_heatmap))

heatmap_final <- pheatmap(
  expr_heatmap[, samples_order],
  color = color_gradient,
  cluster_cols = F,
  annotation_colors = colors_heatmap,
  annotation_col = metadata_selection,
  cutree_rows = 2,
  show_rownames = F,
  show_colnames = F,
  gaps_col = (which(metadata_selection$Condition == "Control")[1]-1), # Cuts the heatmap a position before the start the first Control.
  treeheight_row = 0 # To hide the dendrogram
)

heatmap_final

# Export in SVG format with Cairo

library(Cairo)
Cairo::CairoSVG(file = "./res/illumina_heatmap_1500.svg", width = 14, height = 7)
heatmap_final
dev.off()

# PDF export

pdf(file = "./res/illumina_heatmap_1500.pdf", width = 14, height = 7)
heatmap_final
dev.off()
