# 📦 Load packages
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(here)
library(heatmaply)

set.seed(1) 

# 📂 Load data
loaded_data <- read_csv2(here("data", "onGoingWork_Tableau_recap_fluidigm_noRef_GenoStade.csv"))


# 🧮 Build a matrix that keeps all 3 mice per geno_stade
geno_labels <- make.unique(as.character(loaded_data$Geno_Stade))

matrix_expr_3_mice <- loaded_data[, !(names(loaded_data) %in% "Geno_Stade")]
rownames(matrix_expr_3_mice) <- geno_labels
matrix_expr_3_mice <- as.matrix(matrix_expr_3_mice)
matrix_expr_3_mice <- t(matrix_expr_3_mice)
colnames(matrix_expr_3_mice) <- geno_labels
matrix_expr_3_mice <- log2(matrix_expr_3_mice)

# 🧮 Average replicates by condition
averaged_data <- loaded_data %>%
  group_by(Geno_Stade) %>% # so to perform operations within each group
  summarise(across(everything(), mean), .groups = "drop") %>% # Applies the mean function to all columns within each group
  column_to_rownames("Geno_Stade") %>%
  as.matrix() %>%
  t() %>%
  as_tibble(rownames = "Genes")

# ✨ Reorder columns according to desired stage/genotype order
col_order <- c("Genes", "S1_E18.5", "M1_E18.5", "S1_P4", "M1_P4", 
               "S1_P11", "M1_P11", "S1_P14", "M1_P14", "S1_Adult", "M1_Adult")

ordered_data <- averaged_data[, col_order]

# 📏 Compute fold changes (M1/S1)
fold_changes <- ordered_data %>%
  transmute(
    Genes,
    E18.5 = M1_E18.5 / S1_E18.5,
    P4    = M1_P4    / S1_P4,
    P11   = M1_P11   / S1_P11,
    P14   = M1_P14   / S1_P14,
    Adult = M1_Adult / S1_Adult
  )


# 🧬 Prepare heatmap matrix
matrix_expr <- ordered_data %>%
  column_to_rownames("Genes") %>%
  as.matrix()

min_value_expr <- min(matrix_expr, na.rm = TRUE)  
max_value_expr <- max(matrix_expr, na.rm = TRUE)  

# 🧬 Prepare FC matrix (fold change, log2)
matrix_fc_log2 <- fold_changes %>%
  column_to_rownames("Genes") %>%
  as.matrix() %>%
  log2()

min_value_fc_log2 <- min(matrix_fc_log2, na.rm = TRUE)  # Minimum (ignores NA)
max_value_fc_log2 <- max(matrix_fc_log2, na.rm = TRUE)  # Maximum (ignores NA)

# 🎨 Heatmap color scale
breaks_list_expr <- seq(0, 8, by = 0.5)
breaks_list_fc <- seq(-3, 3, by = 0.5)
color_expr <- colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(length(breaks_list_expr))
color_fc   <- colorRampPalette(rev(brewer.pal(7, "RdBu")))(length(breaks_list_fc))

# 🔥 Expression heatmap
pheatmap(
  matrix_expr,
  scale = "none",
  cluster_cols = TRUE,
  cutree_cols = 1,
  cluster_rows = FALSE,
  cutree_rows = 3,
  treeheight_row = 40,
  treeheight_col = 40,
  fontsize = 8,
  color = color_expr,
  breaks = breaks_list_expr,
  main = "Expression heatmap"
)


# 🔥 Interactive plot expression heatmap

heatmaply(
  matrix_expr,
  scale = "none",
  colors = color_expr,
  k_col = 1,   # column clustering
  k_row = 3,   # row clustering
  dendrogram = "none",  # column, row or none
  main = "Expression heatmap"
)

# 🔥 Plot fold-change heatmap
PFClog <- pheatmap(
  matrix_fc_log2,
  scale = "none",
  cluster_cols = FALSE,
  cutree_cols = 1,
  cluster_rows = TRUE,
  cutree_rows = 3,
  treeheight_row = 20,
  treeheight_col = 0,
  fontsize = 8,
  color = color_fc,
  breaks = breaks_list_fc,
  main = "Log2(FC)  KO / WT"
)

print(PFClog)

# 🔥 Interactive plot fold-change heatmap

heatmaply(
  matrix_fc_log2,
  scale = "none",
  colors = color_fc,
  k_col = 0,   # column clustering
  k_row = 3,   # row clustering
  dendrogram = "row",  # column, row or none
  main = "log2(fold-change) heatmap"
)

# 😬 Interactive plot data

heatmaply(
  matrix_expr_3_mice,
  scale = "none",
  # colors = color_fc,
  # k_col = 0,   # column clustering
  k_row = 4,   # row clustering
  dendrogram = "row",  # both, column, row or none
  main = "log2(fold-change) heatmap"
)