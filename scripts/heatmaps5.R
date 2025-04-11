# April 2025
#  Draw heatmap from Eulalie's fluidigm results

rm(list = ls(all.names = TRUE))

library(pheatmap)
library(readr)
library(tidyverse)
library(dplyr)
library(RColorBrewer)
set.seed(1)

dataWD <-
  "D:/smyl/Labo/Mouse/KO_NCL/Transcriptomics/fluidigm/Fluidigm_results/Analyse_Fab/data"
resultsWD <-
  "D:/smyl/Labo/Mouse/KO_NCL/Transcriptomics/fluidigm/Fluidigm_results/Analyse_Fab/results"

setwd(dataWD)

loaded_data <- read_csv2("onGoingWork_Tableau_recap_fluidigm_noRef_GenoStade.csv")

aver_data_tibble <- loaded_data %>%                           # take the mean...
                    group_by(Geno_Stade) %>%                  #... per stade of the 3 values for each gene tested 
                    summarise_all(mean)  %>%                  # or summarise_all(list(mean=~mean(.))) 
                    column_to_rownames(.,"Geno_Stade") %>%    # otherwise matrix incorporates names in numerical data
                    as.matrix(.) %>%                          # ugly way to transpose the tibble
                    t(.) %>%                                  # idem
                    as_tibble (., rownames="Genes")           # idem

# write.csv2(loaded_dataset, file = "D:/smyl/Labo/Mouse/KO_NCL/Transcriptomics/fluidigm/Fluidigm_results/Analyse_Fab/results/loaded_dataset.csv")
# write.csv2(aver_dataset, file = "D:/smyl/Labo/Mouse/KO_NCL/Transcriptomics/fluidigm/Fluidigm_results/Analyse_Fab/results/averaged_dataset.csv")

# re-order cols using predifined list of column names
namesOfCol <- c("Genes", "S1_E18.5", "M1_E18.5", "S1_P4", "M1_P4", "S1_P11", "M1_P11",  "S1_P14","M1_P14", "S1_Adult","M1_Adult")
dataset_noFunc <- aver_data_tibble[,namesOfCol]

fold_changes <- transmute (dataset_noFunc, 
                          Genes = Genes,
                          E18.5 = M1_E18.5 - S1_E18.5, 
                          P4 = M1_P4 - S1_P4, 
                          P11 = M1_P11 - S1_P11,
                          P14 = M1_P14 - S1_P14,
                          Adult = M1_Adult - S1_Adult)


##############################" choice of functions to normalise rows
#cal_z_score <- function(x){
#(x - mean(x)) / sd(x)
#}
 
# or
#cal_z_score <- function(x){
#log2(x)
#}

# or
cal_z_score <- function(x){
x
}

# matrix_Func <- as.matrix(apply(dataset_noFunc, 1, cal_z_score)) # t(matrix) to return the transposed of matrix

matrix_Func <-  dataset_noFunc %>%
                column_to_rownames(var = "Genes") %>%
                as.matrix() %>%
                apply(1, cal_z_score) %>%
                t()

matrix_Func_FC <-  fold_changes %>%
  column_to_rownames(var = "Genes") %>%
  as.matrix() %>%
  apply(1, cal_z_score) %>%
  t()

paletteLength <- 10
breaksList <- c(seq(min(matrix_Func_FC), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(matrix_Func_FC)/paletteLength, max(matrix_Func_FC), length.out=floor(paletteLength/2)))
#myColor <- colorRampPalette(c("blue","green","w","yellow", "orange","red"))(paletteLength)
myColor <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

       
# pheatmap(matrix_Func,
#          scale ="none",   # possible values are "row", "column" and "none"
#          cluster_cols=0, cutree_cols=0,
#          cluster_rows=1, cutree_rows=3 ,
#          treeheight_row=0,
#          treeheight_col=0,
#          fontsize=8,
#          color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
#          breaks = breaksList,
#          main = "Hemizygous KO / Wild Type")

pheatmap(matrix_Func_FC,
         scale ="none",   # possible values are "row", "column" and "none"
         cluster_cols=0, cutree_cols=1,
         cluster_rows=1, cutree_rows=1 ,
         treeheight_row=0,
         treeheight_col=0,
         color = myColor,
         breaks = breaksList,
         fontsize=8,
         main = "Difference Hemizygous KO vs Wild Type")

minMat <- min(matrix_Func_FC)
maxMat <- max(matrix_Func_FC)

