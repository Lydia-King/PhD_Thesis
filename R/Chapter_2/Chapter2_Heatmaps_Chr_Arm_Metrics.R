# Chapter 2: Heatmaps of Chromosome Arm metrics
## Load up libraries
library(ComplexHeatmap)
library(tidyverse)
library(circlize)

## Load up data
CNA_Arm_Score_Metrics <-  read.delim(
    "../../data/Processed_Data/Chr_Arm_CNA_Score.txt",
    sep = "\t",
    na.strings = c("", " ", "NA")
)

CNA_Arm_Burden_Metrics <- read.delim(
    "../../data/Processed_Data/Chr_Arm_CNA_Burden.txt",
    sep = "\t",
    na.strings = c("", " ", "NA")
)

## Setup
CNA_Score_metrics_grep <-
    c("CNA_Score",
      "CNA_Amp",
      "CNA_Del",
      "CNA_Difference",
      "CNA_Per_Amp",
      "CNA_Per_Del")

Names_Score <-
    c(
        "CNA \nScore",
        "Amp\nScore",
        "Del\nScore",
        "Difference \nScore",
        "% Amplification \nScore",
        "% Deletion \nScore"
    )

CNA_Burden_metrics_grep <-
    c("CNA_Burden",
      "CNA_Amp",
      "CNA_Del",
      "CNA_Difference",
      "CNA_Per_Amp",
      "CNA_Per_Del")

Names_Burden <-
    c(
        "CNA \nBurden",
        "Amp\nBurden",
        "Del\nBurden",
        "Difference \nBurden",
        "% Amplification \nBurden",
        "% Deletion \nBurden"
    )

CNA_Score_Labels <-
    c(
        "(B) Heatmap of CNA Score on each Chromosome Arm",
        "(B) Heatmap of CNA Amp Score on each Chromosome Arm",
        "(B) Heatmap of CNA Del Score on each Chromosome Arm",
        "(B) Heatmap of CNA Difference Score on each Chromosome Arm",
        "(B) Heatmap of Percentage CNA Amplification Score on each Chromosome Arm",
        "(B) Heatmap of Percentage CNA Deletion Score on each Chromosome Arm"
    )

CNA_Score_Labels_2 <-
    c(
        "(A) Heatmap of CNA Score on each Chromosome Arm",
        "(A) Heatmap of CNA Amp Score on each Chromosome Arm",
        "(A) Heatmap of CNA Del Score on each Chromosome Arm",
        "(A) Heatmap of CNA Difference Score on each Chromosome Arm",
        "(A) Heatmap of Percentage CNA Amplification Score on each Chromosome Arm",
        "(A) Heatmap of Percentage CNA Deletion Score on each Chromosome Arm"
    )


CNA_Burden_Labels <-
    c(
        "(B) Heatmap of CNA Burden on each Chromosome Arm",
        "(B) Heatmap of CNA Amplification Burden on each Chromosome Arm",
        "(B) Heatmap of CNA Deletion Burden on each Chromosome Arm",
        "(B) Heatmap of CNA Difference Burden on each Chromosome Arm",
        "(B) Heatmap of Percentage CNA Amplification Burden on each Chromosome Arm",
        "(B) Heatmap of Percentage CNA Deletion Burden on each Chromosome Arm"
    )

CNA_Burden_Labels_2 <-
    c(
        "(A) Heatmap of CNA Burden on each Chromosome Arm",
        "(A) Heatmap of CNA Amplification Burden on each Chromosome Arm",
        "(A) Heatmap of CNA Deletion Burden on each Chromosome Arm",
        "(A) Heatmap of CNA Difference Burden on each Chromosome Arm",
        "(A) Heatmap of Percentage CNA Amplification Burden on each Chromosome Arm",
        "(A) Heatmap of Percentage CNA Deletion Burden on each Chromosome Arm"
    )

setwd("../../figures/Chapter_2/")

## CNA Score
CNA_Score_Heatmap <-
    CNA_Arm_Score_Metrics %>% dplyr::select(grep("All", names(CNA_Arm_Score_Metrics)))

for (i in 1:length(CNA_Score_metrics_grep)) {
    CNA_Score_Heatmap_1 <-
        CNA_Score_Heatmap %>% dplyr::select(grep(CNA_Score_metrics_grep[i], names(CNA_Score_Heatmap)))
    colnames(CNA_Score_Heatmap_1) <-
        gsub('.*\\_', 'Chr ', colnames(CNA_Score_Heatmap_1))
    
    if (i == 1) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1)), c("white", "purple"))
    } else if (i == 2) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1)), c("white", "red"))
    } else if (i == 3) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1)), c("white", "blue"))
    } else if (i == 4) {
        col_fun <-
            colorRamp2(c(
                min(CNA_Score_Heatmap_1),
                0,
                max(CNA_Score_Heatmap_1)
            ), c("blue", "white", "red"))
    } else if (i == 5) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1)), c("white", "red"))
    } else {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1)), c("white", "blue"))
    }
    
    ## Produce Heatmaps - Add in Common fragile site markers
    filename <-
        paste(CNA_Score_metrics_grep[i], "_Score_Heatmap_All.png", sep = "")
    
    ht11 = Heatmap(
        as.matrix(CNA_Score_Heatmap_1),
        column_title = CNA_Score_Labels[i],
        row_title = "Patients",
        show_column_names = T,
        show_row_names = F,
        name = Names_Score[i],
        column_title_side = "top",
        column_order = colnames(CNA_Score_Heatmap_1),
        row_gap = unit(5, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 16),
        column_names_rot = 45,
        col = col_fun
    )
    
    png(
        filename,
        width = 27,
        height =  16,
        units = "cm",
        res = 400
    )
    draw(ht11)
    dev.off()
    
}

## CNA Burden
CNA_Burden_Heatmap <-
    CNA_Arm_Burden_Metrics %>% dplyr::select(grep("All", names(CNA_Arm_Burden_Metrics)))

for (i in 1:length(CNA_Burden_metrics_grep)) {
    CNA_Burden_Heatmap_1 <-
        CNA_Burden_Heatmap %>% dplyr::select(grep(CNA_Burden_metrics_grep[i], names(CNA_Burden_Heatmap)))
    colnames(CNA_Burden_Heatmap_1) <-
        gsub('.*\\_', 'Chr ', colnames(CNA_Burden_Heatmap_1))
    
    if (i == 1) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1)), c("white", "purple"))
    } else if (i == 2) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1)), c("white", "red"))
    } else if (i == 3) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1)), c("white", "blue"))
    } else if (i == 4) {
        col_fun <-
            colorRamp2(c(
                min(CNA_Burden_Heatmap_1),
                0,
                max(CNA_Burden_Heatmap_1)
            ), c("blue", "white", "red"))
    } else if (i == 5) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1)), c("white", "red"))
    } else {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1)), c("white", "blue"))
    }
    
    ## Produce Heatmaps - Add in Common fragile site markers
    filename <-
        paste(CNA_Burden_metrics_grep[i], "_Burden_Heatmap_All.png", sep = "")
    ht11 = Heatmap(
        as.matrix(CNA_Burden_Heatmap_1),
        column_title = CNA_Burden_Labels[i],
        show_column_names = T,
        show_row_names = F,
        name = Names_Burden[i],
        row_title = "Patients",
        column_title_side = "top",
        column_order = colnames(CNA_Burden_Heatmap_1),
        row_gap = unit(5, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 16),
        column_names_rot = 45,
        col = col_fun
    )
    
    png(
        filename,
        width = 27,
        height =  16,
        units = "cm",
        res = 400
    )
    draw(ht11)
    dev.off()
    
}


## CCA
## CNA Score
CNA_Score_Heatmap <-
    CNA_Arm_Score_Metrics %>% dplyr::select(grep("CCA", names(CNA_Arm_Score_Metrics)))

for (i in 1:6) {
    CNA_Score_Heatmap_1 <-
        CNA_Score_Heatmap %>% dplyr::select(grep(CNA_Score_metrics_grep[i], names(CNA_Score_Heatmap)))
    colnames(CNA_Score_Heatmap_1) <-
        gsub('.*\\_', 'Chr ', colnames(CNA_Score_Heatmap_1))
    
    if (i == 1) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1, na.rm = T)), c("white", "purple"))
    } else if (i == 2) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1, na.rm = T)), c("white", "red"))
    } else if (i == 3) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1, na.rm = T)), c("white", "blue"))
    } else if (i == 4) {
        col_fun <-
            colorRamp2(c(
                min(CNA_Score_Heatmap_1, na.rm = T),
                0,
                max(CNA_Score_Heatmap_1, na.rm = T)
            ), c("blue", "white", "red"))
    } else if (i == 5) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1, na.rm = T)), c("white", "red"))
    } else {
        col_fun <-
            colorRamp2(c(0, max(CNA_Score_Heatmap_1, na.rm = T)), c("white", "blue"))
    }
    
    ## Produce Heatmaps - Add in Common fragile site markers
    filename <-
        paste(CNA_Score_metrics_grep[i], "_Score_Heatmap_CCA.png", sep = "")
    ht11 = Heatmap(
        as.matrix(CNA_Score_Heatmap_1),
        column_title = CNA_Score_Labels_2[i],
        row_title = "Patients",
        show_column_names = T,
        show_row_names = F,
        name = Names_Score[i],
        column_title_side = "top",
        column_order = colnames(CNA_Score_Heatmap_1),
        row_gap = unit(5, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 16),
        column_names_rot = 45,
        col = col_fun,
        na_col = "black"
    )
    
    png(
        filename,
        width = 27,
        height =  16,
        units = "cm",
        res = 400
    )
    draw(ht11)
    dev.off()
    
}

## CNA Burden
CNA_Burden_Heatmap <-
    CNA_Arm_Burden_Metrics %>% dplyr::select(grep("CCA", names(CNA_Arm_Burden_Metrics)))


for (i in 1:length(CNA_Burden_metrics_grep)) {
    CNA_Burden_Heatmap_1 <-
        CNA_Burden_Heatmap %>% dplyr::select(grep(CNA_Burden_metrics_grep[i], names(CNA_Burden_Heatmap)))
    colnames(CNA_Burden_Heatmap_1) <-
        gsub('.*\\_', 'Chr ', colnames(CNA_Burden_Heatmap_1))
    
    if (i == 1) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1, na.rm = T)), c("white", "purple"))
    } else if (i == 2) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1, na.rm = T)), c("white", "red"))
    } else if (i == 3) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1, na.rm = T)), c("white", "blue"))
    } else if (i == 4) {
        col_fun <-
            colorRamp2(c(
                min(CNA_Burden_Heatmap_1, na.rm = T),
                0,
                max(CNA_Burden_Heatmap_1, na.rm = T)
            ), c("blue", "white", "red"))
    } else if (i == 5) {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1, na.rm = T)), c("white", "red"))
    } else {
        col_fun <-
            colorRamp2(c(0, max(CNA_Burden_Heatmap_1, na.rm = T)), c("white", "blue"))
    }
    
    ## Produce Heatmaps - Add in Common fragile site markers
    filename <-
        paste(CNA_Burden_metrics_grep[i], "_Burden_Heatmap_CCA.png", sep = "")
    ht11 = Heatmap(
        as.matrix(CNA_Burden_Heatmap_1),
        column_title = CNA_Burden_Labels_2[i],
        show_column_names = T,
        show_row_names = F,
        name = Names_Burden[i],
        row_title = "Patients",
        column_title_side = "top",
        column_order = colnames(CNA_Burden_Heatmap_1),
        row_gap = unit(5, "mm"),
        column_names_gp = grid::gpar(fontsize = 10),
        row_names_gp = grid::gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 16),
        column_names_rot = 45,
        col = col_fun,
        na_col = "black"
    )
    
    png(
        filename,
        width = 27,
        height =  16,
        units = "cm",
        res = 400
    )
    draw(ht11)
    dev.off()
    
}