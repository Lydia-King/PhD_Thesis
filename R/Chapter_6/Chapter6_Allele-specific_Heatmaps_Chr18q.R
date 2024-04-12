# Chapter 6: Allele-specific Heatmaps - Chromosome 18q

## Load up Libraries 
library(tidyverse)

## Load up Data 
### CNA Location 
CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg19.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

CNA <- CNA_Loc %>% filter(Chr_Loc == "18q")

Data_Chr18q <- CNA %>% select(Hugo_Symbol, Start, End) %>% filter(!is.na(Start)) %>% 
    arrange(Start) 

### ASCAT Data 
ASCAT3_Data_1 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_1.segments.txt", sep="\t")
ASCAT3_Data_2 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_2.segments.txt", sep="\t")
ASCAT3_Data_3 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_3.segments.txt", sep="\t")
ASCAT3_Data_4 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_4.segments.txt", sep="\t")
ASCAT3_Data_5 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_5.segments.txt", sep="\t")
ASCAT3_Data_6 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_6.segments.txt", sep="\t")
ASCAT3_Data_7 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_7.segments.txt", sep="\t")
ASCAT3_Data_8 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_8.segments.txt", sep="\t")

Total_ASCAT18 <- rbind.data.frame(ASCAT3_Data_1, ASCAT3_Data_2, ASCAT3_Data_3, 
                                 ASCAT3_Data_4, ASCAT3_Data_5, ASCAT3_Data_6, 
                                 ASCAT3_Data_7, ASCAT3_Data_8)

### Chromosome 18q
Total_ASCAT18_Chr18 <- Total_ASCAT18 %>% filter(chr == "18")

# Initialize an empty list to store data frames
result_list <- list()

unique_samples <- unique(Total_ASCAT18$sample)

# Iterate over unique samples
for (pat in unique_samples) {
    test <- Total_ASCAT18_Chr18 %>% filter(sample == pat)
    
    # Initialize empty data frames for each sample
    sample_df <- data.frame("Sample" = character(0), 
                            "Gene" = character(0), 
                            "Major" = character(0), 
                            "Minor" = character(0))
    
    # Iterate over rows of Data_Chr18q
    for (i in 1:nrow(Data_Chr18q)) {
        row <- Data_Chr18q[i, ]
        indices <- which(as.numeric(row["Start"]) >= test$startpos & as.numeric(row["End"]) <= test$endpos)
        
        if (length(indices) > 0) {
            # If conditions are met, append the row to the sample_df
            sample_df <- rbind(sample_df, data.frame("Sample" = pat, "Gene" = row[["Hugo_Symbol"]], 
                                                     "Major" = test[indices, "nMajor"],
                                                     "Minor" = test[indices, "nMinor"]))
        } else {
            # If conditions are not met, append NA values
            sample_df <- rbind(sample_df, data.frame("Sample" = pat, "Gene" = row[["Hugo_Symbol"]], 
                                                     "Major" = "NA", "Minor" = "NA"))
        }
    }
    
    # Append the sample_df to the result_list
    result_list[[pat]] <- sample_df
}

# Combine the list of data frames into a single data frame
Chr18q_DF <- do.call(rbind, result_list)

Chr18q_DF <- Chr18q_DF %>% mutate(Major = as.numeric(Major)) %>% mutate(Major = case_when(Major == 0 ~ -1, 
                                                  Major == 1 ~ 0, 
                                                  Major == 2 ~ 1, 
                                                  Major == 3 ~ 1, 
                                                  Major > 3 ~ 2, 
                                                  Major == "NA" ~ NA), 
                                Minor = case_when(Minor == 0 ~ -1, 
                                                  Minor == 1 ~ 0, 
                                                  Minor == 2 ~ 1, 
                                                  Minor == 3 ~ 1, 
                                                  Minor > 3 ~ 2, 
                                                  Minor == "NA" ~ NA))

## Add Nodes to dataset
data2 <-  readRDS(file = "../../data/Chapter_3/Chap3_List_of_PA_Data.rds")
data2 <- data2[[2]]
data2$PATIENT_ID <- gsub("\\-", ".", data2$PATIENT_ID)

Chr18q_DF <- merge(Chr18q_DF, data2, by.x = "Sample", by.y ="PATIENT_ID")
Chr18q_DF$Node <- paste("Node", Chr18q_DF$Node, sep = " ")

Chr18q_DF_HM1_Major <- pivot_wider(
    Chr18q_DF[,c(1,2,3,5)],
    id_cols = c("Sample", "Node"),
    names_from = Gene,
    values_from = Major,
    names_prefix = ""
)

vec <- c("Sample", "Node", Data_Chr18q$Hugo_Symbol)
Chr18q_DF_HM1_Major <- Chr18q_DF_HM1_Major[vec]

Chr18q_DF_HM1_Minor <- pivot_wider(
    Chr18q_DF[,c(1,2,4,5)],
    id_cols = c("Sample", "Node"),
    names_from = Gene,
    values_from = Minor,
    names_prefix = ""
)

vec <- c("Sample", "Node", Data_Chr18q$Hugo_Symbol)
Chr18q_DF_HM1_Minor <- Chr18q_DF_HM1_Minor[vec]

# 3 heatmaps, 1 major, 1 minor and total (mutate to match)
library(circlize)
library(ComplexHeatmap)

## Major
col_fun = colorRamp2(c(min(Chr18q_DF_HM1_Major[,-c(1, 2)], na.rm = T), 0, max(Chr18q_DF_HM1_Major[,-c(1, 2)], na.rm = T)), c("blue", "white", "red"))

ht11 = Heatmap(as.matrix(Chr18q_DF_HM1_Major[,-c(1, 2)]), 
               column_order = colnames(Chr18q_DF_HM1_Major[,-c(1, 2)]), 
               na_col = "black", cluster_rows = T, 
               col = col_fun, column_title = paste("Heatmap of CNAs on the Major Allele of Chromosome 18q"), 
               name ="CNA",  column_names_gp = grid::gpar(fontsize = 6), row_title_gp = grid::gpar(fontsize = 13), 
               column_title_gp = gpar(fontsize = 16),  column_title_side = "top", column_names_rot = 45, row_split = Chr18q_DF_HM1_Major[,c("Node")], 
               show_column_names = F, show_row_names = F, 
               cluster_row_slices = FALSE)

png("../../figures/Chapter_6/Heatmap_Chr18q_Genes_Major.png", width = 23, height =  13, units = "cm", res = 400)
draw(ht11, padding = unit(c(8,2,2,2), "mm"))
dev.off()

## Minor
col_fun = colorRamp2(c(min(Chr18q_DF_HM1_Minor[,-c(1, 2)], na.rm = T), 0, max(Chr18q_DF_HM1_Minor[,-c(1, 2)], na.rm = T)), c("blue", "white", "red"))

ht11 = Heatmap(as.matrix(Chr18q_DF_HM1_Minor[,-c(1, 2)]), 
               column_order = colnames(Chr18q_DF_HM1_Minor[,-c(1, 2)]), 
               na_col = "black", cluster_rows = T, 
               col = col_fun, column_title = paste("Heatmap of CNAs on the Minor Allele of Chromosome 18q"), 
               name ="CNA",  column_names_gp = grid::gpar(fontsize = 6), row_names_gp = grid::gpar(fontsize = 14), row_title_gp = grid::gpar(fontsize = 13),
               column_title_gp = gpar(fontsize = 16),  column_title_side = "top", column_names_rot = 45, row_split = Chr18q_DF_HM1_Minor[,c("Node")], show_column_names = F, show_row_names = F, 
               cluster_row_slices = FALSE)

png("../../figures/Chapter_6/Heatmap_Chr18q_Genes_Minor.png", width = 23, height =  13, units = "cm", res = 400)
draw(ht11, padding = unit(c(8,2,2,2), "mm"))
dev.off()

## Total 
Chr18q_DF <- do.call(rbind, result_list)

Chr18q_DF$Total <- as.numeric(Chr18q_DF$Major) + as.numeric(Chr18q_DF$Minor)

Chr18q_DF <- Chr18q_DF %>% mutate(Total = case_when(Total == 0 ~ -2, 
                                                  Total == 1 ~ -1, 
                                                  Total == 2 ~ 0, 
                                                  Total == 3 ~ 1,
                                                  Total == 4 ~ 1, 
                                                  Total == 5 ~ 1,
                                                  Total > 5 ~ 2, 
                                                  Total == "NA" ~ NA))

## Add Nodes to dataset
Chr18q_DF <- merge(Chr18q_DF, data2, by.x = "Sample", by.y ="PATIENT_ID")

Chr18q_DF_HM1_Total <- pivot_wider(
    Chr18q_DF[,c(1,2,5,6)],
    id_cols = c("Sample", "Node"),
    names_from = Gene,
    values_from = Total,
    names_prefix = ""
)

vec <- c("Sample", "Node", Data_Chr18q$Hugo_Symbol)
Chr18q_DF_HM1_Total <- Chr18q_DF_HM1_Total[vec]

col_fun = colorRamp2(c(min(Chr18q_DF_HM1_Total[,-c(1, 2)], na.rm = T), 0, max(Chr18q_DF_HM1_Total[,-c(1, 2)], na.rm = T)), c("blue", "white", "red"))

ht11 = Heatmap(as.matrix(Chr18q_DF_HM1_Total[,-c(1,2)]), 
               column_order = colnames(Chr18q_DF_HM1_Total[,-c(1,2)]), 
               na_col = "black", cluster_rows = T, 
               col = col_fun, column_title = paste("Heatmap of CNA on Both Alleles of Chromosome 18q"), row_title_gp = grid::gpar(fontsize = 13),
               name ="CNA",  column_names_gp = grid::gpar(fontsize = 6), row_names_gp = grid::gpar(fontsize = 14), 
               column_title_gp = gpar(fontsize = 16),  column_title_side = "top", column_names_rot = 45, show_column_names = F, show_row_names = F, row_split = Chr18q_DF_HM1_Minor[,c("Node")], 
               cluster_row_slices = FALSE)

png("../../figures/Chapter_6/Heatmap_Chr18q_Genes_Both_Alleles.png", width = 23, height =  13, units = "cm", res = 400)
draw(ht11, padding = unit(c(8,2,2,2), "mm"))
dev.off()