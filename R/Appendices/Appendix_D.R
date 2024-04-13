# Comparative Study - IntClust Genes Missing 

## Libraries 
library(tidyverse)
library(readxl)
library(xtable)
library(geneSynonym)

## Load Up CNA and GE Data (2023)
CNA_2023 <- read.delim("../../data/METABRIC_2023/data_cna.txt", sep="\t") 
CNA_2020 <- read.delim("../../data/METABRIC_2019/data_CNA.txt", sep="\t") 
GE_2023 <- read.delim("../../data/METABRIC_2023/data_mrna_illumina_microarray.txt", sep="\t")

## Combine Gene Expression and CNA data (based on Hugo_Symbol)
CNA_GE_Data <- merge(CNA_2023[,c(1,2)], GE_2023[,c(1,2)], by = "Hugo_Symbol")

## Load up IntClust Genes
IC <- read_xls("../../data/IntClust/table_S30.xls")
IC_Sort <- IC[order(IC$All_ANOVA.pval.adj),]

Top_1000_Genes <- IC_Sort[1:1165, ]
Top_1000_Genes_NoDup <- Top_1000_Genes %>% distinct(Gene, .keep_all = T)
Top_1000_Genes_NoDup$Gene <- gsub("ORF", "orf", Top_1000_Genes_NoDup$Gene)

## Analyse overlap and missingness 
Missing_Genes_1 <- Top_1000_Genes_NoDup$Gene[which(Top_1000_Genes_NoDup$Gene %!in% CNA_GE_Data$Hugo_Symbol)]

# Missing_Genes_1 <- Top_1000_Genes$Gene[which(Top_1000_Genes$Gene %!in% GE_Loc_Common$Hugo_Symbol_CNA)]
Old_Gene_Name <- c()
New_Gene_Name <- c()

for(i in 1:length(Missing_Genes_1)){
    
    syn <- geneSynonym(c(Missing_Genes_1[i]), tax = 9606)
    
    Old_Gene_Name <- c(Old_Gene_Name, Missing_Genes_1[i])
    p <- paste(syn[[1]][[1]][which(syn[[1]][[1]] %in% CNA_GE_Data$Hugo_Symbol)], sep="|")
    New_Gene_Name <- c(New_Gene_Name, ifelse(is_empty(p), NA, p))
}

Gene_Swap <- data.frame(Old_Gene_Name, New_Gene_Name)

for(j in 1:length(Missing_Genes_1)){
    if(!is.na(New_Gene_Name)[j]){
        Top_1000_Genes_NoDup[Top_1000_Genes_NoDup$Gene == Old_Gene_Name[j], 3] <- New_Gene_Name[j] 
    }
}

Missing_Genes_2 <- Top_1000_Genes_NoDup$Gene[which(Top_1000_Genes_NoDup$Gene %!in% CNA_GE_Data$Hugo_Symbol)]

print(xtable(data.frame("Gene" = Missing_Genes_2), align = "|c|c|", caption = c("78 IntClust genes not present in CNA and/or gene expression data.")), include.rownames = F)


noquote(paste(noquote(Missing_Genes_2), ", ", sep=""))

