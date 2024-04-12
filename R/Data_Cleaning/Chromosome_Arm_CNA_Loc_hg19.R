# Chapter 2: Get Chromosome Arm Location for Genes in CNA file

## Load up libraries
library(tidyverse)
library(GeneBreak)
library(operator.tools)
library(biomaRt)

### CNA METABRIC 2021
CNA <-
    read.delim(
        "../../data/METABRIC_2021/data_CNA.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

### Annotation Data
data(ens.gene.ann.hg19)

CNA_Loc <- CNA %>% select(1:2)
CNA_Loc <- merge(CNA_Loc, ens.gene.ann.hg19 , by.x="Hugo_Symbol", by.y = "Gene") %>% 
    distinct(Hugo_Symbol, .keep_all = T)

### What genes are missing
Missing_Genes_1 <-
    CNA$Hugo_Symbol[which(CNA$Hugo_Symbol %!in% ens.gene.ann.hg19$Gene)]

#### Biomart 
ensembl <- useMart("ensembl")
ensembl <- useDataset("hsapiens_gene_ensembl", mart = ensembl)

Biomart_Results <-
    getBM(
        attributes = c(
            "hgnc_symbol",
            "chromosome_name",
            "start_position",
            "end_position",
            "band"
        ),
        filters = c("hgnc_symbol"),
        values = c(Missing_Genes_1),
        mart = ensembl
    )

### Combine 
Total_1 <- CNA[CNA$Hugo_Symbol %in% Missing_Genes_1,c(1:2)]
Total_1 <- merge(Total_1, Biomart_Results, by.x="Hugo_Symbol", by.y = "hgnc_symbol")
colnames(Total_1) <- c("Hugo_Symbol", "Entrez_Gene_Id", "Chromosome", "Start", "End", "band")

Total_CNA_2 <- rbind.data.frame(Total_1, CNA_Loc[,-c(3,8)])

### Tidy up and combine with CNA data 
Total_CNA_2$Chromosome <- gsub("chr", "\\1", Total_CNA_2$Chromosome)
Total_CNA_2$band <- gsub("(p).*", "\\1", Total_CNA_2$band)
Total_CNA_2$band <- gsub("(q).*", "\\1", Total_CNA_2$band)

Total_CNA_2 <-
    Total_CNA_2 %>% mutate(Chr_Loc = paste(Chromosome, band, sep = ""))

Total_CNA_2 <- Total_CNA_2 %>% distinct(Hugo_Symbol, .keep_all = T)

### Combine with CNA Data
CNA_with_Loc <-
    merge(
        Total_CNA_2,
        CNA %>% dplyr::select(-Entrez_Gene_Id),
        by.x = "Hugo_Symbol",
        by.y = "Hugo_Symbol"
    )

write.table(CNA_with_Loc,
            file = "../../data/Processed_Data/data_CNA_Loc_hg19.txt",
            sep = "\t",
            row.names = FALSE)
