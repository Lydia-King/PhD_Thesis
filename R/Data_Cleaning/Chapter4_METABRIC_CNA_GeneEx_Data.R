# Chapter 4: Setup CNA and Gene Expression Data

## Load up Libraries
library(tidyverse)
library(operator.tools)

## Load up Data (CNA and GE)
CNA_data <-
    read.delim("../../data/METABRIC_2023/data_cna.txt", sep = "\t")

GE <-
    read.delim("../../data/METABRIC_2023/data_mrna_illumina_microarray.txt",
               sep = "\t")

## Remove Duplicates (CNA)
CNA_data <- CNA_data %>%
    mutate(Hugo_Symbol = ifelse(Hugo_Symbol == "PALM2-AKAP2", "PALM2AKAP2", Hugo_Symbol))

CNA_data_NoDup <-
    CNA_data %>% distinct(Hugo_Symbol, Entrez_Gene_Id, .keep_all = T)

## Remove Duplicates (Gene Expression)
GE_NoDup <-
    GE %>% distinct(Hugo_Symbol, Entrez_Gene_Id, MB.0362, .keep_all = T)

## Which genes are in both datasets (common - using Entrez_Gene_Id)
Common_Genes_EGI <-
    merge(CNA_data_NoDup[, c(1, 2)], GE_NoDup, by = "Entrez_Gene_Id")

missing_CNA <-
    CNA_data_NoDup[, c(1, 2)][which(CNA_data_NoDup$Hugo_Symbol %!in% Common_Genes_EGI$Hugo_Symbol.x),]

missing_GE <-
    GE_NoDup[which(GE_NoDup$Hugo_Symbol %!in% Common_Genes_EGI$Hugo_Symbol.y),]

Common_Genes_EGI  <- Common_Genes_EGI %>%
    dplyr::rename("Hugo_Symbol_CNA" = Hugo_Symbol.x, "Hugo_Symbol_GE" = Hugo_Symbol.y) %>%
    dplyr::select(-Hugo_Symbol_GE,-Entrez_Gene_Id)

## Which missing genes are in both datasets (common - using Hugo_Symbol)
Common_Genes_HS <-
    merge(missing_CNA %>% select(-2), missing_GE[, -c(2)], by = "Hugo_Symbol") %>%
    dplyr::rename("Hugo_Symbol_CNA" = Hugo_Symbol)

Common_Genes <- rbind.data.frame(Common_Genes_HS, Common_Genes_EGI)
Common_Genes %>% distinct(Hugo_Symbol_CNA, .keep_all = T)

## Add Location Information (To Gene Expression)
CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg18.txt",
        sep = ",",
        na.strings = c("", " ", "NA")
    )

CNA_Loc <- CNA_Loc %>% mutate(Chr_Loc = paste(Chromosome, Arm, sep="")) %>% select("Hugo_Symbol", "Chr_Loc")

Common_Genes_Loc <-
    merge(CNA_Loc, Common_Genes, by.x = "Hugo_Symbol", by.y = "Hugo_Symbol_CNA")

## Write Out
write.table(Common_Genes_Loc,
            "../../data/Processed_Data/GE_Common_Genes.txt",
            sep = "\t")
