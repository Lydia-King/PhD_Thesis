# Chapter 2: Get Chromosome Arm Location for Genes in CNA file

## Load up libraries
library(tidyverse)
library(biomaRt)
library(geneSynonym)
library(operator.tools)
library(rCGH)

## Load up Data and Tidy
### CNA METABRIC 2021
CNA <-
    read.delim(
        "../../data/METABRIC_2021/data_CNA.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

### biomaRt Annotation
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
        values = c(CNA$Hugo_Symbol),
        mart = ensembl
    )

Biomart_Results_1 <-
    Biomart_Results[Biomart_Results$band != "", ] %>% distinct(., hgnc_symbol, .keep_all = TRUE)

### What genes are missing
Missing_Genes_1 <-
    CNA$Hugo_Symbol[which(CNA$Hugo_Symbol %!in% Biomart_Results_1$hgnc_symbol)]

### Get alternative names?
Old_Gene_Name <- c()
New_Gene_Name <- c()

for (i in 1:length(Missing_Genes_1)) {
    syn <- geneSynonym(c(Missing_Genes_1[i]), tax = 9606)
    
    Old_Gene_Name <- c(Old_Gene_Name, Missing_Genes_1[i])
    p <- syn[[1]][[1]][[1]]
    New_Gene_Name <- c(New_Gene_Name, ifelse(is_empty(p), NA, p))
}

Biomart_Results_2  <-
    getBM(
        attributes = c(
            "hgnc_symbol",
            "chromosome_name",
            "start_position",
            "end_position",
            "band"
        ),
        filters = c("hgnc_symbol"),
        values = c(New_Gene_Name),
        mart = ensembl
    )


Biomart_Results_2 <-
    merge(
        data.frame("CNA_Hugo_Symbol" = Old_Gene_Name, "Hugo_Symbol_Alt" = New_Gene_Name),
        Biomart_Results_2,
        by.x = "Hugo_Symbol_Alt",
        by.y = "hgnc_symbol"
    )

Biomart_Results_2 <-
    Biomart_Results_2[Biomart_Results_2$band != "", ] %>% 
    distinct(., Hugo_Symbol_Alt, .keep_all = TRUE)

### Combine Two result df
Biomart_Results_1 <-
    Biomart_Results_1 %>% mutate("CNA_Hugo_Symbol" = hgnc_symbol) %>%
    dplyr::select(hgnc_symbol,
                  CNA_Hugo_Symbol,
                  chromosome_name,
                  start_position,
                  end_position,
                  band)
Biomart_Results_2 <-
    Biomart_Results_2 %>% dplyr::rename("hgnc_symbol" = Hugo_Symbol_Alt)

Total_CNA <- rbind.data.frame(Biomart_Results_1, Biomart_Results_2)

### What genes are missing
Missing_Genes_2 <-
    CNA$Hugo_Symbol[which(CNA$Hugo_Symbol %!in% Total_CNA$CNA_Hugo_Symbol)]

### Annotate with Hg18
Hg18 <- read.delim("../../data/Annotation_Data/Hg18", sep = "\t")
Hg18 <- Hg18[Hg18$name2 %in% Missing_Genes_2, ]

Centromere_Loc <- hg18

Results_3 <-
    data.frame(
        "hgnc_symbol" = character(),
        "CNA_Hugo_Symbol" = character(),
        "chromosome_name" = character(),
        "start_position" = numeric(),
        "end_position" =  numeric(),
        "band" = character()
    )

for (i in 1:23) {
    Cent <- Centromere_Loc[i, ]
    
    if (i == 23) {
        dat <-
            Hg18 %>% filter(chrom == "chrX") %>% 
            mutate(band = ifelse(txStart < Cent$centromerStart, "p", "q"))
    } else {
        dat <-
            Hg18 %>% filter(chrom == paste("chr", i, sep = ""))  %>% 
            mutate(band = ifelse(txStart < Cent$centromerStart, "p", "q"))
    }
    
    Results_3 <-
        rbind.data.frame(
            Results_3,
            data.frame(
                "hgnc_symbol" = dat$name2,
                "CNA_Hugo_Symbol" = dat$name2,
                "chromosome_name" = dat$chrom,
                "start_position" = dat$txStart,
                "end_position" =  dat$txEnd,
                "band" = dat$band
            )
        )
    
}

Total_CNA_1 <-
    rbind.data.frame(
        Total_CNA,
        Results_3 %>% distinct(., CNA_Hugo_Symbol, chromosome_name, band, .keep_all = T)
    )

### What genes are missing
Missing_Genes_3 <-
    CNA$Hugo_Symbol[which(CNA$Hugo_Symbol %!in% Total_CNA_1$CNA_Hugo_Symbol)] # 38

### Manually Annotate
hgnc_symbol <-
    c(
        "SAXO5",
        "CIMIP2C",
        "RNF32-DT",
        "SPMIP7",
        "CRIPAK",
        "ERCC6",
        "WASHC2C",
        "LOC124903857",
        "FAM231B",
        "FAM231C",
        "FAM231D",
        "H2BW4P",
        "H2BC4",
        "H2BC5",
        "H2BC6",
        "H2BC7",
        "H2BC10",
        "SPMIP11",
        "MIR570HG",
        "RGS2-AS1",
        "MIR1254-1",
        "MIR1254-2",
        "MIR1273A",
        "MIR1273F",
        "MIR1273G",
        "MIR3656",
        "MIR4417",
        "MIR4419A",
        "MIR4419B",
        "MIR4459",
        "MIR4461",
        "MIR4532",
        "MIR4792",
        "MIR566",
        "OCLM",
        "PRAMEF16",
        "PRR21",
        "SPHAR"
    )
CNA_Hugo_Symbol <- Missing_Genes_3
chromosome_name <-
    c(
        "19",
        "2",
        "7",
        "7",
        "4",
        "10",
        "10",
        "1",
        "1",
        "1",
        "1",
        "X",
        "6",
        "6",
        "6",
        "6",
        "6",
        "12",
        "3",
        "1",
        "10",
        "10",
        "8",
        "1",
        "1",
        "11",
        "1",
        "1",
        "12",
        "5",
        "5",
        "20",
        "3",
        "3",
        "1",
        "1",
        "2",
        "1"
    )
start_position <- rep(NA, 38)
end_position <- rep(NA, 38)
band <-
    c(
        "p",
        "p",
        "q",
        "p",
        "p",
        "q",
        "q",
        "p",
        "p",
        "p",
        "q",
        "q",
        "p",
        "p",
        "p",
        "p",
        "p",
        "q",
        "q",
        "q",
        "p",
        "p",
        "q",
        "p",
        "p",
        "q",
        "p",
        "p",
        "q",
        "q",
        "q",
        "q",
        "p",
        "p",
        "q",
        "p",
        "q",
        "q"
    )

Total_CNA_2 <-
    rbind.data.frame(
        Total_CNA_1,
        data.frame(
            hgnc_symbol,
            CNA_Hugo_Symbol,
            chromosome_name,
            start_position,
            end_position,
            band
        )
    )

Missing_Genes_4 <-
    CNA$Hugo_Symbol[which(CNA$Hugo_Symbol %!in% Total_CNA_2$CNA_Hugo_Symbol)]

### Tidy up
Total_CNA_2$chromosome_name <- gsub("chr", "\\1", Total_CNA_2$chromosome_name)
Total_CNA_2$band <- gsub("(p).*", "\\1", Total_CNA_2$band)
Total_CNA_2$band <- gsub("(q).*", "\\1", Total_CNA_2$band)

Total_CNA_2 <-
    Total_CNA_2 %>% mutate(Chr_Loc = paste(chromosome_name, band, sep = ""))

### Combine with CNA Data
CNA_with_Loc <-
    merge(
        Total_CNA_2,
        CNA %>% dplyr::select(-Entrez_Gene_Id),
        by.x = "CNA_Hugo_Symbol",
        by.y = "Hugo_Symbol"
    )

write.table(CNA_with_Loc,
            file = "../../data/Processed_Data/data_CNA_Loc.txt",
            sep = "\t",
            row.names = FALSE)