# Libraries 
library(tidyverse)
library(rCGH)
library(GeneBreak)

# Data
setwd("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/ASCAT_Data/Formatted_Data/")
ASCAT_Data_3Step_NoNeut <- read.csv2("Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

ASCAT_Data_3Step_NoNeut_1 <- read.delim("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/PhD_Thesis/data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
data(ens.gene.ann.hg19)

CNA_Loc <- read.delim("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/data/METABRIC_Data_Common/CNA_Data_Loc_Inc.txt", sep="\t")
CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 3, band == "p") %>% select(1:7)

## 
result_list_1 <- list()
unique_samples <- unique(ASCAT_Data_3Step_NoNeut$Sample)
CHR <- c(1:22, "X")

Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "3") 
Total_ASCAT3_Chr3$Breakpoint <- as.numeric(Total_ASCAT3_Chr3$Breakpoint)

# Iterate over unique samples
sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(Total_ASCAT3_Chr3$Breakpoint >= CNA_Loc_3p[i,]$Start & Total_ASCAT3_Chr3$Breakpoint <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)
    
sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))


## Barplot?
ggplot(sample_df_2, aes(x=Gene, y=n, color = Allele)) + 
    geom_bar(stat = "identity")  + facet_grid(~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 3p by Allele") + ylab("Frequency") + 
    xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                          axis.title.x = element_text(vjust = +4)) + 
    theme(legend.position = "none")
    

ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_3p_Barplot_Allele.png", last_plot(), width = 10, height = 5)

as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                            Category == "Neutral/Deletion" ~ "Neut/Del",
                                            Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                            Category == "Deletion/Neutral" ~ "Del/Neut", 
                                            Category == "Deletion/Amplification" ~ "Del/Amp", 
                                            Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                 "Neut/Amp", "Neut/Del", 
                                                                                                                                                 "Amp/Del", "Del/Amp"))) %>% 

    ggplot(aes(x=Gene, y=n, color = interaction(Category, Allele, sep=':'))) + 
    geom_bar(stat = "identity") + facet_grid(Category~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 3p by Category and Allele") + 
    ylab("Frequency") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "none")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_3p_Barplot_AlleleCate.png", last_plot(), width = 10, height = 7)

## Node 
CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 3, band == "p") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "3") 

setwd("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/Draft_Corrections/Updated_Chapter3/data/")
data2 <- readRDS("Chap3_List_of_PA_Data.rds")
data2 <- data2[[1]]
data2$PATIENT_ID <- gsub("\\-", ".", data2$PATIENT_ID)
data2$Node <-paste("Node", data2$Node, sep = " ")

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, data2, by.x = "Sample", by.y = "PATIENT_ID")

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "Node" = character(0))


for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele, Node) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)

sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))

as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                     Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                     Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                     Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                     Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                     Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                                          "Neut/Amp", "Neut/Del", 
                                                                                                                                                                          "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=Gene, y=n, color = Allele)) + 
    geom_bar(stat = "identity",  width=0.5) + facet_grid(Category~ Node) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 3p by Category and Allele") + 
    ylab("Frequency") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "none")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_3p_Barplot_Node.png", last_plot(), width = 10, height = 7)

## Chromosome 18q 

CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 18, band == "q") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "18") 

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)

sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))


## Barplot?
ggplot(sample_df_2, aes(x=Gene, y=n, color = Allele, fill = Allele)) + 
    geom_bar(stat = "identity")  + facet_grid(~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 18q by Allele") + ylab("Frequency") + 
    xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                          axis.title.x = element_text(vjust = +4)) + 
    theme(legend.position = "none")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_18q_Barplot_Allele.png", last_plot(), width = 10, height = 5)

as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                     Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                     Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                     Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                     Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                     Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                                          "Neut/Amp", "Neut/Del", 
                                                                                                                                                                          "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=Gene, y=n, color = interaction(Category, Allele, sep=':'), fill =  interaction(Category, Allele, sep=':'))) + 
    geom_bar(stat = "identity") + facet_grid(Category~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 18q by Category and Allele") + 
    ylab("Frequency") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "none")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_18q_Barplot_AlleleCate.png", last_plot(), width = 10, height = 7)




## Chromosome 11p 
CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 11, band == "p") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "11") 

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)

sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))


## Barplot?
ggplot(sample_df_2, aes(x=Gene, y=n, color = Allele, fill = Allele)) + 
    geom_bar(stat = "identity")  + facet_grid(~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 11p by Allele") + ylab("Frequency") + 
    xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                          axis.title.x = element_text(vjust = +4)) + 
    theme(legend.position = "none")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_11p_Barplot_Allele.png", last_plot(), width = 10, height = 5)

as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                     Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                     Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                     Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                     Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                     Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                                          "Neut/Amp", "Neut/Del", 
                                                                                                                                                                          "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=Gene, y=n, color = interaction(Category, Allele, sep=':'), fill =  interaction(Category, Allele, sep=':'))) + 
    geom_bar(stat = "identity") + facet_grid(Category~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 11p by Category and Allele") + 
    ylab("Frequency") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "none")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_11p_Barplot_AlleleCate.png", last_plot(), width = 10, height = 7)


## PAM50 
## Chr 3p 
CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 3, band == "p") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "3") 

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "PAM50" = character(0))

Clinical <- read.delim("/home/lydia/Downloads/brca_metabric/data_clinical_patient_2.txt", comment.char = "#")
Clinical$PATIENT_ID <- gsub("\\-", ".", Clinical$PATIENT_ID)
Clinical <- Clinical %>% filter(CLAUDIN_SUBTYPE %!in% c("", " ", "NC"))

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, Clinical[,c("PATIENT_ID", "CLAUDIN_SUBTYPE")], by.x = "Sample", by.y = "PATIENT_ID")


for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele, CLAUDIN_SUBTYPE) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)

sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))

# PAM50
as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                     Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                     Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                     Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                     Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                     Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                                          "Neut/Amp", "Neut/Del", 
                                                                                                                                                                          "Amp/Del", "Del/Amp"))) %>% 
    
ggplot(aes(x=Gene, y=n, color = Category, fill = Category)) + 
    geom_bar(stat = "identity")  + facet_grid(CLAUDIN_SUBTYPE~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 3p by Allele and PAM50") + ylab("Frequency") + 
    xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                          axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "top")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_3p_Barplot_Allele_PAM50.png", last_plot(), width = 12, height = 13)






## Chr 18q
CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 18, band == "q") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "18") 

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "PAM50" = character(0))

Clinical <- read.delim("/home/lydia/Downloads/brca_metabric/data_clinical_patient_2.txt", comment.char = "#")
Clinical$PATIENT_ID <- gsub("\\-", ".", Clinical$PATIENT_ID)
Clinical <- Clinical %>% filter(CLAUDIN_SUBTYPE %!in% c("", " ", "NC"))

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, Clinical[,c("PATIENT_ID", "CLAUDIN_SUBTYPE")], by.x = "Sample", by.y = "PATIENT_ID")


for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele, CLAUDIN_SUBTYPE) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)

sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))

# PAM50
as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                     Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                     Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                     Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                     Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                     Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                                          "Neut/Amp", "Neut/Del", 
                                                                                                                                                                          "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=Gene, y=n, color = Category, fill = Category)) + 
    geom_bar(stat = "identity")  + facet_grid(CLAUDIN_SUBTYPE~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 18q by Allele and PAM50") + ylab("Frequency") + 
    xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                          axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "top")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_18q_Barplot_Allele_PAM50.png", last_plot(), width = 12, height = 13)



## Chr 11p
CNA_Loc_3p <- CNA_Loc %>% filter(Chromosome == 11, band == "p") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "11") 

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "PAM50" = character(0))

Clinical <- read.delim("/home/lydia/Downloads/brca_metabric/data_clinical_patient_2.txt", comment.char = "#")
Clinical$PATIENT_ID <- gsub("\\-", ".", Clinical$PATIENT_ID)
Clinical <- Clinical %>% filter(CLAUDIN_SUBTYPE %!in% c("", " ", "NC"))

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, Clinical[,c("PATIENT_ID", "CLAUDIN_SUBTYPE")], by.x = "Sample", by.y = "PATIENT_ID")


for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc_3p[i,]$End) %>% 
        group_by(Category, Allele, CLAUDIN_SUBTYPE) %>% summarise(n = n()) %>% mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
    
    sample_df_2 <- rbind.data.frame(sample_df_2, test)
}


CNA3p_Order <- CNA_Loc_3p %>% arrange(Start)

sample_df <- sample_df[order(match(sample_df[,1],CNA3p_Order$Hugo_Symbol)),]
sample_df_2$Gene <- factor(sample_df_2$Gene, levels = c(CNA3p_Order$Hugo_Symbol))

# PAM50
as_tibble(sample_df_2) %>% dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                     Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                     Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                     Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                     Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                     Category == "Amplification/Deletion" ~ "Amp/Del")) %>% mutate(Category = factor(Category, levels = c("Amp/Neut", "Del/Neut", 
                                                                                                                                                                          "Neut/Amp", "Neut/Del", 
                                                                                                                                                                          "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=Gene, y=n, color = Category, fill = Category)) + 
    geom_bar(stat = "identity")  + facet_grid(CLAUDIN_SUBTYPE~Allele) + ggtitle("Frequency of Breakpoints in Genes on Chromosome 11p by Allele and PAM50") + ylab("Frequency") + 
    xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                          axis.title.x = element_text(vjust = +2)) + 
    theme(legend.position = "top")


ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/Chromosome_11p_Barplot_Allele_PAM50.png", last_plot(), width = 12, height = 13)

## All Genes 

setwd("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/ASCAT_Data/Formatted_Data/")
ASCAT_Data_3Step_NoNeut<- read.csv2("Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
data(ens.gene.ann.hg19)

CNA_Loc <- read.delim("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/data/METABRIC_Data_Common/CNA_Data_Loc_Inc.txt", sep="\t")
CNA_Loc_3p <- CNA_Loc 

## 
result_list_1 <- list()
unique_samples <- unique(ASCAT_Data_3Step_NoNeut$Sample)
CHR <- c(1:22, "X")

Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut
Total_ASCAT3_Chr3$Breakpoint <- as.numeric(Total_ASCAT3_Chr3$Breakpoint)
CNA_Loc <- CNA_Loc %>% mutate(Chromosome = ifelse(Chromosome == 23, "X", Chromosome))

## Barplot?
# Plot of WHole Gene Genome

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "Chr" = character(0))

for (i in 1:nrow(CNA_Loc)) {
    
    test <- Total_ASCAT3_Chr3 %>% mutate(Category = factor(Category), Allele = factor(Allele))
    
    test <- test %>% 
        filter(as.numeric(Total_ASCAT3_Chr3$Breakpoint) >= CNA_Loc[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Breakpoint) <= CNA_Loc[i,]$End & Chr == CNA_Loc[i,]$Chromosome) 
    
    df_to_fill <- data.frame("Category" = rep(levels(test$Category), 2), "Allele" = rep(levels(test$Allele), 7), "Chr" = CNA_Loc[i,]$Chromosome, "n" = 0, "Gene" = CNA_Loc[i,]$Hugo_Symbol)
        
    if(nrow(test) == 0){
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
        
    } else {
        test <- test %>% 
            group_by(Category, Allele, Chr) %>% summarise(n = n()) %>%
            mutate(Gene = CNA_Loc[i,]$Hugo_Symbol,)
        
        for(j in 1:nrow(test)){
            df_to_fill$n[df_to_fill$Category == test$Category[j] & df_to_fill$Allele == test$Allele[j] & df_to_fill$Gene == test$Gene[j] & df_to_fill$Chr == test$Chr[j]] <- test$n[j]
        }
    
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
    } 
}


sample_df_2 <- read.delim("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/ASCAT_Data/Formatted_Data/Frequency_BP_Genes.txt", sep="\t")

Final <- merge(sample_df_2, CNA_Loc[,c(1,4,5)], by.x = "Gene" , by.y = "Hugo_Symbol")
sample_df_2_Order <- Final %>% arrange(Chr, Start)

sample_df_2_Order <- sample_df_2_Order %>% mutate(Chromosome = ifelse(Chromosome == 23, "X", Chromosome))
Test <- CNA_Loc %>% select(Hugo_Symbol, Chromosome, Start, band) %>% group_by(Chromosome) %>% arrange(Chromosome, Start)

Test1 <- Test %>% filter(band == "q") %>% group_by(Chromosome) %>% top_n(n=-1, wt=Start) %>% select(Chromosome, Hugo_Symbol)
Test1 <- Test1 %>% mutate(Chromosome = as.factor(Chromosome))

#Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
#Cent_Data <- merge(Test1, sample_df_2, by.x = "Chromosome", by.y = "Chr")
##Cent_Data <- Cent_Data %>% mutate(chrom = factor(chrom, levels = c(as.character(1:22), "X"))) %>% 
# #   dplyr::rename(., Chromosome = chrom)
#
#Cent_Data <- as_tibble(Cent_Data) %>% mutate(Chromosome = paste("Chr", Chromosome, sep= " ")) %>% 
#    mutate(Chromosome = factor(Chromosome, levels = c(paste("Chr", 1:22, sep= " "), "Chr X"))) 
#
#Cent_Data <- Cent_Data %>% mutate(Gene = as.character(Gene))
#Cent_Data <- merge(Cent_Data, CNA_Loc[,c("Hugo_Symbol", "Start")], by.x = "Gene", by.y = "Hugo_Symbol")
#Cent_Data <- Cent_Data %>% group_by(Chromosome) %>% arrange(Start)

sample_df_2_Order <- sample_df_2_Order %>% mutate(Chromosome = paste("Chr", Chromosome, sep= " ")) %>% mutate(Chromosome = factor(Chromosome, levels = c(paste("Chr", 1:22, sep= " "), "Chr X"))) 

Test1 <- Test1 %>% mutate(Chromosome = paste("Chr", Chromosome, sep= " ")) %>% mutate(Chromosome = factor(Chromosome, levels = c(paste("Chr", 1:22, sep= " "), "Chr X"))) 
Test1[23,"Chromosome"] <- "Chr X"

ggplot(sample_df_2_Order, aes(x=factor(Gene, levels = c(unique(sample_df_2_Order$Gene))), y=n)) + 
    geom_bar(stat = "identity")  +  facet_wrap(~Chromosome, scales = "free_x", ncol = 3) + 
    ggtitle("Frequency of Breakpoints in Genes across Chromosomes") + 
    ylab("") + 
    xlab("") + theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(legend.position = "none") + geom_vline(data = Test1, aes(xintercept = c(Hugo_Symbol)), color = "red", size = 0.2) + 
    theme(legend.position = "top") + scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=TRUE)


ggsave("/home/lydia/PhD_Insync/PhD_Project/PhD_Thesis_Submission/figures/Chapter_6/All_Genes_Barplot.png", last_plot(), width = 8,
       height = 10, units = c("in"), dpi = 700)

unique(sample_df_2_Order$Gene)
### Whole Genome (Barplots)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut 

## Combine Centromere info with Output
Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
Cent_Data <- merge(Centromere_Loc, Total_ASCAT3_Chr3, by.x = "chrom", by.y = "Chr")
Cent_Data <- Cent_Data %>% mutate(chrom = factor(chrom, levels = c(as.character(1:22), "X"))) %>% dplyr::rename(., Chromosome = chrom)

# Plot 
Cent_Data %>% ggplot(aes(x=as.numeric(Breakpoint), fill = Category, color = Category)) + geom_histogram(binwidth = 1000000) + 
    facet_wrap(~Chromosome, scales = "free_x") + geom_vline(data = Cent_Data, aes(xintercept = c(centromerStart)), color = "red") + 
    geom_vline(data = Cent_Data, aes(xintercept = c(centromerEnd)),  color = "red") + theme(legend.position = "top")

ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/scaled.png", last_plot(), width = 8, height = 10, units = c("in"), dpi = 700)

Cent_Data %>% ggplot(aes(x=as.numeric(Breakpoint), fill = Category, color = Category)) + geom_histogram(binwidth = 1000000) + 
    facet_wrap(~Chromosome, scales = "free") + geom_vline(data = Cent_Data, aes(xintercept = c(centromerStart)), color = "red") + 
    geom_vline(data = Cent_Data, aes(xintercept = c(centromerEnd)),  color = "red") + theme(legend.position = "top")

ggsave("/home/lydia/PhD_Insync/PhD_Project/GH_Thesis/figures/unscaled.png", width = 10, height = 6)