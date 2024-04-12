library(tidyverse)
library(rCGH)
library(GeneBreak)
library(survival)
library(survminer)

# Data
ASCAT_Data_1 <- read.csv2("../../data/ASCAT_Data/3Step/ASCAT_Part_1.segments.txt", sep="\t")

ASCAT_Data_3Step_NoNeut <- read.csv2("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
data(ens.gene.ann.hg19)

Chr_Info <- as.data.frame(ASCAT_Data_1 %>% group_by(chr) %>% summarise(start = min(startpos), end = max(endpos))) %>% 
    mutate(chr = factor(chr, levels = c(1:22, "X")))

Chr_Info <- Chr_Info[order(Chr_Info$chr),]

## Segment 27 Chromosome 1 
Chr_K <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == 1)
seq_segment <- unique(c(seq(Chr_Info[1, 2], Chr_Info[1, 3], 5000000),  Chr_Info[1, 3]))

start <- seq_segment[27] 
end <- seq_segment[27+1]

Pat <- Chr_K[which(as.numeric(Chr_K$Changepoint) >= start & as.numeric(Chr_K$Changepoint) <= end & Chr_K$Category == "Neut/Amp" & Chr_K$Allele == "Major"),] 

MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t")

MB_Data$PATIENT_ID <- gsub("\\-", ".", MB_Data$PATIENT_ID)
MB_Data <- MB_Data %>% mutate(Chr1_Seg27 = ifelse(PATIENT_ID %in% unique(Pat$Sample), "Changepoint", "NoChangepoint"))   

fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ Chr1_Seg27, data = MB_Data)

png("../../figures/Chapter_6/survplot_Chr1Seg27.png", width = 7, height = 4, units = "in", res = 600)
ggsurvplot(
    fit,
    censor.shape = "+",
    xlab = "Survival time (months)",
    ylab = "Survival probability",
    data = MB_Data,
    size = 1,
    conf.int = F,
    pval = T,
    risk.table = F,
    legend = c(0.84, 0.80),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(A) DSS", "for Neut/Amp Changepoints in Chromosome 1 Segment 27", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)
dev.off()

## Segment 31 Chromosome 1 
Chr_K <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == 1)
seq_segment <- unique(c(seq(Chr_Info[1, 2], Chr_Info[1, 3], 5000000),  Chr_Info[1, 3]))

start <- seq_segment[31] 
end <- seq_segment[31+1]

Pat <- Chr_K[which(as.numeric(Chr_K$Changepoint) >= start & as.numeric(Chr_K$Changepoint) <= end & Chr_K$Category == "Del/Amp" & Chr_K$Allele == "Major"),] 

MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t")

MB_Data$PATIENT_ID <- gsub("\\-", ".", MB_Data$PATIENT_ID)
MB_Data <- MB_Data %>% mutate(Chr1_Seg31 = ifelse(PATIENT_ID %in% unique(Pat$Sample), "Del/Amp Changepoint", "No Del/Amp Changepoint"))   

fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ Chr1_Seg31, data = MB_Data)

png("../../figures/Chapter_6/survplot_Chr1Seg31.png", width = 7, height = 4, units = "in", res = 600)
ggsurvplot(
    fit,
    censor.shape = "+",
    xlab = "Survival time (months)",
    ylab = "Survival probability",
    data = MB_Data,
    size = 1,
    conf.int = F,
    pval = T,
    risk.table = F,
    legend = c(0.79, 0.83),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(B) DSS", "for Del/Amp Changepoints in Chromosome 1 Segment 31 (Major)", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)
dev.off()

## Segment 9 Chromosome 16 
Chr_K <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == 16)
seq_segment <- unique(c(seq(Chr_Info[16, 2], Chr_Info[16, 3], 5000000),  Chr_Info[16, 3]))

start <- seq_segment[9] 
end <- seq_segment[9+1]

Pat <- Chr_K[which(as.numeric(Chr_K$Changepoint) >= start & as.numeric(Chr_K$Changepoint) <= end & Chr_K$Category == "Neut/Del" & Chr_K$Allele == "Minor"),] 

# Clinical Data 
MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t")

MB_Data$PATIENT_ID <- gsub("\\-", ".", MB_Data$PATIENT_ID)
MB_Data <- MB_Data %>% mutate(Chr16_Seg9 = ifelse(PATIENT_ID %in% unique(Pat$Sample), "Changepoint", "NoChangepoint"))    

## Plot KM Curves

fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ Chr16_Seg9, data = MB_Data)

png("../../figures/Chapter_6/survplot_Chr16Seg9.png", width = 7, height = 4, units = "in", res = 600)
ggsurvplot(
    fit,
    censor.shape = "+",
    xlab = "Survival time (months)",
    ylab = "Survival probability",
    data = MB_Data,
    size = 1,
    conf.int = F,
    pval = T,
    risk.table = F,
    legend = c(0.84, 0.80),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(B) DSS", "for Neut/Del Changepoints in Chromosome 16 Segment 9", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)
dev.off()

