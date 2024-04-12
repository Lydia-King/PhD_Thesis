# Libraries 
library(tidyverse)
library(rCGH)
library(GeneBreak)
library(survival)
library(survminer)

# Data
ASCAT_Data_3Step_NoNeut<- read.csv2("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

Centromere_Loc <- hg19 %>% 
    mutate(chrom = ifelse(chrom == 23, "X", chrom))

data(ens.gene.ann.hg19)

CNA_Loc <- read.delim("../../data/Processed_Data/data_CNA_Loc_hg19.txt", sep="\t")
Gene_Of_Int <- CNA_Loc %>% filter(Hugo_Symbol == "OR52N1")

Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% 
    filter(Chr == "11") 

Pat <- Total_ASCAT3_Chr3[which(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= Gene_Of_Int[1,]$Start & as.numeric(Total_ASCAT3_Chr3$Changepoint) <= Gene_Of_Int[1,]$End & Total_ASCAT3_Chr3$Category == "Del/Amp" & Total_ASCAT3_Chr3$Allele == "Major"),]

# Clinical Data 
MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t")

MB_Data$PATIENT_ID <- gsub("\\-", ".", MB_Data$PATIENT_ID)
MB_Data <- MB_Data %>% mutate(OR52N1 = ifelse(PATIENT_ID %in% unique(Pat$Sample), "Del/Amp Changepoint", "No Del/Amp Changepoint"))    

## Plot KM Curves
fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ OR52N1, data = MB_Data)

survp <- ggsurvplot(
    fit,
    censor.shape = "+",
    xlab = "Survival time (months)",
    ylab = "Survival probability",
    data = MB_Data,
    size = 1,
    conf.int = F,
    pval = T,
    risk.table = F,
    legend = c(0.8, 0.80),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(A) DSS", "for Del/Amp Changepoint in OR52N1 (Major Allele) ", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)

png("../../figures/Chapter_6/survplot_OR52N1.png", width = 6, height = 4, units = "in", res = 600)
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
    legend = c(0.76, 0.82),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(A) DSS", "for Del/Amp Changepoint in OR52N1 (Major Allele) ", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)
dev.off()

# (B)
Gene_Of_Int <- CNA_Loc %>% filter(Hugo_Symbol == "TRIM5")

Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "11") 

Pat <- Total_ASCAT3_Chr3[which(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= Gene_Of_Int[1,]$Start & as.numeric(Total_ASCAT3_Chr3$Changepoint) <= Gene_Of_Int[1,]$End & Total_ASCAT3_Chr3$Category == "Del/Amp" & Total_ASCAT3_Chr3$Allele == "Major"),]

# Clinical Data 
MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t")

MB_Data$PATIENT_ID <- gsub("\\-", ".", MB_Data$PATIENT_ID)
MB_Data <- MB_Data %>% mutate(TRIM5 = ifelse(PATIENT_ID %in% unique(Pat$Sample), "Del/Amp Changepoint", "No Del/Amp Changepoint"))    

## Plot KM Curves
fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ TRIM5, data = MB_Data)

survp <- ggsurvplot(
    fit,
    censor.shape = "+",
    xlab = "Survival time (months)",
    ylab = "Survival probability",
    data = MB_Data,
    size = 1,
    conf.int = F,
    pval = T,
    risk.table = F,
    legend = c(0.8, 0.80),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(B) DSS", "for Del/Amp Changepoint in TRIM5 (Major Allele) ", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)

png("../../figures/Chapter_6/survplot_TRIM5.png", width = 6, height = 4, units = "in", res = 600)
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
    legend = c(0.76, 0.82),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(B) DSS", "for Del/Amp Changepoint in TRIM5 (Major Allele) ", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)
dev.off()

# (C)
Gene_Of_Int <- CNA_Loc %>% filter(Hugo_Symbol == "ALG1L2")

Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "3") 

Pat <- Total_ASCAT3_Chr3[which(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= Gene_Of_Int[1,]$Start & as.numeric(Total_ASCAT3_Chr3$Changepoint) <= Gene_Of_Int[1,]$End & Total_ASCAT3_Chr3$Category == "Del/Amp"),]

# Clinical Data 
MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t")

MB_Data$PATIENT_ID <- gsub("\\-", ".", MB_Data$PATIENT_ID)
MB_Data <- MB_Data %>% mutate(ALG1L2 = ifelse(PATIENT_ID %in% unique(Pat$Sample), "Del/Amp Changepoint", "No Del/Amp Changepoint"))    

## Plot KM Curves
fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ ALG1L2, data = MB_Data)

png("../../figures/Chapter_6/survplot_ALG.png", width = 6, height = 4, units = "in", res = 600)
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
    legend = c(0.76, 0.82),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("(C) DSS", "for Del/Amp Changepoint in ALG1L2 ", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)
dev.off()
