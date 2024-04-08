# Chapter 2: CNA Score and Burden Metric Calculation (All Data and Complete Case)
## Load up necessary libraries
library(dplyr)
library(operator.tools)

## Load up data (CNA file from cBioPortal 2021)
CNA <-
    read.delim(
        "../../data/METABRIC_2021/data_CNA.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

# CNA Score
## Create dataframe to store CNA Score for all data and CC data
CNA_Score_Metrics_All <-
    data.frame("PATIENT_ID" = colnames(CNA[, 3:ncol(CNA)]))
CNA_Score_Metrics_CCA <-
    data.frame("PATIENT_ID" = colnames(CNA[, 3:ncol(CNA)]))

## All Data
### Absolute CNA Score
CNA_Score_Metrics_All$CNA_Score_All <-
    colSums(abs(CNA[, 3:ncol(CNA)]), na.rm = T)

### CNA Amplification Score
CNA_Score_Metrics_All$Amp_Score_All <-
    apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x[x > 0], na.rm = T))

### CNA Deletion Score
CNA_Score_Metrics_All$Del_Score_All <-
    apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        abs(sum(x[x < 0], na.rm = T)))

### Difference Score
CNA_Score_Metrics_All <- CNA_Score_Metrics_All %>%
    mutate(Difference_Score_All = Amp_Score_All - Del_Score_All)

### Percentage of Alterations Amplified/Deleted
CNA_Score_Metrics_All <- CNA_Score_Metrics_All %>%
    mutate(Percentage_Score_Amp_All = (Amp_Score_All / CNA_Score_All) *
               100) %>%
    mutate(Percentage_Score_Del_All = (Del_Score_All / CNA_Score_All) *
               100)
CNA_Score_Metrics_All[CNA_Score_Metrics_All == "NaN"] <- 0

## Complete Case (CC) Data
### Absolute CNA Score
CNA_Score_Metrics_CCA$CNA_Score_CCA <-
    colSums(abs(CNA[, 3:ncol(CNA)]))
CNA_Score_Metrics_CCA <- na.omit(CNA_Score_Metrics_CCA)

### CNA Amplification Score
CNA_Score_Metrics_CCA$Amp_Score_CCA <-
    na.omit(apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x[x > 0])))

### CNA Deletion Score
CNA_Score_Metrics_CCA$Del_Score_CCA <-
    na.omit(apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        abs(sum(x[x < 0]))))

### Difference Score
CNA_Score_Metrics_CCA <- CNA_Score_Metrics_CCA %>%
    mutate(Difference_Score_CCA = Amp_Score_CCA - Del_Score_CCA)

### Percentage of Alterations Amplified/Deleted
CNA_Score_Metrics_CCA <- CNA_Score_Metrics_CCA  %>%
    mutate(Percentage_Score_Amp_CCA = (Amp_Score_CCA / CNA_Score_CCA) *
               100) %>%
    mutate(Percentage_Score_Del_CCA = (Del_Score_CCA / CNA_Score_CCA) *
               100)
CNA_Score_Metrics_CCA[CNA_Score_Metrics_CCA == "NaN"] <- 0


# CNA Burden
## Create dataframe to store CNA Burden for all data and CC data
CNA_Burden_Metrics_All <-
    data.frame("PATIENT_ID" = colnames(CNA[, 3:ncol(CNA)]))
CNA_Burden_Metrics_CCA <-
    data.frame("PATIENT_ID" = colnames(CNA[, 3:ncol(CNA)]))

## All Data
### CNA Burden
Number_Genes_All <-
    sapply(X = CNA[, 3:ncol(CNA)], function(y)
        sum(length(which(!is.na(
            y
        )))))
Affected_Genes_All <-
    apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x != 0, na.rm = T))
CNA_Burden_Metrics_All$CNA_Burden_All <-
    (Affected_Genes_All / Number_Genes_All) * 100

### CNA Amplification Burden
Amp_Genes_All <-
    apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x > 0, na.rm = T))
CNA_Burden_Metrics_All$Amp_Burden_All = (Amp_Genes_All / Number_Genes_All) *
    100

### CNA Deletion Burden
Del_Genes_All <-
    apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x < 0, na.rm = T))
CNA_Burden_Metrics_All$Del_Burden_All = (Del_Genes_All / Number_Genes_All) *
    100

### Difference Burden
CNA_Burden_Metrics_All <- CNA_Burden_Metrics_All %>%
    mutate(Difference_Burden_All = Amp_Burden_All - Del_Burden_All)

### Percentage of Alterations Amplified/Deleted
CNA_Burden_Metrics_All <- CNA_Burden_Metrics_All %>%
    mutate(Percentage_Burden_Amp_All = (Amp_Burden_All / CNA_Burden_All) *
               100) %>%
    mutate(Percentage_Burden_Del_All = (Del_Burden_All / CNA_Burden_All) *
               100)
CNA_Burden_Metrics_All[CNA_Burden_Metrics_All == "NaN"] <- 0

## Complete Case (CC) Data
### CNA Burden
Number_Genes_CCA <- nrow(CNA)
Affected_Genes_CCA <-
    apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x != 0))
CNA_Burden_Metrics_CCA$CNA_Burden_CCA = (Affected_Genes_CCA / Number_Genes_CCA) *
    100
CNA_Burden_Metrics_CCA <- na.omit(CNA_Burden_Metrics_CCA)

### CNA Amplification Burden
Amp_Genes_CCA <-
    na.omit(apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x > 0)))
CNA_Burden_Metrics_CCA$Amp_Burden_CCA = (Amp_Genes_CCA / Number_Genes_CCA) *
    100

### CNA Deletion Burden
Del_Genes_CCA <-
    na.omit(apply(X = CNA[, 3:ncol(CNA)], MARGIN = 2, function(x)
        sum(x < 0)))
CNA_Burden_Metrics_CCA$Del_Burden_CCA = (Del_Genes_CCA / Number_Genes_CCA) *
    100

# Difference Burden
CNA_Burden_Metrics_CCA <- CNA_Burden_Metrics_CCA %>%
    mutate(Difference_Burden_CCA = Amp_Burden_CCA - Del_Burden_CCA)

# Percentage of Alterations Amplified/Deleted (CCA)
CNA_Burden_Metrics_CCA <- CNA_Burden_Metrics_CCA  %>%
    mutate(Percentage_Burden_Amp_CCA = (Amp_Burden_CCA / CNA_Burden_CCA) *
               100) %>%
    mutate(Percentage_Burden_Del_CCA = (Del_Burden_CCA / CNA_Burden_CCA) *
               100)
CNA_Burden_Metrics_CCA[CNA_Burden_Metrics_CCA == "NaN"] <- 0

# Save CNA Score and Burden Metrics
write.table(
    CNA_Score_Metrics_All,
    file = "../../data/Processed_Data/CNA_Score_All.txt",
    sep = "\t",
    row.names = FALSE
)
write.table(
    CNA_Score_Metrics_CCA,
    file = "../../data/Processed_Data/CNA_Score_CCA.txt",
    sep = "\t",
    row.names = FALSE
)
write.table(
    CNA_Burden_Metrics_All,
    file = "../../data/Processed_Data/CNA_Burden_All.txt",
    sep = "\t",
    row.names = FALSE
)
write.table(
    CNA_Burden_Metrics_CCA,
    file = "../../data/Processed_Data/CNA_Burden_CCA.txt",
    sep = "\t",
    row.names = FALSE
)
