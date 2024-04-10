# Chapter 5: Setup Example Dataset

## Load up libraries
library(tidyverse)
library(ds4psy)
library(rbenchmark)
library(operator.tools)
library(truncnorm)
library(gt)

## Load up functions to simulate data
source("Chapter5_Simulate_Data_Functions.R")

## Create Datasets - No Neut
### Profile (A) 
MU_SD_Neut <-
    matrix(
        c(27641, 35854, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Neut) <-
    c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <-
    matrix(
        c(
            700,
            1,
            15249,
            28815,
            700,
            1,
            100,
            1,
            22777,
            35235,
            100,
            1,
            100,
            1,
            100,
            1
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Amp) <-
    c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <-
    matrix(
        c(50, 1, 50, 1, 50, 1, 25, 1, 25, 1, 25, 1, 25, 1, 9769, 19739),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Del) <-
    c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

dataset1_Major_A <-
    Simulate_Data_NoNeut(
        Num_Samples = 4,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(100),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = "Major",
        seed = 2344,
        Chr_Start = 1,
        Chr_End = 250000
    )

dataset1_Minor_A <-
    Simulate_Data_NoNeut(
        Num_Samples = 4,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(100),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = "Minor",
        seed = 2344,
        Chr_Start = 1,
        Chr_End = 250000
    )

Total_Dataset_A <-
    lapply(seq_along(dataset1_Major_A), function(x)
        rbind(dataset1_Major_A[[x]], dataset1_Minor_A[[x]]))
Total_Dataset_A <- lapply(Total_Dataset_A, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

### Profile (B)
dataset1_Major_B <-
    Simulate_Data_NoNeut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(100),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 2,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 1,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 0,
        Allele = c("Minor"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 4,
        Chr_Start = 1,
        Chr_End = 250000
    )

dataset1_Minor_B <-
    Simulate_Data_NoNeut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(0),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(100),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = c("Major"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 4,
        Chr_Start = 1,
        Chr_End = 250000
    )

Total_Dataset_B <-
    lapply(seq_along(dataset1_Major_B), function(x)
        rbind(dataset1_Major_B[[x]], dataset1_Minor_B[[x]]))
Total_Dataset_B <- lapply(Total_Dataset_B, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)), ])

### Profile (C)
MU_SD_Neut <-
    matrix(
        c(
            27641,
            35854,
            8997,
            18675,
            20,
            10,
            20,
            10,
            31129,
            38125,
            20,
            10,
            20,
            10,
            20,
            10
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Neut) <-
    c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <-
    matrix(
        c(
            700,
            1,
            15249,
            28815,
            700,
            1,
            100,
            1,
            68331,
            35235,
            100,
            1,
            100,
            1,
            100,
            1
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Amp) <-
    c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <-
    matrix(
        c(
            8997,
            18675,
            8997,
            18675,
            50,
            1,
            25,
            1,
            25,
            1,
            25,
            1,
            25,
            1,
            3 * 9769,
            19739
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Del) <-
    c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

dataset1_Major_C <-
    Simulate_Data_NoNeut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(0),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(100),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = c("Major"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 12,
        Chr_Start = 1,
        Chr_End = 250000
    )


dataset1_Minor_C <-
    Simulate_Data_NoNeut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(0),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(100),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 8,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 1,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 0,
        Allele = c("Minor"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 12,
        Chr_Start = 1,
        Chr_End = 250000
    )


Total_Dataset_C <-
    lapply(seq_along(dataset1_Major_C), function(x)
        rbind(dataset1_Major_C[[x]], dataset1_Minor_C[[x]]))

Total_Dataset_C <- lapply(Total_Dataset_C, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

## Combine all
Total_Dataset_3 <-
    lapply(seq_along(Total_Dataset_A), function(x)
        rbind(Total_Dataset_A[[x]], Total_Dataset_B[[x]]))
Total_Dataset_3 <-
    lapply(seq_along(Total_Dataset_3), function(x)
        rbind(Total_Dataset_3[[x]], Total_Dataset_C[[x]]))

Total_Dataset_3 <- lapply(Total_Dataset_3, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

DF3 <-
    within(Total_Dataset_3[[1]],
           Category <- relevel(factor(Type), ref = "NoChangepoint"))

DF3 <- DF3 %>% mutate("TS" = TS, "TE" = TE)

## Create Datasets - Neut
### Profile (A) 
MU_SD_Neut <-
    matrix(
        c(27641, 35854, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Neut) <-
    c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <-
    matrix(
        c(
            700,
            1,
            15249,
            28815,
            700,
            1,
            100,
            1,
            22777,
            35235,
            100,
            1,
            100,
            1,
            100,
            1
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Amp) <-
    c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <-
    matrix(
        c(50, 1, 50, 1, 50, 1, 25, 1, 25, 1, 25, 1, 25, 1, 9769, 19739),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Del) <-
    c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

dataset1_Major_A <-
    Simulate_Data_Neut(
        Num_Samples = 4,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(100),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = "Major",
        seed = 2344,
        Chr_Start = 1,
        Chr_End = 250000
    )

dataset1_Minor_A <-
    Simulate_Data_Neut(
        Num_Samples = 4,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(100),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = "Minor",
        seed = 2344,
        Chr_Start = 1,
        Chr_End = 250000
    )

Total_Dataset_A <-
    lapply(seq_along(dataset1_Major_A), function(x)
        rbind(dataset1_Major_A[[x]], dataset1_Minor_A[[x]]))
Total_Dataset_A <- lapply(Total_Dataset_A, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

### Profile (B)
dataset1_Major_B <-
    Simulate_Data_Neut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(100),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 2,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 1,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 0,
        Allele = c("Minor"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 4,
        Chr_Start = 1,
        Chr_End = 250000
    )

dataset1_Minor_B <-
    Simulate_Data_Neut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(0),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(100),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = c("Major"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 4,
        Chr_Start = 1,
        Chr_End = 250000
    )

Total_Dataset_B <-
    lapply(seq_along(dataset1_Major_B), function(x)
        rbind(dataset1_Major_B[[x]], dataset1_Minor_B[[x]]))
Total_Dataset_B <- lapply(Total_Dataset_B, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

### Profile (C)
MU_SD_Neut <-
    matrix(
        c(
            27641,
            35854,
            8997,
            18675,
            20,
            10,
            20,
            10,
            31129,
            38125,
            20,
            10,
            20,
            10,
            20,
            10
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Neut) <-
    c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <-
    matrix(
        c(
            700,
            1,
            15249,
            28815,
            700,
            1,
            100,
            1,
            68331,
            35235,
            100,
            1,
            100,
            1,
            100,
            1
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Amp) <-
    c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <-
    matrix(
        c(
            8997,
            18675,
            8997,
            18675,
            50,
            1,
            25,
            1,
            25,
            1,
            25,
            1,
            25,
            1,
            3 * 9769,
            19739
        ),
        nrow = 8,
        ncol = 2,
        byrow = T
    )

rownames(MU_SD_Del) <-
    c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

dataset1_Major_C <-
    Simulate_Data_Neut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(0),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(100),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(0),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 3,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 0,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 1,
        Allele = c("Major"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 12,
        Chr_Start = 1,
        Chr_End = 250000
    )


dataset1_Minor_C <-
    Simulate_Data_Neut(
        Num_Samples = 8,
        Lengths_Matrix = Total_Mu_Sd1,
        Percent_Neutral = c(0),
        Percent_Amp_Neutral = c(0),
        Percent_Neutral_Amp = c(0),
        Percent_Del_Neutral = c(0),
        Percent_Neutral_Del = c(100),
        Percent_Amp_Del_FP = c(0),
        Percent_Del_Amp_FP = c(0),
        Max_Num_CP = 8,
        Prob_Amp_Neut_to_Neut_Amp = 0,
        Prob_Neut_Amp_to_Amp_Neut = 0,
        Prob_Del_Neut_to_Neut_Amp = 0,
        Prob_Neut_Del_to_Del_Neut = 1,
        Prob_Del_Amp_to_Amp_Del = 0,
        Prob_Amp_Del_to_Del_Neut = 0,
        Allele = c("Minor"),
        seed = 2344,
        Sample_Same = T,
        samp_start = 12,
        Chr_Start = 1,
        Chr_End = 250000
    )


Total_Dataset_C <-
    lapply(seq_along(dataset1_Major_C), function(x)
        rbind(dataset1_Major_C[[x]], dataset1_Minor_C[[x]]))
Total_Dataset_C <- lapply(Total_Dataset_C, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

## Combine all
Total_Dataset_3 <-
    lapply(seq_along(Total_Dataset_A), function(x)
        rbind(Total_Dataset_A[[x]], Total_Dataset_B[[x]]))
Total_Dataset_3 <-
    lapply(seq_along(Total_Dataset_3), function(x)
        rbind(Total_Dataset_3[[x]], Total_Dataset_C[[x]]))

Total_Dataset_3 <- lapply(Total_Dataset_3, function(x)
    x[with(x, order(factor(Sample, levels = c(
        paste("Sample", 1:20, sep = " ")
    )), Allele)),])

DF3_Neut <-
    within(Total_Dataset_3[[1]],
           Category <- relevel(factor(Type), ref = "NoChangepoint"))
DF3_Neut <-
    DF3_Neut %>% mutate("TS" = TS, "TE" = TE)

## Write Out Data
write.table(DF3, file = "../../data/Chapter_5/Example_Data_1_NoNeut.txt", sep = "\t", row.names = FALSE)
write.table(DF3_Neut, file = "../../data/Chapter_5/Example_Data_1_Neut.txt", sep = "\t", row.names = FALSE)
