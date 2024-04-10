# Chapter 5: Simulation Study Data

## Load up Libraries
library(tidyverse)
library(ds4psy)
library(rbenchmark)
library(operator.tools)
library(truncnorm)

## Load up functions to simulate data
source("Chapter5_Simulate_Data_Functions.R")

## Scenario 1
# Matrix of Distributions
MU_SD_Neut <- matrix(c(27641, 35854, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10), nrow = 8, ncol = 2, byrow = T)
rownames(MU_SD_Neut) <- c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <- matrix(c(700, 1, 15249 , 28815, 700, 1, 100, 1, 100, 1, 100, 1, 100, 1, 100, 1), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Amp) <- c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <- matrix(c(50, 1, 50, 1, 50, 1, 25, 1, 25, 1, 25, 1, 25, 1, 25, 1), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Del) <- c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

Percentage_Neut_Amp <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)
Sample_Size <- c(20, 50, 80, 100, 200, 500)

Same_Structure_list <- list()
Total_Dataset1_list <- list()
Samp_list <- list()

set.seed(2433)
for(i in 1:length(Percentage_Neut_Amp)){
    for(Samp in c(Sample_Size)){
        dataset1_Major <- Simulate_Data_NoNeut(Num_Samples = Samp, Lengths_Matrix = Total_Mu_Sd1,
                                        Percent_Neutral = c(100-Percentage_Neut_Amp[i]), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(Percentage_Neut_Amp[i]),
                                        Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 2,
                                        Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 1,
                                        Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                        Prob_Amp_Del_to_Del_Neut = 0, Allele = "Major")
        
        dataset1_Minor <- Simulate_Data_NoNeut(Num_Samples = Samp,Lengths_Matrix = Total_Mu_Sd1,
                                        Percent_Neutral = c(100), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                        Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 0,
                                        Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                        Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                        Prob_Amp_Del_to_Del_Neut = 0,  Allele = "Minor")
        
        Total_Dataset1 <- lapply(seq_along(dataset1_Major), function(x) rbind(dataset1_Major[[x]], dataset1_Minor[[x]]))
        Total_Dataset1 <- lapply(Total_Dataset1, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:20, sep = " "))), Allele)),
        ])
        
        Samp_list[[paste("Samplesize_", Samp, sep="")]] <- Total_Dataset1
    }
    Total_Dataset1_list[[paste("Percentage_", Percentage_Neut_Amp[i], sep="")]] <- Samp_list
}

Multiple_Dataset_Scen1 <- Total_Dataset1_list

## Scenario 2
# Matrix of Distributions
MU_SD_Neut <- matrix(c(27641, 35854, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Neut) <- c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <- matrix(c(700, 1, 15249 , 28815, 700, 1, 100, 1, 22777, 35235, 100, 1, 100, 1, 100, 1), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Amp) <- c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <- matrix(c(50, 1, 50, 1, 50, 1, 25, 1, 25, 1, 25, 1, 25, 1, 9769, 19739), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Del) <- c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

Percentage_A <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)
Percentage_B <- c(90, 80, 70, 60, 50, 40, 30, 20, 10)
Sample_Size <- c(20, 50, 80, 100, 200, 500)

Same_Structure_list <- list()
Total_Dataset1_list <- list()
Samp_list <- list()

set.seed(2433)
for(i in 1:length(Percentage_Neut_Amp)){
    for(Samp in c(Sample_Size)){
        
        ## A
        Samp_Count_A <- Samp*Percentage_A[i]/100
        dataset1_Major_Samp_A1 <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_A, Lengths_Matrix = Total_Mu_Sd1,
                                                Percent_Neutral = c(100), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                                Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 2,
                                                Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 1,
                                                Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                                Prob_Amp_Del_to_Del_Neut = 0, Allele = c("Minor"), Sample_Same = T)
        
        dataset1_Minor_Samp_A2 <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_A, Lengths_Matrix = Total_Mu_Sd1,
                                                Percent_Neutral = c(0), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(100),
                                                Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 3,
                                                Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                                Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                                Prob_Amp_Del_to_Del_Neut = 1, Allele = c("Major"), Sample_Same = T)
        
        Total_Dataset1_A <- lapply(seq_along(dataset1_Major_Samp_A1), function(x) rbind(dataset1_Major_Samp_A1[[x]], dataset1_Minor_Samp_A2[[x]]))
        Total_Dataset1_A <- lapply(Total_Dataset1_A, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:Samp_Count_A, sep = " "))), Allele)),
        ])
        
        ## B
        Samp_Count_B <- Samp*Percentage_B[i]/100
        dataset1_Major_Samp_B1 <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_B, Lengths_Matrix = Total_Mu_Sd1,
                                                Percent_Neutral = c(0), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(100),
                                                Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 2,
                                                Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 1,
                                                Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                                Prob_Amp_Del_to_Del_Neut = 0, Allele = c("Major"), Sample_Same = T, samp_start = Samp_Count_A)
        
        dataset1_Minor_Samp_B2 <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_B, Lengths_Matrix = Total_Mu_Sd1,
                                                Percent_Neutral = c(100), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                                Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 2,
                                                Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 1,
                                                Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                                Prob_Amp_Del_to_Del_Neut = 0, Allele = c("Minor"), Sample_Same = T, samp_start = Samp_Count_A)
        
        Total_Dataset1_B <- lapply(seq_along(dataset1_Major_Samp_B1), function(x) rbind(dataset1_Major_Samp_B1[[x]], dataset1_Minor_Samp_B2[[x]]))
        Total_Dataset1_B <- lapply(Total_Dataset1_B, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:800, sep = " "))), Allele)),
        ])
        
        Total_Dataset_All <- lapply(seq_along(Total_Dataset1_A), function(x) rbind(Total_Dataset1_A[[x]], Total_Dataset1_B[[x]]))
        Total_Dataset_All <- lapply(Total_Dataset_All, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:800, sep = " "))), Allele)),
        ])
        
        Samp_list[[paste("Samplesize_", Samp, sep="")]] <- Total_Dataset_All
    }
    Total_Dataset1_list[[paste("Percentage_", Percentage_A[i], sep="")]] <- Samp_list
}

Multiple_Dataset_Scen2 <- Total_Dataset1_list

## Scenario 3
## (A) Samp1 - Samp4
MU_SD_Neut <- matrix(c(27641, 35854, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10, 20, 10), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Neut) <- c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
colnames(MU_SD_Neut) <- c("Mu", "Sd")

MU_SD_Amp <- matrix(c(700, 1, 15249 , 28815, 700, 1, 100, 1, 22777, 35235, 100, 1, 100, 1, 100, 1), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Amp) <- c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
colnames(MU_SD_Amp) <- c("Mu", "Sd")

MU_SD_Del <- matrix(c(50, 1, 50, 1, 50, 1, 25, 1, 25, 1, 25, 1, 25, 1, 9769, 19739), nrow = 8, ncol = 2, byrow = T)

rownames(MU_SD_Del) <- c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
colnames(MU_SD_Del) <- c("Mu", "Sd")

Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)

Percentage_A <- c(20, 20, 20, 20, 20, 20, 20) 
Percentage_B <- c(10, 20, 30, 40, 50, 60, 70)
Percentage_C <- c(70, 60, 50, 40, 30, 20, 10)
Sample_Size <- c(20, 50, 80, 100, 200, 500)

Same_Structure_list <- list()
Total_Dataset1_list <- list()
Samp_list <- list()

set.seed(2433)
for(i in 1:length(Percentage_A)){
    for(Samp in c(Sample_Size)){
        
        ## A
        Samp_Count_A <- Samp*Percentage_A[i]/100
        dataset1_Major_A <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_A, Lengths_Matrix = Total_Mu_Sd1,
                                          Percent_Neutral = c(100), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                          Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 3,
                                          Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                          Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                          Prob_Amp_Del_to_Del_Neut = 1,  Allele = "Major")
        
        dataset1_Minor_A <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_A, Lengths_Matrix = Total_Mu_Sd1,
                                          Percent_Neutral = c(100), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                          Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 3,
                                          Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                          Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                          Prob_Amp_Del_to_Del_Neut = 1,  Allele = "Minor")
        
        Total_Dataset_A <- lapply(seq_along(dataset1_Major_A), function(x) rbind(dataset1_Major_A[[x]], dataset1_Minor_A[[x]]))
        Total_Dataset_A <- lapply(Total_Dataset_A, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:800, sep = " "))), Allele)),
        ])
        
        
        ## B
        Samp_Count_B <- Samp*Percentage_B[i]/100
        dataset1_Major_B <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_B, Lengths_Matrix = Total_Mu_Sd1,
                                          Percent_Neutral = c(100), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                          Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 2,
                                          Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 1,
                                          Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                          Prob_Amp_Del_to_Del_Neut = 0, Allele = c("Minor"), Sample_Same = T, samp_start = Samp_Count_A)
        
        dataset1_Minor_B <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_B, Lengths_Matrix = Total_Mu_Sd1,
                                          Percent_Neutral = c(0), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(100),
                                          Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 3,
                                          Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                          Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                          Prob_Amp_Del_to_Del_Neut = 1, Allele = c("Major"), Sample_Same = T, samp_start = Samp_Count_A)
        
        Total_Dataset_B <- lapply(seq_along(dataset1_Major_B), function(x) rbind(dataset1_Major_B[[x]], dataset1_Minor_B[[x]]))
        Total_Dataset_B <- lapply(Total_Dataset_B, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:800, sep = " "))), Allele)),
        ])
        
        ## C
        ## Major
        MU_SD_Neut <- matrix(c(27641, 35854, 8997, 18675, 20, 10, 20, 10, 31129, 38125, 20, 10, 20, 10, 20, 10), nrow = 8, ncol = 2, byrow = T)

        rownames(MU_SD_Neut) <- c("SNNA", "SNND", "ANNA", "ANND", "DNND", "DNNA", "ANNE", "DNNE")
        colnames(MU_SD_Neut) <- c("Mu", "Sd")
        
        MU_SD_Amp <- matrix(c(700, 1, 15249 , 28815, 700, 1, 100, 1, 68331, 35235, 100, 1, 100, 1, 100, 1), nrow = 8, ncol = 2, byrow = T)
        
        rownames(MU_SD_Amp) <- c("SAAN", "NAAN", "NAAE", "SAAD", "NAAD", "DAAD", "DAAE", "DAAN")
        colnames(MU_SD_Amp) <- c("Mu", "Sd")
        
        MU_SD_Del <- matrix(c(8997, 18675, 8997, 18675, 50, 1, 25, 1, 25, 1, 25, 1, 25, 1, 3*9769, 19739), nrow = 8, ncol = 2, byrow = T)
        
        rownames(MU_SD_Del) <- c("SDDN", "NDDN", "NDDE", "SDDA", "NDDA", "ADDA", "ADDE", "ADDN")
        colnames(MU_SD_Del) <- c("Mu", "Sd")
        
        Total_Mu_Sd1 <- rbind.data.frame(MU_SD_Neut, MU_SD_Amp, MU_SD_Del)
        
        Samp_Count_C <- Samp*Percentage_C[i]/100
        
        dataset1_Major_C <-Simulate_Data_NoNeut(Num_Samples = Samp_Count_C, Lengths_Matrix = Total_Mu_Sd1,
                                         Percent_Neutral = c(0), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(100),
                                         Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(0), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 3,
                                         Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                         Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 0, Prob_Del_Amp_to_Amp_Del = 0,
                                         Prob_Amp_Del_to_Del_Neut = 1, Allele = c("Major"), Sample_Same = T, samp_start = Samp_Count_A+Samp_Count_B)
        
        
        ## Minor
        dataset1_Minor_C <- Simulate_Data_NoNeut(Num_Samples = Samp_Count_C, Lengths_Matrix = Total_Mu_Sd1,
                                          Percent_Neutral = c(0), Percent_Amp_Neutral = c(0), Percent_Neutral_Amp = c(0),
                                          Percent_Del_Neutral = c(0), Percent_Neutral_Del = c(100), Percent_Amp_Del_FP = c(0), Percent_Del_Amp_FP = c(0), Max_Num_CP = 8,
                                          Prob_Amp_Neut_to_Neut_Amp = 0, Prob_Neut_Amp_to_Amp_Neut = 0,
                                          Prob_Del_Neut_to_Neut_Amp = 0, Prob_Neut_Del_to_Del_Neut = 1, Prob_Del_Amp_to_Amp_Del = 0,
                                          Prob_Amp_Del_to_Del_Neut = 0, Allele = c("Minor"), Sample_Same = T, samp_start = Samp_Count_A+Samp_Count_B)
        
        
        Total_Dataset_C <- lapply(seq_along(dataset1_Major_C), function(x) rbind(dataset1_Major_C[[x]], dataset1_Minor_C[[x]]))
        Total_Dataset_C <- lapply(Total_Dataset_C, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:800, sep = " "))), Allele)),
        ])
        
        Total_Dataset_3 <- lapply(seq_along(Total_Dataset_A), function(x) rbind(Total_Dataset_A[[x]], Total_Dataset_B[[x]]))
        Total_Dataset_3 <- lapply(seq_along(Total_Dataset_3), function(x) rbind(Total_Dataset_3[[x]], Total_Dataset_C[[x]]))
        
        Total_Dataset_3 <- lapply(Total_Dataset_3, function(x) x[
            with(x, order(factor(Sample, levels = c(paste("Sample", 1:800, sep = " "))), Allele)),
        ])
        
        Samp_list[[paste("Samplesize_", Samp, sep="")]] <- Total_Dataset_3
    }
    Total_Dataset1_list[[paste("Percentage_", Percentage_B[i], sep="")]] <- Samp_list
}

Multiple_Dataset_Scen3 <- Total_Dataset1_list

## Write Out
Simulation_Data <- list(Multiple_Dataset_Scen1, Multiple_Dataset_Scen2, Multiple_Dataset_Scen3)
saveRDS(Simulation_Data, file = "../../data/Chapter_5/Simulation_Study_Data.rds")
