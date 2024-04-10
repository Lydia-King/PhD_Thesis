# Chapter 5: Run Simulation Study

## Load Libraries 
library(tidyverse)
library(MCMCglmm)

## Source functions 
source("Chapter5_Simulation_Study_Functions.R")

## Load up datasets 
Multiple_Dataset_Scen <- readRDS(file = "../../data/Chapter_5/Simulation_Study_Data.rds")

## Set up datasets
Scenario1 <- Multiple_Dataset_Scen[[1]]
Scenario2 <- Multiple_Dataset_Scen[[2]]
Scenario3 <- Multiple_Dataset_Scen[[3]]

# Univariate AD model - lm with changepoint
Test_Results_Scenario1_Uni_lm <- Test_lm_MCMC_Results_Uni(Scenario1, tests = "lm")
Test_Results_Scenario2_Uni_lm <- Test_lm_MCMC_Results_Uni(Scenario2, tests = "lm")
Test_Results_Scenario3_Uni_lm <- Test_lm_MCMC_Results_Uni(Scenario3, tests = "lm")

Uni_Final_Scenario_1_lm <- Unlist_function(Summarise_Test(Test_Results_Scenario1_Uni_lm), Percentage_Neut_Amp)
Uni_Final_Scenario_2_lm <- Unlist_function(Summarise_Test(Test_Results_Scenario2_Uni_lm), Percentage_Neut_Amp)
Uni_Final_Scenario_3_lm <- Unlist_function(Summarise_Test(Test_Results_Scenario3_Uni_lm), Percentage_B)

Uni_Final_Scenario_1_lm <- Uni_Final_Scenario_1_lm %>% mutate(Scen = "Scenario 1")
Uni_Final_Scenario_2_lm <- Uni_Final_Scenario_2_lm %>% mutate(Scen = "Scenario 2")
Uni_Final_Scenario_3_lm <- Uni_Final_Scenario_3_lm %>% mutate(Scen = "Scenario 3")
write.table(Uni_Final_Scenario_1_lm, "../../data/Chapter_5/Univariate_7_lm_Model_2_scenario_1.txt", sep = "\t")
write.table(Uni_Final_Scenario_2_lm, "../../data/Chapter_5/Univariate_7_lm_Model_2_scenario_2.txt", sep = "\t")
write.table(Uni_Final_Scenario_3_lm, "../../data/Chapter_5/Univariate_7_lm_Model_2_scenario_3.txt", sep = "\t")

# Univariate AD model - MCMC with changepoint
Test_Results_Scenario1_Uni_MCMC <- Test_lm_MCMC_Results_Uni(Scenario1, tests = "MCMC")
Test_Results_Scenario2_Uni_MCMC <- Test_lm_MCMC_Results_Uni(Scenario2, tests = "MCMC")
Test_Results_Scenario3_Uni_MCMC <- Test_lm_MCMC_Results_Uni(Scenario3, tests = "MCMC")

Uni_Final_Scenario_1_MCMC <- Unlist_function(Summarise_Test(Test_Results_Scenario1_Uni_MCMC), Percentage_Neut_Amp)
Uni_Final_Scenario_2_MCMC <- Unlist_function(Summarise_Test(Test_Results_Scenario2_Uni_MCMC), Percentage_Neut_Amp)
Uni_Final_Scenario_3_MCMC <- Unlist_function(Summarise_Test(Test_Results_Scenario3_Uni_MCMC), Percentage_B)

Uni_Final_Scenario_1_MCMC <- Uni_Final_Scenario_1_MCMC %>% mutate(Scen = "Scenario 1")
Uni_Final_Scenario_2_MCMC <- Uni_Final_Scenario_2_MCMC %>% mutate(Scen = "Scenario 2")
Uni_Final_Scenario_3_MCMC <- Uni_Final_Scenario_3_MCMC %>% mutate(Scen = "Scenario 3")

write.table(Uni_Final_Scenario_1_MCMC, "../../data/Chapter_5/Univariate_7_MCMC_Model_2_scenario_1.txt", sep = "\t")
write.table(Uni_Final_Scenario_2_MCMC, "../../data/Chapter_5/Univariate_7_MCMC_Model_2_scenario_2.txt", sep = "\t")
write.table(Uni_Final_Scenario_3_MCMC, "../../data/Chapter_5/Univariate_7_MCMC_Model_2_scenario_3.txt", sep = "\t")

# Multivariate ADIM model - lm with changepoint
Test_Results_Scenario1_Multi_lm <- Test_lm_MCMC_Results_Multi(Scenario1, tests = "lm")
Test_Results_Scenario2_Multi_lm <- Test_lm_MCMC_Results_Multi(Scenario2, tests = "lm")
Test_Results_Scenario3_Multi_lm <- Test_lm_MCMC_Results_Multi(Scenario3, tests = "lm")

Multi_Final_Scenario_1_lm <- Unlist_function(Summarise_Test(Test_Results_Scenario1_Multi_lm, Multi_lm = T), Percentage_Neut_Amp)
Multi_Final_Scenario_2_lm <- Unlist_function(Summarise_Test(Test_Results_Scenario2_Multi_lm, Multi_lm = T), Percentage_Neut_Amp)
Multi_Final_Scenario_3_lm <- Unlist_function(Summarise_Test(Test_Results_Scenario3_Multi_lm, Multi_lm = T), Percentage_B)

Multi_Final_Scenario_1_lm <- Multi_Final_Scenario_1_lm %>% mutate(Scen = "Scenario 1")
Multi_Final_Scenario_2_lm <- Multi_Final_Scenario_2_lm %>% mutate(Scen = "Scenario 2")
Multi_Final_Scenario_3_lm <- Multi_Final_Scenario_3_lm %>% mutate(Scen = "Scenario 3")

write.table(Multi_Final_Scenario_1_lm, "../../data/Chapter_5/Multivariate_7_lm_Model_2_scenario_1.txt", sep = "\t")
write.table(Multi_Final_Scenario_2_lm, "../../data/Chapter_5/Multivariate_7_lm_Model_2_scenario_2.txt", sep = "\t")
write.table(Multi_Final_Scenario_3_lm, "../../data/Chapter_5/Multivariate_7_lm_Model_2_scenario_3.txt", sep = "\t")

# Multivariate ADIM model - MCMC with changepoint
Test_Results_Scenario1_Multi_MCMC <- Test_lm_MCMC_Results_Multi(Scenario1, tests = "MCMC")
Test_Results_Scenario2_Multi_MCMC <- Test_lm_MCMC_Results_Multi(Scenario2, tests = "MCMC")
Test_Results_Scenario3_Multi_MCMC <- Test_lm_MCMC_Results_Multi(Scenario3, tests = "MCMC")

Multi_Final_Scenario_1_MCMC <- Unlist_function(Summarise_Test(Test_Results_Scenario1_Multi_MCMC), Percentage_Neut_Amp)
Multi_Final_Scenario_2_MCMC <- Unlist_function(Summarise_Test(Test_Results_Scenario2_Multi_MCMC), Percentage_Neut_Amp)
Multi_Final_Scenario_3_MCMC <- Unlist_function(Summarise_Test(Test_Results_Scenario3_Multi_MCMC), Percentage_B)

Multi_Final_Scenario_1_MCMC <- Multi_Final_Scenario_1_MCMC %>% mutate(Scen = "Scenario 1")
Multi_Final_Scenario_2_MCMC <- Multi_Final_Scenario_2_MCMC %>% mutate(Scen = "Scenario 2")
Multi_Final_Scenario_3_MCMC <- Multi_Final_Scenario_3_MCMC %>% mutate(Scen = "Scenario 3")

write.table(Multi_Final_Scenario_1_MCMC, "../../data/Chapter_5/Multivariate_7_MCMC_Model_2_scenario_1.txt", sep = "\t")
write.table(Multi_Final_Scenario_2_MCMC, "../../data/Chapter_5/Multivariate_7_MCMC_Model_2_scenario_2.txt", sep = "\t")
write.table(Multi_Final_Scenario_3_MCMC, "../../data/Chapter_5/Multivariate_7_MCMC_Model_2_scenario_3.txt", sep = "\t")

