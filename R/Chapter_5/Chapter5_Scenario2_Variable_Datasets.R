# Chapter 5: Scenario 2 Variable Datasets

## Source functions 
source("Chapter5_Simulation_Study_Functions.R")

## Load up Libraries 
library(gt)
library(tidyverse)

## Load up datasets 
Multiple_Dataset_Scen <- readRDS(file = "../../data/Chapter_5/Simulation_Study_Data.rds")

## Scenario 2 
Scenario2 <- Multiple_Dataset_Scen[[2]]
Scenario2_Sub_1 <- list("Percentage_90" = list("Samplesize_20" = Scenario2[[9]][[1]]))
Scenario2_Sub_2 <- list("Percentage_10" = list("Samplesize_20" = Scenario2[[1]][[1]]))

## Run Tests again
Test_Results_Scenario2_Uni_lm_1 <- Test_lm_MCMC_Results_Uni(Scenario2_Sub_1, tests = "lm")
Test_Results_Scenario2_Uni_lm_2 <- Test_lm_MCMC_Results_Uni(Scenario2_Sub_2, tests = "lm")

Uni_lm_7_1 <- bind_rows(Test_Results_Scenario2_Uni_lm_1[[1]][[1]], .id = "id")
Uni_lm_7_2 <- bind_rows(Test_Results_Scenario2_Uni_lm_2[[1]][[1]], .id = "id")
Uni_lm_7_1  <- Uni_lm_7_1 %>% filter(Allele == "Major", Category == "Amp/Neut", Direction == "TS")
Uni_lm_7_2  <- Uni_lm_7_2 %>% filter(Allele == "Major", Category == "Amp/Del", Direction == "TE")

## Get Table 
Uni_lm_7_1 |>
    select(-c(Chr)) |> select(id, Category, Allele, Direction, Samplesize, Pred, LB, UB) |> 
    rename(n = Samplesize, Dir = Direction, Dataset = id, Fit = Pred) |>
    gt() |>
    tab_header(
        title =  md("**(A) Univariate TS Parameter Estimates**")
    )  |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    tab_style(
        style = "padding-top:8px;padding-bottom:8px;padding-left:4px;padding-right:4px",
        locations = cells_body()
    ) |>
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
    fmt_scientific(
        columns = c(),
        decimals = 3
    ) |>
    opt_table_outline() %>% gtsave("Model_Scenario_2_Closer_Look_1.png", 
                                   path = "../../tables/Chapter_5/", vwidth = 1900, vheight = 1000)

Uni_lm_7_2 |>
    select(-c(Chr)) |> select(id, Category, Allele, Direction, Samplesize, Pred, LB, UB) |> 
    rename(n = Samplesize, Dir = Direction, Dataset = id, Fit = Pred) |>
    gt() |>
    tab_header(
        title =  md("**(B) Univariate TE Parameter Estimates**")
    )  |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    tab_style(
        style = "padding-top:8px;padding-bottom:8px;padding-left:4px;padding-right:4px",
        locations = cells_body()
    ) |>
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
    fmt_scientific(
        columns = c(),
        decimals = 3
    ) |>
    opt_table_outline() %>% gtsave("Model_Scenario_2_Closer_Look_2.png", 
                                   path = "../../tables/Chapter_5/", vwidth = 1900, vheight = 1000)


sum(Uni_lm_7_1$LB <= 0)/20
sum(Uni_lm_7_2$LB <= 0)/20
