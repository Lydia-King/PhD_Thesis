# Chapter 5: Clean ASCAT Data

## Load up Libraries
library(tidyverse)
library(gt)

## Load up data
ASCAT3_Data_1 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_1.segments.txt", sep="\t")
ASCAT3_Data_2 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_2.segments.txt", sep="\t")
ASCAT3_Data_3 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_3.segments.txt", sep="\t")
ASCAT3_Data_4 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_4.segments.txt", sep="\t")
ASCAT3_Data_5 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_5.segments.txt", sep="\t")
ASCAT3_Data_6 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_6.segments.txt", sep="\t")
ASCAT3_Data_7 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_7.segments.txt", sep="\t")
ASCAT3_Data_8 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_8.segments.txt", sep="\t")

Total_ASCAT3 <- rbind.data.frame(ASCAT3_Data_1, ASCAT3_Data_2, 
                                 ASCAT3_Data_3, ASCAT3_Data_4, 
                                 ASCAT3_Data_5, ASCAT3_Data_6, 
                                 ASCAT3_Data_7, ASCAT3_Data_8)

rm(list=setdiff(ls(), c("Total_ASCAT4", "Total_ASCAT3")))

## Reformat Changepoints
ASCAT_Data <- Total_ASCAT3 %>% mutate(nMajor_Record = nMajor) %>%
    mutate(nMinor_Record = nMinor) %>% mutate(nMajor = ifelse(nMajor > 2, 2, nMajor)) %>%
    mutate(nMinor = ifelse(nMinor > 2, 2, nMinor))

Start <- c()
End <- c()

for (i in c(1:22, "X")) {
    Start <-
        c(Start, min(Total_ASCAT3 %>% filter(chr == i) %>% select(startpos)))
    End <-
        c(End, max(Total_ASCAT3 %>% filter(chr == i) %>% select(endpos)))
}

Reformat_ASCAT_Output <-
    function(Data, Allele = c("nMajor", "nMinor")) {
        Data <- Data %>% select(sample, chr, startpos, endpos,!!Allele)
        unique_samples <- unique(Data$sample)
        unique_chromosomes <- unique(Data$chr)
        
        
        Output_ASCAT <-
            data.frame(
                "Sample" = character(),
                "Chr" = character(),
                "Changepoint" = numeric(),
                "TS" = numeric(),
                "TE" = numeric(),
                "Category" = character()
            )
        
        for (sample_ID in unique_samples) {
            Sample_Data <- Data %>% filter(sample == sample_ID)
            for (chromosome in unique_chromosomes) {
                Chromosome_Data <- Sample_Data %>% filter(chr == chromosome)
                n_rows <- nrow(Chromosome_Data)
                
                if (n_rows > 1 &&
                    length(unique(Chromosome_Data[, 5])) != 1) {
                    for (Row in 2:(n_rows)) {
                        if (Chromosome_Data[Row - 1, 5] != Chromosome_Data[Row, 5]) {
                            if (isTRUE(Chromosome_Data[Row - 1, 5] == 1 &&
                                       Chromosome_Data[Row, 5] == 2)) {
                                Type_Cate <- "Neut/Amp"
                            } else if (isTRUE(Chromosome_Data[Row - 1, 5] == 1 &&
                                              Chromosome_Data[Row, 5] == 0)) {
                                Type_Cate <- "Neut/Del"
                            } else if (isTRUE(Chromosome_Data[Row - 1, 5] == 0 &&
                                              Chromosome_Data[Row, 5] == 2)) {
                                Type_Cate <- "Del/Amp"
                            } else if (isTRUE(Chromosome_Data[Row - 1, 5] == 0 &&
                                              Chromosome_Data[Row, 5] == 1)) {
                                Type_Cate <- "Del/Neut"
                            } else if (isTRUE(Chromosome_Data[Row - 1, 5] == 2 &&
                                              Chromosome_Data[Row, 5] == 0)) {
                                Type_Cate <- "Amp/Del"
                            } else if (isTRUE(Chromosome_Data[Row - 1, 5] == 2 &&
                                              Chromosome_Data[Row, 5] == 1)) {
                                Type_Cate <- "Amp/Neut"
                            }
                            
                            New_df <-
                                data.frame(
                                    "Sample" = sample_ID,
                                    "Chr" = chromosome,
                                    "Changepoint" = (Chromosome_Data[Row -
                                                                         1, 4] + Chromosome_Data[Row, 3]) / 2,
                                    "TS" = NA,
                                    "TE" = NA,
                                    "Category" = Type_Cate
                                )
                            Output_ASCAT <-
                                bind_rows(Output_ASCAT, New_df)
                        }
                    }
                    
                } else {
                    Type_Cate <- "NoChangepoint"
                    New_df <-
                        data.frame(
                            "Sample" = sample_ID,
                            "Chr" = chromosome,
                            "Changepoint" = NA,
                            "TS" = 0,
                            "TE" = 0,
                            "Category" = Type_Cate
                        )
                    Output_ASCAT <- bind_rows(Output_ASCAT, New_df)
                }
            }
        }
        
        Data <- Output_ASCAT
        Output_ASCAT <-
            data.frame(
                "Sample" = character(),
                "Chr" = character(),
                "Changepoint" = numeric(),
                "TS" = numeric(),
                "TE" = numeric(),
                "Category" = character()
            )
        
        for (sample_ID in unique_samples) {
            Sample_Data <- Data %>% filter(Sample == sample_ID)
            chr_num <- 1
            for (chromosome in unique_chromosomes) {
                Chromosome_Data <- Sample_Data %>% filter(Chr == chromosome)
                n_rows <- nrow(Chromosome_Data)
                # Take each chromosome
                if (n_rows > 1) {
                    for (Row in 1:n_rows) {
                        if (Row == 1) {
                            Chromosome_Data[1, "TS"] = as.numeric(Chromosome_Data[1, "Changepoint"]) - Start[chr_num]
                            Chromosome_Data[1, "TE"] =  as.numeric(Chromosome_Data[2, "Changepoint"]) - as.numeric(Chromosome_Data[1, "Changepoint"])
                        } else if (Row == n_rows) {
                            Chromosome_Data[n_rows, "TS"] =  as.numeric(Chromosome_Data[n_rows, "Changepoint"]) - as.numeric(Chromosome_Data[(n_rows -
                                                                                                                                                  1), "Changepoint"])
                            Chromosome_Data[n_rows, "TE"] = End[chr_num] - as.numeric(Chromosome_Data[n_rows, "Changepoint"])
                        } else {
                            Chromosome_Data[Row, "TS"] =  as.numeric(Chromosome_Data[Row, "Changepoint"]) - as.numeric(Chromosome_Data[(Row -
                                                                                                                                            1), "Changepoint"])
                            Chromosome_Data[Row, "TE"] =  as.numeric(Chromosome_Data[Row + 1, "Changepoint"]) - as.numeric(Chromosome_Data[Row, "Changepoint"])
                        }
                    }
                } else if (n_rows == 1 &
                           Chromosome_Data$Category != "NoChangepoint") {
                    Chromosome_Data[1, "TS"] = as.numeric(Chromosome_Data[1, "Changepoint"]) - Start[chr_num]
                    Chromosome_Data[1, "TE"] = End[chr_num] - as.numeric(Chromosome_Data[n_rows, "Changepoint"])
                } else {
                    Chromosome_Data[1, "TS"] = 0
                    Chromosome_Data[1, "TE"] =  0
                }
                
                Output_ASCAT <-
                    rbind.data.frame(Output_ASCAT,  Chromosome_Data)
                chr_num <- chr_num + 1
            }
        }
        return(Output_ASCAT)
    }

## Save Data (writetable)
data_out_Major <- Reformat_ASCAT_Output(ASCAT_Data, Allele = "nMajor")
data_out_Minor <- Reformat_ASCAT_Output(ASCAT_Data, Allele = "nMinor")

data_out_Major_NoNeut <- data_out_Major %>% mutate(TS = ifelse(Category %in% c("Neut/Amp", "Neut/Del"), 0, TS)) %>% 
    mutate(TE = ifelse(Category %in% c("Amp/Neut", "Del/Neut"), 0, TE)) %>% mutate(Allele = "Major")

data_out_Minor_NoNeut <- data_out_Minor %>% mutate(TS = ifelse(Category %in% c("Neut/Amp", "Neut/Del"), 0, TS)) %>% 
    mutate(TE = ifelse(Category %in% c("Amp/Neut", "Del/Neut"), 0, TE)) %>% mutate(Allele = "Minor")

data_out_Major_Neut <- data_out_Major
data_out_Minor_Neut <- data_out_Minor

Total_ASCAT_Data_NoNeut <- rbind.data.frame(data_out_Major_NoNeut, data_out_Minor_NoNeut)
Total_ASCAT_Data_Neut <- rbind.data.frame(data_out_Major_Neut %>% 
                                              mutate(Allele = "Major"), data_out_Minor_Neut  %>% 
                                              mutate(Allele = "Minor"))

## Write Out
write.table(Total_ASCAT_Data_NoNeut, "../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")
write.table(Total_ASCAT_Data_Neut, "../../data/Chapter_5/Total_ASCAT_Data_3Step_Neut.txt", sep="\t")

Total_ASCAT_Data_NoNeut <- Total_ASCAT_Data_NoNeut %>% mutate(TS = TS/1000, TE = TE/1000)
Total_ASCAT_Data_Neut <- Total_ASCAT_Data_Neut %>% mutate(TS = TS/1000, TE = TE/1000)

write.table(Total_ASCAT_Data_NoNeut, "../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")
write.table(Total_ASCAT_Data_Neut, "../../data/Chapter_5/Total_ASCAT_Data_3Step_Neut.txt", sep="\t")

## Create Summary Statistic Tables
Total_ASCAT_Data_NoNeut |>
    mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp",
                                                  "Neut/Del", "Amp/Neut", 
                                                  "Del/Neut", "Amp/Del", "Del/Amp"))) |>
    dplyr::group_by(Category) |> 
    dplyr::summarise("n" = n(), "mean" = mean(TS), "median" = median(TS), 
                     "sd"  = sd(TS), "mean " = mean(TE), "median " = median(TE), 
                     "sd " = sd(TE)) |>
    gt() |>
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_header(
        title =  md("**Summary Statistics of Segment Kilobase Lengths**")
    ) |> 
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_spanner(
        label = md("**TS**"),
        columns = c(3,4,5)
    ) |> 
    tab_spanner(
        label = md("**TE**"),
        columns = c(6,7,8)
    ) |> 
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
        style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
        locations = cells_body()
    ) |>
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels(columns=c(1))) |>
    fmt_number(
        columns = c(3,4,5,6,7,8),
        decimals = 2,
        use_seps = TRUE
    ) |> fmt_number(
        columns = c(2),
        decimals = 0,
        use_seps = TRUE
    ) |>
    opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("ASCAT_Summary_Table_All_Zero.png", 
           path = "../../tables/Chapter_5/", vwidth = 1900, vheight = 1000)

Total_ASCAT_Data_Neut |>
    mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp",
                                                  "Neut/Del", "Amp/Neut", 
                                                  "Del/Neut", "Amp/Del", "Del/Amp"))) |>
    dplyr::group_by(Category) |> 
    dplyr::summarise("n" = n(), "mean" = mean(TS), "median" = median(TS), 
                     "sd"  = sd(TS), "mean " = mean(TE), "median " = median(TE), 
                     "sd " = sd(TE)) |>
    gt() |>
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_header(
        title =  md("**Summary Statistics of Segment Kilobase Lengths**")
    ) |> 
    tab_spanner(
        label = md("**TS**"),
        columns = c(3,4,5)
    ) |> 
    tab_spanner(
        label = md("**TE**"),
        columns = c(6,7,8)
    ) |> 
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
        style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
        locations = cells_body()
    ) |>
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels(columns=c(1))) |>
    fmt_number(
        columns = c(3,4,5,6,7,8),
        decimals = 2,
        use_seps = TRUE
    ) |> fmt_number(
        columns = c(2),
        decimals = 0,
        use_seps = TRUE
    ) |>
    opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("ASCAT_Summary_Table_All_NoZero.png", 
           path = "../../tables/Chapter_5/", vwidth = 1900, vheight = 1000)