# Chapter 5: Example Dataset Summary Tables

## Load up libraries
library(tidyverse)
library(gt)

## Load Up Data
DF3 <-
    read.delim("../../data/Chapter_5/Example_Data_1_NoNeut.txt",
               sep = "\t")


DF3_Neut <-
    read.delim("../../data/Chapter_5/Example_Data_1_Neut.txt",
               sep = "\t")

## Tables
### Structure of Data
DF3 %>% filter(Sample %in% c("Sample 1", "Sample 5", "Sample 13")) %>%
    select(Sample, Category, Allele, CP, TS, TE) |>
    mutate(Category = as.character(Category)) |>
    rename("Changepoint" = CP) |>
    gt()  |>
    tab_header(title =  md("**Example of Possible Simulated Samples**")) |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:12px;padding-right:12px",
              locations = cells_body()) |>
    fmt_number(columns = c(3, 4, 5, 6),
               decimals = 0,
               use_seps = TRUE) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5, 6))) |>
    opt_table_outline() |> tab_options(table.width = pct(38)) %>% gtsave(
        "Indv_Simulated_Example_NoNeut.png",
        path = "../../tables/Chapter_5",
        vwidth = 1900,
        vheight = 1000
    )

DF3_Neut %>% filter(Sample %in% c("Sample 1", "Sample 5", "Sample 13")) %>%
    select(Sample, Category, Allele, CP, TS, TE) |>
    mutate(CP = as.numeric(ifelse(is.na(CP), "-", CP))) |>
    rename("Changepoint" = CP) |>
    gt()  |>
    tab_header(title =  md("**Example of Possible Simulated Samples**")) |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:12px;padding-right:12px",
              locations = cells_body()) |>
    fmt_number(columns = c(3, 4, 5, 6),
               decimals = 0,
               use_seps = TRUE) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5, 6))) |>
    opt_table_outline() |> tab_options(table.width = pct(38)) %>% gtsave(
        "Indv_Simulated_Example_Neut.png",
        path = "../../tables/Chapter_5",
        vwidth = 1900,
        vheight = 1000
    )

## Summary Statistics
DF3 |>
    mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Del/Neut", "Amp/Del"))) |>
    dplyr::group_by(Category) |>
    dplyr::summarise(
        "n" = n(),
        "mean" = mean(TS),
        "median" = median(TS),
        "sd"  = sd(TS),
        "mean " = mean(TE),
        "median " = median(TE),
        "sd " = sd(TE)
    ) |>
    gt() |>
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_header(title =  md("**Summary Statistics for Simulated Dataset by Category**")) |>
    tab_spanner(label = md("**TS**"),
                columns = c(3, 4, 5)) |>
    tab_spanner(label = md("**TE**"),
                columns = c(6, 7, 8)) |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1))) |>
    fmt_number(
        columns = c(3, 4, 5, 6, 7, 8),
        decimals = 2,
        use_seps = TRUE
    ) |> fmt_number(columns = c(2),
                    decimals = 0,
                    use_seps = TRUE) |>
    opt_table_outline() |> tab_options(table.width = pct(60)) %>% gtsave(
        "Indv_Simulated_Example_Summary_Table_NoNeut.png",
        path = "../../tables/Chapter_5",
        vwidth = 1900,
        vheight = 1000
    )

DF3_Neut |>
    mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Del/Neut", "Amp/Del"))) |>
    dplyr::group_by(Category) |>
    dplyr::summarise(
        "n" = n(),
        "mean" = mean(TS),
        "median" = median(TS),
        "sd"  = sd(TS),
        "mean " = mean(TE),
        "median " = median(TE),
        "sd " = sd(TE)
    ) |>
    gt() |>
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_header(title =  md("**Summary Statistics for Simulated Dataset by Category**")) |>
    tab_spanner(label = md("**TS**"),
                columns = c(3, 4, 5)) |>
    tab_spanner(label = md("**TE**"),
                columns = c(6, 7, 8)) |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1))) |>
    fmt_number(
        columns = c(3, 4, 5, 6, 7, 8),
        decimals = 2,
        use_seps = TRUE
    ) |> fmt_number(columns = c(2),
                    decimals = 0,
                    use_seps = TRUE) |>
    opt_table_outline() |> tab_options(table.width = pct(60)) %>% gtsave(
        "Indv_Simulated_Example_Summary_Table_Neut.png",
        path = "../../tables/Chapter_5",
        vwidth = 1900,
        vheight = 1000
    )

## Summary Stats by Allele
DF3 |>
    mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Del/Neut", "Amp/Del"))) |>
    dplyr::group_by(Allele, Category) |>
    dplyr::summarise(
        "n" = n(),
        "mean" = mean(TS),
        "median" = median(TS),
        "sd"  = sd(TS),
        "mean " = mean(TE),
        "median " = median(TE),
        "sd " = sd(TE)
    ) |>
    gt() |>
    cols_align(
        align = "left",
        columns = Category
    ) |>
    tab_header(title =  md("**Summary Statistics for Simulated Dataset by Allele and Category**")) |>
    tab_spanner(label = md("**TS**"),
                columns = c(4,5,6)) |>
    tab_spanner(label = md("**TE**"),
                columns = c(7,8,9)) |>
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) |>
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1))) |>
    tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels(columns=c(2))) |>
    fmt_number(
        columns = c(4,5,6,7,8,9),
        decimals = 2,
        use_seps = TRUE
    ) |>
    fmt_number(
        columns = c(2,3),
        decimals = 0,
        use_seps = TRUE
    ) |>
    opt_table_outline() |>
    cols_label(
        Category = md("")) |>
    tab_row_group(
        label = md("**Category (Major Allele)**"),
        rows = c(1:4)
    ) |> 
    tab_row_group(
        label = md("**Category (Minor Allele)**"),
        rows = c(5:7)
    ) |>
    opt_table_outline() |> tab_options(table.width = pct(60)) %>% gtsave(
        "Indv_Simulated_Example_Allele_Summary_Table.png",
        path = "../../tables/Chapter_5",
        vwidth = 1900,
        vheight = 1000
    )