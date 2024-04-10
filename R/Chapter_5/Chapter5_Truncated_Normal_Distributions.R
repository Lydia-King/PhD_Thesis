# Chapter 5: Truncated Normal Distribution Tables 

## Load up Libraries 
library("kableExtra")
library(gt)
library(tidyverse)

## Create Table 
text_tbl <- data.frame(
    Scenarios = c(
        "Scenario 1",
        "Scenario 1",
        "Scenario 2",
        "Scenario 2",
        "Scenario 3",
        "Scenario 3",
        "Scenario 3"
    ),
    Profiles = c(
        "Profile A",
        "Profile B",
        "Profile B",
        "Profile C",
        "Profile A",
        "Profile C",
        "Profile D"
    ),
    Major_Allele = c(
        "No Breakpoint = 0",
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 15,249, σ = 28,815, a = 1, b = 250,000)",
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 15,249, σ = 28,815, a = 1, b = 250,000)",
        
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 22,777, σ = 35,235, a = 1, b = 250,000)
            <br> Del ~ TN(μ = 9,769, σ = 19,739, a = 1, b = 250,000)",
        
        "No Breakpoint = 0",
        
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 22,777, σ = 35,235, a = 1, b = 250,000)
            <br> Del ~ TN(μ = 9,769, σ = 19,739, a = 1, b = 250,000)",
        
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 68,331, σ = 35,235, a = 1, b = 250,000)
            <br> Del ~ TN(μ = 29,307, σ = 19,739, a = 1, b = 250,000)"
    ),
    Minor_Allele = c(
        "No Breakpoint = 0",
        "No Breakpoint = 0 ",
        "No Breakpoint = 0",
        "No Breakpoint = 0",
        
        "No Breakpoint = 0",
        "No Breakpoint = 0",
        "Neutral ~ TN(μ = 31,129, σ = 38,125, a = 1, b = 250,000) <br> Del ~ TN(μ = 8,997, σ = 18,675, a = 1, b = 250,000)"
    ),
    Properties = c(
        "P = 10%, 20%, ..., 90%",
        "n = 20, 50, 80, 100, 200, 500",
        "P = 10%, 20%, ..., 90%",
        "n = 20, 50, 80, 100, 200, 500",
        "P = 20%",
        "P<sub>1</sub>= 10%, 20%, .., 70% <br> P<sub>2</sub> = 100% - P - P<sub>1</sub>",
        "n = 20, 50, 80, 100, 200, 500"
    )
)

## Save Table 
text_tbl %>%
    gt(groupname_col = "Scenarios") %>% fmt_markdown(columns = everything()) %>%
    tab_options(
        table.font.color.light = "#000000",
        column_labels.font.size = 20,
        table.font.size = 18,
        heading.title.font.size = 24
    ) |> tab_header(title =  md("**Scenario Distribution Parameters and Properties**"))  |>
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
              locations = cells_body())  |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5))) |>
    cols_align(align = c("center"),
               columns = everything()) |>  cols_label(Profiles = "",
                                                      Major_Allele = "Major Allele",
                                                      Minor_Allele = "Minor Allele") |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_row_groups())  |>
    opt_table_outline() |>
    tab_options(table.width = pct(70)) %>% gtsave(
        "TN_Distribution.png",
        path = "../../tables/Chapter_5/",
        vwidth = 2280,
        vheight = 1000
    )


## Create Table (Smaller)
text_tbl_2 <- data.frame(
    Profiles = c(
        "**Profile A**",
        "**Profile C**",
        "**Profile D**"
    ),
    Major_Allele = c(
        "No Breakpoint = 0",
        
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 22,777, σ = 35,235, a = 1, b = 250,000)
            <br> Del ~ TN(μ = 9,769, σ = 19,739, a = 1, b = 250,000)",
        
        "Neutral ~ TN(μ = 27,641, σ = 35,854, a = 1, b = 250,000) <br> Amp ~ TN(μ = 68,331, σ = 35,235, a = 1, b = 250,000)
            <br> Del ~ TN(μ = 29,307, σ = 19,739, a = 1, b = 250,000)"
    ),
    Minor_Allele = c(
        "No Breakpoint = 0",
        "No Breakpoint = 0",
        "Neutral ~ TN(μ = 31,129, σ = 38,125, a = 1, b = 250,000) <br> Del ~ TN(μ = 8,997, σ = 18,675, a = 1, b = 250,000)"
    ),
    Properties = c(
        "P<sub>A</sub>= 20%",
        "P<sub>C</sub> = 40%",
        "P<sub>D</sub> = 40%" 
    )
)

## Save Table 
text_tbl_2 %>%
    gt(groupname_col = "Profile") %>% fmt_markdown(columns = everything()) %>%
    tab_options(
        table.font.color.light = "#000000",
        column_labels.font.size = 20,
        table.font.size = 18,
        heading.title.font.size = 24
    ) |> tab_header(title =  md("**Truncated Normal Distribution Parameters and Properties**"))  |>
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
              locations = cells_body())  |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4))) |>
    cols_align(align = c("center"),
               columns = everything()) |>  cols_label(Profiles = "",
                                                      Major_Allele = "Major Allele",
                                                      Minor_Allele = "Minor Allele") |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_row_groups())  |>
    opt_table_outline() |>
    tab_options(table.width = pct(63)) %>% gtsave(
        "TN_Distribution_2.png",
        path = "../../tables/Chapter_5/",
        vwidth = 2280,
        vheight = 1000
    )

