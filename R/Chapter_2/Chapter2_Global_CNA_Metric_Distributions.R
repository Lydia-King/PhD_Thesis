# Chapter 2: Summary Stats and Density Plots of Global Metrics
## Load up necessary libraries
library(gt)
library(tidyverse)
library(reshape2)
library(overlapping)
library(kableExtra)
library(xtable)

## Load up data
CNA_Score_Metrics_All <-
    read.delim(
        "../../data/Processed_Data/CNA_Score_All.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )
CNA_Score_Metrics_CCA <-
    read.delim(
        "../../data/Processed_Data/CNA_Score_CCA.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )
CNA_Burden_Metrics_All <-
    read.delim(
        "../../data/Processed_Data/CNA_Burden_All.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )
CNA_Burden_Metrics_CCA <-
    read.delim(
        "../../data/Processed_Data/CNA_Burden_CCA.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

## Tables
### CNA Score (All)
CNA_Score_Metrics_All %>%
    rename(
        "Absolute CNA Score" = CNA_Score_All,
        "CNA Amp Score" = Amp_Score_All,
        "CNA Del Score" = Del_Score_All,
        "Difference Score" = Difference_Score_All,
        "Percentage Score Amp" = Percentage_Score_Amp_All,
        "Percentage Score Del" = Percentage_Score_Del_All
    ) %>%
    melt() %>%
    rename("CNA Score Metric" = variable) %>%
    group_by(`CNA Score Metric`) %>%
    dplyr::summarise(
        "n" = n(),
        "min" = min(value),
        "mean" = mean(value),
        "median" = median(value),
        "max" = max(value),
        "sd"  = sd(value)
    ) %>%
    gt() %>%
    tab_header(title =  md("**Summary Statistics of CNA Score Metrics (All)**")) %>%  
    tab_stubhead(label = md("**CNA Score Metric**")) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5, 6, 7))) %>%
    fmt_number(
        columns = c(3, 4, 5, 6, 7),
        decimals = 2,
        use_seps = TRUE
    ) %>%
    fmt_number(columns = c(2),
               decimals = 0,
               use_seps = TRUE) %>% cols_align(align = "left",
                                               columns = `CNA Score Metric`) |>
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% 
    gtsave("Global_CNA_Score_Metric_All_Summary.png", path = "../../tables/Chapter_2/")


### CNA Score (CCA)
CNA_Score_Metrics_CCA %>%
    rename(
        "Absolute CNA Score" = CNA_Score_CCA,
        "CNA Amp Score" = Amp_Score_CCA,
        "CNA Del Score" = Del_Score_CCA,
        "Difference Score" = Difference_Score_CCA,
        "Percentage Score Amp" = Percentage_Score_Amp_CCA,
        "Percentage Score Del" = Percentage_Score_Del_CCA
    ) %>%
    melt() %>%
    rename("CNA Score Metric" = variable) %>%
    group_by(`CNA Score Metric`) %>%
    dplyr::summarise(
        "n" = n(),
        "min" = min(value),
        "mean" = mean(value),
        "median" = median(value),
        "max" = max(value),
        "sd"  = sd(value)
    ) %>%
    gt() %>%
    tab_header(title =  md("**Summary Statistics of CNA Score Metrics (CC)**")) %>%  
    tab_stubhead(label = md("**CNA Score Metric**")) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5, 6, 7))) %>%
    fmt_number(
        columns = c(3, 4, 5, 6, 7),
        decimals = 2,
        use_seps = TRUE
    ) %>%
    fmt_number(columns = c(2),
               decimals = 0,
               use_seps = TRUE) %>% cols_align(align = "left",
                                               columns = `CNA Score Metric`) |>
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% 
    gtsave("Global_CNA_Score_Metric_CCA_Summary.png", path = "../../tables/Chapter_2/")


### CNA Burden (All)
CNA_Burden_Metrics_All %>%
    rename(
        "CNA Burden" = CNA_Burden_All,
        "CNA Amp Burden" = Amp_Burden_All,
        "CNA Del Burden" = Del_Burden_All,
        "Difference Burden" = Difference_Burden_All,
        "Percentage Burden Amp" = Percentage_Burden_Amp_All,
        "Percentage Burden Del" = Percentage_Burden_Del_All
    ) %>%
    melt() %>%
    rename("CNA Burden Metric" = variable) %>%
    group_by(`CNA Burden Metric`) %>%
    dplyr::summarise(
        "n" = n(),
        "min" = min(value),
        "mean" = mean(value),
        "median" = median(value),
        "max" = max(value),
        "sd"  = sd(value)
    ) %>%
    gt() %>%
    tab_header(title =  md("**Summary Statistics of CNA Burden Metrics (All)**")) %>%  
    tab_stubhead(label = md("**CNA Burden Metric**")) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5, 6, 7))) %>%
    fmt_number(
        columns = c(3, 4, 5, 6, 7),
        decimals = 2,
        use_seps = TRUE
    ) %>%
    fmt_number(columns = c(2),
               decimals = 0,
               use_seps = TRUE) %>% cols_align(align = "left",
                                               columns = `CNA Burden Metric`) |>
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% 
    gtsave("Global_CNA_Burden_Metric_All_Summary.png", path = "../../tables/Chapter_2/")


### CNA Burden (CCA)
CNA_Burden_Metrics_CCA %>%
    rename(
        "CNA Burden" = CNA_Burden_CCA,
        "CNA Amp Burden" = Amp_Burden_CCA,
        "CNA Del Burden" = Del_Burden_CCA,
        "Difference Burden" = Difference_Burden_CCA,
        "Percentage Burden Amp" = Percentage_Burden_Amp_CCA,
        "Percentage Burden Del" = Percentage_Burden_Del_CCA
    ) %>%
    melt() %>%
    rename("CNA Burden Metric" = variable) %>%
    group_by(`CNA Burden Metric`) %>%
    dplyr::summarise(
        "n" = n(),
        "min" = min(value),
        "mean" = mean(value),
        "median" = median(value),
        "max" = max(value),
        "sd"  = sd(value)
    ) %>%
    gt() %>%
    tab_header(title =  md("**Summary Statistics of CNA Burden Metrics (CC)**")) %>%  
    tab_stubhead(label = md("**CNA Burden Metric**")) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    tab_style(style = "padding-top:7px;padding-bottom:7px;padding-left:7px;padding-right:7px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5, 6, 7))) %>%
    fmt_number(
        columns = c(3, 4, 5, 6, 7),
        decimals = 2,
        use_seps = TRUE
    ) %>%
    fmt_number(columns = c(2),
               decimals = 0,
               use_seps = TRUE) %>%  cols_align(align = "left",
                                                columns = `CNA Burden Metric`) |>
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% 
    gtsave("Global_CNA_Burden_Metric_CCA_Summary.png", path = "../../tables/Chapter_2/")


## Density Plots
### CNA Score Metric Comparative Density Plot

pal_2 <-
    c(
        "#1E7900" ,
        "#6FDC6D",
        "#005A49",
        "#76F1BE",
        "#332288",
        "#4EB9EF",
        "#9200D0",
        "#FD7CE0",
        "#9C0000",
        "#F16827",
        "#792C00",
        "#DEBD17"
    )

Plot_CNA_Score <-
    merge(CNA_Score_Metrics_All,
          CNA_Score_Metrics_CCA,
          by = "PATIENT_ID",
          all = TRUE) %>%
    melt(.) %>%
    mutate(., State = variable) %>%
    mutate(
        .,
        State = recode(
            State,
            "CNA_Score_All" = "Absolute CNA Score",
            "CNA_Score_CCA" = "Absolute CNA Score",
            "Amp_Score_All" = "Amplification Score",
            "Amp_Score_CCA" = "Amplification Score",
            "Del_Score_All" = "Deletion Score",
            "Del_Score_CCA" = "Deletion Score",
            "Difference_Score_All" = "Difference",
            "Difference_Score_CCA" = "Difference",
            "Percentage_Score_Amp_All" = "% CNA Score Amplified",
            "Percentage_Score_Amp_CCA" = "% CNA Score Amplified",
            "Percentage_Score_Del_All" = "% CNA Score Deleted",
            "Percentage_Score_Del_CCA" = "% CNA Score Deleted"
        )
    ) %>%
    mutate(., State = factor(
        State,
        levels = c(
            'Absolute CNA Score',
            'Difference',
            'Amplification Score',
            'Deletion Score',
            '% CNA Score Amplified',
            "% CNA Score Deleted"
        )
    )) %>%
    mutate(., variable = factor(
        variable,
        levels = c(
            "CNA_Score_All",
            "CNA_Score_CCA",
            "Difference_Score_All",
            "Difference_Score_CCA",
            "Amp_Score_All",
            "Amp_Score_CCA",
            "Del_Score_All",
            "Del_Score_CCA",
            "Percentage_Score_Amp_All",
            "Percentage_Score_Amp_CCA",
            "Percentage_Score_Del_All",
            "Percentage_Score_Del_CCA"
        )
    ))


Plot_CNA_Score %>% ggplot() +
    geom_histogram(
        aes(x = value, y = ..density.., fill = variable),
        position = "identity",
        alpha = 0.45,
        color = "grey30"
    ) +
    geom_density(aes(x = value, color = variable)) +
    facet_wrap( ~ State, scales = "free", ncol = 2) +
    ggtitle("Global CNA Score Distribution") +
    ylab("Density") +
    xlab("CNA Score Metric") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(hjust = 0.5, size = 20),
        axis.title.y = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(
            colour = "black",
            size = 14,
            face = "bold"
        ),
        legend.text = element_text(size = 14),
        legend.position = "top",
        strip.text = element_text(size = 16)
    ) +
    scale_color_manual(values = pal_2) +
    scale_fill_manual(
        values = pal_2,
        name = "Score Metric:",
        labels = c(
            "Score All",
            "Score CC",
            "Diff All",
            "Diff CC",
            "Amp All",
            "Amp CC",
            "Del All",
            "Del CC",
            "% Amp All",
            "% Amp CC",
            "% Del All",
            "% Del CC"
        ),
        guide = guide_legend()
    ) +
    guides(color = "none") +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE))

ggsave(
    "../../figures/Chapter_2/Global_CNA_Score_Comparative_Density.png",
    width = 14,
    height = 10
)

### CNA Burden Metric Comparative Density Plot

Plot_CNA_Burden <-
    merge(CNA_Burden_Metrics_All,
          CNA_Burden_Metrics_CCA,
          by = "PATIENT_ID",
          all = TRUE) %>%
    melt(.) %>%
    mutate(., State = variable) %>%
    mutate(
        .,
        State = recode(
            State,
            "CNA_Burden_All" = "CNA Burden",
            "CNA_Burden_CCA" = "CNA Burden",
            "Amp_Burden_All" = "Amplification Burden",
            "Amp_Burden_CCA" = "Amplification Burden",
            "Del_Burden_All" = "Deletion Burden",
            "Del_Burden_CCA" = "Deletion Burden",
            "Difference_Burden_All" = "Difference",
            "Difference_Burden_CCA" = "Difference",
            "Percentage_Burden_Amp_All" = "% CNA Burden Amplified",
            "Percentage_Burden_Amp_CCA" = "% CNA Burden Amplified",
            "Percentage_Burden_Del_All" = "% CNA Burden Deleted",
            "Percentage_Burden_Del_CCA" = "% CNA Burden Deleted"
        )
    ) %>%
    mutate(., State = factor(
        State,
        levels = c(
            'CNA Burden',
            'Difference',
            'Amplification Burden',
            'Deletion Burden',
            '% CNA Burden Amplified',
            "% CNA Burden Deleted"
        )
    )) %>%
    mutate(., variable = factor(
        variable,
        levels = c(
            "CNA_Burden_All",
            "CNA_Burden_CCA",
            "Difference_Burden_All",
            "Difference_Burden_CCA",
            "Amp_Burden_All",
            "Amp_Burden_CCA",
            "Del_Burden_All",
            "Del_Burden_CCA",
            "Percentage_Burden_Amp_All",
            "Percentage_Burden_Amp_CCA",
            "Percentage_Burden_Del_All",
            "Percentage_Burden_Del_CCA"
        )
    ))

ggplot(Plot_CNA_Burden, aes(x = value)) +
    geom_histogram(
        aes(x = value, y = ..density.., fill = variable),
        position = "identity",
        alpha = 0.45,
        color = "grey30"
    ) +
    geom_density(aes(x = value, color = variable)) +
    facet_wrap( ~ State, scales = "free", ncol = 2) +
    ggtitle("Global CNA Burden Distribution") +
    ylab("Density") +
    xlab("CNA Burden Metric") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 22),
        axis.title.x = element_text(hjust = 0.5, size = 20),
        axis.title.y = element_text(hjust = 0.5, size = 20),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        legend.title = element_text(
            colour = "black",
            size = 14,
            face = "bold"
        ),
        legend.text = element_text(size = 14),
        legend.position = "top",
        strip.text = element_text(size = 16)
    ) +
    scale_fill_manual(
        values = pal_2,
        name = "Burden Metric:",
        labels = c(
            "Burden All",
            "Burden CC",
            "Diff All",
            "Diff CC",
            "Amp All",
            "Amp CC",
            "Del All",
            "Del CC",
            "% Amp All",
            "% Amp CC",
            "% Del All",
            "% Del CC"
        ),
        guide = guide_legend()
    ) +
    guides(color = "none") +
    guides(fill = guide_legend(nrow = 2, byrow = FALSE)) +
    scale_color_manual(values = pal_2)

ggsave(
    "../../figures/Chapter_2/Global_CNA_Burden_Comparative_Density.png",
    width = 14,
    height = 10
)


## CNA Score and Burden Metric Table of Overlap:
### Set up dataframe to store results
Overlap_CNA_Score <-
    data.frame("Metric" = character(), "Overlap" = numeric())
Overlap_CNA_Burden <-
    data.frame("Metric" = character(), "Overlap" = numeric())

### Calculate % overlap for CNA Score metrics
for (i in c(2:ncol(CNA_Score_Metrics_All))) {
    Score_All <- CNA_Score_Metrics_All[, i]
    Score_CCA <- na.omit(CNA_Score_Metrics_CCA[, i])
    x <- list(X1 = Score_All , X2 = Score_CCA)
    y <- overlap(x, plot = F, type = "2")
    Overlap_CNA_Score[i - 1, 1] <- names(CNA_Score_Metrics_All[i])
    Overlap_CNA_Score[i - 1, 2] <- round(y$OV * 100, digits = 2)
}

### Calculate % overlap for CNA Burden metrics
for (i in c(2:ncol(CNA_Burden_Metrics_All))) {
    Burden_All <- CNA_Burden_Metrics_All[, i]
    Burden_CCA <- na.omit(CNA_Burden_Metrics_CCA[, i])
    x <- list(X1 = Burden_All , X2 = Burden_CCA)
    y <- overlap(x, plot = F, type = "2")
    Overlap_CNA_Burden[i - 1, 1] <- names(CNA_Burden_Metrics_All[i])
    Overlap_CNA_Burden[i - 1, 2] <- round(y$OV * 100, digits = 2)
}

### Order
Overlap_CNA_Score <-
    Overlap_CNA_Score[order(Overlap_CNA_Score$Overlap), ]
Overlap_CNA_Burden <-
    Overlap_CNA_Burden[order(Overlap_CNA_Burden$Overlap), ]

### Change/Format metric names contained in column 1
Overlap_CNA_Score[, 1] <-
    c(
        "Absolute CNA Score",
        "CNA Deletion Score",
        "CNA Amplification Score",
        "% CNA Score Amplified",
        "% CNA Score Deleted",
        "Difference Score"
    )

Overlap_CNA_Burden[, 1] <- c(
    "CNA Burden",
    "CNA Deletion Burden",
    "CNA Amplification Burden",
    "% CNA Burden Amplified",
    "% CNA Burden Deleted",
    "Difference Burden"
)

### Remove row names
rownames(Overlap_CNA_Score) <- NULL
rownames(Overlap_CNA_Burden) <- NULL

### Save
write.table(
    Overlap_CNA_Score,
    file = "../../data/Processed_Data/Global_CNA_Score_Metrics_Overlap.txt",
    sep = "\t",
    row.names = FALSE
)
write.table(
    Overlap_CNA_Burden,
    file = "../../data/Processed_Data/Global_CNA_Burden_Metrics_Overlap.txt",
    sep = "\t",
    row.names = FALSE
)
