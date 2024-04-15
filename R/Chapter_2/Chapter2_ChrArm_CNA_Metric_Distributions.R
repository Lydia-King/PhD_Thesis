# Chapter 2: Summary and Density Plots of Chr Arm Metrics
## Load up necessary libraries
library(gt)
library(tidyverse)
library(reshape2)
library(overlapping)
library(kableExtra)
library(xtable)
library(ggpubr)

## Load up data
PA_CNA_Score_Metrics <-
    read.delim(
        "../../data/Processed_Data/Chr_Arm_CNA_Score.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )
PA_CNA_Burden_Metrics <-
    read.delim(
        "../../data/Processed_Data/Chr_Arm_CNA_Burden.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

## Separate out CNA Score metrics
CNA_Arm_Score_All <-
    PA_CNA_Score_Metrics[, grepl("All" , names(PA_CNA_Score_Metrics)) |
                             grepl("PATIENT_ID" , names(PA_CNA_Score_Metrics))]
CNA_Arm_Score_CCA <-
    PA_CNA_Score_Metrics[, grepl("CCA" , names(PA_CNA_Score_Metrics)) |
                             grepl("PATIENT_ID" , names(PA_CNA_Score_Metrics))]
CNA_Arm_Score_CCA <- na.omit(CNA_Arm_Score_CCA)

## Separate out CNA Burden metrics
CNA_Arm_Burden_All <-
    PA_CNA_Burden_Metrics[, grepl("All" , names(PA_CNA_Burden_Metrics)) |
                              grepl("PATIENT_ID" , names(PA_CNA_Burden_Metrics))]
CNA_Arm_Burden_CCA <-
    PA_CNA_Burden_Metrics[, grepl("CCA" , names(PA_CNA_Burden_Metrics)) |
                              grepl("PATIENT_ID" , names(PA_CNA_Burden_Metrics))]
CNA_Arm_Burden_CCA <- na.omit(CNA_Arm_Burden_CCA)

# Chromosome 1q
## Summary Table (Score)
CNA_Arm_Score_All %>% select_if(grepl("_1q", names(.))) %>%
    dplyr::rename(
        "CNA Score" = CNA_Score_All_1q,
        "CNA Amp Score" = CNA_Amp_All_1q,
        "CNA Del Score" = CNA_Del_All_1q,
        "Difference Score" = CNA_Difference_All_1q,
        "Percentage Amp Score" = CNA_Per_Amp_All_1q,
        "Percentage Del Score" = CNA_Per_Del_All_1q
    ) %>%
    melt() %>%
    dplyr::rename("CNA Score Metric" = variable) %>%
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
    tab_header(title =  md(
        "**Summary Statistics of CNA Score Metrics on Chromosome 1q (All)**"
    )) %>%  tab_stubhead(label = md("**CNA Score Metric**")) %>%
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
                                                columns = `CNA Score Metric`) |>
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% 
    gtsave("CNA_Score_Metric_All_Chr1q_Summary.png",
           path = "../../tables/Chapter_2/")

CNA_Arm_Score_CCA %>% select_if(grepl("_1q", names(.))) %>%
    dplyr::rename(
        "CNA Score" = CNA_Score_CCA_1q,
        "CNA Amp Score" = CNA_Amp_CCA_1q,
        "CNA Del Score" = CNA_Del_CCA_1q,
        "Difference Score" = CNA_Difference_CCA_1q,
        "Percentage Amp Score" = CNA_Per_Amp_CCA_1q,
        "Percentage Del Score" = CNA_Per_Del_CCA_1q
    ) %>%
    melt() %>%
    dplyr::rename("CNA Score Metric" = variable) %>%
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
    tab_header(title =  md("**Summary Statistics of CNA Score Metrics on Chromosome 1q (CC)**")) %>%  
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
               use_seps = TRUE) %>%  cols_align(align = "left",
                                                columns = `CNA Score Metric`) |>
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>%
    gtsave("CNA_Score_Metric_CCA_Chr1q_Summary.png",
           path = "../../tables/Chapter_2/")


## Summary Table (Burden)
CNA_Arm_Burden_All %>% select_if(grepl("_1q", names(.))) %>%
    dplyr::rename(
        "CNA Burden" = CNA_Burden_All_1q,
        "CNA Amp Burden" = CNA_Amp_All_1q,
        "CNA Del Burden" = CNA_Del_All_1q,
        "Difference Burden" = CNA_Difference_All_1q,
        "Percentage Amp Burden" = CNA_Per_Amp_All_1q,
        "Percentage Del Burden" = CNA_Per_Del_All_1q
    ) %>%
    melt() %>%
    dplyr::rename("CNA Burden Metric" = variable) %>%
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
    tab_header(title =  md(
        "**Summary Statistics of CNA Burden Metrics on Chromosome 1q (All)**"
    )) %>%  tab_stubhead(label = md("**CNA Burden Metric**")) %>%
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
    tab_options(table.width = pct(100)) %>% gtsave("CNA_Burden_Metric_All_Chr1q_Summary.png",
                                                   path = "../../tables/Chapter_2/")


CNA_Arm_Burden_CCA %>% select_if(grepl("_1q", names(.))) %>%
    dplyr::rename(
        "CNA Burden" = CNA_Burden_CCA_1q,
        "CNA Amp Burden" = CNA_Amp_CCA_1q,
        "CNA Del Burden" = CNA_Del_CCA_1q,
        "Difference Burden" = CNA_Difference_CCA_1q,
        "Percentage Amp Burden" = CNA_Per_Amp_CCA_1q,
        "Percentage Del Burden" = CNA_Per_Del_CCA_1q
    ) %>%
    melt() %>%
    dplyr::rename("CNA Burden Metric" = variable) %>%
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
    tab_header(title =  md(
        "**Summary Statistics of CNA Burden Metrics on Chromosome 1q (CC)**"
    )) %>%  tab_stubhead(label = md("**CNA Burden Metric**")) %>%
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
    tab_options(table.width = pct(100)) %>% gtsave("CNA_Burden_Metric_CCA_Chr1q_Summary.png",
                                                   path = "../../tables/Chapter_2/")

## Density Plot
### CNA Score Metric Comparative Density Plot
#### Setup
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


# write for loop:
chrom <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q",
           "5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q",
           "9p", "9q", "10p", "10q", "11p", "11q", "12p",
           "12q", "13q", "14q", "15q", "16p", "16q", "17p",
           "17q", "18p", "18q", "19p", "19q", "20p", "20q",
           "21p", "21q", "22q", "Xp", "Xq")

for(i in 1:length(chrom)){
    
    Chr_1q_Score <-
        PA_CNA_Score_Metrics %>% select_if(grepl(paste0("_", chrom[i]), names(.)))
    
    Chr_1q_Burden <-
        PA_CNA_Burden_Metrics %>% select_if(grepl(paste0("_", chrom[i]), names(.)))
    
    colnames(Chr_1q_Score) <- gsub(pattern = paste0("_", chrom[i]), replacement = "", x = colnames(Chr_1q_Score))
    colnames(Chr_1q_Burden) <- gsub(pattern = paste0("_", chrom[i]), replacement = "", x = colnames(Chr_1q_Burden))
    
    Plot_CNA_Score <- Chr_1q_Score %>%
        melt(.) %>%
        mutate(., State = variable) %>%
        mutate(
            .,
            State = recode(
                State,
                "CNA_Score_All" = "Absolute CNA Score",
                "CNA_Score_CCA" = "Absolute CNA Score",
                "CNA_Amp_All" = "CNA Amp Score",
                "CNA_Amp_CCA" = "CNA Amp Score",
                "CNA_Del_All" = "CNA Del Score",
                "CNA_Del_CCA" = "CNA Del Score",
                "CNA_Difference_All" = "Difference",
                "CNA_Difference_CCA" = "Difference",
                "CNA_Per_Amp_All" = "% Amp Score",
                "CNA_Per_Amp_CCA" = "% Amp Score",
                "CNA_Per_Del_All" = "% Del Score",
                "CNA_Per_Del_CCA" = "% Del Score"
            )
        ) %>%
        mutate(., State = factor(
            State,
            levels = c(
                'Absolute CNA Score',
                'Difference',
                'CNA Amp Score',
                'CNA Del Score',
                '% Amp Score',
                "% Del Score"
            )
        )) %>%
        mutate(., variable = factor(
            variable,
            levels = c(
                "CNA_Score_All",
                "CNA_Score_CCA",
                "CNA_Difference_All",
                "CNA_Difference_CCA",
                "CNA_Amp_All",
                "CNA_Amp_CCA",
                "CNA_Del_All",
                "CNA_Del_CCA",
                "CNA_Per_Amp_All",
                "CNA_Per_Amp_CCA",
                "CNA_Per_Del_All",
                "CNA_Per_Del_CCA"
            )
        ))
    
    #### Plot density plots
    Plot_CNA_Score %>% ggplot() +
        geom_histogram(
            aes(x = value, y = ..density.., fill = variable),
            position = "identity",
            alpha = 0.45,
            color = "grey30"
        ) +
        geom_density(aes(x = value, color = variable)) +
        facet_wrap( ~ State, scales = "free", ncol = 2) +
        ggtitle(paste("CNA Score Distribution on Chromosome", chrom[i])) +
        ylab("Density") +
        xlab("") +
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
    
    ggsave(paste0("../../figures/Chapter_2/CNA_Score_Comparative_Density_Chr", chrom[i], ".png"),
           width = 14,
           height = 10
    )
    
    Plot_CNA_Burden <- Chr_1q_Burden %>%
        melt(.) %>%
        mutate(., State = variable) %>%
        mutate(
            .,
            State = recode(
                State,
                "CNA_Burden_All" = "Absolute CNA Burden",
                "CNA_Burden_CCA" = "Absolute CNA Burden",
                "CNA_Amp_All" = "CNA Amp Burden",
                "CNA_Amp_CCA" = "CNA Amp Burden",
                "CNA_Del_All" = "CNA Del Burden",
                "CNA_Del_CCA" = "CNA Del Burden",
                "CNA_Difference_All" = "Difference",
                "CNA_Difference_CCA" = "Difference",
                "CNA_Per_Amp_All" = "% Amp Burden",
                "CNA_Per_Amp_CCA" = "% Amp Burden",
                "CNA_Per_Del_All" = "% Del Burden",
                "CNA_Per_Del_CCA" = "% Del Burden"
            )
        ) %>%
        mutate(., State = factor(
            State,
            levels = c(
                'Absolute CNA Burden',
                'Difference',
                'CNA Amp Burden',
                'CNA Del Burden',
                '% Amp Burden',
                "% Del Burden"
            )
        )) %>%
        mutate(., variable = factor(
            variable,
            levels = c(
                "CNA_Burden_All",
                "CNA_Burden_CCA",
                "CNA_Difference_All",
                "CNA_Difference_CCA",
                "CNA_Amp_All",
                "CNA_Amp_CCA",
                "CNA_Del_All",
                "CNA_Del_CCA",
                "CNA_Per_Amp_All",
                "CNA_Per_Amp_CCA",
                "CNA_Per_Del_All",
                "CNA_Per_Del_CCA"
            )
        ))
    
    #### Plot density plots
    Plot_CNA_Burden %>% ggplot() +
        geom_histogram(
            aes(x = value, y = ..density.., fill = variable),
            position = "identity",
            alpha = 0.45,
            color = "grey30"
        ) +
        geom_density(aes(x = value, color = variable)) +
        facet_wrap( ~ State, scales = "free", ncol = 2) +
        ggtitle(paste("CNA Burden Distribution on Chromosome", chrom[i])) +
        ylab("Density") +
        xlab("") +
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
        guides(fill = guide_legend(nrow = 2, byrow = FALSE))
    
    ggsave(
        paste0("../../figures/Chapter_2/CNA_Burden_Comparative_Density_Chr", chrom[i], ".png"),
        width = 14,
        height = 10
    )
    
}



# Per Chromosome Arm CNA Score and Burden Metric Table of Overlap
## Set up dataframe to store results
Overlap_CNA_Arm_Score <-
    data.frame("Metric" = character(), "Overlap" = numeric())
Overlap_CNA_Arm_Burden <-
    data.frame("Metric" = character(), "Overlap" = numeric())

## Calculate % overlap for CNA Score metrics
for (i in c(2:ncol(CNA_Arm_Score_All))) {
    Score_All <- CNA_Arm_Score_All[, i]
    Score_CCA <- na.omit(CNA_Arm_Score_CCA[, i])
    x <- list(X1 = Score_All , X2 = Score_CCA)
    y <- overlap(x, plot = F, type = "2")
    Overlap_CNA_Arm_Score[i - 1, 1] <- names(CNA_Arm_Score_All[i])
    Overlap_CNA_Arm_Score[i - 1, 2] <- round(y$OV * 100, digits = 2)
}

## Calculate % overlap for CNA Burden metrics
for (i in c(2:ncol(CNA_Arm_Burden_All))) {
    Burden_All <- CNA_Arm_Burden_All[, i]
    Burden_CCA <- na.omit(CNA_Arm_Burden_CCA[, i])
    x <- list(X1 = Burden_All , X2 = Burden_CCA)
    y <- overlap(x, plot = F, type = "2")
    Overlap_CNA_Arm_Burden[i - 1, 1] <- names(CNA_Arm_Burden_All[i])
    Overlap_CNA_Arm_Burden[i - 1, 2] <- round(y$OV * 100, digits = 2)
}

## Order
Overlap_CNA_Arm_Score <-
    Overlap_CNA_Arm_Score[order(Overlap_CNA_Arm_Score$Overlap), ]
Overlap_CNA_Arm_Burden <-
    Overlap_CNA_Arm_Burden[order(Overlap_CNA_Arm_Burden$Overlap), ]


## Change/Format metric names for top 10 entries contained in column 1
Overlap_CNA_Arm_Score[1:10,1] <- c("% CNA Score Amplified 9p", "CNA Amplification Score 9p",
                                   "% CNA Score Deleted 7p", "% CNA Score Amplified 18q",
                                   "Difference 18p", "% CNA Score Deleted 8q",
                                   "% CNA Score Amplified 22q",   "Difference Xq",
                                   "Difference 18q", "Difference 19p")

Overlap_CNA_Arm_Burden[1:10,1] <- c("% CNA Burden Amplified 9p", "CNA Amplification Burden 9p",
                                    "% CNA Burden Deleted 7p", "% CNA Burden Amplified 18q",
                                    "CNA Amplification Burden Xq", "Difference 19p", 
                                    "% CNA Burden Amplified 22q",
                                    "Difference 18p", "Difference Xq", "Difference 18q")

## Remove row names
rownames(Overlap_CNA_Arm_Score) <- NULL
rownames(Overlap_CNA_Arm_Burden) <- NULL

## Save
write.table(
    Overlap_CNA_Arm_Score[1:10,],
    file = "../../data/Processed_Data/CNA_Arm_Score_Metrics_Overlap.txt",
    sep = "\t",
    row.names = FALSE
)
write.table(
    Overlap_CNA_Arm_Burden[1:10,],
    file = "../../data/Processed_Data/CNA_Arm_Burden_Metrics_Overlap.txt",
    sep = "\t",
    row.names = FALSE
)

# Plot Density Plots for 9p and 7p
## Chromosome 9p Comparative Density Plot
p1 <- ggplot() +
    geom_density(aes(
        CNA_Arm_Score_All$CNA_Per_Amp_All_9p,
        fill = "data1",
        color = "data1"
    ),
    alpha = .4) +
    geom_density(aes(
        CNA_Arm_Score_CCA$CNA_Per_Amp_CCA_9p,
        fill = "data2",
        color = "data2"
    ),
    alpha = .4) +
    scale_fill_manual(
        name = "Score Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Amp All", "% Amp CC")
    ) +
    scale_color_manual(
        name = "Score Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Amp All", "% Amp CC")
    ) +
    ggtitle("% Amp Score on Chromosome 9p") +
    ylab("Density") +
    xlab("% Amp Score") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(hjust = 0.5, size = 28),
        axis.title.y = element_text(hjust = 0.5, size = 28),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.title = element_text(
            colour = "black",
            size = 22,
            face = "bold"
        ),
        legend.text = element_text(size = 22),
        legend.position = "top",
        strip.text = element_text(size = 24),
        legend.margin = margin(0, 0, 0, 0)
    ) +
    coord_cartesian(xlim = c(0, 100))

p2 <- ggplot() +
    geom_density(aes(
        CNA_Arm_Burden_All$CNA_Per_Amp_All_9p,
        fill = "data1",
        color = "data1"
    ),
    alpha = .4) +
    geom_density(aes(
        CNA_Arm_Burden_CCA$CNA_Per_Amp_CCA_9p,
        fill = "data2",
        color = "data2"
    ),
    alpha = .4) +
    scale_fill_manual(
        name = "Burden Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Amp All", "% Amp CC")
    ) +
    scale_color_manual(
        name = "Burden Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Amp All", "% Amp CC")
    ) +
    ggtitle("% Amp Burden on Chromosome 9p") +
    ylab("Density") + xlab("% Amp Burden") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(hjust = 0.5, size = 28),
        axis.title.y = element_text(hjust = 0.5, size = 28),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.title = element_text(
            colour = "black",
            size = 22,
            face = "bold"
        ),
        legend.text = element_text(size = 22),
        legend.position = "top",
        strip.text = element_text(size = 24),
        legend.margin = margin(0, 0, 0, 0)
    )


ggsave("../../figures/Chapter_2/CNA_Score_9p.png",
       p1,
       width = 11,
       height = 7.5)
ggsave("../../figures/Chapter_2/CNA_Burden_9p.png",
       p2,
       width = 11,
       height = 7.5)


## Chromosome 7p Comparative Density Plot
p1 <- ggplot() +
    geom_density(aes(
        CNA_Arm_Score_All$CNA_Per_Del_All_7p,
        fill = "data1",
        color = "data1"
    ),
    alpha = .4) +
    geom_density(aes(
        CNA_Arm_Score_CCA$CNA_Per_Del_CCA_7p,
        fill = "data2",
        color = "data2"
    ),
    alpha = .4) +
    scale_fill_manual(
        name = "Score Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Del All", "% Del CC")
    ) +
    scale_color_manual(
        name = "Score Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Del All", "% Del CC")
    ) +
    ggtitle("% Del Score on Chromosome 7p") +
    ylab("Density") +
    xlab("% Del Score") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(hjust = 0.5, size = 28),
        axis.title.y = element_text(hjust = 0.5, size = 28),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.title = element_text(
            colour = "black",
            size = 22,
            face = "bold"
        ),
        legend.text = element_text(size = 22),
        legend.position = "top",
        strip.text = element_text(size = 24),
        legend.margin = margin(0, 0, 0, 0)
    ) +
    coord_cartesian(xlim = c(0, 100))

p2 <- ggplot() +
    geom_density(aes(
        CNA_Arm_Burden_All$CNA_Per_Del_All_7p,
        fill = "data1",
        color = "data1"
    ),
    alpha = .4) +
    geom_density(aes(
        CNA_Arm_Burden_CCA$CNA_Per_Del_CCA_7p,
        fill = "data2",
        color = "data2"
    ),
    alpha = .4) +
    scale_fill_manual(
        name = "Burden Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Del All", "% Del CC")
    ) +
    scale_color_manual(
        name = "Burden Metric:",
        values = c(data1 = "#002d74",
                   data2 = "#942092"),
        labels = c("% Del All", "% Del CC")
    ) +
    ggtitle("% Del Burden on Chromosome 7p") +
    ylab("Density") + xlab("% Del Burden") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 30),
        axis.title.x = element_text(hjust = 0.5, size = 28),
        axis.title.y = element_text(hjust = 0.5, size = 28),
        axis.text.x = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.title = element_text(
            colour = "black",
            size = 22,
            face = "bold"
        ),
        legend.text = element_text(size = 22),
        legend.position = "top",
        strip.text = element_text(size = 24),
        legend.margin = margin(0, 0, 0, 0)
    )


ggsave("../../figures/Chapter_2/CNA_Score_7p.png",
       p1,
       width = 11,
       height = 7.5)
ggsave("../../figures/Chapter_2/CNA_Burden_7p.png",
       p2,
       width = 11,
       height = 7.5)
