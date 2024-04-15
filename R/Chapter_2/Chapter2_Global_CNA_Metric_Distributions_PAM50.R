# Chapter 2: Distribution of CNA Scores across PAM50 Subtypes

library(tidyverse)
library(DescTools)
library(reshape2)
library(ggpubr)
library(dunn.test)
library(gt)

## Load up data
CNA_Score_Metrics_All <-
    read.delim("../../data/Processed_Data/CNA_Score_All.txt", sep="\t", na.strings=c(""," ","NA"))
CNA_Burden_Metrics_All <-
    read.delim("../../data/Processed_Data/CNA_Burden_All.txt", sep="\t", na.strings=c(""," ","NA"))


Clinical <-
    read.delim(
        "../../data/METABRIC_2021/data_clinical_patient.txt",
        sep = "\t",
        na.strings = c("", " ", "NA"),
        skip = 4
    )

Clinical$PATIENT_ID <- gsub("\\-", ".", Clinical$PATIENT_ID)
Clinical <- Clinical %>%
    mutate(CLAUDIN_SUBTYPE = case_when(CLAUDIN_SUBTYPE  == "Basal" ~ "Basal",
                                       CLAUDIN_SUBTYPE  == "claudin-low" ~ "Claudin-low",
                                       CLAUDIN_SUBTYPE  == "Her2" ~ "HER2",
                                       CLAUDIN_SUBTYPE  == "LumA" ~ "LumA",
                                       CLAUDIN_SUBTYPE  == "LumB" ~ "LumB",
                                       CLAUDIN_SUBTYPE  == "Normal" ~ "Normal",
                                       CLAUDIN_SUBTYPE  == "NC" ~ NA) 
    )

Clinical$INTCLUST <-
    factor(Clinical$INTCLUST,
           levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10"))

CNA_Score_Metrics_With_MC <-
    merge(CNA_Score_Metrics_All, Clinical[, c("PATIENT_ID", "CLAUDIN_SUBTYPE", "INTCLUST")], by =
              "PATIENT_ID")
CNA_Score_Metrics_With_MC <- CNA_Score_Metrics_With_MC %>%
    mutate(Del_Score_All = abs(Del_Score_All))
CNA_Score_Metrics_With_MC_Melted <-
    melt(CNA_Score_Metrics_With_MC) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Score_All" = "Absolute CNA Score",
            "Amp_Score_All" = "CNA Amp Score",
            "Del_Score_All" = "CNA Del Score",
            "Difference_Score_All" = "Difference",
            "Percentage_Score_Amp_All" = "% Amp Score",
            "Percentage_Score_Del_All" = "% Del Score"
        )
    )

ggplot(
    data = CNA_Score_Metrics_With_MC_Melted,
    mapping = aes(x = CLAUDIN_SUBTYPE, y = value, color = CLAUDIN_SUBTYPE)
) +
    geom_boxplot(varwidth = FALSE,
                 aes(fill = CLAUDIN_SUBTYPE),
                 alpha = 0.3) +
    facet_wrap( ~ variable, scales = "free", ncol = 2) +
    geom_jitter(alpha = 0.15) +
    ggtitle("Global CNA Score Metrics by PAM50 Subtype") +
    ylab("CNA Score Metrics") +
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
        strip.text = element_text(size = 16),
        panel.spacing = unit(1, "lines")
    ) +
    scale_x_discrete(na.translate = TRUE) +
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Score_Metrics_With_MC_Melted$CLAUDIN_SUBTYPE
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Score_Metrics_With_MC %>%
                                        group_by(CLAUDIN_SUBTYPE) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Score_Metrics_With_MC_Melted$CLAUDIN_SUBTYPE
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Score_Metrics_With_MC %>%
                                          group_by(CLAUDIN_SUBTYPE) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) +
    stat_kruskal_test(p.adjust.method = "BH", label = "Kruskal-Wallis, italic(p) = {p.adj.format}")

ggsave(
    "../../figures/Chapter_2/Global_CNA_Score_Metrics_Across_PAM50.png",
    width = 15,
    height = 13,
    dpi = 120
)

# Distribution of CNA Burden across PAM50 Subtypes
CNA_Burden_Metrics_With_MC <-
    merge(CNA_Burden_Metrics_All, Clinical[, c("PATIENT_ID", "CLAUDIN_SUBTYPE", "INTCLUST")], by =
              "PATIENT_ID")
CNA_Burden_Metrics_With_MC_Melted <-
    melt(CNA_Burden_Metrics_With_MC) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Burden_All" = "CNA Burden",
            "Amp_Burden_All" = "CNA Amp Burden",
            "Del_Burden_All" = "CNA Del Burden",
            "Difference_Burden_All" = "Difference",
            "Percentage_Burden_Amp_All" = "% Amp Burden",
            "Percentage_Burden_Del_All" = "% Del Burden"
        )
    )

ggplot(
    data = CNA_Burden_Metrics_With_MC_Melted,
    mapping = aes(x = CLAUDIN_SUBTYPE, y = value, color = CLAUDIN_SUBTYPE)
) +
    geom_boxplot(varwidth = FALSE,
                 aes(fill = CLAUDIN_SUBTYPE),
                 alpha = 0.3) +
    facet_wrap( ~ variable, scales = "free", ncol = 2) +
    geom_jitter(alpha = 0.15) +
    ggtitle("Global CNA Burden Metrics by PAM50 Subtype") +
    ylab("CNA Burden Metrics") +
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
        strip.text = element_text(size = 16),
        panel.spacing = unit(1, "lines")
    ) +
    scale_x_discrete(na.translate = TRUE) +
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Burden_Metrics_With_MC_Melted$CLAUDIN_SUBTYPE
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Burden_Metrics_With_MC %>%
                                        group_by(CLAUDIN_SUBTYPE) %>%
                                        summarise(number = n()) %>% select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Burden_Metrics_With_MC_Melted$CLAUDIN_SUBTYPE
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Score_Metrics_With_MC %>%
                                          group_by(CLAUDIN_SUBTYPE) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) +
    stat_kruskal_test(p.adjust.method = "BH", label = "Kruskal-Wallis, italic(p) = {p.adj.format}")

ggsave(
    "../../figures/Chapter_2/Global_CNA_Burden_Metrics_Across_PAM50.png",
    width = 15,
    height = 13,
    dpi = 120
)


# Distribution of CNA Amp/Del Scores across PAM50 Subtypes
CNA_Score_Metrics_Melt <-
    melt(CNA_Score_Metrics_With_MC[, c("Amp_Score_All",
                                       "Del_Score_All",
                                       "CLAUDIN_SUBTYPE",
                                       "INTCLUST")])
significance <-
    list(
        cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
        symbols = c("****", "***", "**", "*", "ns")
    )

ggplot(data = CNA_Score_Metrics_Melt,
       mapping = aes(
           x = CLAUDIN_SUBTYPE,
           y = abs(value),
           color = variable
       )) +
    geom_boxplot(aes(fill = variable), alpha = 0.3) +
    geom_point(position = position_jitterdodge(), alpha = 0.2)  +
    ggtitle("Boxplot of Amplification and Deletion Score Metrics by PAM50 Subtype") +
    ylab("CNA Score Metrics") + xlab("Subtype") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(
            colour = "black",
            size = 10,
            face = "bold"
        ),
        legend.text = element_text(size = 10),
        legend.position = "top"
    ) +
    scale_x_discrete(na.translate = TRUE) +
    theme(legend.position = "top") +
    scale_fill_discrete(name = "Score Metric:",
                        labels = c("Amplification Score", "Deletion Score")) +
    scale_color_discrete(name = "Score Metric:",
                         labels = c("Amplification Score", "Deletion Score")) +
    stat_kruskal_test(p.adjust.method = "BH", label = "p = {p.adj.format}")

ggsave(
    "../../figures/Chapter_2/Global_CNA_Score_AmpDel_Across_PAM50.png",
    width = 10,
    height = 6,
    dpi = 200
)

# Distribution of CNA Amp/Del Metrics across PAM50 Subtypes
CNA_Burden_Metrics_Melt <-
    melt(CNA_Burden_Metrics_With_MC[, c("Amp_Burden_All",
                                        "Del_Burden_All",
                                        "CLAUDIN_SUBTYPE",
                                        "INTCLUST")])

ggplot(data = CNA_Burden_Metrics_Melt,
       mapping = aes(
           x = CLAUDIN_SUBTYPE,
           y = abs(value),
           color = variable
       )) +
    geom_boxplot(aes(fill = variable), alpha = 0.3) +
    geom_point(position = position_jitterdodge(), alpha = 0.2)  +
    ggtitle("Boxplot of Amplification and Deletion Burden Metrics by PAM50 Subtype") +
    ylab("CNA Burden Metrics") + xlab("Subtype") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(
            colour = "black",
            size = 10,
            face = "bold"
        ),
        legend.text = element_text(size = 10),
        legend.position = "top"
    ) +
    scale_x_discrete(na.translate = TRUE) +
    theme(legend.position = "top") +
    scale_fill_discrete(name = "Burden Metric:",
                        labels = c("Amplification Burden", "Deletion Burden")) +
    scale_color_discrete(name = "Burden Metric:",
                         labels = c("Amplification Burden", "Deletion Burden")) +
    stat_kruskal_test(p.adjust.method = "BH", label = "p = {p.adj.format}")

ggsave(
    "../../figures/Chapter_2/Global_CNA_Burden_AmpDel_Across_PAM50.png",
    width = 10,
    height = 6,
    dpi = 200
)

# Dunns Test
## CNA Score
DT_1 <- dunn.test(
    CNA_Score_Metrics_With_MC$CNA_Score_All,
    CNA_Score_Metrics_With_MC$CLAUDIN_SUBTYPE,
    method = "bh",
    list = T
)

DT_2 <- dunn.test(
    CNA_Score_Metrics_With_MC$Amp_Score_All,
    CNA_Score_Metrics_With_MC$CLAUDIN_SUBTYPE,
    method = "bh",
    list = T
)

DT_3 <- dunn.test(
    CNA_Score_Metrics_With_MC$Del_Score_All,
    CNA_Score_Metrics_With_MC$CLAUDIN_SUBTYPE,
    method = "bh",
    list = T
)

pvals1 <- round(DT_1$P.adjusted, digits = 2)
pvals1[pvals1 < 0.0001] <- "<0.0001"

pvals2 <- round(DT_2$P.adjusted, digits = 2)
pvals2[pvals2 < 0.0001] <- "<0.0001"

pvals3 <- round(DT_3$P.adjusted, digits = 2)
pvals3[pvals3 < 0.0001] <- "<0.0001"

CNA_Score_DF <- data.frame(
    "Comparisons" = DT_1$comparisons,
    "AbsoluteCNAScore" =
        paste(round(DT_1$Z, digits = 2),
              " (", pvals1, ")",
              sep = ""),
    "CNAAmplificationScore" =
        paste(round(DT_2$Z, digits = 2),
              " (", pvals2, ")",
              sep = ""),
    "CNADeletionScore" =
        paste(round(DT_3$Z, digits = 2),
              " (", pvals3, ")",
              sep = ""),
    check.names = FALSE
)

z1 <- DT_1$Z
CNA_Score_DF <- CNA_Score_DF[order(-z1),]

CNA_Score_DF %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Score Metrics by PAM50 Subtype**")) %>%
    cols_align(align = "right",
               columns = c(2, 3, 4)) %>% cols_align(align = "left",
                                                    columns = c(1)) |>
    cols_label(
        AbsoluteCNAScore = "Absolute CNA Score <br> Z (adj p-value)",
        CNAAmplificationScore = "CNA Amp Score <br> Z (adj p-value)",
        CNADeletionScore = "CNA Del Score <br> Z (adj p-value)",
        .fn = md
    ) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    tab_style(style = "padding-top:10px;padding-bottom:10px;padding-left:12px;padding-right:12px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4))) %>%
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% gtsave("Global_CNA_Score_Metric_Comparisons.png", path = "../../tables/Chapter_2/")

## CNA Burden
DT_1 <- dunn.test(
    CNA_Burden_Metrics_With_MC$CNA_Burden_All,
    CNA_Burden_Metrics_With_MC$CLAUDIN_SUBTYPE,
    method = "bh",
    list = T
)

DT_2 <- dunn.test(
    CNA_Burden_Metrics_With_MC$Amp_Burden_All,
    CNA_Burden_Metrics_With_MC$CLAUDIN_SUBTYPE,
    method = "bh",
    list = T
)

DT_3 <- dunn.test(
    CNA_Burden_Metrics_With_MC$Del_Burden_All,
    CNA_Burden_Metrics_With_MC$CLAUDIN_SUBTYPE,
    method = "bh",
    list = T
)

pvals1 <- round(DT_1$P.adjusted, digits = 2)
pvals1[pvals1 < 0.0001] <- "<0.0001"

pvals2 <- round(DT_2$P.adjusted, digits = 2)
pvals2[pvals2 < 0.0001] <- "<0.0001"

pvals3 <- round(DT_3$P.adjusted, digits = 2)
pvals3[pvals3 < 0.0001] <- "<0.0001"

CNA_Burden_DF <- data.frame(
    "Comparisons" = DT_1$comparisons,
    "AbsoluteCNABurden" =
        paste(round(DT_1$Z, digits = 2),
              " (", pvals1, ")",
              sep = ""),
    "CNAAmplificationBurden" =
        paste(round(DT_2$Z, digits = 2),
              " (", pvals2, ")",
              sep = ""),
    "CNADeletionBurden" =
        paste(round(DT_3$Z, digits = 2),
              " (", pvals3, ")",
              sep = ""),
    check.names = FALSE
)

z1 <- DT_1$Z
CNA_Burden_DF <- CNA_Burden_DF[order(-z1),]

CNA_Burden_DF %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Burden Metrics by PAM50 Subtype**")) %>%
    cols_align(align = "right",
               columns = c(2, 3, 4)) %>% cols_align(align = "left",
                                                    columns = c(1)) |>
    cols_label(
        AbsoluteCNABurden = "Absolute CNA Burden <br> Z (adj p-value)",
        CNAAmplificationBurden = "CNA Amp Burden <br> Z (adj p-value)",
        CNADeletionBurden = "CNA Del Burden <br> Z (adj p-value)",
        .fn = md
    ) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    tab_style(style = "padding-top:10px;padding-bottom:10px;padding-left:12px;padding-right:12px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4))) %>%
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% gtsave("Global_CNA_Burden_Metric_Comparisons.png", path = "../../tables/Chapter_2/")
