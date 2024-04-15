# Chapter 2: Distribution of Chr Arm CNA Score across PAM50 Subtypes and IntClust
## Load up libraries
library(tidyverse)
library(stringr)
library(ggridges)
library(viridis)
library(ggh4x)
library(reshape2)
library(dunn.test)
library(gt)
library(gtExtras)

## Load up data
CNA_Arm_Score_All <-
    read.delim(
        "../../data/Processed_Data/Chr_Arm_CNA_Score.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )


CNA_Arm_Score_All  <-
    cbind.data.frame(data.frame("PATIENT_ID" = CNA_Arm_Score_All[, str_detect(colnames(CNA_Arm_Score_All), "PATIENT_ID")]),
                     CNA_Arm_Score_All[, str_detect(colnames(CNA_Arm_Score_All), "All")])
CNA_Arm_Burden_All <-
    read.delim(
        "../../data/Processed_Data/Chr_Arm_CNA_Burden.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )


CNA_Arm_Burden_All  <-
    cbind.data.frame(data.frame("PATIENT_ID" = CNA_Arm_Burden_All[, str_detect(colnames(CNA_Arm_Burden_All), "PATIENT_ID")]),
                     CNA_Arm_Burden_All[, str_detect(colnames(CNA_Arm_Burden_All), "All")])


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

CNA_Arm_Score_All_With_MC <- merge(CNA_Arm_Score_All, Clinical[,c("PATIENT_ID", "CLAUDIN_SUBTYPE", "INTCLUST")], by="PATIENT_ID")
CNA_Arm_Burden_All_With_MC <- merge(CNA_Arm_Burden_All, Clinical[,c("PATIENT_ID", "CLAUDIN_SUBTYPE", "INTCLUST")], by="PATIENT_ID")

# All Chromosomes (PAM50)
chrom <- c("1p", "1q", "2p", "2q", "3p", "3q", "4p", "4q",
           "5p", "5q", "6p", "6q", "7p", "7q", "8p", "8q",
           "9p", "9q", "10p", "10q", "11p", "11q", "12p",
           "12q", "13q", "14q", "15q", "16p", "16q", "17p",
           "17q", "18p", "18q", "19p", "19q", "20p", "20q",
           "21p", "21q", "22q", "Xp", "Xq")

for(i in chrom) {
    CNA_Arm_Score_All_With_MC_Melted <-
        melt(CNA_Arm_Score_All_With_MC) %>%
        filter(str_detect(variable, paste("_", i, sep = ""))) %>%
        mutate(value = ifelse(variable == str_detect(variable, "Del"), abs(value), value))
    
    CNA_Arm_Score_All_With_MC_Melted %>%
        ggplot(aes(y = CLAUDIN_SUBTYPE, x = value, fill = CLAUDIN_SUBTYPE)) +
        geom_density_ridges(alpha = 0.6) +
        scale_fill_viridis(discrete = TRUE) +
        facet_wrap( ~ variable, ncol = 2, scales = "free") +
        scale_color_viridis(discrete = TRUE) +
        ggtitle(paste(
            "Density Plot of CNA Score Metrics by PAM50 Subtype on Chromosome",
            i,
            sep = " "
        )) +
        ylab("Subtype") +
        xlab("CNA Score Metrics") +
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
        scale_fill_discrete(name = "Subtype:",
                            labels = c(paste(
                                c(levels(
                                    factor(CNA_Arm_Score_All_With_MC_Melted$CLAUDIN_SUBTYPE)
                                ), "NA"),
                                paste0(" (n = ", c(
                                    unlist(
                                        CNA_Arm_Score_All_With_MC_Melted %>%
                                            group_by(CLAUDIN_SUBTYPE) %>%
                                            summarise(number = n() /
                                                          6) %>%
                                            select(number)
                                    )
                                ), ")")
                            ))) +
        scale_colour_discrete(name = "Subtype:",
                              labels = c(paste(
                                  c(levels(
                                      factor(CNA_Arm_Score_All_With_MC_Melted$CLAUDIN_SUBTYPE)
                                  ), "NA"),
                                  paste0(" (n = ", c(
                                      unlist(
                                          CNA_Arm_Score_All_With_MC_Melted %>%
                                              group_by(CLAUDIN_SUBTYPE) %>%
                                              summarise(number = n() /
                                                            6) %>%
                                              select(number)
                                      )
                                  ), ")")
                              )))
    
    ggsave(
        paste(
            "../../figures/Chapter_2/ChrArm_DensityPlots/Score/CNA_Score_Metrics_Across_PAM50_",
            i  ,
            ".png"
        ),
        width = 16,
        height = 17
    )
}

for (i in chrom) {
    CNA_Arm_Burden_All_With_MC_Melted <-
        melt(CNA_Arm_Burden_All_With_MC) %>%
        filter(str_detect(variable, paste("_", i, sep = ""))) %>%
        mutate(value = ifelse(variable == str_detect(variable, "Del"), abs(value), value))
    
    CNA_Arm_Burden_All_With_MC_Melted %>%
        ggplot(aes(y = CLAUDIN_SUBTYPE, x = value, fill = CLAUDIN_SUBTYPE)) +
        geom_density_ridges(alpha = 0.6) +
        scale_fill_viridis(discrete = TRUE) +
        facet_wrap( ~ variable, ncol = 2, scales = "free") +
        scale_color_viridis(discrete = TRUE) +
        ggtitle(paste(
            "Density Plot of CNA Burden Metrics by PAM50 Subtype on Chromosome",
            i,
            sep = " "
        )) +
        ylab("Subtype") +
        xlab("CNA Burden Metrics") +
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
        scale_fill_discrete(name = "Subtype:",
                            labels = c(paste(
                                c(levels(
                                    factor(CNA_Arm_Burden_All_With_MC_Melted$CLAUDIN_SUBTYPE)
                                ), "NA"),
                                paste0(" (n = ", c(
                                    unlist(
                                        CNA_Arm_Burden_All_With_MC_Melted %>%
                                            group_by(CLAUDIN_SUBTYPE) %>%
                                            summarise(number = n() /
                                                          6) %>%
                                            select(number)
                                    )
                                ), ")")
                            ))) +
        scale_colour_discrete(name = "Subtype:",
                              labels = c(paste(
                                  c(levels(
                                      factor(CNA_Arm_Burden_All_With_MC_Melted$CLAUDIN_SUBTYPE)
                                  ), "NA"),
                                  paste0(" (n = ", c(
                                      unlist(
                                          CNA_Arm_Burden_All_With_MC_Melted %>%
                                              group_by(CLAUDIN_SUBTYPE) %>%
                                              summarise(number = n() /
                                                            6) %>%
                                              select(number)
                                      )
                                  ), ")")
                              )))
    
    ggsave(
        paste(
            "../../figures/Chapter_2/ChrArm_DensityPlots/Burden/CNA_Burden_Metrics_Across_PAM50_",
            i  ,
            ".png"
        ),
        width = 16,
        height = 17
    )
}


## PAM50 Individual Figures
### Basal

CNA_Arm_Burden_All_With_MC_Melted <-
    melt(CNA_Arm_Burden_All_With_MC) %>%
    filter(
        str_detect(variable, "CNA_Del_All_3p") |
            str_detect(variable, "CNA_Del_All_4p") |
            str_detect(variable, "CNA_Del_All_5q") |
            str_detect(variable, "CNA_Del_All_15q") |
            str_detect(variable, "CNA_Amp_All_3q") |
            str_detect(variable, "CNA_Amp_All_10p")
    ) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Del_All_3p" = "CNA Del Burden 3p",
            "CNA_Del_All_4p" = "CNA Del Burden 4p",
            "CNA_Del_All_5q" = "CNA Del Burden 5q",
            "CNA_Del_All_15q" = "CNA Del Burden 15q",
            "CNA_Amp_All_3q" = "CNA Amp Burden 3q",
            "CNA_Amp_All_10p" = "CNA Amp Burden 10p"
            
        )
    ) %>%
    mutate(value = ifelse(
        variable %in% c(
            "CNA Del Burden 3p",
            "CNA Del Burden 4p",
            "CNA Del Burden 5q",
            "CNA Del Burden 14q",
            "CNA Del Burden 15q",
            "CNA Amp Burden 3q",
            "CNA Amp Burden 10p"
        ),
        abs(value),
        value
    ))

CNA_Arm_Burden_All_With_MC_Melted %>%
    ggplot(aes(y = CLAUDIN_SUBTYPE, x = value, fill = CLAUDIN_SUBTYPE)) +
    geom_density_ridges(alpha = 0.6) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap( ~ variable, ncol = 2, scales = "free") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Density Plot of Selected CNA Burden Metrics by PAM50 Subtype (Basal)") +
    ylab("Subtype") +
    xlab("CNA Burden Metrics") +
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
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Arm_Burden_All_With_MC %>%
                                        group_by(CLAUDIN_SUBTYPE) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Arm_Burden_All_With_MC %>%
                                          group_by(CLAUDIN_SUBTYPE) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) + annotate("text", label = "Kruskal-Wallis, p<0.0001", size = 5, x = 90, y = 8)

ggsave(
    "../../figures/Chapter_2/ChrArm_CNA_Burden_Metrics_Across_PAM50_Basal_Burden.png",
    width = 16,
    height = 15
)

Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))
Pval <- c()
for(j in 1:length(Name)){
    
    t <-  kruskal.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$value, 
                       CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$CLAUDIN_SUBTYPE)
    
    Pval <- c(Pval, t$p.value)
}

Pval <- Pval*6

## Dunn's Test
Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))

DT_1 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_2 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_3 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[3], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[3], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_4 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[4], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[4], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_5 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[5], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[5], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_6 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[6], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[6], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)


pvals1 <- round(DT_1$P.adjusted, digits = 2)
pvals1[pvals1 < 0.0001] <- "<0.0001"
pvals2 <- round(DT_2$P.adjusted, digits = 2)
pvals2[pvals2 < 0.0001] <- "<0.0001"
pvals3 <- round(DT_3$P.adjusted, digits = 2)
pvals3[pvals3 < 0.0001] <- "<0.0001"
pvals4 <- round(DT_4$P.adjusted, digits = 2)
pvals4[pvals4 < 0.0001] <- "<0.0001"
pvals5 <- round(DT_5$P.adjusted, digits = 2)
pvals5[pvals5 < 0.0001] <- "<0.0001"
pvals6 <- round(DT_6$P.adjusted, digits = 2)
pvals6[pvals6 < 0.0001] <- "<0.0001"

CNA_Burden_DF_Basal <- data.frame(
    "Comparisons" = DT_1$comparisons,
    "CNA Del Burden 3p Z (adj p-value)" =
        paste(round(DT_1$Z, digits = 2),
              " (", pvals1, ")",
              sep = ""),
    "CNA Amp Burden 3q Z (adj p-value)" =
        paste(round(DT_2$Z, digits = 2),
              " (", pvals2, ")",
              sep = ""),
    "CNA Del Burden 4p Z (adj p-value)" =
        paste(round(DT_3$Z, digits = 2),
              " (", pvals3, ")",
              sep = ""),
    "CNA Del Burden 5q Z (adj p-value)" =
        paste(round(DT_4$Z, digits = 2),
              " (", pvals4, ")",
              sep = ""),
    "CNA Amp Burden 10p Z (adj p-value)" =
        paste(round(DT_5$Z, digits = 2),
              " (", pvals5, ")",
              sep = ""),
    "CNA Del Burden 15q Z (adj p-value)" =
        paste(round(DT_6$Z, digits = 2),
              " (", pvals6, ")",
              sep = ""),
    check.names = FALSE
)


z1 <- DT_1$Z
CNA_Burden_DF_Basal <- CNA_Burden_DF_Basal[order(-z1),]

CNA_Burden_DF_Basal %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Burden Metrics by PAM50 Subtype (Basal)**")) %>%
    cols_align(align = "right",
               columns = c(2, 3, 4,5,6,7)) %>% cols_align(align = "left",
                                                          columns = c(1)) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    gt_highlight_rows(
        rows = c(1,3,4,5,8),
        fill = "lightgrey",
        bold_target_only = TRUE
    ) %>%
    tab_style(style = "padding-top:10px;padding-bottom:10px;padding-left:12px;padding-right:12px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4,5,6,7))) %>%
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% gtsave("ChrArm_CNA_Burden_Metric_Comparisons_Basal.png", 
                                                   path = "../../tables/Chapter_2/", vwidth = 1370, vheight = 1500)


### HER2 Patients (Deletions/Amp)
CNA_Arm_Burden_All_With_MC_Melted <-
    melt(CNA_Arm_Burden_All_With_MC) %>%
    filter(
        str_detect(variable, "CNA_Amp_All_1q") |
            str_detect(variable, "CNA_Amp_All_17q") |
            str_detect(variable, "CNA_Del_All_8p") |
            str_detect(variable, "CNA_Amp_All_8q") |
            str_detect(variable, "CNA_Del_All_17p") |
            str_detect(variable, "CNA_Del_All_17q")
    ) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Amp_All_1q" = "CNA Amp Burden 1q",
            "CNA_Amp_All_17q" = "CNA Amp Burden 17q",
            "CNA_Amp_All_8q" = "CNA Amp Burden 8q",
            "CNA_Del_All_8p" = "CNA Del Burden 8p",
            "CNA_Del_All_17p" = "CNA Del Burden 17p",
            "CNA_Del_All_17q" = "CNA Del Burden 17q"
        )
    ) %>%
    mutate(value = ifelse(
        variable %in% c("CNA Del Burden 8p", "CNA Del Burden 17p",
                        "CNA Del Burden 17q"),
        abs(value),
        value
    ))

CNA_Arm_Burden_All_With_MC_Melted %>%
    ggplot(aes(y = CLAUDIN_SUBTYPE, x = value, fill = CLAUDIN_SUBTYPE)) +
    geom_density_ridges(alpha = 0.6) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap( ~ variable, ncol = 2, scales = "free") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Density Plot of CNA Burden Metrics by PAM50 Subtype (HER2)") +
    ylab("Subtype") +
    xlab("CNA Burden Metrics") +
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
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Arm_Burden_All_With_MC %>%
                                        group_by(CLAUDIN_SUBTYPE) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Arm_Burden_All_With_MC %>%
                                          group_by(CLAUDIN_SUBTYPE) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) + annotate("text", label = "Kruskal-Wallis, p<0.0001", size = 5, x = 90, y = 8)

ggsave(
    "../../figures/Chapter_2/ChrArm_CNA_Burden_Metrics_Across_PAM50_HER2_Burden.png",
    width = 16,
    height = 15
)


Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))
Pval <- c()
for(j in 1:length(Name)){
    
    t <-  kruskal.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$value, 
                       CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$CLAUDIN_SUBTYPE)
    
    Pval <- c(Pval, t$p.value)
}

Pval <- Pval*6

## Dunn's Test
Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))

DT_1 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_2 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_3 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[3], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[3], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_4 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[4], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[4], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_5 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[5], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[5], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_6 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[6], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[6], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)


pvals1 <- round(DT_1$P.adjusted, digits = 2)
pvals1[pvals1 < 0.0001] <- "<0.0001"
pvals2 <- round(DT_2$P.adjusted, digits = 2)
pvals2[pvals2 < 0.0001] <- "<0.0001"
pvals3 <- round(DT_3$P.adjusted, digits = 2)
pvals3[pvals3 < 0.0001] <- "<0.0001"
pvals4 <- round(DT_4$P.adjusted, digits = 2)
pvals4[pvals4 < 0.0001] <- "<0.0001"
pvals5 <- round(DT_5$P.adjusted, digits = 2)
pvals5[pvals5 < 0.0001] <- "<0.0001"
pvals6 <- round(DT_6$P.adjusted, digits = 2)
pvals6[pvals6 < 0.0001] <- "<0.0001"

Name
CNA_Burden_DF_HER2 <- data.frame(
    "Comparisons" = DT_1$comparisons,
    "CNA Amp Burden 1q Z (adj p-value)" =
        paste(round(DT_1$Z, digits = 2),
              " (", pvals1, ")",
              sep = ""),
    "CNA Del Burden 8p Z (adj p-value)" =
        paste(round(DT_2$Z, digits = 2),
              " (", pvals2, ")",
              sep = ""),
    " CNA Amp Burden 8q Z (adj p-value)" =
        paste(round(DT_3$Z, digits = 2),
              " (", pvals3, ")",
              sep = ""),
    "CNA Del Burden 17p Z (adj p-value)" =
        paste(round(DT_4$Z, digits = 2),
              " (", pvals4, ")",
              sep = ""),
    "CNA Amp Burden 17q Z (adj p-value)" =
        paste(round(DT_5$Z, digits = 2),
              " (", pvals5, ")",
              sep = ""),
    "CNA Del Burden 17q Z (adj p-value)" =
        paste(round(DT_6$Z, digits = 2),
              " (", pvals6, ")",
              sep = ""),
    check.names = FALSE
)

z1 <- DT_1$Z
CNA_Burden_DF_HER2 <- CNA_Burden_DF_HER2[order(-z1),]

CNA_Burden_DF_HER2 %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Burden Metrics by PAM50 Subtype (HER2)**")) %>%
    cols_align(align = "right",
               columns = c(2, 3, 4,5,6,7)) %>% cols_align(align = "left",
                                                          columns = c(1)) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    gt_highlight_rows(
        rows = c(6,7,10,12,13),
        fill = "lightgrey",
        bold_target_only = TRUE
    ) %>%
    tab_style(style = "padding-top:10px;padding-bottom:10px;padding-left:12px;padding-right:12px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4,5,6,7))) %>%
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% gtsave("ChrArm_CNA_Burden_Metric_Comparisons_HER2.png", 
                                                   path = "../../tables/Chapter_2/", vwidth = 1370, vheight = 1500)



### Luminal Patients (Deletions/Amp)
CNA_Arm_Burden_All_With_MC_Melted<-
    melt(CNA_Arm_Burden_All_With_MC) %>%
    filter(
        str_detect(variable, "CNA_Amp_All_16p") |
            str_detect(variable, "CNA_Del_All_16q") |
            str_detect(variable, "CNA_Del_All_11q") |
            str_detect(variable,  "CNA_Del_All_13q")
    ) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Amp_All_16p" = "CNA Amp Burden 16p",
            "CNA_Del_All_16q" = "CNA Del Burden 16q",
            "CNA_Del_All_11q" = "CNA Del Burden 11q",
            "CNA_Del_All_13q" = "CNA Del Burden 13q",
        )
    ) %>%
    mutate(value = ifelse(
        variable %in% c(
            "CNA Del Burden 8p",
            "CNA Del Burden 16q",
            "CNA Del Burden 11q",
            "CNA Del Burden 13q"
        ),
        abs(value),
        value
    ))

CNA_Arm_Burden_All_With_MC_Melted %>%
    ggplot(aes(y = CLAUDIN_SUBTYPE, x = value, fill = CLAUDIN_SUBTYPE)) +
    geom_density_ridges(alpha = 0.6) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap( ~ variable, ncol = 2, scales = "free") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Density Plot of CNA Burden Metrics by PAM50 Subtype (Luminal)") +
    ylab("Subtype") +
    xlab("CNA Burden Metrics") +
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
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Arm_Burden_All_With_MC %>%
                                        group_by(CLAUDIN_SUBTYPE) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Arm_Burden_All_With_MC %>%
                                          group_by(CLAUDIN_SUBTYPE) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) + annotate("text", label = "Kruskal-Wallis, p<0.0001", size = 5, x = 90, y = 7.9)

ggsave(
    "../../figures/Chapter_2/ChrArm_CNA_Burden_Metrics_Across_PAM50_Luminal_Burden.png",
    width = 15,
    height = 9
)

Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))
Pval <- c()
for(j in 1:length(Name)){
    
    t <-  kruskal.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$value, 
                       CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$CLAUDIN_SUBTYPE)
    
    Pval <- c(Pval, t$p.value)
}

Pval <- Pval*6

## Dunn's Test
Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))

DT_1 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_2 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_3 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[3], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[3], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

DT_4 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[4], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[4], ]$CLAUDIN_SUBTYPE,
                  method = "bh",
                  list = T)

pvals1 <- round(DT_1$P.adjusted, digits = 2)
pvals1[pvals1 < 0.0001] <- "<0.0001"
pvals2 <- round(DT_2$P.adjusted, digits = 2)
pvals2[pvals2 < 0.0001] <- "<0.0001"
pvals3 <- round(DT_3$P.adjusted, digits = 2)
pvals3[pvals3 < 0.0001] <- "<0.0001"
pvals4 <- round(DT_4$P.adjusted, digits = 2)
pvals4[pvals4 < 0.0001] <- "<0.0001"

Name
CNA_Burden_DF_Lum <- data.frame(
    "Comparisons" = DT_1$comparisons,
    "CNA Del Burden 11q Z (adj p-value)" =
        paste(round(DT_1$Z, digits = 2),
              " (", pvals1, ")",
              sep = ""),
    "CNA Del Burden 13q Z (adj p-value)" =
        paste(round(DT_2$Z, digits = 2),
              " (", pvals2, ")",
              sep = ""),
    "CNA Amp Burden 16p Z (adj p-value)" =
        paste(round(DT_3$Z, digits = 2),
              " (", pvals3, ")",
              sep = ""),
    "CNA Del Burden 16q Z (adj p-value)" =
        paste(round(DT_4$Z, digits = 2),
              " (", pvals4, ")",
              sep = ""),
    check.names = FALSE
)

z1 <- DT_1$Z
CNA_Burden_DF_Lum <- CNA_Burden_DF_Lum[order(-z1),]

CNA_Burden_DF_Lum %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Burden Metrics by PAM50 Subtype (Luminal)**")) %>%
    cols_align(align = "right",
               columns = c(2, 3, 4,5)) %>% cols_align(align = "left",
                                                      columns = c(1)) %>%
    tab_style(
        style = cell_borders(
            sides = c("all"),
            color = "lightgrey",
            weight = px(1),
            style = "solid"
        ),
        locations = cells_body()
    ) %>%
    gt_highlight_rows(
        rows = c(1,3,4,8,10,11,13,14,15),
        fill = "lightgrey",
        bold_target_only = TRUE
    ) %>%
    tab_style(style = "padding-top:10px;padding-bottom:10px;padding-left:12px;padding-right:12px",
              locations = cells_body()) %>%
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4,5))) %>%
    opt_table_outline() %>%
    tab_options(table.width = pct(100)) %>% gtsave("ChrArm_CNA_Burden_Metric_Comparisons_Luminal.png", 
                                                   path = "../../tables/Chapter_2/", vwidth = 1000, vheight = 1200)

## IntClust Density Plots
# All Chromosomes (IC)
for (i in chrom) {
    CNA_Arm_Score_All_With_MC_Melted <-
        melt(CNA_Arm_Score_All_With_MC) %>%
        filter(str_detect(variable, paste("_", i, sep = ""))) %>%
        mutate(value = ifelse(variable == str_detect(variable, "Del"), abs(value), value))
    
    CNA_Arm_Score_All_With_MC_Melted %>%
        ggplot(aes(y = INTCLUST, x = value, fill = INTCLUST)) +
        geom_density_ridges(alpha = 0.6) +
        scale_fill_viridis(discrete = TRUE) +
        facet_wrap(~ variable, ncol = 2, scales = "free") +
        scale_color_viridis(discrete = TRUE) +
        ggtitle(paste(
            "Density Plot of CNA Score Metrics by IntClust on Chromosome",
            i,
            sep = " "
        )) +
        ylab("Subtype") +
        xlab("CNA Score Metrics") +
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
        scale_fill_discrete(name = "Subtype:",
                            labels = c(paste(
                                c(levels(
                                    factor(CNA_Arm_Score_All_With_MC_Melted$INTCLUST)
                                ), "NA"),
                                paste0(" (n = ", c(
                                    unlist(
                                        CNA_Arm_Score_All_With_MC_Melted %>%
                                            group_by(INTCLUST) %>%
                                            summarise(number = n() /
                                                          6) %>%
                                            select(number)
                                    )
                                ), ")")
                            ))) +
        scale_colour_discrete(name = "Subtype:",
                              labels = c(paste(
                                  c(levels(
                                      factor(CNA_Arm_Score_All_With_MC_Melted$INTCLUST)
                                  ), "NA"),
                                  paste0(" (n = ", c(
                                      unlist(
                                          CNA_Arm_Score_All_With_MC_Melted %>%
                                              group_by(INTCLUST) %>%
                                              summarise(number = n() /
                                                            6) %>%
                                              select(number)
                                      )
                                  ), ")")
                              )))
    
    ggsave(
        paste(
            "../../figures/Chapter_2/ChrArm_DensityPlots/Score/CNA_Score_Metrics_Across_IC_",
            i  ,
            ".png"
        ),
        width = 16,
        height = 17
    )
}

for (i in chrom) {
    CNA_Arm_Burden_All_With_MC_Melted <-
        melt(CNA_Arm_Burden_All_With_MC) %>%
        filter(str_detect(variable, paste("_", i, sep = ""))) %>%
        mutate(value = ifelse(variable == str_detect(variable, "Del"), abs(value), value))
    
    CNA_Arm_Burden_All_With_MC_Melted %>%
        ggplot(aes(y = INTCLUST, x = value, fill = INTCLUST)) +
        geom_density_ridges(alpha = 0.6) +
        scale_fill_viridis(discrete = TRUE) +
        facet_wrap(~ variable, ncol = 2, scales = "free") +
        scale_color_viridis(discrete = TRUE) +
        ggtitle(paste(
            "Density Plot of CNA Burden Metrics by IntClust on Chromosome",
            i,
            sep = " "
        )) +
        ylab("Subtype") +
        xlab("CNA Burden Metrics") +
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
        scale_fill_discrete(name = "Subtype:",
                            labels = c(paste(
                                c(levels(
                                    factor(CNA_Arm_Burden_All_With_MC_Melted$INTCLUST)
                                ), "NA"),
                                paste0(" (n = ", c(
                                    unlist(
                                        CNA_Arm_Burden_All_With_MC_Melted %>%
                                            group_by(INTCLUST) %>%
                                            summarise(number = n() /
                                                          6) %>%
                                            select(number)
                                    )
                                ), ")")
                            ))) +
        scale_colour_discrete(name = "Subtype:",
                              labels = c(paste(
                                  c(levels(
                                      factor(CNA_Arm_Burden_All_With_MC_Melted$INTCLUST)
                                  ), "NA"),
                                  paste0(" (n = ", c(
                                      unlist(
                                          CNA_Arm_Burden_All_With_MC_Melted %>%
                                              group_by(INTCLUST) %>%
                                              summarise(number = n() /
                                                            6) %>%
                                              select(number)
                                      )
                                  ), ")")
                              )))
    
    ggsave(
        paste(
            "../../figures/Chapter_2/ChrArm_DensityPlots/Burden/CNA_Burden_Metrics_Across_IC_",
            i  ,
            ".png"
        ),
        width = 16,
        height = 17
    )
}


## IC Patients (Amplification Score/Burden)
CNA_Arm_Burden_All_With_MC_Melted <-
    melt(CNA_Arm_Burden_All_With_MC) %>%
    filter(str_detect(variable, "CNA_Del_All_3p") |
               str_detect(variable, "CNA_Del_All_4p")) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Del_All_3p" = "CNA Del Burden 3p",
            "CNA_Del_All_4p" = "CNA Del Burden 4p"
        )
    )

CNA_Arm_Burden_All_With_MC_Melted %>%
    ggplot(aes(y = INTCLUST, x = value, fill = INTCLUST)) +
    geom_density_ridges(alpha = 0.6) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap( ~ variable, ncol = 2, scales = "free") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Density Plot of CNA Burden Metrics by IntClust") +
    ylab("Subtype") +
    xlab("CNA Burden Metrics") +
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
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(CNA_Arm_Burden_All_With_MC$INTCLUST)), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Arm_Burden_All_With_MC %>%
                                        group_by(INTCLUST) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(CNA_Arm_Burden_All_With_MC$INTCLUST)), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Arm_Burden_All_With_MC %>%
                                          group_by(INTCLUST) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) + annotate("text", label = "Kruskal-Wallis, p<0.0001", size = 5, x = 90, y = 12.8)

Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))
Pval <- c()
for(j in 1:length(Name)){
    
    t <-  kruskal.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$value, 
                       CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                             Name[j], ]$CLAUDIN_SUBTYPE)
    
    Pval <- c(Pval, t$p.value)
}

Pval <- Pval*2

ggsave(
    "../../figures/Chapter_2/ChrArm_CNA_Burden_Metrics_Across_IC.png",
    width = 16,
    height = 7
)


Name <- c(unique(CNA_Arm_Burden_All_With_MC_Melted$variable))

DT_1 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[1], ]$INTCLUST,
                  method = "bh",
                  list = T)

DT_2 <- dunn.test(CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$value, 
                  CNA_Arm_Burden_All_With_MC_Melted[CNA_Arm_Burden_All_With_MC_Melted$variable == 
                                                        Name[2], ]$INTCLUST,
                  method = "bh",
                  list = T)

pvals1 <- round(DT_1$P.adjusted, digits = 2)
pvals1[pvals1 < 0.0001] <- "<0.0001"
pvals2 <- round(DT_2$P.adjusted, digits = 2)
pvals2[pvals2 < 0.0001] <- "<0.0001"


CNA_Burden_DF_IC <- data.frame(
    "Comparisons" = DT_1$comparisons,
    "CNA Del Burden 3p Z (adj p-value)" =
        paste(round(DT_1$Z, digits = 2),
              " (", pvals1, ")",
              sep = ""),
    "CNA Del Burden 4p Z (adj p-value)" =
        paste(round(DT_2$Z, digits = 2),
              " (", pvals2, ")",
              sep = ""),
    check.names = FALSE
)

z1 <- DT_1$Z
CNA_Burden_DF_IC <- CNA_Burden_DF_IC[order(-z1),]

CNA_Burden_DF_IC[1:27, ] %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Burden Metrics by IntClust**")) %>%
    cols_align(align = "right",
               columns = c(2, 3)) %>% cols_align(align = "left",
                                                 columns = c(1)) %>%
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
              locations = cells_column_labels(columns = c(1, 2, 3))) %>%
    opt_table_outline() %>% 
    tab_options(table.width = pct(60)) %>% gtsave("ChrArm_CNA_Burden_Metric_Comparisons_IC_1.png", 
                                                  path = "../../tables/Chapter_2/", vwidth = 890, vheight = 2000)


CNA_Burden_DF_IC[28:55, ] %>%
    gt() %>%
    tab_header(title =  md("**Comparison of CNA Burden Metrics by IntClust**")) %>%
    cols_align(align = "right",
               columns = c(2, 3)) %>% cols_align(align = "left",
                                                 columns = c(1)) %>%
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
              locations = cells_column_labels(columns = c(1, 2, 3))) %>%
    opt_table_outline() %>% 
    tab_options(table.width = pct(60)) %>% gtsave("ChrArm_CNA_Burden_Metric_Comparisons_IC_2.png", 
                                                  path = "../../tables/Chapter_2/", vwidth = 890, vheight = 2000)




## Presentation 
CNA_Arm_Burden_All_With_MC_Melted <-
    melt(CNA_Arm_Burden_All_With_MC) %>%
    filter(str_detect(variable, "CNA_Amp_All_1q") |
               str_detect(variable, "CNA_Del_All_4p") |
               str_detect(variable, "CNA_Del_All_5q") |
               str_detect(variable, "CNA_Amp_All_16p") |
               str_detect(variable, "CNA_Del_All_16q") |
               str_detect(variable, "CNA_Amp_All_17q")
    ) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Amp_All_4p" = "CNA Amp Burden 1q",
            "CNA_Del_All_4p" = "CNA Del Burden 4p",
            "CNA_Del_All_5q" = "CNA Del Burden 5q",
            "CNA_Amp_All_16p" = "CNA Amp Burden 16p",
            "CNA_Del_All_16q" = "CNA Del Burden 16q",
            "CNA_Amp_All_17p" = "CNA Amp Burden 17q",
        )
    ) %>%
    mutate(value = ifelse(
        variable %in% c(
            "CNA Amp Burden 1q",
            "CNA Del Burden 4p",
            "CNA Del Burden 5q",
            "CNA Amp Burden 16p",
            "CNA Del Burden 16q",
            "CNA Amp Burden 17q"
        ),
        abs(value),
        value
    ))

CNA_Arm_Burden_All_With_MC_Melted %>%
    ggplot(aes(y = CLAUDIN_SUBTYPE, x = value, fill = CLAUDIN_SUBTYPE)) +
    geom_density_ridges(alpha = 0.6) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap( ~ variable, ncol = 2, scales = "free") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Density Plot of Selected CNA Burden Metrics by PAM50 Subtype") +
    ylab("Subtype") +
    xlab("CNA Burden Metrics") +
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
    scale_fill_discrete(name = "Subtype:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Arm_Burden_All_With_MC %>%
                                        group_by(CLAUDIN_SUBTYPE) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "Subtype:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Arm_Burden_All_With_MC$CLAUDIN_SUBTYPE
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Arm_Burden_All_With_MC %>%
                                          group_by(CLAUDIN_SUBTYPE) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) + annotate("text", label = "Kruskal-Wallis, p<0.0001", size = 5, x = 90, y = 8)

ggsave(
    "../../figures/Chapter_2/ChrArm_CNA_Burden_Metrics_Across_PAM50_Deletions_Burden.png",
    width = 16,
    height = 15
)

## IntClust
CNA_Arm_Burden_All_With_MC_Melted <-
    melt(CNA_Arm_Burden_All_With_MC) %>%
    filter(str_detect(variable, "CNA_Amp_All_1q") |
               str_detect(variable, "CNA_Del_All_5q") |
               str_detect(variable, "CNA_Amp_All_16p") |
               str_detect(variable, "CNA_Amp_All_17q") 
    ) %>%
    mutate(
        variable = recode(
            variable,
            "CNA_Amp_All_1q" = "CNA Amp Burden 1q",
            "CNA_Del_All_8q" = "CNA Del Burden 5q",
            "CNA_Amp_All_16p" = "CNA Amp Burden 16p",
            "CNA_Amp_All_17q" = "CNA Amp Burden 17q",
        )
    ) %>%
    mutate(value = ifelse(
        variable %in% c(
            
            "CNA Amp Burden 1q",
            "CNA Del Burden 5q",
            "CNA Amp Burden 16p",
            "CNA Amp Burden 17q"
            
        ),
        abs(value),
        value
    ))

CNA_Arm_Burden_All_With_MC_Melted %>%
    ggplot(aes(y = INTCLUST, x = value, fill = INTCLUST)) +
    geom_density_ridges(alpha = 0.6) +
    scale_fill_viridis(discrete = TRUE) +
    facet_wrap( ~ variable, ncol = 2, scales = "free") +
    scale_color_viridis(discrete = TRUE) +
    ggtitle("Density Plot of Selected CNA Burden Metrics by IntClust") +
    ylab("IntClust") +
    xlab("CNA Burden Metrics") +
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
    scale_fill_discrete(name = "IntClust:",
                        labels = c(paste(c(
                            levels(factor(
                                CNA_Arm_Burden_All_With_MC$INTCLUST
                            )), "NA"
                        ),
                        paste0(
                            " (n = ", c(
                                unlist(
                                    CNA_Arm_Burden_All_With_MC %>%
                                        group_by(INTCLUST) %>%
                                        summarise(number = n()) %>%
                                        select(number)
                                )
                            ), ")"
                        )))) +
    scale_colour_discrete(name = "IntClust:",
                          labels = c(paste(c(
                              levels(factor(
                                  CNA_Arm_Burden_All_With_MC$INTCLUST
                              )), "NA"
                          ),
                          paste0(
                              " (n = ", c(
                                  unlist(
                                      CNA_Arm_Burden_All_With_MC %>%
                                          group_by(INTCLUST) %>%
                                          summarise(number = n()) %>%
                                          select(number)
                                  )
                              ), ")"
                          )))) + annotate("text", label = "Kruskal-Wallis, p<0.0001", size = 5, x = 90, y = 12.5)

ggsave(
    "../../figures/Chapter_2/ChrArm_CNA_Burden_Metrics_Across_PAM50_Amplifications_Burden.png",
    width = 16,
    height = 12
)

