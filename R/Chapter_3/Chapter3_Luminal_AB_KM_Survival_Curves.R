# Chapter 3: Survival Analysis Luminal A/B

## Load up necessary libraries
library(tidyverse)
library(RColorBrewer)
library(survival)
library(survminer)

## Load up data 
Luminal_Data <- read.delim("../../data/Processed_Data/LuminalAB_Data.txt", 
                           sep = "\t")

## Overall Survival
datafit <-
    survfit(Surv(OS_MONTHS, OS) ~ Subset_Quartile, data = Luminal_Data)

plot_km_OS <-
    ggsurvplot(
        datafit,
        censor.shape = "",
        data = Luminal_Data,
        size = 1,
        conf.int = F,
        pval = T,
        risk.table = T,
        xlab = "Overall Survival time (months)",
        ylab = "Survival probability",
        legend.title = "CNA Quartile",
        legend.labs = paste("Q", 1:4, sep = ""),
        pval.size = 3,
        risk.table.height = 0.28,
        ggtheme = theme_gray() +
            theme(plot.title = element_text(size = 12, hjust = 0.5)) +
            theme(legend.title = element_text(
                colour = "black",
                size = 11,
                face = "bold"
            )),
        break.time.by = 50,
        risk.table.y.text.col = T,
        risk.table.y.text = FALSE,
        title = "Luminal Breast Cancer Patients in METABRIC Data",
        font.main = c(12, "plain", "black"),
        font.x = c(11, "plain", "black"),
        font.y = c(11, "plain", "black"),
        font.legend = c(8, "plain", "black"),
        font.tickslab = c(8, "plain", "black"),
        pval.coord = c(250, 0.85),
        fontsize = 2,
        palette = c("#3cb082", "#9ccb85", "#efb47a", "#e98472"),
        legend = c(0.1, 0.35)
    )

plot_km_OS$table <-
    plot_km_OS$table + theme(
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9)
    )

png(
    file = "../../figures/Chapter_3/Luminal_AB_Score_OS.png",
    units = "cm",
    width = 15.24,
    height = 9.8, 
    res = 720
)
plot_km_OS
dev.off()


## Disease-Specific Survival
datafit <-
    survfit(Surv(OS_MONTHS, DSS) ~ Subset_Quartile, data = Luminal_Data)

plot_km_DSS <-
    ggsurvplot(
        datafit,
        censor.shape = "",
        data = Luminal_Data,
        size = 1,
        conf.int = F,
        pval = T,
        risk.table = T,
        xlab = "Disease-Specific Survival time (months)",
        ylab = "Survival probability",
        legend.title = "CNA Quartile",
        legend.labs = paste("Q", 1:4, sep = ""),
        pval.size = 3,
        risk.table.height = 0.28,
        ggtheme = theme_gray() +
            theme(plot.title = element_text(size = 12, hjust = 0.5)) +
            theme(legend.title = element_text(
                colour = "black",
                size = 11,
                face = "bold"
            )),
        break.time.by = 50,
        risk.table.y.text.col = T,
        risk.table.y.text = FALSE,
        title = "Luminal Breast Cancer Patients in METABRIC Data",
        font.main = c(12, "plain", "black"),
        font.x = c(11, "plain", "black"),
        font.y = c(11, "plain", "black"),
        font.legend = c(8, "plain", "black"),
        font.tickslab = c(8, "plain", "black"),
        pval.coord = c(250, 0.85),
        fontsize = 2,
        palette = c("#3cb082", "#9ccb85", "#efb47a", "#e98472"),
        legend = c(0.1, 0.35)
    )

plot_km_DSS$table <-
    plot_km_DSS$table + theme(
        plot.title = element_text(size = 10),
        axis.title.x = element_text(size = 9),
        axis.title.y = element_text(size = 9)
    )

png(
    file = "../../figures/Chapter_3/Luminal_AB_Score_DSS.png",
    units = "cm",
    width = 15.24,
    height = 9.8, 
    res = 720
)
plot_km_DSS
dev.off()