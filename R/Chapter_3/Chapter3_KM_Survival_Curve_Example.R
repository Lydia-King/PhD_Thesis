# Chapter 3: Example of KM survival curve

## Load up necessary libraries
library(survival)
library(survminer)
library(ggplot2)

## Load up data (Clinical patient and sample file from cBioPortal)
MB_Data <- read.delim("../../data/Processed_Data/METABRIC_PAM50_Data.txt", 
                      sep = "\t")

## Plot KM Curves
fit <-
    survfit(Surv(OS_MONTHS, DSS) ~ Her2_Status, data = MB_Data)

fit

ggsurvplot(
    fit,
    censor.shape = "+",
    xlab = "Survival time (months)",
    ylab = "Survival probability",
    data = MB_Data,
    surv.median.line = "hv",
    size = 1,
    conf.int = F,
    pval = T,
    risk.table = F,
    legend = c(0.88, 0.80),
    legend.title = "HER2 Status:",
    legend.labs = c("Negative", "Positive"),
    pval.size = 3,
    risk.table.height = 0.25,
    ggtheme = theme_gray()  + theme(plot.title = element_text(hjust = 0.5)),
    break.time.by = 50,
    risk.table.y.text.col = T,
    risk.table.y.text = FALSE,
    title = (
        paste("DSS", "for breast cancer patients in METABRIC data", sep = " ")
    ),
    font.main = c(12, "plain", "black"),
    font.x = c(10, "plain", "black"),
    font.y = c(10, "plain", "black"),
    font.legend = c(8, "plain", "black"),
    font.tickslab = c(8, "plain", "black")
)

ggsave("../../figures/Chapter_3/Example_Survival_Curve.png", width = 9, height = 5)
