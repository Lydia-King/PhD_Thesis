# Chapter 3: Adjusted Survival Curves Luminal A/B

## Load up necessary libraries
library(tidyverse)
library(survival)
library(survminer)

## Load up data
Luminal_Data <-
    read.delim("../../data/Processed_Data/LuminalAB_Data.txt", sep = "\t")  %>%
    mutate_at(c("Subset_Quartile",
                "GRADE",
                "TUMOR_STAGE"),
              funs(as.factor(as.character(.))))

## Final Cox model
res.cox_DSS <-
    coxph(
        Surv(OS_MONTHS, DSS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            TUMOR_SIZE +
            AGE_AT_DIAGNOSIS +
            Her2_Status +
            GRADE +
            CLAUDIN_SUBTYPE *
            Subset_Quartile,
        data = Luminal_Data
    )

## Set up 'new' dataframe
getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
}

NewData_Lum <-
    cbind("Subset_Quartile" = rep(1:4, 2),
          "CLAUDIN_SUBTYPE" = c(rep("LumA", 4), rep("LumB", 4)))

NewData_Lum <- as.data.frame(NewData_Lum)

## Hold everything else constant
NewData_Lum <- NewData_Lum %>% mutate(
    "LYMPH_NODES_EXAMINED_POSITIVE" = median(Luminal_Data$LYMPH_NODES_EXAMINED_POSITIVE, na.rm = T),
    "AGE_AT_DIAGNOSIS" =  median(Luminal_Data$AGE_AT_DIAGNOSIS, na.rm = T),
    "TUMOR_SIZE" = median(Luminal_Data$TUMOR_SIZE, na.rm = T),
    "GRADE" = getmode(Luminal_Data$GRADE),
    "Her2_Status" = getmode(Luminal_Data$Her2_Status)
)

NewData_Lum1 <- NewData_Lum %>% mutate(
    "LYMPH_NODES_EXAMINED_POSITIVE" = mean(Luminal_Data$LYMPH_NODES_EXAMINED_POSITIVE, na.rm = T),
    "AGE_AT_DIAGNOSIS" =  mean(Luminal_Data$AGE_AT_DIAGNOSIS, na.rm = T),
    "TUMOR_SIZE" = mean(Luminal_Data$TUMOR_SIZE, na.rm = T),
    "GRADE" = getmode(Luminal_Data$GRADE),
    "Her2_Status" = getmode(Luminal_Data$Her2_Status)
)

## Luminal A
fit <-
    survfit(res.cox_DSS, newdata = NewData_Lum1[NewData_Lum1$CLAUDIN_SUBTYPE == "LumA",])

plot_A <-
    ggsurvplot(
        fit,
        censor.shape = "",
        data = NewData_Lum,
        size = 1,
        conf.int = F,
        xlab = "Disease-Specific Survival time (months)",
        ylab = "Survival probability",
        legend.title = "CNA Quartile",
        legend.labs = paste("Q", 1:4, sep = ""),
        ggtheme = theme_gray() +
            theme(plot.title = element_text(size = 15, hjust = 0.5)) +
            theme(legend.title = element_text(
                colour = "black",
                size = 14,
                face = "bold"
            )),
        break.time.by = 50,
        title = "Adjusted Survival Curves for DSS (LumA)",
        font.main = c(15, "plain", "black"),
        font.x = c(14, "plain", "black"),
        font.y = c(14, "plain", "black"),
        font.legend = c(11, "plain", "black"),
        font.tickslab = c(11, "plain", "black"),
        palette = c("#3cb082", "#9ccb85", "#efb47a", "#e98472"),
        legend = c(0.15, 0.35)
    )

png(
    filename = "../../figures/Chapter_3/LuminalA_Adsurvplot.png",
    width = 440,
    height = 330,
    units = "px",
    bg = "white"
)
plot_A
invisible(dev.off())

## Luminal B
fit <-
    survfit(res.cox_DSS, newdata = NewData_Lum1[NewData_Lum1$CLAUDIN_SUBTYPE == "LumB",])

plot_B <-
    ggsurvplot(
        fit,
        censor.shape = "",
        data = NewData_Lum,
        size = 1,
        conf.int = F,
        xlab = "Disease-Specific Survival time (months)",
        ylab = "Survival probability",
        legend.title = "CNA Quartile",
        legend.labs = paste("Q", 1:4, sep = ""),
        ggtheme = theme_gray() +
            theme(plot.title = element_text(size = 15, hjust = 0.5)) +
            theme(legend.title = element_text(
                colour = "black",
                size = 14,
                face = "bold"
            )),
        break.time.by = 50,
        title = "Adjusted Survival Curves for DSS (LumB)",
        font.main = c(15, "plain", "black"),
        font.x = c(14, "plain", "black"),
        font.y = c(14, "plain", "black"),
        font.legend = c(11, "plain", "black"),
        font.tickslab = c(11, "plain", "black"),
        palette = c("#3cb082", "#9ccb85", "#efb47a", "#e98472"),
        legend = c(0.15, 0.35)
    )

png(
    filename = "../../figures/Chapter_3/LuminalB_Adsurvplot.png",
    width = 440,
    height = 330,
    units = "px",
    bg = "white"
)
print(plot_B)
invisible(dev.off())
