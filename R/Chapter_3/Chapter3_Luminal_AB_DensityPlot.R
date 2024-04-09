# Chapter 3: CNA Score Metric Density Plot for Luminal Patients

## Load up necessary libraries
library(RColorBrewer)
library(ggplot2)

## Load up data 
Luminal_Data <- read.delim("../../data/Processed_Data/LuminalAB_Data.txt", 
                           sep = "\t")

## Density Plot
colourCount = 4
Palette = c("#3cb082", "#9ccb85", "#efb47a", "#e98472")
lab <- as.character(1:4)

dt <-
    data.frame(x = c(1:length(Luminal_Data$CNA_Score)),
               y = Luminal_Data$CNA_Score)
dt <- na.omit(dt)
dens <- density(dt$y)
df <- data.frame(x = dens$x, y = dens$y)
probs1 = c(0:4 / 4)
probs <- probs1[-c(1, length(probs1))]
quantiles <- quantile(dt$y, prob = probs)
df$quant <- factor(findInterval(df$x, quantiles))

## Plot
ggplot(df, aes(x, y)) +
    geom_line() +
    geom_ribbon(aes(ymin = 0, ymax = y, fill = quant)) +
    scale_x_continuous(breaks = quantiles) +
    xlab("CNA Score") +
    ylab("Density") +
    ggtitle("Segmented Density Plot of CNA Scores") +
    theme(legend.position = c(0.9, 0.7)) +
    theme(legend.key.size = unit(0.9, "cm")) +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    theme(axis.title.x = element_text(hjust = 0.5, size = 18)) +
    theme(axis.title.y = element_text(hjust = 0.5, size = 18)) +
    theme(axis.text.x = element_text(size = 15)) +
    theme(axis.text.y = element_text(size = 15)) +
    theme(legend.title = element_text(
        colour = "black",
        size = 15,
        face = "bold"
    )) +
    theme(strip.text = element_text(size = 15)) +
    theme(legend.text = element_text(size = 15)) +
    scale_fill_manual(
        values = Palette,
        labels = paste("CNA Q", lab, sep = ""),
        name = "Legend"
    )


ggsave("../../figures/Chapter_3/Luminal_AB_Score_Density.png", width = 12, height = 6)