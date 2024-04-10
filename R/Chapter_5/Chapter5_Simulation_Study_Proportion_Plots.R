# Chapter 5: Simulation Study Proportion Plot

## Load Libraries 
library(ggh4x)
library(ggpubr)
library(tidyverse)

## Load up datasets 
## Univariate lm
Scenario_1_lm_uni <- read.delim("../../data/Chapter_5/Univariate_7_lm_Model_2_scenario_1.txt")
Scenario_2_lm_uni <- read.delim("../../data/Chapter_5/Univariate_7_lm_Model_2_scenario_2.txt")
Scenario_3_lm_uni <- read.delim("../../data/Chapter_5/Univariate_7_lm_Model_2_scenario_3.txt")

## Univariate MCMC
Scenario_1_mcmc_uni <- read.delim("../../data/Chapter_5/Univariate_7_MCMC_Model_2_scenario_1.txt")
Scenario_2_mcmc_uni <- read.delim("../../data/Chapter_5/Univariate_7_MCMC_Model_2_scenario_2.txt")
Scenario_3_mcmc_uni <- read.delim("../../data/Chapter_5/Univariate_7_MCMC_Model_2_scenario_3.txt")

## Multivariate lm
Scenario_1_lm_multi <- read.delim("../../data/Chapter_5/Multivariate_7_lm_Model_2_scenario_1.txt")
Scenario_2_lm_multi <- read.delim("../../data/Chapter_5/Multivariate_7_lm_Model_2_scenario_2.txt")
Scenario_3_lm_multi <- read.delim("../../data/Chapter_5/Multivariate_7_lm_Model_2_scenario_3.txt")

## Multivariate MCMC
Scenario_1_mcmc_multi <- read.delim("../../data/Chapter_5/Multivariate_7_MCMC_Model_2_scenario_1.txt")
Scenario_2_mcmc_multi <- read.delim("../../data/Chapter_5/Multivariate_7_MCMC_Model_2_scenario_2.txt")
Scenario_3_mcmc_multi <- read.delim("../../data/Chapter_5/Multivariate_7_MCMC_Model_2_scenario_3.txt")


## Plot Function
Plot_Function <- function(Data1, Data2, Data3){
    scen1_plot <- ggplot(Data1 %>% filter(SampleSize == 20 | SampleSize == 50 | SampleSize == 100 | SampleSize == 200) %>% mutate(Direction = factor(Direction, levels = c("TS", "TE"))), aes(x=as.factor(Percentage), y=Prop, group=Category)) +
        geom_point(aes(color=Category, shape=Category), alpha = .5, size = 2) + 
        facet_nested(Allele ~ SampleSize + Direction) + ggtitle("Scenario 1") +xlab("Simulated percentage") + ylab("Proportion datasets significant") +
        geom_line(aes(color=Category), alpha = .5, linewidth=1)+
        theme(plot.title = element_text(hjust = 0.5, size = 20, face="bold")) +
        theme(axis.text.y = element_text(size = 11)) +  theme(axis.text.x = element_text(size = 9)) +
        theme(legend.text = element_text(size = 12))+ theme(strip.text.x = element_text(size = 16)) + theme(legend.title = element_text(size=13)) +
        theme(axis.title.x = element_text(size = 14)) + theme(plot.subtitle=element_text(hjust=0.5, size= 21, face="bold")) +
        theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=14), legend.spacing.x = unit(0.2, 'cm')) +
        theme(strip.text = element_text(size = 16)) + theme(axis.title.y = element_text(size=14)) + guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(color=guide_legend(nrow=2,byrow=TRUE)) + 
        theme(plot.margin = unit(c(1,0,2,0), 'lines'))
    
    scen2_plot <-  ggplot(Data2 %>% filter(SampleSize == 20 | SampleSize == 50 | SampleSize == 100 | SampleSize == 200) %>% mutate(Direction = factor(Direction, levels = c("TS", "TE"))), aes(x=as.factor(Percentage), y=Prop, group=Category)) +
        geom_point(aes(color=Category, shape=Category), alpha = .5, size = 2) + facet_nested(Allele ~ SampleSize + Direction) + ggtitle("Scenario 2") +xlab("Simulated percentage") + ylab("Proportion datasets significant") +
        geom_line(aes(color=Category), alpha = .5, linewidth=1)+
        theme(plot.title = element_text(hjust = 0.5, size = 20, face="bold")) +
        theme(axis.text.y = element_text(size = 11)) +  theme(axis.text.x = element_text(size = 9)) +
        theme(legend.text = element_text(size = 12))+ theme(strip.text.x = element_text(size = 16)) + theme(legend.title = element_text(size=13)) +
        theme(axis.title.x = element_text(size = 14)) + theme(plot.subtitle=element_text(hjust=0.5, size= 21, face="bold")) +
        theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=14), legend.spacing.x = unit(0.2, 'cm')) +
        theme(strip.text = element_text(size = 16)) + theme(axis.title.y = element_text(size=14)) + guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(color=guide_legend(nrow=2,byrow=TRUE)) +
        theme(plot.margin = unit(c(1,0,2,0), 'lines'))
    
    scen3_plot <-  ggplot(Data3 %>% filter(SampleSize == 20 | SampleSize == 50 | SampleSize == 100 | SampleSize == 200) %>% mutate(Direction = factor(Direction, levels = c("TS", "TE"))), aes(x=as.factor(Percentage), y=Prop, group=Category)) +
        geom_point(aes(color=Category, shape=Category), alpha = .5, size = 2) + facet_nested(Allele ~ SampleSize + Direction) + ggtitle("Scenario 3") +xlab("Simulated percentage") + ylab("Proportion datasets significant") +
        geom_line(aes(color=Category), alpha = .5, linewidth=1)+
        theme(plot.title = element_text(hjust = 0.5, size = 20, face="bold")) +
        theme(axis.text.y = element_text(size = 11)) +  theme(axis.text.x = element_text(size = 9)) +
        theme(legend.text = element_text(size = 12))+ theme(strip.text.x = element_text(size = 16)) + theme(legend.title = element_text(size=13)) +
        theme(axis.title.x = element_text(size = 14)) + theme(plot.subtitle=element_text(hjust=0.5, size= 21, face="bold")) +
        theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=14), legend.spacing.x = unit(0.2, 'cm')) +
        theme(strip.text = element_text(size = 16)) + theme(axis.title.y = element_text(size=14)) + guides(shape=guide_legend(nrow=2,byrow=TRUE)) + guides(color=guide_legend(nrow=2,byrow=TRUE)) +
        theme(plot.margin = unit(c(1,0,2,0), 'lines'))
    
    return(list(scen1_plot, scen2_plot, scen3_plot))
}

## Plot Outputs
### lm Univariate Plot
Uni_lm <- Plot_Function(Scenario_1_lm_uni, Scenario_2_lm_uni, Scenario_3_lm_uni)
png(filename = paste("../../figures/Chapter_5/lm_uni_sim_plot", ".png", sep=""), width = 900, height = 1150)

plot <- ggarrange(Uni_lm[[1]], Uni_lm[[2]], Uni_lm[[3]], ncol = 1)

annotate_figure(plot, top = text_grob("Univariate lm ADIM", color = "black", face = "bold", size = 26))
dev.off()

### MCMC Univariate Plot
Uni_mcmc <- Plot_Function(Scenario_1_mcmc_uni, Scenario_2_mcmc_uni, Scenario_3_mcmc_uni)
png(filename = paste("../../figures/Chapter_5/mcmc_uni_sim_plot", ".png", sep=""), width = 900, height = 1150)

plot <- ggarrange(Uni_mcmc[[1]], Uni_mcmc[[2]], Uni_mcmc[[3]], ncol = 1)

annotate_figure(plot, top = text_grob("Univariate MCMCglmm ADIM", color = "black", face = "bold", size = 26))
dev.off()

### lm Multivariate Plot
Multi_lm <- Plot_Function(Scenario_1_lm_multi, Scenario_2_lm_multi, Scenario_3_lm_multi)
png(filename = paste("../../figures/Chapter_5/lm_multi_sim_plot", ".png", sep=""), width = 900, height = 1150)

plot <- ggarrange(Multi_lm[[1]], Multi_lm[[2]], Multi_lm[[3]], ncol = 1)

annotate_figure(plot, top = text_grob("Multivariate lm ADIM", color = "black", face = "bold", size = 26))
dev.off()

### MCMC Multivariate Plot
Multi_mcmc <- Plot_Function(Scenario_1_mcmc_multi, Scenario_2_mcmc_multi, Scenario_3_mcmc_multi)
png(filename = paste("../../figures/Chapter_5/mcmc_multi_sim_plot", ".png", sep=""), width = 900, height = 1150)

plot <- ggarrange(Multi_mcmc[[1]], Multi_mcmc[[2]], Multi_mcmc[[3]], ncol = 1)

annotate_figure(plot, top = text_grob("Multivariate MCMCglmm ADIM", color = "black", face = "bold", size = 26))
dev.off()