# Chapter 5: Dotplot for selected scenario

## Load Libraries 
library(tidyverse)
library(MCMCglmm)

## Source functions 
source("Chapter5_Simulation_Study_Functions.R")

## Load up datasets 
Multiple_Dataset_Scen <- readRDS(file = "../../data/Chapter_5/Simulation_Study_Data.rds")

## Select Dataset to focus on 
Dataset <- Multiple_Dataset_Scen[[1]][[2]][[2]]

## Run tests
Test_Results_Scenario_Uni_lm <- Test_lm_MCMC_Results_Uni_DP(Dataset, tests = "lm")
Test_Results_Scenario_Uni_MCMC <- Test_lm_MCMC_Results_Uni_DP(Dataset, tests = "MCMC")
Test_Results_Scenario_Multi_lm <- Test_lm_MCMC_Results_Multi_DP(Dataset, tests = "lm")
Test_Results_Scenario_Multi_MCMC <- Test_lm_MCMC_Results_Multi_DP(Dataset, tests = "MCMC")

## Unlist Results 
Test_Results_Scenario_Uni_lm <- bind_rows(Test_Results_Scenario_Uni_lm, .id = "Dataset")
Test_Results_Scenario_Uni_MCMC <- bind_rows(Test_Results_Scenario_Uni_MCMC, .id = "Dataset")
Test_Results_Scenario_Multi_lm <- bind_rows(Test_Results_Scenario_Multi_lm, .id = "Dataset")
Test_Results_Scenario_Multi_MCMC <- bind_rows(Test_Results_Scenario_Multi_MCMC, .id = "Dataset")

Results_List <- list(Test_Results_Scenario_Uni_lm, Test_Results_Scenario_Uni_MCMC, 
                     Test_Results_Scenario_Multi_lm, Test_Results_Scenario_Multi_MCMC)

Names_List <- list("Uni_lm", "Uni_MCMC", "Multi_lm", "Multi_MCMC")

## Plot Interval Plots
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 8
cols = gg_color_hue(n)

for(i in 1:length(Results_List)){
    m1 <- 
        Results_List[[1]] %>% 
        mutate(Both = Category) %>%     #paste(Category, " (N = ", Samplesize,")", sep="")) %>%
        mutate(Direction = factor(Direction, levels = c("TS", "TE"))) %>% 
        filter(Samplesize > 0)
    
    Colors <-setNames(c(cols), c(unique(m1$Both)))
    
    png(filename = paste("../../figures/Chapter_5/", Names_List[[i]], "_Prediction_Simulation.png", sep=""), 
        width = 900, height = 500)
    
    if(i  %in% c(1,3)){
        print(m1 %>% ggplot(aes(x = Both, y = Pred, color = Both)) +
                  geom_point(size = 5) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  ggh4x::facet_grid2(Allele~Direction, scales = "free", independent = "x") + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + 
                  ggtitle("(A) Plot of Parameter Estimates from 20 Simulated Datasets"))
        dev.off()
    } else {
        print(m1 %>% ggplot(aes(x = Both, y = Pred, color = Both)) +
                  geom_point(size = 5) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  ggh4x::facet_grid2(Allele~Direction, scales = "free", independent = "x")  + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + 
                  ggtitle("(B) Plot of Parameter Estimates from 20 Simulated Datasets"))
        dev.off() 
    }
}
    

