## Library
library(tidyverse)

## Load up two datasets
Model4_Genes_lm <- read.delim("../../data/Chapter_6/Model4_Genes_func1_LM.txt")
Model4_Genes_MCMCglmm <- read.delim("../../data/Chapter_6/Model4_Genes_func1_MCMC.txt")

## How many didnt get fit at all?
nrow(Model4_Genes_MCMCglmm %>% filter(is.na(Category)))
nrow(Model4_Genes_lm %>% filter(is.na(Category)))

## Remove 10,934 genes that could be fit 
Model4_lm_Genes_Filtered <- Model4_Genes_lm[!is.na(Model4_Genes_lm$Category),]

Model4_lm_Genes_Filtered$Gene <- gsub("chr_\\d+_Gene_", "", Model4_lm_Genes_Filtered$Gene)
Model4_lm_Genes_Filtered$Gene <- gsub("chr_X_Gene_", "", Model4_lm_Genes_Filtered$Gene)

Model4_Genes_MCMCglmm$Gene <- gsub("chr_\\d+_Gene_", "", Model4_Genes_MCMCglmm$Gene)
Model4_Genes_MCMCglmm$Gene <- gsub("chr_X_Gene_", "", Model4_Genes_MCMCglmm$Gene)

Model4_mcmc_Genes_Filtered <- Model4_Genes_MCMCglmm[Model4_Genes_MCMCglmm$Gene %in% unique(Model4_lm_Genes_Filtered$Gene),]

## Lood up Location Data
CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg19.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

CNA_Loc <- CNA_Loc %>% select(1, 3,4)

CNA_Loc_Order <- CNA_Loc %>% arrange(Chromosome, Start)

# Tile Plots (lm model)
data <-  Model4_lm_Genes_Filtered

# Need to Order them!
data$Gene <- gsub("chr_\\d+_Gene_", "", data$Gene)
data$Gene <- gsub("chr_X_Gene_", "", data$Gene)

data <- data %>% filter(Samplesize != 0)

data <- data[order(match(data[,"Gene"],CNA_Loc_Order[,"Hugo_Symbol"])),]

## LM
## 0KB

grouped_data <- data %>% dplyr::mutate(Category = ifelse(Category == "NoChangepoint", "NoCP", Category)) %>%
    mutate(Category = factor(Category,  levels = rev(c("NoCP", "Amp/Neut", "Del/Neut", 
                                                       "Neut/Amp", "Neut/Del", 
                                                       "Amp/Del", "Del/Amp")))) %>%  mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
mutate(Chr = factor(Chr, levels = c(paste0("Chr_", c(1:22, "X"))))) %>%
    group_by(Chr, Gene, Category, Allele) %>%
    mutate(Sig = case_when(
        any(Direction == "TS" & LB > 10) & any(Direction == "TE" & LB > 10) ~ "Sig",
        any(Direction == "TS" & LB > 10) & any(Direction == "TE" & LB <= 10) ~ "Sig_TS",
        any(Direction == "TS" & LB <= 10) & any(Direction == "TE" & LB > 10) ~ "Sig_TE",
        TRUE ~ "NSig"
    )) %>%
    group_by(Chr, Gene, Category, Allele, Sig) %>%
    summarise(
        n_new = sum(Samplesize),
        n_Sig = sum(ifelse(Sig == "Sig", Samplesize, 0)),
        n_Sig_TS = ifelse(Sig == "Sig_TS", sum(Samplesize[Direction == "TS"]), 0),
        n_Sig_TE = ifelse(Sig == "Sig_TE", sum(Samplesize[Direction == "TE"]), 0)) %>% mutate("Final_n" = case_when(Sig == "Sig" ~ n_new,
                                                                                                                    Sig == "Sig_TS" ~ n_Sig_TS,
                                                                                                                    Sig == "Sig_TE" ~ n_Sig_TE,
                                                                                                                    Sig == "NSig" ~ NA))
grouped_data <- grouped_data[order(match(grouped_data$Gene, CNA_Loc_Order$Hugo_Symbol)),]
grouped_data$Gene <- factor(grouped_data$Gene, levels = unique(CNA_Loc_Order$Hugo_Symbol))

ggplot(grouped_data,
       aes(x = Gene, y = Category, fill = Sig)) + geom_tile()  + scale_fill_manual(
           values = c("white", "white", "white", "white"),
           na.value = "lightgrey",
           drop = F,
           guide = "none") + geom_point(aes(size = Final_n, colour = Sig),
                                        shape = 16) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())  + ggtitle("Changepoints of Significant Length across Genes (LB > 10kb)") + 
    theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(
        name = "Sig & Sample Size",
        labels = c("250", "500", "750"),
        values = c(16, 16, 16)
    ) +  scale_size_continuous(name = "Sample\nSize") +  facet_wrap(Chr~Allele, scales = "free_x", ncol = 6) + xlab("") + ylab("") +  theme(plot.title = element_text(size=18))

ggsave(paste0("PerGene_LM_Zero_Thesis.png"), 
       plot = last_plot(),
       path = "../../figures/Chapter_6/",
       width = 13,
       height = 15,
       units = c("in"),
       dpi = 300
)

## 10,000 kb
grouped_data <- data %>%  dplyr::mutate(Category = ifelse(Category == "NoChangepoint", "NoCP", Category)) %>%
    mutate(Category = factor(Category,  levels = rev(c("NoCP", "Amp/Neut", "Del/Neut", 
                                                       "Neut/Amp", "Neut/Del", 
                                                       "Amp/Del", "Del/Amp")))) %>%  
    mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
    mutate(Chr = factor(Chr, levels = c(paste0("Chr_", c(1:22, "X"))))) %>%
    group_by(Chr, Gene, Category, Allele) %>%
    mutate(Sig = case_when(
        any(Direction == "TS" & LB > 10000) & any(Direction == "TE" & LB > 10000) ~ "Sig",
        any(Direction == "TS" & LB > 10000) & any(Direction == "TE" & LB <= 10000) ~ "Sig_TS",
        any(Direction == "TS" & LB <= 10000) & any(Direction == "TE" & LB > 10000) ~ "Sig_TE",
        TRUE ~ "NSig"
    )) %>%
    group_by(Chr, Gene, Category, Allele, Sig) %>%
    summarise(
        n_new = sum(Samplesize),
        n_Sig = sum(ifelse(Sig == "Sig", Samplesize, 0)),
        n_Sig_TS = ifelse(Sig == "Sig_TS", sum(Samplesize[Direction == "TS"]), 0),
        n_Sig_TE = ifelse(Sig == "Sig_TE", sum(Samplesize[Direction == "TE"]), 0)) %>% mutate("Final_n" = case_when(Sig == "Sig" ~ n_new,
                                                                                                                    Sig == "Sig_TS" ~ n_Sig_TS,
                                                                                                                    Sig == "Sig_TE" ~ n_Sig_TE,
                                                                                                                    Sig == "NSig" ~ NA))
grouped_data <- grouped_data[order(match(grouped_data$Gene, CNA_Loc_Order$Hugo_Symbol)),]
grouped_data$Gene <- factor(grouped_data$Gene, levels = unique(CNA_Loc_Order$Hugo_Symbol))

ggplot(grouped_data,
       aes(x = Gene, y = Category, fill = Sig)) + geom_tile()  + scale_fill_manual(
    values = c("white", "white", "white", "white"),
    na.value = "lightgrey",
    drop = F,
    guide = "none") + geom_point(aes(size = Final_n, colour = Sig),
               shape = 16) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())  + ggtitle("Changepoints of Significant Length across Genes (LB > 10000kb)") + 
    theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(
        name = "Sig & Sample Size",
        labels = c("250", "500", "750"),
        values = c(16, 16, 16)
    ) +  scale_size_continuous(name = "Sample\nSize") +  facet_wrap(Chr~Allele, scales = "free_x", ncol = 6) + xlab("") + ylab("") +  theme(plot.title = element_text(size=18))

ggsave(paste0("PerGene_LM_10000_Thesis.png"), 
       plot = last_plot(),
       path = "../../figures/Chapter_6/",
       width = 13,
       height = 15,
       units = c("in"),
       dpi = 300
)



## MCMC 
data <- Model4_mcmc_Genes_Filtered
data$Gene <- gsub("chr_X_Gene_", "", data$Gene)
data$Gene <- gsub("chr_\\d+_Gene_", "", data$Gene)

data <- data %>% complete(Gene, Category, Allele) %>% group_by(Gene) %>% mutate(Chr = na.omit(unique(Chr)))
data <- data %>% filter(!is.na(Category)) %>% filter(!is.na(Allele))

data <- data[order(match(data$Gene,CNA_Loc_Order[,"Hugo_Symbol"])),]

## 10,000 kb
grouped_data <- data %>%  dplyr::mutate(Category = ifelse(Category == "NoChangepoint", "NoCP", Category)) %>%
    mutate(Category = factor(Category,  levels = rev(c("NoCP", "Amp/Neut", "Del/Neut", 
                                                       "Neut/Amp", "Neut/Del", 
                                                       "Amp/Del", "Del/Amp")))) %>%  
    mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
    mutate(Chr = factor(Chr, levels = c(paste0("Chr_", c(1:22, "X"))))) %>%
    group_by(Chr, Gene, Category, Allele) %>%
    mutate(Sig = case_when(
        is.na(Direction) & is.na(LB) ~ NA,
        any(Direction == "TS" & LB > 10000) & any(Direction == "TE" & LB > 10000) ~ "Sig",
        any(Direction == "TS" & LB > 10000) & any(Direction == "TE" & LB <= 10000) ~ "Sig_TS",
        any(Direction == "TS" & LB <= 10000) & any(Direction == "TE" & LB > 10000) ~ "Sig_TE",
        TRUE ~ "NSig"
    )) %>%
    group_by(Chr, Gene, Category, Allele, Sig) %>%
    summarise(
        n_new = sum(n),
        n_Sig = sum(ifelse(Sig == "Sig", n, 0)),
        n_Sig_TS = ifelse(Sig == "Sig_TS", sum(n[Direction == "TS"]), 0),
        n_Sig_TE = ifelse(Sig == "Sig_TE", sum(n[Direction == "TE"]), 0)) %>% mutate("Final_n" = case_when(Sig == "Sig" ~ n_new,
                                                                                                           Sig == "Sig_TS" ~ n_Sig_TS,
                                                                                                           Sig == "Sig_TE" ~ n_Sig_TE,
                                                                                                           Sig == "NSig" ~ NA))

grouped_data <- grouped_data[order(match(grouped_data$Gene, CNA_Loc_Order$Hugo_Symbol)),]
grouped_data$Gene <- factor(grouped_data$Gene, levels = unique(CNA_Loc_Order$Hugo_Symbol))

ggplot(grouped_data,
       aes(x = Gene, y = Category, fill = Sig)) + geom_tile()  + scale_fill_manual(
           values = c("white", "white", "white", "white"),
           na.value = "gray91",
           drop = F,
           guide = "none") + geom_point(aes(size = Final_n, colour = Sig),
                                        shape = 16) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())  + ggtitle("Changepoints of Significant Length across Genes (LB > 10000kb)") + 
    theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(
        name = "Sig & Sample Size",
        labels = c("250", "500", "750"),
        values = c(16, 16, 16)
    ) +  scale_size_continuous(name = "Sample\nSize") +  facet_wrap(Chr~Allele, scales = "free_x", ncol = 6) + xlab("") + ylab("") +  theme(plot.title = element_text(size=18))

ggsave(paste0("PerGene_MCMC_10000_Thesis.png"), 
       plot = last_plot(),
       path = "../../figures/Chapter_6/",
       width = 13,
       height = 15,
       units = c("in"),
       dpi = 300
)

## 0 kb
grouped_data <- data %>%  dplyr::mutate(Category = ifelse(Category == "NoChangepoint", "NoCP", Category)) %>%
    mutate(Category = factor(Category,  levels = rev(c("NoCP", "Amp/Neut", "Del/Neut", 
                                                       "Neut/Amp", "Neut/Del", 
                                                       "Amp/Del", "Del/Amp")))) %>%   
    mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
    mutate(Chr = factor(Chr, levels = c(paste0("Chr_", c(1:22, "X"))))) %>%
    group_by(Chr, Gene, Category, Allele) %>%
    mutate(Sig = case_when(
        is.na(Direction) & is.na(LB) ~ NA,
        any(Direction == "TS" & LB > 10) & any(Direction == "TE" & LB > 10) ~ "Sig",
        any(Direction == "TS" & LB > 10) & any(Direction == "TE" & LB <= 10) ~ "Sig_TS",
        any(Direction == "TS" & LB <= 10) & any(Direction == "TE" & LB > 10) ~ "Sig_TE",
        TRUE ~ "NSig"
    )) %>%
    group_by(Chr, Gene, Category, Allele, Sig) %>%
    summarise(
        n_new = sum(n),
        n_Sig = sum(ifelse(Sig == "Sig", n, 0)),
        n_Sig_TS = ifelse(Sig == "Sig_TS", sum(n[Direction == "TS"]), 0),
        n_Sig_TE = ifelse(Sig == "Sig_TE", sum(n[Direction == "TE"]), 0)) %>% mutate("Final_n" = case_when(Sig == "Sig" ~ n_new,
                                                                                                           Sig == "Sig_TS" ~ n_Sig_TS,
                                                                                                           Sig == "Sig_TE" ~ n_Sig_TE,
                                                                                                           Sig == "NSig" ~ NA))

grouped_data <- grouped_data[order(match(grouped_data$Gene, CNA_Loc_Order$Hugo_Symbol)),]
grouped_data$Gene <- factor(grouped_data$Gene, levels = unique(CNA_Loc_Order$Hugo_Symbol))

ggplot(grouped_data,
       aes(x = Gene, y = Category, fill = Sig)) + geom_tile()  + scale_fill_manual(
           values = c("white", "white", "white", "white"),
           na.value = "gray91",
           drop = F,
           guide = "none") + geom_point(aes(size = Final_n, colour = Sig),
                                        shape = 16) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())  + ggtitle("Changepoints of Significant Length across Genes (LB > 10kb)") + 
    theme(plot.title = element_text(hjust = 0.5)) +  scale_shape_manual(
        name = "Sig & Sample Size",
        labels = c("250", "500", "750"),
        values = c(16, 16, 16)
    ) +  scale_size_continuous(name = "Sample\nSize") +  facet_wrap(Chr~Allele, scales = "free_x", ncol = 6) + xlab("") + ylab("") +  theme(plot.title = element_text(size=18))

ggsave(paste0("PerGene_MCMC_Zero_Thesis.png"), 
       plot = last_plot(),
       path = "../../figures/Chapter_6/",
       width = 13,
       height = 15,
       units = c("in"),
       dpi = 300
)


## Tables
## LM

library(gt)
Model4_lm_Genes_Filtered$Gene <- gsub("chr_\\d+_Gene_", "", Model4_lm_Genes_Filtered$Gene)

Tab_1 <- Model4_lm_Genes_Filtered %>% select(Gene, Chr, Allele, Category, Samplesize, Direction, Pred, LB, UB) %>% 
    dplyr::mutate(Category = ifelse(Category == "NoChangepoint", "NoCP", Category)) %>% 
    dplyr::rename(n = Samplesize, Fit = Pred) %>%
    filter(LB > 10000, n > 30) %>% arrange(desc(Fit)) %>% 
    gt()  %>%
    tab_header(
        title = md("**(B) Genes Containing Changepoints with Large CNAs* (lm)*"), 
    )

Model4_mcmc_Genes_Filtered$Gene <- gsub("chr_\\d+_Gene_", "", Model4_mcmc_Genes_Filtered$Gene)

Tab_2 <- Model4_mcmc_Genes_Filtered %>% select(Gene, Chr, Allele, Category, n, Direction, Fit, LB, UB) %>% 
    dplyr::mutate(Category = ifelse(Category == "NoChangepoint", "NoCP", Category)) %>% 
    filter(LB > 10000, n > 30) %>% arrange(desc(Fit)) %>% 
    gt()  %>%
    tab_header(
        title = md("**(A) Genes Containing Changepoints with Large CNAs (MCMCglmm)**"), 
    )


Tab_1 |> gtsave(filename = "../../tables/Chapter_6/Gene_LM_1_Thesis.png", vwidth = 3000, vheight = 2500)
Tab_2 |> gtsave(filename = "../../tables/Chapter_6/Gene_MCMC_1_Thesis.png", vwidth = 3000, vheight = 2500)
