# Chapter 6: Barplots by Nodes

## Load up Libraries 
library(tidyverse)
library(rCGH)
library(gt)

## Load up Data 
ASCAT_Data_3Step_NoNeut<- read.delim("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))

CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg19.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

## Node 3p
CNA_Loc_3p <- CNA_Loc %>% filter(Chr_Loc == "3p") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "3") 

## Trees
data2 <-  readRDS(file = "../../data/Chapter_3/Chap3_List_of_PA_Data.rds")
data2 <- data2[[1]]

# data2 %>% group_by(Node) %>% count()
data2$PATIENT_ID <- gsub("\\-", ".", data2$PATIENT_ID)
data2$Node <-paste("Node", data2$Node, sep = " ")

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, data2, by.x = "Sample", by.y = "PATIENT_ID")
# length(unique(Total_ASCAT3_Chr3$Sample))

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "Node" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% mutate(Category = factor(Category), Allele = factor(Allele), Node = factor(Node))
    
    test <- test %>% 
        filter(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Changepoint) <= CNA_Loc_3p[i,]$End) 
    
    df_to_fill <- data.frame("Category" = rep(levels(test$Category), 6), "Allele" = rep(levels(test$Allele), 21), "n" = 0, 
                             "Gene" = CNA_Loc_3p[i,]$Hugo_Symbol, "Node" = rep(levels(test$Node), 14), Chromosome = 3)
    
    if(nrow(test) == 0){
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
        
    } else {
        test <- test %>% 
            group_by(Category, Allele, Chr, Node) %>% summarise(n = n()) %>%
            mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
        
        for(j in 1:nrow(test)){
            df_to_fill$n[df_to_fill$Category == test$Category[j] & df_to_fill$Allele == test$Allele[j] & df_to_fill$Gene == test$Gene[j] & df_to_fill$Chr == test$Chr[j] & 
                             df_to_fill$Node == test$Node[j]] <- test$n[j]
        }
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
    } 
}

Final <- merge(sample_df_2, CNA_Loc[,c(1,4)], by.x = "Gene" , by.y = "Hugo_Symbol")
sample_df_2_Order <- Final %>% arrange(Chromosome, Start)
sample_df_2_Order$Gene <- as.character(sample_df_2_Order$Gene)

p <- sample_df_2_Order  %>%
    group_by(Gene, Allele, Node) %>%
    summarise(Category = 'Total', 
              n = sum(n), 
              Node = Node, 
              Chromosome = Chromosome, 
              Start = Start) %>%
    bind_rows(sample_df_2_Order)

p <- p %>%   group_by(Gene, Category, Allele, Node) %>% distinct()

test <- p %>% select(Gene, Allele, Node, n) 

test1 <- pivot_wider(
    test,
    id_cols = c("Gene", "Node", "Category"),
    names_from = Allele,
    values_from = n,
    names_prefix = ""
) %>% mutate(Total = Minor + Major)

as.data.frame(test1) %>% filter(Node == "Node 2",  Category == "Total") %>% select(-Node, -Category) %>% arrange(desc(Total)) %>% dplyr::slice(1:10) %>% gt() |>
    tab_header(
        title =  md("**(A) Frequency of Changepoints within Genes (Node 2)**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Chr3_Node_2.png", path = "../../tables/Chapter_6/", vwidth = 500, vheight = 1000)

as.data.frame(test1) %>% filter(Node == "Node 4", Category == "Total") %>% select(-Node, -Category) %>% arrange(desc(Total)) %>% dplyr::slice(1:10) %>% gt() |>
    tab_header(
        title =  md("**(B) Frequency of Changepoints within Genes (Node 4)**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Chr3_Node_4.png", path = "../../tables/Chapter_6/", vwidth = 500, vheight = 1000)

as.data.frame(test1) %>% filter(Node == "Node 5",  Category == "Total") %>% select(-Node, -Category) %>% arrange(desc(Total)) %>% dplyr::slice(1:10) %>% gt() |>
    tab_header(
        title =  md("**(C) Frequency of Changepoints within Genes (Node 5)**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Chr3_Node_5.png", path = "../../tables/Chapter_6/", vwidth = 500, vheight = 1000)

Total_ASCAT3_Chr3 %>% distinct(Sample, Node, .keep_all = T) %>% group_by(Node) %>% summarise(n = n())

p <- p %>% mutate(Node = case_when(Node == "Node 2" ~ "Node 2 (1044)", 
                                   Node == "Node 4" ~ "Node 4 (790)", 
                                   Node == "Node 5" ~ "Node 5 (127)"))

as_tibble(p) %>% filter(Category != "NoChangepoint") %>% 
    mutate(Category = factor(Category,
                             levels = c("Total","Amp/Neut", "Del/Neut", "Neut/Amp", "Neut/Del", "Amp/Del", "Del/Amp"))) %>% 
    ggplot(aes(x=factor(Gene, levels = c(unique(sample_df_2_Order$Gene))), y=n, fill = Allele)) + 
    geom_bar(stat = "identity",  width=1) + facet_grid(Category~Node) + ggtitle("Frequency of Changepoints in Genes on Chromosome 3p by Node and Category") + 
    ylab("Frequency of Changepoints") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +1)) + theme(panel.spacing = unit(0.4, "cm")) + 
    theme(legend.position = "none")

library(patchwork)
design <- 
"
AA
BB
"
ploto <- as_tibble(p) %>% filter(Category != "NoChangepoint") %>% mutate(Category = factor(Category, levels = c("Total","Amp/Neut", "Del/Neut", 
                                                                                                                "Neut/Amp", "Neut/Del", 
                                                                                                                "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=factor(Gene, levels = c(unique(sample_df_2_Order$Gene))), y=n, fill = Allele)) + 
    geom_bar(stat = "identity",  width=1) + facet_grid(Category~Node) + ggtitle("Frequency of Changepoints in Genes on Chromosome 3p by Node and Category") + 
    ylab("Frequency of Changepoints") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                                              axis.title.x = element_text(vjust = +1)) + theme(panel.spacing = unit(0.4, "cm")) + 
    theme(legend.position = "none")

data2 <-  readRDS(file = "../../data/Chapter_3/Chap3_List_of_PA_Trees.rds")
data2 <- data2[[1]]

layout.matrix <- matrix(c(1, 2, 0, 0), nrow = 2, ncol = 1)

layout(mat = layout.matrix,
       heights = c(1, 2), # Heights of the two rows
       widths = c(2, 2)) # Widths of the two columns

layout.show(2)

par(mar = c(5, 4, 0, 0))
plot(data2)

par(mar = c(0, 4, 0, 0))
ploto
ggsave("../../figures/Chapter_6/Chromosome_3p_Barplot_Node.png", last_plot(), width = 9, height = 7)

## Tables 

## Node 18q
CNA_Loc_3p <- CNA_Loc %>% filter(Chr_Loc == "18q") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "18") 

data2 <-  readRDS(file = "../../data/Chapter_3/Chap3_List_of_PA_Data.rds")
data2 <- data2[[2]]
data2$PATIENT_ID <- gsub("\\-", ".", data2$PATIENT_ID)
data2$Node <-paste("Node", data2$Node, sep = " ")

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, data2, by.x = "Sample", by.y = "PATIENT_ID")

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "Node" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% mutate(Category = factor(Category), Allele = factor(Allele), Node = factor(Node))
    
    test <- test %>% 
        filter(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Changepoint) <= CNA_Loc_3p[i,]$End) 
    
    df_to_fill <- data.frame("Category" = rep(levels(test$Category), 6), "Allele" = rep(levels(test$Allele), 21), "n" = 0, 
                             "Gene" = CNA_Loc_3p[i,]$Hugo_Symbol, "Node" = rep(levels(test$Node), 14), Chromosome = 18)
    
    if(nrow(test) == 0){
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
        
    } else {
        test <- test %>% 
            group_by(Category, Allele, Chr, Node) %>% summarise(n = n()) %>%
            mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
        
        for(j in 1:nrow(test)){
            df_to_fill$n[df_to_fill$Category == test$Category[j] & df_to_fill$Allele == test$Allele[j] & df_to_fill$Gene == test$Gene[j] & df_to_fill$Chr == test$Chr[j] & 
                             df_to_fill$Node == test$Node[j]] <- test$n[j]
        }
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
    } 
}

Final <- merge(sample_df_2, CNA_Loc[,c(1,4)], by.x = "Gene" , by.y = "Hugo_Symbol")
sample_df_2_Order <- Final %>% arrange(Chromosome, Start)
sample_df_2_Order$Gene <- as.character(sample_df_2_Order$Gene)

p <- sample_df_2_Order  %>%
    group_by(Gene, Allele, Node) %>%
    summarise(Category = 'Total', 
              n = sum(n), 
              Node = Node, 
              Chromosome = Chromosome, 
              Start = Start) %>%
    bind_rows(sample_df_2_Order)

Total_ASCAT3_Chr3 %>% distinct(Sample, Node, .keep_all = T) %>% group_by(Node) %>% summarise(n = n())

p <- p %>%   group_by(Gene, Category, Allele, Node) %>% distinct()
p <- p %>% mutate(Node = case_when(Node == "Node 3" ~ "Node 3 (909)", 
                                   Node == "Node 4" ~ "Node 4 (125)", 
                                   Node == "Node 5" ~ "Node 5 (926)"))


as_tibble(p) %>% filter(Category != "NoChangepoint") %>% mutate(Category = factor(Category, levels = c("Total","Amp/Neut", "Del/Neut", 
                                                                                                                                                                                    "Neut/Amp", "Neut/Del", 
                                                                                                                                                                                    "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=factor(Gene, levels = c(unique(sample_df_2_Order$Gene))), y=n, fill = Allele)) + 
    geom_bar(stat = "identity",  width=1) + facet_grid(Category~Node) + ggtitle("Frequency of Changepoints in Genes on Chromosome 18q by Node and Category") + 
    ylab("Frequency of Changepoints") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +1)) + theme(panel.spacing = unit(0.4, "cm")) + 
    theme(legend.position = "none")


ggsave("../../figures/Chapter_6/Chromosome_18q_Barplot_Node.png", last_plot(), width = 9, height = 7)


## 11p 
CNA_Loc_3p <- CNA_Loc %>% filter(Chr_Loc == "11p") %>% select(1:7)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == "11") 

data2 <-  readRDS(file = "../../data/Chapter_3/Chap3_List_of_PA_Data.rds")
data2 <- data2[[3]]
data2$PATIENT_ID <- gsub("\\-", ".", data2$PATIENT_ID)
data2$Node <-paste("Node", data2$Node, sep = " ")

Total_ASCAT3_Chr3 <- merge(Total_ASCAT3_Chr3, data2, by.x = "Sample", by.y = "PATIENT_ID")

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "Node" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% mutate(Category = factor(Category), Allele = factor(Allele), Node = factor(Node))
    
    test <- test %>% 
        filter(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= CNA_Loc_3p[i,]$Start & as.numeric(Total_ASCAT3_Chr3$Changepoint) <= CNA_Loc_3p[i,]$End) 
    
    df_to_fill <- data.frame("Category" = rep(levels(test$Category), 10), "Allele" = rep(levels(test$Allele), 35), "n" = 0, 
                             "Gene" = CNA_Loc_3p[i,]$Hugo_Symbol, "Node" = rep(levels(test$Node), 14), Chromosome = 11)
    
    if(nrow(test) == 0){
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
        
    } else {
        test <- test %>% 
            group_by(Category, Allele, Chr, Node) %>% summarise(n = n()) %>%
            mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
        
        for(j in 1:nrow(test)){
            df_to_fill$n[df_to_fill$Category == test$Category[j] & df_to_fill$Allele == test$Allele[j] & df_to_fill$Gene == test$Gene[j] & df_to_fill$Chr == test$Chr[j] & 
                             df_to_fill$Node == test$Node[j]] <- test$n[j]
        }
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
    } 
}

Final <- merge(sample_df_2, CNA_Loc[,c(1,4)], by.x = "Gene" , by.y = "Hugo_Symbol")
sample_df_2_Order <- Final %>% arrange(Chromosome, Start)
sample_df_2_Order$Gene <- as.character(sample_df_2_Order$Gene)

p <- sample_df_2_Order  %>%
    group_by(Gene, Allele, Node) %>%
    summarise(Category = 'Total', 
              n = sum(n), 
              Node = Node, 
              Chromosome = Chromosome, 
              Start = Start) %>%
    bind_rows(sample_df_2_Order)

Total_ASCAT3_Chr3 %>% distinct(Sample, Node, .keep_all = T) %>% group_by(Node) %>% summarise(n = n())

p <- p %>%   group_by(Gene, Category, Allele, Node) %>% distinct()

test <- p %>% select(Gene, Allele, Node, n) 

test1 <- pivot_wider(
    test,
    id_cols = c("Gene", "Node", "Category"),
    names_from = Allele,
    values_from = n,
    names_prefix = ""
) %>% mutate(Total = Minor + Major)

as.data.frame(test1) %>% filter(Node == "Node 3",  Category == "Total") %>% select(-Node, -Category) %>% arrange(desc(Total)) %>% dplyr::slice(1:10) %>% gt() |>
    tab_header(
        title =  md("**(A) Frequency of Changepoints within Genes (Node 3)**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Chr11_Node_3.png", path = "../../tables/Chapter_6/", vwidth = 500, vheight = 1000)

as.data.frame(test1) %>% filter(Node == "Node 7",  Category == "Total") %>% select(-Node, -Category) %>% arrange(desc(Total)) %>% dplyr::slice(1:10) %>% gt() |>
    tab_header(
        title =  md("**(B) Frequency of Changepoints within Genes (Node 7)**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Chr11_Node_7.png", path = "../../tables/Chapter_6/", vwidth = 500, vheight = 1000)

p <- p %>% mutate(Node = case_when(Node == "Node 3" ~ "Node 3 (602)", 
                                   Node == "Node 4" ~ "Node 4 (95)", 
                                   Node == "Node 7" ~ "Node 7 (614)", 
                                   Node == "Node 8" ~ "Node 8 (222)", 
                                   Node == "Node 9" ~ "Node 9 (427)"))


as_tibble(p) %>% filter(Category != "NoChangepoint")  %>% mutate(Category = factor(Category, levels = c("Total","Amp/Neut", "Del/Neut", 
                                                                                                                                                                                    "Neut/Amp", "Neut/Del", 
                                                                                                                                                                                    "Amp/Del", "Del/Amp"))) %>% 
    
    ggplot(aes(x=factor(Gene, levels = c(unique(sample_df_2_Order$Gene))), y=n, fill = Allele)) + 
    geom_bar(stat = "identity",  width=1) + facet_grid(Category~Node) + ggtitle("Frequency of Changepoints in Genes on Chromosome 11p by Node and Category") + 
    ylab("Frequency of Changepoints") + xlab("Genes") + theme(plot.title = element_text(hjust = 0.5)) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), 
                                                                                              axis.title.x = element_text(vjust = +1)) + theme(panel.spacing = unit(0.4, "cm")) + 
    theme(legend.position = "none")


ggsave("../../figures/Chapter_6/Chromosome_11p_Barplot_Node.png", last_plot(), width = 9, height = 7)