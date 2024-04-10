# Chapter 4: Global CNA Volcano Plots (Based on CNA Burden Tree Nodes)

## Load Libraries 
library(fastDummies)
library(data.table)
library(ggplot2)
library(ggrepel)
library(cowplot)
library(patchwork)
library(operator.tools)
library(tidyverse)
library(limma)

## Load up Trees 
List_of_Trees <- 
    readRDS(file = "../../data/Chapter_3/Chap3_List_of_Tree_Global.rds")

List_of_Data <- 
    readRDS(file = "../../data/Chapter_3/Chap3_List_of_Data_Global.rds")

## Volcano Plot Function
Volcano_Plot_Function <- function(de, title, save_path_name, width, height, x_1, x_2){
    
    # Label Genes as Up or Down Reg
    de$diffexpressed <- "Not Sig"
    de$diffexpressed[de$logFC > 0.58 & de$adj.P.Val < 0.05] <- "Up-reg"
    de$diffexpressed[de$logFC < -0.58 & de$adj.P.Val < 0.05] <- "Down-reg"
    de$AlogFC <- abs(de$logFC)
    de <- de[order(-de$AlogFC, de$adj.P.Val), ]
    de$rank <- 1:nrow(de)
    de$delabel_1 <- NA
    de$gene_symbol <- rownames(de)
    de$delabel_1[de$diffexpressed != "Not Sig" & de$rank %in% c(1:15)] <- 
        de$gene_symbol[de$diffexpressed != "Not Sig" & de$rank %in% c(1:15)]
    de$diffexpressed <- as.factor(de$diffexpressed)
    levels(de$diffexpressed) <- c(levels(de$diffexpressed), "Down-reg")
    
    # Plot
    p <- ggplot(data = de,
                aes(
                    x = logFC,
                    y = -log10(adj.P.Val),
                    col = diffexpressed,
                    label = delabel_1
                )) +
        geom_point() +
        theme_minimal() +
        geom_text_repel(size = 3, nudge_x = 0.15) +
        scale_color_manual(
            "",
            values = c(
                "Down-reg" = "blue",
                "Not Sig" = "grey",
                "Up-reg" = "red"
            ),
            labels = c("Down", "Not Sig", "Up")
        ) +
        geom_vline(
            xintercept = c(-0.58, 0.58),
            col = "black",
            linetype = 3
        ) +
        geom_hline(
            yintercept = -log10(0.05),
            col = "black",
            linetype = 3
        ) + ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5)) + 
        theme(title = element_text(size = 20)) + 
        xlab("Log Fold Change") + 
        ylab("-Log10 (Adjusted P-value)") + 
        theme_grey(base_size = 12) + 
        theme(plot.title = element_text(hjust = 0.5)) +
        xlim(c(x_1, x_2)) + 
        theme(legend.position = "right")
    
    
    ggsave(save_path_name, p, width = width, height = height)
}

## Load up Data
GE_Loc <- read.delim("../../data/Processed_Data/GE_Common_Genes.txt", sep = "\t")

## Gene Expression Analysis - Data setup
### Set up Gene Expression Data 
GE_Loc$Amean <- rowMeans(GE_Loc[,3:1982])
o <- order(GE_Loc$Amean, decreasing=TRUE)
data.fit2 <- GE_Loc[o,]
d <- duplicated(data.fit2$Hugo_Symbol)

data.fit3 <- data.fit2[!d,]

Gene_Log_1 <- data.fit3 %>% select(-Amean)

### Get Common Patients (Patients who have CNA and Gene Expression Data)
colnames(GE_Loc) <- gsub("\\-", ".", colnames(GE_Loc))
PATIENT_ID <- unique(c(List_of_Data[[1]]$PATIENT_ID, List_of_Data[[2]]$PATIENT_ID, List_of_Data[[3]]$PATIENT_ID, List_of_Data[[4]]$PATIENT_ID))
Common_Patients <- intersect(PATIENT_ID, colnames(Gene_Log_1))
Gene_Log_1_Sub <- Gene_Log_1 %>% dplyr::select(c(colnames(Gene_Log_1[colnames(Gene_Log_1) %in% Common_Patients]), "Hugo_Symbol"))

### Set up Design List 
Fit_List <- List_of_Trees
Data_List <- List_of_Data

Nodes_List <- list(c("Node_2", "Node_4", "Node_5"), 
                   c("Node_2", "Node_5", "Node_6", "Node_8", "Node_9"), 
                   c("Node_3", "Node_4", "Node_5"), 
                   c("Node_3", "Node_4", "Node_6", "Node_7"))

Genes_List <- list()
Design_List <- list()

for(i in 1:length(Fit_List)) {
    Genes_of_Interest <- Gene_Log_1_Sub$Hugo_Symbol
    
    Gene_Log_of_Interest <-
        Gene_Log_1_Sub %>% filter(Hugo_Symbol%in% Genes_of_Interest)
    
    row.names(Gene_Log_of_Interest) <-
        Gene_Log_of_Interest$Hugo_Symbol
    
    Gene_Log_of_Interest <-
        Gene_Log_of_Interest %>% dplyr::select(-c(Hugo_Symbol))
    
    Data_Ctree <- Data_List[[i]]
    
    Design <-
        Data_Ctree %>% filter(PATIENT_ID %in% c(colnames(Gene_Log_of_Interest)))
    
    Design_Dummy <- dummy_cols(Design, select_columns = "Node")
    
    rownames(Design_Dummy) <- Design$PATIENT_ID
    
    colnames(Design_Dummy) <-
        c("PATIENT_ID", "Node", c(Nodes_List[[i]]))
    
    Design_Dummy <-
        Design_Dummy %>% dplyr::select(-c(Node, PATIENT_ID))
    
    # Make column order gene_log_of_interest are in same order as rows in design
    k <- rownames(Design_Dummy)
    
    Gene_Log_of_Interest <-
        setcolorder(Gene_Log_of_Interest, as.character(k))
    
    # Data for glm function
    Design <-
        Data_Ctree %>% 
        filter(PATIENT_ID %in% c(colnames(Gene_Log_of_Interest)))
    
    data_expr <- t(as.data.frame(Gene_Log_of_Interest))
    data <- as.data.frame(data_expr)
    data$PATIENT_ID <- rownames(data)
    data <- merge(data, Design, by = "PATIENT_ID")
    
    Genes_List[[i]] <- Gene_Log_of_Interest
    Design_List[[i]] <- data
}

## Gene Expression Analysis - Fit Volcano Plots
## Volcano Plot 1
data <- Design_List[[1]]
de2 <- model.matrix(~ 0 + data$Node)
colnames(de2) <- c("Node2", "Node4", "Node5")
aw <- arrayWeights(Genes_List[[1]], method = "reml", design = de2)

fit <- lmFit(Genes_List[[1]], design=de2, weights = aw)

contrast.matrix <- makeContrasts(
    BVsT2 = Node5 - Node4,
    levels=colnames(de2))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
de <- topTable(fit2, n=nrow(fit2), adjust.method = "BH", coef = 1) 

nrow(de %>% filter(abs(logFC) > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC < -0.58 & de$adj.P.Val < 0.05))

Volcano_Plot_Function(de, 
                      title = "(A) Volcano Plot for Node 5 versus Node 4", 
                      save_path_name = "../../figures/Chapter_4/Volcano_1.png", 
                      width = 9, height = 5, x_1 = -1.8, x_2 = 1.8)

## Volcano Plot 2
data <- Design_List[[2]]
de2 <- model.matrix(~ 0 + data$Node)
colnames(de2) <- c("Node2", "Node5", "Node6", "Node8", "Node9")

aw <- arrayWeights(Genes_List[[2]], method = "reml", design = de2)

fit <- lmFit(Genes_List[[2]], design=de2, weights = aw)

contrast.matrix <- makeContrasts(
    BVsT = Node9 - Node8,
    levels=colnames(de2))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

de <- topTable(fit2, n=nrow(fit2), adjust.method = "BH") 

nrow(de %>% filter(abs(logFC) > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC < -0.58 & de$adj.P.Val < 0.05))

Volcano_Plot_Function(de, 
                      title = "(B) Volcano Plot for Node 9 versus Node 8", 
                      save_path_name = "../../figures/Chapter_4/Volcano_2.png", 
                      width = 9, height = 5, x_1 = -1.1, x_2 = 1.1)

## Volcano Plot 3
data <- Design_List[[3]]
de2 <- model.matrix(~ 0 + data$Node)
colnames(de2) <- c("Node3", "Node4", "Node5")

Genes_List[[3]] <- as.data.frame(Genes_List[[3]]) %>% 
    select(which(colnames(as.data.frame(Genes_List[[3]])) %in% data$PATIENT_ID))

aw <- arrayWeights(Genes_List[[3]], method = "reml", design = de2)

fit <- lmFit(Genes_List[[3]], design=de2, weights = aw)

contrast.matrix <- makeContrasts(
    BVsT = Node4 - Node3 ,
    levels=colnames(de2))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

de <- topTable(fit2, n=nrow(fit2), adjust.method = "BH") 

nrow(de %>% filter(abs(logFC) > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC < -0.58 & de$adj.P.Val < 0.05))

Volcano_Plot_Function(de, 
                      title = "(A) Volcano Plot for Node 4 versus Node 3", 
                      save_path_name = "../../figures/Chapter_4/Volcano_3.png", 
                      width = 9, height = 5, x_1 = -1, x_2 = 1)

## Volcano Plot 4 
data <- Design_List[[4]]
de2 <- model.matrix(~ 0 + data$Node)
colnames(de2) <- c("Node3", "Node4", "Node6", "Node7")

Genes_List[[4]] <- as.data.frame(Genes_List[[4]]) %>% 
    select(which(colnames(as.data.frame(Genes_List[[4]])) %in% data$PATIENT_ID))

aw <- arrayWeights(Genes_List[[4]], method = "reml", design = de2)

fit <- lmFit(Genes_List[[4]], design=de2, weigths =aw)

contrast.matrix <- makeContrasts(
    BVsT = Node7 - Node6,
    levels=colnames(de2))

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

de <- topTable(fit2, n=nrow(fit2), adjust.method = "BH") 

nrow(de %>% filter(abs(logFC) > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC > 0.58 & de$adj.P.Val < 0.05))
nrow(de %>% filter(logFC < -0.58 & de$adj.P.Val < 0.05))

Volcano_Plot_Function(de, 
                      title = "(B) Volcano Plot for Node 7 versus Node 6", 
                      save_path_name = "../../figures/Chapter_4/Volcano_4.png", 
                      width = 9, height = 5, x_1 = -1.1, x_2 = 1.1)