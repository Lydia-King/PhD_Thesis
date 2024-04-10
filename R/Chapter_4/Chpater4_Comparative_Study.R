# Chapter 4: Comparative Study

## Load Libraries 
library(tidyverse)
library(operator.tools)
library(ggVennDiagram)
library(ggplot2)
library(sf)
library(genefu)
library(readxl)
library(geneSynonym)
library(limma)

## Load up Data (for comparison)
GE_Loc <- read.delim("../../data/Processed_Data/GE_Common_Genes.txt", sep = "\t")
# CNA_Loc <- read.delim("../../data/Processed_Data/data_CNA_Loc_hg18.txt", sep="\t")
GE_Raw <- read.delim("../../data/METABRIC_2023/data_mrna_illumina_microarray.txt", sep = "\t")

## Oncotype dx 
data("sig.oncotypedx")
OncoType_Genes <- sig.oncotypedx$symbol

length(OncoType_Genes[which(OncoType_Genes %!in% GE_Loc$Hugo_Symbol)])
OncoType_Genes[which(OncoType_Genes %!in% GE_Loc$Hugo_Symbol)]
OncoType_Genes[OncoType_Genes == "CTSL2"] <- "CTSV"
length(OncoType_Genes[which(OncoType_Genes %in% GE_Loc$Hugo_Symbol)]) # 21/21

## BCI 
BCI <- c("HOXB13", "IL17RB", "BUB1B", "CENPA", "NEK2", "RACGAP1", 
         "RRM2", "ACTB", "HMBS", "SDHA", "UBC")

length(BCI[which(BCI %in% GE_Loc$Hugo_Symbol)]) # 11/11

## PAM50 
data(pam50)
Data_PAM <- pam50[["centroids.map"]]

length(Data_PAM$probe[which(Data_PAM$probe %in% GE_Loc$Hugo_Symbol)]) # 47

Data_PAM$probe[which(Data_PAM$probe %!in% GE_Loc$Hugo_Symbol)]

Data_PAM[Data_PAM$probe == "ORC6L", 1] <- "ORC6"
Data_PAM[Data_PAM$probe == "CDCA1", 1] <- "NUF2"
Data_PAM[Data_PAM$probe == "KNTC2", 1] <- "NDC80"

length(Data_PAM$probe[which(Data_PAM$probe %in% GE_Loc$Hugo_Symbol)]) # 50/50

### Housekeeping genes (8)
sum(c("ACTB", "GUSB", "MRPL19", "PSMC4", "PUM1", "RPLP0", "SF3A1", "TFRC") %in% GE_Loc$Hugo_Symbol)

PAM50_Genes <- c(Data_PAM$probe, "ACTB", "GUSB", "MRPL19", "PSMC4", "PUM1", "RPLP0", "SF3A1", "TFRC")

## Mammaprint
data("sig.gene70")

MP_Genes_Sig <- sig.gene70
MP_Genes_Sig <- MP_Genes_Sig[,c(1,4,5,6,8)]

My_MP_Genes <- c("ALDH4A1", "FGF18", "BBC3", "EBF4", "SCUBE2", "RUNDC1", "WISP1", 
                 "GSTM3", "ZNF385B", "RTN4RL1", "ECI2", "TGFB3", "STK32B", "ECI2", 
                 "MS4A7", "AP2B1", "DHX58", "TMEM74B", "ESM1", "CCNE2", "EGLN1", "CENPA", "PRC1", 
                 "PALM2AKAP2", "NMU", "IGFBP5", "PITRM1", "HRASLS", "IGFBP5",  "MSANTD3", 
                 "MCM6", "CDCA7", "RFC4", "ORC6", "SLC2A3", "GPR126", "DCK", "DTL", "COL4A2", "MELK", 
                 "MTDH", "UCHL5", "RAB6B", "GPR180", "LPCAT1", "SERF1A", "CDC42BPA", "NDC80", 
                 "GMPS", "ECT2", "MMP9", "OXCT1", "GNAZ", "FLT1", "EXT1", "CMC2", "DIAPH3", 
                 "DIAPH3", "QSOX2", "NUSAP1", "DIAPH3", "TSPYL5", "AA555029_RC", "RASSF7",  "LIN9",  
                 "RECQL5", "LOC100288906", "LOC730018", "LOC100131053", "KDM7A")

length(MP_Genes_Sig$HUGO.gene.symbol[which(MP_Genes_Sig$HUGO.gene.symbol %in% GE_Loc$Hugo_Symbol)])

length(My_MP_Genes[which(My_MP_Genes %in% GE_Loc$Hugo_Symbol)])

length(unique(My_MP_Genes))
My_MP_Genes[which(My_MP_Genes %!in% GE_Loc$Hugo_Symbol)]
My_MP_Genes[which(duplicated(My_MP_Genes))] # 3 duplicated = 4

## IntClust 
IC <- read_xls("../../data/IntClust/table_S30.xls")

## Top 1000 
IC_Sort <- IC[order(IC$All_ANOVA.pval.adj),]

Top_1000_Genes <- IC_Sort[1:1166, ]
Top_1000_Genes_NoDup <- Top_1000_Genes %>% distinct(Gene, .keep_all = T)
Top_1000_Genes_NoDup$Gene <- gsub("ORF", "orf", Top_1000_Genes_NoDup$Gene)

sum(Top_1000_Genes_NoDup$Gene %in% GE_Loc$Hugo_Symbol) # 849

# What's missing?
Missing <- Top_1000_Genes_NoDup$Gene[which(Top_1000_Genes_NoDup$Gene %!in% GE_Loc$Hugo_Symbol)] # 151

Old_Gene_Name <- c()
New_Gene_Name <- c()

for(i in 1:length(Missing)){
    
    syn <- geneSynonym(c(Missing[i]), tax = 9606)
    
    Old_Gene_Name <- c(Old_Gene_Name, Missing[i])
    p <- paste(syn[[1]][[1]][which(syn[[1]][[1]] %in% GE_Loc$Hugo_Symbol)], sep="|")
    New_Gene_Name <- c(New_Gene_Name, ifelse(is_empty(p), NA, p))
}

Gene_Swap <- data.frame(Old_Gene_Name, New_Gene_Name)

for(j in 1:length(Missing)){
    if(!is.na(New_Gene_Name)[j]){
        Top_1000_Genes_NoDup[Top_1000_Genes_NoDup$Gene == Old_Gene_Name[j], 3] <- New_Gene_Name[j] 
    }
}

sum(Top_1000_Genes_NoDup$Gene %in% GE_Loc$Hugo_Symbol) # 959

Missing2 <- Top_1000_Genes_NoDup$Gene[which(Top_1000_Genes_NoDup$Gene %!in% GE_Loc$Hugo_Symbol)] #41

## Differentially expressed genes 
load('../../data/Chapter_4/CNA_State_Environment.RData')

## Pval + LogFC
Gain_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 7) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Gain_Genes_Small$Genes))) 

Amp_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 4) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Amp_Genes_Small$Genes))) 

HetDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 9) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(HetDel_Genes_Small$Genes))) 

HomoDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 10) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(HomoDel_Genes_Small$Genes))) 

Amp_Genes <- rownames(Amp_Both) 
Gain_Genes <- rownames(Gain_Both) 
HetDel_Genes <- rownames(HetDel_Both) 
HomoDel_Genes <- rownames(HomoDel_Both) 

Genes_2022_CNA_5_Both_Complete<- unique(c(rownames(Amp_Both), rownames(HetDel_Both), rownames(Gain_Both), rownames(HomoDel_Both)))

## Pval Only 
Gain_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 7) %>% filter(adj.P.Val < 0.05) %>% 
    filter(rownames(.) %!in% c(unique(Gain_Genes_Small$Genes))) 
Amp_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 4) %>% filter(adj.P.Val < 0.05) %>% 
    filter(rownames(.) %!in% c(unique(Amp_Genes_Small$Genes))) 

HetDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 9) %>% filter(adj.P.Val < 0.05) %>% 
    filter(rownames(.) %!in% c(unique(HetDel_Genes_Small$Genes))) 
HomoDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 10) %>% filter(adj.P.Val < 0.05) %>% 
    filter(rownames(.) %!in% c(unique(HomoDel_Genes_Small$Genes))) 

Amp_Genes <- rownames(Amp_Both) 
Gain_Genes <- rownames(Gain_Both) 
HetDel_Genes <- rownames(HetDel_Both) 
HomoDel_Genes <- rownames(HomoDel_Both) 

Genes_2022_CNA_5_Both_Complete_Pval <- unique(c(rownames(Amp_Both), rownames(HetDel_Both), rownames(Gain_Both), rownames(HomoDel_Both)))


sum(Top_1000_Genes_NoDup$Gene %in%  Genes_2022_CNA_5_Both_Complete_Pval)
# Plot Overlap in Genesets 
## Setup 
plot_venn <- function (x, show_intersect, set_color, set_size, label, label_geom, 
                       label_alpha, label_color, label_size, label_percent_digit, 
                       label_txtWidth, edge_lty, edge_size, ...)  {
    venn <- Venn(x)
    data <- process_data(venn)
    p <- ggplot() + geom_sf(aes_string(fill = "count"), data = data@region) + 
        geom_sf(aes_string(color = "name"), data = data@setEdge, 
                show.legend = F, lty = edge_lty, size = edge_size, color = set_color) + 
        geom_sf_text(aes_string(label = "name"), data = data@setLabel, 
                     size = set_size, color = set_color) + theme_void()
    if (label != "none" & show_intersect == FALSE) {
        region_label <- data@region %>% dplyr::filter(.data$component == 
                                                          "region") %>% dplyr::mutate(percent = paste(round(.data$count * 
                                                                                                                100/sum(.data$count), digits = label_percent_digit), 
                                                                                                      "%", sep = "")) %>% dplyr::mutate(both = paste(.data$count, 
                                                                                                                                                     paste0("(", .data$percent, ")"), sep = "\n"))
        if (label_geom == "label") {
            p <- p + geom_sf_label(aes_string(label = label), 
                                   data = region_label, alpha = label_alpha, color = label_color, 
                                   size = label_size, lineheight = 0.85, label.size = NA)
        }
        if (label_geom == "text") {
            p <- p + geom_sf_text(aes_string(label = label), 
                                  data = region_label, alpha = label_alpha, color = label_color, 
                                  size = label_size, lineheight = 0.85)
        }
    }
    if (show_intersect == TRUE) {
        items <- data@region %>% dplyr::rowwise() %>% dplyr::mutate(text = stringr::str_wrap(paste0(.data$item, 
                                                                                                    collapse = " "), width = label_txtWidth)) %>% sf::st_as_sf()
        label_coord = sf::st_centroid(items$geometry) %>% sf::st_coordinates()
        p <- ggplot(items) + geom_sf(aes_string(fill = "count")) + 
            geom_sf_text(aes_string(label = "name"), data = data@setLabel, 
                         inherit.aes = F) + geom_text(aes_string(label = "count", 
                                                                 text = "text"), x = label_coord[, 1], y = label_coord[, 
                                                                                                                       2], show.legend = FALSE) + theme_void()
        ax <- list(showline = FALSE)
        p <- plotly::ggplotly(p, tooltip = c("text")) %>% plotly::layout(xaxis = ax, 
                                                                         yaxis = ax)
    }
    p
}

## Replace the plot_venn function with the modified version
assignInNamespace(x="plot_venn", value=plot_venn, ns="ggVennDiagram") 

ggVennDiagram(x = list(OncoType_Genes, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("OncoType \n DX ","(", length(OncoType_Genes),")"), 
                                                                                           paste0("DEGs Pval & LogFC \n","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_Onco_Ours.png", last_plot(), width = 10, height = 8)

ggVennDiagram(x = list(My_MP_Genes, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("MammaPrint ","(", length(My_MP_Genes),") \n unique (66)"), 
                                                                                       paste0("DEGs Pval & LogFC \n","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_Mamma_Ours.png", last_plot(), width = 10, height = 8)

ggVennDiagram(x = list(PAM50_Genes, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("Prosigna/PAM50 ","(", length(PAM50_Genes),")"), 
                                                                                    paste0("DEGs Pval & LogFC \n","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_Prosig_PAM5_Ours.png", last_plot(), width = 10, height = 8)

ggVennDiagram(x = list(BCI, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("BCI ","(", length(BCI),")"), 
                                                                                paste0("DEGs Pval & LogFC \n","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_BCI_Ours.png", last_plot(), width = 10, height = 8)

Top_1000_Genes_NoDup <- Top_1000_Genes_NoDup %>% distinct(Gene, .keep_all = T)
Top_1000_Genes_NoDup$Gene[which(duplicated(Top_1000_Genes_NoDup$Gene))]

ggVennDiagram(x = list(Top_1000_Genes_NoDup$Gene, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("IntClust ","(", length(Top_1000_Genes_NoDup$Gene),")"), 
                                                                                                           paste0("DEGs Pval & LogFC \n","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_IC_Ours.png", last_plot(), width = 10, height = 8)

# Check 1,000  
Gain_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 7) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Gain_Genes_Small$Genes))) 

Amp_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 4) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Amp_Genes_Small$Genes))) 

HetDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 9) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(HetDel_Genes_Small$Genes))) 

HomoDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 10) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(HomoDel_Genes_Small$Genes))) 

TopTable_1000 <- rbind.data.frame(Gain_Both, Amp_Both, HomoDel_Both, HetDel_Both)
TopTable_1000 <- TopTable_1000 %>% arrange(adj.P.Val)  

TopTable_1000 <- TopTable_1000[1:1000, ]
TopTable_1000 <- TopTable_1000 %>% distinct(rownames(.), .keep_all = T)

Genes_2022_CNA_5_Both_Complete <- rownames(TopTable_1000)

Top_1000_Genes_NoDup <- Top_1000_Genes_NoDup %>% distinct(Gene, .keep_all = T)
Top_1000_Genes_NoDup$Gene[which(duplicated(Top_1000_Genes_NoDup$Gene))]

ggVennDiagram(x = list(Top_1000_Genes_NoDup$Gene, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("IntClust ","(", length(Top_1000_Genes_NoDup$Gene),")"), 
                                                                                                      paste0("DEGs Pval & LogFC \n","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_IC_Ours_1000.png", last_plot(), width = 10, height = 8)

## Overlap between both 
Del_Both <- topTable(eBayes_CNA_3_NInt_2022, n=nrow(eBayes_CNA_3_NInt_2022), adjust.method = "BH", coef = 3) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Del_3_Small$Genes))) 

Amp_Both <- topTable(eBayes_CNA_3_NInt_2022, n=nrow(eBayes_CNA_3_NInt_2022), adjust.method = "BH", coef = 2) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Amp_3_Small$Genes))) 

Amp_Genes <- rownames(Amp_Both) 
Del_Genes <- rownames(Del_Both) 

Genes_2022_CNA_3_Both_Complete<- unique(c(rownames(Amp_Both), rownames(Del_Both)))

Gain_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 7) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Gain_Genes_Small$Genes))) 

Amp_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 4) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(Amp_Genes_Small$Genes))) 

HetDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 9) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(HetDel_Genes_Small$Genes))) 

HomoDel_Both <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 10) %>% 
    filter(abs(logFC) > 0.58, adj.P.Val < 0.05) %>% filter(rownames(.) %!in% c(unique(HomoDel_Genes_Small$Genes))) 

Amp_Genes <- rownames(Amp_Both) 
Gain_Genes <- rownames(Gain_Both) 
HetDel_Genes <- rownames(HetDel_Both) 
HomoDel_Genes <- rownames(HomoDel_Both) 

Genes_2022_CNA_5_Both_Complete<- unique(c(rownames(Amp_Both), rownames(HetDel_Both), rownames(Gain_Both), rownames(HomoDel_Both)))


ggVennDiagram(x = list(Genes_2022_CNA_3_Both_Complete, Genes_2022_CNA_5_Both_Complete), category.names = c(paste0("ModLim3 DEGs Pval \n & LogFC ","(", length(Genes_2022_CNA_3_Both_Complete),")"), 
                                                                                                           paste0("ModLim5 DEGs Pval & \n LogFC ","(", length(Genes_2022_CNA_5_Both_Complete),")")), 
              label = "count", label_size = 7, set_size = 7, set_color = c("navy","purple"))  + scale_fill_gradient(low="#002D74", high = "#942092")  +
    theme(legend.key.height= unit(1, 'cm'),
          legend.key.width= unit(1.5, 'cm')) +
    theme(legend.text = element_text(size=16))  +
    theme(legend.title = element_text(size=16)) + scale_x_continuous(expand = expansion(mult = .1)) + theme(legend.position = "bottom")

ggsave("../../figures/Chapter_4/Venn_3State_5State.png", last_plot(), width = 10, height = 8)

## Appendix D
library(xtable)
print(xtable(Top_1000_Genes_NoDup %>% filter(Gene %in% Missing2) %>% 
                 dplyr::select(ProbeID, Gene, Gene_description) %>% 
                 distinct(Gene, .keep_all = T)), include.rownames = F)