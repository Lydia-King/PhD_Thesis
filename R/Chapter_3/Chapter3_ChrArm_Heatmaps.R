# Chapter 3: Heatmaps of Chromosome 3p, 18q and 11p (Pre-selected)

## Load up libraries 
library(stringr)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(partykit)
library(survival)
library(survminer)
library(RColorBrewer)

## Complete Function
completeFun <- function(data, desiredCols) {
    completeVec <- complete.cases(data[, desiredCols])
    return(data[completeVec,])
}

## Load Data
CNA_Arm_Score_Metrics <-
    read.delim("../../data/Processed_Data/Chr_Arm_CNA_Score.txt", sep = "\t")

CNA_Arm_Burden_Metrics <-
    read.delim("../../data/Processed_Data/Chr_Arm_CNA_Burden.txt", sep = "\t")

Clinical_Sample <-
    read.delim(
        "../../data/Processed_Data/Chapter3_METABRIC_PAM50_Data.txt",
        sep = "\t",
        na.strings = c("", " ", "NA"),
    )

Clinical_Sample$INTCLUST <-
    factor(
        Clinical_Sample$INTCLUST,
        levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10")
    )

Clinical_Sample$PAM50 <-
    as.factor(as.character(Clinical_Sample$PAM50))

Clinical_Sample_MC <-
    Clinical_Sample %>% select(
        PATIENT_ID,
        PAM50,
        INTCLUST,
        OS_MONTHS,
        OS,
        DSS,
        FiveYearTimeOS,
        FiveYearOS,
        FiveYearTimeDSS,
        FiveYearDSS,
        TenYearTimeOS,
        TenYearOS,
        TenYearTimeDSS,
        TenYearDSS
    )

Clinical_Sample_MC <- Clinical_Sample_MC %>% mutate_if(~ any(is.character(.x)),~ as.factor(.)) 

CNA_Clin_Data <-
    merge(CNA_Arm_Burden_Metrics, Clinical_Sample_MC, by = "PATIENT_ID")

colnames(CNA_Clin_Data) <-
    gsub("_All", "", colnames(CNA_Clin_Data))

colnames(CNA_Clin_Data) <-
    gsub("Per_Amp", "PerAmp", colnames(CNA_Clin_Data))

colnames(CNA_Clin_Data) <-
    gsub("Per_Del", "PerDel", colnames(CNA_Clin_Data))

## Setup
StatusHeat <- c("^DSS$", "^FiveYearDSS$", "^TenYearDSS$")
Time <- c("OS_MONTHS", "FiveYearTimeDSS", "TenYearTimeDSS")
Label <- c("DSS", "5-year DSS", "10-year DSS")
Classification <- list("PAM50", "INTCLUST")
Classification_Name <- c("PAM50", "INTCLUST")
Status <- c("DSS", "FiveYearDSS", "TenYearDSS")

Fit_Burden_Rpart <- list()
Data_Burden_Rpart <- list()
Fit_Burden_Ctree <- list()
Data_Burden_Ctree <- list()

set.seed(123)

## Survival Trees 
for(j in 1:length(Classification)) {

    for (i in 1:length(Label)) {

        Names_All <-
            CNA_Clin_Data %>% 
            dplyr::select(!grep("CCA|PATIENT_ID", names(CNA_Arm_Burden_Metrics)))
        
        Names_All <-
            Names_All %>% 
            dplyr::select(!grep("OS|DSS", colnames(Names_All)))
        
        Names_All <-
            Names_All %>% 
            dplyr::select(!grep("PAM50|INTCLUST", colnames(Names_All)))
        Formula <-
            noquote(str_c(c(
                colnames(Names_All), c(Classification[[j]])
            ), collapse = ' + '))
        
        CNA_Clin_Data_Filter <-
            CNA_Clin_Data %>% filter(OS_MONTHS > 0)
        
        CNA_Clin_Data_Filter <-
            completeFun(CNA_Clin_Data_Filter, "DSS")
        
        pfit <-
            rpart(
                paste("Surv(", Time[i], ",", Status[i], ") ~ ", 
                      Formula, sep = ""),
                data = CNA_Clin_Data_Filter,
                method =  "exp",
                model = TRUE
            )
        
        ## PartyKit
        filename <-
            paste(
                "PA_PartyKit_Survival_Burden_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        
        Fit_Burden_Rpart[[filename]] <- pfit
        
        ## Store Corresponding Data
        if (is.null(pfit$na.action)) {
            data1 <-
                CNA_Clin_Data_Filter %>% mutate(Node = as.factor(pfit$where))
        } else {
            data1 <-
                CNA_Clin_Data_Filter[-c(as.numeric(pfit$na.action)), ] %>% 
                mutate(Node = as.factor(pfit$where))
        } # Get node information
        
        data <- data1 %>% select(PATIENT_ID, Node)
        Data_Burden_Rpart[[filename]] <- data
        
        ## PartyKit KM Curve
        data1 <-
            CNA_Clin_Data_Filter %>% mutate(Node = as.factor(pfit$where))
        data1 <-
            data.frame(Time = data1[, c(Time[i])],
                       Node = data1[, c("Node")],
                       Cen  =  data1[, c(Status[i])])
        
        CNA_Clin_Data <- completeFun(CNA_Clin_Data, Status[i])
        
        filename <-
            paste(
                "PA_Ctree_Survival_Burden_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        tree <-
            partykit::ctree(
                as.formula(
                    paste("Surv(", Time[i], ",", Status[i], ") ~ ", 
                          Formula, sep = "")
                ),
                data = CNA_Clin_Data,
                control = ctree_control(
                    minsplit = 30L,
                    minbucket = 30L,
                    minprob = 0.02
                )
            )
        Fit_Burden_Ctree[[filename]] <- tree
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Burden_Ctree[[filename]] <- data
    }
}

## Which Trees?
Burden_Tree_Fits <- append(Fit_Burden_Ctree, Fit_Burden_Rpart)
Fit_List_Burden <- Burden_Tree_Fits[c(1,10,8)] # 3p, 18q, 11p
List_For_Each_Tree_Burden <- c("3p", "18q", "11p")

## Fit Heatmaps 
# 3p, 11p and 18q
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}
cols_Subtype = gg_color_hue(6)
IC_Col = brewer.pal(n = 11, name = "Paired")
colours1 <-
    list(
        'PAM50' = c(
            'Basal' = cols_Subtype[1],
            'HER2' = cols_Subtype[3],
            'Claudin-low' = cols_Subtype[2],
            'Normal' = cols_Subtype[6],
            'LumA' = cols_Subtype[4],
            'LumB' = cols_Subtype[5]
        ),
        'INTCLUST' = c(
            '1' = IC_Col[1],
            '2' = IC_Col[2],
            '3' = IC_Col[3],
            '4ER+' = IC_Col[4],
            '4ER-' = IC_Col[5],
            '5' = IC_Col[6],
            '6' = IC_Col[7],
            '7' = IC_Col[8],
            '8' = IC_Col[9],
            '9' = IC_Col[10] ,
            '10' = IC_Col[11]
        )
    )

## Chromosome Arm Location for Genes
CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg18.txt",
        sep = ",",
        na.strings = c("", " ", "NA")
    )

CNA_Loc <- CNA_Loc %>% mutate(Location = paste(Chromosome, Arm, sep="")) %>%
    select("Location","txStart",7:ncol(CNA_Loc))

CNA_Metrics_Burden <-
    merge(CNA_Arm_Burden_Metrics, Clinical_Sample[, c(
        "PATIENT_ID",
        "PAM50",
        "INTCLUST",
        "OS_MONTHS",
        "OS",
        "DSS",
        "FiveYearTimeOS",
        "FiveYearOS",
        "TenYearTimeOS",
        "TenYearOS",
        "FiveYearTimeDSS",
        "FiveYearDSS",
        "TenYearTimeDSS",
        "TenYearDSS"
    )], by = "PATIENT_ID")

CNA_Metrics_Burden <-
    CNA_Metrics_Burden %>% mutate(INTCLUST = as.factor(INTCLUST), PAM50 = as.factor(PAM50))

names(CNA_Metrics_Burden) <-
    gsub(" ", "_", names(CNA_Metrics_Burden), fixed = TRUE)

# CNA Burden
Datasets_Burden_CNA_All <- list()

for(i in 1:length(Fit_List_Burden)) {
    # For each tree
    Chr_Heat <- List_For_Each_Tree_Burden[i]
    
    ## CNA Data (All)
    CNA_Loc_Chr <- CNA_Loc %>% filter(Location == Chr_Heat) 

    CNA_Loc_Chr <-
        CNA_Loc_Chr[order(CNA_Loc_Chr$txStart),] # Order based on Gene Location
    
    CNA_Loc_Chr <-
        CNA_Loc_Chr %>% select(-txStart, -Location) # Remove txStart column (only after ordering)
    
    rownames(CNA_Loc_Chr) <-
        CNA_Loc_Chr$Hugo_Symbol # Alter rownames
    
    CNA_Loc_Chr$Hugo_Symbol <- NULL
    
    CNA_T_Chr <- as.data.frame(t(CNA_Loc_Chr)) # Transpose dataframe
    
    CNA_T_Chr <-
        CNA_T_Chr %>% mutate(PATIENT_ID = rownames(CNA_T_Chr))
    
    CNA_T_Chr <-
        merge(Clinical_Sample[, c("PATIENT_ID", "PAM50", "INTCLUST")], CNA_T_Chr, by = "PATIENT_ID")
    
    ## Complete Data
    DATA <- completeFun(CNA_Metrics_Burden, c("DSS"))
    
    ## Add node information
    if (i %in% c(2, 3)) {
        data2 <-
            DATA %>% filter(OS_MONTHS > 0) %>% 
            mutate(Node = paste("Node", as.factor(Fit_List_Burden[[i]]$where), 
                                sep = " "))
        CNA_Sub <-
            merge(data2[, c("PATIENT_ID", "Node")], CNA_T_Chr, BY = "PATIENT_ID")
    } else {
        data2 <-
            DATA %>% mutate(Node = as.factor(predict(Fit_List_Burden[[i]], type = "node")))
        
        CNA_Sub <-
            merge(data2[, c("PATIENT_ID", "Node")], CNA_T_Chr, BY = "PATIENT_ID")
        
        CNA_Sub <-
            CNA_Sub %>% mutate(Node = paste("Node", Node, sep = " "))
    }
    
    ## Store in dataset list
    name <- paste0(names(Fit_List_Burden[i]), "_", Chr_Heat, "_All")
    Datasets_Burden_CNA_All[[name]] <- CNA_Sub
}

## CNA Burden Heatmaps
## List for each chromosome arm used to split patients in each tree using CNA Burden metrics 
setwd("../../figures/Chapter_3/")

## Set up lists to store across chromosome heatmaps 
CNA_Burden_Heatmaps_Fig_CNA_1 <- list()
library(ComplexHeatmap)
Letter <- c("(A)", "(B)")

## CNA Burden heatmaps
for (i in 1:length(c(Datasets_Burden_CNA_All))) {
    Chr_Heat <- List_For_Each_Tree_Burden[i]
    
    Fit_CNA_All <- Datasets_Burden_CNA_All[[i]]
    
    Fit_CNA_All <-
        Fit_CNA_All %>%  mutate(Node = factor(Node, level  = c(paste(
            "Node", 1:10, sep = " "
        ))))
    
    ## CNA All - Across Chromosome with PAM50/INTCLUST (and NAs in order for Gene Expression)
    row_ha = rowAnnotation(
        df = Fit_CNA_All[, c(3, 4)],
        simple_anno_size = unit(2, "cm"),
        col = colours1,
        gap = unit(1, 'mm')
    )
    
    ## Produce Heatmaps - Add in Common fragile site markers
    filename <-
        paste(names(Datasets_Burden_CNA_All)[[i]], "_Heatmap.png", sep = "")
    
    ht11 = Heatmap(
        as.matrix(Fit_CNA_All[, -c(1, 2, 3, 4)]),
        column_title = paste("Heatmap of CNAs across Chromosome", Chr_Heat),
        cluster_row_slices = FALSE,
        show_column_names = F,
        show_row_names = F,
        name = "CNA",
        column_title_side = "top",
        row_split = Fit_CNA_All[, c("Node")],
        row_gap = unit(5, "mm"),
        column_order = colnames(Fit_CNA_All[, -c(1, 2, 3, 4)]),
        column_names_gp = grid::gpar(fontsize = 14),
        row_names_gp = grid::gpar(fontsize = 14),
        column_title_gp = gpar(fontsize = 16),
        right_annotation = row_ha
    )
    
    ## Store heatmap in list
    CNA_Burden_Heatmaps_Fig_CNA_1[[filename]] <- ht11
    
    png(
        filename,
        width = 27,
        height =  16,
        units = "cm",
        res = 400
    )
    draw(ht11)
    dev.off()
}