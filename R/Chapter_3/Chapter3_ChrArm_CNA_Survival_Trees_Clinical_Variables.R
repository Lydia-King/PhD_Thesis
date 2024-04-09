# Chapter 3: ChrArm CNA Metric Survival Trees with Clinical Variables

## Load libraries
library(stringr)
library(tidyverse)
library(partykit)
library(rpart)
library(rpart.plot)
library(survival)
library(survminer)

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
        TenYearDSS,
        "LYMPH_NODES_EXAMINED_POSITIVE",
        "NPI",
        "PR_STATUS",
        "ER_STATUS",
        "HER2_STATUS",
        "AGE_AT_DIAGNOSIS",
        "TUMOR_SIZE",
        "TUMOR_STAGE",
        "GRADE",
        "CANCER_TYPE_DETAILED"
    )

Clinical_Sample_MC <-
    Clinical_Sample_MC %>% mutate_if( ~ any(is.character(.x)),  ~ as.factor(.))

## Survival Trees - CNA Score in Combination with Classifications with Clinical Vars
## Setup
## Naming Conventions (For for loop)
Status <- c("DSS", "FiveYearDSS", "TenYearDSS")
StatusHeat <- c("^DSS$", "^FiveYearDSS$", "^TenYearDSS$")
Time <- c("OS_MONTHS", "FiveYearTimeDSS", "TenYearTimeDSS")
Label <- c("DSS", "5-year DSS", "10-year DSS")
Classification <- list("PAM50", "INTCLUST")
Classification_Name <- c("PAM50", "INTCLUST")

setwd("../../figures/Chapter_3/")
set.seed(123)

Fit_Burden_Rpart <- list()
Data_Burden_Rpart <- list()
Fit_Burden_Ctree <- list()
Data_Burden_Ctree <- list()

## Survival Trees - CNA Burden in Combination with Classifications 
## Setup
CNA_Clin_Data <-
    merge(CNA_Arm_Burden_Metrics, Clinical_Sample_MC, by = "PATIENT_ID")

colnames(CNA_Clin_Data) <-
    gsub("_All", "", colnames(CNA_Clin_Data))

colnames(CNA_Clin_Data) <-
    gsub("Per_Amp", "PerAmp", colnames(CNA_Clin_Data))

colnames(CNA_Clin_Data) <-
    gsub("Per_Del", "PerDel", colnames(CNA_Clin_Data))

Fit_Burden_Rpart <- list()
Data_Burden_Rpart <- list()
Fit_Burden_Ctree <- list()
Data_Burden_Ctree <- list()

set.seed(123)

CNA_Clin_Data <- CNA_Clin_Data %>% dplyr::rename(
    "Histologic_Grade" = "GRADE",
    "HER2_Status" =  "HER2_STATUS",
    "PR_Status" = "PR_STATUS",
    "Clinical_Stage" = "TUMOR_STAGE",
    "Age" = "AGE_AT_DIAGNOSIS",
    "Positive_Lymph_Nodes" = "LYMPH_NODES_EXAMINED_POSITIVE",
    "NPI" = "NPI",
    "Tumour_Size" = "TUMOR_SIZE",
    "ER_Status" = "ER_STATUS",
    "Cancer_Type_Detailed" =  "CANCER_TYPE_DETAILED"
)

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
                "Clin_PA_PartyKit_Survival_Burden_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        Fit_Burden_Rpart[[filename]] <- pfit
        png(
            filename,
            width = 22,
            height =  12,
            units = "cm",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 10),
                edge_panel = edge_simple(as.party(pfit), digits = 2, abbreviate = FALSE, justmin = 10)
            )
        myplot
        dev.off()
        
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
                "Clin_PA_Ctree_Survival_Burden_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        png(
            filename,
            width = 25,
            height =  13,
            units = "cm",
            res = 400
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
        plot(tree, gp = gpar(fontsize = 9), edge_panel = edge_simple(tree, digits = 2, justmin = 10))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) 
        # Get node information
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Burden_Ctree[[filename]] <- data
    }
}