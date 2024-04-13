# Appendix B - OS Survival Trees

# Global CNA Metric Survival Trees

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
    return(data[completeVec, ])
}

## Load Data
CNA_Score_Metrics_All <- read.delim("../../data/Processed_Data/CNA_Score_All.txt", sep = "\t")
CNA_Burden_Metrics_All <- read.delim("../../data/Processed_Data/CNA_Burden_All.txt", sep = "\t")

Clinical_Sample <-
    read.delim(
        "../../data/Processed_Data/Chapter3_METABRIC_PAM50_Data.txt",
        sep = "\t",
        na.strings = c("", " ", "NA"),
    )

Clinical_Sample$INTCLUST <- 
    factor(Clinical_Sample$INTCLUST, 
           levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10"))

Clinical_Sample$PAM50<-
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

CNA_Clin_Data <-
    merge(CNA_Score_Metrics_All, Clinical_Sample_MC, by = "PATIENT_ID")

## Setup
CNA_Clin_Data$CNA_Score <- CNA_Clin_Data$CNA_Score_All
CNA_Clin_Data$Deletion_Score <- abs(CNA_Clin_Data$Del_Score_All)
CNA_Clin_Data$Amplification_Score <- CNA_Clin_Data$Amp_Score_All
CNA_Clin_Data$Difference_Score <- CNA_Clin_Data$Difference_Score_All
CNA_Clin_Data$Percentage_CNAS_Amplified <- CNA_Clin_Data$Percentage_Score_Amp_All
CNA_Clin_Data$Percentage_CNAS_Deleted <- CNA_Clin_Data$Percentage_Score_Del_All

## Naming Conventions (For for loop)
Status <- c("OS", "FiveYearOS", "TenYearOS")
StatusHeat <- c("^OS$", "^FiveYearOS$", "^TenYearOS$")
Time <- c("OS_MONTHS", "FiveYearTimeOS", "TenYearTimeOS")
Label <- c("OS", "5-year OS", "10-year OS")
Classification <- list("PAM50", "INTCLUST")
Classification_Name <- c("PAM50", "INTCLUST")

setwd("../../figures/Appendices/Appendix_B/")

Fit_Score_Rpart <- list()
Fit_Score_Ctree <- list()
Data_Score_Rpart <- list()
Data_Score_Ctree <- list()

set.seed(189)

## Survival Trees - Individual Classification
for (j in 1:length(Classification)) {
    for (i in 1) {
        Names_All <- c(Classification[[j]])
        
        Formula <-  noquote(str_c(c(Names_All), collapse = ' + '))
        
        CNA_Clin_Data_Filter <-
            CNA_Clin_Data %>% filter(OS_MONTHS > 0)
        
        CNA_Clin_Data_Filter <-
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
                "Ind_PartyKit_Survival_Score_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        png(
            filename,
            width = 10,
            height =  8,
            units = "in",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 18),
                ep_args = list(justmin = 12)
            )
        myplot
        dev.off()
        
        CNA_Clin_Data <- completeFun(CNA_Clin_Data, Status[i])
        
        filename <-
            paste(
                "Ind_Ctree_Survival_Score_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        png(
            filename,
            width = 14,
            height =  9,
            units = "in",
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
        
        plot(tree, gp = gpar(fontsize = 21))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Score_Ctree[[filename]] <- data
    }
}

## Survival Trees - CNA Score in Combination with Classifications
Fit_Score_Rpart <- list()
Fit_Score_Ctree <- list()
Data_Score_Rpart <- list()
Data_Score_Ctree <- list()

for(j in 1:length(Classification)) {
    # PAM50 First, then INTCLUST
    for (i in 1:length(Label)) {
        
        Names_All <-
            c(
                "CNA_Score",
                "Amplification_Score",
                "Deletion_Score",
                "Difference_Score",
                "Percentage_CNAS_Amplified",
                "Percentage_CNAS_Deleted",
                Classification[[j]]
            )
        
        Formula <-  noquote(str_c(c(Names_All), collapse = ' + '))
        CNA_Clin_Data_Filter <-
            CNA_Clin_Data %>% filter(OS_MONTHS > 0)
        CNA_Clin_Data_Filter <-
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
                "PartyKit_Survival_Score_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        Fit_Score_Rpart[[filename]] <- pfit
        png(
            filename,
            width = 20,
            height =  12,
            units = "cm",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 10),
                ep_args = list(justmin = 10)
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
        Data_Score_Rpart[[filename]] <- data
        
        ## PartyKit KM Curve
        data1 <-
            CNA_Clin_Data_Filter %>% mutate(Node = as.factor(pfit$where))
        
        data1 <-
            data.frame(Time = data1[, c(Time[i])],
                       Node = data1[, c("Node")],
                       Cen  =  data1[, c(Status[i])])
        
        CNA_Clin_Data <- completeFun(CNA_Clin_Data, Status[i])
        
        filename <-
            paste("Ctree_Survival_Score_",
                  Status[i],
                  "_",
                  Classification_Name[j],
                  ".png",
                  sep = "")
        png(
            filename,
            width = 20,
            height =  12,
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
        Fit_Score_Ctree[[filename]] <- tree
        plot(tree, gp = gpar(fontsize = 10))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Score_Ctree[[filename]] <- data
    }
}

## Survival Trees - CNA Burden in Combination with Classifications
## Setup
CNA_Clin_Data <- merge(CNA_Burden_Metrics_All, 
                       Clinical_Sample_MC, by="PATIENT_ID")

## Check how many missing values
sapply(CNA_Clin_Data, function(x) sum(is.na(x)))

CNA_Clin_Data$CNA_Burden <- CNA_Clin_Data$CNA_Burden_All
CNA_Clin_Data$Deletion_Burden <- abs(CNA_Clin_Data$Del_Burden_All)
CNA_Clin_Data$Amplification_Burden <- CNA_Clin_Data$Amp_Burden_All
CNA_Clin_Data$Difference_Burden <- CNA_Clin_Data$Difference_Burden_All
CNA_Clin_Data$Percentage_CNAS_Amplified <- CNA_Clin_Data$Percentage_Burden_Amp_All
CNA_Clin_Data$Percentage_CNAS_Deleted <- CNA_Clin_Data$Percentage_Burden_Del_All

Fit_Burden_Rpart <- list()
Data_Burden_Rpart <- list()
Fit_Burden_Ctree <- list()
Data_Burden_Ctree <- list()

for (j in 1:length(Classification)) {
    
    for (i in 1:length(Label)) {
        
        Names_All <-
            c(
                "CNA_Burden",
                "Amplification_Burden",
                "Deletion_Burden",
                "Difference_Burden",
                "Percentage_CNAS_Amplified",
                "Percentage_CNAS_Deleted",
                Classification[[j]]
            )
        
        Formula <-  noquote(str_c(c(Names_All), collapse = ' + '))
        
        CNA_Clin_Data_Filter <-
            CNA_Clin_Data %>% filter(OS_MONTHS > 0)
        
        CNA_Clin_Data_Filter <-
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
                "PartyKit_Survival_Burden_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        Fit_Burden_Rpart[[filename]] <- pfit
        png(
            filename,
            width = 20,
            height =  12,
            units = "cm",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 10),
                ep_args = list(justmin = 10)
            )
        myplot
        dev.off()
        
        ## PartyKit KM Curve
        data1 <-
            CNA_Clin_Data_Filter %>% mutate(Node = as.factor(pfit$where))
        data <- data1 %>% select(PATIENT_ID, Node)
        Data_Burden_Rpart[[filename]] <- data
        
        data1 <-
            data.frame(Time = data1[, c(Time[i])],
                       Node = data1[, c("Node")],
                       Cen  =  data1[, c(Status[i])])
        
        CNA_Clin_Data <- completeFun(CNA_Clin_Data, Status[i])
        
        filename <-
            paste("Ctree_Survival_Burden_",
                  Status[i],
                  "_",
                  Classification_Name[j],
                  ".png",
                  sep = "")
        png(
            filename,
            width = 20,
            height =  12,
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
        plot(tree, gp = gpar(fontsize = 10))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Burden_Ctree[[filename]] <- data
    }
}

# Global CNA Metric Survival Trees with Clinical Variables
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

## Survival Trees - CNA Burden in Combination with Classifications and Clinical
CNA_Clin_Data <- merge(CNA_Burden_Metrics_All, Clinical_Sample_MC, by="PATIENT_ID")

CNA_Clin_Data$CNA_Burden <- CNA_Clin_Data$CNA_Burden_All
CNA_Clin_Data$Deletion_Burden <- abs(CNA_Clin_Data$Del_Burden_All)
CNA_Clin_Data$Amplification_Burden <- CNA_Clin_Data$Amp_Burden_All
CNA_Clin_Data$Difference_Burden <- CNA_Clin_Data$Difference_Burden_All
CNA_Clin_Data$Percentage_CNAS_Amplified <- CNA_Clin_Data$Percentage_Burden_Amp_All
CNA_Clin_Data$Percentage_CNAS_Deleted <- CNA_Clin_Data$Percentage_Burden_Del_All

Fit_Burden_Rpart <- list()
Data_Burden_Rpart <- list()
Fit_Burden_Ctree <- list()
Data_Burden_Ctree <- list()

set.seed(189)

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


for (j in 1:length(Classification)) {
    
    for (i in 1:length(Label)) {
        
        Names_All <-
            c(
                "CNA_Burden",
                "Amplification_Burden",
                "Deletion_Burden",
                "Difference_Burden",
                "Percentage_CNAS_Amplified",
                "Percentage_CNAS_Deleted",
                Classification[[j]],
                "Positive_Lymph_Nodes",
                "NPI",
                "PR_Status",
                "ER_Status",
                "HER2_Status",
                "Age",
                "Tumour_Size",
                "Clinical_Stage",
                "Histologic_Grade",
                "Cancer_Type_Detailed"
            )
        
        Formula <-  noquote(str_c(c(Names_All), collapse = ' + '))
        
        CNA_Clin_Data_Filter <-
            CNA_Clin_Data %>% filter(OS_MONTHS > 0)
        CNA_Clin_Data_Filter <-
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
                "Clin_PartyKit_Survival_Burden_",
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
                ep_args = list(justmin = 6)
            )
        myplot
        dev.off()
        
        ## PartyKit KM Curve
        data1 <-
            CNA_Clin_Data_Filter %>% mutate(Node = as.factor(pfit$where))
        data <- data1 %>% select(PATIENT_ID, Node)
        Data_Burden_Rpart[[filename]] <- data
        
        data1 <-
            data.frame(Time = data1[, c(Time[i])],
                       Node = data1[, c("Node")],
                       Cen  =  data1[, c(Status[i])])
        
        CNA_Clin_Data <- completeFun(CNA_Clin_Data, Status[i])
        
        filename <-
            paste("Clin_Ctree_Survival_Burden_",
                  Status[i],
                  "_",
                  Classification_Name[j],
                  ".png",
                  sep = "")
        png(
            filename,
            width = 22,
            height =  12,
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
                    minprob = 0.02,
                    maxdepth = 4
                )
            )
        Fit_Burden_Ctree[[filename]] <- tree
        plot(tree, gp = gpar(fontsize = 8))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Burden_Ctree[[filename]] <- data
    }
}

# Chapter 3: ChrArm CNA Metric Survival Trees 
## Load Data
setwd("../../../R/Appendices/")

CNA_Arm_Score_Metrics <-
    read.delim("../../data/Processed_Data/Chr_Arm_CNA_Score_hg18.txt", sep = "\t")

CNA_Arm_Burden_Metrics <-
    read.delim("../../data/Processed_Data/Chr_Arm_CNA_Burden_hg18.txt", sep = "\t")

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

CNA_Clin_Data <-
    merge(CNA_Arm_Score_Metrics, Clinical_Sample_MC, by = "PATIENT_ID")

colnames(CNA_Clin_Data) <-
    gsub("_All", "", colnames(CNA_Clin_Data))

colnames(CNA_Clin_Data) <-
    gsub("Per_Amp", "PerAmp", colnames(CNA_Clin_Data))

colnames(CNA_Clin_Data) <-
    gsub("Per_Del", "PerDel", colnames(CNA_Clin_Data))

## Survival Trees - CNA Score in Combination with Classifications

setwd("../../figures/Appendices/Appendix_B/")
set.seed(123)

Fit_Score_Rpart <- list()
Fit_Score_Ctree <- list()
Data_Score_Rpart <- list()
Data_Score_Ctree <- list()

for (j in 1:length(Classification)) {
    
    for (i in 1:length(Label)) {
        
        Names_All <-
            CNA_Clin_Data %>% 
            dplyr::select(!grep("CCA|PATIENT_ID", names(CNA_Arm_Score_Metrics)))
        
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
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
                "PA_PartyKit_Survival_Score_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        Fit_Score_Rpart[[filename]] <- pfit
        png(
            filename,
            width = 20,
            height =  12,
            units = "cm",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 10),
                ep_args = list(justmin = 10)
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
        Data_Score_Rpart[[filename]] <- data
        
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
                "PA_Ctree_Survival_Score_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        png(
            filename,
            width = 20,
            height =  12,
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
        Fit_Score_Ctree[[filename]] <- tree
        plot(tree, gp = gpar(fontsize = 10))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Score_Ctree[[filename]] <- data
    }
}

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
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
        png(
            filename,
            width = 20,
            height =  12,
            units = "cm",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 10),
                ep_args = list(justmin = 10)
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
                "PA_Ctree_Survival_Burden_",
                Status[i],
                "_",
                Classification_Name[j],
                ".png",
                sep = ""
            )
        png(
            filename,
            width = 20,
            height =  12,
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
        plot(tree, gp = gpar(fontsize = 10))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) # Get node information
        
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Burden_Ctree[[filename]] <- data
    }
}

# ChrArm CNA Metric Survival Trees with Clinical Variables
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

set.seed(123)

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
            completeFun(CNA_Clin_Data_Filter, "OS")
        
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
            width = 20,
            height =  12,
            units = "cm",
            res = 400
        )
        myplot <-
            plot(
                as.party(pfit),
                gp = gpar(fontsize = 10),
                ep_args = list(justmin = 10)
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
            width = 23,
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
        plot(tree, gp = gpar(fontsize = 9))
        dev.off()
        
        data2 <-
            CNA_Clin_Data %>% 
            mutate(Node = as.factor(predict(tree, type = "node"))) 
        # Get node information
        data <- data2 %>% select(PATIENT_ID, Node)
        Data_Burden_Ctree[[filename]] <- data
    }
}