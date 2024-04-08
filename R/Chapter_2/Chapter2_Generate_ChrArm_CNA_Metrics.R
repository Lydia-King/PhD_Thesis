# Chapter 2: Chromosome Arm CNA Score Metric Calculation (All Data and Complete Case)
## Load up necessary libraries
library(reshape)
library(dplyr)

## Load up data (CNA data with gene locations)
CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg18.txt",
        sep = ",",
        na.strings = c("", " ", "NA")
    )

CNA_Loc <- CNA_Loc %>% mutate(Location = paste(Chromosome, Arm, sep="")) %>%
    select(Location, 7:ncol(CNA_Loc)) %>% mutate(Location = ifelse(Location == "23p", "Xp", 
                                                                   ifelse(Location == "23q", "Xq", 
                                                                          Location)))

## Desired order of chromosome arms
variable_order <-
    c(
        "PATIENT_ID",
        "1p",
        "1q",
        "2p",
        "2q",
        "3p",
        "3q",
        "4p",
        "4q",
        "5p",
        "5q",
        "6p",
        "6q",
        "7p",
        "7q",
        "8p",
        "8q",
        "9p",
        "9q",
        "10p",
        "10q",
        "11p",
        "11q",
        "12p",
        "12q",
        "13q",
        "14q",
        "15q",
        "16p",
        "16q",
        "17p",
        "17q",
        "18p",
        "18q",
        "19p",
        "19q",
        "20p",
        "20q",
        "21p",
        "21q",
        "22q",
        "Xp",
        "Xq"
    )

# Chromosome Arm CNA Score Metrics
Per_Arm_CNA_Score <- function(Data) {
    ## Absolute CNA Score
    CNA_Score_All <-
        Data %>% group_by(Location) %>% select(2:ncol(.)) %>%
        summarise_each(funs(sum(abs(.), na.rm = TRUE)))
    CNA_Score_CCA <-
        Data %>% group_by(Location) %>% select(2:ncol(.)) %>%
        summarise_each(funs(sum(abs(.), na.rm = FALSE)))
    
    ## Amplification Score
    CNA_Score_Amp_All <- Data %>% group_by(Location) %>%
        summarise_each(funs(sum(.[. > 0], na.rm = TRUE)))
    CNA_Score_Amp_CCA <- Data %>% group_by(Location) %>%
        summarise_each(funs(sum(.[. > 0], na.rm = FALSE)))
    
    ## Deletion Score
    CNA_Score_Del_All <- Data %>% group_by(Location) %>%
        summarise_each(funs(sum(abs(.[. < 0]), na.rm = TRUE)))
    CNA_Score_Del_CCA <- Data %>% group_by(Location) %>%
        summarise_each(funs(sum(abs(.[. < 0]), na.rm = FALSE)))
    
    ## Difference Score
    CNA_Score_Diff_All <-
        CNA_Score_Amp_All[, 2:ncol(CNA_Score_Amp_All)] -
        CNA_Score_Del_All[, 2:ncol(CNA_Score_Del_All)]
    CNA_Score_Diff_All <- CNA_Score_Diff_All %>%
        mutate(Location = CNA_Score_All$Location) %>%
        select(Location, everything())
    
    CNA_Score_Diff_CCA <-
        CNA_Score_Amp_CCA[, 2:ncol(CNA_Score_Amp_CCA)] -
        CNA_Score_Del_CCA[, 2:ncol(CNA_Score_Del_CCA)]
    CNA_Score_Diff_CCA <- CNA_Score_Diff_CCA %>%
        mutate(Location = CNA_Score_CCA$Location) %>%
        select(Location, everything())
    
    ## Percentage of Alterations Amplification Score
    CNA_Score_Per_Amp_All <-
        (CNA_Score_Amp_All[, 2:ncol(CNA_Score_Amp_All)] /
             CNA_Score_All[, 2:ncol(CNA_Score_All)]) *
        100
    CNA_Score_Per_Amp_All[CNA_Score_Per_Amp_All == "NaN"] <- 0
    CNA_Score_Per_Amp_All <- CNA_Score_Per_Amp_All %>%
        mutate(Location = CNA_Score_All$Location) %>%
        select(Location, everything())
    
    CNA_Score_Per_Amp_CCA <-
        (CNA_Score_Amp_CCA[, 2:ncol(CNA_Score_Amp_CCA)] /
             CNA_Score_CCA[, 2:ncol(CNA_Score_CCA)]) *
        100
    CNA_Score_Per_Amp_CCA[CNA_Score_Per_Amp_CCA == "NaN"] <- 0
    CNA_Score_Per_Amp_CCA <- CNA_Score_Per_Amp_CCA %>%
        mutate(Location = CNA_Score_CCA$Location) %>%
        select(Location, everything())
    
    ## Percentage of Alterations Deletion Score
    CNA_Score_Per_Del_All  <-
        (CNA_Score_Del_All[, 2:ncol(CNA_Score_Del_All)] /
             CNA_Score_All[, 2:ncol(CNA_Score_All)]) *
        100
    CNA_Score_Per_Del_All[CNA_Score_Per_Del_All == "NaN"] <- 0
    CNA_Score_Per_Del_All <- CNA_Score_Per_Del_All %>%
        mutate(Location = CNA_Score_All$Location) %>%
        select(Location, everything())
    
    CNA_Score_Per_Del_CCA  <-
        (CNA_Score_Del_CCA[, 2:ncol(CNA_Score_Del_CCA)] /
             CNA_Score_CCA[, 2:ncol(CNA_Score_CCA)]) *
        100
    CNA_Score_Per_Del_CCA[CNA_Score_Per_Del_CCA == "NaN"] <- 0
    CNA_Score_Per_Del_CCA <- CNA_Score_Per_Del_CCA %>%
        mutate(Location = CNA_Score_CCA$Location) %>%
        select(Location, everything())
    
    list_Scores <-
        list(
            CNA_Score_All,
            CNA_Score_Amp_All,
            CNA_Score_Del_All,
            CNA_Score_Diff_All,
            CNA_Score_Per_Amp_All,
            CNA_Score_Per_Del_All,
            CNA_Score_CCA,
            CNA_Score_Amp_CCA,
            CNA_Score_Del_CCA,
            CNA_Score_Diff_CCA,
            CNA_Score_Per_Amp_CCA,
            CNA_Score_Per_Del_CCA
        )
    
    names_Scores <-
        c(
            "CNA_Score_All",
            "CNA_Amp_All",
            "CNA_Del_All",
            "CNA_Difference_All",
            "CNA_Per_Amp_All",
            "CNA_Per_Del_All",
            "CNA_Score_CCA",
            "CNA_Amp_CCA",
            "CNA_Del_CCA",
            "CNA_Difference_CCA",
            "CNA_Per_Amp_CCA",
            "CNA_Per_Del_CCA"
        )
    
    for (i in 1:length(list_Scores)) {
        list_Scores[[i]] <-
            as.data.frame(t(list_Scores[[i]]), stringsAsFactors = FALSE)
        names(list_Scores[[i]]) <- list_Scores[[i]][1, ]
        list_Scores[[i]] <- list_Scores[[i]][-1, ]
        list_Scores[[i]] <-
            lapply(list_Scores[[i]], type.convert, as.is = TRUE)
        list_Scores[[i]] <-
            cbind.data.frame(PATIENT_ID = colnames(Data[, -c(1)]),
                             list_Scores[[i]])
        list_Scores[[i]] <-
            list_Scores[[i]] %>% select(c(variable_order))
        
        for (j in 2:ncol(list_Scores[[i]])) {
            colnames(list_Scores[[i]])[j] <-
                as.character(paste0(names_Scores[i], "_", colnames(list_Scores[[i]])[j]))
        }
    }
    
    names(list_Scores) <- names_Scores
    return(list_Scores)
}

## Run function
CNA_Arm_Score_Metrics_List <- Per_Arm_CNA_Score(CNA_Loc)

## Merge all CNA Score metrics
CNA_Arm_Score_Metrics <-
    merge_recurse(CNA_Arm_Score_Metrics_List, by = "PATIENT_ID")

## Order columns
Order_Function <- function(Data) {
    new_order <- c()
    
    for (i in 2:ncol(Data[[1]])) {
        a <- colnames(Data[[1]][i])
        b <- colnames(Data[[7]][i])
        c <- colnames(Data[[2]][i])
        d <- colnames(Data[[8]][i])
        e <- colnames(Data[[3]][i])
        f <- colnames(Data[[9]][i])
        g <- colnames(Data[[4]][i])
        h <- colnames(Data[[10]][i])
        j <- colnames(Data[[5]][i])
        k <- colnames(Data[[11]][i])
        l <- colnames(Data[[6]][i])
        m <- colnames(Data[[12]][i])
        new_order <- c(new_order, a, c, e, g, j, l, b, d, f, h, k, m)
    }
    new_order <<- new_order
}

## Run function
Order_Function(CNA_Arm_Score_Metrics_List)

## Apply new order
CNA_Arm_Score_Metrics  <- CNA_Arm_Score_Metrics %>%
    select(c("PATIENT_ID", all_of(new_order)))

# Chromosome Arm CNA Burden
Per_Arm_CNA_Burden <- function(Data) {
    ## Absolute CNA Burden
    CNA_Burden_All <-
        Data %>% group_by(Location) %>% select(2:ncol(.)) %>%
        summarise_each(funs(sum(. != 0, na.rm = TRUE)))
    Number_Genes_All <-
        Data %>% group_by(Location) %>% select(2:ncol(.)) %>%
        summarise_each(funs(sum(length(which(
            !is.na(.)
        )))))
    CNA_Burden_All <-
        (CNA_Burden_All %>% select(2:ncol(.)) / Number_Genes_All %>%
             select(2:ncol(.))) * 100
    CNA_Burden_All <-
        CNA_Burden_All %>% mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    CNA_Burden_CCA <- Data %>% group_by(Location) %>%
        summarise_each(funs((
            sum(. != 0, na.rm = FALSE) / length(Location)
        ) * 100))
    
    ## CNA Amplification Burden
    CNA_Burden_Amp_All <- Data %>% group_by(Location) %>%
        summarise_each(funs(sum(. > 0, na.rm = TRUE)))
    CNA_Burden_Amp_All <-
        (CNA_Burden_Amp_All[, 2:ncol(CNA_Burden_Amp_All)] /
             Number_Genes_All[, 2:ncol(CNA_Burden_Amp_All)]) *
        100
    CNA_Burden_Amp_All <- CNA_Burden_Amp_All %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    CNA_Burden_Amp_CCA <- Data %>% group_by(Location) %>%
        summarise_each(funs((
            sum(. > 0, na.rm = FALSE) / length(Location)
        ) * 100))
    
    ## CNA Deletion Burden
    CNA_Burden_Del_All <- Data %>% group_by(Location) %>%
        summarise_each(funs(sum(. < 0, na.rm = TRUE)))
    CNA_Burden_Del_All <-
        (CNA_Burden_Del_All[, 2:ncol(CNA_Burden_Del_All)] /
             Number_Genes_All[, 2:ncol(CNA_Burden_Del_All)]) *
        100
    CNA_Burden_Del_All <- CNA_Burden_Del_All %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    CNA_Burden_Del_CCA <- Data %>% group_by(Location) %>%
        summarise_each(funs((
            sum(. < 0, na.rm = FALSE) / length(Location)
        ) * 100))
    
    ## Difference Burden
    CNA_Burden_Diff_All <-
        CNA_Burden_Amp_All[, 2:ncol(CNA_Burden_Amp_All)] -
        CNA_Burden_Del_All[, 2:ncol(CNA_Burden_Del_All)]
    CNA_Burden_Diff_All <- CNA_Burden_Diff_All %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    CNA_Burden_Diff_CCA <-
        CNA_Burden_Amp_CCA[, 2:ncol(CNA_Burden_Amp_CCA)] -
        CNA_Burden_Del_CCA[, 2:ncol(CNA_Burden_Del_CCA)]
    CNA_Burden_Diff_CCA <- CNA_Burden_Diff_CCA %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    ## Percentage of Alterations Amplification Burden
    CNA_Burden_Per_Amp_All  <-
        (abs(CNA_Burden_Amp_All[, 2:ncol(CNA_Burden_Amp_All)]) /
             CNA_Burden_All[, 2:ncol(CNA_Burden_All)]) *
        100
    CNA_Burden_Per_Amp_All[CNA_Burden_Per_Amp_All == "NaN"] <- 0
    CNA_Burden_Per_Amp_All <- CNA_Burden_Per_Amp_All %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    CNA_Burden_Per_Amp_CCA  <-
        (abs(CNA_Burden_Amp_CCA[, 2:ncol(CNA_Burden_Amp_CCA)]) /
             CNA_Burden_CCA[, 2:ncol(CNA_Burden_CCA)]) *
        100
    CNA_Burden_Per_Amp_CCA[CNA_Burden_Per_Amp_CCA == "NaN"] <- 0
    CNA_Burden_Per_Amp_CCA <- CNA_Burden_Per_Amp_CCA %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    ## Percentage of Alterations Deletion Burden
    CNA_Burden_Per_Del_All  <-
        (abs(CNA_Burden_Del_All[, 2:ncol(CNA_Burden_Del_All)]) /
             CNA_Burden_All[, 2:ncol(CNA_Burden_All)]) *
        100
    CNA_Burden_Per_Del_All[CNA_Burden_Per_Del_All == "NaN"] <- 0
    CNA_Burden_Per_Del_All <- CNA_Burden_Per_Del_All %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    CNA_Burden_Per_Del_CCA  <-
        (abs(CNA_Burden_Del_CCA[, 2:ncol(CNA_Burden_Del_CCA)]) /
             CNA_Burden_CCA[, 2:ncol(CNA_Burden_CCA)]) *
        100
    CNA_Burden_Per_Del_CCA[CNA_Burden_Per_Del_CCA == "NaN"] <- 0
    CNA_Burden_Per_Del_CCA <- CNA_Burden_Per_Del_CCA %>%
        mutate(Location = Number_Genes_All$Location) %>%
        select(Location, everything())
    
    list_Burden <-
        list(
            CNA_Burden_All,
            CNA_Burden_Amp_All,
            CNA_Burden_Del_All,
            CNA_Burden_Diff_All,
            CNA_Burden_Per_Amp_All,
            CNA_Burden_Per_Del_All,
            CNA_Burden_CCA,
            CNA_Burden_Amp_CCA,
            CNA_Burden_Del_CCA,
            CNA_Burden_Diff_CCA,
            CNA_Burden_Per_Amp_CCA,
            CNA_Burden_Per_Del_CCA
        )
    
    names_Burden <-
        c(
            "CNA_Burden_All",
            "CNA_Amp_All",
            "CNA_Del_All",
            "CNA_Difference_All",
            "CNA_Per_Amp_All",
            "CNA_Per_Del_All",
            "CNA_Burden_CCA",
            "CNA_Amp_CCA",
            "CNA_Del_CCA",
            "CNA_Difference_CCA",
            "CNA_Per_Amp_CCA",
            "CNA_Per_Del_CCA"
        )
    
    for (i in 1:length(list_Burden)) {
        list_Burden[[i]] <-
            as.data.frame(t(list_Burden[[i]]), stringsAsFactors = FALSE)
        names(list_Burden[[i]]) <- list_Burden[[i]][1, ]
        list_Burden[[i]] <- list_Burden[[i]][-1, ]
        list_Burden[[i]] <-
            lapply(list_Burden[[i]], type.convert, as.is = TRUE)
        list_Burden[[i]] <-
            cbind.data.frame(PATIENT_ID = colnames(Data[, -c(1)]),
                             list_Burden[[i]])
        list_Burden[[i]] <-
            list_Burden[[i]] %>% select(c(variable_order))
        
        for (j in 2:ncol(list_Burden[[i]])) {
            colnames(list_Burden[[i]])[j] <-
                as.character(paste0(names_Burden[i], "_", colnames(list_Burden[[i]])[j]))
        }
    }
    
    names(list_Burden) <- names_Burden
    return(list_Burden)
}

## Run function
CNA_Arm_Burden_Metrics_List <- Per_Arm_CNA_Burden(CNA_Loc)

## Merge all CNA Burden metrics
CNA_Arm_Burden_Metrics <-
    merge_recurse(CNA_Arm_Burden_Metrics_List, by = "PATIENT_ID")

## Run function
Order_Function(CNA_Arm_Burden_Metrics_List)

## Apply new order
CNA_Arm_Burden_Metrics  <- CNA_Arm_Burden_Metrics %>%
    select(c("PATIENT_ID", all_of(new_order)))

# Save Per Arm CNA Score and Burden Metrics
write.table(
    CNA_Arm_Score_Metrics,
    file = "../../data/Processed_Data/Chr_Arm_CNA_Score.txt",
    sep = "\t",
    row.names = FALSE
)
write.table(
    CNA_Arm_Burden_Metrics,
    file = "../../data/Processed_Data/Chr_Arm_CNA_Burden.txt",
    sep = "\t",
    row.names = FALSE
)