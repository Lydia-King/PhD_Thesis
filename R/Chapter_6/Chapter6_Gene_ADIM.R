### Load up Library 
library(tidyverse)
library(dplyr)
library(MCMCglmm)

### 2) Load up data
### Load up and Format Data (ASCAT outputs)
ASCAT_Data_1 <- read.delim("../../data/ASCAT_Data/3Step/ASCAT_Part_1.segments.txt", sep="\t")

ASCAT_Data_3Step_NoNeut<- read.delim("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg19.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

CNA_Loc <- CNA_Loc %>% select(1:7)

### 4) Segment function
## function 1 

Segment_function <- function(data, CNA_Loc){
    seq_list <- list()
    
    for(i in 1:nrow(CNA_Loc)){

        data_in <- data %>% filter(as.numeric(Changepoint) >= CNA_Loc[i,]$Start & as.numeric(Changepoint) <= CNA_Loc[i,]$End) 
        data_out <- data %>% filter(as.numeric(Changepoint) < CNA_Loc[i,]$Start | as.numeric(Changepoint) > CNA_Loc[i,]$End | Category == "NoChangepoint")
        
        data_in <- data_in %>% mutate(TS_Seg = TS, TE_Seg = TE, Category_Seg = Category, Changepoint_Seg = Changepoint) 
        data_out <- data_out %>% mutate(TS_Seg = 0, TE_Seg = 0, Category_Seg = "NoChangepoint", Changepoint_Seg = NA)
        
        data_out <- data_out[!duplicated(data_out[c("Sample", "Chr", "Allele")]),]
        
        int_com <- intersect(data_in[,c("Sample", "Allele")], data_out[,c("Sample", "Allele")])
        
        if(!is_empty(int_com) & nrow(int_com) > 0){
            for(k in 1:nrow(int_com)){
                data_out <- data_out[-c(which(data_out$Sample == int_com$Sample[k] & data_out$Allele == int_com$Allele[k])), ]
            }
        }
        
        tot <- rbind.data.frame(data_in, data_out)     
        
        seq_list[[paste("chr_", unique(data$Chr), "_Gene_", CNA_Loc[i, "Hugo_Symbol"], sep="")]] <- tot
    }
    
    return(seq_list)
}

## function
### 5) Segment Data (Into Different Size Chunks) - Neut 

WG_List_20million_func1 <- list()

for(i in c(1:22, "X")){
    Chr_K <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == i)
    CNA_Loc_Chr <- CNA_Loc %>% mutate(Chromosome = ifelse(Chromosome == 23, "X", Chromosome)) %>% filter(Chromosome == i) %>% arrange(Start)
    WG_List_20million_func1[[paste("Chr", i)]] <- Segment_function(Chr_K, CNA_Loc_Chr)
}



## MCMC multi
Table_function_Pred_Allele_Multi <- function(Dataset){
    
    Res_list <- list()
    
    DF <- as.data.frame(Dataset)
    DF <- within(DF, Category_Seg <- relevel(factor(Category_Seg), ref = "NoChangepoint"))
    DF <- within(DF, Allele <- relevel(factor(Allele), ref = "Major"))
    DF <- DF %>% mutate(TS_Seg = as.numeric(TS_Seg), TE_Seg = as.numeric(TE_Seg)) 
    
    tryCatch(
        {
            mod_TS <- MCMCglmm(cbind(TS_Seg, TE_Seg) ~ trait:(Category_Seg*Allele),  data = DF, family = c("gaussian", "gaussian"), 
                               nitt=60000, thin=10, burnin=15000, rcov = ~ us(trait):units)
            
            SS <- DF %>% mutate(Category_Seg = as.factor(Category_Seg), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category_Seg, .drop = F)
            SS1 <- SS %>% filter(n != 0) %>% arrange(desc(Category_Seg), Allele)
            
            sum_TS <- summary(mod_TS)
            
            
            newdata = data.frame("Category_Seg" = c(unique(DF[,c('Category_Seg', "Allele")])[1]),
                                 "TS_Seg" = c(rep(0, nrow(unique(DF[,c('Category_Seg', "Allele")])[1]))),
                                 "TE_Seg" =  c(rep(0, nrow(unique(DF[,c('Category_Seg', "Allele")])[1]))),
                                 "Allele" = c(unique(DF[,c('Category_Seg', "Allele")])[2]), 
                                 "Sample" = c(paste("Sample", 1:nrow(unique(DF[,c('Category_Seg', "Allele")])[1]))))
            
            pred_TS <- MCMCglmm::predict.MCMCglmm(mod_TS, newdata = newdata, interval = "confidence")
            
            
            Model_DF <- rbind.data.frame(data.frame("Coefficients" = rownames(sum_TS$solutions), "Beta" = sum_TS$solutions[,1], "P" = sum_TS$solutions[,5]))
            
            newdata$id  <- 1:nrow(newdata)
            out <- merge(newdata, SS1, by = c("Category_Seg", "Allele"))
            out <- out[order(out$id), ]
            
            Pred_DF <- data.frame("Category" = rep(unlist(unique(DF[,c('Category_Seg', "Allele")])[1]), 2),
                                  "n" = c(rep(out$n, 2)), 
                                  "Direction" = rep(c("TS", "TE"), each = nrow(unique(DF[,c('Category_Seg', "Allele")])[1])),
                                  "Allele" = rep(c(unlist(unique(DF[,c('Category_Seg', "Allele")])[2])), 2),
                                  "Fit" = c(pred_TS[,1]), 
                                  "LB" = c(pred_TS[,2]), 
                                  "UB" = c(pred_TS[,3]))
            Pred_DF
        },
        error=function(e) {
            Pred_DF1 <- data.frame("Category" = NA,
                                   "n" = NA, 
                                   "Direction" = NA,
                                   "Allele" = NA,
                                   "Fit" = NA, 
                                   "LB" = NA, 
                                   "UB" = NA)
            
            Pred_DF1
            
            
        }
    )
}

## Table for Plot (Multi)
Test_Results_Multi <- function(dataset){
    List1 <- list()
    
    for(chr in 1:length(dataset)){
        List2 <- list()
        
        Dataset <- dataset[[chr]]
        
        print(length(Dataset))
        for(seg in 1:length(Dataset)){
            Dataset1 <- Dataset[[seg]]
            print(paste("Current chr:", chr, sep = " "))
            print(paste("Current Gene:", names(Dataset)[seg], sep = " "))
            List2[[names(Dataset)[seg]]] <- Table_function_Pred_Allele_Multi(Dataset1)
        }
        
        List1[[paste("Chr_", chr, sep="")]] <- List2
    }
    
    return(List1)
}

Test_Results_20million_Multi <- Test_Results_Multi(WG_List_20million_func1)

df3 <- purrr::map(Test_Results_20million_Multi, .f = ~bind_rows(., .id = "Gene"))
df4 <- bind_rows(df3, .id = "Chr")

df4$Gene <- gsub("chr_._Gene_", "", df4$Gene)

write.table(df4, file='../../data/Chapter_6/Model4_Genes_func1_MCMC.txt', quote=FALSE, sep='\t', row.names = F)

### 8) Run Tests - Model 4
# Functions

## LM Uni
Univariate_lm_Model_2 <- function(Dataset, NoChangepoint_Rem = F){
    
    tryCatch(
        {
            
            if(isTRUE(NoChangepoint_Rem)){
                DF <- Dataset #%>% mutate(Category = Type) %>% dplyr::rename(TS = TS, TE = TE)
                DF <- within(DF, Category_Seg <- relevel(factor(Category_Seg), ref = "NoChangepoint"))
                DF <- DF %>% filter(Category_Seg != "NoChangepoint")
                mod_TS <- lm(cbind(TS_Seg) ~ 0 + Category_Seg*Allele, data =  DF)
                mod_TE <- lm(cbind(TE_Seg) ~ 0 + Category_Seg*Allele, data =  DF)
            } else {
                DF <- Dataset# %>% mutate(Category = Type) %>% dplyr::rename(TS = TS, TE = TE)
                DF <- within(DF, Category_Seg <- relevel(factor(Category_Seg), ref = "NoChangepoint"))
                mod_TS <- lm(cbind(TS_Seg) ~ Category_Seg*Allele, data =  DF)
                mod_TE <- lm(cbind(TE_Seg) ~ Category_Seg*Allele, data =  DF)
            }
            
            SS <- DF %>% mutate(Category_Seg = as.factor(Category_Seg), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category_Seg, .drop = F)
            
            if(isTRUE(NoChangepoint_Rem)){
                SS <- SS %>% filter(Category_Seg != "NoChangepoint")
            } 
            
            
            sum_TS <- summary(mod_TS)
            sum_TE <- summary(mod_TE)
            
            pval_TS <- c()
            j <- 1
            for(i in 1:length(names(mod_TS$coefficients))){
                if(rownames(sum_TS$coefficients)[i] %in% c(names(mod_TS$coefficients))){
                    pval_TS <- c(pval_TS, sum_TS$coefficients[j,4])
                    j <- j + 1
                } else {
                    pval_TS <- c(pval_TS, NA)
                }
            }
            
            pval_TE <- c()
            j <- 1
            for(i in 1:length(names(mod_TE$coefficients))){
                if(rownames(sum_TE$coefficients)[i] %in% c(names(mod_TE$coefficients))){
                    pval_TE <- c(pval_TE, sum_TE$coefficients[j,4])
                    j <- j + 1
                } else {
                    pval_TE <- c(pval_TE, NA)
                }
            }
            
            cate <- names(mod_TS$coefficients)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub(":AlleleMinor", "", cate)
            cate <- gsub("AlleleMinor", names(mod_TS$coefficients)[1], cate)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub("Category_Seg", "", cate)
            pred_TS <- predict(mod_TS, newdata = data.frame("Category_Seg" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
            pred_TE <- predict(mod_TE, newdata = data.frame("Category_Seg" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
            
            Model_DF1 <- rbind.data.frame(data.frame("Category" = cate, "Chr" = rep(1, length(mod_TS$coefficients)),
                                                     "Direction" = rep("TS", length(mod_TS$coefficients)), "Samplesize" = SS$n,
                                                     "Beta" = unname(mod_TS$coefficients), "P" = pval_TS,
                                                     "LB" = pred_TS[,2], "UB" = pred_TS[,3], "Pred" = pred_TS[,1],
                                                     "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),
                                          data.frame("Category" = cate, "Chr" = rep(1, length(mod_TE$coefficients)), "Direction" = rep("TE", length(mod_TE$coefficients)),
                                                     "Samplesize" = SS$n,  "Beta" = unname(mod_TE$coefficients), "P" = pval_TE,
                                                     "LB" = pred_TE[,2], "UB" = pred_TE[,3], "Pred" = pred_TE[,1],
                                                     "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))))
            
            Data_Chr1 <- Model_DF1  %>% mutate(Category = gsub("Category_Seg", "", Category)) %>% mutate(Category = gsub("Category", "", Category))
            Data_Chr1 <- Data_Chr1 %>% filter(Samplesize > 0)
            Data_Chr1$Category <- gsub("\\(Intercept\\)", "NoChangepoint", Data_Chr1$Category)
            Data_Chr1 <- Data_Chr1 %>% mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Amp/Neut", "Del/Neut", "Amp/Del", "Del/Amp")))
            
            return(Data_Chr1)
            
        },
        error=function(e) {
            Pred_DF1 <- data.frame("Category" = NA,
                                   "Chr" = 1,
                                   "Direction" = NA,
                                   "Samplesize" = NA,
                                   "Beta" = NA, 
                                   "P" = NA,
                                   "LB" = NA, 
                                   "UB" = NA, 
                                   "Pred" = NA,
                                   "Allele" = NA)
            
            Pred_DF1
            
        }
    )
    
} 

Test_Results_Multi <- function(dataset){
    List1 <- list()
    
    for(chr in 1:length(dataset)){
        List2 <- list()
        
        Dataset <- dataset[[chr]]
    
        for(seg in 1:length(Dataset)){
            Dataset1 <- Dataset[[seg]]
            print(paste("Current chr:", chr, sep = " "))
            print(paste("Current Gene:", names(Dataset)[seg], sep = " "))
            List2[[names(Dataset)[seg]]] <- Univariate_lm_Model_2(Dataset1)
        }
        
        List1[[paste("Chr_", chr, sep="")]] <- List2
    }
    
    return(List1)
}

Test_Results_5million_Multi <- Test_Results_Multi(WG_List_20million_func1)


df1 <- purrr::map(Test_Results_5million_Multi, .f = ~bind_rows(., .id = "Gene"))
df2 <- bind_rows(df1, .id = "Chr")

write.table(df2, file='../../data/Chapter_6/Model4_Genes_func1_LM.txt', quote=FALSE, sep='\t', row.names = F)
