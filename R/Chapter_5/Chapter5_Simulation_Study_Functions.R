# Chapter 5: Simulation Study Functions

## Allele Dependant Models 
# Univariate linear model
Univariate_lm_Model_2 <- function(Dataset, NoChangepoint_Rem = F){
    
    if(isTRUE(NoChangepoint_Rem)){
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        DF <- DF %>% filter(Category != "NoChangepoint")
        mod_TS <- lm(cbind(Towards.Start) ~ 0 + Category*Allele, data =  DF)
        mod_TE <- lm(cbind(Towards.End) ~ 0 + Category*Allele, data =  DF)
    } else {
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        mod_TS <- lm(cbind(Towards.Start) ~ Category*Allele, data =  DF)
        mod_TE <- lm(cbind(Towards.End) ~ Category*Allele, data =  DF)
    }
    
    SS <- DF %>% mutate(Category = as.factor(Category), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category, .drop = F)
    
    if(isTRUE(NoChangepoint_Rem)){
        SS <- SS %>% filter(Category != "NoChangepoint")
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
    cate <- gsub("Category", "", cate)
    pred_TS <- predict(mod_TS, newdata = data.frame("Category" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
    pred_TE <- predict(mod_TE, newdata = data.frame("Category" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
    
    Model_DF1 <- rbind.data.frame(data.frame("Category" = cate, "Chr" = rep(1, length(mod_TS$coefficients)),
                                             "Direction" = rep("TS", length(mod_TS$coefficients)), "Samplesize" = SS$n,
                                             "Beta" = unname(mod_TS$coefficients), "P" = pval_TS,
                                             "LB" = pred_TS[,2], "UB" = pred_TS[,3], "Pred" = pred_TS[,1],
                                             "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),
                                  data.frame("Category" = cate, "Chr" = rep(1, length(mod_TE$coefficients)), "Direction" = rep("TE", length(mod_TE$coefficients)),
                                             "Samplesize" = SS$n,  "Beta" = unname(mod_TE$coefficients), "P" = pval_TE,
                                             "LB" = pred_TE[,2], "UB" = pred_TE[,3], "Pred" = pred_TE[,1],
                                             "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))))
    
    Data_Chr1 <- Model_DF1  %>% mutate(Category = gsub("Category", "", Category)) %>% mutate(Category = gsub("Category", "", Category))
    Data_Chr1 <- Data_Chr1 %>% filter(Samplesize > 0)
    Data_Chr1$Category <- gsub("\\(Intercept\\)", "NoChangepoint", Data_Chr1$Category)
    Data_Chr1 <- Data_Chr1 %>% mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Amp/Neut", "Del/Neut", "Amp/Del", "Del/Amp")))
    
    return(Data_Chr1)
    
} 

# Univariate MCMC 
Univariate_MCMC_Model_2 <- function(Dataset, family = c("gaussian"), NoChangepoint_Rem = F){
    
    if(isTRUE(NoChangepoint_Rem)){
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        DF <- DF %>% filter(Category != "NoChangepoint")
        mod_TS <- MCMCglmm(cbind(Towards.Start) ~ 0 + Category*Allele, data = DF, family = family, nitt=60000, thin=10, burnin=15000)
        mod_TE <- MCMCglmm(cbind(Towards.End) ~ 0 + Category*Allele,  data = DF, family = family, nitt=60000, thin=10, burnin=15000)
    } else {
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        mod_TS <- MCMCglmm(cbind(Towards.Start) ~ Category*Allele, data = DF, family = family, nitt=60000, thin=10, burnin=15000)
        mod_TE <- MCMCglmm(cbind(Towards.End) ~ Category*Allele,  data = DF, family = family, nitt=60000, thin=10, burnin=15000)
    }
    
    
    SS <- DF %>% mutate(Category = as.factor(Category), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category, .drop = F)
    SS1 <- SS %>% filter(n != 0)
    
    sum_TS <- summary(mod_TS)
    sum_TE <- summary(mod_TE)
    
    pred_TS <- MCMCglmm::predict.MCMCglmm(mod_TS, newdata = data.frame("Category" = SS1$Category, "Allele" = SS1$Allele,
                                                                       "Towards.Start" = c(rep(0, (length(SS1$Allele)))),
                                                                       "Sample" = c(paste("Sample", 1:((length(SS1$Allele))))))
                                          , interval = "confidence")
    
    pred_TE <- MCMCglmm::predict.MCMCglmm(mod_TE, newdata = data.frame("Category" = SS1$Category, "Allele" = SS1$Allele,
                                                                       "Towards.End" =  c(rep(0, (length(SS1$Allele)))),
                                                                       "Sample" = c(paste("Sample", 1:((length(SS1$Allele))))))
                                          , interval = "confidence")
    
    Model_DF <- rbind.data.frame(data.frame("Category" = SS1$Category, "Chr" = rep(1, length(SS1$Category)),
                                            "Direction" = rep("TS", length(SS1$Category)), "Samplesize" = SS1$n,  "Beta" = sum_TS$solutions[,1], "Pred" = pred_TS[,1],
                                            "LB" = pred_TS[,2], "UB" = pred_TS[,3], "P" = sum_TS$solutions[,5], "Allele" = SS1$Allele), 
                                 data.frame("Category" = SS1$Category, "Chr" = rep(1, length(SS1$Category)),
                                            "Direction" = rep("TE", length(SS1$Category)), "Samplesize" = SS1$n,  "Beta" = sum_TE$solutions[,1], "Pred" = pred_TE[,1],
                                            "LB" = pred_TE[,2], "UB" = pred_TE[,3], "P" = sum_TE$solutions[,5], "Allele" = SS1$Allele))
    
    Data_Chr <- Model_DF  %>% mutate(Category = gsub("Category", "", Category)) %>% mutate(Category = gsub("Category", "", Category))
    Data_Chr$Category <- gsub(":AlleleMinor", "", Data_Chr$Category)
    Data_Chr$Category <- gsub("AlleleMinor", "NoChangepoint", Data_Chr$Category)
    
    Data_Chr <- Data_Chr%>% mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Amp/Neut", "Del/Neut", "Amp/Del", "Del/Amp")))
    rownames(Data_Chr) <- NULL
    return(Data_Chr)
}

# Multivariate linear model
Multivariate_lm_Model_2 <-  function(Dataset, NoChangepoint_Rem = F){
    
    if(isTRUE(NoChangepoint_Rem)){
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        DF <- DF %>% filter(Category != "NoChangepoint")
        mod_TS <- lm(cbind(Towards.Start, Towards.End) ~ 0 + Category*Allele, data =  DF)
    } else {
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        mod_TS <- lm(cbind(Towards.Start, Towards.End) ~ Category*Allele, data =  DF)
    }
    
    SS <- DF %>% mutate(Category = as.factor(as.character(Category)), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category, .drop = F)
    
    sum_TS <- summary(mod_TS)
    
    pval_TS <- c()
    j <- 1
    for(i in 1:length(rownames(mod_TS$coefficients))){
        if(c(rownames(mod_TS$coefficients))[i] %in% c(rownames(sum_TS$`Response Towards.Start`$coefficients))){
            pval_TS <- c(pval_TS, sum_TS$`Response Towards.Start`$coefficients[j, 4])
            j <- j + 1
        } else {
            pval_TS <- c(pval_TS, NA)
        }
    }
    
    pval_TE <- c()
    j <- 1
    for(i in 1:length(rownames(mod_TS$coefficients))){
        if(c(rownames(mod_TS$coefficients))[i] %in% c(rownames(sum_TS$`Response Towards.End`$coefficients))){
            pval_TE <- c(pval_TE, sum_TS$`Response Towards.End`$coefficients[j, 4])
            j <- j + 1
        } else {
            pval_TE <- c(pval_TE, NA)
        }
    }
    
    cate <- rownames(mod_TS$coefficients)
    cate <- gsub(":AlleleMinor", "", cate)
    cate <- gsub("AlleleMinor", "NoChangepoint", cate)
    cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
    cate <- gsub("Category", "", cate)
    pred_TS <- predict(mod_TS, newdata = data.frame("Category" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
    
    SS <- SS %>% arrange(Allele, factor(Category, levels = c(unique(cate))))
    
    Model_DF1 <- data.frame("Category" = rep(cate, 2), "Chr" = rep(1, length(mod_TS$coefficients)),
                            "Direction" = c(rep("TS", length(mod_TS$coefficients)/2), rep("TE", length(mod_TS$coefficients)/2)), "Samplesize" = rep(SS$n, 2),
                            "Beta" = c(unname(mod_TS$coefficients)[,1], unname(mod_TS$coefficients)[,2]), "P" = c(pval_TS, pval_TE),
                            "LB" = rep(NA, length(mod_TS$coefficients)), "UB" = rep(NA, length(mod_TS$coefficients)), "Pred" = c(pred_TS[,1], pred_TS[,2]),
                            "Allele" = c(rep(c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor")))), 2)))
    
    Data_Chr1 <- Model_DF1  %>% mutate(Category = gsub("Category", "", Category)) %>% mutate(Category = gsub("Category", "", Category))
    Data_Chr1$Category <- gsub("\\(Intercept\\)", "NoChangepoint", Data_Chr1$Category)
    Data_Chr1 <- Data_Chr1 %>% filter(Samplesize > 0)
    Data_Chr1 <- Data_Chr1 %>% mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Amp/Neut", "Del/Neut", "Amp/Del", "Del/Amp")))
    
    return(Data_Chr1)
}

# Multivariate MCMC model
Multivariate_MCMC_Model_2 <-  function(Dataset, NoChangepoint_Rem = F){
    
    if(isTRUE(NoChangepoint_Rem)){
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        DF <- within(DF, Allele <- relevel(factor(Allele), ref = "Major"))
        DF <- DF %>% filter(Category != "NoChangepoint")
        mod_TS <- MCMCglmm(cbind(Towards.Start, Towards.End) ~ 0 + trait:(Category*Allele),  data = DF, family = c("gaussian", "gaussian"), nitt=60000, thin=10, burnin=15000, rcov = ~ us(trait):units)
    } else {
        DF <- Dataset %>% mutate(Category = Type) %>% dplyr::rename(Towards.Start = TS, Towards.End = TE)
        DF <- within(DF, Category <- relevel(factor(Category), ref = "NoChangepoint"))
        DF <- within(DF, Allele <- relevel(factor(Allele), ref = "Major"))
        mod_TS <- MCMCglmm(cbind(Towards.Start, Towards.End) ~ trait:(Category*Allele),  data = DF, family = c("gaussian", "gaussian"), nitt=60000, thin=10, burnin=15000, rcov = ~ us(trait):units)
    }
    
    SS <- DF %>% mutate(Category = as.factor(Category), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category, .drop = F)
    SS1 <- SS %>% filter(n != 0) %>% arrange(Category, Allele)
    
    sum_TS <- summary(mod_TS)
    
    newdata = data.frame("Category" = c(unique(DF[,c('Category', "Allele")])[1]),
                         "Towards.Start" = c(rep(0, nrow(unique(DF[,c('Category', "Allele")])[1]))),
                         "Towards.End" =  c(rep(0, nrow(unique(DF[,c('Category', "Allele")])[1]))),
                         "Allele" = c(unique(DF[,c('Category', "Allele")])[2]), 
                         "Sample" = c(paste("Sample", 1:nrow(unique(DF[,c('Category', "Allele")])[1]))))
    
    pred_TS <- MCMCglmm::predict.MCMCglmm(mod_TS, newdata = newdata, interval = "confidence")
    
    Model_DF <- rbind.data.frame(data.frame("Coefficients" = rownames(sum_TS$solutions), "Beta" = sum_TS$solutions[,1], "P" = sum_TS$solutions[,5]))
    
    Pred_DF <- data.frame("Category" = rep(unlist(unique(DF[,c('Category', "Allele")])[1]), 2),
                          "n" = c(rep(SS1$n, 2)), 
                          "Direction" = rep(c("TS", "TE"), each = nrow(unique(DF[,c('Category', "Allele")])[1])),
                          "Allele" = rep(c(unlist(unique(DF[,c('Category', "Allele")])[2])), 2),
                          "Pred" = c(pred_TS[,1]), 
                          "LB" = c(pred_TS[,2]), 
                          "UB" = c(pred_TS[,3]))
    
    
    
    Data_Chr <- Model_DF
    rownames(Data_Chr) <- NULL
    
    return(Pred_DF)
    
}

# Test Results (Dot Plot of selected Data)
Test_lm_MCMC_Results_Uni_DP <- function(dataset, tests = c("lm", "MCMC"), NoChangepoint_Rem = F){
    List1 <- list()
    
    for(k in 1:length(dataset)){
        Dataset_Test <- dataset[k]
        if(tests == "lm"){
            List1[[paste("Dataset_", k, sep="")]] <- Univariate_lm_Model_2(Dataset_Test[[1]], 
                                                                           NoChangepoint_Rem = NoChangepoint_Rem)
        } else {
            List1[[paste("Dataset_", k, sep="")]] <- Univariate_MCMC_Model_2(Dataset_Test[[1]], 
                                                                             NoChangepoint_Rem = NoChangepoint_Rem)
        }
    }
    
    return(List1)
}
Test_lm_MCMC_Results_Multi_DP <- function(dataset, tests = c("lm", "MCMC"), NoChangepoint_Rem = F){
    List1 <- list()
    
    for(k in 1:length(dataset)){
        Dataset_Test <- dataset[k]
        if(tests == "lm"){
            List1[[paste("Dataset_", k, sep="")]] <- Multivariate_lm_Model_2(Dataset_Test[[1]], 
                                                                             NoChangepoint_Rem = NoChangepoint_Rem)
        } else {
            List1[[paste("Dataset_", k, sep="")]] <- Multivariate_MCMC_Model_2(Dataset_Test[[1]], 
                                                                               NoChangepoint_Rem = NoChangepoint_Rem)
        }
    }
    
    return(List1)
}

# Test Results for Simulation Study
Test_lm_MCMC_Results_Uni <- function(dataset, tests = c("lm", "MCMC"), NoChangepoint_Rem = F){
    List1 <- list()
    List2 <- list()
    List3 <- list()
    
    for(i in 1:length(dataset)){
        Dataset <- dataset[[i]]
        
        for(j in 1:length(Dataset)){
            Dataset1 <- Dataset[[j]]
            
            for(k in 1:length(Dataset1)){
                Dataset_Test <- Dataset1[k]
                if(tests == "lm"){
                    List3[[paste("Dataset_", k, sep="")]] <- Univariate_lm_Model_2(Dataset_Test[[1]], 
                                                                                   NoChangepoint_Rem = NoChangepoint_Rem)
                } else {
                    List3[[paste("Dataset_", k, sep="")]] <- Univariate_MCMC_Model_2(Dataset_Test[[1]], 
                                                                                     NoChangepoint_Rem = NoChangepoint_Rem)
                }
            }
            
            List2[[paste("Samplesize_", j, sep="")]] <- List3
        }
        
        List1[[paste("Percentage_", i, sep="")]] <- List2
    }
    
    return(List1)
}
Test_lm_MCMC_Results_Multi <- function(dataset, tests = c("lm", "MCMC"), NoChangepoint_Rem = F){
    List1 <- list()
    List2 <- list()
    List3 <- list()
    
    for(i in 1:length(dataset)){
        Dataset <- dataset[[i]]
        
        for(j in 1:length(Dataset)){
            Dataset1 <- Dataset[[j]]
            
            for(k in 1:length(Dataset1)){
                Dataset_Test <- Dataset1[k]
                if(tests == "lm"){
                    List3[[paste("Dataset_", k, sep="")]] <- Multivariate_lm_Model_2(Dataset_Test[[1]], 
                                                                                     NoChangepoint_Rem = NoChangepoint_Rem)
                } else {
                    List3[[paste("Dataset_", k, sep="")]] <- Multivariate_MCMC_Model_2(Dataset_Test[[1]], 
                                                                                       NoChangepoint_Rem = NoChangepoint_Rem)
                }
            }
            
            List2[[paste("Samplesize_", j, sep="")]] <- List3
        }
        
        List1[[paste("Percentage_", i, sep="")]] <- List2
    }
    
    return(List1)
}

# Summary Table for Plots
Summarise_Test <- function(dataset, Multi_lm = F){
    Sample_Size_List <- list()
    Trial_List <- list()
    
    for(i in 1:length(dataset)){
        Dataset <- dataset[[i]]
        
        for(j in 1:length(Dataset)){
            Dataset1 <- Dataset[[j]]
            results_1 <- list()
            
            for(k in 1:length(Dataset1)){
                if(isTRUE(Multi_lm)){
                    results_1[[k]] <- Dataset1[[k]] %>% mutate("Res" = ifelse(.$Pred > 0, "Sig", "NSig"))
                } else {
                    results_1[[k]] <- Dataset1[[k]] %>% mutate("Res" = ifelse(.$LB > 0, "Sig", "NSig"))
                }
            }
            
            results_1 <- bind_rows(results_1, .id = "group_id")
            results_2 <- results_1 %>% group_by(Allele, Category, Direction, Res) %>% summarise(n = n(), Prop = ifelse(Res == "NSig", 1 - n()/20, n()/20))
            Sample_Size_List[[paste("Samplesize_", j, sep="")]] <- unique(results_2)
        }
        
        Trial_List[[paste("Percentage_", i, sep="")]] <- Sample_Size_List
    }
    return(Trial_List)
}

# Unlist function (unlist results)
Sample_Size <- c(20, 50, 80, 100, 200, 500)
Percentage_Neut_Amp <- c(10, 20, 30, 40, 50, 60, 70, 80, 90)
Percentage_B <- c(10, 20, 30, 40, 50, 60, 70)

Unlist_function <- function(dataset, per){
    Samp_List <- list()
    PC_List <- list()
    
    for(i in 1:length(dataset)){
        Dataset <- dataset[[i]]
        
        for(j in 1:length(Dataset)){
            Dataset1 <- Dataset[[j]]
            Dataset1 <- Dataset1 %>% mutate("SampleSize" = Sample_Size[j]) %>% mutate("Percentage" = per[i])
            Samp_List[[j]] <- Dataset1
        }
        
        PC_List[[i]] <- Samp_List
    }
    final <- bind_rows(PC_List, .id = "group_id")
    return(final)
}
