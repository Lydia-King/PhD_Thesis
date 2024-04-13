# Chapter 5: Allele-Dependent (AD) Models

## Load up Libraries
library(gt)
library(tidyverse)
library(MCMCglmm)

## Load Up Data
DF3 <-
    read.delim("../../data/Chapter_5/Example_Data_1_NoNeut.txt",
               sep = "\t")


# Allele-Dependent (AD) Models
## Univariate and Multivariate Models Function
Table_function_Pred <- function(Dataset,
                                model = c("lm", "mcmc"),
                                model_type = c("univariate", "multivariate"),
                                category_spec = c("Intercept", "No_Intercept")) {
    if (category_spec == "No_Intercept") {
        DF <- Dataset
        DF <- DF %>% filter(Category != "NoChangepoint")
    } else {
        DF <- Dataset
        DF <-
            within(DF, Category <-
                       relevel(factor(Category), ref = "NoChangepoint"))
    }
    
    if (model_type == "univariate" & category_spec == "Intercept") {
        if (model == "lm") {
            mod_TS <- lm(cbind(TS) ~ Category * Allele, data =  DF)
            mod_TE <-
                lm(cbind(TE) ~ Category * Allele, data =  DF)
        } else {
            mod_TS <-
                MCMCglmm(
                    cbind(TS) ~ Category * Allele,
                    data = DF,
                    family = c("gaussian"),
                    nitt = 60000,
                    thin = 10,
                    burnin = 15000
                )
            mod_TE <-
                MCMCglmm(
                    cbind(TE) ~ Category * Allele,
                    data = DF,
                    family = c("gaussian"),
                    nitt = 60000,
                    thin = 10,
                    burnin = 15000
                )
        }
    } else if (model_type == "univariate" &
               category_spec == "No_Intercept") {
        if (model == "lm") {
            mod_TS <- lm(cbind(TS) ~ 0 + Category * Allele, data =  DF)
            mod_TE <-
                lm(cbind(TE) ~ 0 + Category * Allele, data =  DF)
        } else {
            mod_TS <-
                MCMCglmm(
                    cbind(TS) ~ 0 + Category * Allele,
                    data = DF,
                    family = c("gaussian"),
                    nitt = 60000,
                    thin = 10,
                    burnin = 15000
                )
            mod_TE <-
                MCMCglmm(
                    cbind(TE) ~ 0 + Category * Allele,
                    data = DF,
                    family = c("gaussian"),
                    nitt = 60000,
                    thin = 10,
                    burnin = 15000
                )
        }
    } else if (model_type == "multivariate" &
               category_spec == "Intercept") {
        if (model == "lm") {
            mod_TS <-
                lm(cbind(TS, TE) ~ Category * Allele,
                   data =  DF)
        } else {
            mod_TS <-
                MCMCglmm(
                    cbind(TS, TE) ~ trait:(Category * Allele),
                    data = DF3,
                    family = c("gaussian", "gaussian"),
                    nitt = 60000,
                    thin = 10,
                    burnin = 15000,
                    rcov = ~ us(trait):units
                )
        }
    } else if (model_type == "multivariate" &
               category_spec == "No_Intercept") {
        if (model == "lm") {
            mod_TS <-
                lm(cbind(TS, TE) ~ 0 + Category * Allele,
                   data =  DF)
        } else {
            mod_TS <-
                MCMCglmm(
                    cbind(TS, TE) ~ 0 + trait:(Category * Allele),
                    data = DF,
                    family = c("gaussian", "gaussian"),
                    nitt = 60000,
                    thin = 10,
                    burnin = 15000,
                    rcov = ~ us(trait):units
                )
        }
    }
    
    SS <-
        DF %>% mutate(Category = as.factor(Category),
                      Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category, .drop = F)
    SS1 <- SS %>% filter(n != 0) %>% arrange(Category, Allele)
    
    if (model_type == "univariate") {
        sum_TS <- summary(mod_TS)
        sum_TE <- summary(mod_TE)
    } else {
        sum_TS <- summary(mod_TS)
    }
    
    if (model_type == "univariate" & category_spec == "Intercept") {
        if (model == "lm") {
            pval_TS <- c()
            j <- 1
            for (i in 1:length(names(mod_TS$coefficients))) {
                if (rownames(sum_TS$coefficients)[i] %in% c(names(mod_TS$coefficients))) {
                    pval_TS <- c(pval_TS, sum_TS$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TS <- c(pval_TS, NA)
                }
            }
            
            pval_TE <- c()
            j <- 1
            for (i in 1:length(names(mod_TE$coefficients))) {
                if (rownames(sum_TE$coefficients)[i] %in% c(names(mod_TE$coefficients))) {
                    pval_TE <- c(pval_TE, sum_TE$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TE <- c(pval_TE, NA)
                }
            }
            
            cate <- names(mod_TS$coefficients)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub(":AlleleMinor", "", cate)
            cate <- gsub("AlleleMinor", "NoChangepoint", cate)
            cate <- gsub("Category", "", cate)
            pred_TS <-
                predict(
                    mod_TS,
                    newdata = data.frame(
                        "Category" = cate,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    interval = "confidence"
                )
            pred_TE <-
                predict(
                    mod_TE,
                    newdata = data.frame(
                        "Category" = cate,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    interval = "confidence"
                )
            
            Model_DF1 <-
                rbind.data.frame(
                    data.frame(
                        "Coefficients" = names(mod_TS$coefficients),
                        "Chr" = rep(1, length(mod_TS$coefficients)),
                        "Direction" = rep("TS", length(mod_TS$coefficients)),
                        "Samplesize" = SS$n,
                        "Beta" = unname(mod_TS$coefficients),
                        "P" = pval_TS,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    data.frame(
                        "Coefficients" = names(mod_TS$coefficients),
                        "Chr" = rep(1, length(mod_TE$coefficients)),
                        "Direction" = rep("TE", length(mod_TE$coefficients)),
                        "Samplesize" = SS$n,
                        "Beta" = unname(mod_TE$coefficients),
                        "P" = pval_TE,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    )
                )
            
            Pred_DF <- data.frame(
                "Category" = rep(c(cate), 2),
                "n" = c(rep(SS$n, 2)),
                "Direction" = c(rep(c(
                    "TS", "TE"
                ), each = length(
                    c(cate)
                ))),
                "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                "UB" = c(pred_TS[, 3], pred_TE[, 3]),
                "Allele" = c(rep("Major", length(
                    which(SS$Allele == "Major")
                )), rep("Minor", length(
                    which(SS$Allele == "Minor")
                )))
            )
            
            
            Data_Chr <- Model_DF1
            rownames(Data_Chr) <- NULL
        } else {
            pred_TS <-
                MCMCglmm::predict.MCMCglmm(
                    mod_TS,
                    newdata = data.frame(
                        "Category" = SS1$Category,
                        "Allele" = SS1$Allele,
                        "TS" = c(rep(0, (
                            length(SS1$Allele)
                        ))),
                        "Sample" = c(paste("Sample", 1:((length(SS1$Allele))
                        )))
                    )
                    ,
                    interval = "confidence"
                )
            
            pred_TE <-
                MCMCglmm::predict.MCMCglmm(
                    mod_TE,
                    newdata = data.frame(
                        "Category" = SS1$Category,
                        "Allele" = SS1$Allele,
                        "TE" =  c(rep(0, (
                            length(SS1$Allele)
                        ))),
                        "Sample" = c(paste("Sample", 1:((length(SS1$Allele))
                        )))
                    )
                    ,
                    interval = "confidence"
                )
            
            Model_DF <-
                rbind.data.frame(
                    data.frame(
                        "Coefficients" = SS1$Category,
                        "Chr" = rep(1, length(SS1$Category)),
                        "Direction" = rep("TS", length(SS1$Category)),
                        "Samplesize" = SS1$n,
                        "Beta" = sum_TS$solutions[, 1],
                        "P" = sum_TS$solutions[, 5],
                        "Allele" = SS1$Allele
                    ),
                    data.frame(
                        "Coefficients" = SS1$Category,
                        "Chr" = rep(1, length(SS1$Category)),
                        "Direction" = rep("TE", length(SS1$Category)),
                        "Samplesize" = SS1$n,
                        "Beta" = sum_TE$solutions[, 1],
                        "P" = sum_TE$solutions[, 5],
                        "Allele" = SS1$Allele
                    )
                )
            
            
            Pred_DF <-
                data.frame(
                    "Category" = c(SS1$Category, SS1$Category),
                    "n" = c(rep(SS1$n, 2)),
                    "Direction" = c(rep("TS", length(
                        SS1$Category
                    )),  rep("TE", length(
                        SS1$Category
                    ))),
                    "Allele" = c(SS1$Allele, SS1$Allele),
                    "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                    "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                    "UB" = c(pred_TS[, 3], pred_TE[, 3])
                )
            
            
            Data_Chr <- Model_DF
            rownames(Data_Chr) <- NULL
        }
    } else if (model_type == "univariate" &
               category_spec == "No_Intercept") {
        if (model == "lm") {
            pval_TS <- c()
            j <- 1
            for (i in 1:length(names(mod_TS$coefficients))) {
                if (rownames(sum_TS$coefficients)[i] %in% c(names(mod_TS$coefficients))) {
                    pval_TS <- c(pval_TS, sum_TS$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TS <- c(pval_TS, NA)
                }
            }
            
            pval_TE <- c()
            j <- 1
            for (i in 1:length(names(mod_TE$coefficients))) {
                if (rownames(sum_TE$coefficients)[i] %in% c(names(mod_TE$coefficients))) {
                    pval_TE <- c(pval_TE, sum_TE$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TE <- c(pval_TE, NA)
                }
            }
            
            cate <- names(mod_TS$coefficients)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub(":AlleleMinor", "", cate)
            cate <- gsub("AlleleMinor", cate[1], cate)
            cate <- gsub("Category", "", cate)
            SS <- SS %>% filter(Category != "NoChangepoint")
            pred_TS <-
                predict(
                    mod_TS,
                    newdata = data.frame(
                        "Category" = cate,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    interval = "confidence"
                )
            pred_TE <-
                predict(
                    mod_TE,
                    newdata = data.frame(
                        "Category" = cate,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    interval = "confidence"
                )
            
            Model_DF1 <-
                rbind.data.frame(
                    data.frame(
                        "Coefficients" = names(mod_TS$coefficients),
                        "Chr" = rep(1, length(mod_TS$coefficients)),
                        "Direction" = rep("TS", length(mod_TS$coefficients)),
                        "Samplesize" = SS$n,
                        "Beta" = unname(mod_TS$coefficients),
                        "P" = pval_TS,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    data.frame(
                        "Coefficients" = names(mod_TS$coefficients),
                        "Chr" = rep(1, length(mod_TE$coefficients)),
                        "Direction" = rep("TE", length(mod_TE$coefficients)),
                        "Samplesize" = SS$n,
                        "Beta" = unname(mod_TE$coefficients),
                        "P" = pval_TE,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    )
                )
            
            Pred_DF <- data.frame(
                "Category" = rep(c(cate), 2),
                "n" = c(rep(SS$n, 2)),
                "Direction" = c(rep(c(
                    "TS", "TE"
                ), each = length(
                    c(cate)
                ))),
                "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                "UB" = c(pred_TS[, 3], pred_TE[, 3]),
                "Allele" = c(rep("Major", length(
                    which(SS$Allele == "Major")
                )), rep("Minor", length(
                    which(SS$Allele == "Minor")
                )))
            )
            
            
            Data_Chr <- Model_DF1
            rownames(Data_Chr) <- NULL
            
        } else {
            pred_TS <-
                MCMCglmm::predict.MCMCglmm(
                    mod_TS,
                    newdata = data.frame(
                        "Category" = SS1$Category,
                        "Allele" = SS1$Allele,
                        "TS" = c(rep(0, (
                            length(SS1$Allele)
                        ))),
                        "Sample" = c(paste("Sample", 1:((length(SS1$Allele))
                        )))
                    )
                    ,
                    interval = "confidence"
                )
            
            pred_TE <-
                MCMCglmm::predict.MCMCglmm(
                    mod_TE,
                    newdata = data.frame(
                        "Category" = SS1$Category,
                        "Allele" = SS1$Allele,
                        "TE" =  c(rep(0, (
                            length(SS1$Allele)
                        ))),
                        "Sample" = c(paste("Sample", 1:((length(SS1$Allele))
                        )))
                    )
                    ,
                    interval = "confidence"
                )
            
            Model_DF <-
                rbind.data.frame(
                    data.frame(
                        "Coefficients" = SS1$Category,
                        "Chr" = rep(1, length(SS1$Category)),
                        "Direction" = rep("TS", length(SS1$Category)),
                        "Samplesize" = SS1$n,
                        "Beta" = sum_TS$solutions[, 1],
                        "P" = sum_TS$solutions[, 5],
                        "Allele" = SS1$Allele
                    ),
                    data.frame(
                        "Coefficients" = SS1$Category,
                        "Chr" = rep(1, length(SS1$Category)),
                        "Direction" = rep("TE", length(SS1$Category)),
                        "Samplesize" = SS1$n,
                        "Beta" = sum_TE$solutions[, 1],
                        "P" = sum_TE$solutions[, 5],
                        "Allele" = SS1$Allele
                    )
                )
            
            
            Pred_DF <-
                data.frame(
                    "Category" = c(SS1$Category, SS1$Category),
                    "n" = c(rep(SS1$n, 2)),
                    "Direction" = c(rep("TS", length(
                        SS1$Category
                    )),  rep("TE", length(
                        SS1$Category
                    ))),
                    "Allele" = c(SS1$Allele, SS1$Allele),
                    "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                    "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                    "UB" = c(pred_TS[, 3], pred_TE[, 3])
                )
            
            
            Data_Chr <- Model_DF
            rownames(Data_Chr) <- NULL
        }
    } else if (model_type == "multivariate" &
               category_spec == "Intercept") {
        if (model == "lm") {
            pval_TS <- c()
            j <- 1
            for (i in 1:length(rownames(mod_TS$coefficients))) {
                if (c(rownames(mod_TS$coefficients))[i] %in% c(rownames(
                    sum_TS$`Response TS`$coefficients
                ))) {
                    pval_TS <-
                        c(pval_TS,
                          sum_TS$`Response TS`$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TS <- c(pval_TS, NA)
                }
            }
            
            pval_TE <- c()
            j <- 1
            for (i in 1:length(rownames(mod_TS$coefficients))) {
                if (c(rownames(mod_TS$coefficients))[i] %in% c(rownames(sum_TS$`Response TE`$coefficients))) {
                    pval_TE <-
                        c(pval_TE,
                          sum_TS$`Response TE`$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TE <- c(pval_TE, NA)
                }
            }
            
            cate <- rownames(mod_TS$coefficients)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub(":AlleleMinor", "", cate)
            cate <- gsub("AlleleMinor", "NoChangepoint", cate)
            cate <- gsub("Category", "", cate)
            pred_TS <-
                predict(
                    mod_TS,
                    newdata = data.frame(
                        "Category" = cate,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    interval = "confidence"
                )
            
            Model_DF1 <-
                data.frame(
                    "Coefficients" = rownames(mod_TS$coefficients),
                    "Chr" = rep(1, length(mod_TS$coefficients)),
                    "Direction" = c(rep(
                        "TS", length(mod_TS$coefficients) / 2
                    ), rep(
                        "TE", length(mod_TS$coefficients) / 2
                    )),
                    "Samplesize" = rep(SS$n, 2),
                    "Beta" = c(
                        unname(mod_TS$coefficients)[, 1],
                        unname(mod_TS$coefficients)[, 2]
                    ),
                    "P" = c(pval_TS, pval_TE),
                    "Allele" = c(rep("Major", length(
                        which(SS$Allele == "Major")
                    )), rep("Minor", length(
                        which(SS$Allele == "Minor")
                    )))
                )
            
            Pred_DF <- data.frame(
                "Category" = rep(c(cate), 2),
                "n" = c(rep(SS$n, 2)),
                "Direction" = c(rep(c(
                    "TS", "TE"
                ), each = length(
                    c(cate)
                ))),
                "Fit" = c(pred_TS[, 1], pred_TS[, 2]),
                "LB" = NA,
                "UB" = NA,
                "Allele" = c(rep("Major", length(
                    which(SS$Allele == "Major")
                )), rep("Minor", length(
                    which(SS$Allele == "Minor")
                )))
            )
            
            
            Data_Chr <- Model_DF1
            
        } else {
            newdata = data.frame(
                "Category" = c(unique(DF3[, c('Category', "Allele")])[1]),
                "TS" = c(rep(0, nrow(
                    unique(DF3[, c('Category', "Allele")])[1]
                ))),
                "TE" =  c(rep(0, nrow(
                    unique(DF3[, c('Category', "Allele")])[1]
                ))),
                "Allele" = c(unique(DF3[, c('Category', "Allele")])[2]),
                "Sample" = c(paste("Sample", 1:nrow(
                    unique(DF3[, c('Category', "Allele")])[1]
                )))
            )
            
            pred_TS <-
                MCMCglmm::predict.MCMCglmm(mod_TS, newdata = newdata, interval = "confidence")
            
            Model_DF <-
                rbind.data.frame(
                    data.frame(
                        "Coefficients" = rownames(sum_TS$solutions),
                        "Samplesize" = c(rep(SS1$n, 2)),
                        "Chr" = rep(1, length(rownames(
                            sum_TS$solutions
                        ))),
                        "Beta" = sum_TS$solutions[, 1],
                        "P" = sum_TS$solutions[, 5],
                        "Allele" = ifelse(grepl("AlleleMinor", rownames(sum_TS$solutions), fixed = TRUE), "Minor", "Major")
                    )
                )
            
            Pred_DF <-
                data.frame(
                    "Category" = rep(unlist(unique(DF[, c('Category', "Allele")])[1]), 2),
                    "n" = c(rep(SS1$n, 2)),
                    "Direction" = c(rep(c(
                        "TS", "TE"
                    ), each = length(SS1$n))),
                    "Allele" = rep(c(unlist(
                        unique(DF[, c('Category', "Allele")])[2]
                    )), 2),
                    "Fit" = c(pred_TS[, 1]),
                    "LB" = c(pred_TS[, 2]),
                    "UB" = c(pred_TS[, 3])
                )
            
            
            
            Data_Chr <- Model_DF
            rownames(Data_Chr) <- NULL
        }
    } else if (model_type == "multivariate" &
               category_spec == "No_Intercept") {
        if (model == "lm") {
            pval_TS <- c()
            j <- 1
            for (i in 1:length(rownames(mod_TS$coefficients))) {
                if (c(rownames(mod_TS$coefficients))[i] %in% c(rownames(
                    sum_TS$`Response TS`$coefficients
                ))) {
                    pval_TS <-
                        c(pval_TS,
                          sum_TS$`Response TS`$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TS <- c(pval_TS, NA)
                }
            }
            
            pval_TE <- c()
            j <- 1
            for (i in 1:length(rownames(mod_TS$coefficients))) {
                if (c(rownames(mod_TS$coefficients))[i] %in% c(rownames(sum_TS$`Response TE`$coefficients))) {
                    pval_TE <-
                        c(pval_TE,
                          sum_TS$`Response TE`$coefficients[j, 4])
                    j <- j + 1
                } else {
                    pval_TE <- c(pval_TE, NA)
                }
            }
            
            cate <- rownames(mod_TS$coefficients)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub(":AlleleMinor", "", cate)
            cate <- gsub("AlleleMinor", cate[1], cate)
            cate <- gsub("Category", "", cate)
            SS <- SS %>% filter(Category != "NoChangepoint")
            pred_TS <-
                predict(
                    mod_TS,
                    newdata = data.frame(
                        "Category" = cate,
                        "Allele" = c(rep("Major", length(
                            which(SS$Allele == "Major")
                        )), rep("Minor", length(
                            which(SS$Allele == "Minor")
                        )))
                    ),
                    interval = "confidence"
                )
            
            Model_DF1 <-
                data.frame(
                    "Coefficients" = rownames(mod_TS$coefficients),
                    "Chr" = rep(1, length(mod_TS$coefficients)),
                    "Direction" = c(rep(
                        "TS", length(mod_TS$coefficients) / 2
                    ), rep(
                        "TE", length(mod_TS$coefficients) / 2
                    )),
                    "Samplesize" = rep(SS$n, 2),
                    "Beta" = c(
                        unname(mod_TS$coefficients)[, 1],
                        unname(mod_TS$coefficients)[, 2]
                    ),
                    "P" = c(pval_TS, pval_TE),
                    "Allele" = c(rep("Major", length(
                        which(SS$Allele == "Major")
                    )), rep("Minor", length(
                        which(SS$Allele == "Minor")
                    )))
                )
            
            Pred_DF <- data.frame(
                "Category" = rep(c(cate), 2),
                "n" = c(rep(SS$n, 2)),
                "Direction" = c(rep(c(
                    "TS", "TE"
                ), each = length(
                    c(cate)
                ))),
                "Fit" = c(pred_TS[, 1], pred_TS[, 2]),
                "LB" = NA,
                "UB" = NA,
                "Allele" = c(rep("Major", length(
                    which(SS$Allele == "Major")
                )), rep("Minor", length(
                    which(SS$Allele == "Minor")
                )))
            )
            
            
            Data_Chr <- Model_DF1
            
        } else {
            newdata = data.frame(
                "Category" = c(unique(DF[, c('Category', "Allele")])[1]),
                "TS" = c(rep(0, nrow(
                    unique(DF[, c('Category', "Allele")])[1]
                ))),
                "TE" =  c(rep(0, nrow(
                    unique(DF[, c('Category', "Allele")])[1]
                ))),
                "Allele" = c(unique(DF[, c('Category', "Allele")])[2]),
                "Sample" = c(paste("Sample", 1:nrow(
                    unique(DF[, c('Category', "Allele")])[1]
                )))
            )
            
            pred_TS <-
                MCMCglmm::predict.MCMCglmm(mod_TS, newdata = newdata, interval = "confidence")
            
            
            Model_DF <-
                rbind.data.frame(
                    data.frame(
                        "Coefficients" = rownames(sum_TS$solutions),
                        "Chr" = rep(1, length(rownames(
                            sum_TS$solutions
                        ))),
                        "Samplesize" = c(rep(SS1$n, 2)),
                        "Beta" = sum_TS$solutions[, 1],
                        "P" = sum_TS$solutions[, 5],
                        "Allele" = ifelse(grepl("AlleleMinor", rownames(sum_TS$solutions), fixed = TRUE), "Minor", "Major")
                    )
                )
            
            Pred_DF <-
                data.frame(
                    "Category" = rep(unlist(unique(DF[, c('Category', "Allele")])[1]), 2),
                    "n" = c(rep(SS1$n, 2)),
                    "Direction" = c(rep(c(
                        "TS", "TE"
                    ), each = length(SS1$n))),
                    "Allele" = rep(c(unlist(
                        unique(DF[, c('Category', "Allele")])[2]
                    )), 2),
                    "Fit" = c(pred_TS[, 1]),
                    "LB" = c(pred_TS[, 2]),
                    "UB" = c(pred_TS[, 3])
                )
            
            
            
            Data_Chr <- Model_DF
            rownames(Data_Chr) <- NULL
        }
    }
    
    return(list(Data_Chr, Pred_DF))
}

## Apply Models
### Univariate
Uni_lm_7 <-
    Table_function_Pred(DF3,
                        model = "lm",
                        model_type = "univariate",
                        category_spec = "Intercept")
Uni_mcmc_7 <-
    Table_function_Pred(DF3,
                        model = "mcmc",
                        model_type = "univariate",
                        category_spec = "Intercept")
Uni_lm_6 <-
    Table_function_Pred(DF3,
                        model = "lm",
                        model_type = "univariate",
                        category_spec = "No_Intercept")
Uni_mcmc_6 <-
    Table_function_Pred(DF3,
                        model = "mcmc",
                        model_type = "univariate",
                        category_spec = "No_Intercept")

### Multivariate
Multi_lm_7 <-
    Table_function_Pred(DF3,
                        model = "lm",
                        model_type = "multivariate",
                        category_spec = "Intercept")
Multi_mcmc_7 <-
    Table_function_Pred(DF3,
                        model = "mcmc",
                        model_type = "multivariate",
                        category_spec = "Intercept")
Multi_lm_6 <-
    Table_function_Pred(DF3,
                        model = "lm",
                        model_type = "multivariate",
                        category_spec = "No_Intercept")
Multi_mcmc_6 <-
    Table_function_Pred(DF3,
                        model = "mcmc",
                        model_type = "multivariate",
                        category_spec = "No_Intercept")

### Save Tables + Figures
Results_List_NoNeut <-
    list(
        Uni_lm_7,
        Uni_mcmc_7,
        Uni_lm_6,
        Uni_mcmc_6,
        Multi_lm_7,
        Multi_mcmc_7,
        Multi_lm_6,
        Multi_mcmc_6
    )

Names_List_NoNeut <-
    c(
        "Univariate_lm_7",
        "Univariate_mcmc_7",
        "Univariate_lm_6",
        "Univariate_mcmc_6",
        "Multivariate_lm_7",
        "Multivariate_mcmc_7",
        "Multivariate_lm_6",
        "Multivariate_mcmc_6"
    )


gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 8
cols = gg_color_hue(n)

## Model - NoNeutral
for(i in 1:length(Results_List_NoNeut)){
    
    dat <- cbind.data.frame(Results_List_NoNeut[[i]][[1]], 
                            Results_List_NoNeut[[i]][[2]] %>% select(-c(n, Direction, Allele)))
    
    dat |>
        select(-c(Chr)) |> select(Coefficients, Allele, everything()) |> rename(n = Samplesize) |> filter(n > 0) |>
        gt() |>
        cols_align(
            align = "left",
            columns = c(Category, Coefficients)
        ) |>
        tab_header(
            title =  md("**Parameter Estimates and Confidence Intervals**")
        )  |>
        tab_style(
            style = cell_borders(
                sides = c("all"),
                color = "lightgrey",
                weight = px(1),
                style = "solid"
            ),
            locations = cells_body()
        ) |>
        tab_style(
            style = "padding-top:8px;padding-bottom:8px;padding-left:4px;padding-right:4px",
            locations = cells_body()
        ) |>
        tab_style(style = cell_text(weight = "bold"), locations = cells_column_labels()) |>
        fmt_scientific(
            columns = c("Beta", "P", "Fit", "LB", "UB"),
            decimals = 3
        ) |>
        opt_table_outline() %>% gtsave(paste(Names_List_NoNeut[[i]],  "_AD_Model_Pred.png", sep=""), 
                                       path = "../../tables/Chapter_5/", 
                                       vwidth = 1900, vheight = 1000)
    
    # Interval Plots 
    m1 <- Results_List_NoNeut[[i]][[2]] %>% mutate(Both = paste(Category, " (N = ", n,")", sep="")) %>% 
        mutate(Direction = factor(Direction, levels = c("TS", "TE"))) %>% filter(n > 0)
    
    Colors <-setNames( c(cols), c(unique(m1$Both)))
    
    png(filename = paste("../../figures/Chapter_5/", Names_List_NoNeut[[i]], "_AD_Interval.png", sep=""), width = 900, height = 500)
    
    if(i %in% c(5,7)){
       
        print(m1 %>% ggplot(aes(x=Both, y= Fit, color = Both)) +
                  geom_point(size = 5) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  facet_grid(Allele~Direction,  scales = "free")  + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + ggtitle("Plot of Parameter Point Estimates"))
        dev.off()
    } else {
        print(m1 %>% ggplot(aes(x=Both, y= Fit, ymin=LB, ymax=UB, color = Both)) +
                  geom_pointrange(position=position_dodge(width=0.6),  size = 1) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  facet_grid(Allele~Direction,  scales = "free")  + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + ggtitle("Interval Plot of Parameter Estimates"))
        dev.off()
    }
}