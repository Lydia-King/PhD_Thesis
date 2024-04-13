# Chapter 5: Allele-Independent (AI) Models

## Load up Libraries
library(gt)
library(tidyverse)
library(MCMCglmm)

## Load Up Data
DF3 <-
    read.delim("../../data/Chapter_5/Example_Data_1_NoNeut.txt",
               sep = "\t")


DF3_Neut <-
    read.delim("../../data/Chapter_5/Example_Data_1_Neut.txt",
               sep = "\t")

## Univariate and Multivariate Models Function
Table_function_Pred <-
    function(Dataset,
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
                mod_TS <- lm(cbind(TS) ~ Category, data =  DF)
                mod_TE <- lm(cbind(TE) ~ Category, data =  DF)
                fit_model <- mod_TS
            } else {
                mod_TS <-
                    MCMCglmm(
                        cbind(TS) ~ Category,
                        data = DF,
                        family = c("gaussian"),
                        nitt = 60000,
                        thin = 10,
                        burnin = 15000
                    )
                mod_TE <-
                    MCMCglmm(
                        cbind(TE) ~ Category,
                        data = DF,
                        family = c("gaussian"),
                        nitt = 60000,
                        thin = 10,
                        burnin = 15000
                    )
                fit_model <- mod_TS
            }
        } else if (model_type == "univariate" &
                   category_spec == "No_Intercept") {
            if (model == "lm") {
                mod_TS <- lm(cbind(TS) ~ 0 + Category, data =  DF)
                mod_TE <-
                    lm(cbind(TE) ~ 0 + Category, data =  DF)
                fit_model <- mod_TS
            } else {
                mod_TS <-
                    MCMCglmm(
                        cbind(TS) ~ 0 + Category,
                        data = DF,
                        family = c("gaussian"),
                        nitt = 60000,
                        thin = 10,
                        burnin = 15000
                    )
                mod_TE <-
                    MCMCglmm(
                        cbind(TE) ~ 0 + Category,
                        data = DF,
                        family = c("gaussian"),
                        nitt = 60000,
                        thin = 10,
                        burnin = 15000
                    )
                fit_model <- mod_TS
            }
        } else if (model_type == "multivariate" &
                   category_spec == "Intercept") {
            if (model == "lm") {
                mod_TS <-
                    lm(cbind(TS, TE) ~ Category,
                       data =  DF)
                fit_model <- mod_TS
            } else {
                mod_TS <-
                    MCMCglmm(
                        cbind(TS, TE) ~ trait:Category,
                        data = DF,
                        family = c("gaussian", "gaussian"),
                        rcov = ~ us(trait):units,
                        nitt = 60000,
                        thin = 10,
                        burnin = 15000
                    )
                fit_model <- mod_TS
            }
        } else if (model_type == "multivariate" &
                   category_spec == "No_Intercept") {
            if (model == "lm") {
                mod_TS <-
                    lm(cbind(TS, TE) ~ 0 + Category,
                       data =  DF)
                fit_model <- mod_TS
            } else {
                mod_TS <-
                    MCMCglmm(
                        cbind(TS, TE) ~ 0 + trait:Category,
                        data = DF,
                        family = c("gaussian", "gaussian"),
                        rcov = ~ us(trait):units,
                        nitt = 60000,
                        thin = 10,
                        burnin = 15000
                    )
                fit_model <- mod_TS
            }
        }
        
        SS <- DF %>%  group_by(Category) %>% count()
        
        if (model_type == "univariate") {
            sum_TS <- summary(mod_TS)
            sum_TE <- summary(mod_TE)
        } else {
            sum_TS <- summary(mod_TS)
        }
        
        if (model_type == "univariate" & category_spec == "Intercept") {
            if (model == "lm") {
                cate <- rownames(sum_TS$coefficients)
                cate <- na.omit(gsub("Category", "", cate)[2:7])
                pred_TS <-
                    predict(mod_TS,
                            newdata = data.frame("Category" = c("NoChangepoint", cate)),
                            interval = "confidence")
                pred_TE <-
                    predict(mod_TE,
                            newdata = data.frame("Category" = c("NoChangepoint", cate)),
                            interval = "confidence")
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$coefficients),
                            "Chr" = rep(1, nrow(sum_TS$coefficients)),
                            "Direction" = rep("TS", nrow(sum_TS$coefficients)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$coefficients[, 1],
                            "P" = sum_TS$coefficients[, 4]
                        ),
                        data.frame(
                            "Coefficients" = rownames(sum_TE$coefficients),
                            "Chr" = rep(1, nrow(sum_TE$coefficients)),
                            "Direction" = rep("TE", nrow(sum_TE$coefficients)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TE$coefficients[, 1],
                            "P" = sum_TE$coefficients[, 4]
                        )
                    )
                
                Pred_DF <-
                    data.frame(
                        "Category" = rep(c("NoChangepoint", cate), 2),
                        "n" = c(rep(SS$n, 2)),
                        "Direction" = c(rep(c(
                            "TS", "TE"
                        ), each = length(
                            c("NoChangepoint", cate)
                        ))),
                        "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                        "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                        "UB" = c(pred_TS[, 3], pred_TE[, 3])
                    )
                
                Data_Chr <- Model_DF
                rownames(Data_Chr) <- NULL
            } else {
                pred_TS <-
                    MCMCglmm::predict.MCMCglmm(
                        mod_TS,
                        newdata = data.frame(
                            "Category" = c(rownames(sum_TS$solutions)) %>% gsub("Category", "", .),
                            "TS" = c(rep(0, length(
                                rownames(sum_TS$solutions)
                            )))
                        ),
                        interval = "confidence"
                    )
                pred_TE <-
                    MCMCglmm::predict.MCMCglmm(
                        mod_TE,
                        newdata = data.frame(
                            "Category" = c(rownames(sum_TE$solutions)) %>% gsub("Category", "", .),
                            "TE" = c(rep(0, length(
                                rownames(sum_TS$solutions)
                            )))
                        ),
                        interval = "confidence"
                    )
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$solutions),
                            "Chr" = rep(1, nrow(sum_TS$solutions)),
                            "Direction" = rep("TS", nrow(sum_TS$solutions)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$solutions[, 1],
                            "P" = sum_TS$solutions[, 5]
                        ),
                        data.frame(
                            "Coefficients" = rownames(sum_TE$solutions),
                            "Chr" = rep(1, nrow(sum_TE$solutions)),
                            "Direction" = rep("TE", nrow(sum_TE$solutions)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TE$solutions[, 1],
                            "P" = sum_TE$solutions[, 5]
                        )
                    )
                
                
                Pred_DF <-
                    data.frame(
                        "Category" = rep(c(rownames(
                            sum_TS$solutions
                        )) %>% gsub("Category", "", .), 2),
                        "n" = c(rep(SS$n, 2)),
                        "Direction" = c(rep(
                            c("TS", "TE"), each = nrow(sum_TS$solutions)
                        )),
                        "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                        "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                        "UB" = c(pred_TS[, 3], pred_TE[, 3])
                    )
                
                Data_Chr <- Model_DF
                Pred_DF$Category <-
                    gsub("\\(Intercept\\)",
                         "NoChangepoint",
                         Pred_DF$Category)
                Pred_DF <-
                    Pred_DF %>% mutate(Category = gsub("Category", "", Category)) %>% mutate(Category = factor(
                        Category,
                        levels = c(
                            "NoChangepoint",
                            "Neut/Amp",
                            "Neut/Del",
                            "Amp/Neut",
                            "Del/Neut",
                            "Amp/Del",
                            "Del/Amp"
                        )
                    ))
                rownames(Data_Chr) <- NULL
            }
        } else if (model_type == "univariate" &
                   category_spec == "No_Intercept") {
            if (model == "lm") {
                cate <- rownames(sum_TS$coefficients)
                cate <- gsub("Category", "", cate)
                pred_TS <-
                    predict(mod_TS,
                            newdata = data.frame("Category" = c(cate)),
                            interval = "confidence")
                pred_TE <-
                    predict(mod_TE,
                            newdata = data.frame("Category" = c(cate)),
                            interval = "confidence")
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$coefficients),
                            "Chr" = rep(1, nrow(sum_TS$coefficients)),
                            "Direction" = rep("TS", nrow(sum_TS$coefficients)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$coefficients[, 1],
                            "P" = sum_TS$coefficients[, 4]
                        ),
                        data.frame(
                            "Coefficients" = rownames(sum_TE$coefficients),
                            "Chr" = rep(1, nrow(sum_TE$coefficients)),
                            "Direction" = rep("TE", nrow(sum_TE$coefficients)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TE$coefficients[, 1],
                            "P" = sum_TE$coefficients[, 4]
                        )
                    )
                
                Pred_DF <- data.frame(
                    "Category" = rep(c(cate), 2),
                    "n" = c(rep(SS$n, 2)),
                    "Direction" = c(rep(c(
                        "TS", "TE"
                    ), each = 4)),
                    "Fit" = c(pred_TS[, 1], pred_TE[, 1]),
                    "LB" = c(pred_TS[, 2], pred_TE[, 2]),
                    "UB" = c(pred_TS[, 3], pred_TE[, 3])
                )
                
                Data_Chr <- Model_DF
                rownames(Data_Chr) <- NULL
                
            } else {
                pred_TS <-
                    MCMCglmm::predict.MCMCglmm(
                        mod_TS,
                        newdata = data.frame(
                            "Category" = c(rownames(sum_TS$solutions)) %>% gsub("Category", "", .),
                            "TS" = c(rep(0, length(
                                rownames(sum_TS$solutions)
                            )))
                        ),
                        interval = "confidence"
                    )
                pred_TE <-
                    MCMCglmm::predict.MCMCglmm(
                        mod_TE,
                        newdata = data.frame(
                            "Category" = c(rownames(sum_TE$solutions)) %>% gsub("Category", "", .),
                            "TE" = c(rep(0, length(
                                rownames(sum_TS$solutions)
                            )))
                        ),
                        interval = "confidence"
                    )
                
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$solutions),
                            "Chr" = rep(1, nrow(sum_TS$solutions)),
                            "Direction" = rep("TS", nrow(sum_TS$solutions)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$solutions[, 1],
                            "P" = sum_TS$solutions[, 5]
                        ),
                        data.frame(
                            "Coefficients" = rownames(sum_TE$solutions),
                            "Chr" = rep(1, nrow(sum_TE$solutions)),
                            "Direction" = rep("TE", nrow(sum_TE$solutions)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TE$solutions[, 1],
                            "P" = sum_TE$solutions[, 5]
                        )
                    )
                
                
                Pred_DF <-
                    data.frame(
                        "Category" = rep(c(rownames(
                            sum_TS$solutions
                        )) %>% gsub("Category", "", .), 2),
                        "n" = c(rep(SS$n, 2)),
                        "Direction" = c(rep(
                            c("TS", "TE"), each = nrow(sum_TS$solutions)
                        )),
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
                cate <- rownames(sum_TS$`Response TS`$coefficients)
                cate <- na.omit(gsub("Category", "", cate)[2:7])
                pred_TS <-
                    predict(mod_TS,
                            newdata = data.frame("Category" = c("NoChangepoint", cate)),
                            interval = "confidence")
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$`Response TS`$coefficients),
                            "Chr" = rep(
                                1,
                                nrow(sum_TS$`Response TS`$coefficients)
                            ),
                            "Direction" = rep(
                                "TS",
                                nrow(sum_TS$`Response TS`$coefficients)
                            ),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$`Response TS`$coefficients[, 1],
                            "P" = sum_TS$`Response TS`$coefficients[, 4]
                        ),
                        data.frame(
                            "Coefficients" = rownames(sum_TS$`Response TE`$coefficients),
                            "Chr" = rep(
                                1,
                                nrow(sum_TS$`Response TE`$coefficients)
                            ),
                            "Direction" = rep(
                                "TE",
                                nrow(sum_TS$`Response TE`$coefficients)
                            ),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$`Response TE`$coefficients[, 1],
                            "P" = sum_TS$`Response TE`$coefficients[, 4]
                        )
                    )
                
                
                Pred_DF <-
                    data.frame(
                        "Category" = rep(c("NoChangepoint", cate), 2),
                        "n" = c(rep(SS$n, 2)),
                        "Direction" = c(rep(c(
                            "TS", "TE"
                        ), each = length(
                            c("NoChangepoint", cate)
                        ))),
                        "Fit" = c(pred_TS[, 1], pred_TS[, 2]),
                        "LB" = NA,
                        "UB" = NA
                    )
                
                Data_Chr <- Model_DF
                rownames(Data_Chr) <- NULL
            } else {
                newdata = data.frame(
                    "Category" = c(unique(DF[, c('Category')])),
                    "TS" = as.numeric(c(rep(
                        0, length(unique(DF[, c('Category')]))
                    ))),
                    "TE" =  as.numeric(c(rep(
                        0, length(unique(DF[, c('Category')]))
                    )))
                )
                
                pred_TS <-
                    MCMCglmm::predict.MCMCglmm(mod_TS, newdata = newdata, interval = "confidence")
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$solutions),
                            "Chr" = rep(1, nrow(sum_TS$solutions)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$solutions[, 1],
                            "P" = sum_TS$solutions[, 5]
                        )
                    )
                
                newdata$id  <- 1:nrow(newdata)
                out <- merge(newdata, SS, by = c("Category"))
                out <- out[order(out$id),]
                
                Pred_DF <-
                    data.frame(
                        "Category" = rep(unlist(unique(DF[, c('Category')])), 2),
                        "n" = c(rep(out$n, 2)),
                        "Direction" = rep(c("TS", "TE"), each = nrow(unique(DF[, c('Category')]))),
                        "Fit" = c(pred_TS[, 1]),
                        "LB" = c(pred_TS[, 2]),
                        "UB" = c(pred_TS[, 3])
                    )
                
                Data_Chr <- Model_DF
                
                Pred_DF <-
                    Pred_DF %>% mutate(Category = gsub("Category", "", Category)) %>% mutate(
                        Category = gsub("Category", "", Category) %>% gsub("traitTE:", "", .) %>% gsub("traitTS:", "", .)
                    ) %>%
                    mutate(Category = factor(
                        Category,
                        levels = c(
                            "NoChangepoint",
                            "Neut/Amp",
                            "Neut/Del",
                            "Amp/Neut",
                            "Del/Neut",
                            "Amp/Del",
                            "Del/Amp"
                        )
                    ))
                rownames(Data_Chr) <- NULL
            }
        } else if (model_type == "multivariate" &
                   category_spec == "No_Intercept") {
            if (model == "lm") {
                cate <- rownames(sum_TS$`Response TS`$coefficients)
                cate <- gsub("Category", "", cate)
                pred_TS <-
                    predict(mod_TS,
                            newdata = data.frame("Category" = c(cate)),
                            interval = "confidence")
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$`Response TS`$coefficients),
                            "Chr" = rep(
                                1,
                                nrow(sum_TS$`Response TS`$coefficients)
                            ),
                            "Direction" = rep(
                                "TS",
                                nrow(sum_TS$`Response TS`$coefficients)
                            ),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$`Response TS`$coefficients[, 1],
                            "P" = sum_TS$`Response TS`$coefficients[, 4]
                        ),
                        data.frame(
                            "Coefficients" = rownames(sum_TS$`Response TE`$coefficients),
                            "Chr" = rep(
                                1,
                                nrow(sum_TS$`Response TE`$coefficients)
                            ),
                            "Direction" = rep(
                                "TE",
                                nrow(sum_TS$`Response TE`$coefficients)
                            ),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$`Response TE`$coefficients[, 1],
                            "P" = sum_TS$`Response TE`$coefficients[, 4]
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
                    "Fit" = c(pred_TS[, 1], pred_TS[, 2]),
                    "LB" = NA,
                    "UB" = NA
                )
                
                Data_Chr <- Model_DF
                rownames(Data_Chr) <- NULL
            } else {
                newdata = data.frame(
                    "Category" = c(unique(DF[, c('Category')])),
                    "TS" = as.numeric(c(rep(
                        0, length(unique(DF[, c('Category')]))
                    ))),
                    "TE" =  as.numeric(c(rep(
                        0, length(unique(DF[, c('Category')]))
                    )))
                )
                
                pred_TS <-
                    MCMCglmm::predict.MCMCglmm(mod_TS, newdata = newdata, interval = "confidence")
                
                Model_DF <-
                    rbind.data.frame(
                        data.frame(
                            "Coefficients" = rownames(sum_TS$solutions),
                            "Chr" = rep(1, nrow(sum_TS$solutions)),
                            "Samplesize" = SS$n,
                            "Beta" = sum_TS$solutions[, 1],
                            "P" = sum_TS$solutions[, 5]
                        )
                    )
                
                newdata$id  <- 1:nrow(newdata)
                out <- merge(newdata, SS, by = c("Category"))
                out <- out[order(out$id),]
                
                Pred_DF <-
                    data.frame(
                        "Category" = rep(unlist(unique(DF[, c('Category')])), 2),
                        "n" = c(rep(out$n, 2)),
                        "Direction" = rep(c("TS", "TE"), each = length(unique(DF[, c('Category')]))),
                        "Fit" = c(pred_TS[, 1]),
                        "LB" = c(pred_TS[, 2]),
                        "UB" = c(pred_TS[, 3])
                    )
                
                Data_Chr <- Model_DF
                
                Pred_DF <-
                    Pred_DF %>% mutate(Category = gsub("Category", "", Category)) %>% mutate(
                        Category = gsub("Category", "", Category) %>% gsub("traitTE:", "", .) %>% gsub("traitTS:", "", .)
                    ) %>%
                    mutate(Category = factor(
                        Category,
                        levels = c(
                            "NoChangepoint",
                            "Neut/Amp",
                            "Neut/Del",
                            "Amp/Neut",
                            "Del/Neut",
                            "Amp/Del",
                            "Del/Amp"
                        )
                    ))
                rownames(Data_Chr) <- NULL
            }
        }
        
        return(list(Data_Chr, Pred_DF, fit_model))
    }


## Apply Models
Uni_lm_7 <- Table_function_Pred(DF3, model = "lm", model_type = "univariate", category_spec = "Intercept")
Uni_mcmc_7 <- Table_function_Pred(DF3, model = "mcmc", model_type = "univariate", category_spec = "Intercept")
Uni_lm_6 <- Table_function_Pred(DF3, model = "lm", model_type = "univariate", category_spec = "No_Intercept")
Uni_mcmc_6 <- Table_function_Pred(DF3, model = "mcmc", model_type = "univariate", category_spec = "No_Intercept")

Uni_lm_7_Neut <- Table_function_Pred(DF3_Neut, model = "lm", model_type = "univariate", category_spec = "Intercept")
Uni_mcmc_7_Neut <- Table_function_Pred(DF3_Neut, model = "mcmc", model_type = "univariate", category_spec = "Intercept")
Uni_lm_6_Neut <- Table_function_Pred(DF3_Neut, model = "lm", model_type = "univariate", category_spec = "No_Intercept")
Uni_mcmc_6_Neut <- Table_function_Pred(DF3_Neut, model = "mcmc", model_type = "univariate", category_spec = "No_Intercept")

Multi_lm_7 <- Table_function_Pred(DF3, model = "lm", model_type = "multivariate", category_spec = "Intercept")
Multi_mcmc_7 <- Table_function_Pred(DF3, model = "mcmc", model_type = "multivariate", category_spec = "Intercept")
Multi_lm_6 <- Table_function_Pred(DF3, model = "lm", model_type = "multivariate", category_spec = "No_Intercept")
Multi_mcmc_6 <- Table_function_Pred(DF3, model = "mcmc", model_type = "multivariate", category_spec = "No_Intercept")

Multi_lm_7_Neut <- Table_function_Pred(DF3_Neut, model = "lm", model_type = "multivariate", category_spec = "Intercept")
Multi_mcmc_7_Neut <- Table_function_Pred(DF3_Neut, model = "mcmc", model_type = "multivariate", category_spec = "Intercept")
Multi_lm_6_Neut <- Table_function_Pred(DF3_Neut, model = "lm", model_type = "multivariate", category_spec = "No_Intercept")
Multi_mcmc_6_Neut <- Table_function_Pred(DF3_Neut, model = "mcmc", model_type = "multivariate", category_spec = "No_Intercept")

## Save Tables + Figures
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


Results_List_Neut <-
    list(
        Uni_lm_7_Neut,
        Uni_mcmc_7_Neut,
        Uni_lm_6_Neut,
        Uni_mcmc_6_Neut,
        Multi_lm_7_Neut,
        Multi_mcmc_7_Neut,
        Multi_lm_6_Neut,
        Multi_mcmc_6_Neut
    )

Names_List_Neut <-
    c(
        "Univariate_lm_7_Neut",
        "Univariate_mcmc_7_Neut",
        "Univariate_lm_6_Neut",
        "Univariate_mcmc_6_Neut",
        "Multivariate_lm_7_Neut",
        "Multivariate_mcmc_7_Neut",
        "Multivariate_lm_6_Neut",
        "Multivariate_mcmc_6_Neut"
    )

## Outputs
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

n = 8
cols = gg_color_hue(n)

## Model - NoNeutral
for(i in 1:length(Results_List_NoNeut)){

    Results_List_NoNeut[[i]][[1]] |>
        select(-c(Chr)) |> rename(n = Samplesize) |>
        gt() |>
        tab_header(
            title =  md("**(A) Parameter Estimates**")
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
            columns = c("Beta", "P"),
            decimals = 3
        ) |>
        opt_table_outline() %>% gtsave(paste(Names_List_NoNeut[[i]],  "_AI_Model.png", sep=""), 
                                       path = "../../tables/Chapter_5/", 
                                       vwidth = 1900, vheight = 1000)
    
    ## Pred
    Results_List_NoNeut[[i]][[2]] |>
        gt() |>
        tab_header(
            title =  md("**(A) Parameter Estimates and Confidence Intervals**")
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
            columns = c(4,5,6),
            decimals = 3
        ) |> fmt_number(
            columns = c(3),
            decimals = 0,
            use_seps = TRUE
        )  |>
        opt_table_outline() %>% tab_options(table.width = pct(42)) %>% 
        gtsave(paste(Names_List_NoNeut[[i]],  "_AI_Pred.png", sep=""), 
               path = "../../tables/Chapter_5/", 
               vwidth = 1900, vheight = 1000)
    
    
    # Interval Plots 
    
    m1 <- Results_List_NoNeut[[i]][[2]] %>% mutate(Both = paste(Category, " (N = ", n,")", sep="")) %>% 
        mutate(Direction = factor(Direction, levels = c("TS", "TE")))
    Colors <-setNames( c(cols), c(unique(m1$Both)))

        if(i %in% c(5,7)){
            png(filename = paste("../../figures/Chapter_5/", Names_List_NoNeut[[i]], "_AI_Interval.png", sep=""), width = 800, height = 400)
            print(m1 %>% ggplot(aes(x=Both, y= Fit, color = Both)) +
                      geom_point(size = 5) +
                      geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                      coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                      scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                      theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                      theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                      theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                      facet_grid(~Direction,  scales = "free") + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + ggtitle("(A) Plot of Parameter Point Estimates"))
            dev.off()
        } else {
            png(filename = paste("../../figures/Chapter_5/", Names_List_NoNeut[[i]], "_AI_Interval.png", sep=""), width = 800, height = 400)
        print(m1 %>% ggplot(aes(x=Both, y= Fit, ymin=LB, ymax=UB, color = Both)) +
                  geom_pointrange(position=position_dodge(width=0.6),  size = 1) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  facet_grid(~Direction,  scales = "free") + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + ggtitle("(A) Interval Plot of Parameter Estimates"))
        dev.off()
        }
}


## Model - Neutral
for(i in 1:length(Results_List_Neut)){
    
    Results_List_Neut[[i]][[1]] |>
        select(-c(Chr)) |> rename(n = Samplesize) |>
        gt() |>
        tab_header(
            title =  md("**(B) Parameter Estimates**")
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
            columns = c("Beta", "P"),
            decimals = 3
        ) |>
        opt_table_outline() %>% gtsave(paste(Names_List_Neut[[i]],  "_AI_Model.png", sep=""), 
                                       path = "../../tables/Chapter_5/", 
                                       vwidth = 1900, vheight = 1000)
    
    ## Pred
    Results_List_Neut[[i]][[2]] |>
        gt() |>
        cols_align(
            align = "left",
            columns = Category
        ) |>
        tab_header(
            title =  md("**(B) Parameter Estimates and Confidence Intervals**")
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
            columns = c(4,5,6),
            decimals = 3
        ) |> fmt_number(
            columns = c(3),
            decimals = 0,
            use_seps = TRUE
        )  |>
        opt_table_outline() %>% tab_options(table.width = pct(42)) %>% 
        gtsave(paste(Names_List_Neut[[i]],  "_AI_Pred.png", sep=""), 
               path = "../../tables/Chapter_5/", vwidth = 1900, vheight = 1000)
    
    
    # Interval Plots 
    
    m1 <- Results_List_Neut[[i]][[2]] %>% mutate(Both = paste(Category, " (N = ", n,")", sep="")) %>% 
        mutate(Direction = factor(Direction, levels = c("TS", "TE")))
    Colors <-setNames( c(cols), c(unique(m1$Both)))
    
    if(i %in% c(5,7)){
        png(filename = paste("../../figures/Chapter_5/", Names_List_Neut[[i]], "_AI_Interval.png", sep=""), width = 800, height = 400)
        print(m1 %>% ggplot(aes(x=Both, y= Fit, color = Both)) +
                  geom_point(size = 5) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  facet_grid(~Direction,  scales = "free") + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + ggtitle("(B) Plot of Parameter Point Estimates"))
        dev.off()
    } else {
        png(filename = paste("../../figures/Chapter_5/", Names_List_Neut[[i]], "_AI_Interval.png", sep=""), width = 800, height = 400)
        print(m1 %>% ggplot(aes(x=Both, y= Fit, ymin=LB, ymax=UB, color = Both)) +
                  geom_pointrange(position=position_dodge(width=0.6),  size = 1) +
                  geom_hline(yintercept = 0, linetype=2, size = 0.8)+
                  coord_flip() +  xlab('') + ylab("Estimates") + theme(plot.title = element_text(hjust = 0.5, size = 24, face="bold")) +
                  scale_color_manual(values = Colors) + theme(axis.text.y = element_text(size = 15)) +  theme(axis.text.x = element_text(size = 13)) +
                  theme(legend.text = element_text(size = 25))+ theme(strip.text.x = element_text(size = 20)) + theme(legend.title = element_text(size=22)) +
                  theme(axis.title.x = element_text(size = 18)) + theme(plot.subtitle=element_text(hjust=0.5, size= 25, face="bold")) +
                  theme(legend.position = c("top"), legend.direction = "horizontal") + theme(axis.title.x = element_text(size=18), legend.spacing.x = unit(0.6, 'cm')) +
                  facet_grid(~Direction,  scales = "free") + theme(strip.text = element_text(size = 20)) + theme(legend.position = "none") + ggtitle("(B) Interval Plot of Parameter Estimates"))
        dev.off()
    }
}
