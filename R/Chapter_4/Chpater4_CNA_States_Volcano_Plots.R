# Chapter 4: CNA State Volcano Plots

## Load Libraries 
library(limma)
library(tidyverse)
library(dplyr)
library(operator.tools)
library(ggrepel)
library(patchwork)

## Create functions for alternative limma
glmfit_model_list <-
    function(y,
             formula = "y[,i] ~ 0",
             model_matrix_formula = NULL,
             CNA = NULL,
             data = Test,
             family = 'gaussian("identity")',
             weights = NULL) {
        # Set up vectors/lists
        fit <-
            Coef <-
            rank <-
            df.residual <-
            sigma <-
            Amean <-
            stdev.un <- std.error <- cov.coefficients <- design <- list()
        
        y <- as.data.frame(y)
        gene_name <- colnames(y)
        
        for (i in 1:ncol(y)) {
            # Coefficients
            if (is.null(weights)) {
                model <-
                    glm(as.formula(formula),
                        data = data,
                        family = eval(parse(text = family)))
            } else {
                model <-
                    glm(
                        as.formula(formula),
                        data = data,
                        family = eval(parse(text = family)),
                        weights = weights[i, ]
                    )
            } # Option to include weights from array weights but has to be in matrix form
            
            Coef_1 <- c(as.numeric(model$coefficients))
            names(Coef_1) <- names(model$coefficients)
            names(Coef_1) <- gsub("\\[, i]", "_", names(Coef_1))
            
            
            # Rank
            rank <- model$rank
            
            # df Residual
            df.residual <- model$df.residual
            
            # Sigma
            sigma <- sigma(model)
            
            # Amean (mean of expression for each gene)
            Amean <- mean(y[, i], na.rm = T)
            
            # Covariance matrix (Do both just to compare)
            cov.coefficients <- summary(model)$cov.unscaled
            colnames(cov.coefficients) <-
                gsub("\\[, i]", "_", colnames(cov.coefficients))
            rownames(cov.coefficients) <-
                gsub("\\[, i]", "_", rownames(cov.coefficients))
            
            # Standard deviation unscaled
            stdev.un <- sqrt(diag(cov.coefficients))
            # stdev.un_formula <- sqrt(diag(cov.coefficients_formula))
            stdev.unscaled <-
                matrix(stdev.un, 1, length(names(Coef_1)), byrow = TRUE)
            colnames(stdev.unscaled) <- names(Coef_1)
            colnames(stdev.unscaled) <-
                gsub("\\[, i]", "_", colnames(stdev.unscaled))
            
            # Standard Errors
            stderror <- sigma(model) * sqrt(diag(cov.coefficients))
            names(stderror) <- c(names(model$coefficients))
            names(stderror) <- gsub("\\[, i]", "_", names(stderror))
            
            stderror2 <- summary(model)$coefficients[, 2]
            names(stderror2) <- gsub("\\[, i]", "_", names(stderror2))
            
            # Design matrix
            if (!is.null(model_matrix_formula)) {
                design <- model.matrix(eval(parse(text = model_matrix_formula)))
            }
            
            # Fit
            fit1 <- list()
            fit1[["coefficients"]] <- as.matrix(t(Coef_1))
            fit1[["rank"]] <- rank
            fit1[["df.residual"]] <- df.residual
            fit1[["sigma"]] <- sigma
            fit1[["Amean"]] <- Amean
            fit1[["cov.coefficients"]] <- cov.coefficients
            fit1[["stdev.unscaled"]] <- stdev.unscaled
            fit1[["stderror"]] <- stderror
            fit1[["stderror2"]] <- stderror2
            
            if (!is.null(model_matrix_formula)) {
                fit1[["design"]] <- design
            }
            
            fit[[gene_name[[i]]]] <- fit1
            
        }
        
        # Fit
        fit
    }

expand.grid.unique <- function(x, y, include.equals = FALSE) {
    x <- unique(x)
    y <- unique(y)
    
    g <- function(i) {
        z <- setdiff(y, x[seq_len(i - include.equals)])
        
        if (length(z))
            cbind(x[i], z, deparse.level = 0)
    }
    
    do.call(rbind, lapply(seq_along(x), g))
}

Create_Fit_Contrasts <- function(model_fit_list) {
    cm_list <- fit_contrast_list <- list()
    
    for (i in 1:length(model_fit_list)) {
        con <-
            apply(expand.grid.unique(
                colnames(model_fit_list[[i]]$coefficients),
                colnames(model_fit_list[[i]]$coefficients)
            ),
            1,
            paste,
            collapse = " - ")
        
        cm_list[[i]] <-
            makeContrasts(contrasts = con,
                          levels = colnames(model_fit_list[[i]]$coefficients))
        
        fit_contrast_list[[names(model_fit_list[i])]] <-
            contrasts.fit(new("MArrayLM", model_fit_list[[i]]), cm_list[[i]])
    }
    
    fit_contrast_list
}

restructure_new <- function(fit, con){
    
    coefficients <- data.frame(matrix(ncol = length(con), nrow = 0))
    colnames(coefficients) <- con
    
    for(i in 1:length(fit)){
        coefficients <- plyr::rbind.fill(coefficients, as.data.frame(fit[[i]]$coefficients))}
    
    rownames(coefficients) <- names(fit)
    
    coefficients_1 <- c()
    
    for(i in 1:length(fit)){
        coefficients_1  <- rbind(coefficients_1, as.double(coefficients[i,]))
    }
    
    rownames(coefficients_1) <- names(fit)
    colnames(coefficients_1) <- con
    
    ## Standard Deviation Unscaled
    stdev.unscaled <- data.frame(matrix(ncol = length(con), nrow = 0))
    colnames(stdev.unscaled) <- con
    
    for(i in 1:length(fit)){
        stdev.unscaled <- plyr::rbind.fill(stdev.unscaled, 
                                           as.data.frame(fit[[i]]$stdev.unscaled))}
    
    rownames(stdev.unscaled) <- names(fit)
    
    stdev.unscaled_1 <- c()
    
    for(i in 1:length(fit)){
        stdev.unscaled_1  <- rbind(stdev.unscaled_1, 
                                   as.double(stdev.unscaled[i,]))
    }
    
    rownames(stdev.unscaled_1) <- names(fit)
    colnames(stdev.unscaled_1) <- con
    
    ## Sigma
    sigma <- c()
    for(i in 1:length(fit)){
        sigma <- c(sigma, fit[[i]]$sigma)
    }
    
    names(sigma) <- names(fit)
    
    ## Residual DF
    df.residual <- c()
    for(i in 1:length(fit)){
        df.residual <- c(df.residual, fit[[i]]$df.residual)
    }
    
    ## Amean 
    Amean <- c()
    for(i in 1:length(fit)){
        Amean <- c(Amean, fit[[i]]$Amean)
    }
    
    ## Cov.coefficients 
    allSame <- function(x) length(unique(x)) == 1
    
    cov.co <- list()
    
    for(i in 1:length(fit)){
        cov.co[[i]] <- fit[[i]]$cov.coefficients
    }
    
    restructured <- list()
    
    # Also need to consider design matrix
    if(!is.null(fit[[1]]$design)){
        design_list <- list()
        
        for(i in 1:length(fit)){
            design_list[[i]] <- fit[[i]]$design
        }
        
        vec <- unlist(lapply(design_list, is.fullrank))
        
        if(all(vec)){
            restructured$design <- design_list[[1]]
        }
    }
    
    
    restructured$coefficients <- coefficients_1
    restructured$stdev.unscaled <- stdev.unscaled_1
    restructured$sigma <- sigma
    restructured$df.residual <- df.residual
    restructured$Amean <- Amean 
    
    if(allSame(cov.co)){
        restructured$cov.coefficients <- cov.co[[1]]
    }
    
    # Convert to MArrayLM structure
    restructured <- new("MArrayLM", restructured)
    restructured
}

## Volcano Plot Function
Volcano_Plot_Function <- function(de, title){
    
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
        theme(legend.position = "right")
    
 #   ggsave(save_path_name, p, width = width, height = height)
    
    return(p)
}

## Load up Data
### GE data
GE_Loc <- read.delim("../../data/Processed_Data/GE_Common_Genes.txt", sep = "\t")

GE_Loc_NoDup <- GE_Loc
GE_Loc_NoDup$Amean <- rowMeans(GE_Loc_NoDup[,3:1982])
o <- order(GE_Loc_NoDup$Amean, decreasing=TRUE)
data.fit2 <- GE_Loc_NoDup[o,]
d <- duplicated(data.fit2$Hugo_Symbol)
data.fit2 <- data.fit2[!d,]

GE_Loc_NoDup <- data.fit2 %>% dplyr::select(-Amean)

## CNA data
CNA_Loc <- read.delim("../../data/Processed_Data/data_CNA_Loc_hg18.txt", sep=",")
CNA_Loc_NoDup <- CNA_Loc %>% filter(Hugo_Symbol %in% GE_Loc_NoDup$Hugo_Symbol)
CNA_Loc_NoDup <- CNA_Loc_NoDup %>% distinct(Hugo_Symbol, .keep_all = T)

CNA_Matrix_NoDup <- CNA_Loc_NoDup
rownames(CNA_Matrix_NoDup) <- CNA_Matrix_NoDup$Hugo_Symbol

CNA_Matrix_NoDup <- CNA_Matrix_NoDup %>% 
    dplyr::select(-X, -Hugo_Symbol, -Chromosome, -Arm, -txStart, -Entrez_Gene_Id) 

CNA_Matrix_NoDup <- as.data.frame(t(CNA_Matrix_NoDup))

## 3 level
CNA_Matrix_NoDup_3 <- CNA_Matrix_NoDup
CNA_Matrix_NoDup_3[CNA_Matrix_NoDup_3 == "-2"] <- "Del"
CNA_Matrix_NoDup_3[CNA_Matrix_NoDup_3 == "-1"] <- "Del"
CNA_Matrix_NoDup_3[CNA_Matrix_NoDup_3 == "0"] <- "Neut"
CNA_Matrix_NoDup_3[CNA_Matrix_NoDup_3 == "1"] <- "Amp"
CNA_Matrix_NoDup_3[CNA_Matrix_NoDup_3 == "2"] <- "Amp"

# 5 level
CNA_Matrix_NoDup_5 <- CNA_Matrix_NoDup
CNA_Matrix_NoDup_5[CNA_Matrix_NoDup_5 == "-2"] <- "HomoDel"
CNA_Matrix_NoDup_5[CNA_Matrix_NoDup_5 == "-1"] <- "HetDel"
CNA_Matrix_NoDup_5[CNA_Matrix_NoDup_5 == "0"] <- "Neut"
CNA_Matrix_NoDup_5[CNA_Matrix_NoDup_5 == "1"] <- "Gain"
CNA_Matrix_NoDup_5[CNA_Matrix_NoDup_5 == "2"] <- "Amp"

CNA_Matrix_NoDup_Chr_All <- CNA_Matrix_NoDup
CNA_Matrix_NoDup_3_Chr_All <- CNA_Matrix_NoDup_3
CNA_Matrix_NoDup_5_Chr_All <- CNA_Matrix_NoDup_5

Common_Genes_All <- intersect(colnames(CNA_Matrix_NoDup_Chr_All),
              GE_Loc_NoDup$Hugo_Symbol)

Common_Patients <- intersect(rownames(CNA_Matrix_NoDup_Chr_All), 
                             colnames(GE_Loc_NoDup))

CNA_Matrix_NoDup_3_Chr_All <- CNA_Matrix_NoDup_3_Chr_All %>%
    select(all_of(Common_Genes_All)) %>% filter(rownames(.) %in% Common_Patients)

CNA_Matrix_NoDup_5_Chr_All <- CNA_Matrix_NoDup_5_Chr_All %>%
    select(all_of(Common_Genes_All)) %>% filter(rownames(.) %in% Common_Patients)

## GE Data
Gene_Log_Chr_All <- GE_Loc_NoDup %>% 
    filter(Hugo_Symbol %in% Common_Genes_All) %>% 
    select(Hugo_Symbol, which(c(colnames(GE_Loc_NoDup) %in% Common_Patients)))

rownames(Gene_Log_Chr_All) <- Gene_Log_Chr_All$Hugo_Symbol
Gene_Log_Chr_All <- Gene_Log_Chr_All %>% select(-Hugo_Symbol)
Gene_Log_Chr_All <- as.data.frame(t(Gene_Log_Chr_All))

## Now Order the same 
### Now make sure they are ordered correctly 
Gene_Log_Chr_All_3 <-
    Gene_Log_Chr_All[match(rownames(CNA_Matrix_NoDup_3_Chr_All),
                           rownames(Gene_Log_Chr_All)), ]
Gene_Log_Chr_All_3 <-
    Gene_Log_Chr_All_3[, match(colnames(CNA_Matrix_NoDup_3_Chr_All),
                               colnames(Gene_Log_Chr_All_3))]

Gene_Log_Chr_All_5 <-
    Gene_Log_Chr_All[match(rownames(CNA_Matrix_NoDup_5_Chr_All),
                           rownames(Gene_Log_Chr_All)), ]
Gene_Log_Chr_All_5 <-
    Gene_Log_Chr_All_5[, match(colnames(CNA_Matrix_NoDup_5_Chr_All),
                               colnames(Gene_Log_Chr_All_5))]

# Model Fits
## Fit Models - 3 Levels
fit_CNA_3_2022_NInt <-  glmfit_model_list(Gene_Log_Chr_All_3, 
                                          CNA = CNA_Matrix_NoDup_3_Chr_All, 
                                          formula = "y[,i] ~ CNA[,i] + 0", 
                                          family='gaussian("identity")', 
                                          data = Gene_Log_Chr_All_3)

## Make and Fit Contrasts
CC_CNA_3_2022_NInt <- Create_Fit_Contrasts(fit_CNA_3_2022_NInt)

## Restructure 
Re_CNA_3_NInt_2022 <- restructure_new(CC_CNA_3_2022_NInt,
                                      con = c("CNA_Amp - CNA_Del", 
                                              "CNA_Amp - CNA_Neut", 
                                              "CNA_Del - CNA_Neut"))

## eBayes 
eBayes_CNA_3_NInt_2022 <- eBayes(Re_CNA_3_NInt_2022)

## Fit Models - 5 Levels
fit_CNA_5_2022_NInt <-  glmfit_model_list(
    Gene_Log_Chr_All_5,
    CNA = CNA_Matrix_NoDup_5_Chr_All,
    formula = "y[,i] ~ CNA[,i] + 0",
    family = 'gaussian("identity")',
    data = Gene_Log_Chr_All_5
)

## Make and Fit Contrasts
CC_CNA_5_2022_NInt <- Create_Fit_Contrasts(fit_CNA_5_2022_NInt)

## Restructure
Re_CNA_5_NInt_2022 <- restructure_new(
    CC_CNA_5_2022_NInt,
    con = c(
        "CNA_Amp - CNA_Gain",
        "CNA_Amp - CNA_HetDel",
        "CNA_Amp - CNA_HomoDel",
        "CNA_Amp - CNA_Neut",
        "CNA_Gain - CNA_HetDel",
        "CNA_Gain - CNA_HomoDel",
        "CNA_Gain - CNA_Neut",
        "CNA_HetDel - CNA_HomoDel",
        "CNA_HetDel - CNA_Neut",
        "CNA_HomoDel - CNA_Neut"
    )
)

## eBayes 
eBayes_CNA_5_NInt_2022 <- eBayes(Re_CNA_5_NInt_2022)

# Identify Small Sample Size Genes
## Small Sample Size Genes (3 state)
## Amp 
Genes <- c()
State <- c()
for(i in 1:ncol(CNA_Matrix_NoDup_3_Chr_All)){
    CNA_State <- CNA_Matrix_NoDup_3_Chr_All[,i]
    n <- sum(CNA_State == "Amp", na.rm = T)
    prop <- n/length(CNA_State)
    if(prop < 0.01){
        Genes<- c(Genes, colnames(CNA_Matrix_NoDup_3_Chr_All)[i])
        State <- c(State, "Amp")
    }
}

Amp_Genes_Small <- cbind.data.frame(Genes, State)

## Del 
Genes <- c()
State <- c()
for(i in 1:ncol(CNA_Matrix_NoDup_3_Chr_All)){
    CNA_State <- CNA_Matrix_NoDup_3_Chr_All[,i]
    n <- sum(CNA_State == "Del", na.rm = T)
    prop <- n/length(CNA_State)
    if(prop < 0.01){
        Genes<- c(Genes, colnames(CNA_Matrix_NoDup_3_Chr_All)[i])
        State <- c(State, "Del")
    }
}

Del_Genes_Small <- cbind.data.frame(Genes, State)

Amp_3_Small <- Amp_Genes_Small
Del_3_Small <- Del_Genes_Small

## Small Sample Size Genes (5 state)
## Amp 
Genes <- c()
State <- c()
for(i in 1:ncol(CNA_Matrix_NoDup_5_Chr_All)){
    CNA_State <- CNA_Matrix_NoDup_5_Chr_All[,i]
    n <- sum(CNA_State == "Amp", na.rm = T)
    prop <- n/length(CNA_State)
    if(prop < 0.01){
        Genes<- c(Genes, colnames(CNA_Matrix_NoDup_5_Chr_All)[i])
        State <- c(State, "Amp")
    }
}

Amp_Genes_Small<- cbind.data.frame(Genes, State)

## Gain
Genes <- c()
State <- c()
for(i in 1:ncol(CNA_Matrix_NoDup_5_Chr_All)){
    CNA_State <- CNA_Matrix_NoDup_5_Chr_All[,i]
    n <- sum(CNA_State == "Gain", na.rm = T)
    prop <- n/length(CNA_State)
    if(prop < 0.01){
        Genes<- c(Genes, colnames(CNA_Matrix_NoDup_5_Chr_All)[i])
        State <- c(State, "Gain")
    }
}

Gain_Genes_Small <- cbind.data.frame(Genes, State)

## HemiDel 
Genes <- c()
State <- c()
for(i in 1:ncol(CNA_Matrix_NoDup_5_Chr_All)){
    CNA_State <- CNA_Matrix_NoDup_5_Chr_All[,i]
    n <- sum(CNA_State == "HetDel", na.rm = T)
    prop <- n/length(CNA_State)
    if(prop < 0.01){
        Genes<- c(Genes, colnames(CNA_Matrix_NoDup_5_Chr_All)[i])
        State <- c(State, "HetDel")
    }
}

HetDel_Genes_Small <- cbind.data.frame(Genes, State)

## HomoDel 
Genes <- c()
State <- c()
for(i in 1:ncol(CNA_Matrix_NoDup_5_Chr_All)){
    CNA_State <- CNA_Matrix_NoDup_5_Chr_All[,i]
    n <- sum(CNA_State == "HomoDel", na.rm = T)
    prop <- n/length(CNA_State)
    if(prop < 0.01){
        Genes<- c(Genes, colnames(CNA_Matrix_NoDup_5_Chr_All)[i])
        State <- c(State, "HomoDel")
    }
}
HomoDel_Genes_Small <- cbind.data.frame(Genes, State)

# Volcano Plots (3 State)
de_Amp_Neut <- topTable(eBayes_CNA_3_NInt_2022, n=nrow(eBayes_CNA_3_NInt_2022), 
                        adjust.method = "BH", coef = 2) 

de_Amp_Neut_Small <- topTable(eBayes_CNA_3_NInt_2022, n=nrow(eBayes_CNA_3_NInt_2022), 
                              adjust.method = "BH", coef = 2) %>% 
    filter(rownames(.) %!in% c(unique(Amp_3_Small$Genes))) 

nrow(de_Amp_Neut_Small  %>% filter(abs(logFC) > 0.58 & de_Amp_Neut_Small$adj.P.Val < 0.05))
nrow(de_Amp_Neut_Small  %>% filter(logFC > 0.58 & de_Amp_Neut_Small$adj.P.Val < 0.05))
nrow(de_Amp_Neut_Small  %>% filter(logFC < -0.58 & de_Amp_Neut_Small$adj.P.Val < 0.05))

de_Del_Neut <- topTable(eBayes_CNA_3_NInt_2022, n=nrow(eBayes_CNA_3_NInt_2022), 
                        adjust.method = "BH", coef = 3)

de_Del_Neut_Small <- topTable(eBayes_CNA_3_NInt_2022, n=nrow(eBayes_CNA_3_NInt_2022), 
                              adjust.method = "BH", coef = 3) %>% 
    filter(rownames(.) %!in% c(unique(Del_3_Small$Genes))) 

nrow(de_Del_Neut_Small  %>% filter(abs(logFC) > 0.58 & de_Del_Neut_Small$adj.P.Val < 0.05))
nrow(de_Del_Neut_Small  %>% filter(logFC > 0.58 & de_Del_Neut_Small$adj.P.Val < 0.05))
nrow(de_Del_Neut_Small  %>% filter(logFC < -0.58 & de_Del_Neut_Small$adj.P.Val < 0.05))


write.table(de_Amp_Neut, "../../supplementary_information/DE_Genes/Amp_vs_Neut_3State.txt", row.names = TRUE)
write.table(de_Amp_Neut_Small, "../../supplementary_information/DE_Genes/Amp_vs_Neut_3State_Reduced.txt", row.names = TRUE)
write.table(de_Del_Neut, "../../supplementary_information/DE_Genes/Del_vs_Neut_3State.txt", row.names = TRUE)
write.table(de_Del_Neut_Small, "../../supplementary_information/DE_Genes/Del_vs_Neut_3State_Reduced.txt", row.names = TRUE)

vol_1 <- Volcano_Plot_Function(de_Amp_Neut, 
                           title = "(A) Volcano Plot for CNA Amp versus CNA Neutral")

vol_2 <- Volcano_Plot_Function(de_Amp_Neut_Small, 
                               title = "(B) Volcano Plot for CNA Amp versus CNA Neutral (Reduced)")

vol_3 <- Volcano_Plot_Function(de_Del_Neut, 
                               title = "(C) Volcano Plot for CNA Del versus CNA Neutral")

vol_4 <- Volcano_Plot_Function(de_Del_Neut_Small, 
                               title = "(D) Volcano Plot for CNA Del versus CNA Neutral (Reduced)")

layout <- "
AABB
CCDD
"

vol_1 + vol_2 + vol_3 + vol_4 + plot_layout(design = layout)  
ggsave("../../figures/Chapter_4/Volcano_CNA_3State.png", last_plot(), width = 17, height = 10)

# Volcano Plots (5 State)
de_Amp_Neut <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 4) 

de_Amp_Neut_Small <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 4) %>% 
    filter(rownames(.) %!in% c(unique(Amp_Genes_Small$Genes))) 

de_HetDel_Neut <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 9)

de_HetDel_Neut_Small <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 9) %>% 
    filter(rownames(.) %!in% c(unique(HetDel_Genes_Small$Genes))) 

de_Gain_Neut <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 7) 

de_Gain_Neut_Small <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 7) %>% 
    filter(rownames(.) %!in% c(unique(Gain_Genes_Small$Genes))) 

de_HomDel_Neut <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 10)

de_HomDel_Neut_Small <- topTable(eBayes_CNA_5_NInt_2022, n=nrow(eBayes_CNA_5_NInt_2022), adjust.method = "BH", coef = 10) %>% 
    filter(rownames(.) %!in% c(unique(HomoDel_Genes_Small$Genes))) 

write.table(de_Amp_Neut, "../../supplementary_information/DE_Genes/Amp_vs_Neut_5State.txt", row.names = TRUE)
write.table(de_Amp_Neut_Small, "../../supplementary_information/DE_Genes/Amp_vs_Neut_5State_Reduced.txt", row.names = TRUE)
write.table(de_HomDel_Neut, "../../supplementary_information/DE_Genes/HomoDel_vs_Neut_5State.txt", row.names = TRUE)
write.table(de_HomDel_Neut_Small, "../../supplementary_information/DE_Genes/HomoDel_vs_Neut_5State_Reduced.txt", row.names = TRUE)
write.table(de_Gain_Neut, "../../supplementary_information/DE_Genes/Gain_vs_Neut_5State.txt", row.names = TRUE)
write.table(de_Gain_Neut_Small, "../../supplementary_information/DE_Genes/Gain_vs_Neut_5State_Reduced.txt", row.names = TRUE)
write.table(de_HetDel_Neut, "../../supplementary_information/DE_Genes/HetDel_vs_Neut_5State.txt", row.names = TRUE)
write.table(de_HetDel_Neut_Small, "../../supplementary_information/DE_Genes/HetDel_vs_Neut_5State_Reduced.txt", row.names = TRUE)

vol_Gain <- Volcano_Plot_Function(de_Gain_Neut, 
                                  title = "(A) Volcano Plot for CNA Gain versus CNA Neutral")

vol_Gain_Small <- Volcano_Plot_Function(de_Gain_Neut_Small, 
                                        title = "(B) Volcano Plot for CNA Gain versus CNA Neutral (Reduced)")

vol_Amp <- Volcano_Plot_Function(de_Amp_Neut, 
                               title = "(C) Volcano Plot for CNA Amp versus CNA Neutral")

vol_Amp_Small <- Volcano_Plot_Function(de_Amp_Neut_Small, 
                               title = "(D) Volcano Plot for CNA Amp versus CNA Neutral (Reduced)")

vol_HetDel <- Volcano_Plot_Function(de_HetDel_Neut, 
                                  title = "(E) Volcano Plot for CNA HetDel versus CNA Neutral")

vol_HetDel_Small <- Volcano_Plot_Function(de_HetDel_Neut_Small, 
                                        title = "(F) Volcano Plot for CNA HetDel versus CNA Neutral (Reduced)")

vol_HomDel <- Volcano_Plot_Function(de_HomDel_Neut, 
                                 title = "(G) Volcano Plot for CNA HomoDel versus CNA Neutral")

vol_HomDel_Small <- Volcano_Plot_Function(de_HomDel_Neut_Small, 
                                       title = "(H) Volcano Plot for CNA HomoDel versus CNA Neutral (Reduced)")


layout <- "
AABB
CCDD
EEFF
GGHH
"

vol_Gain + vol_Gain_Small + vol_Amp + vol_Amp_Small + vol_HetDel + vol_HetDel_Small + vol_HomDel + vol_HomDel_Small + plot_layout(design = layout)
ggsave("../../figures/Chapter_4/Volcano_CNA_5State.png", last_plot(), width = 18, height = 18)

# Export Gene Sets 
rm(list=setdiff(ls(), c("eBayes_CNA_5_NInt_2022", 
                        "eBayes_CNA_3_NInt_2022", 
                        "Gain_Genes_Small", 
                        "Amp_Genes_Small", 
                        "HetDel_Genes_Small",
                        "HomoDel_Genes_Small", 
                        "Amp_3_Small",
                        "Del_3_Small")))

save.image(file='../../data/Chapter_4/CNA_State_Environment.RData')
