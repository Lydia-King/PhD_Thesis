# Chapter 3: Survival Association Luminal A/B 9Cox)

## Load up necessary libraries
library(tidyverse)
library(fabricatr)
library(RColorBrewer)
library(survival)
library(operator.tools)
library(gt)

## Load up data
Luminal_Data <-
    read.delim("../../data/Processed_Data/LuminalAB_Data.txt", sep = "\t")

Luminal_Data <- Luminal_Data %>% mutate_at(c("Subset_Quartile",
                                             "GRADE",
                                             "TUMOR_STAGE"),
                                           funs(as.factor(as.character(.))))

# Cox Models for Quartiles
## Overall Survival
res.cox <-
    coxph(Surv(OS_MONTHS, OS) ~ Subset_Quartile, data = Luminal_Data)
summary(res.cox)

## Disease Specific Survival
res.cox <-
    coxph(Surv(OS_MONTHS, DSS) ~ Subset_Quartile, data = Luminal_Data)
summary(res.cox)


# Cox Models for Score
## Overall Survival
res.cox <-
    coxph(Surv(OS_MONTHS, OS) ~ CNA_Score_All, data = Luminal_Data)

summary(res.cox)


## Disease Specific Survival
res.cox <-
    coxph(Surv(OS_MONTHS, DSS) ~ CNA_Score_All, data = Luminal_Data)
summary(res.cox)


# Identify Variables Associated with Survival
colnames(Luminal_Data)
UniVarList <-
    c(
        "LYMPH_NODES_EXAMINED_POSITIVE",
        "NPI",
        "CELLULARITY",
        "CHEMOTHERAPY",
        "ER_IHC",
        "HER2_SNP6",
        "HORMONE_THERAPY",
        "INFERRED_MENOPAUSAL_STATE",
        "INTCLUST",
        "AGE_AT_DIAGNOSIS",
        "CLAUDIN_SUBTYPE",
        "THREEGENE",
        "LATERALITY",
        "RADIO_THERAPY",
        "HISTOLOGICAL_SUBTYPE",
        "BREAST_SURGERY",
        "CANCER_TYPE_DETAILED",
        "ER_STATUS",
        "Her2_Status",
        "GRADE",
        "PR_STATUS",
        "TUMOR_SIZE",
        "TUMOR_STAGE"
    )

## OS
univ_formulas <-
    sapply(UniVarList, function(x)
        as.formula(paste('Surv(OS_MONTHS, OS)~', x)))

univ_models <-
    lapply(univ_formulas, function(x) {
        summary(coxph(x, data = Luminal_Data))
    })

Variable <- P.val_LRT <- P.val_Wald <- c()

for (i in 1:23) {
    res.cox <- univ_models[i]
    n <- res.cox[1]
    m <-  n[[UniVarList[i]]]
    pval_LRT <- m$logtest[[3]]
    pval_Wald <- m$waldtest[[3]]
    Variable <- c(Variable, UniVarList[i])
    P.val_LRT <- c(P.val_LRT, pval_LRT)
    P.val_Wald <- c(P.val_Wald, pval_Wald)
}

UV_OS_All <- cbind.data.frame(Variable, P.val_LRT, P.val_Wald)

UV_OS_All$P.val.Adj_LRT <-
    p.adjust(UV_OS_All$P.val_LRT,
             method = "BH",
             n = length(UV_OS_All$P.val_LRT))

UV_OS_All$P.val.Adj_Wald <-
    p.adjust(UV_OS_All$P.val_Wald,
             method = "BH",
             n = length(UV_OS_All$P.val_Wald))

UV_OS_All <- UV_OS_All[order(-UV_OS_All$P.val.Adj_Wald), ]
UV_OS_All$C <- "Overall Survival"

Association_OS_All <-
    UV_OS_All[UV_OS_All$P.val.Adj_LRT < 0.05 |
                  UV_OS_All$P.val.Adj_Wald < 0.05, ]

A_OS_All <- as.character(Association_OS_All$Variable)

UniVarList[which(UniVarList %!in% Association_OS_All$Variable)]

## Not Associated with survival: ER_STATUS, CELLULARITY, LATERALITY, 
## CHEMOTHERAPY, CANCER_TYPE_DETAILED

## Create Table
UV_OS_All <-  dplyr::rename(
    UV_OS_All,
    "Clinical Variable" = Variable,
    "LRT" = P.val_LRT,
    "Wald Test" = P.val_Wald,
    "Adjusted LRT" = P.val.Adj_LRT,
    "Adjusted Wald Test" = P.val.Adj_Wald
)

UV_OS_All <- UV_OS_All %>%
    mutate_if(is.numeric, round, digits = 2)

UV_OS_All[UV_OS_All < 0.0001] <- "<0.0001"

UV_OS_All$`Clinical Variable` <-
    dplyr::recode(
        as.character(UV_OS_All$`Clinical Variable`),
        "BREAST_SURGERY" = "Breast Surgery",
        "CHEMOTHERAPY" = "Chemotherapy",
        "CLAUDIN_SUBTYPE" = "PAM50",
        "ER_IHC" = "ER Immunohistochemistry",
        "GRADE" = "Histological Grade",
        "HER2_SNP6" = "HER2 SNP6",
        "HER2_STATUS" = "HER2 Status",
        "HISTOLOGICAL_SUBTYPE" = "Histological Subtype",
        "HORMONE_THERAPY" = "Hormone Therapy",
        "INFERRED_MENOPAUSAL_STATE" = "Inferred Menopausal State",
        "INTCLUST" = "Integrative Cluster",
        "PR_STATUS" = "PR Status",
        "THREEGENE" = "Three Gene Classification",
        "TUMOR_STAGE" = "Clinical Stage",
        "AGE_AT_DIAGNOSIS" = "Age at Diagnosis",
        "LYMPH_NODES_EXAMINED_POSITIVE" = "Positive Lymph Nodes",
        "NPI" = "NPI",
        "TUMOR_SIZE" = "Tumour Size",
        "RADIO_THERAPY" = "Radiotherapy",
        "ER_STATUS" = "ER Status",
        "CANCER_TYPE_DETAILED" = "Cancer Type Detailed",
        "CELLULARITY" = "Cellularity",
        "LATERALITY" = "Laterality",
        "Her2_Status" =  "HER2 Status"
    )

UV_OS_All |> select(-c(C)) |>
    gt() |>
    tab_header(
        title =  md(
            "**Univariate Cox models for each clinical variable for OS, within the Luminal METABRIC cohort**"
        )
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
    tab_style(style = "padding-top:11px;padding-bottom:11px;padding-left:7px;padding-right:7px",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5))) |>
    fmt_scientific(columns = c(2, 3, 4, 5),
                   decimals = 3) |>
    cols_align(align = "right",
               columns = c(2, 3, 4, 5)) |>
    opt_table_outline() %>% tab_options(table.width = pct(30)) %>% gtsave(
        "Luminal_Survival_Association_OS.png",
        path = "../../tables/Chapter_3/",
        vwidth = 3800,
        vheight = 2000
    )

## DSS
univ_formulas <-
    sapply(UniVarList, function(x)
        as.formula(paste('Surv(OS_MONTHS, DSS)~', x)))

univ_models <-
    lapply(univ_formulas, function(x) {
        summary(coxph(x, data = Luminal_Data))
    })

Variable <- P.val_LRT <- P.val_Wald <- c()

for (i in 1:23) {
    res.cox <- univ_models[i]
    n <- res.cox[1]
    m <-  n[[UniVarList[i]]]
    pval_LRT <- m$logtest[[3]]
    pval_Wald <- m$waldtest[[3]]
    Variable <- c(Variable, UniVarList[i])
    P.val_LRT <- c(P.val_LRT, pval_LRT)
    P.val_Wald <- c(P.val_Wald, pval_Wald)
}

UV_DSS_All <- cbind.data.frame(Variable, P.val_LRT, P.val_Wald)

UV_DSS_All$P.val.Adj_LRT <-
    p.adjust(UV_DSS_All$P.val_LRT,
             method = "BH",
             n = length(UV_DSS_All$P.val_LRT))

UV_DSS_All$P.val.Adj_Wald <-
    p.adjust(UV_DSS_All$P.val_Wald,
             method = "BH",
             n = length(UV_DSS_All$P.val_Wald))

UV_DSS_All <- UV_DSS_All[order(-UV_DSS_All$P.val.Adj_Wald), ]

UV_DSS_All$C <- "Disease-specific Survival"


Association_DSS_All <-
    UV_DSS_All[UV_DSS_All$P.val.Adj_LRT < 0.05 |
                   UV_DSS_All$P.val.Adj_Wald < 0.05, ]

A_DSS_All <- as.character(Association_DSS_All$Variable)

UniVarList[which(UniVarList %!in% Association_DSS_All$Variable)]

## Not Associated with survival: CELLULARITY, LATERALITY, RADIO_THERAPY, 
## CANCER_TYPE_DETAILED, ER_STATUS

UV_DSS_All <-  dplyr::rename(
    UV_DSS_All,
    "Clinical Variable" = Variable,
    "LRT" = P.val_LRT,
    "Wald Test" = P.val_Wald,
    "Adjusted LRT" = P.val.Adj_LRT,
    "Adjusted Wald Test" = P.val.Adj_Wald
)

UV_DSS_All <- UV_DSS_All %>%
    mutate_if(is.numeric, round, digits = 2)

UV_DSS_All[UV_DSS_All < 0.0001] <- "<0.0001"


## Create Table
UV_DSS_All$`Clinical Variable` <-
    dplyr::recode(
        as.character(UV_DSS_All$`Clinical Variable`),
        "BREAST_SURGERY" = "Breast Surgery",
        "CHEMOTHERAPY" = "Chemotherapy",
        "CLAUDIN_SUBTYPE" = "PAM50",
        "ER_IHC" = "ER Immunohistochemistry",
        "GRADE" = "Histological Grade",
        "HER2_SNP6" = "HER2 SNP6",
        "HER2_STATUS" = "HER2 Status",
        "HISTOLOGICAL_SUBTYPE" = "Histological Subtype",
        "HORMONE_THERAPY" = "Hormone Therapy",
        "INFERRED_MENOPAUSAL_STATE" = "Inferred Menopausal State",
        "INTCLUST" = "Integrative Cluster",
        "PR_STATUS" = "PR Status",
        "THREEGENE" = "Three Gene Classification",
        "TUMOR_STAGE" = "Clinical Stage",
        "AGE_AT_DIAGNOSIS" = "Age at Diagnosis",
        "LYMPH_NODES_EXAMINED_POSITIVE" = "Positive Lymph Nodes",
        "NPI" = "NPI",
        "TUMOR_SIZE" = "Tumour Size",
        "RADIO_THERAPY" = "Radiotherapy",
        "ER_STATUS" = "ER Status",
        "CANCER_TYPE_DETAILED" = "Cancer Type Detailed",
        "CELLULARITY" = "Cellularity",
        "LATERALITY" = "Laterality",
        "Her2_Status" =  "HER2 Status"
        
    )

UV_DSS_All |> select(-c(C)) |>
    gt() |>
    tab_header(
        title =  md(
            "**Univariate Cox models for each clinical variable for DSS, within the Luminal METABRIC cohort**"
        )
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
    tab_style(style = "padding-top:11px;padding-bottom:11px;padding-left:7px;padding-right:7px",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5))) |>
    fmt_scientific(columns = c(2, 3, 4, 5),
                   decimals = 3) |>  cols_align(align = "right",
                                                columns = c(2, 3, 4, 5)) |>
    opt_table_outline() %>% tab_options(table.width = pct(30)) %>% gtsave(
        "Luminal_Survival_Association_DSS.png",
        path = "../../tables/Chapter_3/",
        vwidth = 3800,
        vheight = 2000
    )
