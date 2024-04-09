# Chapter 3: Association Analysis Luminal A/B

## Load up necessary libraries
library(tidyverse)
library(gt)

## Load up data
Luminal_Data <-
    read.delim("../../data/Processed_Data/LuminalAB_Data.txt", sep = "\t")

Luminal_Data <- Luminal_Data %>% mutate_at(c("Subset_Quartile",
                                             "GRADE",
                                             "TUMOR_STAGE"),
                                           funs(as.factor(as.character(.))))

# Association Tests
## Chi-squared
Chi <-
    c(
        "ER_IHC",
        "HER2_SNP6",
        "HORMONE_THERAPY",
        "INFERRED_MENOPAUSAL_STATE",
        "INTCLUST",
        "CLAUDIN_SUBTYPE",
        "THREEGENE",
        "RADIO_THERAPY",
        "HISTOLOGICAL_SUBTYPE",
        "BREAST_SURGERY",
        "CANCER_TYPE_DETAILED",
        "Her2_Status",
        "GRADE",
        "PR_STATUS",
        "TUMOR_STAGE",
        "CHEMOTHERAPY"
    )

Variable <- P.val_Chi <- c()


for (i in 1:length(Chi)) {
    C <- chisq.test(Luminal_Data[, Chi[i]],
                    Luminal_Data$Subset_Quartile,
                    correct = FALSE)
    Variable <- c(Variable, Chi[i])
    pval <- C$p.value
    P.val_Chi <- c(P.val_Chi, pval)
}

Chi_S <- cbind.data.frame(Variable, P.val_Chi)

## Fishers Exact Test
Fish_name <-
    c(
        "ER_IHC",
        "HORMONE_THERAPY",
        "INFERRED_MENOPAUSAL_STATE",
        "RADIO_THERAPY",
        "BREAST_SURGERY",
        "Her2_Status",
        "PR_STATUS",
        "CHEMOTHERAPY"
    )

Fish_Sim <-
    c(
        "CLAUDIN_SUBTYPE",
        "INTCLUST",
        "GRADE",
        "TUMOR_STAGE",
        "HER2_SNP6",
        "THREEGENE",
        "HISTOLOGICAL_SUBTYPE",
        "CANCER_TYPE_DETAILED"
    )

Variable <- P.val_FE <- c()

for (i in 1:length(Fish_name)) {
    C <- fisher.test(Luminal_Data[, Fish_name[i]],
                     Luminal_Data$Subset_Quartile)
    Variable <- c(Variable, Fish_name[i])
    pval <- C$p.value
    P.val_FE <- c(P.val_FE, pval)
}

Fish_name <- cbind.data.frame(Variable, P.val_FE)

Variable <- P.val_FE <- c()

for (i in 1:length(Fish_Sim)) {
    C <- fisher.test(Luminal_Data[, Fish_Sim[i]],
                     Luminal_Data$Subset_Quartile,
                     simulate.p.value = T)
    Variable <- c(Variable, Fish_Sim[i])
    pval <- C$p.value
    P.val_FE <- c(P.val_FE, pval)
}

Fish_Sim <- cbind.data.frame(Variable, P.val_FE)
Fish_All <- rbind.data.frame(Fish_name, Fish_Sim)

## Kruskal-Wallis
KW <- c("LYMPH_NODES_EXAMINED_POSITIVE",
        "NPI",
        "AGE_AT_DIAGNOSIS",
        "TUMOR_SIZE")

Variable <- P.val_KW <- c()

for (i in 1:length(KW)) {
    K <-
        kruskal.test(Luminal_Data[, KW[i]]  ~ Luminal_Data$Subset_Quartile)
    Variable <- c(Variable, KW[i])
    pval <- K$p.value
    P.val_KW <- c(P.val_KW, pval)
}

KW1 <- cbind.data.frame(Variable, P.val_KW)

## Create Dataframe with all Results
Tot <- merge(Chi_S, Fish_All, by = "Variable")

Tot$P.val_KW <- c("-")
KW1$P.val_FE <- c("-")
KW1$P.val_Chi <- c("-")

Total_Association_All <- rbind.data.frame(Tot, KW1)

Total_Association_All$P.val.adj <-
    p.adjust(c(
        Total_Association_All$P.val_Chi[1:16],
        Total_Association_All$P.val_KW[17:20]
    ),
    "BH",
    n = 20)

Total_Association_All <-  dplyr::rename(
    Total_Association_All,
    "Clinical Variable" = Variable,
    "Chi-Squared Test" = P.val_Chi,
    "Fisher's Exact Test" = P.val_FE,
    "Kruskal-Wallis Test" = P.val_KW,
    "Adjusted P-value" = P.val.adj
)

## Create Table
Total_Association_All$`Clinical Variable` <-
    dplyr::recode(
        as.character(Total_Association_All$`Clinical Variable`),
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

Total_Association_All$`Chi-Squared Test` <-
    as.numeric(Total_Association_All$`Chi-Squared Test`)

Total_Association_All$`Fisher's Exact Test` <-
    as.numeric(Total_Association_All$`Fisher's Exact Test`)

Total_Association_All$`Kruskal-Wallis Test` <-
    as.numeric(Total_Association_All$`Kruskal-Wallis Test`)

Total_Association_All$`Adjusted P-value` <-
    as.numeric(Total_Association_All$`Adjusted P-value`)

Total_Association_All <-
    Total_Association_All [order(-Total_Association_All$`Adjusted P-value`), ]

Total_Association_All <- Total_Association_All %>%
    mutate_if(is.numeric, round, digits = 2)

Total_Association_All[Total_Association_All < 0.0001] <- "<0.0001"

Total_Association_All |>
    gt() |>
    tab_header(
        title =  md(
            "**Association between categorised CNA Score metric quartiles
                    and Clinical variables, within the Luminal METABRIC cohort**"
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
    tab_style(style = "padding-top:11px;padding-bottom:11px;padding-left:7px;padding-right:7",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4, 5))) |>
    fmt_scientific(columns = c(2, 3, 4, 5),
                   decimals = 3) |> cols_align(align = "right",
                                               columns = c(2, 3, 4, 5)) |>
    opt_table_outline() %>% tab_options(table.width = pct(30)) %>% gtsave(
        "CNA_Quartile_Clin_Association.png",
        path = "../../tables/Chapter_3/",
        vwidth = 3800,
        vheight = 2000
    )

    
## CNA Score
## Kruskal-Wallis
KW <- c(
    "ER_IHC",
    "HER2_SNP6",
    "HORMONE_THERAPY",
    "INFERRED_MENOPAUSAL_STATE",
    "INTCLUST",
    "CLAUDIN_SUBTYPE",
    "THREEGENE",
    "RADIO_THERAPY",
    "HISTOLOGICAL_SUBTYPE",
    "BREAST_SURGERY",
    "CANCER_TYPE_DETAILED",
    "Her2_Status",
    "GRADE",
    "PR_STATUS",
    "TUMOR_STAGE",
    "CHEMOTHERAPY"
)

kruskal.test(Luminal_Data$CNA_Score_All ~ Luminal_Data$CELLULARITY)
kruskal.test(Luminal_Data$CELLULARITY ~ Luminal_Data$CNA_Score_All)

? kruskal.test()

Variable <- P.val_KW <- c()

for (i in 1:length(KW)) {
    K <-
        kruskal.test(Luminal_Data$CNA_Score_All ~ Luminal_Data[, KW[i]])
    Variable <- c(Variable, KW[i])
    pval <- K$p.value
    P.val_KW <- c(P.val_KW, pval)
}

KW1 <- cbind.data.frame(Variable, P.val_KW)

# Correlation
cor_var <- c("LYMPH_NODES_EXAMINED_POSITIVE",
             "NPI",
             "AGE_AT_DIAGNOSIS",
             "TUMOR_SIZE")

Variable <- P.val_cor <- c()

for (i in 1:length(cor_var)) {
    print(cor_var[i])
    K <-
        cor.test(Luminal_Data$CNA_Score_All, Luminal_Data[, cor_var[i]],
                 use = 'complete.obs')
    Variable <- c(Variable, cor_var[i])
    pval <- K$p.value
    P.val_cor <- c(P.val_cor, pval)
}

cor1 <- cbind.data.frame(Variable, P.val_cor)
cor1$P.val_KW <- ("-")

## Create Dataframe with all Results
Tot <- KW1

Tot$P.val_cor <- c("-")

Total_Association_All <- rbind.data.frame(Tot, cor1)

Total_Association_All$P.val.adj <-
    c(p.adjust(
        c(
            Total_Association_All$P.val_KW[1:16],
            Total_Association_All$P.val_cor[17:20]
        ),
        "BH",
        n = 20
    ))

Total_Association_All <-  dplyr::rename(
    Total_Association_All,
    "Clinical Variable" = Variable,
    "Kruskal-Wallis Test" = P.val_KW,
    "Pearson's Correlation" = P.val_cor,
    "Adjusted P-value" = P.val.adj
)

## Create Table
Total_Association_All$`Clinical Variable` <-
    dplyr::recode(
        as.character(Total_Association_All$`Clinical Variable`),
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

Total_Association_All$`Kruskal-Wallis Test` <-
    as.numeric(Total_Association_All$`Kruskal-Wallis Test`)

Total_Association_All$`Pearson's Correlation` <-
    as.numeric(Total_Association_All$`Pearson's Correlation`)

Total_Association_All$`Adjusted P-value` <-
    as.numeric(Total_Association_All$`Adjusted P-value`)

Total_Association_All <-
    Total_Association_All [order(-Total_Association_All$`Adjusted P-value`), ]

Total_Association_All <- Total_Association_All %>%
    mutate_if(is.numeric, round, digits = 2)

Total_Association_All[Total_Association_All < 0.0001] <- "<0.0001"

Total_Association_All |>
    gt() |>
    tab_header(
        title =  md(
            "**Association between CNA Score metric and Clinical variables,
                    within the Luminal METABRIC cohort**"
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
    tab_style(style = "padding-top:11px;padding-bottom:11px;padding-left:7px;padding-right:7",
              locations = cells_body()) |>
    tab_style(style = cell_text(weight = "bold"),
              locations = cells_column_labels(columns = c(1, 2, 3, 4))) |>
    fmt_scientific(columns = c(2, 3, 4),
                   decimals = 3) |> cols_align(align = "right",
                                               columns = c(2, 3, 4)) |>
    opt_table_outline() %>% tab_options(table.width = pct(30)) %>%
    gtsave(
        "CNA_Score_Clin_Association.png",
        path = "../../tables/Chapter_3/",
        vwidth = 3800,
        vheight = 2000
    )
