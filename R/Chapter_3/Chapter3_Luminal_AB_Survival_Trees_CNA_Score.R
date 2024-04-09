# Chapter 3: CNA Score Survival Trees Luminal A/B

## Load up necessary libraries
library(tidyverse)
library(survival)
library(rpart)
library(partykit)
library(rpart.plot)

## Load up data
Luminal_Data <-
    read.delim("../../data/Processed_Data/LuminalAB_Data.txt", sep = "\t")  %>%
    mutate_at(
        c(
            "Subset_Quartile",
            "GRADE",
            "TUMOR_STAGE",
            "CLAUDIN_SUBTYPE",
            "Her2_Status"
        ),
        funs(as.factor(as.character(.)))
    ) %>% dplyr::rename(
        "LNEP" = "LYMPH_NODES_EXAMINED_POSITIVE",
        "Age" = "AGE_AT_DIAGNOSIS",
        "Subtype" = "CLAUDIN_SUBTYPE",
        "Tumour_Size" = "TUMOR_SIZE",
        "Her2" = "Her2_Status",
        "CNA_Score" = "CNA_Score_All",
        "Tumour_Grade" = "GRADE",
        "CNA_Quartile" = "Subset_Quartile"
    )

## CNA Score with Clinical variables OS and DSS
### Rpart Trees
pfit_OS <-
    rpart(
        Surv(OS_MONTHS, OS) ~  
            LNEP + 
            Age + 
            Her2 + 
            Tumour_Size + 
            Tumour_Grade + 
            Subtype + 
            CNA_Score,
        data =  Luminal_Data,
        method = "exp"
    )

pfit_DSS <-
    rpart(
        Surv(OS_MONTHS, DSS) ~  
            LNEP + 
            Age + 
            Her2 + 
            Tumour_Size + 
            Tumour_Grade + 
            Subtype + 
            CNA_Score,
        data =  Luminal_Data,
        method = "exp"
    )

tfit_OS <- as.party(pfit_OS)
tfit_DSS <- as.party(pfit_DSS)

#png(
#    filename = "../../figures/Chapter_3/LuminalAB_Rpart_OS_Score.png",
#    width = 440,
#    height = 330,
#    units = "px",
#    bg = "white"
#)
#plot(tfit_OS)
#invisible(dev.off())

png(
    filename = "../../figures/Appendices/Appendix_A/LuminalAB_Rpart_DSS_Score.png",
    width = 650,
    height = 450,
    units = "px",
    bg = "white"
)
plot(tfit_DSS)
invisible(dev.off())

### Ctree Trees
set.seed(1004)
Ctree_OS <-
    ctree(
        Surv(OS_MONTHS, OS) ~ 
            LNEP + 
            Age + 
            Her2 + 
            Tumour_Size + 
            Tumour_Grade + 
            Subtype + 
            CNA_Score,
        data = Luminal_Data
    )

set.seed(1004)
Ctree_DSS <-
    ctree(
        Surv(OS_MONTHS, DSS) ~ 
            LNEP + 
            Age + 
            Her2 + 
            Tumour_Size + 
            Tumour_Grade + 
            Subtype + 
            CNA_Score,
        data = Luminal_Data
    )

#png(
#    filename = "../../figures/Chapter_3/LuminalAB_Ctree_OS_Score.png",
#    width = 440,
#    height = 500,
#    units = "px",
#    bg = "white"
#)
#plot(Ctree_OS)
#invisible(dev.off())

png(
    filename = "../../figures/Chapter_3/LuminalAB_Ctree_DSS_Score.png",
    width = 1125,
    height = 600,
    units = "px",
    bg = "white"
)
plot(Ctree_DSS, gp = gpar(fontsize = 16))
invisible(dev.off())
