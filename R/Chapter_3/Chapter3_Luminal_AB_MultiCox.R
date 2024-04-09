# Chapter 3: Multivariable Cox Models Luminal A/B

## Load up necessary libraries
library(tidyverse)
library(survival)

## Load up data
Luminal_Data <-
    read.delim("../../data/Processed_Data/LuminalAB_Data.txt", sep = "\t")  %>% 
    mutate_at(c("Subset_Quartile",
                "GRADE",
                "TUMOR_STAGE"),
              funs(as.factor(as.character(.))))

# Overall Survival
## Final Model with LNEP, Age, Claudin_Subtype, HER2, Size:
res.cox_OS <-
    coxph(
        Surv(OS_MONTHS, OS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            AGE_AT_DIAGNOSIS +
            CLAUDIN_SUBTYPE +
            Her2_Status +
            TUMOR_SIZE,
        data =  Luminal_Data
    )

summary(res.cox_OS)

## Add CNA Score
res.cox_OS <-
    coxph(
        Surv(OS_MONTHS, OS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            AGE_AT_DIAGNOSIS +
            CLAUDIN_SUBTYPE +
            Her2_Status +
            TUMOR_SIZE +
            CNA_Score_All,
        data =  Luminal_Data
    )

summary(res.cox_OS)

## Add CNA Quartiles
res.cox_OS <-
    coxph(
        Surv(OS_MONTHS, OS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            AGE_AT_DIAGNOSIS +
            CLAUDIN_SUBTYPE +
            Her2_Status +
            TUMOR_SIZE +
            Subset_Quartile,
        data =  Luminal_Data
    )

summary(res.cox_OS)


# Disease-specific Survival
# Final Model with LNEP, Age, Claudin_Subtype, HER2, Size:
res.cox_DSS <-
    coxph(
        Surv(OS_MONTHS, DSS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            AGE_AT_DIAGNOSIS +
            CLAUDIN_SUBTYPE +
            Her2_Status +
            TUMOR_SIZE,
        data =  Luminal_Data
    )

summary(res.cox_DSS)

## Add CNA Score
res.cox_DSS <-
    coxph(
        Surv(OS_MONTHS, DSS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            AGE_AT_DIAGNOSIS +
            CLAUDIN_SUBTYPE +
            Her2_Status +
            GRADE +
            TUMOR_SIZE +
            CNA_Score_All,
        data =  Luminal_Data
    )

summary(res.cox_DSS)

## Add CNA Quartiles
res.cox_DSS <-
    coxph(
        Surv(OS_MONTHS, DSS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            AGE_AT_DIAGNOSIS +
            CLAUDIN_SUBTYPE +
            Her2_Status +
            GRADE +
            TUMOR_SIZE +
            Subset_Quartile,
        data =  Luminal_Data
    )

summary(res.cox_DSS)

## Add interaction term
Luminal_Data <-
    within(Luminal_Data,
           CLAUDIN_SUBTYPE <-
               relevel(as.factor(CLAUDIN_SUBTYPE), ref = "LumB"))

res.cox_DSS <-
    coxph(
        Surv(OS_MONTHS, DSS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            TUMOR_SIZE +
            AGE_AT_DIAGNOSIS +
            Her2_Status +
            GRADE +
            CLAUDIN_SUBTYPE *
            CNA_Score_All,
        data = Luminal_Data
    )

summary(res.cox_DSS)
cox.zph(res.cox_DSS)


res.cox_DSS <-
    coxph(
        Surv(OS_MONTHS, DSS) ~
            LYMPH_NODES_EXAMINED_POSITIVE +
            TUMOR_SIZE +
            AGE_AT_DIAGNOSIS +
            Her2_Status +
            GRADE +
            CLAUDIN_SUBTYPE *
            Subset_Quartile,
        data = Luminal_Data
    )

summary(res.cox_DSS)
cox.zph(res.cox_DSS)