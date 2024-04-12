# Chapter 3: Create dataset to be used for downstream analysis

## Load up necessary libraries
library(dplyr)
library(fabricatr)

## Load up clinical data (clinical patient and sample file from cBioPortal 2021)
Patient <-
    read.delim(
        "../../data/METABRIC_2021/data_clinical_patient.txt",
        sep = "\t",
        na.strings = c("", " ", "NA"),
        comment.char = "#"
    )

Sample <-
    read.delim(
        "../../data/METABRIC_2021/data_clinical_sample.txt",
        sep = "\t",
        na.strings = c("", " ", "NA"),
        comment.char = "#"
    )

Clinical_Sample <-
    merge(Patient, Sample, by.x = "PATIENT_ID", by.y = "PATIENT_ID")

Clinical_Sample <-
    Clinical_Sample %>% mutate_at(
        c(
            "CELLULARITY",
            "CHEMOTHERAPY",
            "ER_IHC",
            "HER2_SNP6",
            "HORMONE_THERAPY",
            "INFERRED_MENOPAUSAL_STATE",
            "INTCLUST",
            "CLAUDIN_SUBTYPE",
            "THREEGENE",
            "LATERALITY",
            "RADIO_THERAPY",
            "HISTOLOGICAL_SUBTYPE",
            "BREAST_SURGERY",
            "CANCER_TYPE_DETAILED",
            "ER_STATUS",
            "HER2_STATUS",
            "GRADE",
            "PR_STATUS",
            "TUMOR_STAGE"
        ),
        funs(as.factor(as.character(.)))
    ) %>% mutate(
        OS = ifelse(OS_STATUS == "0:LIVING", 0, 1),
        DSS = ifelse(VITAL_STATUS == "Died of Disease", 1, 0)
    ) %>% dplyr::rename(., Her2_Status = HER2_STATUS)

## Write Out 
write.table(Clinical_Sample, file = "../../data/Processed_Data/METABRIC_PAM50_Data.txt", sep = "\t", row.names = FALSE)