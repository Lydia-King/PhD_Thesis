# Create dataset to be used for downstream analysis
## Load up necessary libraries
library(dplyr)
library(fabricatr)

## Load up clinical data (clinical patient and sample file from cBioPortal 2019)
Patient <-
    read.delim(
        "../../data/METABRIC_2019/data_clinical_patient.txt",
        sep = "\t",
        na.strings = c("", " ", "NA"),
        comment.char = "#"
    )

Sample <-
    read.delim(
        "../../data/METABRIC_2019/data_clinical_sample.txt",
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
        OS = ifelse(OS_STATUS == "LIVING", 0, 1),
        DSS = ifelse(VITAL_STATUS == "Died of Disease", 1, 0)
    ) %>% dplyr::rename(., Her2_Status = HER2_STATUS)


Clinical_Sample_Luminal <-
    Clinical_Sample %>%
    filter(CLAUDIN_SUBTYPE %in% c("LumA", "LumB"))

## Load up CNA Data and get CNA Score for Luminal Patients
CNA <- read.delim("../../data/METABRIC_2019/data_CNA.txt", sep="\t", na.strings=c(""," ","NA"))

CNA_Score_Metrics_All <- data.frame("PATIENT_ID" = colnames(CNA[,3:ncol(CNA)]))
CNA_Score_Metrics_All$CNA_Score_All <- colSums(abs(CNA[,3:ncol(CNA)]), na.rm = T)

CNA_Score_Metrics_All$PATIENT_ID <-
    gsub("\\.", "-", CNA_Score_Metrics_All$PATIENT_ID)

## Calculate CNA Score Quartiles for Luminal A and B patients (n = 1175)
Luminal_Data <-
    merge(Clinical_Sample_Luminal, CNA_Score_Metrics_All, BY = "PATIENT_ID")

Luminal_Data$Subset_Quartile <-
    split_quantile(Luminal_Data$CNA_Score, type = 4)

## Write Out 
write.table(Luminal_Data, file = "../../data/Processed_Data/LuminalAB_Data.txt", sep = "\t", row.names = FALSE)
