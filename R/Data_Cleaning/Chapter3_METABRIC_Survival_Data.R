# Chapter 3: Data Preprocessing

## Indicate where data is located
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


## Combine Clinical and Sample into one dataset
Clinical_Sample <-
    merge(Patient, Sample, by.x = "PATIENT_ID", by.y = "PATIENT_ID")

Clinical_Sample$PATIENT_ID <-
    gsub("\\-", ".", Clinical_Sample$PATIENT_ID)


Clinical_Sample$CLAUDIN_SUBTYPE[Clinical_Sample$CLAUDIN_SUBTYPE == "NC"] <- NA
Clinical_Sample$CLAUDIN_SUBTYPE[Clinical_Sample$CLAUDIN_SUBTYPE == "Her2"] <- "HER2"
Clinical_Sample$CLAUDIN_SUBTYPE[Clinical_Sample$CLAUDIN_SUBTYPE == "claudin-low"] <- "Claudin-low"

Clinical_Sample$CLAUDIN_SUBTYPE <-
    as.factor(as.character(Clinical_Sample$CLAUDIN_SUBTYPE))

names(Clinical_Sample)[16] <- "PAM50"

Clinical_Sample$INTCLUST <-
    factor(
        Clinical_Sample$INTCLUST,
        levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10")
    )

Clinical_Sample <-
    Clinical_Sample %>% mutate(
        OS = ifelse(OS_STATUS == "0:LIVING", 0, 1),
        DSS = ifelse(VITAL_STATUS == "Died of Disease", 1, 0)
    ) %>% mutate(FiveYearTimeOS = OS_MONTHS, FiveYearOS = OS_STATUS) %>% mutate(
        FiveYearTimeOS = case_when(
            OS_STATUS == "0:LIVING" & OS_MONTHS > 60 ~ 60,
            OS_STATUS == "1:DECEASED" & OS_MONTHS > 60 ~ 60,
            OS_STATUS == "0:LIVING" & OS_MONTHS <= 60 ~ OS_MONTHS,
            OS_STATUS == "1:DECEASED" & OS_MONTHS <= 60 ~ OS_MONTHS,
        )
    ) %>%  mutate(
        FiveYearOS = case_when(
            OS_STATUS == "0:LIVING" & OS_MONTHS > 60 ~ "LIVING",
            OS_STATUS == "1:DECEASED" & OS_MONTHS > 60 ~ "LIVING",
            OS_STATUS == "0:LIVING" & OS_MONTHS <= 60 ~ OS_STATUS,
            OS_STATUS == "1:DECEASED" &
                OS_MONTHS <= 60 ~ OS_STATUS
        )
    ) %>% mutate(FiveYearOS = ifelse(FiveYearOS == "LIVING", 0, 1)) %>%
    mutate(FiveYearTimeDSS = OS_MONTHS, FiveYearDSS = VITAL_STATUS) %>%
    mutate(
        FiveYearTimeDSS = case_when(
            VITAL_STATUS == "Living" & OS_MONTHS > 60 ~ 60,
            VITAL_STATUS == "Died of Disease" & OS_MONTHS > 60 ~ 60,
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS > 60 ~ 60,
            VITAL_STATUS == "Living" & OS_MONTHS <= 60 ~ OS_MONTHS,
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS <= 60 ~ OS_MONTHS,
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS <= 60 ~ OS_MONTHS
        )
    ) %>% mutate(
        FiveYearDSS = case_when(
            VITAL_STATUS == "Living" & OS_MONTHS > 60 ~ "Living",
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS > 60 ~ "Living",
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS > 60 ~ "Living",
            VITAL_STATUS == "Living" &
                OS_MONTHS <= 60 ~ VITAL_STATUS,
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS <= 60 ~ VITAL_STATUS,
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS <= 60 ~ VITAL_STATUS
        )
    ) %>% mutate(FiveYearDSS = ifelse(FiveYearDSS == "Died of Disease", 1, 0)) %>%
    mutate(TenYearTimeOS = OS_MONTHS, TenYearOS = OS_STATUS) %>% mutate(
        TenYearTimeOS = case_when(
            OS_STATUS == "0:LIVING" & OS_MONTHS > 120 ~ 120,
            OS_STATUS == "1:DECEASED" & OS_MONTHS > 120 ~ 120,
            OS_STATUS == "0:LIVING" &
                OS_MONTHS <= 120 ~ OS_MONTHS,
            OS_STATUS == "1:DECEASED" &
                OS_MONTHS <= 120 ~ OS_MONTHS,
        )
    ) %>%  mutate(
        TenYearOS = case_when(
            OS_STATUS == "0:LIVING" & OS_MONTHS > 120 ~ "LIVING",
            OS_STATUS == "1:DECEASED" &
                OS_MONTHS > 120 ~ "LIVING",
            OS_STATUS == "0:LIVING" &
                OS_MONTHS <= 120 ~ OS_STATUS,
            OS_STATUS == "1:DECEASED" &
                OS_MONTHS <= 120 ~ OS_STATUS
        )
    ) %>% mutate(TenYearOS = ifelse(TenYearOS == "LIVING", 0, 1)) %>%
    mutate(TenYearTimeDSS = OS_MONTHS, TenYearDSS = VITAL_STATUS) %>%
    mutate(
        TenYearTimeDSS = case_when(
            VITAL_STATUS == "Living" & OS_MONTHS > 120 ~ 120,
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS > 120 ~ 120,
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS > 120 ~ 120,
            VITAL_STATUS == "Living" &
                OS_MONTHS <= 120 ~ OS_MONTHS,
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS <= 120 ~ OS_MONTHS,
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS <= 120 ~ OS_MONTHS
        )
    ) %>% mutate(
        TenYearDSS = case_when(
            VITAL_STATUS == "Living" & OS_MONTHS > 120 ~ "Living",
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS > 120 ~ "Living",
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS > 120 ~ "Living",
            VITAL_STATUS == "Living" &
                OS_MONTHS <= 120 ~ VITAL_STATUS,
            VITAL_STATUS == "Died of Disease" &
                OS_MONTHS <= 120 ~ VITAL_STATUS,
            VITAL_STATUS == "Died of Other Causes" &
                OS_MONTHS <= 120 ~ VITAL_STATUS
        )
    ) %>% mutate(TenYearDSS = ifelse(TenYearDSS == "Died of Disease", 1, 0))

## Write Out Data
write.table(Clinical_Sample,
            file = "../../data/Processed_Data/Chapter3_METABRIC_PAM50_Data.txt",
            sep = "\t",
            row.names = FALSE)
