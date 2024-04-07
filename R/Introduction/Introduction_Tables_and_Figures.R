# Chapter 1: Introduction 
## Libraries
library(tidyverse)
library(gt)
library(gtsummary)

## Load up clinical and sample data
patient_data <-
    read.delim(
        "../../data/METABRIC_2021/data_clinical_patient.txt",
        sep = "\t",
        comment.char = "#",
        na.strings = c("", "NA")
    )
sample_data <-
    read.delim(
        "../../data/METABRIC_2021/data_clinical_sample.txt",
        sep = "\t",
        comment.char = "#",
        na.strings = c("", "NA")
    )
clinical_data <- merge(patient_data, sample_data, by = "PATIENT_ID")

# Tidy up clinical data
clinical_data <- clinical_data %>% mutate(CLAUDIN_SUBTYPE = ifelse(CLAUDIN_SUBTYPE == "NC", NA, CLAUDIN_SUBTYPE))

## Table 1 - Clinical Characteristics
Table_1 <- clinical_data %>%
    select(AGE_AT_DIAGNOSIS, NPI, LYMPH_NODES_EXAMINED_POSITIVE, TUMOR_SIZE, ER_STATUS,
           PR_STATUS, HER2_STATUS, TUMOR_STAGE, GRADE, CLAUDIN_SUBTYPE, INTCLUST) %>%
    mutate(CLAUDIN_SUBTYPE = case_when(CLAUDIN_SUBTYPE  == "Basal" ~ "Basal",
                                       CLAUDIN_SUBTYPE  == "claudin-low" ~ "Claudin-low",
                                       CLAUDIN_SUBTYPE  == "Her2" ~ "HER2",
                                       CLAUDIN_SUBTYPE  == "LumA" ~ "Luminal A",
                                       CLAUDIN_SUBTYPE  == "LumB" ~ "Luminal B",
                                       CLAUDIN_SUBTYPE  == "Normal" ~ "Normal")) %>%
    mutate(INTCLUST = factor(INTCLUST, levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10"))) %>%
    rename("Histological Grade" = GRADE,
           "Tumour Stage" = TUMOR_STAGE,
           "Lymph Nodes Positive" = LYMPH_NODES_EXAMINED_POSITIVE,
           "ER Status" = ER_STATUS,
           "PR Status" = PR_STATUS,
           "Age" = AGE_AT_DIAGNOSIS,
           "Tumour Size" = TUMOR_SIZE,
           "HER2 Status" = HER2_STATUS,
           "PAM50" = CLAUDIN_SUBTYPE,
           "IntClust" = INTCLUST) %>%
    tbl_summary(missing_text = "NA") %>%
    bold_labels() %>%
    modify_header(label ~ "**Clinical Characteristics**") %>%
    tbl_split(., c("Tumour Stage"))

as_gt(Table_1[[1]]) %>% gtsave("Table1_Clin_pt1.png", path = "../../tables/Introduction/")
as_gt(Table_1[[2]]) %>% gtsave("Table1_Clin_pt2.png", path = "../../tables/Introduction/")

## Table 2 - Treatment Characteristics
treatment_data <- clinical_data %>%
    select(HORMONE_THERAPY, CHEMOTHERAPY, RADIO_THERAPY, BREAST_SURGERY) %>%
    mutate("HORMONE_THERAPY" = ifelse(HORMONE_THERAPY == "YES", "Yes", ifelse(HORMONE_THERAPY == "NO", "No", HORMONE_THERAPY))) %>%
    mutate("CHEMOTHERAPY" = ifelse(CHEMOTHERAPY == "YES", "Yes", ifelse(CHEMOTHERAPY == "NO", "No", CHEMOTHERAPY))) %>%
    mutate("RADIO_THERAPY" = ifelse(RADIO_THERAPY == "YES", "Yes", ifelse(RADIO_THERAPY == "NO", "No", RADIO_THERAPY))) %>%
    mutate("BREAST_SURGERY" = ifelse(BREAST_SURGERY == "BREAST CONSERVING", "Breast Conserving",
                                     ifelse(BREAST_SURGERY == "MASTECTOMY", "Mastectomy", BREAST_SURGERY))) %>%
    rename("Hormone Therapy" = HORMONE_THERAPY,
           "Chemotherapy" = CHEMOTHERAPY,
           "Radiotherapy" = RADIO_THERAPY,
           "Breast Surgery" = BREAST_SURGERY) %>%
    tbl_summary(missing_text = "NA",  type = all_dichotomous() ~ "categorical") %>%
    bold_labels() %>%
    modify_header(label ~ "**Treatment Characteristics**")

as_gt(treatment_data) %>% gtsave("Table2_Treatment.png", path = "../../tables/Introduction/")

## Table 3 - Survival Characteristics
surv_data <- clinical_data %>%
    select(OS_MONTHS, OS_STATUS, VITAL_STATUS, RFS_MONTHS, RFS_STATUS) %>%
    mutate("OS_STATUS" = ifelse(OS_STATUS == "0:LIVING", "Living", "Deceased")) %>%
    mutate(RFS_STATUS = ifelse(RFS_STATUS == "0:Not Recurred", "Not Recurred", "Recurred")) %>%
    rename("Survival Time (Months)" = OS_MONTHS,
           "Overall Survival" = OS_STATUS,
           "Disease-specific Survival" = VITAL_STATUS,
           "Recurrence-free Survival" = RFS_STATUS,
           "Recurrence Time (Months)" = RFS_MONTHS) %>%
    tbl_summary(missing_text = "NA",  type = all_dichotomous() ~ "categorical" ) %>%
    bold_labels() %>%
    modify_header(label ~ "**Survival Characteristics**")

as_gt(surv_data) %>% gtsave("Table3_Survival.png", path = "../../tables/Introduction/")


# Plot IntClust on x-axis and display percentage of IntClust group that is made up of each PAM50 subtype
percentDataPam50 <-
    clinical_data %>%
    select(CLAUDIN_SUBTYPE, INTCLUST) %>%
    mutate(CLAUDIN_SUBTYPE = case_when(CLAUDIN_SUBTYPE  == "Basal" ~ "Basal",
                                       CLAUDIN_SUBTYPE  == "claudin-low" ~ "Claudin-low",
                                       CLAUDIN_SUBTYPE  == "Her2" ~ "HER2",
                                       CLAUDIN_SUBTYPE  == "LumA" ~ "LumA",
                                       CLAUDIN_SUBTYPE  == "LumB" ~ "LumB",
                                       CLAUDIN_SUBTYPE  == "Normal" ~ "Normal", 
                                       CLAUDIN_SUBTYPE == "NC" ~ NA)) %>%
    mutate(INTCLUST = factor(INTCLUST, levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10"))) %>%
    group_by(INTCLUST) %>%
    count(CLAUDIN_SUBTYPE) %>%
    mutate(ratio = scales::percent(n / sum(n)))

ggplot(clinical_data %>%
           select(CLAUDIN_SUBTYPE, INTCLUST) %>%
           mutate(CLAUDIN_SUBTYPE = case_when(CLAUDIN_SUBTYPE  == "Basal" ~ "Basal",
                                              CLAUDIN_SUBTYPE  == "claudin-low" ~ "Claudin-low",
                                              CLAUDIN_SUBTYPE  == "Her2" ~ "HER2",
                                              CLAUDIN_SUBTYPE  == "LumA" ~ "LumA",
                                              CLAUDIN_SUBTYPE  == "LumB" ~ "LumB",
                                              CLAUDIN_SUBTYPE  == "Normal" ~ "Normal",
                                              CLAUDIN_SUBTYPE == "NC" ~ NA)) %>%
           mutate(INTCLUST = factor(INTCLUST, levels = c("1", "2", "3", "4ER-", "4ER+", "5", "6", "7", "8", "9", "10"))) ,
       aes(x = INTCLUST, fill = CLAUDIN_SUBTYPE)) +
    geom_bar(position = "fill") +
    geom_text(
        data = percentDataPam50,
        aes(y = n, label = ratio),
        position = position_fill(vjust = 0.5),
        size = 3
    ) +
    ggtitle("Stacked Barplot of Integrative Cluster Composition") +
    ylab("Percentage") +
    xlab("IntClust") +
    theme(
        plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(
            colour = "black",
            size = 10,
            face = "bold"
        ),
        legend.text = element_text(size = 10),
        legend.position = "top"
    ) +
    scale_fill_discrete(
        name = "Subtype:",
        labels = c("Basal", "Claudin-low", "HER2", "LumA", "LumB", "Normal", "NA")
    ) +
    scale_colour_discrete(
        name = "Subtype:",
        labels = c("Basal", "Claudin-low", "HER2", "LumA", "LumB", "Normal", "NA")
    ) +
    scale_y_continuous(labels = scales::percent) +
    guides(fill = guide_legend(ncol = 7))

ggsave(
    "../../figures/Introduction/IntClust_Composition.png",
    width = 10,
    height = 7
)
