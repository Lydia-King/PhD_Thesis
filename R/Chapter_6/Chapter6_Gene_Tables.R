# Chapter 6: Gene Table

## Load up Libraries 
library(tidyverse)
library(GeneBreak)
library(rCGH)
library(gt)

ASCAT_Data_3Step_NoNeut <- read.delim("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

colnames(ASCAT_Data_3Step_NoNeut) <- c("Sample", "Chr", "Changepoint", "TS", "TE", "Category", "Allele")
Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))

CNA_Loc <-
    read.delim(
        "../../data/Processed_Data/data_CNA_Loc_hg19.txt",
        sep = "\t",
        na.strings = c("", " ", "NA")
    )

CNA_Loc_3p <- CNA_Loc

CNA_Loc_3p <- CNA_Loc_3p  %>% mutate(Chromosome = ifelse(Chromosome == 23, "X", Chromosome))

Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut  

sample_df_2 <- data.frame("Gene" = character(0), "n" = numeric(0), "Allele" = character(0), "Category" = character(0), "Node" = character(0))

for (i in 1:nrow(CNA_Loc_3p)) {
    
    test <- Total_ASCAT3_Chr3 %>% mutate(Category = factor(Category), Allele = factor(Allele))
    
    test <- test %>% 
        filter(as.numeric(Total_ASCAT3_Chr3$Changepoint) >= CNA_Loc_3p[i,]$Start& as.numeric(Total_ASCAT3_Chr3$Changepoint) <= CNA_Loc_3p[i,]$End &
                   Total_ASCAT3_Chr3$Chr == CNA_Loc_3p[i,]$Chromosome) 
    
    df_to_fill <- data.frame("Category" = rep(levels(test$Category), 2), "Allele" = rep(levels(test$Allele), 7), "n" = 0, 
                             "Gene" = CNA_Loc_3p[i,]$Hugo_Symbol, Chromosome = CNA_Loc_3p[i,]$Chromosome)
    
    if(nrow(test) == 0){
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
        
    } else {
        test <- test %>% 
            group_by(Category, Allele, Chr) %>% summarise(n = n()) %>%
            mutate(Gene = CNA_Loc_3p[i,]$Hugo_Symbol,)
        
        for(j in 1:nrow(test)){
            df_to_fill$n[df_to_fill$Category == test$Category[j] & df_to_fill$Allele == test$Allele[j] & df_to_fill$Gene == test$Gene[j] & df_to_fill$Chr == test$Chr[j]] <- test$n[j]
        }
        
        sample_df_2 <- rbind.data.frame(sample_df_2, df_to_fill)
    } 
}

Final <- merge(sample_df_2, CNA_Loc[,c(1,4)], by.x = "Gene" , by.y = "Hugo_Symbol")
sample_df_2_Order <- Final %>% arrange(Chromosome, Start)
sample_df_2_Order$Gene <- as.character(sample_df_2_Order$Gene)

sample_df_2 <- write.table(sample_df_2_Order, "../../data/Chapter_6/Frequency_BP_Genes.txt", sep="\t")

p <- sample_df_2_Order  %>%
    group_by(Gene, Allele) %>%
    summarise(Category = 'Total', 
              n = sum(n), 
              Chromosome = Chromosome, 
              Start= Start) %>%
    bind_rows(sample_df_2_Order)

p <- p %>%   group_by(Gene, Category, Allele) %>% distinct()

test <- p %>% select(Gene, Chromosome, Allele, n) 

test1 <- pivot_wider(
    test,
    id_cols = c("Gene", "Chromosome", "Category"),
    names_from = Allele,
    values_from = n,
    names_prefix = ""
) %>% mutate(Total = Minor + Major) %>% dplyr::rename("Chr" = Chromosome)

as.data.frame(test1) %>% arrange(desc(Total)) %>% select(-Category) %>% dplyr::slice(1:10) %>% gt() |>
    tab_header(
        title =  md("**Frequency of Changepoints within Genes**"),
        subtitle = md("**Rows 1 to 10**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Gene_Table_Freq_1.png", path = "../../tables/Chapter_6/", vwidth = 700, vheight = 1000)

as.data.frame(test1) %>% arrange(desc(Total)) %>% select(-Category) %>% dplyr::slice(11:20) %>% gt() |>
    tab_header(
        title =  md("**Frequency of Changepoints within Genes**"),
        subtitle = md("**Rows 11 to 20**")
    ) |> tab_style(
        style = cell_text(weight = "bold"),
        locations = cells_column_labels()
    ) %>% opt_table_outline() |> tab_options(table.width = pct(70)) %>% 
    gtsave("Gene_Table_Freq_2.png", path = "../../tables/Chapter_6/", vwidth = 700, vheight = 1000)

Final <- merge(sample_df_2_Order %>% dplyr::select(-Start), CNA_Loc[,c(1,4)], by.x = "Gene" , by.y = "Hugo_Symbol")
sample_df_2_Order <- Final %>% arrange(Chromosome, Start)

sample_df_2_Order <- sample_df_2_Order %>% mutate(Chromosome = ifelse(Chromosome == 23, "X", Chromosome))

Test <- CNA_Loc %>% select(Hugo_Symbol, Chromosome, Start, band) %>% group_by(Chromosome) %>% arrange(Chromosome, Start)

Test$band <- gsub("(p|q).*", "\\1", Test$band)
Test1 <- Test %>% filter(band == "q") %>% group_by(Chromosome) %>% top_n(n=-1, wt=Start) %>% select(Chromosome, Hugo_Symbol)
Test1 <- Test1 %>% mutate(Chromosome = as.factor(Chromosome))

Sample_df_2_Order <- sample_df_2_Order %>% mutate(Chromosome = paste("Chr", Chromosome, sep= " ")) %>% mutate(Chromosome = factor(Chromosome, levels = c(paste("Chr", 1:22, sep= " "), "Chr X"))) 

Test1 <- Test1 %>% mutate(Chromosome = paste("Chr", Chromosome, sep= " ")) %>% mutate(Chromosome = factor(Chromosome, levels = c(paste("Chr", 1:23 , sep= " "), "Chr X"))) 

Test1[Test1$Chromosome == "Chr 23","Chromosome"] <- "Chr X"
Test1$Chromosome <- as.factor(as.character(Test1$Chromosome))
levels(Sample_df_2_Order$Chromosome)

ggplot(Sample_df_2_Order, aes(x=factor(Gene, levels = c(unique(sample_df_2_Order$Gene))), y=n)) + 
    geom_bar(stat = "identity")  +  facet_wrap(~Chromosome + Allele, scales = "free_x", ncol = 4) + 
    ggtitle("Frequency of Changepoints in Genes across Chromosomes") + 
    ylab("") + 
    xlab("") + theme(plot.title = element_text(hjust = 0.5)) + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
    theme(legend.position = "none") + geom_vline(data = Test1, aes(xintercept = c(Hugo_Symbol)), color = "red", size = 0.2) + 
    theme(legend.position = "top") + scale_fill_discrete(drop=FALSE) + scale_x_discrete(drop=TRUE)


ggsave("../../figures/Chapter_6/All_Genes_Barplot.png", last_plot(), width = 10,
       height = 16, units = c("in"), dpi = 400)
