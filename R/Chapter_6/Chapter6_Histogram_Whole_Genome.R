# Libraries 
library(tidyverse)
library(rCGH)
library(GeneBreak)

# Data
ASCAT_Data_3Step_NoNeut<- read.delim("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")
Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))

### Whole Genome (Barplots)
Total_ASCAT3_Chr3 <- ASCAT_Data_3Step_NoNeut 

## Combine Centromere info with Output
Centromere_Loc <- hg19 %>% mutate(chrom = ifelse(chrom == 23, "X", chrom))
Cent_Data <- merge(Centromere_Loc, Total_ASCAT3_Chr3, by.x = "chrom", by.y = "Chr")
Cent_Data <- Cent_Data %>% mutate(chrom = factor(chrom, levels = c(as.character(1:22), "X"))) %>% dplyr::rename(., Chromosome = chrom)
Cent_Data <- Cent_Data %>% mutate(Centre = (centromerStart+centromerEnd)/2) %>% mutate(Chromosome = paste("Chr", Chromosome, sep=" ")) %>% 
    mutate(Chromosome = factor(Chromosome, levels = c(paste("Chr", 1:22, sep=" "), "Chr X")))

# Plot 
Cent_Data %>% ggplot(aes(x=as.numeric(Changepoint))) + geom_histogram(binwidth = 1000000) + 
    facet_wrap(~Chromosome + Allele, scales = "free_x", ncol = 4) + geom_vline(data = Cent_Data, aes(xintercept = c(Centre)), color = "red") + 
    theme(legend.position = "top") + ylab("Frequency") + xlab("Genomic Location") + theme(strip.text.x = element_text(size = 10)) + 
    theme(axis.text.x = element_text(size = 7)) + ggtitle("Frequency of Changepoints across Whole Genome") +  theme(plot.title = element_text(hjust = 0.5)) + 
    theme(plot.title = element_text(size=18))

ggsave("../../figures/Chapter_6/Scaled_Frequency_Genomic_Seg.png", width = 10, height = 16, dpi = 400)

Cent_Data %>% ggplot(aes(x=as.numeric(Changepoint))) + geom_histogram(binwidth = 1000000) + 
    facet_wrap(~Chromosome + Allele, scales = "free", ncol = 4) + geom_vline(data = Cent_Data, aes(xintercept = c(Centre)), color = "red") + 
    theme(legend.position = "top") + ylab("Frequency") + xlab("Genomic Location") + theme(strip.text.x = element_text(size = 10)) + 
    theme(axis.text.x = element_text(size = 7)) +  ggtitle("Frequency of Changepoints across Whole Genome") +  theme(plot.title = element_text(hjust = 0.5)) +
    theme(plot.title = element_text(size=18))

ggsave("../../figures/Chapter_6/Unscaled_Frequency_Genomic_Seg.png", width = 10, height = 16)