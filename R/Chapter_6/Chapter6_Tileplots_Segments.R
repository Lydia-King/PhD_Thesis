# Load up Data
Model4_5mil_MCMCglmm <- read.delim("../../data/ASCAT_Data/Model4_Segments/Model4_5million_func1.txt", sep="\t")

# Libraries 
library(gt)
library(tidyverse)

## LM Tile Plots
### Segment function
ASCAT_Data_1 <- read.csv2("../../data/ASCAT_Data/3Step/ASCAT_Part_1.segments.txt", sep="\t")

ASCAT_Data_3Step_NoNeut <- read.csv2("../../data/Chapter_5/Total_ASCAT_Data_3Step_NoNeut.txt", sep="\t")

### 4) Segment function
## function 1 
Segment_function <- function(data, chr_start, chr_end, bin.size){
    seq_list <- list()
    seq_segment <- unique(c(seq(from = chr_start, to = chr_end, by = bin.size), chr_end))
    
    for(i in 1:(length(seq_segment)-1)){
        start <- seq_segment[i] 
        end <- seq_segment[i+1]
        
        data_in <- data %>% filter(c(as.numeric(Changepoint) > start & as.numeric(Changepoint) <= end) | Category %in% c("NoChangepoint_Amp", "NoChangepoint_Del"))
        data_out <- data %>% filter(as.numeric(Changepoint) < start | as.numeric(Changepoint) > end | Category == "NoChangepoint")
        
        data_in <- data_in %>% mutate(TS_Seg = TS, TE_Seg = TE, Category_Seg = Category, Changepoint_Seg = Changepoint) 
        data_out <- data_out %>% mutate(TS_Seg = 0, TE_Seg = 0, Category_Seg = "NoChangepoint", Changepoint_Seg = NA)
        
        data_out <- data_out[!duplicated(data_out[c("Sample", "Chr", "Allele")]),]
        
        int_com <- as.data.frame(intersect(data_in[,c("Sample", "Allele")], data_out[,c("Sample", "Allele")]))
        
        if(nrow(int_com) > 0){
            for(k in 1:nrow(int_com)){
                data_out <- data_out[-c(which(data_out$Sample == int_com$Sample[k] & data_out$Allele == int_com$Allele[k])), ]
            }
        }
        
        tot <- rbind.data.frame(data_in, data_out)     
        
        seq_list[[paste("chr_", unique(data$Chr), "_Segment_", i, sep="")]] <- tot
    }
    
    return(seq_list)
}

Chr_Info <- as.data.frame(ASCAT_Data_1 %>% group_by(chr) %>% summarise(start = min(startpos), end = max(endpos))) %>% 
    mutate(chr = factor(chr, levels = c(1:22, "X")))

Chr_Info <- Chr_Info[order(Chr_Info$chr),]

### 5) Segment Data (Into Different Size Chunks) - Neut 
ASCAT_Data_3Step_NoNeut <- ASCAT_Data_3Step_NoNeut %>% 
    mutate(TS = as.numeric(TS), TE = as.numeric(TE))

Univariate_lm_Model_2 <- function(Dataset, NoChangepoint_Rem = F){
    
    tryCatch(
        {
            
            if(isTRUE(NoChangepoint_Rem)){
                DF <- Dataset 
                DF <- within(DF, Category_Seg <- relevel(factor(Category_Seg), ref = "NoChangepoint"))
                DF <- DF %>% filter(Category_Seg != "NoChangepoint")
                mod_TS <- lm(cbind(TS_Seg) ~ 0 + Category_Seg*Allele, data =  DF)
                mod_TE <- lm(cbind(TE_Seg) ~ 0 + Category_Seg*Allele, data =  DF)
            } else {
                DF <- Dataset
                DF <- within(DF, Category_Seg <- relevel(factor(Category_Seg), ref = "NoChangepoint"))
                mod_TS <- lm(cbind(TS_Seg) ~ Category_Seg*Allele, data =  DF)
                mod_TE <- lm(cbind(TE_Seg) ~ Category_Seg*Allele, data =  DF)
            }
            
            SS <- DF %>% mutate(Category_Seg = as.factor(Category_Seg), Allele = as.factor(Allele)) %>% group_by(Allele) %>% count(Category_Seg, .drop = F)
            
            if(isTRUE(NoChangepoint_Rem)){
                SS <- SS %>% filter(Category_Seg != "NoChangepoint")
            } 
            
            
            sum_TS <- summary(mod_TS)
            sum_TE <- summary(mod_TE)
            
            pval_TS <- c()
            j <- 1
            for(i in 1:length(names(mod_TS$coefficients))){
                if(rownames(sum_TS$coefficients)[i] %in% c(names(mod_TS$coefficients))){
                    pval_TS <- c(pval_TS, sum_TS$coefficients[j,4])
                    j <- j + 1
                } else {
                    pval_TS <- c(pval_TS, NA)
                }
            }
            
            pval_TE <- c()
            j <- 1
            for(i in 1:length(names(mod_TE$coefficients))){
                if(rownames(sum_TE$coefficients)[i] %in% c(names(mod_TE$coefficients))){
                    pval_TE <- c(pval_TE, sum_TE$coefficients[j,4])
                    j <- j + 1
                } else {
                    pval_TE <- c(pval_TE, NA)
                }
            }
            
            cate <- names(mod_TS$coefficients)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub(":AlleleMinor", "", cate)
            cate <- gsub("AlleleMinor", names(mod_TS$coefficients)[1], cate)
            cate <- gsub("\\(Intercept\\)", "NoChangepoint", cate)
            cate <- gsub("Category_Seg", "", cate)
            pred_TS <- predict(mod_TS, newdata = data.frame("Category_Seg" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
            pred_TE <- predict(mod_TE, newdata = data.frame("Category_Seg" = cate, "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),  interval = "confidence")
            
            Model_DF1 <- rbind.data.frame(data.frame("Category" = cate, "Chr" = rep(1, length(mod_TS$coefficients)),
                                                     "Direction" = rep("TS", length(mod_TS$coefficients)), "n" = SS$n,
                                                     "Beta" = unname(mod_TS$coefficients), "P" = pval_TS,
                                                     "LB" = pred_TS[,2], "UB" = pred_TS[,3], "Pred" = pred_TS[,1],
                                                     "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))),
                                          data.frame("Category" = cate, "Chr" = rep(1, length(mod_TE$coefficients)), "Direction" = rep("TE", length(mod_TE$coefficients)),
                                                     "n" = SS$n,  "Beta" = unname(mod_TE$coefficients), "P" = pval_TE,
                                                     "LB" = pred_TE[,2], "UB" = pred_TE[,3], "Pred" = pred_TE[,1],
                                                     "Allele" = c(rep("Major", length(which(SS$Allele == "Major"))), rep("Minor", length(which(SS$Allele == "Minor"))))))
            
            Data_Chr1 <- Model_DF1  %>% mutate(Category = gsub("Category_Seg", "", Category)) %>% mutate(Category = gsub("Category", "", Category))
            Data_Chr1 <- Data_Chr1 %>% filter(n > 0)
            Data_Chr1$Category <- gsub("\\(Intercept\\)", "NoChangepoint", Data_Chr1$Category)
            Data_Chr1 <- Data_Chr1 %>% mutate(Category = factor(Category, levels = c("NoChangepoint", "Neut/Amp", "Neut/Del", "Amp/Neut", "Del/Neut", "Amp/Del", "Del/Amp")))
            
            return(Data_Chr1)
            
        },
        error=function(e) {
            Pred_DF1 <- data.frame("Category" = NA,
                                   "Chr" = 1,
                                   "Direction" = NA,
                                   "n" = NA,
                                   "Beta" = NA, 
                                   "P" = NA,
                                   "LB" = NA, 
                                   "UB" = NA, 
                                   "Pred" = NA,
                                   "Allele" = NA)
            
            Pred_DF1
            
        }
    )
    
} 

Test_Results_Multi <- function(dataset){
    List1 <- list()
    
    for(chr in 1:length(dataset)){
        Dataset <- dataset[[chr]]
        List2 <- list()
        
        for(seg in 1:length(Dataset)){
            Dataset1 <- Dataset[[seg]]
            print(paste("Current chr:", chr, sep = " "))
            print(paste("Current seg:", seg, sep = " "))
            List2[[paste("Seg_", seg, sep="")]] <- Univariate_lm_Model_2(Dataset1)
        }
        
        List1[[paste("Chr_", chr, sep="")]] <- List2
    }
    
    return(List1)
}

Segment_Function <- function(segment_width = 1000000){
    
    WG_List_5million_func1 <- list()
    
    chr <- 1
    for(i in c(1:22, "X")){
        Chr_K <- ASCAT_Data_3Step_NoNeut %>% filter(Chr == i)
        WG_List_5million_func1[[paste("Chr", i)]] <- Segment_function(Chr_K, Chr_Info[chr, 2], Chr_Info[chr, 3], segment_width)
        chr <- chr + 1
    }
    
    Test_Results_5million_Multi <- Test_Results_Multi(WG_List_5million_func1)
    
    df1 <- purrr::map(Test_Results_5million_Multi, .f = ~bind_rows(., .id = "Segment"))
    df2 <- bind_rows(df1, .id = "Chr")
    
    return(df2)
}

Model4_5mil_LM <- Segment_Function(segment_width = 5000000)

Name <- c("Five_Million")

TilePlot_Func <- function(model = c("LM", "MCMC"), thres = 10, Data, Name, width){
    
    
    Dataframe_Empty <- data.frame("Chr" = character(), "Segment" = character(),  "Category" = character(),  "n" = numeric(), 
                                  "Direction" = character(),  "Allele" = character(),    "Fit" = numeric(), "LB" = numeric(),        
                                  "UB" = numeric())
    
    for(i in 1:23){
        Data_1 <- Data %>% filter(Chr == paste("Chr_", i, sep="")) %>% 
            complete(Segment, Chr, Category, Allele) %>% filter(!is.na(Allele)) %>% 
            filter(!is.na(Category))
        
        Dataframe_Empty <- rbind.data.frame(Dataframe_Empty, Data_1)
    }
    
    grouped_data <- Dataframe_Empty %>%  dplyr::mutate(Category = dplyr::case_when(Category == "NoChangepoint" ~ "NoCP",
                                                                                   Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                                                   Category == "Neutral/Deletion" ~ "Neut/Del",
                                                                                   Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                                                   Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                                                   Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                                                   Category == "Amplification/Deletion" ~ "Amp/Del", 
                                                                                   .default = Category)) %>%
        mutate(Category = factor(Category,  levels = rev(c("NoCP", "Amp/Neut", "Del/Neut", 
                                                           "Neut/Amp", "Neut/Del", 
                                                           "Amp/Del", "Del/Amp")))) %>%  
        mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
        mutate(Chr = factor(Chr, levels = c(paste0("Chr_", c(1:22, "X"))))) %>%
        group_by(Chr, Segment, Category, Allele) %>%
        mutate(Sig = case_when(
            is.na(Direction) & is.na(LB) ~ NA,
            any(Direction == "TS" & LB > thres) & any(Direction == "TE" & LB > thres) ~ "Sig",
            any(Direction == "TS" & LB > thres) & any(Direction == "TE" & LB <= thres) ~ "Sig_TS",
            any(Direction == "TS" & LB <= thres) & any(Direction == "TE" & LB > thres) ~ "Sig_TE",
            TRUE ~ "NSig"
        )) %>%
        group_by(Chr, Segment, Category, Allele, Sig) %>%
        summarise(
            n_new = sum(n),
            n_Sig = sum(ifelse(Sig == "Sig", n, 0)),
            n_Sig_TS = ifelse(Sig == "Sig_TS", sum(n[Direction == "TS"]), 0),
            n_Sig_TE = ifelse(Sig == "Sig_TE", sum(n[Direction == "TE"]), 0)) %>% mutate("Final_n" = case_when(Sig == "Sig" ~ n_new,
                                                                                                               Sig == "Sig_TS" ~ n_Sig_TS,
                                                                                                               Sig == "Sig_TE" ~ n_Sig_TE,
                                                                                                               Sig == "NSig" ~ NA))
    
    grouped_data$Segment<- factor(grouped_data$Segment, levels = c(paste("Seg_", 1:300, sep= "")))
    grouped_data <- grouped_data %>% distinct()
    
    ggplot(grouped_data,
           aes(x = Segment, y = Category, fill = Sig)) + geom_tile()  + scale_fill_manual(
               values = c("white", "white", "white", "white"),
               na.value = "gray91",
               drop = F,
               guide = "none") + geom_point(aes(size = Final_n, colour = Sig),
                                            shape = 16) + 
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()) +  theme(axis.text.x = element_blank()) + theme(axis.ticks.x = element_blank())  + ggtitle(paste("Changepoints of Significant Length across", width, "Segments (LB > ", thres, "kb)", sep="")) + 
        theme(plot.title = element_text(hjust = 0.5)) +  
        scale_color_discrete(na.translate=FALSE) +
        scale_shape_manual(
            name = "Sig & Sample Size",
            labels = c("250", "500", "750"),
            values = c(16, 16, 16)
        ) +  scale_size_continuous(name = "Sample\nSize") +  facet_wrap(Chr~Allele, scales = "free_x", ncol = 6) + xlab("") + ylab("") +  theme(plot.title = element_text(size=18))
    
    ggsave(paste0("PerSegments_", model, "_", thres, "kb_", Name ,"_Thesis.png"), 
           plot = last_plot(),
           path = "../../figures/Chapter_6/",
           width = 13,
           height = 15,
           units = c("in"),
           dpi = 300
    )
}

# Data cleaning 
Model4_5mil_LM <- Model4_5mil_LM[!is.na(Model4_5mil_LM$Category),]

TilePlot_Func(model = c("LM"), thres = 10, Model4_5mil_LM %>% filter(n!= 0), Name = "fivemillion", width = " 5 million base ")
TilePlot_Func(model = c("LM"), thres = 10000, Model4_5mil_LM %>% filter(n!= 0), Name = "fivemillion", width = " 5 million base ")

## MCMC
Chr_Seg_List_5 <- unique(paste0(Model4_5mil_LM$Chr, "_", Model4_5mil_LM$Segment))
Model4_5mil_MCMCglmm$Ext <- paste0(Model4_5mil_MCMCglmm$Chr, "_",  Model4_5mil_MCMCglmm$Segment)
Model4_5mil_MCMCglmm <-  Model4_5mil_MCMCglmm[ Model4_5mil_MCMCglmm$Ext %in% Chr_Seg_List_5,]

## MCMC Tile Plots 
TilePlot_Func(model = c("MCMC"), thres = 10, Model4_5mil_MCMCglmm, Name = "fivemillion", width = " 5 million base ")
TilePlot_Func(model = c("MCMC"), thres = 10000, Model4_5mil_MCMCglmm, Name = "fivemillion", width = " 5 million base ")



# MCMC and LM Tables
Table_Func <- function(thres = 10000, Data_LM, Data_MCMC, Name){
    Tab_1 <- Data_LM %>% select(Chr, Segment, Allele, Category, n, Direction, Pred, LB, UB) %>% 
        dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                  Category == "Neutral/Deletion" ~ "Neut/Del",
                                                  Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                  Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                  Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                  Category == "Amplification/Deletion" ~ "Amp/Del", 
                                                  Category == "NoChangepoint" ~ "NoCP")) %>% 
        dplyr::rename(Fit = Pred) %>%
        filter(LB > thres, n >200) %>% arrange(desc(Fit)) %>% 
        gt()  %>%
        tab_header(
            title = md("**(A) Segments Containing Changepoints with Large CNAs**"), 
        )
    
    
    Tab_2 <- Data_MCMC %>% select(Chr, Segment, Allele, Category, n, Direction, Fit, LB, UB) %>% 
        dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                                  Category == "Neutral/Deletion" ~ "Neut/Del",
                                                  Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                                  Category == "Deletion/Neutral" ~ "Del/Neut", 
                                                  Category == "Deletion/Amplification" ~ "Del/Amp", 
                                                  Category == "Amplification/Deletion" ~ "Amp/Del", 
                                                  Category == "NoChangepoint" ~ "NoCP")) %>% 
        filter(LB > thres, n > 200) %>% arrange(desc(Fit)) %>% 
        gt()  %>%
        tab_header(
            title = md("**(B) Segments Containing Changepoints with Large CNAs**"), 
        )
    
    Tab_1 |> gtsave(filename = paste("../../tables/Chapter_6/Segment_LM_", thres, "_", Name, "_Thesis.png", sep=""), vwidth = 3000, vheight = 2500)
    Tab_2 |> gtsave(filename = paste("../../tables/Chapter_6/Segment_MCMC_", thres,  "_", Name, "_Thesis.png", sep=""), vwidth = 3000, vheight = 2500)
    
}

Table_Func(thres = 10000, Model4_5mil_LM, Model4_5mil_MCMCglmm, Name = "fivemillion")

### Tables 
## N >200, 10kb

Tab_1 <- Model4_5mil_LM %>% 
    mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
    select(Chr, Segment, Allele, Category, n, Direction, Pred, LB, UB) %>% 
    dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                              Category == "Neutral/Deletion" ~ "Neut/Del",
                                              Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                              Category == "Deletion/Neutral" ~ "Del/Neut", 
                                              Category == "Deletion/Amplification" ~ "Del/Amp", 
                                              Category == "Amplification/Deletion" ~ "Amp/Del", 
                                              Category == "NoChangepoint" ~ "NoCP")) %>% 
    dplyr::rename(Fit = Pred) %>%
    filter(LB > 10, n >200) %>% arrange(desc(Fit)) %>% 
    gt()  %>%
    tab_header(
        title = md("**(A) Segments Containing n > 200 Changepoints with CNAs > 10kb**"), 
    )


Tab_2 <- Model4_5mil_MCMCglmm %>% 
    mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
    select(Chr, Segment, Allele, Category, n, Direction, Fit, LB, UB) %>% 
    dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                              Category == "Neutral/Deletion" ~ "Neut/Del",
                                              Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                              Category == "Deletion/Neutral" ~ "Del/Neut", 
                                              Category == "Deletion/Amplification" ~ "Del/Amp", 
                                              Category == "Amplification/Deletion" ~ "Amp/Del", 
                                              Category == "NoChangepoint" ~ "NoCP")) %>% 
    filter(LB > 10, n > 200) %>% arrange(desc(Fit)) %>% 
    gt()  %>%
    tab_header(
        title = md("**(B) Segments Containing n > 200 Changepoints with CNAs > 10kb**"), 
    )

Tab_1 |> gtsave(filename = paste("../../tables/Chapter_6/Segment_LM_", 10, "_", "Five", "_Thesis.png", sep=""), vwidth = 3000, vheight = 2500)
Tab_2 |> gtsave(filename = paste("../../tables/Chapter_6/Segment_MCMC_", 10,  "_", "Five", "_Thesis.png", sep=""), vwidth = 3000, vheight = 2500)

## n > 200, 10000kb 

Tab_1 <- Model4_5mil_LM %>% select(Chr, Segment, Allele, Category, n, Direction, Pred, LB, UB) %>% 
    dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                              Category == "Neutral/Deletion" ~ "Neut/Del",
                                              Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                              Category == "Deletion/Neutral" ~ "Del/Neut", 
                                              Category == "Deletion/Amplification" ~ "Del/Amp", 
                                              Category == "Amplification/Deletion" ~ "Amp/Del", 
                                              Category == "NoChangepoint" ~ "NoCP")) %>% 
    dplyr::rename(Fit = Pred) %>%
    filter(LB > 10000, n >200) %>% arrange(desc(Fit)) %>% 
    gt()  %>%
    tab_header(
        title = md("**(A) Segments Containing n > 200 Changepoints with CNAs > 10,000kb**"), 
    )


Tab_2 <- Model4_5mil_MCMCglmm %>% 
    mutate(Chr = ifelse(Chr == "Chr_23", "Chr_X", Chr)) %>%
    select(Chr, Segment, Allele, Category, n, Direction, Fit, LB, UB) %>% 
    dplyr::mutate(Category = dplyr::case_when(Category == "Neutral/Amplification" ~ "Neut/Amp",
                                              Category == "Neutral/Deletion" ~ "Neut/Del",
                                              Category == "Amplification/Neutral" ~ "Amp/Neut", 
                                              Category == "Deletion/Neutral" ~ "Del/Neut", 
                                              Category == "Deletion/Amplification" ~ "Del/Amp", 
                                              Category == "Amplification/Deletion" ~ "Amp/Del", 
                                              Category == "NoChangepoint" ~ "NoCP")) %>% 
    filter(LB > 10000, n > 200) %>% arrange(desc(Fit))  %>%
    gt()  %>%
    tab_header(
        title = md("**Segments Containing n > 200 Changepoints with CNAs > 10,000kb**"), 
    )

Tab_1 |> gtsave(filename = paste("../../tables/Chapter_6/Segment_LM_", 10000, "_", "Five", "_Thesis.png", sep=""), vwidth = 3000, vheight = 2500)
Tab_2 |> gtsave(filename = paste("../../tables/Chapter_6/Segment_MCMC_", 10000,  "_", "Five", "_Thesis.png", sep=""), vwidth = 3000, vheight = 2500)
