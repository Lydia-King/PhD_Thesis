# Chapter 5: Simulate Data Functions

## Functions
Simulate_Data_NoNeut <-
    function(Num_Samples = 10,
             Percent_Neutral = c(90),
             Percent_Amp_Neutral = c(10),
             Percent_Neutral_Amp = c(0),
             Percent_Del_Neutral = c(0),
             Percent_Neutral_Del = c(0),
             Percent_Amp_Del_FP = c(0),
             Percent_Del_Amp_FP = c(0),
             Prob_CP_Occurs = c(1),
             Max_Num_CP = c(5),
             Prob_Amp_Neut_to_Neut_Amp = 0.5,
             Prob_Neut_Amp_to_Amp_Neut = 0.8,
             Prob_Del_Neut_to_Neut_Amp = 0.5,
             Prob_Neut_Del_to_Del_Neut = 0.5,
             Prob_Del_Amp_to_Amp_Del = 0.5,
             Prob_Amp_Del_to_Del_Neut = 0.5,
             Sample_Same = T,
             Chr_Start = 1,
             Chr_End = 250000000,
             Allele,
             Shuffle = F,
             SimData = 20,
             Lengths_Matrix,
             seed,
             samp_start = 0) {
        if (!missing(seed))
            set.seed(seed)
        
        # Create vector of numbers of samples in each category
        # Make sure number of individuals are whole numbers and add up to sample size
        times <-
            (
                Num_Samples * c(
                    Percent_Neutral,
                    Percent_Amp_Neutral,
                    Percent_Neutral_Amp,
                    Percent_Del_Neutral,
                    Percent_Neutral_Del,
                    Percent_Amp_Del_FP,
                    Percent_Del_Amp_FP
                )
            ) / 100 # Get number of individuals
        round_up <-
            c(sample(which(!is_wholenumber(times)), length(which(
                !is_wholenumber(times)
            )) / 2)) # elements to round up
        round_down <-
            which(!is_wholenumber(times))[which(!is_wholenumber(times)) %!in% round_up] #  elements to round down
        
        times[round_up] <- ceiling(times[round_up])
        times[round_down] <- floor(times[round_down])
        #  print(times)
        
        ## Create Matrix/df with Probabilities of transition between each state pair
        Prob_Matrix <-
            matrix(
                c(
                    Prob_Amp_Neut_to_Neut_Amp,
                    1 - Prob_Amp_Neut_to_Neut_Amp,
                    Prob_Neut_Amp_to_Amp_Neut,
                    1 - Prob_Neut_Amp_to_Amp_Neut,
                    Prob_Del_Neut_to_Neut_Amp,
                    1 - Prob_Del_Neut_to_Neut_Amp,
                    Prob_Neut_Del_to_Del_Neut,
                    1 -  Prob_Neut_Del_to_Del_Neut,
                    Prob_Del_Amp_to_Amp_Del,
                    1 - Prob_Del_Amp_to_Amp_Del,
                    Prob_Amp_Del_to_Del_Neut,
                    1 - Prob_Amp_Del_to_Del_Neut
                )
                ,
                nrow = 12,
                ncol = 1,
                byrow = T
            )
        
        rownames(Prob_Matrix) <-
            c(
                "ANNA",
                "ANND",
                "NAAN",
                "NAAD",
                "DNNA",
                "DNND",
                "NDDN",
                "NDDA",
                "DAAD",
                "DAAN",
                "ADDN",
                "ADDA"
            )
        colnames(Prob_Matrix) <- c("Probability")
        
        States_Vec <-
            c(
                "NA",
                "Amp/Neut",
                "Neut/Amp",
                "Del/Neut",
                "Neut/Del",
                "Amp/Del",
                "Del/Amp"
            )
        Start_Vec <-
            c("NA", "SAAN", "SNNA", "SDDN", "SNND", "SAAD", "SDDA")
        
        Store <- list()
        
        # Create Simulated Data
        for (NumSim in 1:SimData) {
            ## No Breakpoint
            if (times[1] == 0) {
                Output4 <-
                    data.frame(
                        "Sample" = character(),
                        "Chr" = numeric(),
                        "Type" = character(),
                        "Combine" = character(),
                        "CumLength" = numeric(),
                        "CP" = numeric(),
                        "TS" = numeric(),
                        "TE" = numeric(),
                        "Allele" = character()
                    )
            } else {
                if (length(Allele) == 1) {
                    Output4 <-
                        data.frame(
                            "Sample" = paste("Sample", c((1 + samp_start):(times[1] + samp_start)
                            )),
                            "Chr" = 1,
                            "Type" = "NoChangepoint",
                            "Combine" = NA,
                            "CumLength" = 0,
                            "CP" = NA,
                            "TS" = 0,
                            "TE" = 0,
                            "Allele" = Allele
                        )
                } else {
                    Output4 <-
                        data.frame(
                            "Sample" = rep(paste("Sample", 1:times[1]), each = 2),
                            "Chr" = rep(1, 2),
                            "Type" = rep("NoChangepoint", 2),
                            "Combine" = rep(NA, 2),
                            "CumLength" = rep(0, 2),
                            "CP" = rep(NA, 2),
                            "TS" = rep(0, 2),
                            "TE" = rep(0, 2),
                            "Allele" = c("Major", "Minor")
                        )
                }
            }
            
            DF <- Output4
            ## All Other Breakpoints
            Samp_Num <-
                times[1] + samp_start ## Keep Track of Sample Number
            ## Create DF to store outputs
            DF_CNA_CP <-
                data.frame(
                    "Sample" = character(),
                    "Chr" = numeric(),
                    "Type" = character(),
                    "Combine" = character(),
                    "CumLength" = numeric(),
                    "CP" = numeric(),
                    "TS" = numeric(),
                    "TE" = numeric(),
                    "Allele" = character()
                )
            
            Total_Mu_Sd <- Lengths_Matrix
            
            ## Create Sequence of CP for each patient (Same or Different)
            for (Cate_Type in c(which(times != 0)[!which(times != 0) %in% 1])) {
                ifelse(Sample_Same == T, patient1 <-
                           1, patient1 <- times[Cate_Type])
                
                for (patient in 1:patient1) {
                    Allele_DF <-
                        data.frame(
                            "Sample" = character(),
                            "Chr" = numeric(),
                            "Allele" = character(),
                            "Type" = character(),
                            "Combine" =  character(),
                            "CP" = numeric(),
                            "TS" = numeric(),
                            "TE" = numeric(),
                            "CumLength" = numeric()
                        )
                    
                    # ifelse(Sample_Same == T, Allele_Try <- 1,  Allele_Try <- length(Allele))
                    
                    for (allele in 1:length(Allele)) {
                        ## Number of mutations
                        Binom_CNA_CP <-
                            rbinom(Max_Num_CP[1] - 1, 1, Prob_CP_Occurs[1])
                        
                        ## Get Starting Category first
                        DF_CNA_CP_Temp <-
                            data.frame(
                                "Sample" = paste("Sample", Samp_Num + patient),
                                "Chr" = 1,
                                "Allele" = Allele[allele],
                                "Type_1" = States_Vec[Cate_Type],
                                "Combine_1" = Start_Vec[Cate_Type],
                                "CP_1" = NA ,
                                "TS_1" = NA,
                                "TE_1" = NA,
                                "CumLength_1" = NA
                            )
                        
                        ## Get rest of categories for each patient
                        State <- States_Vec[Cate_Type]
                        num_CNA_CP <- 1
                        
                        while (Binom_CNA_CP[num_CNA_CP] == 1 &&
                               num_CNA_CP <= length(Binom_CNA_CP)) {
                            ## Append Columns
                            if (State == "Neut/Amp") {
                                State <-
                                    sample(
                                        c("Amp/Neut", "Amp/Del"),
                                        1,
                                        prob = c(Prob_Matrix["NAAN", 1], Prob_Matrix["NAAD", 1])
                                    )
                                if (State == "Amp/Neut") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NAAN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NAAD",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Neut/Del") {
                                State <-
                                    sample(
                                        c("Del/Neut", "Del/Amp"),
                                        1,
                                        prob = c(Prob_Matrix["NDDN", 1], Prob_Matrix["NDDA", 1])
                                    )
                                if (State == "Del/Neut") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NDDN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NDDA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Amp/Neut") {
                                State <-
                                    sample(
                                        c("Neut/Amp", "Neut/Del"),
                                        1,
                                        prob = c(Prob_Matrix["ANNA", 1], Prob_Matrix["ANND", 1])
                                    )
                                if (State == "Neut/Amp") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ANNA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ANND",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Del/Neut") {
                                State <-
                                    sample(
                                        c("Neut/Amp", "Neut/Del"),
                                        1,
                                        prob = c(Prob_Matrix["DNNA", 1], Prob_Matrix["DNND", 1])
                                    )
                                if (State == "Neut/Amp") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DNNA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DNND",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Amp/Del") {
                                State <-
                                    sample(
                                        c("Del/Amp", "Del/Neut"),
                                        1,
                                        prob = c(Prob_Matrix["ADDA", 1], Prob_Matrix["ADDN", 1])
                                    )
                                if (State == "Del/Amp") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ADDA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ADDN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            }  else if (State == "Del/Amp") {
                                State <-
                                    sample(
                                        c("Amp/Del", "Amp/Neut"),
                                        1,
                                        prob = c(Prob_Matrix["DAAD", 1], Prob_Matrix["DAAN", 1])
                                    )
                                if (State == "Amp/Del") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DAAD",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DAAN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            }
                            
                            #  Temp <- as.data.frame(Temp)
                            colnames(Temp) <-
                                c(
                                    paste("Type_",  num_CNA_CP + 1, sep = ""),
                                    paste("Combine_",  num_CNA_CP + 1, sep = ""),
                                    paste("CP_",  num_CNA_CP + 1, sep = ""),
                                    paste("TS_",  num_CNA_CP + 1, sep = ""),
                                    paste("TE_",  num_CNA_CP + 1, sep = ""),
                                    paste("CumLength_",  num_CNA_CP + 1, sep = "")
                                )
                            DF_CNA_CP_Temp <-
                                cbind.data.frame(DF_CNA_CP_Temp, Temp)
                            num_CNA_CP <- num_CNA_CP + 1
                        }
                        
                        #     if(length(Allele) == 1 && Sample_Same == T){
                        if (Sample_Same == T) {
                            DF_CNA_CP_Temp <- DF_CNA_CP_Temp[rep(1,  times[Cate_Type]),]
                            
                            DF_CNA_CP_Temp[, 1] <-
                                c(paste("Sample", c(
                                    Samp_Num + 1:times[Cate_Type]
                                )))
                            DF_CNA_CP_Temp[, 3] <- Allele[allele]
                        }
                        
                        for (row in 1:nrow(DF_CNA_CP_Temp)) {
                            for (Combination_State in c(seq(5, ncol(DF_CNA_CP_Temp), 6))) {
                                Cate_State <- DF_CNA_CP_Temp[row, Combination_State]
                                
                                if (Combination_State == 5) {
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[Cate_State, 1],
                                                sd = Total_Mu_Sd[Cate_State, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row, Combination_State + 2] <-
                                        len
                                    DF_CNA_CP_Temp[row, Combination_State + 1]  <-
                                        1 + len
                                    
                                } else if (Combination_State == tail(c(seq(
                                    5, ncol(DF_CNA_CP_Temp), 6
                                )), n = 1)) {
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[Cate_State, 1],
                                                sd = Total_Mu_Sd[Cate_State, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row , Combination_State - 3] <-
                                        len
                                    DF_CNA_CP_Temp[row , Combination_State + 2] <-
                                        len
                                    
                                    DF_CNA_CP_Temp[row , Combination_State + 1]  <-
                                        DF_CNA_CP_Temp[row, Combination_State - 3] + DF_CNA_CP_Temp[row, Combination_State - 5]
                                    
                                    mock1 <-
                                        paste0(
                                            substr(Cate_State, 3, 4),
                                            substr(Cate_State, 4, 4),
                                            "E"
                                        )
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[mock1, 1],
                                                sd = Total_Mu_Sd[mock1, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row , Combination_State + 3] <-
                                        len
                                }  else {
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[Cate_State, 1],
                                                sd = Total_Mu_Sd[Cate_State, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row, Combination_State - 3] <-
                                        len
                                    DF_CNA_CP_Temp[row, Combination_State + 2] <-
                                        len
                                    
                                    DF_CNA_CP_Temp[row, Combination_State + 1]  <-
                                        DF_CNA_CP_Temp[row, Combination_State - 3] + DF_CNA_CP_Temp[row, Combination_State - 5]
                                }
                            }
                            
                            for (Combination_State in c(seq(5, ncol(DF_CNA_CP_Temp), 6))) {
                                Cate_State <- DF_CNA_CP_Temp[1, Combination_State]
                                
                                if (Combination_State == 5) {
                                    # Cumulative Length
                                    DF_CNA_CP_Temp[row, Combination_State + 4] <-
                                        DF_CNA_CP_Temp[row, Combination_State + 3] +  DF_CNA_CP_Temp[row, Combination_State + 2]
                                }  else {
                                    DF_CNA_CP_Temp[row, Combination_State + 4] <-
                                        DF_CNA_CP_Temp[row, Combination_State + 3] +  DF_CNA_CP_Temp[row, Combination_State - 2]
                                }
                            }
                        }
                        
                        DF_CNA_CP_Temp <-
                            DF_CNA_CP_Temp %>% pivot_longer(
                                cols = !c(Sample, Chr, Allele),
                                names_to = c(".value", "Samp"),
                                names_sep = "_",
                                values_drop_na = TRUE
                            ) %>% select(Sample,
                                         Chr,
                                         Type,
                                         Combine,
                                         CumLength,
                                         CP,
                                         TS,
                                         TE,
                                         Allele)
                        
                        
                        Data1 <-
                            data.frame(
                                "Sample" = character(),
                                "Chr" = numeric(),
                                "Type" = character(),
                                "Combine" = character(),
                                "CumLength" = numeric(),
                                "CP" = numeric(),
                                "TS" = numeric(),
                                "TE" = numeric(),
                                "Allele" = character()
                            )
                        for (pat in 1:length(unique(DF_CNA_CP_Temp$Sample))) {
                            Data <-
                                DF_CNA_CP_Temp %>% filter(Sample == unique(DF_CNA_CP_Temp$Sample)[pat])
                            if (length(which(Data$CumLength >= Chr_End)) > 0) {
                                Data <- Data[1:which(Data$CumLength >= Chr_End)[1], ]
                                Data[nrow(Data), "CumLength"] <- Chr_End
                                Data[nrow(Data), "TE"] <-
                                    Chr_End - Data[nrow(Data), "CP"]
                            } else {
                                Data[nrow(Data), "CumLength"] <- Chr_End
                                Data[nrow(Data), "TE"] <-
                                    Chr_End - Data[nrow(Data), "CP"]
                            }
                            Data1 <- rbind.data.frame(Data1, Data)
                        }
                        DF_CNA_CP_Temp <- Data1
                        # Allele_DF <- rbind.data.frame(Allele_DF, DF_CNA_CP_Temp)
                        
                        Allele_DF <-
                            rbind.data.frame(Allele_DF, DF_CNA_CP_Temp)
                    }
                    DF_CNA_CP <- rbind.data.frame(DF_CNA_CP,  Allele_DF)
                }
                
                Output3 <- rbind.data.frame(DF, DF_CNA_CP)
                
                Output3 <-
                    Output3 %>% mutate(TS = ifelse(Type %in% c(
                        "Neut/Amp", "Neut/Del"
                    ), 0, TS),
                    TE = ifelse(Type %in% c(
                        "Amp/Neut", "Del/Neut"
                    ), 0, TE))
                
                if (Shuffle == T) {
                    Output3 <- Output3 %>%
                        group_by(Allele) %>%
                        mutate(Sample = sample(Sample)) %>% arrange(Sample, Allele)
                }
                
                Output4 <- Output3
                
                Samp_Num <-
                    Samp_Num + times[Cate_Type] # Update Sample Number
            }
            
            Store[[NumSim]] <- Output4
        }
        
        return(Store)
    }

Simulate_Data_Neut <-
    function(Num_Samples = 10,
             Percent_Neutral = c(90),
             Percent_Amp_Neutral = c(10),
             Percent_Neutral_Amp = c(0),
             Percent_Del_Neutral = c(0),
             Percent_Neutral_Del = c(0),
             Percent_Amp_Del_FP = c(0),
             Percent_Del_Amp_FP = c(0),
             Prob_CP_Occurs = c(1),
             Max_Num_CP = c(5),
             Prob_Amp_Neut_to_Neut_Amp = 0.5,
             Prob_Neut_Amp_to_Amp_Neut = 0.8,
             Prob_Del_Neut_to_Neut_Amp = 0.5,
             Prob_Neut_Del_to_Del_Neut = 0.5,
             Prob_Del_Amp_to_Amp_Del = 0.5,
             Prob_Amp_Del_to_Del_Neut = 0.5,
             Sample_Same = T,
             Chr_Start = 1,
             Chr_End = 250000000,
             Allele,
             Shuffle = F,
             SimData = 20,
             Lengths_Matrix,
             seed,
             samp_start = 0) {
        if (!missing(seed))
            set.seed(seed)
        
        # Create vector of numbers of samples in each category
        # Make sure number of individuals are whole numbers and add up to sample size
        times <-
            (
                Num_Samples * c(
                    Percent_Neutral,
                    Percent_Amp_Neutral,
                    Percent_Neutral_Amp,
                    Percent_Del_Neutral,
                    Percent_Neutral_Del,
                    Percent_Amp_Del_FP,
                    Percent_Del_Amp_FP
                )
            ) / 100 # Get number of individuals
        round_up <-
            c(sample(which(!is_wholenumber(times)), length(which(
                !is_wholenumber(times)
            )) / 2)) # elements to round up
        round_down <-
            which(!is_wholenumber(times))[which(!is_wholenumber(times)) %!in% round_up] #  elements to round down
        
        times[round_up] <- ceiling(times[round_up])
        times[round_down] <- floor(times[round_down])
        #  print(times)
        
        ## Create Matrix/df with Probabilities of transition between each state pair
        Prob_Matrix <-
            matrix(
                c(
                    Prob_Amp_Neut_to_Neut_Amp,
                    1 - Prob_Amp_Neut_to_Neut_Amp,
                    Prob_Neut_Amp_to_Amp_Neut,
                    1 - Prob_Neut_Amp_to_Amp_Neut,
                    Prob_Del_Neut_to_Neut_Amp,
                    1 - Prob_Del_Neut_to_Neut_Amp,
                    Prob_Neut_Del_to_Del_Neut,
                    1 -  Prob_Neut_Del_to_Del_Neut,
                    Prob_Del_Amp_to_Amp_Del,
                    1 - Prob_Del_Amp_to_Amp_Del,
                    Prob_Amp_Del_to_Del_Neut,
                    1 - Prob_Amp_Del_to_Del_Neut
                )
                ,
                nrow = 12,
                ncol = 1,
                byrow = T
            )
        
        rownames(Prob_Matrix) <-
            c(
                "ANNA",
                "ANND",
                "NAAN",
                "NAAD",
                "DNNA",
                "DNND",
                "NDDN",
                "NDDA",
                "DAAD",
                "DAAN",
                "ADDN",
                "ADDA"
            )
        colnames(Prob_Matrix) <- c("Probability")
        
        States_Vec <-
            c(
                "NA",
                "Amp/Neut",
                "Neut/Amp",
                "Del/Neut",
                "Neut/Del",
                "Amp/Del",
                "Del/Amp"
            )
        Start_Vec <-
            c("NA", "SAAN", "SNNA", "SDDN", "SNND", "SAAD", "SDDA")
        
        Store <- list()
        
        # Create Simulated Data
        for (NumSim in 1:SimData) {
            ## No Breakpoint
            if (times[1] == 0) {
                Output4 <-
                    data.frame(
                        "Sample" = character(),
                        "Chr" = numeric(),
                        "Type" = character(),
                        "Combine" = character(),
                        "CumLength" = numeric(),
                        "CP" = numeric(),
                        "TS" = numeric(),
                        "TE" = numeric(),
                        "Allele" = character()
                    )
            } else {
                if (length(Allele) == 1) {
                    Output4 <-
                        data.frame(
                            "Sample" = paste("Sample", c((1 + samp_start):(times[1] + samp_start)
                            )),
                            "Chr" = 1,
                            "Type" = "NoChangepoint",
                            "Combine" = NA,
                            "CumLength" = 0,
                            "CP" = NA,
                            "TS" = 0,
                            "TE" = 0,
                            "Allele" = Allele
                        )
                } else {
                    Output4 <-
                        data.frame(
                            "Sample" = rep(paste("Sample", 1:times[1]), each = 2),
                            "Chr" = rep(1, 2),
                            "Type" = rep("NoChangepoint", 2),
                            "Combine" = rep(NA, 2),
                            "CumLength" = rep(0, 2),
                            "CP" = rep(NA, 2),
                            "TS" = rep(0, 2),
                            "TE" = rep(0, 2),
                            "Allele" = c("Major", "Minor")
                        )
                }
            }
            
            DF <- Output4
            ## All Other Breakpoints
            Samp_Num <-
                times[1] + samp_start ## Keep Track of Sample Number
            ## Create DF to store outputs
            DF_CNA_CP <-
                data.frame(
                    "Sample" = character(),
                    "Chr" = numeric(),
                    "Type" = character(),
                    "Combine" = character(),
                    "CumLength" = numeric(),
                    "CP" = numeric(),
                    "TS" = numeric(),
                    "TE" = numeric(),
                    "Allele" = character()
                )
            
            Total_Mu_Sd <- Lengths_Matrix
            
            ## Create Sequence of CP for each patient (Same or Different)
            for (Cate_Type in c(which(times != 0)[!which(times != 0) %in% 1])) {
                ifelse(Sample_Same == T, patient1 <-
                           1, patient1 <- times[Cate_Type])
                
                for (patient in 1:patient1) {
                    Allele_DF <-
                        data.frame(
                            "Sample" = character(),
                            "Chr" = numeric(),
                            "Allele" = character(),
                            "Type" = character(),
                            "Combine" =  character(),
                            "CP" = numeric(),
                            "TS" = numeric(),
                            "TE" = numeric(),
                            "CumLength" = numeric()
                        )
                    
                    # ifelse(Sample_Same == T, Allele_Try <- 1,  Allele_Try <- length(Allele))
                    
                    for (allele in 1:length(Allele)) {
                        ## Number of mutations
                        Binom_CNA_CP <-
                            rbinom(Max_Num_CP[1] - 1, 1, Prob_CP_Occurs[1])
                        
                        ## Get Starting Category first
                        DF_CNA_CP_Temp <-
                            data.frame(
                                "Sample" = paste("Sample", Samp_Num + patient),
                                "Chr" = 1,
                                "Allele" = Allele[allele],
                                "Type_1" = States_Vec[Cate_Type],
                                "Combine_1" = Start_Vec[Cate_Type],
                                "CP_1" = NA ,
                                "TS_1" = NA,
                                "TE_1" = NA,
                                "CumLength_1" = NA
                            )
                        
                        ## Get rest of categories for each patient
                        State <- States_Vec[Cate_Type]
                        num_CNA_CP <- 1
                        
                        while (Binom_CNA_CP[num_CNA_CP] == 1 &&
                               num_CNA_CP <= length(Binom_CNA_CP)) {
                            ## Append Columns
                            if (State == "Neut/Amp") {
                                State <-
                                    sample(
                                        c("Amp/Neut", "Amp/Del"),
                                        1,
                                        prob = c(Prob_Matrix["NAAN", 1], Prob_Matrix["NAAD", 1])
                                    )
                                if (State == "Amp/Neut") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NAAN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NAAD",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Neut/Del") {
                                State <-
                                    sample(
                                        c("Del/Neut", "Del/Amp"),
                                        1,
                                        prob = c(Prob_Matrix["NDDN", 1], Prob_Matrix["NDDA", 1])
                                    )
                                if (State == "Del/Neut") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NDDN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "NDDA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Amp/Neut") {
                                State <-
                                    sample(
                                        c("Neut/Amp", "Neut/Del"),
                                        1,
                                        prob = c(Prob_Matrix["ANNA", 1], Prob_Matrix["ANND", 1])
                                    )
                                if (State == "Neut/Amp") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ANNA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ANND",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Del/Neut") {
                                State <-
                                    sample(
                                        c("Neut/Amp", "Neut/Del"),
                                        1,
                                        prob = c(Prob_Matrix["DNNA", 1], Prob_Matrix["DNND", 1])
                                    )
                                if (State == "Neut/Amp") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DNNA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DNND",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            } else if (State == "Amp/Del") {
                                State <-
                                    sample(
                                        c("Del/Amp", "Del/Neut"),
                                        1,
                                        prob = c(Prob_Matrix["ADDA", 1], Prob_Matrix["ADDN", 1])
                                    )
                                if (State == "Del/Amp") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ADDA",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "ADDN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            }  else if (State == "Del/Amp") {
                                State <-
                                    sample(
                                        c("Amp/Del", "Amp/Neut"),
                                        1,
                                        prob = c(Prob_Matrix["DAAD", 1], Prob_Matrix["DAAN", 1])
                                    )
                                if (State == "Amp/Del") {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DAAD",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                } else {
                                    Temp <-
                                        data.frame(
                                            "Type" = State,
                                            "Combine" = "DAAN",
                                            "CP" = NA,
                                            "TS" = NA,
                                            "TE" = NA,
                                            "CumLength" = NA
                                        )
                                }
                            }
                            
                            #  Temp <- as.data.frame(Temp)
                            colnames(Temp) <-
                                c(
                                    paste("Type_",  num_CNA_CP + 1, sep = ""),
                                    paste("Combine_",  num_CNA_CP + 1, sep = ""),
                                    paste("CP_",  num_CNA_CP + 1, sep = ""),
                                    paste("TS_",  num_CNA_CP + 1, sep = ""),
                                    paste("TE_",  num_CNA_CP + 1, sep = ""),
                                    paste("CumLength_",  num_CNA_CP + 1, sep = "")
                                )
                            DF_CNA_CP_Temp <-
                                cbind.data.frame(DF_CNA_CP_Temp, Temp)
                            num_CNA_CP <- num_CNA_CP + 1
                        }
                        
                        #     if(length(Allele) == 1 && Sample_Same == T){
                        if (Sample_Same == T) {
                            DF_CNA_CP_Temp <- DF_CNA_CP_Temp[rep(1,  times[Cate_Type]),]
                            
                            DF_CNA_CP_Temp[, 1] <-
                                c(paste("Sample", c(
                                    Samp_Num + 1:times[Cate_Type]
                                )))
                            DF_CNA_CP_Temp[, 3] <- Allele[allele]
                        }
                        
                        for (row in 1:nrow(DF_CNA_CP_Temp)) {
                            for (Combination_State in c(seq(5, ncol(DF_CNA_CP_Temp), 6))) {
                                Cate_State <- DF_CNA_CP_Temp[row, Combination_State]
                                
                                if (Combination_State == 5) {
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[Cate_State, 1],
                                                sd = Total_Mu_Sd[Cate_State, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row, Combination_State + 2] <-
                                        len
                                    DF_CNA_CP_Temp[row, Combination_State + 1]  <-
                                        1 + len
                                    
                                } else if (Combination_State == tail(c(seq(
                                    5, ncol(DF_CNA_CP_Temp), 6
                                )), n = 1)) {
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[Cate_State, 1],
                                                sd = Total_Mu_Sd[Cate_State, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row , Combination_State - 3] <-
                                        len
                                    DF_CNA_CP_Temp[row , Combination_State + 2] <-
                                        len
                                    
                                    DF_CNA_CP_Temp[row , Combination_State + 1]  <-
                                        DF_CNA_CP_Temp[row, Combination_State - 3] + DF_CNA_CP_Temp[row, Combination_State - 5]
                                    
                                    mock1 <-
                                        paste0(
                                            substr(Cate_State, 3, 4),
                                            substr(Cate_State, 4, 4),
                                            "E"
                                        )
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[mock1, 1],
                                                sd = Total_Mu_Sd[mock1, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row , Combination_State + 3] <-
                                        len
                                }  else {
                                    len <-
                                        round(
                                            rtruncnorm(
                                                1,
                                                a = 1,
                                                b = Chr_End,
                                                mean = Total_Mu_Sd[Cate_State, 1],
                                                sd = Total_Mu_Sd[Cate_State, 2]
                                            )
                                        )
                                    DF_CNA_CP_Temp[row, Combination_State - 3] <-
                                        len
                                    DF_CNA_CP_Temp[row, Combination_State + 2] <-
                                        len
                                    
                                    DF_CNA_CP_Temp[row, Combination_State + 1]  <-
                                        DF_CNA_CP_Temp[row, Combination_State - 3] + DF_CNA_CP_Temp[row, Combination_State - 5]
                                }
                            }
                            
                            for (Combination_State in c(seq(5, ncol(DF_CNA_CP_Temp), 6))) {
                                Cate_State <- DF_CNA_CP_Temp[1, Combination_State]
                                
                                if (Combination_State == 5) {
                                    # Cumulative Length
                                    DF_CNA_CP_Temp[row, Combination_State + 4] <-
                                        DF_CNA_CP_Temp[row, Combination_State + 3] +  DF_CNA_CP_Temp[row, Combination_State + 2]
                                }  else {
                                    DF_CNA_CP_Temp[row, Combination_State + 4] <-
                                        DF_CNA_CP_Temp[row, Combination_State + 3] +  DF_CNA_CP_Temp[row, Combination_State - 2]
                                }
                            }
                        }
                        
                        DF_CNA_CP_Temp <-
                            DF_CNA_CP_Temp %>% pivot_longer(
                                cols = !c(Sample, Chr, Allele),
                                names_to = c(".value", "Samp"),
                                names_sep = "_",
                                values_drop_na = TRUE
                            ) %>% select(Sample,
                                         Chr,
                                         Type,
                                         Combine,
                                         CumLength,
                                         CP,
                                         TS,
                                         TE,
                                         Allele)
                        
                        
                        Data1 <-
                            data.frame(
                                "Sample" = character(),
                                "Chr" = numeric(),
                                "Type" = character(),
                                "Combine" = character(),
                                "CumLength" = numeric(),
                                "CP" = numeric(),
                                "TS" = numeric(),
                                "TE" = numeric(),
                                "Allele" = character()
                            )
                        for (pat in 1:length(unique(DF_CNA_CP_Temp$Sample))) {
                            Data <-
                                DF_CNA_CP_Temp %>% filter(Sample == unique(DF_CNA_CP_Temp$Sample)[pat])
                            if (length(which(Data$CumLength >= Chr_End)) > 0) {
                                Data <- Data[1:which(Data$CumLength >= Chr_End)[1], ]
                                Data[nrow(Data), "CumLength"] <- Chr_End
                                Data[nrow(Data), "TE"] <-
                                    Chr_End - Data[nrow(Data), "CP"]
                            } else {
                                Data[nrow(Data), "CumLength"] <- Chr_End
                                Data[nrow(Data), "TE"] <-
                                    Chr_End - Data[nrow(Data), "CP"]
                            }
                            Data1 <- rbind.data.frame(Data1, Data)
                        }
                        DF_CNA_CP_Temp <- Data1
                        # Allele_DF <- rbind.data.frame(Allele_DF, DF_CNA_CP_Temp)
                        
                        Allele_DF <-
                            rbind.data.frame(Allele_DF, DF_CNA_CP_Temp)
                    }
                    DF_CNA_CP <- rbind.data.frame(DF_CNA_CP,  Allele_DF)
                }
                
                Output3 <- rbind.data.frame(DF, DF_CNA_CP)
                
                if (Shuffle == T) {
                    Output3 <- Output3 %>%
                        group_by(Allele) %>%
                        mutate(Sample = sample(Sample)) %>% arrange(Sample, Allele)

                }
                
                Output4 <- Output3
                
                Samp_Num <-
                    Samp_Num + times[Cate_Type] # Update Sample Number
            }
            
            Store[[NumSim]] <- Output4
        }
        
        return(Store)
    }