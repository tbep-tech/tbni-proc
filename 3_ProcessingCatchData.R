#### Script Info ####
# M Schrandt
# April 3, 2020

# Purpose:
#  1. Prep the catch data for calculating the TB Nekton Index metrics and scores
#  2. Calculate all potential nekton index metrics (the reduced set of ~ 30)
#  3. Calculate reference conditions based on 1998-2015
#  4. Calculate Metric Scores for each reference/sample (only for metrics included in final TBNI)
#  5. Calculate TBNI score for each reference/sample (only for metrics included in final TBNI)
#  6. Graph the TBNI

#### Set Up ####
library(tidyverse)
library(lubridate)

# Specify catch file with TBEP segment assignments from script "2_Assign_TBSites_to_TBEPSegs"
CatchWithSegs <- "TampaBay_NektonIndexTBEPSegAssign_2020-04-19.csv"

# Specify the species guild classification file
class <- "TBIndex_spp_codes.csv"

#### 1. Prep Catch Data ####
input_data <- read.csv(CatchWithSegs, header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(NODCCODE = as.character(NODCCODE),
         NODCCODE = case_when(NODCCODE == "9.998e+09" ~ "9998000000",
                              NODCCODE == "9.999e+09" ~ "9999000000",
                              TRUE ~ NODCCODE),
         Sampling_Date = ymd(Sampling_Date),
         Bay = base::substr(.$Reference, 1, 2),
         Year = year(Sampling_Date),
         Month = month(Sampling_Date),
         #apply splitter calculation
         Count = case_when(!is.na(as.numeric(Splittype)) ~ Number*(as.numeric(Splittype)^as.numeric(Splitlevel)),
                           TRUE ~ as.numeric(Number))) %>%
  select(-Splittype, -Splitlevel,-Species_record_id, -Cells) %>%
  #add up the number of individuals by species within a reference (can be multiple records of a single species within 1 reference)
  group_by(Reference, NODCCODE) %>%
  mutate(Total_N = sum(Count)) %>%
  select(-Count, -Number) %>%
  #now just keep one record of the total number of indidivuals of each species within a record
  distinct() %>%
  #remove any possible instances of "no gear set", which is denoted with a particular NODCCODE
  filter(!NODCCODE == "9999000000")

# Read in the species classification and clean before merging with catch data
SppClass <- read.csv(class, header = TRUE, stringsAsFactors = FALSE) %>%
  mutate(NODCCODE = as.character(NODCCODE),
         #Need to fix an issue with how No fish is imported
         NODCCODE = case_when(NODCCODE == "9.998e+09" ~ "9998000000",
                              TRUE ~ NODCCODE))

# Merge catch data with species classifications
Counts <- input_data %>%
  left_join(SppClass, by = "NODCCODE") %>%
  #remove species we do not include in the index
  filter(Include_TB_Index == "Y") %>%
  arrange(Reference, NODCCODE)

#### 2. Calculate the TBNI Metrics ####
#There are 5 metrics in the TBNI small seine index
#1. Number of Taxa (NumTaxa)
#2. Shannon-Weiner diversity index (Shannon)
#3. Number of Selected Taxa (TaxaSelect)
#4. Number of Benthic Taxa (TaxaBenthic)
#5. Number of Feeding Guilds (NumGuilds)
# But we will calculate almost all 27 potential metrics (just leaving out a few sentinel spp)

# Calculate the number of taxa (NumTaxa) per set
NumTaxa <- Counts %>%
  select(Reference, NODCCODE) %>%
  group_by(Reference) %>%
  mutate(NumTaxa1 = n_distinct(NODCCODE)) %>%
  #make no fish caught count as zero taxa
  mutate(NumTaxa = case_when(NODCCODE == "9998000000" ~ 0,
                             TRUE ~ as.numeric(NumTaxa1))) %>%
  ungroup() %>%
  select(-NODCCODE, -NumTaxa1) %>%
  distinct() %>%
  arrange(Reference)

# Calculate the number of individuals per set
NumIndiv <- Counts %>%
  group_by(Reference) %>%
  summarize(NumIndiv = sum(Total_N)) %>%
  arrange(Reference)

# Calculate Diversity Indices - needs number of taxa and number of individuals
Diversity <- merge(NumTaxa, NumIndiv, by = "Reference") %>%
  merge(Counts) %>%
  filter(!NODCCODE == "9998000000" | !is.na(NumTaxa)) %>%
  group_by(Reference) %>%
  mutate(pi = Total_N/NumIndiv,
         pi2 = pi*pi,
         lnpi = log(pi),
         pilnpi = pi*lnpi) %>%
  mutate(sumpi2 = sum(pi2),
         sumpilnpi = sum(pilnpi),
         Shannon = -1*sumpilnpi,
         Simpson = 1/sumpi2,
         Pielou = Shannon/log(NumTaxa)) %>%
  select(Reference, Shannon, Simpson, Pielou) %>%
  distinct() %>%
  arrange(Reference)
#Will need to fix the No fish sets to have 0 taxa and zero for all metrics at the end - right now they are NaN's

# Calculate number of selected taxa per set
TaxaSelect <- Counts %>%
  select(Reference, NODCCODE, Selected_Taxa) %>%
  mutate(Select = case_when(Selected_Taxa == "Y" ~ 1,
                            Selected_Taxa == "N" ~ 0)) %>%
  #get distinct species list because sometimes multiple lines of 1 species
  distinct() %>%
  group_by(Reference) %>%
  summarize(TaxaSelect = sum(Select)) %>%
  arrange(Reference)
#Will need to fix the No fish sets to have 0 taxa

# Calculate the number of individuals of selected taxa
SelectIndiv <- Counts %>%
  select(Reference, Total_N, Selected_Taxa) %>%
  filter(Selected_Taxa == "Y") %>%
  group_by(Reference) %>%
  summarize(SelectIndiv = sum(Total_N))
  
# Calculate number of feeding guilds per set
NumGuilds <- Counts %>%
  select(Reference, NODCCODE, Feeding_Guild) %>%
  distinct() %>%
  group_by(Reference) %>%
  mutate(NumGuilds = n_distinct(Feeding_Guild)) %>%
  select(-NODCCODE, -Feeding_Guild) %>%
  distinct() %>%
  arrange(Reference)
#Will need to fix the No fish sets to have 0 taxa

# Calculate number of species needed to get to 90% of the catch (dominance)
Dom <- Counts %>%
  select(Reference, NODCCODE, Total_N) %>%
  arrange(Reference, desc(Total_N)) %>%
  left_join(NumIndiv, by = "Reference") %>%
  group_by(Reference) %>%
  mutate(CumProp = cumsum(Total_N)/NumIndiv,
         Taxa = row_number()) %>%
  ungroup %>%
  #retain rows that are >90%
  filter(CumProp >= 0.9) %>%
  # make sure only keep the first observation, closest to 90
  group_by(Reference) %>%
  slice(1) %>%
  ungroup %>%
  select(Reference, Taxa90 = Taxa)
  
# Calculate the number of taxa metrics
TaxaCountPrep <- Counts %>%
  group_by(Reference) %>%
  mutate(TSTaxa = case_when(Feeding_Cat == "TS" ~ 1, TRUE ~ 0),
         TGTaxa = case_when(Feeding_Cat == "TG" ~ 1, TRUE ~ 0),
         BenthicTaxa = case_when(Hab_Cat == "B" ~ 1, TRUE ~ 0),
         PelagicTaxa = case_when(Hab_Cat == "P" ~ 1, TRUE ~ 0),
         OblTaxa = case_when(Est_Use == "O" ~ 1, TRUE ~ 0),
         MSTaxa = case_when(Est_Cat == "MS" ~ 1, TRUE ~ 0),
         ESTaxa = case_when(Est_Cat == "ES" ~ 1, TRUE ~ 0))

taxa <- TaxaCountPrep %>%
  select(Reference, TSTaxa, TGTaxa, BenthicTaxa, PelagicTaxa,
         OblTaxa, MSTaxa , ESTaxa) %>%
  summarise_all(sum)

# Calculate number of individuals per functional ecological guild
abundances <- TaxaCountPrep %>%
  mutate(TSAbund = TSTaxa * Total_N,
         TGAbund = TGTaxa * Total_N,
         BenthicAbund = BenthicTaxa * Total_N,
         PelagicAbund = PelagicTaxa * Total_N,
         OblAbund = OblTaxa * Total_N,
         ESAbund = ESTaxa * Total_N,
         MSAbund = MSTaxa * Total_N) %>%
  select(Reference,TSAbund, TGAbund, BenthicAbund, PelagicAbund,
         OblAbund, ESAbund, MSAbund) %>%
  summarise_all(sum)
  
# Calculate abundances of sentinel species         
Num_LR <- Counts %>%
  ungroup() %>%
  select(Reference, Scientificname, Total_N) %>%
  group_by(Reference) %>%
  filter(Scientificname == "Lagodon rhomboides") %>%
  summarise(Num_LR = sum(Total_N))

# Merge metrics into 1 dataframe
# then calculate proportion metrics, and fix the "No fish" sets to have zero value for all metrics ####

Metrics <- Counts %>%
  ungroup() %>%
  #apply the season - based on avearge water temperatures of FIM sampling sites
  mutate(Season = case_when(Month %in% c(12, 1, 2, 3) ~ "Winter",
                            Month %in% c(4, 5) ~ "Spring",
                            Month %in% c(6, 7, 8, 9) ~ "Summer",
                            Month %in% c(10, 11) ~ "Fall")) %>%
  select(Reference, Year, Month, Season, TBEP_seg, Longitude, Latitude) %>%
  #only keep one record of each reference
  distinct() %>%
  left_join(NumTaxa, by = "Reference") %>%
  left_join(NumIndiv, by = "Reference") %>%
  left_join(Diversity, by = "Reference") %>%
  left_join(TaxaSelect, by = "Reference") %>%
  left_join(SelectIndiv, by = "Reference") %>%
  left_join(NumGuilds, by = "Reference") %>%
  left_join(Dom, by = "Reference") %>%
  left_join(taxa, by = "Reference") %>%
  left_join(abundances, by = "Reference") %>%
  left_join(Num_LR, by = "Reference") %>%
  #calculate all the proportion metrics
  mutate(PropTG = TGAbund/NumIndiv,
         PropTS = TSAbund/NumIndiv,
         PropBenthic = BenthicAbund/NumIndiv,
         PropPelagic = PelagicAbund/NumIndiv,
         PropObl = OblAbund/NumIndiv,
         PropMS = MSAbund/NumIndiv,
         PropES = ESAbund/NumIndiv,
         PropSelect = SelectIndiv/NumIndiv) %>%
  mutate(SelectIndiv = replace_na(SelectIndiv, 0),
         Num_LR = replace_na(Num_LR, 0),
         PropSelect = replace_na(PropSelect, 0)) %>%
  #set all metrics to zero if it's a no fish set
  mutate(NumTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ NumTaxa),
         Shannon = case_when(NumIndiv == 0 ~ 0, TRUE ~ Shannon),
         Simpson = case_when(NumIndiv == 0 ~ 0, TRUE ~ Simpson),
         Pielou = case_when(NumIndiv == 0 ~ 0, TRUE ~ Pielou),
         TaxaSelect = case_when(NumIndiv == 0 ~ 0, TRUE ~ TaxaSelect),
         SelectIndiv = case_when(NumIndiv == 0 ~ 0, TRUE ~ SelectIndiv),
         NumGuilds = case_when(NumIndiv == as.numeric(0) ~ 0, TRUE ~ as.numeric(NumGuilds)),
         Taxa90 = case_when(NumIndiv == as.numeric(0) ~ 0, TRUE ~ as.numeric(Taxa90)),
         TSTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ TSTaxa),
         TGTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ TGTaxa),
         BenthicTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ BenthicTaxa),
         PelagicTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ PelagicTaxa),
         OblTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ OblTaxa),
         MSTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ MSTaxa),
         ESTaxa = case_when(NumIndiv == 0 ~ 0, TRUE ~ ESTaxa),
         TSAbund = case_when(NumIndiv == 0 ~ 0, TRUE ~ TSAbund),
         TGAbund = case_when(NumIndiv == 0 ~ 0, TRUE ~ TGAbund),
         BenthicAbund = case_when(NumIndiv == 0 ~ 0, TRUE ~ BenthicAbund),
         PelagicAbund  = case_when(NumIndiv == 0 ~ 0, TRUE ~ PelagicAbund),
         OblAbund = case_when(NumIndiv == 0 ~ 0, TRUE ~ OblAbund),
         ESAbund = case_when(NumIndiv == 0 ~ 0, TRUE ~ ESAbund),
         MSAbund = case_when(NumIndiv == 0 ~ 0, TRUE ~ MSAbund),
         Num_LR = case_when(NumIndiv == 0 ~ 0, TRUE ~ Num_LR),
         PropTG = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropTG),
         PropTS = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropTS),
         PropBenthic = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropBenthic),
         PropPelagic = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropPelagic),
         PropObl = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropObl),
         PropMS = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropMS),
         PropES = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropES),
         PropSelect = case_when(NumIndiv == 0 ~ 0, TRUE ~ PropSelect))
         
# DO NOT RUN # Checking to see that my R calculations match the SAS metrics
# SASmets <- read.csv("SASmetrics.csv", header = TRUE, stringsAsFactors = FALSE) %>%
#   mutate_if(is.integer, as.numeric)
# library(arsenal)
# summary(comparedf(Metrics, SASmets, by = "Reference",
#                   #allow numbers to be off by 0.1 so output only shows those off by my than 0.1
#                   tol.num.val = 0.1))
# Looks like my SAS code is over-counting number of taxa and number of individuals (because apparently the SAS code
#  is still counting ones that we said to remove from the index). So I trust this R script more than the SAS code 
#  especially since I double-checked these NumTaxa and NumIndiv counts with the SQL database.


#### 3. Calculate Reference Conditions ####
# We need reference conditions for each bay segment and season
# We use the 5th and 95th percentiles for years 1998-2015 as the "Best available conditions"
# We only calculate reference conditions for the final variables in the TBNI
#  Number of Taxa, Number of BenthicTaxa, Number of Selected Taxa, Number of Feeding Guilds, S-W Diversity

# The 95th percentile is used as best avail. ref. condition for negative metrics: NumTaxa, BenthicTaxa, TaxaSelect, NumGuilds, Shannon
# All 5 of our metrics are considered negative metrics: we expect them to decline if there is an increase in perturbation/degradation
# Therefore, we use the 95th percentile as the best available condition beause that's the part associated with the best
#  environmental conditions

# But we need both percentiles in order to calculate the metric score.

RefCond <- Metrics %>%
  filter(between(Year, 1998, 2015)) %>%
  select(Season, TBEP_seg, NumTaxa, BenthicTaxa, TaxaSelect, NumGuilds, Shannon) %>%
  group_by(TBEP_seg, Season) %>%
  summarize(NumTaxa_P5 = round(quantile(NumTaxa, probs = 0.05)),
            NumTaxa_P95 = round(quantile(NumTaxa, probs = 0.95)),
            BenthicTaxa_P5 = round(quantile(BenthicTaxa, probs = 0.05)),
            BenthicTaxa_P95 = round(quantile(BenthicTaxa, probs = 0.95)),
            TaxaSelect_P5 = round(quantile(TaxaSelect, probs = 0.05)),
            TaxaSelect_P95 = round(quantile(TaxaSelect, probs = 0.95)),
            NumGuilds_P5 = round(quantile(NumGuilds, probs = 0.05)),
            NumGuilds_P95 = round(quantile(NumGuilds, probs = 0.95)),
            Shannon_P5 = quantile(Shannon, probs = 0.05),
            Shannon_P95 = quantile(Shannon, probs = 0.95))


#### 4. Calculate Metric Scores ####
Metric.Scores <- Metrics %>%
  select(Reference, Year, Month, Season, TBEP_seg, Longitude, Latitude, NumTaxa, BenthicTaxa, TaxaSelect, NumGuilds, Shannon) %>%
  left_join(RefCond, by = c("TBEP_seg", "Season")) %>%
  #calculate metric scores (rounded to nearest whole number)
  mutate(ScoreNumTaxa = case_when(NumTaxa == 0 ~ 0,
                                  NumTaxa > NumTaxa_P95 ~ 10,
                                  TRUE ~ ((NumTaxa-NumTaxa_P5)/(NumTaxa_P95-NumTaxa_P5))*10),
         ScoreShannon = case_when(NumTaxa == 0 ~ 0,
                                  Shannon > Shannon_P95 ~ 10,
                                  TRUE ~ ((Shannon-Shannon_P5)/(Shannon_P95-Shannon_P5))*10),
         ScoreTaxaSelect = case_when(NumTaxa == 0 ~ 0,
                                     TaxaSelect > TaxaSelect_P95 ~ 10,
                                     TRUE ~ ((TaxaSelect-TaxaSelect_P5)/(TaxaSelect_P95-TaxaSelect_P5))*10),
         ScoreTaxaBenthic = case_when(NumTaxa == 0 ~ 0,
                                      BenthicTaxa > BenthicTaxa_P95 ~ 10,
                                      TRUE ~ ((BenthicTaxa-BenthicTaxa_P5)/(BenthicTaxa_P95-BenthicTaxa_P5))*10),
         ScoreNumGuilds = case_when(NumTaxa == 0 ~ 0,
                                    NumGuilds > NumGuilds_P95 ~ 10,
                                    TRUE ~ ((NumGuilds-NumGuilds_P5)/(NumGuilds_P95-NumGuilds_P5))*10)) %>%
  #remove columns that are no-longer needed
  select(-BenthicTaxa_P95, -NumGuilds_P95, -NumTaxa_P95, -TaxaSelect_P95, -Shannon_P95,
         -BenthicTaxa_P5, -NumGuilds_P5, -NumTaxa_P5, -TaxaSelect_P5, -Shannon_P5,
         -NumTaxa, -BenthicTaxa, -TaxaSelect, -NumGuilds, -Shannon) %>%
  #clean up issue with some negative numbers for Benthic score; make them a score of zero
  mutate(ScoreTaxaBenthic = replace(ScoreTaxaBenthic, ScoreTaxaBenthic <0, 0)) %>%
  # now round all scores to nearest integer
  mutate(ScoreNumTaxa = round(ScoreNumTaxa),
         ScoreShannon = round(ScoreShannon),
         ScoreTaxaSelect = round(ScoreTaxaSelect),
         ScoreTaxaBenthic = round(ScoreTaxaBenthic),
         ScoreNumGuilds = round(ScoreNumGuilds))

#### 5. Calculate TBNI Scores ####
TBNI.Scores <- Metric.Scores %>%
  mutate(TBNI_Score = ((ScoreNumTaxa + ScoreShannon + ScoreTaxaSelect + ScoreTaxaBenthic + ScoreNumGuilds)/5)*10)

#save the scores to a file
write.csv(TBNI.Scores, paste0("TBNI_Scores_", min(TBNI.Scores$Year), "-", max(TBNI.Scores$Year), ".csv"), row.names = FALSE)

TBNI_Yrly_Scores <- TBNI.Scores %>%
  group_by(Year) %>%
  summarize(Yr_TBNI = round(mean(TBNI_Score),0),
            SEM = sd(TBNI_Score)/sqrt(n()),
            Lower.SEM = Yr_TBNI - SEM,
            Upper.SEM = Yr_TBNI + SEM)

TBNI_Seg_Scores <- TBNI.Scores %>%
  group_by(TBEP_seg, Year) %>%
  summarize(Segment_TBNI = round(mean(TBNI_Score),0))

#### 6. Graph the TBNI ####
# Calculate percentile for stoplight color communication
CommPerc <- TBNI.Scores %>%
  filter(between(Year, 1998, 2015)) %>%
  summarize(Score_P33 = quantile(TBNI_Score, probs = 0.33),
            Score_P50 = quantile(TBNI_Score, probs = 0.50))
CommPerc
# Apparently I used the 33 percentile in the report and manuscript but said it was the 25th
# So, we're going with actuall 0.33 percentile from here on our

# Assign our communication percentiles to our graph colors
# just to keep it straight in my head!
RedColor <- CommPerc$Score_P33
YellowColor <- CommPerc$Score_P50
GreenColor <- 100

Bay1 <- TBNI_Yrly_Scores %>%
  # Add max value for each color threshold as a column to refer to
  mutate(Red = RedColor,
         Yellow = YellowColor,
         Green = GreenColor)

TB_plot <- ggplot(Bay1) +
  geom_ribbon(aes(x = Year, ymin = 22, ymax= RedColor), fill = "red", alpha = 3/10) +
  geom_ribbon(aes(x = Year, ymin = RedColor, ymax = YellowColor), fill = "yellow", alpha = 3/10) +
  geom_ribbon(aes(x = Year, ymax = 58, ymin = YellowColor), fill = "green", alpha = 3/10) +
  geom_ribbon(aes(x = Year, ymin = Lower.SEM, ymax = Upper.SEM), fill = "gray88") +
  geom_line(aes(x = Year, y = Yr_TBNI), color = "black", size = 1.5) +
  geom_line(aes(x = Year, y = CommPerc$Score_P33), color = "black", linetype = "dotted") +
  geom_line(aes(x = Year, y = CommPerc$Score_P50), color = "black", linetype = "dotted") +
  scale_y_continuous(name = "TBNI score", limits = c(22, 58), breaks = seq(22, 58, 4), expand = c(0,0)) +
  scale_x_continuous(limits = c(1998, max(Bay1$Year)), breaks = seq(1998,max(Bay1$Year), 1), expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5)) +
  theme(axis.title.x = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
TB_plot

Seg_plot <- ggplot(TBNI_Seg_Scores) +
  geom_ribbon(aes(x = Year, ymin = 22, ymax= RedColor), fill = "red", alpha = 3/10) +
  geom_ribbon(aes(x = Year, ymin = RedColor, ymax = YellowColor), fill = "yellow", alpha = 3/10) +
  geom_ribbon(aes(x = Year, ymax = 58, ymin = YellowColor), fill = "green", alpha = 3/10) +
  geom_line(aes(x = Year, y = Segment_TBNI, linetype = TBEP_seg, color = TBEP_seg), size = 1.25) +
  scale_linetype_manual(name = "",
                        breaks = c("OTB", "HB", "MTB", "LTB"),
                        labels = c("OTB", "HB", "MTB", "LTB"),
                        #labels = c("Old Tampa Bay", "Hillsborough Bay",
                        #             "Middle Tampa Bay", "Lower Tampa Bay"),
                        values = c("dashed", "dotdash", "dotted", "solid")) +
  scale_color_manual(name = "",
                     breaks = c("OTB", "HB", "MTB", "LTB"),
                     labels = c("OTB", "HB", "MTB", "LTB"),
                     #labels = c("Old Tampa Bay", "Hillsborough Bay",
                     #             "Middle Tampa Bay", "Lower Tampa Bay"),
                     values = c("black", "black", "gray40", "gray40")) +
  scale_y_continuous(name = "TBNI score", limits = c(22, 58), breaks = seq(22, 58, 4), expand = c(0,0)) +
  scale_x_continuous(limits = c(1998, max(TBNI_Seg_Scores$Year)), breaks = seq(1998,max(TBNI_Seg_Scores$Year), 1), expand = c(0,0)) +
  theme(axis.text.y = element_text(size = 12)) +
  theme(axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5)) +
  geom_line(aes(x = Year, y = CommPerc$Score_P33), color = "black", linetype = "dotted") +
  geom_line(aes(x = Year, y = CommPerc$Score_P50), color = "black", linetype = "dotted") +
  theme(legend.title = element_text(size=12, face="bold")) +
  theme(legend.justification=c(1,0), legend.position=c(1,0),
        legend.key.size = unit(1, 'lines')) +
  #increase the length of the legend lines to be able to distinguish among lines
  theme(legend.key.width = unit(6, "line")) +
  theme(legend.key=element_blank()) +
  theme(legend.background=element_blank())+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
Seg_plot

#### STOP HERE ####
# This is for the manuscript
# Use cowplot package to arrange the graphs into a panel graph
# library(cowplot) 
#  TBNI_plot <- plot_grid(TB_plot, Seg_plot,
#                         nrow = 2,
#                         labels = c("A", "B"),
#                         align = "v",
#                         label_size = 14)
#  TBNI_plot
# 
# #Save the graph as 300 dpi for publication-quality
# #Assign a filename
#  name_as = paste0("TBNI_", max(Bay1$Year), "_", Sys.Date(), ".jpg")
#  save_plot(name_as, TBNI_plot, base_height = 7.2, base_width = 5)
 