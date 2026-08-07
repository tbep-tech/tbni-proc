#### Script Info ####
# Originally developed by: M Schrandt (April 3, 2020)
# Revisions by: M. Schram (August 7, 2026)

# Purpose: query FIM inshore database to essentially data-dump the catch data; no splitter
#          calculations or combining of similar gears. The idea is to have this data as close to "raw" as
#          possible and it is ONLY TO BE USED FOR CALCULATING THE TAMPA BAY NEKTON INDEX!

#### Set Up ####
library(tidyverse)
library(odbc)
library(dbplyr)
library(DBI)
library(here)

#### Connect to SQL database (corporate only) ####
# Note that server credentials have been removed from this script
Connect_FIMSQL <- function(db = "FIMCorpInshore"){
  dbConnect(
    odbc::odbc(),
    driver = "xxxx",
    server = "xxxx",
    database = db,
    uid = "xxxx",
    pwd = "xxxx",
    Encrypt = "Yes",
    TrustServerCertificate = "No"
  )
}


##### Function: Query FIM Data ####
AnnRep_getPhysBio_data <- 
  function(
    yr,
    bay
  ){
    # Name a connection for FIMCorpInshore (need to name it for use in next step(s))
    conn <- Connect_FIMSQL("FIMCorpInshore")
    
    # Physical data for all sets to get sampling effort and list of references
    Phys1 <- 
      tbl(
        conn,
        in_schema("hsdb", "tbl_corp_physical_master")
      ) %>%
      select(
        Reference, 
        Sampling_Date, 
        Project_1, 
        Project_2, 
        Project_3, 
        Gear, 
        Longitude, 
        Latitude, 
        Zone, 
        Grid, 
        Stratum
      ) %>%
      # Decompose Reference value into relevant data fields
      mutate(
        Bay   = substring(Reference, 1, 2),
        Type  = substring(Reference, 3,3),
        Year  = as.numeric(substring(Reference, 4,7)),
        Month = substring(Reference, 8,9),
        .after = Sampling_Date
      ) %>%
      collect() %>%
      filter(
        Bay %in% bay,
        Year %in% yr
      ) %>%
      mutate(
        # Remove white space from Zone column
        Zone = str_trim(
          Zone, 
          side = "both"
        ),
        # pad the Gear column with a zero so that it will sort properly later
        Gear = formatC(
          Gear, 
          width = 3, 
          flag = "0"
        )
      ) %>%
      arrange(
        Sampling_Date,
        Bay,
        Type,
        Reference
      )
    
    # List of relevant References for additional filtering 
    FullRefsList <- 
      Phys1 %>%
      filter(
        Type == "M",
        Gear %in% c("019", "020")
      ) %>%
      pull(Reference)
    
    # Pull in 'weather' data 
    Clouds <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_weather")
      ) %>%
      collect() %>%
      filter(
        # Earlier surveys recorded weather at the beginning and end of the survey.
        # This resulted in duplicate entries which led to merging issues. Since
        # there are at least an order of magnitude more "B" entries, felt
        # justified in dropping those "E" values for consistency. Likewise, this
        # only dropped a total of 4 values where the only record for that
        # Reference was "E", and 3 of those 4 were NA for CloudCover, with the
        # last being "0"
        Beg_end == "B",
      ) %>%
      select(
        Reference,
        CloudCover,
        Tide
      )
    
    # Pull in the effort column for later calculation of fish/100m2
    Gear_effort <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_gear")
      ) %>% 
      select(
        Reference, 
        Dist_tow, 
        Soaktime
      ) %>%
      collect() %>%
      mutate(
        Soaktime = hms::as_hms(Soaktime),
        Soaktime = lubridate::minute(lubridate::hms(Soaktime))
      )
    
    # Pull in gear details
    Gear_data <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_gear")
      ) %>%
      collect() %>%
      select(
        -Flag,
        -Dist_tow,
        -Soaktime
      )
    
    # Combine previous datasets
    Phys <- 
      Phys1 %>%
      left_join(
        Clouds,
        by = "Reference"
      ) %>%
      left_join(
        Gear_effort, 
        by = "Reference"
      ) %>%
      left_join(
        Gear_data,
        by = "Reference"
      ) %>%
      # remove white space from Zone column
      mutate(
        Zone = str_trim(Zone, side = "both")
      ) %>%
      # get effort expressed per 100m2
      mutate(
        effort = case_when(
          # 6.1-m seine
          Gear == "002" ~ 31.17/100,
          # 9.1-m seine
          Gear == "005" ~ 11/100,
          # 21-m offshore seine
          Gear %in% c("019", "020") ~ 140/100,
          # 21-m beach seine
          Gear == "022" ~ 368/100,
          # 21-m boat seine
          Gear == "023" ~ 68/100,
          # 61-m blocknet: No effort denoted; check procedure manual
          Gear == "153" ~ NA, 
          # 183-m haul seine
          Gear == "160" ~ 4120/100,
          # 183-m purse seine
          Gear == "170" ~ 2209/100,
          # 61-m haul seine
          Gear == "180" ~ 465/100,
          # Gillnets (in soaktime, not 100m2)
          Gear == "207" ~ Soaktime/60,
          # 6.1-m otter trawl
          Gear == "300" & !is.na(Dist_tow) ~ (Dist_tow*4*1853)/100, # 1853 m in nautical mile
          Gear == "300" & is.na(Dist_tow) & Soaktime > 1 ~ 0.2*Soaktime,
          # 1-m roving dropnet
          Gear == "350" ~ 1/100,
          .default = NA
        ),
        .after = Gear
      )
    
    rm(Clouds, Gear_effort, Gear_data, Phys1)
    
    # Biology (number of animals) data for all sets to get all the animals
    Bio2 <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_biology_number")
      ) %>%
      select(
        Reference, 
        Species_record_id, 
        Splittype, 
        Splitlevel,
        SC,
        Cells,
        NODCCODE, 
        Number, 
        FHC
      ) %>%
      mutate(
        Bay   = substring(Reference, 1, 2),
        Type  = substring(Reference, 3,3),
        Year  = as.numeric(substring(Reference, 4,7)),
        Month = substring(Reference, 8,9),
        .after = Reference
      ) %>%
      mutate(
        SC = toupper(SC)
      ) %>%
      filter(
        FHC != "D"
      ) %>%
      collect() %>%
      filter(
        Reference %in% FullRefsList
      ) %>%
      left_join(
        select(
          Phys,
          Reference,
          Sampling_Date,
        ),
        by = "Reference"
      ) %>%
      mutate(
        Count = case_when(
          !is.na(as.numeric(Splittype)) ~ Number*(as.numeric(Splittype)^as.numeric(Splitlevel)),
          .default = as.numeric(Number)
        )
      ) %>%
      mutate(
        N = sum(Count),
        .by = c(
          Reference,
          NODCCODE,
          Splittype,
          Splitlevel,
          SC,
        ),
      ) %>%
      select(
        -c(
          Number,
          Count,
          FHC,
          Cells
        )
      ) %>%
      relocate(
        Sampling_Date,
        .after = Reference
      ) %>%
      arrange(
        Sampling_Date,
        Bay,
        Type,
        Reference
      )
    
    # Add in the scientific names and taxa type
    Spp <- 
      tbl(
        conn,
        in_schema("hsdb", "tbl_corp_ref_species_list")
      ) %>%
      select(
        NODCCODE, 
        Scientificname, 
        Commonname
      ) %>%
      collect() %>%
      #assign type of taxa based on NODCCODE
      mutate(
        Taxa_Type = case_when(
          NODCCODE < '2000000000' ~ "Gear",
          NODCCODE > '2000000000' & NODCCODE < '7000000000' ~ "Invert",
          NODCCODE > '7000000000' & NODCCODE < '9000000000' ~ "Fish",
          NODCCODE > '9000000000' & NODCCODE < '9900000000' ~ "Turtle",
          NODCCODE > '9900000000' ~ 'Misc')
      )
    
    Species_List <- 
      tbl(
        conn,
        in_schema("hsdb", "tbl_corp_ref_species_list_selected")
      ) %>%
      collect() %>%
      right_join(
        Spp,
        by = "NODCCODE"
      ) %>%
      replace_na(
        list(
          Inshore_selected_taxa = FALSE,
          Offshore_selected_taxa = FALSE,
          FIM_selected_taxa = FALSE
        )
      ) %>%
      relocate(
        Scientificname:Taxa_Type,
        .after = NODCCODE
      ) %>%
      arrange(
        NODCCODE
      )
    
    Bio1 <- 
      Bio2 %>%
      left_join(
        Spp, 
        by = 'NODCCODE'
      ) %>%
      relocate(
        Taxa_Type,
        .after = NODCCODE
      ) %>%
      arrange(
        Reference, 
        Species_record_id
      )
    
    # Gather length data for the selected References
    Len1 <-
      tbl(
        conn,
        in_schema("hsdb", "tbl_corp_biology_lengths")
      ) %>%
      rename(
        SL = Length
      ) %>%
      select(
        -Flag
      ) %>%
      collect() %>%
      filter(
        Reference %in% FullRefsList,
        !is.na(SL)
      )
    
    # Merge catch and length information
    Len <-
      Len1 %>%
      left_join(
        Bio1,
        by = c(
          "Reference",
          "Species_record_id"
        )
      ) %>%
      relocate(
        c(Species_record_id, Length_record_id),
        .before = NODCCODE
      ) %>%
      relocate(
        SL,
        .after = Commonname
      ) %>%
      # Number of fish where same SL was recorded
      mutate(
        N_at_length = n(),
        .by = c(
          Reference,
          Scientificname,
          Splittype,
          Splitlevel,
          SC,
          SL
        )
      ) %>%
      # Look at the number measured in each species_record_id
      mutate(
        N_measured = n(),
        .by = c(
          Reference,
          Scientificname,
          Splittype,
          Splitlevel,
          SC
        )
      ) %>%
      # Total number of fish (inc. splitter est. from Bio data)
      mutate(
        N_total = sum(unique(N)),
        .by = c(
          Reference,
          Scientificname
        )
      ) %>%
      select(
        -c(
          Splittype,
          Splitlevel,
          Length_record_id
        )
      ) %>%
      # Run distinct so that individuals with same measurement at a length are on same line
      distinct(
        Reference,
        Scientificname,
        SL,
        N_measured,
        N_total,
        N_at_length,
        .keep_all = T
      ) %>%
      filter(
        !is.na(SL)
      ) %>%
      # This now provides the same information as the length data set that is
      # output from the sas program calculate ratio for fish per sz (SL), ceiling
      # will round up the weighting factor to nearest integer so don't have o.5
      mutate(
        # wf is now the adjusted count of fish per record & is column to build LF
        # distribution from
        wf   = floor(0.5 + (N_at_length/N_measured)*N),
      ) %>%
      arrange(
        Sampling_Date,
        Bay,
        Type,
        Reference,
        Scientificname
      )
    
    Bio <- 
      Bio1 %>%
      distinct(
        Reference, 
        NODCCODE,
        Splittype,
        Splitlevel,
        SC,
        .keep_all = TRUE
      ) %>%
      mutate(
        N_Total = sum(N),
        .by = c(Reference, NODCCODE)
      ) %>%
      distinct(
        Reference,
        NODCCODE,
        .keep_all = TRUE
      ) %>%
      select(
        -c(
          Splittype,
          Splitlevel,
          N,
          Taxa_Type
        )
      )
    
    # Remove the Bio1 data set to clean up the environment
    rm("Bio1", "Bio2", "Len1")
    
    Habitat <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_habitat")
      ) %>%
      collect() %>%
      filter(
        Reference %in% FullRefsList
      ) %>%
      mutate(
        Bay   = substring(Reference, 1, 2),
        Type  = substring(Reference, 3,3),
        Year  = as.numeric(substring(Reference, 4,7)),
        Month = substring(Reference, 8,9),
        .after = Reference
      ) %>%
      left_join(
        select(
          Phys,
          Reference,
          Sampling_Date
        ),
        by = "Reference"
      ) %>%
      relocate(
        Sampling_Date,
        .after = Reference
      ) %>%
      select(
        -Flag
      )
    
    Hydro <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_hydrolab")
      ) %>%
      collect() %>%
      filter(
        Reference %in% FullRefsList
      ) %>%
      mutate(
        Bay   = substring(Reference, 1, 2),
        Type  = substring(Reference, 3,3),
        Year  = as.numeric(substring(Reference, 4,7)),
        Month = substring(Reference, 8,9),
        .after = Reference
      ) %>%
      left_join(
        select(
          Phys,
          Reference,
          Sampling_Date
        ),
        by = "Reference"
      ) %>%
      relocate(
        Sampling_Date,
        .after = Reference
      ) %>%
      select(
        -Flag
      ) 
    
    RefCodes <- 
      tbl(
        conn, 
        in_schema("hsdb", "tbl_corp_ref_fim_codes")
      ) %>%
      collect() %>%
      filter(
        Bay %in% c("All", "TB")
      ) %>%
      mutate(
        Description = case_when(
          Bay == "TB" & FieldName == "Zone" ~ "See Tampa Bay FIM Sampling Map for Zone details",
          .default = Description
        )
      )
    
    # Disconnect from SQL database
    dbDisconnect(conn)
    rm("conn")  
    
    conn <- Connect_FIMSQL("FSAv1")
    
    dbDisconnect(conn)
    rm("conn")  
    
    PhysBioList <- list(Phys, Bio, Len, Habitat, Hydro, RefCodes, Species_List)
    
    return(PhysBioList)
  }

#### Query FIM SQL Database ####
beg_year <- 1998
end_year <- 2025
Bays <- c("TB")

dat <- 
  AnnRep_getPhysBio_data(
    yr   = beg_year:end_year,
    bay  = Bays
  ) 

#### Decompose FIM SQL Query output ####
FIM_PhysicalMaster <- dat[[1]]
FIM_BiologyCounts  <- dat[[2]]
FIM_BiologyLengths <- dat[[3]]
FIM_Habitat        <- dat[[4]]
FIM_HydroLab       <- dat[[5]]
FIM_ReferenceCodes <- dat[[6]]
FIM_SpeciesCodes   <- dat[[7]]

#### Export FIM data ####
saveRDS(
  FIM_PhysicalMaster,
  file = here("data/FIM_TB_PhysicalMaster.rds")
)

saveRDS(
  FIM_BiologyCounts,
  file = here("data/FIM_TB_BiologyCounts.rds")
)

saveRDS(
  FIM_BiologyLengths,
  file = here("data/FIM_TB_BiologyLengths.rds")
)

saveRDS(
  FIM_Habitat,
  file = here("data/FIM_TB_Habitat.rds")
)

saveRDS(
  FIM_HydroLab,
  file = here("data/FIM_TB_HydroLab.rds")
)

saveRDS(
  FIM_ReferenceCodes,
  file = here("data/FIM_TB_ReferenceCodes.rds")
)

saveRDS(
  FIM_SpeciesCodes,
  file = here("data/FIM_TB_SpeciesCodes.rds")
)
