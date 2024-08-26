#### Script Info ####
# M Schrandt
# April 3, 2020

# Purpose: query FIM inshore database to essentially data-dump the catch data; no splitter
#          calculations or combining of similar gears. The idea is to have this data as close to "raw" as
#          possible and it is ONLY TO BE USED FOR CALCULATING THE TAMPA BAY NEKTON INDEX!

#### Set Up ####
library(odbc)      #connect to FIM SQL database
library(tidyverse)
library(dbplyr)
library(lubridate) #dates/times
library(RCurl)     #to upload resulting csv to FIM's FTP site

#### Connect to SQL database (corporate only) ####
# Note that server credentials have been removed from this script
connex <- dbConnect(odbc::odbc(),
                    driver = "XXXX",
                    server = "XXXX",
                    database = "FIMCorpInshore",
                    uid = "XXXX",
                    pwd = "XXXX")

##### Collect physical data ####
Phys <- tbl(connex,in_schema("hsdb", "tbl_corp_physical_master")) %>%
  select(Reference, Sampling_Date, Project_1, Project_2, Project_3, Gear, Longitude, Latitude, Zone, Grid, Stratum) %>%
  filter(substring(Reference,1,3) == "TBM" &
        (Project_1 == "AM" | Project_2 == "AM" | Project_3 == "AM") &
         Gear %in% c("19", "20")) %>%
  collect() %>%
  filter(as.Date(Sampling_Date) >= "1998-01-01") %>%
  # add in effort expressed per 100m2
  mutate(effort = 140/100) %>%
  arrange(Reference)

#### Add catch data ####
Catch1 <- tbl(connex,in_schema("hsdb", "tbl_corp_biology_number")) %>%
  select(Reference, Species_record_id, Splittype, Splitlevel, Cells, 
         NODCCODE, Number, FHC) %>%
  filter(FHC != "D") %>%
  collect() %>%
  select(-FHC) %>%
  inner_join(Phys, by = "Reference")

# Add in the scientific names
Catch <- tbl(connex,in_schema("hsdb", "tbl_corp_ref_species_list")) %>%
  select(NODCCODE, Scientificname, Commonname) %>%
  collect() %>%
  inner_join(Catch1, by = 'NODCCODE') %>%
  select(Reference, Project_1, Project_2, Project_3, Sampling_Date, Latitude, Longitude, Zone, Grid, Stratum, effort,
         Species_record_id, NODCCODE, Scientificname, Commonname, Splittype, Splitlevel, Cells, Number) %>%
  arrange(Reference, Species_record_id)

#### Save Full Catch Table ####
write_csv(Catch, paste("TBNI_GitHub/tbni-proc/data/TampaBay_NektonIndexData.csv"))

#### Disconnect from the SQL database ####
dbDisconnect(connex)

#### Update metadata file
#Let's read in the metadata file from FIMaster and subset to the column names found in the nekton index data file
library(xlsx)
FIM_metadata <- read.xlsx("//fwc-spfs1/FIMaster/Data/Metadata/Inshore/metadata_sql_20200902.xlsx", sheetName = "data_fields") %>%
  select(Data.Field_name, Data_type, Field_length, Units.of.Measure, Precision, Description, Date_added, Date_discontinued) %>%
  filter(Data.Field_name %in% c("Reference", "Project_1", "Project_2", "Project_3",
                                "Sampling_Date", "Latitude", "Longitude", "Zone",
                                "Grid", "Stratum", "effort", "Species_record_id",
                                "NODCCODE", "Scientificname", "COmmonname", "Splittype",
                                "Splitlevel", "Cells")) %>%
  distinct() %>%
  add_row(Data.Field_name = "Commonname", Data_type = "Character", Description = "Common name for each species") %>%
  add_row(Data.Field_name = "effort", Data_type = "Decimal", Units.of.Measure = "per 100 m2",
          Description = "effort expressed per 100 m2") %>%
  add_row(Data.Field_name = "Number", Data_type = "Integer", Units.of.Measure = "individuals",
          Description = "Number of individuals collected")
#write to .csv
write_csv(FIM_metadata, paste("TBNI_GitHub/tbni-proc/data/TampaBay_NektonIndex_Metadata.csv"))



