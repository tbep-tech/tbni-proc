library(odbc)      #connect to FIM SQL database
library(tidyverse)
library(dbplyr)

#### Connect to SQL database (corporate only) ####
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

conn <- Connect_FIMSQL("FIMCorpInshore")

# Pull the official species list from the database
spp_list <- 
  collect(
    tbl(
      conn,
      in_schema('hsdb','tbl_corp_ref_species_list')
    )
  )

# Read the TBIndex species list
SppClass <- 
  readRDS(
    here("data/TBIndex_spp_codes.rds")
  )

# Join the two lists, change ScientificName and NODCCODE to match the official list, and write to a new csv
SppClass_New <- 
  left_join(
    SppClass, 
    spp_list, 
    by = 'TSN'
  ) %>%
  mutate(
    ScientificName = case_when(
      Scientificname != ScientificName ~ Scientificname,
      .default = ScientificName
    ),
    NODCCODE = case_when(
      NODCCODE.x != NODCCODE.y ~ NODCCODE.y,
      .default = NODCCODE.x
    )
  ) %>%
  select_at(
    cbind(
      'TSN',
      colnames(SppClass)
    )
  ) %>%
  saveRDS(
    here("data/TBIndex_spp_codes.rds")
  )

