# README

Original processing scripts for Tampa Bay Nekton Index[^1], see [tbeptools](https://tbep-tech.github.io/tbeptools) for working functions.

The data object [fimdata](https://tbep-tech.github.io/tbeptools/reference/fimdata.html)[^2] in tbeptools is created from `data/TampaBay_NektonIndexData.csv`, which is created in `R/1_Query_FIM_Database.R`. This script can only be used by FIM staff with appropriate database credentials and is included here for transparency. 

`1_Query_FIM_Database.R` queries the FIM database and writes numerous .rds files as well as the `TampaBay_NektonIndexData.csv` file. These data files are subsequently transferred to TBEP for use in the Nekton Index and as data deliverables for other TBEP-associated FIM projects.

[^1]: Database query and downstreams scripts updated August 2026 by M. Schram.
[^2]: Any updates to the official FIM species list will not automatically propagate to the TBNI index processing. If the FIM species list is updated (e.g., changes to NODCCodes), an updated species list will need to be provided to TBEP using `R/0_Update_SpeciesList.R`. This script produces an updated version of `TBIndex_spp_codes.rds`, which should be provided to TBEP alongside the other data. Similar to the `1_Query_FIM_Database.R`, the `0_Update_SpeciesList.R` can only be executed by FIM staff with appropriate database credentials. 