# Libraries and source files
# Load geocoded data
library(data.table)
library(dplyr)
library(RecordLinkage)
library(igraph)
#source("./src/mcmcLinkReconcile/00-linkReconcileFunctions.R")

# Load in wvs data, filter columns, and clean

wvsData <- fread("./data/dhs_link_char/working/wvs_geocoded.csv", check.names = TRUE, colClasses = "character")[, .(last_name, first_name, gender, ssn, street_number, route, birth_yr, birth_mo, birth_dy)]
wvsData[, c("last_name", "first_name") := list(toupper(stringr::str_trim(last_name)), toupper(stringr::str_trim(first_name)))]
wvsData[gender == "", gender := NA]
wvsData[last_name == "", last_name := NA]
wvsData[first_name == "", first_name := NA]
wvsData[ssn == "", ssn := NA]

# Remove identical records
wvsData = wvsData[, .N, by = .(last_name, first_name, gender, ssn, street_number, route, birth_yr, birth_mo, birth_dy)][,N := NULL]

# Load in anz data, filter columns, and clean

anzData <- fread("./data/dhs_link_char/working/anz_geocoded.csv", check.names = TRUE, colClasses = "character")[, .(last_name, first_name, gender, ssn, street_number, route, birth_yr, birth_mo, birth_dy)]
anzData[gender == "", gender := NA]
anzData[last_name == "", last_name := NA]
anzData[first_name == "", first_name := NA]
anzData[ssn == "", ssn := NA]

# Remove identical records
anzData = anzData[, .N, by = .(last_name, first_name, gender, ssn, street_number, route, birth_yr, birth_mo, birth_dy)][,N := NULL]

dataTableNames = c("wvsData", "anzData")

demographics = rbindlist(lapply(dataTableNames, function(name){
  #data.table(source = name, get(name))
  data.table(get(name))
}))

rm(dataTableNames, wvsData, anzData)
