
USER <- Sys.getenv("USER") # Linux
if(Sys.info()[["sysname"]]=="Linux") DIR=paste0("/data/", USER, "/R_input_data")
if(Sys.info()[["sysname"]]=="Linux") N_CORES=10
# if(Sys.info()[["sysname"]]=="Windows") N_CORES=3

setwd(DIR)
OUT <- "DIR/R_output/"

if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("rgdal")) { install.packages("rgdal", dependencies = TRUE) ; library(rgdal)}
if (!require("maptools")) { install.packages("maptools", dependencies = TRUE) ; library(maptools)}
if (!require("rgeos")) { install.packages("rgeos", dependencies = TRUE) ; library(rgeos)}
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = TRUE) ; library(doParallel)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = TRUE) ; library(dplyr)}
if (!require("reshape")) { install.packages("reshape", dependencies = TRUE) ; library(reshape)}
if (!require("data.table")) { install.packages("data.table", dependencies = TRUE) ; library(data.table)}
if (!require("viridis")) { install.packages("viridis", dependencies = TRUE) ; library(viridis)}
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)} # scaling
if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("prioritizr")) { install.packages("prioritizr", dependencies = TRUE) ; library(prioritizr)}
if (!require("readr")) { install.packages("readr", dependencies = TRUE) ; library(readr)}
if (!require("stringr")) { install.packages("stringr", dependencies = TRUE) ; library(stringr)}
if (!require("R.utils")) { install.packages("R.utils", dependencies = TRUE) ; library(R.utils)}
if (!require("rgbif")) { install.packages("rgbif", dependencies = TRUE) ; library(rgbif)}
if (!require("taxize")) { install.packages("taxize", dependencies = TRUE) ; library(taxize)} # scaling
if (!require("rfishbase")) { install.packages("rfishbase", dependencies = TRUE) ; library(rfishbase)} # scaling

### Burn in waterways into DEM ###
# load open street waterways (source: )
wtr <- raster(paste0(DIR, "SA_waterways_1km.tif"))
wtr[is.na(wtr)] <- 0
# load digital elevation model (www.EarthEnv.org/topography)
elv <- raster(paste0(DIR, "elevation_1KMmd_GMTEDmd.tif"))
# burn in waterways into DEM
tr_elv <- elv - wtr
writeRaster(tr_elv , paste0(OUT, "elevation.tif"), overwrite=TRUE)

ela <- raster(paste0(OUT, "elevation.tif"))
# crop extent to model domain
ela_subset <- crop(ela, 
                     extent(-81.3, -55.3083333, -12.05, 12.4333333))

x11();plot(ela_subset)
writeRaster(ela_subset , paste0(OUT, "topology/study_area_elevation.tif"), overwrite=TRUE)

### load in species occurrence datasets ###
# load in GBIF data (source: )
fi = read.table(file = (paste0(DIR, "/Gbif/0342829-200613084148143.csv")), 
                header=TRUE, 
                sep = "\t", quote="", 
                fill=TRUE, 
                stringsAsFactors = FALSE, 
                na.strings = c("","NA"))
fi2 = read.table(file = (paste0(DIR, "/Gbif/0342859-200613084148143.csv")), 
                 header=TRUE, 
                 sep = "\t", 
                 quote="", 
                 fill=TRUE, 
                 stringsAsFactors = FALSE, 
                 na.strings = c("","NA"))
# load in fishnet data (source: fishnet.com)
fi3 <- read.table(file = (paste0(DIR, "/Fishnet/SearchResults--8585698942969590597.csv")), 
                  header = TRUE, 
                  sep = ",", 
                  fill = TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = c("", "NA"))
# load in Colombian taxonomists expert data (source:)
fi4 <- read.table(file = (paste0(DIR, "/expert_data/Peces.csv")),
                  sep = ",", 
                  fill=TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = c("","NA"))
# load in Specieslink data (source: specieslink.org.br)
fi5 <- read.table(file = (paste0(DIR, "/Specieslink/speciesLink-20210921055326-0000798.txt")), 
                  header=TRUE, 
                  sep = "\t", 
                  quote="", 
                  fill=TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = c("","NA"))
fi6 <- read.table(file = (paste0(DIR, "/Specieslink/speciesLink-20210921060014-0001860.txt")), 
                  header=TRUE, 
                  sep = "\t", 
                  quote="", 
                  fill=TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = c("","NA"))
fi7 <- read.table(file = (paste0(DIR, "/Specieslink/speciesLink-20210921060105-0002022.txt")), 
                  header=TRUE, 
                  sep = "\t", 
                  quote="", 
                  fill=TRUE, 
                  stringsAsFactors = FALSE, 
                  na.strings = c("","NA"))
# create standard species table
fi5 <- fi5[,c("scientificname", "family", "daycollected", "monthcollected", "yearcollected", "latitude", "longitude", "country", "basisofrecord")]
fi6 <- fi6[,c("scientificname", "family", "daycollected", "monthcollected", "yearcollected", "latitude", "longitude", "country", "basisofrecord")]
fi7 <- fi7[,c("scientificname", "family", "daycollected", "monthcollected", "yearcollected", "latitude", "longitude", "country", "basisofrecord")]
fi5 <- bind_rows(fi5, fi6, fi7)

# load in iDigBio data (source: iDigBio.org)
fi10 <- read.table(file = (paste0(DIR, "/iDigBio/iDigBio1/occurrence.csv")), 
                   header = TRUE, sep = ",", fill = TRUE, stringsAsFactors = FALSE, 
                   na.strings = c("", "NA"))
fi11 <- read.table(file = (paste0(DIR, "/iDigBio/iDigBio2/occurrence.csv")), 
                   header = TRUE, sep = ",", fill = TRUE, stringsAsFactors = FALSE, 
                   na.strings = c("", "NA"))
fi12 <- read.table(file = (paste0(DIR, "/iDigBio/iDigBio3/occurrence.csv")), 
                   header = TRUE, sep = ",", fill = TRUE, stringsAsFactors = FALSE, 
                   na.strings = c("", "NA"))
fi10 <- fi10[, c("dwc.scientificName", "gbif.canonicalName", "dwc.family", "idigbio.geoPoint", "idigbio.eventDate", "dwc.eventDate", "dwc.country", "dwc.basisOfRecord")]
fi11 <- fi11[, c("dwc.scientificName", "gbif.canonicalName", "dwc.family", "idigbio.geoPoint", "idigbio.eventDate", "dwc.eventDate", "dwc.country", "dwc.basisOfRecord")]
fi12 <- fi12[, c("dwc.scientificName", "gbif.canonicalName", "dwc.family", "idigbio.geoPoint", "idigbio.eventDate", "dwc.eventDate", "dwc.country", "dwc.basisOfRecord")]
fi10 <- bind_rows(fi10, fi11, fi12)
# separate coordinates in latitude and longitude columns
fi10 <- tidyr::separate(fi10, idigbio.geoPoint, into = c("latitude", "longitude"), sep = ",")  
fi10$latitude <- str_sub(fi10$latitude, 8, -1)
fi10$longitude <- str_sub(fi10$longitude, 8, -2)
fi10$latitude <- as.numeric(fi10$latitude)
fi10$longitude <- as.numeric(fi10$longitude)
fi10 <- fi10[is.na(fi10$latitude) == FALSE, ]
fi10 <- fi10[is.na(fi10$longitude) == FALSE, ]
# format date column to standard format
fi10$idigbio.eventDate <- substr(fi10$idigbio.eventDate, 0, 10)
fi10$idigbio.eventDate <- as.Date(fi10$idigbio.eventDate, format = "%Y-%m-%d")
fi10 <- fi10 %>% dplyr::mutate(year = lubridate::year(fi10$idigbio.eventDate),
                               month = lubridate::month(fi10$idigbio.eventDate),
                               day = lubridate::day(fi10$idigbio.eventDate))
fi10$dwc.basisOfRecord[fi10$dwc.basisOfRecord %in% c("fossilspecimen")] <- "FOSSIL_SPECIMEN" 
fi10$dwc.scientificName <- capitalize(fi10$dwc.scientificName)

### Data table preparation & aggregation ###
# set first row as header
colnames(fi4) = (fi4[1, ])
# delete first row, which is now column headers
fi4 <- (fi4[-1, ])

# exclude NA values from Latitude and Longitude column 
fi <- fi[is.na(fi$decimalLatitude) == FALSE,]
fi2 <- fi2[is.na(fi2$decimalLatitude) == FALSE,]
fi3 <- fi3[is.na(fi3$Latitude) == FALSE,]
i <- fi4$LOC_Lat == 0
ii <- fi4$LOC_Lon == 0
fi4$LOC_Lat[i] <- NA
fi4$LOC_Lon[ii] <- NA
fi4 <- fi4[is.na(fi4$LOC_Lat) == FALSE,]
fi4 <- fi4[is.na(fi4$LOC_Lon) == FALSE,]
j <- fi5$latitude == 0
jj <- fi5$longitude == 0
fi5$latitude[j] <- NA
fi5$longitude[jj] <- NA
fi5 <- fi5[is.na(fi5$latitude) == FALSE,]
fi5 <- fi5[is.na(fi5$longitude) == FALSE,]

# add column "BasisOfRecords" to enable merging of the different data tables
BasisOfRecord <- as.data.frame(matrix(ncol = 1, nrow = 922))
colnames(BasisOfRecord) <- "BasisOfRecord"
fi4 <- cbind(fi4, BasisOfRecord)

# change coordinates saved as characters to integers
options(digits = 9)
fi4$LOC_Lat <- as.numeric(fi4$LOC_Lat)
fi4$LOC_Lon <- as.numeric(fi4$LOC_Lon)
fi5$latitude <- as.numeric(fi5$latitude)
fi5$longitude <- as.numeric(fi5$longitude)
fi5 <- fi5[is.na(fi5$latitude) == FALSE,]

# order species data columns to the same format
fi <- fi[,c("species","family","day","month","year","decimalLatitude", "decimalLongitude","countryCode","basisOfRecord")]
fi2 <- fi2[,c("species","family","day","month","year","decimalLatitude", "decimalLongitude","countryCode","basisOfRecord")]
fi3 <- fi3[,c("ScientificName", "Family", "DayCollected", "MonthCollected", "YearCollected","Latitude", "Longitude", "Country", "BasisOfRecord")]
fi4 <- fi4[,c("NOMBRECIENTIFICORub", "TAXA_Familia", "DIA_revisado", "MES_revisado", "AÑO_revisado", "LOC_Lat", "LOC_Lon", "PAIS_revisado", "BasisOfRecord")]
fi5 <- fi5[,c("scientificname", "family", "daycollected", "monthcollected", "yearcollected", "latitude", "longitude", "country", "basisofrecord")]
fi10 <- fi10[,c("dwc.scientificName", "dwc.family", "day", "month", "year", "latitude", "longitude", "dwc.country", "dwc.basisOfRecord")]

# remove nonsensical day values and format day values
fi3$DayCollected[fi3$DayCollected %in% c("0", "01-05","07-09", "07-11", "14-15", "67", "78", "1213", "1415", "2629", "2930", "33229", "33310")] <- NA
fi3$DayCollected[fi3$DayCollected %in% c("01", "1")] <- "1"
fi3$DayCollected[fi3$DayCollected %in% c("03", "3")] <- "3"
fi3$DayCollected[fi3$DayCollected %in% c("04", "4")] <- "4"
fi3$DayCollected[fi3$DayCollected %in% c("05", "5")] <- "5"
fi3$DayCollected[fi3$DayCollected %in% c("06", "6")] <- "6"
fi3$DayCollected[fi3$DayCollected %in% c("07", "7")] <- "7"
fi3$DayCollected[fi3$DayCollected %in% c("08", "8")] <- "8"
fi3$DayCollected[fi3$DayCollected %in% c("09", "9")] <- "9"

# convert character column Day, Month and Year columns to integer
fi3 <- within(fi3, ("DayCollected" <- as.numeric(fi3$DayCollected)))
fi4 <- within(fi4, ("DIA_revisado" <- as.numeric(fi4$DIA_revisado )))
fi4 <- within(fi4, ("MES_revisado" <- as.numeric(fi4$MES_revisado)))
fi4 <- within(fi4, ("AÑO_revisado" <- as.numeric(fi4$AÑO_revisado)))

# unify column names among the different tables 
colnames(fi) <- c("Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")
colnames(fi2) <- c("Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")
colnames(fi3) <- c("Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")
colnames(fi4) <- c("Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")
colnames(fi5) <- c("Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")
colnames(fi10) <- c("Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")

fi3 <- within(fi3, ("Month" <- as.numeric(fi3$Month)))
fi3 <- within(fi3, ("Year" <- as.numeric(fi3$Year)))
fi3 <- within(fi3, ("Latitude" <- as.numeric(fi3$Latitude)))
fi3 <- within(fi3, ("Longitude" <- as.numeric(fi3$Longitude)))

# combine the tables into one data table
occ <- bind_rows(fi, fi2, fi3, fi4, fi5, fi10)

# set all year value == 0 to NA
occ$Year[occ$Year %in% c("0","1", "2","21", "101", "197", "700", "966", "988", "1000", "1013", "1070", "1071", "1100", "1111", "1194", "1201","2025", "2030", "2099", "2107", "2140", "2301", "9700", "33241")] <- NA

# set all month value == 0 to NA
occ$Month[occ$Month %in% "0"] <- NA
occ$Month[occ$Month %in% "134"] <- NA

# remove fossils, NA species and NA year values
indx_occ <- which(occ$BasisOfRecord == "FOSSIL_SPECIMEN")
indx_occ2 <- which(is.na(occ$Species))
indx_occ3 <- which(is.na(occ$Year))
occ <- occ[-c(indx_occ, indx_occ2, indx_occ3), ]

# assign consecutive ID values
ID <- 1:nrow(occ)
occ <- cbind(ID=ID, occ)

# reduce the number of colums to the once needed for further analysis
colums_occ <- occ[,c("ID","Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")]

# remove NA coordinates
occ <- occ[is.na(occ$Latitude) == FALSE, ]
occ <- occ[is.na(occ$Longitude) == FALSE, ]

# reduce number of columns to the once needed for further analysis
colums_occ <- occ[,c("ID","Species","Family","Day","Month","Year","Latitude","Longitude","Country","BasisOfRecord")]

# create spatial point occurrences
occ_sp <- SpatialPointsDataFrame(colums_occ[, c( "Longitude","Latitude")], 
                                 data = colums_occ)

# load in basin and sub-catchment rasters generated from grass script (name grass_script) 
basins <- raster(paste0(DIR, "basin_90.tif"))
subbasins <- raster(paste0(DIR, "subbasin_90.tif"))

################################################################################################################################################################################################################# add here code line croping basin extent to study area ########################################
######################################################### and save to disk subbasin study area and model extent ########################################
########################################################################################################################################################

basin <- raster(paste0(DIR, "subbasins_90.tif"))
# create data frame with sub-catchment id column
basin_id <- as.data.frame(basin)
basin_id <- na.omit(basin_id)
basin_id <- as.data.frame(basin_id[!duplicated(basin_id$subbasins_90),])
colnames(basin_id) <- "basin_id"
# model domain
basin_subset <- crop(basin, extent(-81.3, -55.3083333, -12.05, 12.4333333))
# study area
basin_subset <- crop(basin, extent(-81.73828, -66.48926, -4.66184, 12.4333333))

basin_ID <- unique(basin_subset)

indx = which(!(basin[] %in% basin_ID))
basin[indx] = NA
# basin_crop <- trim(basin, 
#                values= NA) 
# indx <- which(basin[] %in% basin_ID)
# basins <- basin[indx]
basins <-raster::trim(basin,
                      values= NA)

x11(); plot(basins)

basin_id <- as.data.frame(basins)
basin_id <- na.omit(basin_id)
basin_id <- as.data.frame(basin_id[!duplicated(basin_id$subbasins_90),])
colnames(basin_id) <- "basin_id"

bas_dat <- basin_id
bas_dat$basins <- basin_id$basin_id
bas_dat$amount <- 1

study_area <- reclassify(basin, bas_dat, right = NA)


raster::writeRaster(study_area, paste0(DIR, "study_crop_basin.tif"), overwrite = TRUE)

# create a stack e.g. for points extraction
sub_overlay <- stack(basins, 
                     subbasins)

# give the new spatial points the same projection/extent as used in the other raster layers of the study area
projection(occ_sp) <- projection(sub_overlay)

# overlay the spatial points with the basin layer to extract the values for spatial points of the respective basin raster
occ_sp_study <- raster::extract(sub_overlay, 
                                occ_sp, 
                                sp = TRUE)

# check if there are cells with NAs in basins but not in sub-catchments and vice versa
occ_study_dt <- as.data.table(occ_sp_study)
NAs_bas_sub <- occ_study_dt[, .N, by = .(basin_90, subbasin_90)]

# exclude all spatial points which do not overlay with a basin (they have a NA value in the previous generated basin colum), hence which are outside the study area (e.g. marine species occuring in the ocean)
occ_sp_study <- occ_sp_study[is.na(occ_sp_study$basin_90) == FALSE, ]
occ_sp_study <- occ_sp_study[is.na(occ_sp_study$subbasin_90) == FALSE, ]

# rewrite the spatial points in a data frame in order to further clean out data 
occ <- as.data.frame(occ_sp_study)

# find species records with the multiple identical coordinates and sampled at the same day/month/year - for Maxent not necessary: there only needed if a species was detected in a catchment
occ_dub_indx <- duplicated(occ[, 2:8])
occ_dub_indx <- which(occ_dub_indx %in% FALSE) 
occ_dub <- occ[-c(occ_dub_indx),]

### cleaning species list ###
# apply taxize to verify species names
occ_va <- validate_names(occ$Species)

# remove data with "invalid" species names
occ_val <- occ_va
occ_val <- na.omit(occ_val)
indx_fish <- which(occ$Species %in% c((occ_val)))
occ_val <- occ[-c(indx_fish), ]
indx_fish <- occ_val$Species
indx_fish <- which(occ$Species %in% c(indx_fish))
occ <- occ[-c(indx_fish), ]

# return information on ecosystem types
eco <- ecosystem(occ_va)

# look up in the column "Salinity" species listed under "saltwater"; index them and remove them from the species list
eco <- eco[, c("Species","Salinity", "EcosystemName", "EcosystemType", "Location"), ]
sal <- which(eco$Salinity %in% "saltwater")
sal <- eco[sal,]
indx_sal <- unique(sal$Species)
indx_sal <- which(occ$Species %in% c(indx_sal))
occ <- occ[-c(indx_sal), ]

# convert data.frame into data.table to receive the total number of occurrences per species
occ_dt <- as.data.table(occ)
Nspecies <- occ_dt[, .N, by = Species]

# index species which occur in less than five sub-catchments; these species will not be modeled in SDMs but used in the spatial prioritization
indx_occ4 <- Nspecies[Nspecies$N <5, ]
indx_occ4 <- indx_occ4$Species
indx_occ4 <- which(occ$Species %in% c(indx_occ4))
# remove species which occur in less than five sub-catchments from the species list
occ <- occ[-c(indx_occ4), ]

### check species list to occur in the Colombian freshwater checklist ###
# load in checklist
###################### put in this step!!! find it in another script!!! ############################################
####################################################################################################################


occ <- occ[, c("Species", "Longitude", "Latitude", "Family", "Day", "Month", "Year", "Country", "BasisOfRecord", "basin_90", "subbasin_90")]

### Preperation environmental predictor variables ###

# load in environmental data
bioclim_list <- list.files(path = (paste0(DIR, "bioclim")), 
                           pattern = "*.tif$", 
                           full.names = TRUE)
landcover_list <- list.files(path = (paste0(DIR, "landcover_data")),
                             pattern = "*.tif$",
                             full.names = TRUE)

topology_list <- list.files(path = (paste0(DIR, "topology_data")),
                             pattern = "*.tif$",
                             full.names = TRUE)

# load in flow accumulation
flow1k <- raster(paste0(DIR, "flow1k_averaged_1960-2015.tif"))
NAvalue(flow1k)
NAvalue(flow1k) <- -1

# define the extent to be the extent of the study area
exte <- extent(-81.3, -55.3083333, -12.05, 12.4333333)

bioclim_stack <- stack(bioclim_list)
landcover_stack <- stack(landcover_list)
topology_stack <- stack(topology_list)
flow1k_stack <- stack(flow1k)

extent(bioclim_stack) <- exte
bioclim_stack <- setExtent(bioclim_stack, exte, keepres = FALSE)
res(bioclim_stack) <- 0.00833333333

extent(landcover_stack) <- exte
landcover_stack <- setExtent(landcover_stack, exte, keepres = FALSE)
res(landcover_stack) <- 0.00833333333

extent(flow1k_stack) <- exte
flow1k_stack <- setExtent(flow1k_stack, exte, keepres = FALSE)
res(flow1k_stack) <- 0.00833333333

### aggregate environmental variables per sub-catchment ###

# do zonal statistics, calculating the mean and range (taking max - min values) for each environmental variables per sub catchtment
bioclim_mean <- zonal(bioclim_stack, subbasins, fun = mean, digits = 0, na.rm = TRUE)
landcover_mean <- zonal(landcover_stack, subbasins, fun = mean, digits = 0, na.rm = TRUE)
topology_mean <- zonal(topology_stack, subbasins, fun = mean, digits = 0, na.rm = TRUE)
flow1k_mean <- zonal(flow1k_stack, subbasins, fun = mean, digits = 1, na.rm = TRUE)

# rewrite as dataframe
bioclim_mean <- as.data.frame(bioclim_mean)
landcover_mean <- as.data.frame(landcover_mean)
topology_mean <- as.data.frame(topology_mean)
flow1k_mean <- as.data.frame(flow1k_mean)

bioclim_max <- zonal(bioclim_stack, subbasins, fun = max, digits = 0, na.rm = TRUE)
landcover_max <- zonal(landcover_stack, subbasins, fun = max, digits = 0, na.rm = TRUE)
topology_max <- zonal(topology_stack, subbasins, fun = max, digits = 0, na.rm = TRUE)
flow1k_max <- zonal(flow1k_stack, subbasins, fun = max, digits = 0, na.rm = TRUE)

bioclim_max <- as.data.frame(bioclim_max)
landcover_max <- as.data.frame(landcover_max)
topology_max <- as.data.frame(topology_max)
flow1k_max <- as.data.frame(flow1k_max)

# calculate environmental variable min values per sub-catchment
bioclim_min <- zonal(bioclim_stack, subbasins, fun = min, digits = 0, na.rm = TRUE)
landcover_min <- zonal(landcover_stack, subbasins, fun = min, digits = 0, na.rm = TRUE)

topology_min <- zonal(topology_stack, subbasins, fun = min, digits = 0, na.rm = TRUE)
flow1k_min <- zonal(flow1k_stack, subbasins, fun = min, digits = 0, na.rm = TRUE)

bioclim_min <- as.data.frame(bioclim_min)
landcover_min <- as.data.frame(landcover_min)
topology_min <- as.data.frame(topology_min)
flow1k_min <- as.data.frame(flow1k_min)

### check which one is correct ###
bioclim_range <- bioclim_max[, 2:5] - bioclim_min[, 2:5]
bioclim_range <- bioclim_max %>%  dplyr::mutate(annual_mean_temp = bioclim_max[, 2] - bioclim_min[, 2], 
                                                isothermality = bioclim_max[, 3] - bioclim_min[, 3],
                                                mean_diurnal_temp = bioclim_max[, 4] - bioclim_min[, 4],
                                                temp_seasonality = bioclim_max[, 5] - bioclim_min[, 5])
bioclim_range <- bioclim_range[, -c(2:5)]

landcover_range <- landcover_max[, 2:5] - landcover_min[, 2:5]
landcover_range <- landcover_max %>% dplyr::mutate(annual_mean_temp = landcover_max[, 2] - landcover_min[, 2], 
                                                isothermality = landcover_max[, 3] - landcover_min[, 3], 
                                                mean_diurnal_temp = landcover_max[, 4] - landcover_min[, 4],                                                        temp_seasonality = landcover_max[, 5] - landcover_min[, 5])
landcover_range <- landcover_range[, -c(2:5)]

topology_range <- topology_max[, 2:5] - topology_min[, 2:5]
topology_range <- topology_max %>%  dplyr::mutate(annual_mean_temp = topology_max[, 2] - topology_min[, 2], 
                                                isothermality = topology_max[, 3] - topology_min[, 3],
                                                mean_diurnal_temp = topology_max[, 4] - topology_min[, 4],
                                                temp_seasonality = topology_max[, 5] - topology_min[, 5])
topology_range <- topology_range[, -c(2:5)]
flow1k_range <- flow1k_max %>%  dplyr::mutate(flow_range = flow1k_max[, 2] - flow1k_min[, 2])
flow1k_range <- flow1k_range[, -2]

### --> rewrite bioclim, landcover and topology data to env_data and cbind them ###
####################### name basin layer ### !!!basin_id!!! ###  and make it the first column #######################

# Merge together
env <- cbind(bioclim_mean, 
             landcover_mean[,2], 
             topology_mean[,2],
             flow1k_mean[,2],
             bioclim_range[,2],
             landcover_range[,2],
             topology_range[,2],
             flow1k_range[,2])

####################################################################################################################
############## test again to cbind environmental stacks dataframes like bioclim per sub_catchment ##################
############################################ or write stacks to disk ###############################################
####################################################################################################################

###--- Get stream order and various stream length values for each basin ----
# https://grass.osgeo.org/grass78/manuals/addons/r.stream.order.html
myfile  <- list.files(paste0(DIR), pattern=".gpkg", all.files=F, full.names=F)
(layers <- ogrListLayers(paste0(DIR, "/", myfile[1]))) # "2" because the 2nd file of the list --> "stream_network_vector_GRASS_5.gpkg"
my_shp <- readOGR(paste0(DIR, "/", myfile[1]), layers[1])

str(my_shp@data)
get_these <- c("cat",  # rename to "basin_id"
               "strahler", # Strahler's stream order
               "length", # stream length
               "stright", # length of stream as stright line
               "sinusoid", # fractal dimension: stream length/stright stream length
               "cum_length", # length of stream from source
               "flow_accum", # upstream contributing area
               "out_dist",  #  distance of current stream init from outlet
               "source_elev", # elevation of stream init
               "outlet_elev", #elevation of stream outlet
               "elev_drop", # difference between source_elev and outlet_elev + drop outlet
               "out_drop", # drop at the outlet of the stream
               "gradient") # drop/length


my_shp <- my_shp[, get_these]
names(my_shp)[1] <- "basin_id"
stream_topo <- as.data.table(my_shp@data)
data.frame(layers=names(stream_topo)) # check names

rm(my_shp, layers, myfile); gc()

### Add stream hydrological data ###
env <- merge(env, stream_topo, by="basin_id")
dim(env)
summary(env)

# calculate log flow accumulation
env$log_flow <- log(env$flow+1)

### Scale and center all data (keep basin_id separated)
env <- cbind(env[, c("basin_id")], scale(env[,-c("basin_id")], center=T))

### Free RAM
rm(list=ls(pattern="env_")); rm(stream_topo); gc()

### correlation analysis ###

env_cor <- cor(env, method = "pearson", use = "complete.obs")
env_cor[lower.tri(env_cor, diag = TRUE)] <- NA
env_cor[env_cor == 1] <- NA
env_cor <- as.data.frame(as.table(env_cor))
env_cor <- na.omit(env_cor)

# remove one variable from highly correlated pairs
env <- env[,c("annual_temp", "diurnal_temp", "temp_seasonality", "barren_veg", "cultivated_veg", "flooded_veg", "elev_range", "elev_drop", "gradient", "flow_accum", "length", "sinusoid", "log_flow" )]

### Load in species data ###
# Species points either from disk or from above
sp_points <- shapefile(paste0(DIR, "species_data/spec_sp.shp"))
sp_points <- occ_sp_study
head(sp_points)
# remove whitespace
sp_points$Species <- gsub(" ", ".", sp_points$Species)  
spec_name <- as.data.frame(sp_points)
spec_name <- spec_name[!duplicated(spec_name$Species),]

names(subbasins) <- "basin_id"

### Run thorugh all species in a loop aggerageting species occurrences to sub-catchments ###
my_species <- c(spec_name$Species)
names(subbasins) <- "basin_id"
sp_dt <- as.data.table(basin_id, na.rm=T)
sp_dt <- sp_dt[!duplicated(sp_dt),]

for (i in 1:length(my_species)) {
  cat("Running species", i, "\n")
  tmp <- as.data.table(extract(basin_id, subset(sp_points, Species==my_species[i]) , df=T))
  tmp <- tmp[!is.na(tmp$basin_id),]
  tmp <- tmp[!duplicated(tmp$basin_id),] # remove duplicates
  tmp$spec <- 1 # write species occurrence
  names(tmp) <- c("ID", "basin_id", my_species[i] )
  sp_dt <- merge(sp_dt, tmp[,c("basin_id",  my_species[i]), with=F] , by="basin_id", all.x=T)
  rm(tmp)
}

### Merge species and environmental data ###
all_data <- merge(env, sp_dt, by="basin_id", all.x=T) # add species data
head(all_data)
tail(all_data)
str(all_data)

### write to disk to use table in ensemble modelling ###
write.csv(all_data, 
          file = (paste0(OUT, "all_data.csv")),
          row.names = FALSE)
