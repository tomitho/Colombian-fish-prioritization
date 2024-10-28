
USER <- Sys.getenv("USER") # Linux
if(Sys.info()[["sysname"]]=="Linux") DIR=paste0("/data/", "/spatial_prioritization")
if(Sys.info()[["sysname"]]=="Linux") N_CORES=10

setwd(DIR)
OUT <- "DIR/output/"

if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = TRUE) ; library(dplyr)}
if (!require("data.table")) { install.packages("data.table", dependencies = TRUE) ; library(data.table)}
if (!require("readr")) { install.packages("readr", dependencies = TRUE) ; library(readr)}

# load in species suitable habitat probabilities
puvspr_dat <- fread(paste0(DIR, "puvspr_dat", ".csv"), h=T)
# load in sub-catchment data of model domain extent
basin <- raster(paste0(DIR, "model_domain", "subbasins_90.tif"))
# create data frame with sub-catchment id column
basin_id <- as.data.frame(basin)
basin_id <- na.omit(basin_id)
basin_id <- as.data.frame(basin_id[!duplicated(basin_id$subbasins_90),])
colnames(basin_id) <- "basin_id"

# group number of species having suitable habitat per sub-catchment
spec_prob <- puvspr_dat %>% group_by(pu) %>% mutate(sum_amount = sum(amount)) %>% summarise_if(is.numeric, sum)
spec_prob <- spec_prob[, c("pu", "amount")]
spec_prob$basin <- spec_prob$pu
spec_prob <- spec_prob[, c("basin", "pu", "amount")]

# index locked in planning units
indx <- which(pu_dat$locked_in %in% TRUE)
locki <- pu_dat[c(indx),]
indx2 <- which(spec_prob$basin %in% locki$id)

# look up rows and add + 1 for the endemic species which were not modeled since they only occur in that one sub-catchment
spec_prob[386,3] <- 7
spec_prob[600,3] <- 20
spec_prob[1144,3] <- 3
spec_prob[5352,3] <- 278
spec_prob[6936,3] <- 62
spec_prob[7369,3] <- 29
spec_prob[15813,3] <- 228

# index sub_catchments which contain suitable habitats
bas_indx <- which(basin_id$basin_id %in% spec_prob$basin)

# identify sub_catchments without suitable habitat
basin_id_mis <- as.data.frame(basin_id[-c(bas_indx),])
colnames(basin_id_mis) <- "basin"
basin_id_mis$pu <- basin_id_mis$basin

# assign not suitable habitat the value 0
basin_id_mis$amount <- 0

spec_prob1 <- rbind(spec_prob, basin_id_mis)
# reclassify raster
spec_prob1 <- reclassify(basin, spec_prob1, right = NA)

# save to disk
writeRaster(spec_prob1,
            paste0(OUT, "stacked_habitat_suitability.tif"),
            overwrite = TRUE)
