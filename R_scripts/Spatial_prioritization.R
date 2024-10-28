# run in linux before using gurobi
source ~/.bashrc
source ~/.profile
source ~/.bash_profile

USER <- Sys.getenv("USER") # Linux
if(Sys.info()[["sysname"]]=="Linux") DIR=paste0("/data/", "/spatial_prioritization")
if(Sys.info()[["sysname"]]=="Linux") N_CORES=10

setwd(DIR)
OUT <- "DIR/output/"

if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("doParallel")) { install.packages("doParallel", dependencies = TRUE) ; library(doParallel)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = TRUE) ; library(dplyr)}
if (!require("reshape")) { install.packages("reshape", dependencies = TRUE) ; library(reshape)}
if (!require("data.table")) { install.packages("data.table", dependencies = TRUE) ; library(data.table)}
if (!require("viridis")) { install.packages("viridis", dependencies = TRUE) ; library(viridis)}
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE) ; library(RColorBrewer)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) ; library(ggplot2)}
if (!require("ggthemes")) { install.packages("ggthemes", dependencies = TRUE) ; library(ggthemes)} # theme_map()
if (!require("caret")) { install.packages("caret", dependencies = TRUE) ; library(caret)} # theme_map()
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)} # scaling
if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("prioritizr")) { install.packages("prioritizr", dependencies = TRUE) ; library(prioritizr)}
if (!require("readr")) { install.packages("readr", dependencies = TRUE) ; library(readr)}
if (!require("lwgeom")) { install.packages("lwgeom", dependencies = TRUE) ; library(lwgeom)}
if (!require("gurobi")) { install.packages("gurobi", dependencies = TRUE) ; library(gurobi)}


# ================================== #
### Creating species data spec_dat ###
# ================================== #

# mention somewhere that "Leporinus.y-ophorus" is renamed to "Leporinus.yophorus" ######################################################################
# load in table containing all species list created in data preperation script
all_data <- fread(paste0(DIR, "all_data", ".csv"), h=T)
MY_SPP <- all_spec$species

# create a matrix/data.frame to be used for the species list needed in Marxan to idendify in the other tables
spec_dat <- matrix(nrow = 1307, ncol =1 , byrow = TRUE) # create a matrix with as many rows as there are species
spec_dat <- as.data.frame(spec_dat)
colnames(spec_dat) <- c("name")
spec_dat$name <- MY_SPP
spec_dat$id <- seq.int(1:nrow(spec_dat))
spec_dat <- spec_dat[, c("id", "name")]

# ====================================== #
### Creating planning unit data pu_dat ###
# ====================================== #
# load in sub-catchment data of study area extent
basin <- raster(paste0(DIR, "study_area", "subbasins_90.tif"))
# create data frame with sub-catchment id as planning unit id
basin_id <- as.data.frame(basin)
basin_id <- na.omit(basin_id)
basin_id <- as.data.frame(basin_id[!duplicated(basin_id$subbasins_90),])
colnames(basin_id) <- "basin_id"

# calculate area per sub-catchment 
basin_area <- area(basin, na.rm=TRUE, weights=FALSE)
basin_dat <- as.data.frame(basin)
area <- as.data.frame(basin_area)
basin_dat <- na.omit(basin_dat)
area <- na.omit(area)
area_basin <- cbind(basin_dat, area)
b <- list(area_basin$subbasins_90)

# combine planning unit id and area per sub-catchment
pu_dat <- aggregate.data.frame(area_basin, by=b, FUN="sum")
pu_dat <-  pu_dat[, -2]
colnames(pu_dat) <- (c("id", "cost"))
# scale pu_dat area cost column
pu_dat$cost <- (pu_dat$cost - min(pu_dat$cost))/(max(pu_dat$cost)-min(pu_dat$cost))

# load in human footprint index data
humanfoot <- raster(paste0(DIR, "humanfoot_cost.tif"))
# scale of human footprint goes until 50, set everything above to NA
humanfoot[humanfoot >= 50] <- NA
# aggregate human footprint index per sub-catchment
humanfoot_mean <- zonal(humanfoot, basin, fun = mean, digits = 0, na.rm = TRUE)
humanfoot_mean <- as.data.frame(humanfoot_mean)

# sub-basins 3422 & 11161 were projected by the human footprint index to be open water which was previsouly set from the value 128 to NA but on the map there is actually land and thus gets an arbitrary number based on visual insepection in QGIS from their surrounding pixels
humanfoot_mean$value[c(3324)] <- 0
humanfoot_mean$value[c(10487)] <- 3
humanfoot_mean$value[c(10552)] <- 0
# set human footprint index to be used as penalty data in spatial prioritization
colnames(humanfoot_mean) <- (c("id", "penalty_data"))
# scale human footprint index
humanfoot_mean$penalty_data <- (humanfoot_mean$penalty_data - min(humanfoot_mean$penalty_data))/(max(humanfoot_mean$penalty_data)-min(humanfoot_mean$penalty_data))
# add human footprint index to pu_dat table
pu_dat$penalty_data <- humanfoot_mean$penalty_data
# check if human footprint penalty data and area cost were changed in the end to be HFI = cost and area = penalty data? ################################
colnames(pu_dat) <- c("id", "penalty_data", "cost")
pu_dat <- pu_dat[, c("id", "cost", "penalty_data")]

########################################################################################################################################################
### add pu_dat locked_in column to set the sub_catchment ids containing species occurring only in a single sub_catchment to be locked in the priority set ####################################################################################################################################################


# ======================= #
### Creating spatial probability data based on SDM probabilities of biological features puvspr_dat ###
# ======================= #
# create empty tmp and puvspr_dat data.frames to be used in the loop
tmp <- NULL
puvspr_dat <- as.data.frame(matrix(nrow = 0, ncol = 4))
colnames(puvspr_dat) <- (c("id", "probs", "amount", "species"))

# loop through ensemble model output based on all predictions data
for (i in 1:length(MY_SPP)) {
  if (file.exists(paste0(getwd(), "spp_ensemble",  MY_SPP[i] ,"/",  MY_SPP[i], "_all_predictions_basin_table.RData"))) {
    # load(paste0(DIR, "/spp_ensemble_output","/", MY_SPP[i],"/", MY_SPP[i], "_all_predictions_basin_table.RData"))
    load(paste0(DIR, "spp_ensemble", "/", MY_SPP[i],"/", MY_SPP[i], "_all_predictions_basin_table.RData"))
    cat("Running species", i, "\n")
    all_predictions <- all_predictions[all_predictions$binary == 1]
    all_predictions$Species <- (MY_SPP[i]) # need to index
    tmp[i] <- list(all_predictions)
    puvspr_dat <- rbindlist(tmp)
    puvspr_dat$species <- puvspr_dat %>% group_indices(Species)
    rm(all_predictions)
    rm(i)
  }
}

# loop with ensemble model output based on weighted distances of suitable sub_basins to occurrence data
tmp <- NULL
puvspr_dat <- as.data.frame(matrix(nrow = 0, ncol = 4))

for (i in 1:length(MY_SPP)) {
  if (file.exists(paste0(getwd(), "spp_ensemble",  MY_SPP[i] ,"/",  MY_SPP[i], "_all_predictions_basin_table.RData"))) {
    load(paste0(DIR, "spp_ensemble", "/", MY_SPP[i],"/", MY_SPP[i], "_all_predictions_basin_table.RData"))
    load(paste0(DIR, "spp_ensemble","/", MY_SPP[i],"/", MY_SPP[i], "_basin_distance_weight_table.RData"))
    cat("Running species", i, "\n")

    basin_indx <- basin_id$basin_id
    basin_indx <- which(basin_for_dist$basin_id %in% basin_indx)
    basin_for_dist <- basin_for_dist[c(basin_indx),]

    basin_indx <- basin_id$basin_id
    basin_indx <- which(all_predictions$basin_id %in% basin_indx)
    all_predictions <- all_predictions[c(basin_indx),]

    basin_for_dist <- basin_for_dist %>% group_by(basin_id) %>% mutate(mean_distance_scaled = mean(distance_scaled)) %>% summarise_if(is.numeric, mean)

    ensemble_eva <- fread(paste0(DIR, "spp_ensemble", "/", MY_SPP[i], "/ensemble_evaluation", ".csv"), h=T)
    ensemble_eva <- ensemble_eva[[1,3]]
    ensemble_eva <- ensemble_eva/1000

    all_predictions$probs_scaled <- all_predictions$probs * basin_for_dist$distance_scaled

    all_predictions$binary_scale <- all_predictions$probs_scaled
    all_predictions$binary_scale[c(all_predictions$binary_scale <= ensemble_eva)] <- 0
    all_predictions$binary_scale[c(all_predictions$binary_scale > ensemble_eva)] <- 1
    all_predictions <- all_predictions[all_predictions$binary_scale == 1]

    all_predictions$Species <- (MY_SPP[i]) # need to index
    tmp[i] <- list(all_predictions)
    puvspr_dat <- rbindlist(tmp)
    puvspr_dat$species <- puvspr_dat %>% group_indices(Species)
    puvspr_dat <- puvspr_dat %>% select(-Species)
    puvspr_dat <- puvspr_dat[, c("species", "basin_id", "probs_scaled", "binary_scale")]
    colnames(puvspr_dat) <- (c("species", "pu", "prob", "amount"))
    rm(basin_for_dist, all_predictions, ensemble_eva)
    rm(i)
  }
}

######################################################################################################################################################################################################################## check if identical to above line ######################################################
########################################################################################################################################################

#' Transform raster to data.table
#'
#' @param x  Rastereco_behav_ object
#' @param row.names  `NULL` or a character vector giving the row names for the data frame. Missing values are not allowed
#' @param optional  logical. If `TRUE`, setting row names and converting column names (to syntactic names: see make.names) is optional
#' @param xy  logical. If `TRUE`, also return the spatial coordinates
#' @param centroids  logical. If TRUE return the centroids instead of all spatial coordinates (only relevant if xy=TRUE)
#' @param sepNA	logical. If TRUE the parts of the spatial objects are separated by lines that are NA (only if xy=TRUE and, for polygons, if centroids=FALSE
#' @param ...	 Additional arguments (none) passed to `raster::as.data.frame`
#'
#' @value returns a data.table object
#' @examples
#' logo <- brick(system.file("external/rlogo.grd", package="raster"))
#' v <- as.data.table(logo)
#' @import
as.data.table.raster <- function(x, row.names = NULL, optional = FALSE, xy=FALSE, inmem = canProcessInMemory(x, 2), ...) {
  stopifnot(require("data.table"))
  if(inmem) {
    v <- as.data.table(raster::as.data.frame(x, row.names=row.names, optional=optional, xy=xy, ...))
  } else {
    tr <- blockSize(x, n=2)
    l <- lapply(1:tr$n, function(i)
      as.data.table(as.data.frame(getValues(x,
                                            row=tr$row[i],
                                            nrows=tr$nrows[i]),
                                  row.names=row.names, optional=optional, xy=xy, ...)))
    v <- rbindlist(l)
  }
  coln <- names(x)
  if(xy) coln <- c("x", "y", coln)
  setnames(v, coln)
  v
}
if (!isGeneric("as.data.table")) {
  setGeneric("as.data.table", function(x, ...)
    standardGeneric("as.data.table"))
}
setMethod('as.data.table', signature(x='data.frame'), data.table::as.data.table)
setMethod('as.data.table', signature(x='Raster'), as.data.table.raster)

sp_dt <- as.data.table(basin, na.rm=T)
sp_dt <- sp_dt[!duplicated(sp_dt),]
tmp <- NULL
names(basin) <- "basin_id"

for (i in 1:length(MY_SPP)) {
  load(paste0(DIR, "spp_ensemble","/", MY_SPP[i],"/", MY_SPP[i], "_all_predictions_basin_table.RData"))
  load(paste0(DIR, "spp_ensemble","/", MY_SPP[i],"/", MY_SPP[i], "_basin_distance_weight_table.RData"))
  cat("Running species", i, "\n")
  basin_for_dist <- basin_for_dist %>% group_by(basin_id) %>% mutate(mean_distance_scaled = mean(distance_scaled)) %>% summarise_if(is.numeric, mean)
  indx_for_dist <- which(is.na(basin_for_dist$basin_id))
  basin_for_dist <- basin_for_dist[-c(indx_for_dist),]

  ensemble_eva <- fread(paste0(DIR,"/Parallel_Run", "/", MY_SPP[i], "/ensemble_evaluation", ".csv"), h=T)
  ensemble_eva <- ensemble_eva[[1,3]]
  ensemble_eva <- ensemble_eva/1000

  all_predictions$probs_scaled <- all_predictions$probs * basin_for_dist$distance_scaled

  all_predictions$binary_scale <- all_predictions$probs_scaled
  all_predictions$binary_scale[c(all_predictions$binary_scale <= ensemble_eva)] <- 0
  all_predictions$binary_scale[c(all_predictions$binary_scale > ensemble_eva)] <- 1


  tmp <- all_predictions[, c("basin_id", "binary_scale")]
  colnames(tmp) <- c("basin_id", MY_SPP[i]) # need to index
  sp_dt <- merge(sp_dt, tmp[,c("basin_id",  MY_SPP[i]), with=F] , by="basin_id", all.x=T)
  rm(basin_for_dist, all_predictions, ensemble_eva)
  rm(i)
}
tmp <- as.data.frame(all_predictions)
tmp <- all_predictions[, c("basin_id", "binary_scale")]
colnames(tmp) <- c("basin_id","Acaronia.nassa")

################################################################################################################################################################################################################################################################################################################


# load in stream routing data to calculate connectivity bmat
load(paste0(DIR, "hydrology/stream_routing_upstream_m_distance.RData"))

# add connectivity matrix
connectivity <- out[, -c("seq_id")]
connect_cat <- which(connectivity$cat %in% basin_id$basin_id)
connect_cat_up <- which(connectivity$cat_up %in% basin_id$basin_id)
connectivity <- connectivity[c(connect_cat, connect_cat_up),]
connectivity_indx <- duplicated(connectivity[, 1:3])
connectivity_indx <- which(connectivity_indx %in% FALSE)
connectivity <- connectivity[-c(connectivity_indx),]
# set 200 km threshold
connectivity <- connectivity[connectivity$dist_m <= 200000]
connectivity$dist_m <- 200000 - connectivity$dist_m
bmat <- connectivity
# set required column names and scale connectivity data
colnames(bmat) <- (c("id1", "id2", "boundary"))
bmat$boundary <- (bmat$boundary - min(bmat$boundary))/(max(bmat$boundary)-min(bmat$boundary))


################################################################################################################################################################################################################################################################################################################

### create planning unit data setting sub_catchments containing at least 50 % protected areas to be locked in the priority area network
# load in protected area data
# Colombia protected areas
pa_col <- fread(paste0(DIR, "/data/tomiczek/Publication/Gurobi/protected_area/protected_tables/col_protected_area_per_basin", ".csv"), h=T)
pa_col <- pa_col[pa_col$area_percent >= 50]
# Brazil protected areas
pa_bra <- fread(paste0(DIR, "bra_protected_area_per_basin", ".csv"), h=T)
pa_bra <- pa_bra[pa_bra$area_percent >= 50]
# Ecuador protected areas
pa_ecu <- fread(paste0(DIR, "ecu_protected_area_per_basin", ".csv"), h=T)
pa_ecu <- pa_ecu[pa_ecu$area_percent >= 50]
# Peru protected areas
pa_per <- fread(paste0(DIR, "per_protect_per_basin", ".csv"), h=T)
pa_per <- pa_per[pa_per$area_percent >= 50]
# Venezuela protected areas
pa_ven <- fread(paste0(DIR, "ven_protected_area_per_basin", ".csv"), h=T)
pa_ven <- pa_ven[pa_ven$area_percent >= 50]

# combine all protected area tables
pa_tot <- rbind(pa_col, pa_bra, pa_ecu, pa_per, pa_ven)

### add how the existing protected areas are locked in the pu_dat to be locked in and call it pu_dat_locked_in

################################################################################################################################################################################################################################################################################################################

###                        ###
### Spatial prioritization ###
###                        ###
# ========================= #
### Add a problem statement ###
# ========================= #

### free choice spatial prioritization ###

# 20 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.2) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)
# 30 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.3) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)
# 40 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.4) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)
# 50 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.5) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)

### locked in spatial prioritization ###

# 20 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.2) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)
# 30 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.3) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)
# 40 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.4) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)
# 50 % conservation target
p <- problem(pu_dat, spec_dat, cost_colum = "cost", rij = puvspr_dat) %>% add_min_set_objective() %>% add_relative_targets(0.5) %>% add_binary_decisions() %>% add_locked_in_constraints(pu_dat$locked_in) %>%  add_linear_penalties(2, data = "penalty_data") %>% add_boundary_penalties(penalty = 0.03, data = bmat) %>% add_gurobi_solver(gap = 0.1)


# solve the problem statement
s <- solve(p)

# save gurobi spatial prioritization to disk
write.csv(s, file = paste0(DIR, "0.2", "/gurobi_solution.csv"), row.names = FALSE)

# ========================================================================= #
### Transform the table output solution to spatially explicit raster format ###
# ========================================================================= #
### reclassify the sub_catchment raster to contain a 1 if it is a priority area ###
# index the prioritized sub-catchments
l1 <- which(s$solution_1 == 1)
# select those rows with the prioritized sub-catchments from the solution table
mat <- s[c(l1),]
# remove the costs and status columns
mat <- mat[, -(2)]
# copy id column
mat2 <- mat$id
mat3 <- cbind(mat, mat2)
mat3 <- mat3[, c("id", "mat2", "solution_1")]

# reclassify sub_catchment raster
basin_rec <- reclassify(basin, mat3, right = NA)

solution_rast <- basin_rec
values(solution_rast)[!values(solution_rast) == 1] = NA

# save reclassified raster showing the priority network to disk
writeRaster(solution_rast, paste0(DIR, "0.2", "/solution.tif"), overwrite = TRUE)