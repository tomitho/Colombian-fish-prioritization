
USER <- Sys.getenv("USER") # Linux
if(Sys.info()[["sysname"]]=="Linux") DIR=paste0("/data/", USER, "parallel_run")
if(Sys.info()[["sysname"]]=="Linux") N_CORES=10

setwd(DIR)
OUT <- "DIR/output/"

if (!require("raster")) { install.packages("raster", dependencies = TRUE) ; library(raster)}
if (!require("foreach")) { install.packages("foreach", dependencies = TRUE) ; library(foreach)}
if (!require("plyr")) { install.packages("plyr", dependencies = TRUE) ; library(plyr)}
if (!require("dplyr")) { install.packages("dplyr", dependencies = TRUE) ; library(dplyr)}
if (!require("data.table")) { install.packages("data.table", dependencies = TRUE) ; library(data.table)}
if (!require("readr")) { install.packages("readr", dependencies = TRUE) ; library(readr)}

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

# load in species data list
all_data <- fread(paste0(DIR, "/all_data", ".csv"), h=T)
MY_SPP <- all_spec$species

# create empty table with an random example species
var_imp <- fread(paste0(DIR, "spp_ensemble", "/",  "Trachydoras.microstomus" ,"/",  "Trachydoras.microstomus", "_variable_importance", ".csv"), h=T)
var_imp <-as.data.frame(var_imp[, var_imp$variable])
sp_dt <- as.data.table(var_imp, na.rm=T)
colnames(sp_dt) <- c("variable")
tmp <- NULL

### loop through all calculated variable importance of all species and add to one table
for (i in 1:length(MY_SPP)) {
  if (file.exists(paste0(getwd(), "spp_ensemble",  MY_SPP[i] ,"/",  MY_SPP[i], "_variable_importance.csv"))) {
    var_imp <- fread(paste0(DIR, "spp_ensemble", "/",  MY_SPP[i] ,"/",  MY_SPP[i], "_variable_importance", ".csv"), h=T)
    cat("Running species", i, "\n")
    
    tmp <- var_imp[, c("variable", "percent_contribution")]
    colnames(tmp) <- c("variable", MY_SPP[i]) # need to index
    sp_dt <- merge(sp_dt, tmp[,c("variable",  MY_SPP[i]), with=F] , by="variable", all.x=T)
    rm(var_imp)
    rm(i)
  }
}

# sum up the variable importance of all species
sp_dt_sum <- rowSums(sp_dt[,2:1307])
sp_dt_sum <- as.data.frame(sp_dt_sum)
sp_dt_sum$sp_dt_sum <- sp_dt_sum$sp_dt_sum/1307
sp_dt_sum$variable <- sp_dt$variable
sp_dt_sum <- sp_dt_sum[, c("variable", "sp_dt_sum")]
colnames(sp_dt_sum) <- c("variable", "var_imp_sum")

sp_dt_sum$var_imp_sum <- sp_dt_sum$var_imp_sum*100
# save to disk
write.csv(sp_dt_sum, file = paste0(OUT, "variable_importance.csv"), row.names = FALSE)