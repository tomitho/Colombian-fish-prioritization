#============================================================================== Biomod Ensemble Run =====================================================================================#

USER <- Sys.getenv("USER") # Linux
if(Sys.info()[["sysname"]]=="Linux") DIR=paste0("/data/", USER, "parallel_run")
if(Sys.info()[["sysname"]]=="Linux") N_CORES=10
dir.create(DIR)
setwd(DIR)
MAXENT_DIR=paste0(DIR, "/MAXENT_DIR")

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
if (!require("biomod2")) { BiocManager::install("biomod2") ; library(biomod2)}
if (!require("viridis")) { install.packages("viridis", dependencies = TRUE) ; library(viridis)}
if (!require("RColorBrewer")) { install.packages("RColorBrewer", dependencies = TRUE) ; library(RColorBrewer)}
if (!require("ggplot2")) { install.packages("ggplot2", dependencies = TRUE) ; library(ggplot2)}
if (!require("ggthemes")) { install.packages("ggthemes", dependencies = TRUE) ; library(ggthemes)} # theme_map()
if (!require("caret")) { install.packages("caret", dependencies = TRUE) ; library(caret)} # theme_map()
if (!require("scales")) { install.packages("scales", dependencies = TRUE) ; library(scales)} # scaling

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


###-------------------------------#
###--- Run an ensemble model -----
###-------------------------------#

all_data <- fread(paste0(DIR, "all_data", ".csv"), h=T)

### Create grid_id layer --> needed later for getting the raster predictions
domain_cells <- as.data.table.raster(basin_id)
domain_cells$seq_id <- seq.int(1:nrow(domain_cells))
setkey(domain_cells, seq_id) # set key

### Download Maxent from within R
download.file(url = "https://biodiversityinformatics.amnh.org/open_source/maxent/maxent.php?op=download", destfile = "maxent.zip", mode='wb')
### Unzip
MAXENT_DIR=paste0(DIR, "/MAXENT_DIR")
unzip("maxent.zip", exdir=MAXENT_DIR, overwrite = T)

###--- Prepare the data for biomod2----
EXP <- all_data[,c("annual_temp", "diurnal_temp", "temp_seasonality", "barren", "cultivated_veg", "flooded_veg", "elevation_range", "log_flow", "length", "sinusoid", "flow_accum","elev_drop", "gradient")]
SPP_NA <- all_data[,..my_species] # get table of species data
COORD <- cbind(all_data[,"basin_id"], all_data[,"basin_id"]); names(COORD) <- c("x", "y")

###--- Select if stream network distance calculations should be done ----
DO_DISTANCE=TRUE
# DO_DISTANCE=FALSE # set to false if it should not be done

### Setup parallel computation
cl <- parallel::makeCluster(N_CORES, outfile = "")
registerDoParallel(cl)

###--- Run in parallel ----  
foreach(i=1:length(SPP_NA), 
        .errorhandling = c("stop"), 
        .verbose = T, 
        .packages = c("raster", "plyr", "dplyr", "data.table", "biomod2", "scales", "viridis")) %dopar% {
          # .errorhandling = c("stop", "remove", "pass")
          
          MY_SPP <- names(SPP_NA[,i,with=F]) # get species name
          MY_SPP_DIR <- paste0(DIR, "/", MY_SPP) # create species folder path
          system(paste0("rm -rf ", MY_SPP_DIR))   # remove folder (previous model run)
          dir.create(MY_SPP_DIR)  # create species folder
          
          ### Copy maxent folder into each species' directory
          file.copy(MAXENT_DIR, MY_SPP_DIR, recursive=TRUE)
          MAXENT_DIR_TMP <- paste0(MY_SPP_DIR, "/MAXENT_DIR")
          
          cat(">>>>>>  Running models for", MY_SPP, "   <<<<<<", "\n")
          
          
          ##-----------------------------------------------------------#
          ####--- Calculate stream-network distance between points ----
          ##-----------------------------------------------------------#
          
          if(DO_DISTANCE==TRUE) { 
            
            ### Load stream network raster
            str_net <- raster(paste0(DIR, "stream_network_GRASS_90", ".tif")) 
            names(str_net) <- "basin_id"
            
            ### Add a "1" is species occurs in the basin that corresponds to the stream reach (no "snapping" needed)
            str_net_tmp <- as.data.table(str_net, na.rm=F) # get data.table
            str_net_tmp$seq_id <- seq.int(1:nrow(str_net_tmp)) # add ID
            setkey(str_net_tmp, seq_id)
            
            ### Prepare species data
            tmp_spp <- sp_dt[,c("basin_id", MY_SPP), with=F]
            tmp_spp <- na.omit(tmp_spp)
            set_these_to_1 <- sort(unique(tmp_spp$basin_id)) # get those basins that have species observations
            
            ### Set the values in the stream network
            str_net_tmp$for_distance <- ifelse(is.na(str_net_tmp$basin_id), NA, 0) # set all stream cells to zero, land remains NA
            
            ### Get index of "set_these_to_1" basins, and repalce those with "1"
            my_index <- which(str_net_tmp$basin_id  %in% set_these_to_1)
            str_net_tmp$for_distance <- replace(str_net_tmp$for_distance, my_index, 1) # write "1" in those basins that have species
            str_net_tmp <- as.data.table(arrange(str_net_tmp, seq_id))
            # subset(str_net_tmp, str_net_tmp$basin_id==8019) # check
            
            ### Insert the new values into the raster (0=no species, 1=species observed in that basin/stream)
            str_net_for_dist <- setValues(basin_id, str_net_tmp$for_distance)
            # writeRaster(str_net_for_dist, paste0(DIR, "/temp_distance_file.tif"), overwrite=T) # check 
            
            ### Calculate the within-stream distance between points
            str_net_for_dist <- gridDistance(str_net_for_dist, 1, omit=NA)
            writeRaster(str_net_for_dist, paste0(DIR, "/", MY_SPP, "/str_net_distance.tif"), overwrite=T)
            
            
            ### Assign the values to the entire basin raster file
            dist_per_bas <- as.data.table(as.data.frame(zonal(str_net_for_dist, basin_id, fun='mean', digits=2, na.rm=T)))
            names(dist_per_bas) <- c("basin_id", "distance")
            
            ### Rescale
            dist_per_bas$distance_scaled <-  rescale(dist_per_bas$distance, to = c(1, 0)) 
            dist_per_bas$distance_scaled <- ifelse(is.na(dist_per_bas$distance_scaled), -1, dist_per_bas$distance_scaled) # -1 if non-connected 
            
            ### Get as raster
            basin_table <- as.data.table(basin_id, na.rm=F) # get basins in a table
            basin_table$seq_id <- seq.int(1:nrow(basin_table))
            setkey(basin_table, seq_id)
            basin_for_dist <- merge(basin_table, dist_per_bas, by="basin_id", all.x=T)
            basin_for_dist <- as.data.table(arrange(basin_for_dist, seq_id))
            basin_for_dist_r <- setValues(basin_id, basin_for_dist$distance_scaled)
            writeRaster(basin_for_dist_r,  paste0(DIR, "/", MY_SPP, "/basin_distance.tif"), overwrite=T)
            
          }
          
          
          ###-------------------#
          ###--- Run models ----
          ###-------------------#
          
          ### Initiate
          myBiomodData <- BIOMOD_FormatingData(resp.var = SPP_NA[,i,with=F] ,
                                               expl.var = EXP,
                                               resp.xy =  COORD, # NULL
                                               resp.name = MY_SPP,
                                               PA.nb.rep = 1,
                                               PA.nb.absences = 10000, 
                                               PA.strategy = 'random')
          
          
          ### Options definition
          # myBiomodOption <- BIOMOD_ModelingOptions() # default options
          myBiomodOption <- BIOMOD_ModelingOptions(
            
            RF = list( do.classif = TRUE,
                       ntree = 1000, # default 500, models more stable with 1000 trees?
                       mtry = 'default',
                       nodesize = 5,
                       maxnodes = NULL),
            # NOTE: Maxent needs x11 connection even without GUI...
            MAXENT.Phillips = list( 
              # path_to_maxent.jar = "./maxent/maxent.jar",
              path_to_maxent.jar = MAXENT_DIR_TMP, # MAXENT_DIR
              # path_to_maxent.jar = "C:\temp\spdep_SDMs\MAXENT_DIR", # in windows
              memory_allocated = 4096, # give more RAM
              background_data_dir = 'default',
              maximumbackground = 'default',
              maximumiterations = 200,
              visible = FALSE,
              linear = TRUE,
              quadratic = TRUE,
              product = TRUE,
              threshold = FALSE, # edited
              hinge = FALSE, # edited
              lq2lqptthreshold = 80,
              l2lqthreshold = 10,
              hingethreshold = 15,
              beta_threshold = -1,
              beta_categorical = -1,
              beta_lqp = -1,
              beta_hinge = -1,
              betamultiplier = 1,
              defaultprevalence = 0.5)
          )
          
          ### Modelling
          myBiomodModelOut <- BIOMOD_Modeling(
            myBiomodData,
            # models = c('SRE','CTA','RF','MARS','FDA','MAXENT.Phillips','GLM','GAM','GBM','ANN'), # ALL ALGORITHMS
            # models = c('CTA','RF','FDA','MAXENT.Phillips','GLM','GAM','GBM','ANN'), # without Mars, SRE
            # models = c('SRE','CTA','RF', 'FDA','MAXENT.Phillips', 'GLM'), # fast ones
            models = c('CTA','RF', 'FDA','MAXENT.Phillips', 'GLM'), # fast ones
            # models = c('RF', 'MAXENT.Phillips'), #  Machine learning subset
            # models = c('RF','CTA', 'FDA', 'MAXENT.Phillips'), #  Machine learning and classification subset
            # models = c('RF', 'MAXENT.Phillips'),
            # models = c('MAXENT.Phillips'), 
            models.options = myBiomodOption,
            NbRunEval=10, # repetitions: use 2 for testing, 10 for the final run
            DataSplit=70,
            Yweights=NULL,
            VarImport=3,
            models.eval.meth = c('TSS','ROC'),
            SaveObj = TRUE,
            rescal.all.models = TRUE)
          
          gc()
          
          ### Ensemble model
          myBiomodEM <- BIOMOD_EnsembleModeling(
            modeling.output = myBiomodModelOut,
            chosen.models = 'all',  # only final model?
            em.by = 'all',
            eval.metric = c('TSS'),
            eval.metric.quality.threshold = c(0.4), # e.g. 0.4
            prob.mean = F,
            prob.cv = F,
            prob.ci = F,
            prob.ci.alpha = 0.05,
            prob.median = F,
            committee.averaging = F,
            prob.mean.weight = T,
            prob.mean.weight.decay = "proportional" ) # 1.6
          
          ###---- Project on present-day data
          myBiomodProj <- BIOMOD_Projection(
            modeling.output = myBiomodModelOut,
            new.env = EXP,
            proj.name = 'present',
            selected.models = 'all',
            binary.meth= 'TSS', # NULL
            compress = 'xz',
            clamping.mask = F)
          
          ### Do ensemble-models projections on present variables
          # ?BIOMOD_EnsembleForecasting
          myBiomodEF <- BIOMOD_EnsembleForecasting(
            projection.output = myBiomodProj,
            EM.output = myBiomodEM,
            binary.meth = 'TSS', # NULL
            total.consensus = TRUE)
          
          gc() 
          ### Save the ensemble evaluation metrics separately:
          x <- get_evaluations(myBiomodEM)
          write.csv(x, paste0(DIR, "/", MY_SPP, "/ensemble_evaluation.csv")); rm(x)
          
          ### Load predictions - probability
          load(paste0(DIR, "/", MY_SPP, "/proj_present/proj_present_", MY_SPP, "_ensemble.RData"))
          tmp_proj_present <- as.data.table(as.data.frame(ef.out)); rm(ef.out)
          names(tmp_proj_present) <- "probs"
          tmp_proj_present <- tmp_proj_present/1000
          tmp_proj_present <- cbind(all_data[,c("basin_id")], tmp_proj_present)
          
          
          ### Prepare map
          tmp_proj_present_domain <- merge(domain_cells, tmp_proj_present, by="basin_id", all.x=T)
          ### Sort data to match the spatial configuration
          tmp_proj_present_domain <- as.data.table(arrange(tmp_proj_present_domain, seq_id))
          ### Insert the values into the raster template
          pred_ensemble_r <- setValues(basin_id, tmp_proj_present_domain$probs)
          names(pred_ensemble_r) <- "pred_ensemble"
          
          
          ### Get species points
          # mypoints <- as.data.table(as.data.frame(subset(sp_points, species==MY_SPP)))
          
          ### Save map as pdf
          pdf(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_probability.pdf"), width=20, height=20)
          plot(pred_ensemble_r, col=viridis(15), zlim=c(0,1))
          points(subset(sp_points, Species==MY_SPP), pch=16, col="red")
          dev.off()
          
          
          ### Load predictions - binary TSS
          load(paste0(DIR, "/", MY_SPP, "/proj_present/proj_present_", MY_SPP, "_ensemble_TSSbin.RData"))
          ef.out_bin <- get(paste0("proj_present_",MY_SPP, "_ensemble_TSSbin"))
          tmp_proj_present_bin <- as.data.table(as.data.frame(ef.out_bin)); rm(ef.out_bin)
          names(tmp_proj_present_bin) <- "binary"
          tmp_proj_present_bin <- tmp_proj_present_bin
          tmp_proj_present_bin <- cbind(all_data[,c("basin_id")], tmp_proj_present_bin)
          
          ### Prepare map
          tmp_proj_present_bin_domain <- merge(domain_cells, tmp_proj_present_bin, by="basin_id", all.x=T)
          ### Sort data to match the spatial configuration
          tmp_proj_present_bin_domain <- as.data.table(arrange(tmp_proj_present_bin_domain, seq_id))
          ### Insert the values into the raster template
          pred_ensemble_bin_r <- setValues(basin_id, tmp_proj_present_bin_domain$binary)
          names(pred_ensemble_bin_r) <- "pred_ensemble_binary"
          
          ### Save map as pdf
          pdf(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_binary.pdf"), width=20, height=20)
          plot(pred_ensemble_bin_r, col=viridis(15), zlim=c(0,1))
          points(subset(sp_points, Species==MY_SPP), pch=16, col="red")
          dev.off()
          
          ### Get variable importance
          var_imp <- as.data.frame(variables_importance(myBiomodEM@em.models[[1]], data=EXP , method="full_rand" , nb_rand=3)$mat)
          var_imp$mean <- rowMeans(var_imp)
          tmp <- subset(var_imp, select=c(mean))
          var_imp$percent_contribution <- round(apply(tmp, 1, function(x) x/colSums(tmp) ),2)  # scale from 0-1
          var_imp$variable <- row.names(var_imp)
          write.csv(var_imp, paste0(DIR, "/", MY_SPP, "_variable_importance.csv"), row.names = F)
          
          ### Write raster to disk
          writeRaster(pred_ensemble_r, paste0(DIR, "/", MY_SPP, "_ensemble.tif"), overwrite=T)
          writeRaster(pred_ensemble_bin_r, paste0(DIR, "/", MY_SPP, "_ensemble_bin.tif"), overwrite=T)
          
          
          gc()
          
          ###--- Check single algorithm results ----
          myCurrentProj <- get_predictions(myBiomodProj)
          ### Get single "full" model maps
          single_full_pred <- as.data.table(as.data.frame(myCurrentProj[,,"Full","PA1"]))
          single_full_pred <- single_full_pred/1000 # rescale to 0-1
          single_full_pred <- cbind(all_data[,c("basin_id")], single_full_pred)
          
          ### Prepare map
          single_full_pred_domain <- merge(domain_cells, single_full_pred, by="basin_id", all.x=T)
          ### Sort data to match the spatial configuration
          single_full_pred_domain <- as.data.table(arrange(single_full_pred_domain, seq_id))
          
          
          ### Create character vector of all possible algorithms
          check_these <- c("SRE", "CTA","RF","MARS","GLM","GAM","GBM","ANN","MAXENT.Phillips","MAXENT.Phillips.2")
          mystack <- stack() # empty stack for storing layers
          
          ### Run through all columns, check if model exists and write into raster stack
          for (k in check_these) {
            if(k %in% colnames(single_full_pred_domain)) { 
              tmp <- setValues(basin_id, single_full_pred_domain[,get(k)])
              names(tmp) <- k 
              mystack <- addLayer(mystack, tmp) 
            }
          }
          
          ### Save single algorighm map as pdf
          pdf(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_single_algorithms.pdf"), width=20, height=20)
          plot(mystack, col=viridis(15), zlim=c(0,1))
          dev.off()
          
          
          ### Export the raw, non-scaled and non-normalized variable importances of each algorthm
          ### careful interpretation is "the higher the value, the more influence it has"
          var_imp_single <- as.data.frame(myBiomodModelOut@variables.importances@val[,,"Full", "PA1"])
          var_imp_single$variable <- row.names(var_imp_single)
          write.csv(var_imp_single, paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_single_algorithm_raw_variable_importance.csv"), row.names = F )
          
          ### Clean up
          rm(mystack, single_full_pred, single_full_pred_domain, var_imp_single)
          gc()
          
          
          ###--- Multiply the predictions with the distance map (clipping)----
          if(DO_DISTANCE==TRUE) {
            pred_ensemble_r_dist <- pred_ensemble_r * basin_for_dist_r
            pred_ensemble_r_dist[pred_ensemble_r_dist<=0] <- 0
            ### Save map as pdf
            pdf(paste0(DIR, "/",MY_SPP, "/", MY_SPP, "_ensemble_distance.pdf"), width=20, height=20)
            plot(pred_ensemble_r_dist, col=viridis(15), zlim=c(0,1))
            points(subset(sp_points, Species==MY_SPP), pch=16, col="red")
            dev.off()
            
            writeRaster(pred_ensemble_r_dist, paste0(MY_SPP_DIR, "/", MY_SPP, "_ensemble_distance.tif"), overwrite=T)
            rm(pred_ensemble_r_dist) # clean up
          }
          
          ### Save all predictions as a table - useful for e.g. Marxan later on..
          all_predictions <- merge(tmp_proj_present, tmp_proj_present_bin, by="basin_id")
          save(all_predictions, file=paste0(MY_SPP_DIR, "/", MY_SPP, "_all_predictions_basin_table.RData")); rm(all_predictions)
          ## if distance_weights needed: load the separate RData file and attach to this table in another loop...
          
          
          ### Remove files for next species run
          rm(pred_ensemble_r, 
             pred_ensemble_bin_r, 
             tmp_proj_present, 
             tmp_proj_present_bin, 
             myBiomodEF, 
             myBiomodProj, 
             myBiomodEM, 
             myBiomodModelOut, 
             myBiomodData, 
             MY_SPP)
          
          system(paste0("rm -rf ", MAXENT_DIR_TMP))   # remove maxent folder
          
          
          
          gc() 
          
        } # end loop

stopCluster(cl) # close the parrallel backend
