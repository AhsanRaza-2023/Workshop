
#### Use this code for creating loops and most impotently split data based on St_ID for train and test data fro models 

#######################################################################################################################
#######################################################################################################################
#######################################################################################################################

base_path <- "C:\\Users\\raza\\Desktop\\rfsi_workshop\\RFSI-master\\RFSI-master\\brandenburg"

RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(42)

wd=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

utm33 <- "+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
#utm33 <- "+proj=longlat +datum=WGS84"
border_wgs84 <- readOGR(paste0(base_path,"\\borders\\osm\\Brandenburg1.shp"))
border <- spTransform(border_wgs84, CRS(utm33))

dir.create("plot")
dir.create("temp_data")
setwd(paste(wd, "/temp_data/", sep = ""))

v = "GWD"
year="2022-22"
covariates <- c("GWD", "Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")
# Loop through years from 2001 to 2022

for (year in "2016-22") {
  file_name <- paste0(base_path, paste0("\\temp_data\\clean_", year, ".csv"))
  bran <- read.csv(file_name)
  bran <- bran[,-1]
  bran <- na.omit(bran)
  
  # Process the data as before
  # bran$DATE <- as.Date(paste0(bran$DATE, "-01"), format = "%Y-%m-%d")    to del
  bran$DATE <- as.Date(paste0(bran$DATE, "-01"), format = "%b-%y-%d")
  
  ob <- bran[, c(5, 6, 10, 18, 19, 20, 21, 22)]
  sta <- bran[, c(5, 3, 4, 7, 8, 11, 12, 13, 14, 15, 16, 17)]
  sta <- sta[!duplicated(sta), ]
  colnames(ob)
  colnames(sta)
  # Create STFDF
  stfdf <- meteo2STFDF(obs = ob,
                       stations = sta,
                       crs = CRS(utm33),
                       obs.staid.time = c(1, 2),
                       stations.staid.lon.lat = c(1, 4, 5)
  )
  
  stfdf@sp@data$staid <- 1:nrow(stfdf@sp)
  
  # Save STFDF as an R data file
  #save(stfdf, file = paste0("stfdf_", year, ".rda"))
}

###### stations plot ######

load("stfdf_2016-22.rda")
#load("grids.rda")
stfdf@sp@data
ss <- as.data.frame(stfdf)

stfdf@sp <- spTransform(stfdf@sp, CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

r <- raster("DEM.tif")

# Transform the raster to WGS84
r <- projectRaster(r, crs = CRS("+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))

# Convert 'r' to SpatialPixelsDataFrame
r <- as(r, "SpatialPixelsDataFrame")

# Create the ggplot plot
sta_dem_plot <- ggplot() +
  geom_raster(data = as.data.frame(r), aes(x = x, y = y, fill = DEM), alpha = 0.8) +
  geom_point(data = as.data.frame(stfdf@sp), aes(x = lon, y = lat), size = 0.8, alpha = 0.5) +
  scale_fill_gradientn(colours = terrain.colors(100), name = "DEM [m]") +
  scale_shape_manual(name = "Stations", values = c("staid" = 17)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()+ theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank()
  )

sta_dem_plot

# tiff("../plot/stations.tiff", width = 100, height = 62, units = 'mm', res = 600, compression = "lzw")
jpeg("../plot/stations.jpeg", width = 150, height = 112, units = 'mm', res = 1800)
sta_dem_plot
dev.off()


his <- ggplot(stfdf@data, aes(x = GWD)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(
  x = "GWD", y = "Frequency") +
  scale_x_continuous(breaks = seq(0, max(stfdf@data$GWD, na.rm = TRUE), by = 5)) +  
  scale_y_continuous(limits = c(0, 30000), breaks = seq(0, 30000, by = 5000)) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18)
  )

# tiff("../plot/histogram.tiff", width = 80, height = 75, units = 'mm', res = 600, compression = "lzw")
ggsave(filename = "../plot/histogram.jpeg", plot = his, width = 25, height = 20, units = "cm", dpi = 100)

############################################
#############################################################  RF 
##########################################

library(doParallel)
library(ranger)
library(caret)

cores <- detectCores()
registerDoParallel(cores = cores)

# Loop through years from 2001 to 2022
year_ranges <- "2016-22"

# Initialize variables for tracking lowest RMSE and iterations since improvement
lowest_rmse <- Inf
iterations_since_improvement <- 0
max_iterations_without_improvement <- 20

{
  
for (year in year_ranges) {
  start_time <- Sys.time()
  cat("Processing year:", year, "\n")
  file_name <- paste0("stfdf_", year, ".rda")
  
  # Load spatiotemporal data frame for the current year
  load(file = file_name)
  
  # Convert the data to a format suitable for modeling
  temp_df <- as.data.frame(stfdf)
  covariates <- c("GWD", "Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")
  temp_df <- temp_df[, c("lon", "lat", "sp.ID", "time", covariates)]
  temp_df = temp_df[complete.cases(temp_df), ]
  
  time = sort(unique(temp_df$time))
  daysNum = length(time)
  days = gsub("-", "", time, fixed = TRUE)
  
  min.node.size <- 2:10
  sample.fraction <- seq(1, 0.832, -0.025) # 0.632 without / 1 with replacement
  ntree <- 250 # 500
  mtry <- 1:9
  
  indices <- CreateSpacetimeFolds(temp_df, spacevar = "sp.ID", k = 5, seed = 42)
  
  hp <- expand.grid(min.node.size = min.node.size, mtry = mtry, sf = sample.fraction)
  hp <- hp[sample(nrow(hp), 50), ]
  hp <- hp[order(hp$mtry), ]
  rmse_hp <- rep(NA, nrow(hp))
  
  for (h in 1:nrow(hp)) {
    comb <- hp[h, ]
    print(paste("combination: ", h, sep = ""))
    print(comb)
    
    fold_obs <- c()
    fold_pred <- c()
    
    for (f in 1:length(indices$index)) {
      print(paste("fold: ", f, sep = ""))
      dev_df1 <- temp_df[indices$index[[f]], covariates]
      val_df1 <- temp_df[indices$indexOut[[f]], covariates]
      
      model <- ranger(GWD ~ ., data = dev_df1, importance = "none", seed = 42,
                      num.trees = ntree, mtry = comb$mtry,
                      splitrule = "extratrees",  #variance  extratrees
                      min.node.size = comb$min.node.size,
                      sample.fraction = comb$sf,
                      oob.error = FALSE)
      fold_obs <- c(fold_obs, val_df1$GWD)
      fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
      
    }
    
    rmse_hp[h] <- round(sqrt(mean((fold_obs - fold_pred)^2, na.rm = TRUE)), 3)
    
    print(paste("rmse: ", rmse_hp[h], sep = ""))
    
    # Check if the current RMSE is lower than the lowest RMSE
    if (rmse_hp[h] < lowest_rmse) {
      lowest_rmse <- rmse_hp[h]
      iterations_since_improvement <- 0
    } else {
      iterations_since_improvement <- iterations_since_improvement + 1
    }
    
    # Check if we've reached the maximum iterations without improvement
    if (iterations_since_improvement >= max_iterations_without_improvement) {
      cat("No improvement in RMSE for", max_iterations_without_improvement, "iterations. Ending the process.\n")
      break
    }
    
  }
  
  if (iterations_since_improvement >= max_iterations_without_improvement) {
    cat("Exiting early due to lack of improvement.\n")
    break
  }
 }


  dev_parameters <- hp[which.min(rmse_hp), ]
  print(dev_parameters)
  
  temp_df <- temp_df[, covariates]
  colnames(temp_df)
  
  # Identify unique combinations of latitude and longitude
  unique_coords <- unique(temp_df[, c("Lat", "Lon")])
  
  # Set a seed for reproducibility
  set.seed(123)
  
  # Create an index for splitting the data
  index <- createDataPartition(1:nrow(unique_coords), p = 0.8, list = FALSE)
  
  # Use the index to split the unique coordinates into training and testing sets
  train_coords <- unique_coords[index, ]
  test_coords <- unique_coords[-index, ]
  
  # Merge the original data with the training and testing coordinates
  train_data <- merge(temp_df, train_coords, by = c("Lat", "Lon"))
  test_data <- merge(temp_df, test_coords, by = c("Lat", "Lon"))
  
  # Select covariates
  covariates <- colnames(temp_df)
  covariates <- covariates[covariates != "GWD"]  # Assuming "GWD" is the response variable
  
  # Subset the data frames with selected covariates
  train_data <- train_data[, c("GWD", covariates)]
  test_data <- test_data[, c("GWD", covariates)]
  
  rf_model <- ranger(GWD ~ ., data = train_data, importance = "impurity", seed = 42,
                     num.trees = ntree, mtry = dev_parameters$mtry,
                     splitrule = "extratrees",  #variance
                     min.node.size = dev_parameters$min.node.size,
                     sample.fraction = dev_parameters$sf,
                     quantreg = TRUE)
  
  # Save the trained model for the current year
  save(rf_model, file = paste0("../models/RF_workshop_1", year, ".rda"))
  
  # Create and save plot_data
  plot_data <- data.frame(Lat = train_data$Lat, Lon = train_data$Lon, Observations = train_data$GWD, Predictions = rf_model$predictions)
  plot_file <- paste0(base_path, (paste0("\\temp_data\\rf_training_workshop_1", year, ".csv")))
  #write.csv(plot_data, plot_file)
  
  # Make predictions on the test data using the trained model
  predictions <- predict(rf_model, data = test_data, type = "response")
  
  # Add the predictions to the test_data dataframe
  test_data$predicted_GWD <- predictions$predictions
  
  plot_data_test <- data.frame(Lat = test_data$Lat, Lon = test_data$Lon, GWD = test_data$GWD, predicted_GWD = test_data$predicted_GWD)
  
  # Save test_data with predictions as a CSV file
  test_file <- paste0(base_path, (paste0("\\temp_data\\rf_test__workshop_1", year, ".csv")))
  #write.csv(plot_data_test, test_file)
  
  end_time <- Sys.time()
  computation_time <- end_time - start_time
  cat("Computation time for year", year, ":", computation_time, "\n")
  
  # Save dev_parameters, RMSE, and computation_time to a text file
  result_file <-  paste0(base_path,"\\computation_times_rf_workshop1.txt")
  write(paste("Year:", year), file = result_file, append = TRUE)
  write(paste("dev_parameters:", toString(dev_parameters)), file = result_file, append = TRUE)
  write(paste("RMSE:", min(rmse_hp, na.rm = TRUE)), file = result_file, append = TRUE)
  write(paste("Computation time:", computation_time), file = result_file, append = TRUE)
  write("\n", file = result_file, append = TRUE)
}
 


############################### create raster RF
# Define directories
model_dir <- paste0(base_path, "\\models\\")
test_data_dir <- paste0(base_path, "\\temp_data\\b_data_for_prediction\\")
output_dir <- paste0(base_path, "\\plot\\del\\")
reference_raster_path <- paste0(base_path, "\\temp_data\\DEM.tif")

# Load the reference raster
reference_raster <- raster(reference_raster_path)

# Define the new CRS
new_crs <- CRS("+proj=longlat +datum=WGS84")

# Project the reference raster
projected_raster <- projectRaster(reference_raster, crs = new_crs)

# Create an empty raster to hold the results
new_raster <- raster(projected_raster)

# Load the model for the current year
model_file <- paste0(model_dir, "RF_2016-22.rda")
load(file = model_file)


# Loop over each year from 2001 to 2022
for (year in 2016:2022) {
  
  
  # Load the test data for the current year
  test_data_file <- paste0(test_data_dir, "data_for_pred_1", year, "_v2.csv")
  test_data <- read.csv(test_data_file)
  test_data <- na.omit(test_data)
  
  # Perform predictions
  predictions <- predict(rf_model, data = test_data)
  
  # Add predictions to the test data
  test_data$Predictions <- predictions$predictions
  
  # Convert the 'Date' column to a Date object if it's not already
  test_data$Date <- as.Date(test_data$Date, format = "%Y-%m-%d")
  
  # List to store the rasters
  monthly_rasters <- list()
  
  # Iterate through each month in the current year
  unique_months <- unique(format(test_data$Date, "%Y-%m"))
  
  for (month in unique_months) {
    # Filter data for the current month
    monthly_data <- test_data %>%
      filter(format(test_data$Date, "%Y-%m") == month)
    
    # Create a data frame for the current month
    data <- data.frame(Longitude = monthly_data$Lon, Latitude = monthly_data$Lat, Predictions = monthly_data$Predictions)
    coordinates(data) <- c("Longitude", "Latitude")
    
    # Rasterize the data
    monthly_raster <- rasterize(data, new_raster, "Predictions")
    
    # Define the output file path
    output_file <- paste0(output_dir, "monthly_rf_1", month, ".tif")
    
    # Save the raster as a GeoTIFF
    writeRaster(monthly_raster, filename = output_file, format = "GTiff", overwrite = TRUE)
  }
  
  print(paste("Completed processing for year:", year))
}



######################################## RFSP
########################################

library(devtools)
#install_github("envirometrix/landmap")
library(landmap)


# Loop through years from 2001 to 2022
#year_ranges <- c("2001-05", "2006-10", "2011-15", "2016-22")
year_ranges <- c("2016-22")
# Loop through years from 2001 to 2022
for (year in year_ranges) {
  # Load spatiotemporal data frame for the current year
  start_time <- Sys.time()
  cat("Processing year:", year, "\n")
  load(file = paste0("stfdf_", year, ".rda"))
  
  covariates <- c("GWD", "Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")
  
  ### Create BBOX ###
  xmin <- min(stfdf@sp@coords[, "lon"]) - 1000
  xmax <- max(stfdf@sp@coords[, "lon"]) + 1000
  ymin <- min(stfdf@sp@coords[, "lat"]) - 1000
  ymax <- max(stfdf@sp@coords[, "lat"]) + 1000
  bbox <- extent(xmin, xmax, ymin, ymax)
  bbox <- as(bbox, "SpatialPolygons")
  bbox@proj4string <- CRS(utm33)  #"+proj=longlat +datum=WGS84"
  fishnet <- raster(bbox)
  res(fishnet) <- c(1000, 1000)  # 1 x 1 km
  fishnet[] <- runif(length(fishnet), -10, 10)
  fishnet <- as(fishnet, "SpatialPixelsDataFrame")
  st <- as(fishnet, "SpatialPoints")
  plot(fishnet)
  #plot(border, add = TRUE)
  extent(fishnet)
  extent(stfdf@sp["staid"])
  class(fishnet)
  
  ## Geographic distances only:
  grid.distP <- landmap::buffer.dist(stfdf@sp["staid"], fishnet[1], as.factor(1:nrow(stfdf@sp)))
  #print(grid.distP)
  
  temp_df <- as.data.frame(stfdf)
  
  temp_df$staid <- as.integer(as.character(temp_df$staid))
  ov_toma <- do.call(cbind, list(stfdf@sp@data, over(stfdf@sp, grid.distP)))
  temp_df <- plyr::join(temp_df, ov_toma, by = "staid")
  
  temp_df <- temp_df[, c(covariates, paste("layer.", sort(unique(as.integer(as.character(temp_df$staid)))), sep = ""))]
  
  temp_df <- temp_df[complete.cases(temp_df), ]
  #colnames(temp_df)
  
  time <- sort(unique(temp_df$time))
  daysNum <- length(time)
  days <- gsub("-", "", time, fixed = TRUE)
  nrow(temp_df)
  
  ntree <- 250  # 500   
  
  
  # tune RFsp as Hengl et al. 2018 did
  rt <- makeRegrTask(data = temp_df[sample(nrow(temp_df), 12000), ], target = "GWD")
  
  estimateTimeTuneRanger(rt, num.threads = detectCores())
  
  # Time-consuming >> do on number cruncher
  set.seed(42)
  t.rfsp <- tuneRanger(rt, num.trees = ntree, iters = 120, build.final.model = FALSE)
  
  dev_parameters <- t.rfsp$recommended.pars
  print(dev_parameters)
  
  # Identify unique combinations of latitude and longitude
  unique_coords <- unique(temp_df[, c("Lat", "Lon")])
  
  # Set a seed for reproducibility
  set.seed(123)
  
  # Create an index for splitting the data
  index <- createDataPartition(1:nrow(unique_coords), p = 0.8, list = FALSE)
  
  # Use the index to split the unique coordinates into training and testing sets
  train_coords <- unique_coords[index, ]
  test_coords <- unique_coords[-index, ]
  
  # Merge the original data with the training and testing coordinates
  train_data <- merge(temp_df, train_coords, by = c("Lat", "Lon"))
  test_data <- merge(temp_df, test_coords, by = c("Lat", "Lon"))
  
  # Select covariates
  covariates <- colnames(temp_df)
  covariates <- covariates[covariates != "GWD"]  # Assuming "GWD" is the response variable
  
  # Subset the data frames with selected covariates
  train_data <- train_data[, c("GWD", covariates)]
  test_data <- test_data[, c("GWD", covariates)]
  
  rfsp_model <- ranger(GWD ~ ., data = train_data, importance = "impurity", seed = 42,
                       num.trees = ntree, mtry = dev_parameters$mtry,
                       splitrule = "variance",
                       min.node.size = dev_parameters$min.node.size,
                       sample.fraction = dev_parameters$sample.fraction,
                       quantreg = TRUE)  ### quantreg???
  
  # Save the trained model for the current year
  save(rfsp_model, file = paste0("../models/RFsp_workshop_1", year, ".rda"))
  
  # Create and save plot_data
  plot_data <- data.frame(Lat = train_data$Lat, Lon = train_data$Lon, Observations = train_data$GWD, Predictions = rfsp_model$predictions)
  plot_file <- paste0(base_path, (paste0("\\temp_data\\rfsp_training_workshop_1", year, ".csv")))
  write.csv(plot_data, plot_file)
  
  # Make predictions on the test data using the trained model
  predictions <- predict(rfsp_model, data = test_data, type = "response")
  
  # Add the predictions to the test_data dataframe
  test_data$predicted_GWD <- predictions$predictions
  
  test_data <- test_data[, c("GWD", "predicted_GWD", "Lat", "Lon")]
  
  # Save test_data with predictions as a CSV file
  test_file <- paste0(base_path, (paste0("\\temp_data\\rfsp_test_workshop_1", year, ".csv")))
  write.csv(test_data, test_file)
  
  end_time <- Sys.time()
  computation_time <- end_time - start_time
  #cat("dev_parameters for year", year, ":", dev_parameters, "\n")
  cat("Computation time for year", year, ":", computation_time, "\n")
  
  # Save dev_parameters, RMSE, and computation_time to a text file
  result_file <- paste0(base_path, "\\temp_data\\computation_times_rfsp_workshop1.txt")
  write(paste("Year:", year), file = result_file, append = TRUE)
  write(paste("dev_parameters:", toString(dev_parameters)), file = result_file, append = TRUE)
  #write(paste("RMSE:", min(rmse_hp, na.rm = TRUE)), file = result_file, append = TRUE)
  write(paste("Computation time:", computation_time), file = result_file, append = TRUE)
  write("\n", file = result_file, append = TRUE)
  
}


###############################  to create raster with rfsp
# Load necessary libraries
# Load necessary libraries
library(doParallel)
library(landmap)
library(raster)
library(data.table)  # For efficient data manipulation
library(dplyr)       # For data manipulation
library(plyr)        # For joining data

# Define directories
model_dir <- paste0(base_path, "\\models\\")
test_data_dir <- paste0(base_path, "\\temp_data\\b_data_for_prediction\\")
output_dir <- paste0(base_path, "\\plot\\del\\")
reference_raster_path <- paste0(base_path, "\\temp_data\\DEM.tif")

# Load the model for the current year
model_file <- paste0(model_dir, "RFsp_2016-22.rda")
load(model_file)

# Load the reference raster and project it
reference_raster <- raster(reference_raster_path)
new_crs <- CRS("+proj=longlat +datum=WGS84")
projected_raster <- projectRaster(reference_raster, crs = new_crs)
new_raster <- raster(projected_raster)

# Loop over each year from 2001 to 2022
for (year in 2016:2022) {
  
  # Load the test data for the current year
  test_data_file <- paste0(test_data_dir, "data_for_pred_1", year, "_v2.csv")
  test_data <- fread(test_data_file) %>% na.omit()
  
  # Add 'staid' column to test data
  test_data[, staid := .I]
  
  # Define bounding box and create fishnet
  bbox <- extent(min(test_data$Lon), max(test_data$Lon), min(test_data$Lat), max(test_data$Lat))
  bbox <- as(bbox, "SpatialPolygons")
  proj4string(bbox) <- new_crs
  fishnet <- raster(bbox)
  res(fishnet) <- c(0.01, 0.01)
  fishnet[] <- runif(ncell(fishnet), -10, 10)
  fishnet <- as(fishnet, "SpatialPixelsDataFrame")
  
  # Convert test data to spatial points
  spatial_points_df <- SpatialPointsDataFrame(coords = test_data[, .(Lon, Lat)], data = test_data, proj4string = new_crs)
  
  # Process each unique date
  unique_dates <- unique(test_data$Date)
  final_results <- list()  # To collect results
  
  for (date in unique_dates) {
    # Filter data for the current date
    date_data <- test_data[Date == date]
    
    # Extract coordinates for the current date
    coords_date <- date_data[, .(Lon, Lat)]
    
    # Reconstruct SpatialPointsDataFrame for the current date
    spatial_points_df_date <- SpatialPointsDataFrame(coords = coords_date, data = date_data, proj4string = new_crs)
    
    # Compute geographic distances using landmap for the current date
    grid.distP1 <- landmap::buffer.dist(spatial_points_df_date["staid"], fishnet[1], as.factor(1:nrow(spatial_points_df_date)))
    
    # Join data based on 'staid'
    ov_toma1 <- cbind(date_data, over(spatial_points_df_date, grid.distP1))
    ov_toma1 <- plyr::join(date_data, ov_toma1, by = "staid")
    
    # Reorder and clean columns
    desired_order <- c("Date", "Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")
    rest_of_columns <- setdiff(names(ov_toma1), desired_order)
    ov_toma1 <- ov_toma1[, c(desired_order, rest_of_columns)]
    ov_toma1 <- ov_toma1[, !names(ov_toma1) %in% c("staid", "Lon.1", "Lat.1"), with = FALSE]
    ov_toma1 <- na.omit(ov_toma1)
    
    
    # Make predictions for the current date data
    predictions <- predict(rfsp_model, data = ov_toma1)
    ov_toma1$Predictions <- predictions$predictions
    
    # Collect results for the current date
    final_results[[date]] <- ov_toma1
  }
  
  # Combine all date-specific results
  combined_results <- rbindlist(final_results)
  
  # Save the predictions to CSV for the entire year
  output_csv <- paste0(output_dir, year, "_pred_rfsp1.csv")
  #write(combined_results[, .(Date, Lat, Lon, Predictions)], output_csv, row.names = FALSE)
  
  # Generate rasters by month
  combined_results[, Date := as.Date(Date, format = "%Y-%m-%d")]
  unique_months <- unique(format(combined_results$Date, "%Y-%m"))
  
  for (month in unique_months) {
    # Filter data for the current month
    monthly_data <- combined_results[format(Date, "%Y-%m") == month]
    
    # Create a data frame for the current month
    data <- data.frame(Longitude = monthly_data$Lon, Latitude = monthly_data$Lat, Predictions = monthly_data$Predictions)
    coordinates(data) <- c("Longitude", "Latitude")
    
    # Rasterize the data
    monthly_raster <- rasterize(data, new_raster, "Predictions")
    
    # Define the output file path
    output_file <- paste0(output_dir, "monthly_rfsp_1", month, ".tif")
    
    # Save the raster as a GeoTIFF
    writeRaster(monthly_raster, filename = output_file, format = "GTiff", overwrite = TRUE)
  }
  
  print(paste("Completed processing for year:", year))
}


###### create RFSI model ###################################################################  (Use this one)
#################################### limit the number of iterations 

year_ranges <- c("2016-22")

# Loop through years from 2001 to 2022
for (year in year_ranges) {
  start_time <- Sys.time()
  cat("Processing year:", year, "\n")
  file_name <- paste0("stfdf_", year, ".rda")
  
  # Load spatiotemporal data frame for the current year
  load(file = file_name)
  
  # Convert the data to a format suitable for modeling
  temp_df <- as.data.frame(stfdf)
  
  covariates <- c("GWD", "Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")
  temp_df <- temp_df[, c("Lon", "Lat", "sp.ID", "time", covariates)]
  temp_df = temp_df[complete.cases(temp_df), ]
  
  time = sort(unique(temp_df$time))
  daysNum = length(time)
  days = gsub("-", "", time, fixed=TRUE)
  
  nos <- 6:14
  min.node.size <- 5:6
  #min.node.size <- 4:6 # 2:6
  sample.fraction <- seq(1, 0.632, -0.05)
  ntree <- 250
  
  
  n_obs = max(nos)
  #mtry <- 3:(9 + 2 * n_obs)
  mtry <- 3
  
  indices <- CreateSpacetimeFolds(temp_df, spacevar = "sp.ID", k = 5, seed = 42)
  
  hp <- expand.grid(min.node.size = min.node.size, mtry = mtry, no = nos, sf = sample.fraction)
  hp <- hp[hp$mtry < (9 + 2 * hp$no - 1), ]
  hp <- hp[sample(nrow(hp), 100), ]
  hp <- hp[order(hp$no), ]
  rmse_hp <- rep(NA, nrow(hp))
  
  lowest_rmse <- Inf  # Initialize with a high value
  iteration_since_lowest <- 0  # Counter for iterations since the lowest RMSE
  
  for (h in 1:nrow(hp)) {
    comb <- hp[h, ]
    print(paste("combination: ", h, sep = ""))
    print(comb)
    
    dev_df <- temp_df
    
    fold_obs <- c()
    fold_pred <- c()
    
    for (f in 1:length(indices$index)) {
      print(paste("fold: ", f, sep = ""))
      cols <- c(covariates, paste("dist", 1:(comb$no), sep = ""), paste("obs", 1:(comb$no), sep = ""))
      dev_df1 <- dev_df[indices$index[[f]], ]
      
      cpus <- detectCores() - 1
      registerDoParallel(cores = cpus)
      nearest_obs <- foreach (t = time) %dopar% {
        dev_day_df <- dev_df1[as.numeric(dev_df1$time) == t, c("Lon", "Lat", "GWD")]
        day_df <- dev_df[as.numeric(dev_df$time) == as.numeric(t), c("Lon", "Lat", "GWD")]
        
        if (nrow(day_df) == 0) {
          return(NULL)
        }
        return(meteo::near.obs(
          locations = day_df,
          observations = dev_day_df,
          obs.col = "GWD",
          n.obs = comb$no
        ))
      }
      
      stopImplicitCluster()
      
      nearest_obs <- do.call("rbind", nearest_obs)
      dev_df <- cbind(dev_df, nearest_obs)
      dev_df = dev_df[complete.cases(dev_df), ]
      
      dev_df1 <- dev_df[indices$index[[f]], cols]
      dev_df1 <- dev_df1[complete.cases(dev_df1), ]
      val_df1 <- dev_df[indices$indexOut[[f]], cols]
      val_df1 <- val_df1[complete.cases(val_df1), ]
      
      model <- ranger(GWD ~ ., data = dev_df1, importance = "none", seed = 42,
                      num.trees = ntree, mtry = comb$mtry,
                      splitrule = "variance",  # variance extratrees
                      min.node.size = comb$min.node.size,
                      sample.fraction = comb$sf,
                      oob.error = FALSE)
      fold_obs <- c(fold_obs, val_df1$GWD)
      fold_pred <- c(fold_pred, predict(model, val_df1)$predictions)
    }
    
    rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm = TRUE))
    print(paste("rmse: ", rmse_hp[h], sep = ""))
    
    if (rmse_hp[h] < lowest_rmse) {
      lowest_rmse <- rmse_hp[h]
      iteration_since_lowest <- 0  # Reset counter
      print(paste("New lowest RMSE:", lowest_rmse))
    } else {
      iteration_since_lowest <- iteration_since_lowest + 1
      print(paste("Iterations since lowest RMSE:", iteration_since_lowest))
    }
    
    # Stop if 5 iterations have passed since the lowest RMSE was found
    if (iteration_since_lowest >= 5) {
      print(paste("Stopping after", iteration_since_lowest, "iterations since the lowest RMSE"))
      break
    }
  }
  
  dev_parameters <- hp[which.min(rmse_hp), ]
  print(dev_parameters)
  
  cpus <- detectCores() - 1
  registerDoParallel(cores = cpus)
  
  nearest_obs <- foreach (t = time) %dopar% {
    day_df <- temp_df[temp_df$time == t, c("Lon", "Lat", "GWD")]
    if (nrow(day_df) == 0) {
      return(NULL)
    }
    return(near.obs(
      locations = day_df,
      observations = day_df,
      obs.col = "GWD",
      n.obs = dev_parameters$no
    ))
  }
  
  stopImplicitCluster()
  
  nearest_obs <- do.call("rbind", nearest_obs)
  temp_df <- cbind(temp_df, nearest_obs)
  
  cols <- c(covariates, paste("dist", 1:(dev_parameters$no), sep = ""), paste("obs", 1:(dev_parameters$no), sep = ""))
  temp_df <- temp_df[, cols]
  colnames(temp_df)
  a <- temp_df
  # Identify unique combinations of latitude and longitude
  unique_coords <- unique(temp_df[, c("Lat", "Lon")])
  
  # Set a seed for reproducibility
  set.seed(123)
  
  # Create an index for splitting the data
  index <- createDataPartition(1:nrow(unique_coords), p = 0.8, list = FALSE)
  
  # Use the index to split the unique coordinates into training and testing sets
  train_coords <- unique_coords[index, ]
  test_coords <- unique_coords[-index, ]
  
  # Merge the original data with the training and testing coordinates
  train_data <- merge(temp_df, train_coords, by = c("Lat", "Lon"))
  test_data <- merge(temp_df, test_coords, by = c("Lat", "Lon"))
  nrow(test_data)
  write_csv(test_data, "test_data_lat_lon.csv")
  # Select covariates
  covariates <- colnames(temp_df)
  covariates <- covariates[covariates != "GWD"]  # Assuming "GWD" is the response variable
  
  # Subset the data frames with selected covariates
  train_data <- train_data[, c("GWD", covariates)]
  test_data <- test_data[, c("GWD", covariates)]
  
  rfsi_model <- ranger(GWD ~ ., data = train_data, importance = "impurity", seed = 42,
                       num.trees = ntree, mtry = dev_parameters$mtry,
                       splitrule = "extratrees",
                       min.node.size = dev_parameters$min.node.size,
                       sample.fraction = dev_parameters$sf,
                       quantreg = TRUE)  ### quantreg???
  
  # save(rfsi_model, file = "../models/RFSI.rda")
  save(rfsi_model, file = paste0("../models/RFSI_workshop_1", year, ".rda"))
  
  # Create and save plot_data
  plot_data <- data.frame(Lat = train_data$Lat, Lon = train_data$Lon, Observations = train_data$GWD, Predictions = rfsi_model$predictions)
  plot_file <-  paste0(base_path, (paste0("\\temp_data\\rfsi_training_workshop_1", year, ".csv")))
  write.csv(plot_data, plot_file)
  
  # Make predictions on the test data using the trained model
  predictions <- predict(rfsi_model, data = test_data, type = "response")
  
  # Add the predictions to the test_data dataframe
  test_data$predicted_GWD <- predictions$predictions
  
  #test_data <- test_data[, c("GWD", "predicted_GWD")]
  
  plot_data_test <- data.frame(Lat = test_data$Lat, Lon = test_data$Lon, GWD = test_data$GWD, predicted_GWD = test_data$predicted_GWD)
  # Save test_data with predictions as a CSV file
  test_file <-  paste0(base_path, (paste0("\\temp_data\\rfsi_test_workshop_1", year, ".csv")))
  write.csv(plot_data_test, test_file)
  
  end_time <- Sys.time()
  computation_time <- end_time - start_time
  #cat("dev_parameters for year", year, ":", dev_parameters, "\n")
  cat("Computation time for year", year, ":", computation_time, "\n")
  
  # Save dev_parameters, RMSE, and computation_time to a text file
  result_file <- paste0(base_path, "\\temp_data\\computation_times_rfsi_workshop1.txt")
  write(paste("Year:", year), file = result_file, append = TRUE)
  write(paste("dev_parameters:", toString(dev_parameters)), file = result_file, append = TRUE)
  write(paste("RMSE:", min(rmse_hp, na.rm = TRUE)), file = result_file, append = TRUE)
  write(paste("Computation time:", computation_time), file = result_file, append = TRUE)
  write("\n", file = result_file, append = TRUE)
}


###########################################
############### Create predictions wtih rfsi
################################

# Load the necessary data and model
load("stfdf_2016-22.rda")
load(file = "../models/RFSI_2016-22.rda")
rfsi_model
covariates <- c("GWD", "Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")

# Extract the data into a dataframe and clean it
temp_df <- as.data.frame(stfdf)
temp_df <- temp_df[, c("sp.ID", "time", covariates)]  # Assuming `covariates` is predefined
temp_df <- temp_df[complete.cases(temp_df), ]

# Extract unique times and sort them
time <- sort(unique(temp_df$time))
daysNum <- length(time)
days <- gsub("-", "", time, fixed=TRUE)

# Convert data to a spatial dataframe for prediction
newdata <- as.data.frame(temp_df)
newdata <- na.omit(newdata)

# Specify the years of interest
years_of_interest <- c(2016:2021)  # Corrected the year range

# Iterate through each year
for (year in years_of_interest) {
  
  # Load the test data for the current year
  test_data_file <- paste0(base_path, sprintf("\\temp_data\\b_data_for_prediction\\data_for_pred_%s_v2.csv", year))
  data <- read.csv(test_data_file)
  
  # Drop rows with NA values
  data <- na.omit(data)
  
  # Convert Lon and Lat to SpatialPointsDataFrame
  spatial_points <- SpatialPoints(coords = data[c("Lon", "Lat")])
  spatial_points_df <- SpatialPointsDataFrame(coords = spatial_points, data = data.frame(data))
  proj4string(spatial_points_df) <- CRS("+proj=longlat +datum=WGS84")
  
  # Reorder columns in desired order
  desired_order <- c("Date","Lat", "Lon", "SLOPE", "ASP", "DEM", "LNF", "TPI", "TWI", "TRI", "PREC", "TEMP", "NDVI", "ET", "EVI")
  rest_of_columns <- setdiff(names(data), desired_order)
  ov_toma1 <- data[, c(desired_order, rest_of_columns)]
  
  # Remove unnecessary columns
  columns_to_delete <- c("staid", "Lon.1", "Lat.1", "X", "Y", "optional")
  ov_toma1 <- ov_toma1[, !names(ov_toma1) %in% columns_to_delete]
  
  # Remove rows with NA values
  ov_toma1 <- ov_toma1[complete.cases(ov_toma1), ]
  
  # Convert the data to a spatial object
  loc_sf <- st_as_sf(ov_toma1, coords = c("Lon", "Lat"), crs = 4326)
  newdata_sf <- st_as_sf(newdata, coords = c("Lon", "Lat"), crs = 4326)
  newdata_sf$time <- as.Date(newdata_sf$time, format = "%Y-%m-%d")
  
  # Initialize a list to store nearest observations
  nearest_obs_list <- list()
  
  # Calculate nearest observations for each date
  unique_dates <- unique(loc_sf$Date)
  
  loc_sf$Date <- as.Date(loc_sf$Date)
  newdata_sf$time <- as.Date(newdata_sf$time)
  
  for (current_date in unique_dates) {
    # Filter data for the current date
    loc_day <- loc_sf[loc_sf$Date == current_date, ]
     
    newdata_day <- newdata_sf[newdata_sf$time == current_date, ]
    
    # Calculate nearest observations
    nearest_obs <- meteo::near.obs(
      locations = loc_day,
      observations = newdata_day,
      obs.col = "GWD",  # Assuming GWD is the observation column
      n.obs = 100  # Adjust the number of nearest neighbors if needed
    )
    
    # Store the nearest observations for the current date
    nearest_obs_list[[as.character(current_date)]] <- nearest_obs
  }
  
  # Combine nearest observations into a single dataframe
  nearest_obs_combined <- do.call(rbind, nearest_obs_list)
  nearest_obs_combined$date <- rownames(nearest_obs_combined)
  
  # Remove the row names as they are now stored in the "date" column
  rownames(nearest_obs_combined) <- NULL
  
  # Optionally move the "date" column to the first position
  nearest_obs_combined <- nearest_obs_combined[, c("date", setdiff(names(nearest_obs_combined), "date"))]
  nearest_obs_combined$date <- gsub("\\.\\d+$", "", nearest_obs_combined$date)
  
  newdata1 <- cbind(loc_sf, nearest_obs_combined)
  
  newdata_sf <- st_as_sf(newdata1, coords = "geometry", crs = 4326)
  
  coordinates <- st_coordinates(newdata_sf)
  newdata1 <- cbind(coordinates, newdata1 %>% select(-geometry))
  colnames(newdata1)[1:2] <- c("Lon", "Lat")
  
  newdata1 <- as.data.frame(newdata1)
  columns_to_remove <- c("Date", "geometry")
  
  # Remove specified columns
  newdata1 <- newdata1[, !names(newdata1) %in% columns_to_remove]
  
  
  # Perform prediction using the RFSI model
  pre <- predict(rfsi_model, newdata1)
  Predictions <- as.data.frame(pre)
  
  # Add predictions as a new column to the original data
  loc_sf$Predictions <- Predictions$prediction
  
  loc_sf11 <- st_as_sf(loc_sf, coords = "geometry", crs = 4326)
  
  coordinates <- st_coordinates(loc_sf11)
  loc_sf11 <- cbind(coordinates, loc_sf11 %>% select(-geometry))
  colnames(loc_sf11)[1:2] <- c("Lon", "Lat")
  
  loc_sf11 <- as.data.frame(loc_sf11)
  columns_to_remove <- c("geometry")
  loc_sf11 <- loc_sf11[, !names(loc_sf11) %in% columns_to_remove]
  
  # Save predictions to CSV file
  output_file <- paste0(base_path, sprintf("\\temp_data\\%s_pred_rfsi1.csv", year))
  #write.csv(loc_sf, file = output_file, row.names = FALSE)
  
  # Load and project the reference raster
  reference_raster <- paste0(base_path, raster("\\temp_data\\DEM.tif"))
  new_crs <- CRS("+proj=longlat +datum=WGS84")
  projected_raster <- projectRaster(reference_raster, crs = new_crs)
  
  # Create and save raster files for each month
  loc_sf11$Date <- as.Date(loc_sf11$Date, format = "%Y-%m-%d")
  unique_months <- unique(format(loc_sf11$Date, "%Y-%m"))
  
  output_dir <- paste0(base_path, "\\plot\\del\\")
  
  for (month in unique_months) {
    monthly_data <- loc_sf11 %>% filter(format(Date, "%Y-%m") == month)
    
    # Rasterize the data
    data <- data.frame(Longitude = monthly_data$Lon, Latitude = monthly_data$Lat, Predictions = monthly_data$Predictions)
    coordinates(data) <- c("Longitude", "Latitude")
    
    monthly_raster <- rasterize(data, projected_raster, field = "Predictions")
    
    # Define the output file path
    output_file <- file.path(output_dir, paste0("monthly_rfsi_1", month, ".tif"))
    
    # Save the raster as a GeoTIFF
    writeRaster(monthly_raster, filename = output_file, format = "GTiff", overwrite = TRUE)
  }
}


