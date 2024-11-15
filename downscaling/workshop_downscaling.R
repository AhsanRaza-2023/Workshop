
####################################### RF for variable importance for complete period 

base_path <- "C:\\Users\\raza\\Desktop\\rfsi_workshop\\downscaling"

csv_dir <- paste0(base_path, "\\var_imp")  #When NDVI and EVI is mean
# List all CSV files in the directory
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
print(csv_files)

# Read, remove first three columns, and combine CSV files
combined_data <- do.call(rbind, lapply(csv_files, function(file) {
  data <- read.csv(file)
  data <- data[, -c(1)]  # Remove first three columns
  return(data)
}))

# Read and combine CSV files
#combined_data <- do.call(rbind, lapply(csv_files, read.csv))
colnames(combined_data)
# Train random forest model
mboston.rf <- randomForest(GWS ~ ., data = combined_data[,-c(1,2)], ntree = 500, importance = TRUE)

# Visualize variable importance
importance <- varImpPlot(mboston.rf)

# Create a data frame for plotting
importance_df <- data.frame(
  Variables = rownames(importance),
  Importance = importance[, "%IncMSE"]
)

# Min-max scaling function
min_max_scale <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

# Scale the variable importance values
importance_df$Importance <- min_max_scale(importance_df$Importance)
# Sort the data frame by importance (optional)
importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
importance_df$Importance <- importance_df$Importance * 100
# Create a ggplot bar plot for variable importance
imp <- ggplot(importance_df, aes(x = reorder(Variables, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  xlab("Variables") +
  ylab("Importance") + 
  theme_bw() +
  ggtitle("Variable Importance Plot") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1, size = 25),
    axis.text.y = element_text(angle = 0, hjust = 1, size = 25), # Increase text size
    axis.title = element_text(size = 27), # Increase axis label size
    plot.title = element_text(size = 27) # Increase title size
  ) +
  coord_flip()

imp

ggsave(paste0(base_path, "\\var_imp\\importance_plot_all.png"), plot = imp, width = 8, height = 10)

##################### Create partial dependence plots 

# First, create a list of all partial dependence data frames
predictors <- c("TEMP", "ET", "LST", "NDVI", "EVI", "PREC", "LNF", "DEM", "SLOPE", "TPI", "ASP")

# Generate all partial dependence data in a loop and combine into one dataframe
pd_data <- lapply(predictors, function(pred) {
  pd <- partial(mboston.rf, pred.var = pred, plot = FALSE)
  pd$Variable <- pred  # Add variable name as a column
  names(pd)[1] <- "Value"  # Rename first column to generic name
  return(pd)
}) %>% bind_rows()

# Create a single plot with facets
ab <- ggplot(pd_data, aes(x = Value, y = yhat)) +
  geom_line(linewidth = 0.5, color = "steelblue") +
  geom_point(size = 0.5, color = "steelblue") +
  facet_wrap(~factor(Variable, c("TEMP", "ET", "LST", "NDVI", "EVI", "PREC", "DEM", "LNF", "SLOPE", "TPI", "ASP")), scales = "free_x", ncol = 4)+
  theme_bw() +
  theme(
    strip.text = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "black"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 10, face = "bold", hjust = 0.5)
  ) +
  labs(
    x = "Predictor Value",
    y = "GWSA (mm)"
  ) +
  coord_cartesian(ylim = c(-250, 250)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", alpha = 0.5)

ab

ggsave(paste0(base_path, "\\var_imp\\pds_contr_all.png"), plot = ab, width = 8, height = 7)


######################################################################################## DOWNASCALING WITH MGWR
########################################################################################

years <- 2015:2020
months <- 01:12
# Initialize the list to store year-month combinations
year_month_combinations <- list()

for (year in years) {
  for (month in months) {
    year_month_combinations <- append(year_month_combinations, list(list(year = as.character(year), month = sprintf("%02d", month))))
  }
}

output_folder <- paste0(base_path,"\\mgwr")


# Iterate over year and month combinations
for (info in year_month_combinations) {
  folder_path <- file.path(base_path,
    info$year,
    info$month
  )
  
  setwd(folder_path)
  
  # Create a raster list of all .tif files in the folder
  rlist <- list.files(pattern = "\\.tif$", full.names = TRUE)
  
  # Stack the raster layers (note all rasters should have same extent and projection)
  rasters_1km <- stack(rlist)
  
  # Load Basin boundary shapefile
  Boundry_indus <- readOGR(paste0(base_path,"\\Germany_shapefile\\DEU_adm0.shp"))
  Boundry_indus <- spTransform(x = Boundry_indus, CRS("+proj=longlat +datum=WGS84"))
  extent(Boundry_indus)
  plot(Boundry_indus)
  
  variable_name <- paste0(info$year, "_", info$month)
  data_file <- file.path(folder_path, paste0(variable_name, ".csv"))
  
  # Check if the CSV file exists
  if (!file.exists(data_file)) {
    message("CSV file missing for: ", info$year, "-", info$month, ". Skipping...")
    next
  }
  
  Data <- read_csv(data_file)
  Longitude <- Data$Longitude
  Latitude <- Data$Latitude
  # Now start to calibrate GWR model on coarse resolution data
  GWRbandwidth <- gwr.sel(GWSA ~ SLOPE+DEM+Prec+Temp+NDVI+LST+	EVI+ET, data=Data, coords=cbind(Longitude, Latitude),adapt=T) 
  #GWRbandwidth <- gwr.sel(GWS ~ LNF+TPI+ASP+SLOPE+DEM+Prec+Temp+NDVI+LST+	EVI+ET, data=Data, coords=cbind(Longitude, Latitude),adapt=T) 
  #run the gwr model
  gwr.model = gwr(GWSA ~ SLOPE+DEM+Prec+Temp+NDVI+LST+	EVI+ET, data=Data, coords=cbind(Latitude, Longitude), adapt=GWRbandwidth, hatmatrix=TRUE, se.fit=TRUE) 
  #print the results of the model
  gwr.model
  results<-as.data.frame(gwr.model$SDF)
  head(results)
  results$GWSA <- Data$GWSA
  
  coordinates(results) = ~Longitude+Latitude 
  proj4string(results) <- CRS("+proj=longlat +datum=WGS84")
  #plot(results, add = TRUE)
  # Define extent of interpolation and grid size (Note: grid size should be same as input data of coarse resolution)
  x.range <- as.numeric(c(5.871619, 15.03811 )) # min/max longitude of the interpolation area
  y.range <- as.numeric(c(47.26986 , 56.05653 )) # min/max latitude of the interpolation area
  grd <- expand.grid(X = seq(from = x.range[1], to = x.range[2], by = 1), Y = seq(from = y.range[1], to = y.range[2], by = 1)) # expand points to grid
  
  coordinates(grd) <- ~ X+Y
  proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")
  gridded(grd) <- TRUE
  #plot(grd, cex=1.5)
  points(results, pch=1, col='red', cex=1)
  
  #Let's creat refrence grid cells for high reesolution prediction (Note: these grid cells should have same resolution at which we are going to downscale our predictant variable)
  grd_resamp <- expand.grid(X = seq(from = x.range[1], to = x.range[2], by = 0.01), Y = seq(from = y.range[1], to = y.range[2], by = 0.01)) # expand points to grid
  coordinates(grd_resamp) <- ~ X+Y
  proj4string(grd_resamp) <- CRS("+proj=longlat +datum=WGS84")
  gridded(grd_resamp) <- TRUE
  r_resample <- raster(grd_resamp)
  
  i = "sum.w"
  as.formula(paste0(i, " ~ 1"))
  #sum.w ~ 1
  timeL<- colnames(data.frame(results))
  m <- c(timeL)
  result = list()
  
  for (i in m) {
    tryCatch({
      f = as.formula(paste0(i, " ~ 1"))
      suppressWarnings({
        variogcloud = variogram(f, locations = results, data = results, cloud = TRUE)
        semivariog = variogram(f, locations = results, data = results)
        model.variog <- vgm(psill = 200, model = "Exp", nugget = 2000, range = 500) #Gau
        fit.variog = fit.variogram(semivariog, model.variog)
        krig = gstat::krige(formula = f, locations = results, newdata = grd, model = model.variog)
      })
      z_int_0.25 = raster(krig)
      #plot(z_int_0.25)
      r_sample_0.01 <- raster::resample(z_int_0.25, r_resample, method = "bilinear")
      masked <- mask(x = r_sample_0.01, mask = Boundry_indus)
      result[[i]] = masked
    }, error = function(e) {
      # Handle or ignore the error (optional)
      # print(paste("Error for variable", i, ":"))
      # print(e)
    })
  }
  
  # Stack all rasters coefficients
  Coeff_stack = stack(result)
  
  variable_name <- paste0(info$year, ".", info$month)
  et_layer_name <- paste0("ET_mean.", variable_name)
  ET <- raster::resample((crop(extend(rasters_1km[[et_layer_name]], Coeff_stack$ET), Coeff_stack$LST)), Coeff_stack$ET, method = "bilinear")
  
  lst_layer_name <- paste0("LST_mean.", variable_name)
  LST <- raster::resample((crop(extend(rasters_1km[[lst_layer_name]], Coeff_stack$LST), Coeff_stack$LST)),Coeff_stack$LST, method ="bilinear")
  
  temp_layer_name <- paste0("Temp_mean.", variable_name)
  TEMP <- raster::resample((crop(extend(rasters_1km[[temp_layer_name]], Coeff_stack$Temp), Coeff_stack$Temp)),Coeff_stack$Temp, method ="bilinear")
  
  evi_layer_name <- paste0("EVI_mean.", variable_name)
  EVI <- raster::resample((crop(extend(rasters_1km[[evi_layer_name]], Coeff_stack$EVI), Coeff_stack$EVI)),Coeff_stack$EVI, method ="bilinear")
  
  prec_layer_name <- paste0("Prec_sum.", variable_name)
  PREC <- raster::resample((crop(extend(rasters_1km[[prec_layer_name]], Coeff_stack$Prec), Coeff_stack$Prec)),Coeff_stack$Prec, method ="bilinear")
  
  ndvi_layer_name <- paste0("NDVI_mean.", variable_name)
  NDVI <- raster::resample((crop(extend(rasters_1km[[ndvi_layer_name]], Coeff_stack$NDVI), Coeff_stack$NDVI)),Coeff_stack$NDVI, method ="bilinear")
  NDVI <- NDVI / 10000
  
  DEM <- raster::resample((crop(extend(rasters_1km$Germany_1000m_dem, Coeff_stack$DEM), Coeff_stack$DEM)),Coeff_stack$DEM, method ="bilinear")
  SLOPE <- raster::resample((crop(extend(rasters_1km$Slope_Germnay_1000m, Coeff_stack$SLOPE), Coeff_stack$SLOPE)),Coeff_stack$SLOPE, method ="bilinear")
  
  # Calculate GWR using the new variable name
  GWR <- Coeff_stack$X.Intercept. + (Coeff_stack$ET * ET) + (Coeff_stack$EVI * EVI) + (Coeff_stack$DEM * DEM) + 
    (Coeff_stack$LST * LST) + (Coeff_stack$NDVI * NDVI) + (Coeff_stack$Prec * PREC) + 
    (Coeff_stack$SLOPE * SLOPE) + (Coeff_stack$Temp * TEMP) + (Coeff_stack$pred - Coeff_stack$GWSA)
  
  # Save raster output with the new variable name
  output_filename <- file.path(output_folder, paste0("GWR_", info$year, "-", info$month, ".tiff"))
  writeRaster(GWR, output_filename, overwrite = TRUE)
  
}


#############################################  RF 2
#############################################
#############################################

library(Metrics) 
library(raster)
library(sp)
library(caret)
library(dplyr)
library(readr)
library(rgdal) # For readOGR function
library(spacetime) # For spatiotemporal analysis
library(ranger)

# Initialize an empty data frame to combine all CSVs
combined_data <- data.frame()

# Loop through years and months
for (year in 2015:2020) {
  for (month in sprintf("%02d", 1:12)) {
    folder_path <- file.path(base_path, year, month)
    
    # Check if the folder exists
    if (dir.exists(folder_path)) {
      # Define the CSV file name
      data_file <- file.path(folder_path, paste0(year, "_", month, ".csv"))
      
      # Check if the file exists
      if (file.exists(data_file)) {
        # Read the CSV file
        df <- read_csv(data_file)
        
        # Add the year column to the data frame
        df$Year <- year
        
        # Bind the current data frame to the combined data frame
        combined_data <- bind_rows(combined_data, df)
        combined_data <- combined_data[!is.na(combined_data$Year), ]
      }
    }
  }
}


combined_data1 <- combined_data

combined_data <- combined_data %>%
  group_by(Latitude, Longitude) %>%
  mutate(sp.ID = as.integer(factor(paste(Latitude, Longitude)))) %>%
  ungroup()

combined_data <- as.data.frame(combined_data)
combined_data <- combined_data %>%
  distinct()

# Set a seed for reproducibility
set.seed(42)

covariates <- c("GWSA", "SLOPE", "DEM", "Prec", "Temp", "NDVI", "LST", "EVI", "ET")
combined_data <- combined_data[, c("Longitude", "Latitude", "sp.ID", covariates)]

combined_data = combined_data[complete.cases(combined_data), ]
min.node.size <- 2:10
sample.fraction <- seq(0.832, 0.663, -0.025) # 0.632 without / 1 with replacement
ntree <- 250 # 500
mtry <- 8

indices <- CreateSpacetimeFolds(combined_data, spacevar = "sp.ID", k=5, seed = 42)

hp <- expand.grid(min.node.size=min.node.size, mtry=mtry, sf=sample.fraction)
hp <- hp[sample(nrow(hp), 50),]
hp <- hp[order(hp$mtry),]
rmse_hp <- rep(NA, nrow(hp))

# Initialize variables to store best model data and RMSE
best_rmse <- Inf
best_dev_df1 <- NULL
best_val_df1 <- NULL
best_dev_preds <- NULL

# Initialize a counter for iterations after the best RMSE is found
no_improvement_counter <- 0
max_no_improvement <- 1  # Number of iterations to check for improvement

# Loop through each hyperparameter setting
for (h in 1:nrow(hp)) {
  fold_obs <- NULL
  fold_pred <- NULL
  
  # Loop through each fold
  for (f in 1:length(indices$index)) {
    print(paste("fold: ", f, sep=""))
    
    # Define development and validation datasets for the current fold
    dev_df1 <- combined_data[indices$index[[f]], covariates]
    val_df1 <- combined_data[indices$indexOut[[f]], covariates]
    
    # Train the model on the development dataset
    model <- ranger(GWSA ~ SLOPE + DEM + Prec + Temp + NDVI + LST + EVI + ET, 
                    data = dev_df1, importance = "none", seed = 42,
                    num.trees = ntree, mtry = hp$mtry[h],
                    splitrule = "extratrees",
                    min.node.size = hp$min.node.size[h],
                    sample.fraction = hp$sf[h],
                    oob.error = FALSE)
    
    # Collect predictions and actual values for validation data
    val_predictions <- predict(model, val_df1)$predictions
    fold_obs <- c(fold_obs, val_df1$GWSA)
    fold_pred <- c(fold_pred, val_predictions)
  }
  
  # Calculate RMSE for the current hyperparameter set
  rmse_hp[h] <- sqrt(mean((fold_obs - fold_pred)^2, na.rm = TRUE))
  rmse_hp[h] <- round(rmse_hp[h], 2)  # Round RMSE to 2 decimal places
  print(paste("rmse: ", rmse_hp[h], sep=""))
  
  # Print the hyperparameter combination used in this iteration
  cat(sprintf("Hyperparameters used: min.node.size = %d, mtry = %d, sample.fraction = %.2f\n", 
              hp$min.node.size[h], hp$mtry[h], hp$sf[h]))
  
  # Check for improvement
  if (rmse_hp[h] < best_rmse) {
    best_rmse <- rmse_hp[h]
    best_dev_df1 <- dev_df1
    best_val_df1 <- val_df1
    best_dev_preds <- predict(model, best_dev_df1)$predictions  # Store predictions for dev_df1
    
    # Reset the no improvement counter since we found a new best RMSE
    no_improvement_counter <- 0
  } else {
    # Increment the no improvement counter
    no_improvement_counter <- no_improvement_counter + 1
  }
  
  # Stop iteration if no improvement has occurred for max_no_improvement iterations
  if (no_improvement_counter >= max_no_improvement) {
    print("Stopping iteration: no improvement in RMSE over the last 10 iterations.")
    break
  }
}

# Add predictions to the best development and validation data frames
best_dev_df1$pred <- best_dev_preds
best_val_df1$pred <- predict(model, best_val_df1)$predictions


# Ensure both data frames have the required columns
# best_dev_df1 <- best_dev_df1[, c("Latitude", "Longitude", "Year", "Month", "GWSA", "pred")]
# best_val_df1 <- best_val_df1[, c("Latitude", "Longitude", "Year", "Month", "GWSA", "pred")]

best_dev_df1 <- best_dev_df1[, c("GWSA", "pred")]
best_val_df1 <- best_val_df1[, c("GWSA", "pred")]

# # Save the best development and validation data frames to CSV
#write.csv(best_dev_df1, paste0(base_path, "\\best_dev_df1.csv"), row.names = FALSE)
#write.csv(best_val_df1, paste0(base_path, "\\best_val_df1.csv"), row.names = FALSE)

# Print the best parameters
dev_parameters  <- hp[which.min(rmse_hp), ]
print(dev_parameters )


model <- lm(pred ~ GWSA, data = best_val_df1)
r2 <- summary(model)$r.squared

mae <- mae(best_val_df1$GWSA, best_val_df1$pred)
rmse <- rmse(best_val_df1$GWSA, best_val_df1$pred)

# Create the 1:1 plot with R², MAE, and RMSE values displayed
ggplot(best_val_df1, aes(x = GWSA, y = pred)) +
  geom_point(color = 'dodgerblue', alpha = 0.1, size = 1) +  # Scatter plot with adjusted color and alpha
  geom_smooth(method = "lm", color = "black", linetype = "solid", size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed', size = 0.5) +  # 1:1 line with thicker line
  labs(
       x = "Actual GWSA (mm)",
       y = "Predicted GWSA (mm)") +
  xlim(-150, 100) +  # Set x-axis limits
  ylim(-150, 100) +  # Set y-axis limits
  theme_bw() +  # Minimal theme with larger base font size
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", color = "darkblue"),  # Center-align title with bold style
    panel.grid.major = element_line(color = "gray80"),  # Lighter grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines for cleaner look
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 12, face = "bold")
  ) +
  annotate("text", x = -120, y = 90, label = paste("R²:", round(r2, 2)), color = "black", size = 5, fontface = "bold") +
  annotate("text", x = -120, y = 70, label = paste("MAE:", round(mae, 2), "mm"), color = "black", size = 5, fontface = "bold") +
  annotate("text", x = -120, y = 50, label = paste("RMSE:", round(rmse, 2), "mm"), color = "black", size = 5, fontface = "bold") +
  scale_color_manual(values = c("Actual" = "dodgerblue", "1:1 Line" = "red"))

ggsave(paste0(base_path,"\\test_Actual_vs_Predicted_Plot.jpeg"), width = 8, height = 5, dpi = 200)


combined_data <- combined_data[, covariates]

rf_model <- ranger(GWSA ~ ., data = combined_data, importance = "impurity", seed = 42,
                   num.trees = ntree, mtry = dev_parameters$mtry,
                   splitrule = "extratrees",
                   min.node.size = dev_parameters$min.node.size,
                   sample.fraction = dev_parameters$sf,
                   quantreg = TRUE) ### quantreg???


save(rf_model, file = paste0(base_path,"\\RF.rda"))
load(file = paste0(base_path,"\\RF.rda"))

# # Define the predictor columns
# predictor_columns <- c("SLOPE", "DEM", "Prec", "Temp", "NDVI", "LST", "EVI", "ET")
# 
# # Select only predictor columns for prediction
# predict_data <- combined_data[, predictor_columns]

# Run prediction using the rf_model on predictor columns only
predictions <- rf_model$predictions


# Combine predictions with actual values and extract year and month
predictions_df <- data.frame(
  Longitude = combined_data1$Longitude,
  Latitude = combined_data1$Latitude,
  Actual = combined_data1$GWSA,
  Predicted = predictions,
  Residual = combined_data1$GWSA - predictions,  # Calculate residuals directly
  Year = combined_data1$Year,      # Extract year
  Month = combined_data1$Month     # Extract month
)

results <- as.data.frame(predictions_df)

model <- lm(Predicted ~ Actual, data = results)
r2 <- summary(model)$r.squared

# Create the 1:1 plot with R² value displayed
ggplot(results, aes(x = Actual, y = Predicted)) +
  geom_point(color = 'blue', alpha = 0.1) +  # Scatter plot
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = 'dashed') +  # 1:1 line
  labs(title = "Actual vs Predicted Values",
       x = "Actual Values",
       y = "Predicted Values") +
  xlim(-150, 100) +  # Set x-axis limits
  ylim(-150, 100) +  # Set y-axis limits
  theme_bw() +  # Minimal theme
  annotate("text", x = -100, y = 90, label = paste("R² =", round(r2, 2)), color = "black", size = 5)

# # Load the reference raster
reference_raster <- raster(paste0(base_path,"\\GWSA_monthly_2020-01.tif"))
Boundary_indus <- readOGR(paste0(base_path,"\\Germany_shapefile\\DEU_adm0.shp"))

# Transform the shapefile to match the raster CRS
Boundary_indus <- spTransform(Boundary_indus, crs(reference_raster))

# Mask the raster with the shapefile
reference_raster <- mask(reference_raster, Boundary_indus)

# Create a function to create and save raster for each unique year and month
create_and_save_raster <- function(predictions_df, reference_raster) {
  unique_year_month <- unique(predictions_df[, c("Year", "Month")])
  
  for (i in 1:nrow(unique_year_month)) {
    year <- unique_year_month$Year[i]
    month <- unique_year_month$Month[i]
    
    # Create an empty raster with the same extent and resolution as the reference raster
    new_raster <- reference_raster
    
    # Fill the raster with predictions
    for (j in 1:nrow(predictions_df)) {
      if (predictions_df$Year[j] == year && predictions_df$Month[j] == month) {
        coords <- c(predictions_df$Longitude[j], predictions_df$Latitude[j])
        
        # Convert coordinates to raster cell index
        cell_index <- cellFromXY(new_raster, coords)
        
        # Assign the predicted value to the corresponding cell
        new_raster[cell_index] <- predictions_df$Predicted[j]
      }
    }
    
    # Define the output filename
    output_filename <- paste0("Predicted_GWSA_", year, "-", sprintf("%02d", month), ".tif")
    output_path <- file.path(paste0(base_path, "\\predictions"), output_filename)
    
    # Save the new raster
    writeRaster(new_raster, filename = output_path, format = "GTiff", overwrite = TRUE)
  }
}

# Call the function with your predictions dataframe
create_and_save_raster(predictions_df, reference_raster)


# Create a SpatialPointsDataFrame for predictions
coordinates(predictions_df) <- ~ Longitude + Latitude
proj4string(predictions_df) <- CRS("+proj=longlat +datum=WGS84")
library(sf)
# Load the basin boundary shapefile
Boundary_indus <- readOGR(paste0(base_path,"\\Germany_shapefile\\DEU_adm0.shp"))
Boundary_indus <- spTransform(x = Boundary_indus, CRS("+proj=longlat +datum=WGS84"))
plot(Boundary_indus)

# Create a reference grid
x.range <- as.numeric(c(5.871619, 15.03811)) # min/max longitude of the interpolation area
y.range <- as.numeric(c(47.26986, 56.05653)) # min/max latitude of the interpolation area
grd <- expand.grid(X = seq(from = x.range[1], to = x.range[2], by = 1), 
                   Y = seq(from = y.range[1], to = y.range[2], by = 1))
coordinates(grd) <- ~ X + Y
proj4string(grd) <- CRS("+proj=longlat +datum=WGS84")
gridded(grd) <- TRUE

# High-resolution grid for resampling
grd_resamp <- expand.grid(X = seq(from = x.range[1], to = x.range[2], by = 0.01), 
                          Y = seq(from = y.range[1], to = y.range[2], by = 0.01))
coordinates(grd_resamp) <- ~ X + Y
proj4string(grd_resamp) <- CRS("+proj=longlat +datum=WGS84")
gridded(grd_resamp) <- TRUE
r_resample <- raster(grd_resamp)

# Loop through each unique Year and Month combination
unique_years_months <- unique(results[, c("Year", "Month")])
for (i in 1:nrow(unique_years_months)) {
   year <- unique_years_months[i, "Year"]
   month <- unique_years_months[i, "Month"]
  
  # Filter results for the specific year and month
  filtered_results <- results %>%
    filter(Year == year & Month == month)
  
  # Convert filtered results to SpatialPointsDataFrame
  sp_results <- SpatialPointsDataFrame(
    coords = filtered_results[, c("Longitude", "Latitude")], 
    data = filtered_results,
    proj4string = CRS("+proj=longlat +datum=WGS84")
  )
  
  # Fit the variogram model for residuals
  residual_variog <- variogram(Residual ~ 1, locations = sp_results)
  residual_model <- vgm(psill = 1, model = "Sph", nugget = 2.5, range = 250)
  fit_residual_variog <- fit.variogram(residual_variog, residual_model)
  
  
  #fit_residual_variog <- fit.variogram(residual_variog, residual_model)
  plot(residual_variog, fit_residual_variog)
  
  # Perform kriging for residuals
  krig_residual <- krige(Residual ~ 1, sp_results, newdata = grd, model = fit_residual_variog)
  z_residual_int <- raster(krig_residual)
  
  # Resample kriged residuals to fine resolution
  r_residual_sample <- raster::resample(z_residual_int, r_resample, method = "bilinear")
  masked_residual  <- mask(r_residual_sample, mask = Boundary_indus)
  plot(masked_residual)
  # Save or plot the masked residual as needed
  output_filename <- file.path(paste0(base_path,"\\residuals"), 
                               paste0("Residuals-", year, "-", sprintf("%02d", as.numeric(month)), ".tiff"))
  writeRaster(masked_residual, output_filename, overwrite = TRUE)
  
  # Optionally plot the results for visualization
  plot(masked_residual, main = paste("Masked Residuals for", year, month))
}



# Define folders for raster input, residuals, and output
raster_folder <- paste0(base_path,"\\all_raster_1km")
residuals_folder <- paste0(base_path,"\\residuals")
output_folder <- paste0(base_path,"\\rf")

# Loop through each year and month combination
for (year in 2015:2020) {
  for (month in sprintf("%02d", 01:12)) {
    variable_name <- paste0(year, "-", month)
    
    # Load the residual raster for the current year and month
    residual_file <- file.path(residuals_folder, paste0("Residuals-", year, "-", month, ".tiff"))
    
    masked_residual <- raster(residual_file)
    
    
    et_layer_name <- paste0("ET_mean-", variable_name)
    et_file_path <- file.path(raster_folder, paste0(et_layer_name, ".tif"))
    ET <- raster(et_file_path)
    ET <- crop(extend(ET, masked_residual), masked_residual)
    
    # Explicitly use raster::resample to avoid conflicts
    ET <- raster::resample(ET, masked_residual, method = "bilinear")
    names(ET) <- "ET"
    
    lst_layer_name <- paste0("LST_mean-", variable_name)
    lst_file_path <- file.path(raster_folder, paste0(lst_layer_name, ".tif"))
    LST <- raster(lst_file_path)
    LST <- crop(extend(LST, masked_residual), masked_residual)
    
    # Explicitly use raster::resample to avoid conflicts
    LST <- raster::resample(LST, masked_residual, method = "bilinear")
    names(LST) <- "LST"
    
    temp_layer_name <- paste0("Temp_mean-", variable_name)
    temp_file_path <- file.path(raster_folder, paste0(temp_layer_name, ".tif"))
    Temp <- raster(temp_file_path)
    Temp <- crop(extend(Temp, masked_residual), masked_residual)
    
    # Explicitly use raster::resample to avoid conflicts
    Temp <- raster::resample(Temp, masked_residual, method = "bilinear")
    names(Temp) <- "Temp"
    
    evi_layer_name <- paste0("EVI_mean-", variable_name)
    evi_file_path <- file.path(raster_folder, paste0(evi_layer_name, ".tif"))
    EVI <- raster(evi_file_path)
    EVI <- crop(extend(EVI, masked_residual), masked_residual)
    
    # Explicitly use raster::resample to avoid conflicts
    EVI <- raster::resample(EVI, masked_residual, method = "bilinear")
    names(EVI) <- "EVI"
    
    prec_layer_name <- paste0("Prec_sum-", variable_name)
    prec_file_path <- file.path(raster_folder, paste0(prec_layer_name, ".tif"))
    Prec <- raster(prec_file_path)
    Prec <- crop(extend(Prec, masked_residual), masked_residual)
    
    # Explicitly use raster::resample to avoid conflicts
    Prec <- raster::resample(Prec, masked_residual, method = "bilinear")
    names(Prec) <- "Prec"
    
    ndvi_layer_name <- paste0("NDVI_mean-", variable_name)
    ndvi_file_path <- file.path(raster_folder, paste0(ndvi_layer_name, ".tif"))
    NDVI <- raster(ndvi_file_path)
    NDVI <- crop(extend(NDVI, masked_residual), masked_residual)
    
    # Explicitly use raster::resample to avoid conflicts
    NDVI <- raster::resample(NDVI, masked_residual, method = "bilinear")
    NDVI <- NDVI / 10000
    names(NDVI) <- "NDVI"
    
    
    DEM <- raster::resample(crop(extend(raster(file.path(raster_folder, "DEM.tif")), masked_residual), masked_residual), masked_residual, method = "bilinear")
    names(DEM) <- "DEM"
    
    SLOPE <- raster::resample(crop(extend(raster(file.path(raster_folder, "SLOPE.tif")), masked_residual), masked_residual), masked_residual, method = "bilinear")
    names(SLOPE) <- "SLOPE"
    
    # Stack and prepare data for prediction
    km <- stack(SLOPE, DEM, Prec, Temp, NDVI, LST,	EVI, ET)
    km_df <- as.data.frame(km, xy = TRUE)
    km_na = km_df
    km_na <- na.omit(km_na)
    
    # Predict using the trained model
    predictions <- predict(rf_model, km_na)  # Use 'data' instead of 'newdata'
    
    p <- data.frame(x = km_na$x, y = km_na$y, predictions = predictions$predictions)
    
    # Create a SpatialPointsDataFrame from the predictions
    spdf <- SpatialPointsDataFrame(coords = p[, c("x", "y")], data = p)
    
    reference_raster <- masked_residual
    new_raster <- rasterize(spdf, reference_raster)
    new_raster
    # Mask with the residuals to get the final output
    mask_raster <- mask(masked_residual, new_raster$predictions)
    plot(mask_raster)
    
    GWR <- new_raster$predictions+mask_raster
    
    # Save raster output
    output_filename <- file.path(output_folder, paste0("GWR_", year, "-", month, ".tiff"))
    writeRaster(GWR, output_filename, overwrite = TRUE)
    
  }
}



