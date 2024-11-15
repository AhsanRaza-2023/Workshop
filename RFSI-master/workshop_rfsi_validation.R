###################################################### (USE this code)

# install.packages("devtools")
# devtools::install_github("LKremer/ggpointdensity")
library(ggpointdensity)
# Set working directory

base_path <- "C:\\Users\\raza\\Desktop\\rfsi_workshop\\RFSI-master\\RFSI-master\\brandenburg"

setwd(paste0(base_path,"\\temp_data\\train_test_data"))

# Function to calculate metrics and create plot
calculate_metrics_and_plot <- function(df, title) {
  rfsp_10f_pred <- df$predicted_GWD
  rfsp_10f_obs <- df$GWD
  
  #rfsp_10f_pred <- df$Predictions
  #rfsp_10f_obs <- df$Observations
  
  rmse <- sqrt(mean((rfsp_10f_obs - rfsp_10f_pred)^2, na.rm = TRUE))
  ccc <- DescTools::CCC(rfsp_10f_obs, rfsp_10f_pred, ci = "z-transform", conf.level = 0.95, na.rm = TRUE)$rho.c
  mae <- mean(abs(rfsp_10f_obs - rfsp_10f_pred), na.rm = TRUE)
  r2 <- sqrt(abs(summary(lm(rfsp_10f_pred ~ rfsp_10f_obs))$r.squared))
  me <- mean((rfsp_10f_obs - rfsp_10f_pred), na.rm = TRUE)
  
  #cat(title, "RMSE:", rmse, "CCC:", ccc, "MAE:", mae, "R2:", r2, "ME:", me, "\n")
  
  ggplot(df, aes(y = predicted_GWD, x = GWD)) +
  #ggplot(df, aes(y = Predictions, x = Observations)) +
    geom_pointdensity(size = 1, alpha = 0.2) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    labs(x = "Observations, m", y = "Predictions, m", title = bquote(.(title))) +
    theme_bw() +
    theme(
      axis.text = element_text(size = 20),
      axis.title = element_text(size = 22),
      plot.title = element_text(size = 22, hjust = 0.5)
    ) +
    scale_x_continuous(
      breaks = seq(0, 50, by = 5),
      limits = c(0, 50),
      labels = scales::number_format(scale = 0.01 * 100)
    ) +
    scale_y_continuous(
      breaks = seq(0, 50, by = 5),
      limits = c(0, 50),
      labels = scales::number_format(scale = 0.01 * 100)
    ) +
    annotate("text", x = 44, y = 7, label = paste("RMSE:", round(rmse, 2)), hjust = 0.5, vjust = 0.5, size = 7) +
    annotate("text", x = 44, y = 4, label = paste("MAE:", round(mae, 2)), hjust = 0.5, vjust = 0.5, size = 7) +
    annotate("text", x = 44, y = 0, label = bquote(paste("R", ":", ~.(round(r2, 2)))), hjust = 0.5, vjust = 0.5, size = 7)
}

# Function to load data, calculate metrics, create plots, and combine them
process_and_plot <- function(prefix) {
  #years <- c("2001-05", "2006-10", "2011-15", "2016-22")
  years <- c("2016-22")
  plots <- list()
  
  for (i in seq_along(years)) {
    file_name <- paste0(prefix, "_test_", years[i], ".csv")
    #file_name <- paste0(prefix, "_training_", years[i], ".csv")
    df <- read.csv(file_name)
    title <- paste(prefix, "-Test_set_", years[i], sep = "")
    #title <- paste(prefix, "-Train_set_", years[i], sep = "")
    plot <- calculate_metrics_and_plot(df, title)
    plots[[i]] <- plot
  }
  
  combined_plot <- grid.arrange(grobs = plots, ncol = 2, nrow = 2)
  output_file <- paste0("test_", prefix, ".jpeg")
  #output_file <- paste0("train_", prefix, ".jpeg")
  ggsave(output_file, combined_plot, dpi = 300, width = 15, height = 12)
  return(plots)
}

# Process and plot for each prefix
prefixes <- c("rf", "rfsp","svm", "rfsi")

# Collect all plots
all_plots <- list()

for (prefix in prefixes) {
  model_plots <- process_and_plot(prefix)
  all_plots <- c(all_plots, model_plots)
}

# Combine all plots into one
#combined_all_plot <- grid.arrange(grobs = all_plots, ncol = 4, nrow = length(all_plots)/4)
combined_all_plot <- grid.arrange(grobs = all_plots, ncol = 2, nrow = 2)
ggsave("all_models_combined_train.jpeg", combined_all_plot, dpi = 300, width = 20, height = 15)

##########################
######## use for mapping the difference between obs and svm,rf,rfsi,rfsp. training set 


# Load required libraries
library(dplyr)
library(readr)
library(ggplot2)
library(purrr)
library(sf)  # For handling spatial data
library(stringr)

# Define the directory containing the CSV files
file_dir <- paste0(base_path,"\\temp_data\\train_test_data\\")

# Define the path to the shapefile
shapefile_path <- paste0(base_path,"\\borders\\osm\\Brandenburg.shp")

# Read the shapefile
polygon_data <- st_read(shapefile_path)

# Function to read and process each file
read_and_process_file <- function(file) {
  # Read the CSV file
  df <- read_csv(file)
  
  # Rename columns to ensure consistency
  df <- df %>%
    rename_with(~ c("Observations", "Predictions", "Lat", "Lon"), 
                dplyr::matches("Observations|Predictions|Lat|Lon")) %>%
    # Remove unnecessary columns if they exist
    select(-matches("^\\.\\.\\.*")) %>%
    dplyr::group_by(Lat, Lon) %>%
    dplyr::summarize(
      mean_observations = mean(Observations, na.rm = TRUE),
      mean_predictions = mean(Predictions, na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    # Calculate the difference between mean observations and predictions
    dplyr::mutate(difference = mean_observations - mean_predictions) %>%
    # Add a new column with the file name to distinguish data from different files
    dplyr::mutate(source_file = basename(file)) %>%
    # Categorize the difference into classes
    dplyr::mutate(class = case_when(
      difference >= -6 & difference < -4 ~ "Class 1 (-6 - -4)",
      difference >= -4 & difference < -2 ~ "Class 2 (-4 - -2)",
      difference >= -2 & difference < 0 ~ "Class 3 (-2 - 0)",
      difference >= 0 & difference < 2 ~ "Class 4 (0 - 2)",
      difference >= 2  & difference < 4 ~ "Class 5 (2 - 4)",
      difference >= 4  & difference < 6 ~ "Class 6 (4 - 6)",
      difference >= 6  & difference <= 8 ~ "Class 7 (6 - 8)",
      difference >= 8  & difference <= 10 ~ "Class 8 (8 - 10)"
      #TRUE ~ "Other"
    )) %>%
    # Convert to sf object
    st_as_sf(coords = c("Lon", "Lat"), crs = 4326)  # Assuming WGS 84
  
  return(df)
}

# Read and combine all CSV files
csv_files <-list.files(file_dir, pattern = "training", full.names = TRUE)

combined_data <- purrr::map_dfr(csv_files, read_and_process_file)

combined_data <- combined_data %>% 
  filter(difference >= -5 & difference <= 10)

combined_data <- combined_data %>%
  mutate(
    method = str_extract(source_file, "^[a-zA-Z]+") %>%
      str_to_upper(),  
    details = source_file %>%
      str_replace("^[a-zA-Z]+_", "") %>%
      str_replace("training_", "") %>%
      str_replace(".csv$", "")
  )

# Check the resulting data frame
print(combined_data)

# Define color scheme for 8 classes
class_colors <- c(
  "Class 1 (-6 - -4)" = "blue",
  "Class 2 (-4 - -2)" = "pink",
  "Class 3 (-2 - 0)" = "orange",
  "Class 4 (0 - 2)" = "brown",
  "Class 5 (2 - 4)" = "purple",
  "Class 6 (4 - 6)" = "darkgrey",
  "Class 7 (6 - 8)" = "green",
  "Class 8 (8 - 10)" = "yellow",
  "Other" = "red"  # Added for completeness
)

# Define the desired order of the methods and details
desired_method_order <- c("RF", "RFSP", "SVM", "RFSI")
#desired_details_order <- c("2001-05", "2006-10", "2011-15", "2016-22")  # Adjust according to your data
desired_details_order <- c("2016-22")  # Adjust according to your data
# Convert method and details columns to factors with specified order
combined_data <- combined_data %>%
  mutate(
    method = factor(method, levels = desired_method_order),
    details = factor(details, levels = desired_details_order)
  )

# Create a combined plot with improvements
combined_plot <- ggplot() +
  geom_sf(data = polygon_data, fill = "white", color = "black") +  # Plot polygons
  geom_sf(data = combined_data, aes(color = class), size = 4) +  # Plot points on top
  scale_color_manual(values = class_colors) +
  labs(
    title = "Difference Plot for All Files",
    x = "Longitude",
    y = "Latitude",
    color = "Difference (m)"
  ) +
  theme_minimal(base_size = 14) +  # Adjust base font size
  theme(
    plot.title = element_blank(),  # Centered and bold title
    axis.title = element_blank(),  # Remove axis titles
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    panel.grid = element_blank(),  # Remove grid lines
    legend.title = element_text(size = 40),
    legend.position = "bottom",
    legend.text = element_text(size = 40),
    #strip.text = element_text(size = 40,  face = "bold")  # Bold facet labels
  ) +
  facet_wrap(method ~ details, nrow = 2, ncol = 2) +
  theme(strip.text.y = element_text(size = 40,  face = "bold", angle = 180))+
  theme(strip.text.x = element_text(size = 40,  face = "bold", angle = 0))
# Get the ggplot grob
#gt <- ggplotGrob(combined_plot)

# Locate the tops of the plot panels
#panels <- grep("panel", gt$layout$name)
#top <- unique(gt$layout$t[panels])

# Remove the rows immediately above the plot panel
#gt = gt[-(top-1), ]

# Save the plot
file_path <- paste0(base_path, "\\temp_data\\combined_plot.jpeg")

# Save the combined plot using ggsave
ggsave(filename = file_path, plot = combined_plot, dpi = 100, width = 35, height = 30)


combined_data_sf <- st_as_sf(combined_data)

# Create bar plots for difference by class, grouped by source_file
bar_plot <- ggplot(combined_data, aes(x = class, fill = source_file)) +
  geom_bar(stat = "count") +  # Count occurrences of each class per source_file
  labs(
    title = "Count of Differences by Class and Source File",
    x = "Class",
    y = "Count",
    fill = "Source File"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 22),
    axis.title = element_text(size = 22)
  ) +
  facet_wrap(~source_file, nrow = 2, ncol = 2)

file_path <- paste0(base_path, "\\temp_data\\difference_by_source_file_and_class.png")

# Save the combined plot using ggsave
ggsave(filename = file_path, plot = bar_plot, dpi = 50, width = 25, height = 20)

########################################################### historic data gap filling graphs

raster_directory <- paste0(base_path,"\\plot\\del")

# List raster files matching the pattern "monthly_svm_yyyy-mm.tif"
raster_files <- list.files(path = raster_directory, 
                           pattern = "monthly_rfsi_[0-9]{4}-[0-9]{2}\\.tif$", 
                           full.names = TRUE)

# Define the coordinates at which to extract values (replace with your coordinates)
coordinates <- data.frame(lon = 13.94207903, lat = 51.96220392)  # Example coordinates  lon = 13.9427, lat = 51.9598) 
sp_points <- SpatialPoints(coordinates, proj4string = CRS("+proj=longlat +datum=WGS84"))


# Initialize an empty vector to store extracted values
extracted_values <- numeric()

# Loop through the raster files and extract the values at the given coordinates
for (raster_file in raster_files) {
  # Load the raster
  raster_data <- raster(raster_file)
  crs(raster_data)
  # Extract the value at the specified coordinates
  value <- extract(raster_data, sp_points)
  
  # Append the extracted value to the vector
  extracted_values <- c(extracted_values, value)
}

# Extract dates from file names
dates <- sub("monthly_rfsi_([0-9]{4}-[0-9]{2}).tif", "\\1", basename(raster_files))
dates <- as.Date(paste0(dates, "-01"), "%Y-%m-%d")

# Create a data frame for the time series
time_series_data <- data.frame(date = dates, value = extracted_values)
time_series_data$value <- time_series_data$value+0.5

time_series_data <- time_series_data %>%
  filter(date >= as.Date("2016-01-01") & date <= as.Date("2022-12-01"))

# time_series_data <- time_series_data %>%
#   filter(!(date >= as.Date("2009-01-01") & date <= as.Date("2009-12-01")))

#head(time_series_data)
#plot(time_series_data$date, time_series_data$value)

# Assume the CSV has columns 'date' and 'GWD' (Ground Water Depth)
observed_data <- read.csv(paste0(base_path, "\\temp_data\\clean_2016-22.csv"))
head(observed_data)

# Filter the observed data for the specific coordinates
filtered_observed_data <- observed_data %>%
  filter(Lon == 13.94207903 & Lat == 51.96220392)    #Lon == 14.30399953 & Lat == 52.70247148 (svm,rf,rfsp)
#Lon == 13.94207903 & Lat == 51.96220392)
#Lon == 13.37630435 & Lat == 52.1047938
#Lon == 12.42743346 & Lat == 53.09865374
# Convert the 'date' column in the filtered observed data to Date format
filtered_observed_data$Date <- as.Date(paste0(filtered_observed_data$Date, "-01"), format = "%b-%y-%d")

filtered_observed_data <- na.omit(filtered_observed_data)

filtered_observed_data <- filtered_observed_data %>%
  select(GWD, Lat, Lon, Date)

# filtered_observed_data <- filtered_observed_data %>%
#    filter(Date >= as.Date("2005-01-01") & Date <= as.Date("2021-12-01"))

head(filtered_observed_data)

a <- ggplot() +
  geom_line(data = time_series_data, aes(x = date, y = value), color = "blue", linetype = "solid") +
  geom_line(data = filtered_observed_data, aes(x = Date, y = GWD), color = "red", linetype = "dashed") +
  geom_rect(aes(xmin = as.Date("2016-01-01"), xmax = as.Date("2021-12-31"), ymin = -Inf, ymax = Inf),
            fill = "grey", alpha = 0.3, inherit.aes = FALSE) +  # Add grey shaded area
  scale_x_date(date_breaks = "2 year", date_labels = "%Y") +
  labs(title = "13.94207903; 51.96220392",
       x = "Date",
       y = "GWD (m)",
       color = "Source") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  theme_classic() + 
  guides(color = guide_legend(nrow = 3), override.aes = list(fill = NA)) +
  theme(
    legend.position = c(.15, .15),  # Adjusted to show legend
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 15),   # Adjust legend text size
    axis.title.x = element_text(size = 16, face = "bold"),  # Make axis titles bold
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(size = 16, angle = 0, hjust = 1),  # Rotate x-axis text for better readability
    axis.text.y = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3) # Add bold title
  ) +
  theme(panel.grid.major.y = element_line(color = "grey",
                                          size = 0.4,
                                          linetype = 2))

a


######################### for validating rfsi with grace (use this one)

library(sp)
library(raster)
library(ggplot2)

# Define the function to process and plot data for a given location
plot_time_series <- function(lat, lon, rfsi_path, grace_path, title_text, output_path) {
  
  # Create a SpatialPoints object for the given location
  location <- SpatialPoints(cbind(lon, lat), proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # ---- First Time Series (RFSI) ----
  files1 <- list.files(path = rfsi_path, pattern = "monthly_rfsi_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)
  rasters1 <- stack(files1)
  extracted_values1 <- extract(rasters1, location)
  column_names1 <- colnames(extracted_values1)
  date_strings1 <- gsub("monthly_rfsi_", "", column_names1)
  dates1 <- as.Date(paste0(date_strings1, "-01"), format = "%Y.%m-%d")
  values1 <- as.numeric(extracted_values1[1, ])
  
  time_series_df1 <- data.frame(date = dates1, value = scale(values1)* (-1), series = "rfsi")
  
  # ---- Second Time Series (GRACE-GWSA) ----
  files2 <- list.files(path = grace_path, pattern = "\\.tif$", full.names = TRUE)
  rasters2 <- stack(files2)
  extracted_values2 <- extract(rasters2, location)
  column_names2 <- colnames(extracted_values2)
  date_strings2 <- gsub("GWR_", "", column_names2)
  date_strings2 <- gsub("_bb", "", date_strings2)
  dates <- as.Date(paste0(date_strings2, "-01"), format = "%Y.%m-%d")
  values2 <- as.numeric(extracted_values2[1, ])
  
  time_series_df2 <- data.frame(date = dates, value = scale(values2) * (1), series = "GRACE")
  
  # ---- Combine Both Time Series Data Frames ----
  combined_time_series_df <- rbind(time_series_df1, time_series_df2)
  combined_time_series_df <- na.omit(combined_time_series_df)
  
  combined_time_series_df1 <- merge(time_series_df1, time_series_df2, by = "date")
  combined_time_series_df1 <- na.omit(combined_time_series_df1)
  combined_time_series_df$date <- as.Date(combined_time_series_df$date)
  
  names(combined_time_series_df1)[names(combined_time_series_df1) == "value.x"] <- combined_time_series_df1$series.x[1]
  names(combined_time_series_df1)[names(combined_time_series_df1) == "value.y"] <- combined_time_series_df1$series.y[1]
  
  # Optionally, you can remove the series columns after renaming
  combined_time_series_df1 <- combined_time_series_df1 %>%
    select(-series.x, -series.y)
  
  rmse <- sqrt(mean((combined_time_series_df1$rfsi - combined_time_series_df1$GRACE)^2, na.rm = TRUE))
  r_squared <- cor(combined_time_series_df1$rfsi , combined_time_series_df1$GRACE)
  bias <- mean(combined_time_series_df1$rfsi - combined_time_series_df1$GRACE)
  nse_value <- 1 - (sum((combined_time_series_df1$rfsi - combined_time_series_df1$GRACE)^2) / sum((combined_time_series_df1$rfsi - mean(combined_time_series_df1$rfsi))^2))
  
  # cat("rmse_value:", rmse_value, "\n")
  # cat("r_squared:", r_squared, "\n")
  # cat("bias:", bias, "\n")
  # cat("nse_value:", nse_value, "\n")
  
  r_value <- cor(combined_time_series_df1$rfsi, combined_time_series_df1$GRACE)
  
  # Create the plot
  p <- ggplot(combined_time_series_df1, aes(x = date)) +
    geom_line(aes(y = rfsi, color = "RFSI"), linetype = "solid") +
    geom_line(aes(y = GRACE, color = "GRACE"), linetype = "solid") +
    scale_x_date(date_breaks = "4 year", date_labels = "%Y") +
    labs(title = title_text,
         x = "Date",
         y = "Standardized Values",
         color = "Series") +
    theme_minimal() +
    theme_classic() +
    guides(color = guide_legend(nrow = 2, override.aes = list(fill = NA))) + 
    theme(
      legend.position = c(.12, .15),
      legend.title = element_blank(),
      legend.text = element_text(size = 16),
      legend.key = element_rect(fill = "transparent", color = NA),  # Removes legend text background
      axis.title.x = element_text(size = 16, face = "bold"),
      axis.title.y = element_text(size = 16, face = "bold"),
      axis.text.x = element_text(size = 16, angle = 0, hjust = 1),
      axis.text.y = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold"),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
      panel.grid.major.y = element_line(color = "grey", size = 0.4, linetype = 2)
    ) +
    annotate("text", x = max(combined_time_series_df1$date) - 100, y = max(c(combined_time_series_df1$rfsi, combined_time_series_df1$GRACE)) + 0.5, 
             label = paste0("R = ", round(r_squared, 3)), size = 6, hjust = 1, color = "black")
  # Save the plot
  ggsave(output_path, p, dpi = 100, width = 6.5, height = 4)
}

rfsi_path = paste0(base_path,"\\plot\\del")
grace_path = "C:\\Users\\raza\\Desktop\\rfsi_workshop\\downscaling\\brandenburg"

# Define locations and paths
locations <- list(
  list(lat = 52.70247148, lon = 14.30399953, title = "14.30399953; 52.70247148",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota1.jpeg")),
  list(lat = 51.96220392, lon = 13.94207903, title = "13.94207903; 51.96220392",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota2.jpeg")),
  list(lat = 52.1047938, lon = 13.37630435, title = "13.37630435; 52.1047938",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota3.jpeg")),
  list(lat = 53.042966, lon = 12.85875225, title = "12.85875225; 53.042966",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota4.jpeg")),
  list(lat = 52.61566363, lon = 12.96818566, title = "12.96818566; 52.61566363",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota5.jpeg")),
  list(lat = 52.79362177, lon = 14.14333931, title = "14.14333931; 52.79362177",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota6.jpeg")),
  list(lat = 53.11452849, lon = 13.50875261, title = "13.50875261; 53.11452849",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota7.jpeg")),
  list(lat = 52.69238305, lon = 13.03130551, title = "13.03130551; 52.69238305",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota8.jpeg")),
  list(lat = 53.03788146, lon = 13.18787977, title = "13.18787977; 53.03788146",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota9.jpeg")),
  list(lat = 53.02569177, lon = 13.88680919, title = "13.88680919; 53.02569177",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota10.jpeg")),
  list(lat = 52.04149066, lon = 12.8508426, title = "12.8508426; 52.04149066",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota11.jpeg")),
  list(lat = 51.79708216, lon = 13.56454914, title = "13.56454914; 51.79708216",
       rfsi_path = rfsi_path,
       grace_path = grace_path,
       output_path = paste0(base_path,"\\temp_data\\plota12.jpeg")),
)

# Loop through locations and create plots
for (loc in locations) {
  plot_time_series(loc$lat, loc$lon, loc$rfsi_path, loc$grace_path, loc$title, loc$output_path)
}


# Extract latitudes and longitudes
lat_lon_df <- do.call(rbind, lapply(locations, function(x) data.frame(lat = x$lat, lon = x$lon)))

# Add row names as "Location"
lat_lon_df$Location <- paste0("Location_", seq_len(nrow(lat_lon_df)))

# Save to an Excel file
#write.csv(lat_lon_df, file = "C:/Users/raza/Downloads/lat_lon_locations.csv", row.names = FALSE)

# Print the data frame to verify
print(lat_lon_df)

########################## corelation between rfsi and grace downscaled brandenburg

# Load necessary libraries
library(raster)
library(ggplot2)
library(dplyr)

# Set your file paths
rfsi_path <- paste0(base_path, "\\plot\\del")
grace_path <- "C:\\Users\\raza\\Desktop\\rfsi_workshop\\downscaling\\brandenburg"

# List all raster files for both datasets
rfsi_files <- list.files(path = rfsi_path, pattern = "monthly_rfsi_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)
grace_files <- list.files(path = grace_path, pattern = "GWR_\\d{4}-\\d{2}_.*\\.tif$", full.names = TRUE)

# Extract years from filenames for RFSI
rfsi_years <- as.numeric(sub("monthly_rfsi_(\\d{4})-\\d{2}\\.tif$", "\\1", basename(rfsi_files)))

# Initialize a list to store RFSI yearly mean values
rfsi_yearly_means <- list()

# Process each year for RFSI
unique_rfsi_years <- unique(rfsi_years)
for (year in unique_rfsi_years) {
  yearly_files <- rfsi_files[rfsi_years == year]
  raster_stack <- stack(yearly_files)
  yearly_mean_raster <- calc(raster_stack, fun = mean)
  mean_value <- cellStats(yearly_mean_raster, stat = "mean")
  rfsi_yearly_means[[as.character(year)]] <- mean_value
}

# Extract years from filenames for GRACE
grace_years <- as.numeric(sub("GWR_(\\d{4})-\\d{2}_.*\\.tif$", "\\1", basename(grace_files)))

# Initialize a list to store GRACE yearly mean values
grace_yearly_means <- list()

# Process each year for GRACE
unique_grace_years <- unique(grace_years)
for (year in unique_grace_years) {
  yearly_files <- grace_files[grace_years == year]
  raster_stack <- stack(yearly_files)
  yearly_mean_raster <- calc(raster_stack, fun = mean)
  mean_value <- cellStats(yearly_mean_raster, stat = "mean")
  grace_yearly_means[[as.character(year)]] <- mean_value
}

# Convert results to data frames
rfsi_results <- data.frame(
  Year = as.numeric(names(rfsi_yearly_means)),
  MeanValue = unlist(rfsi_yearly_means),
  Dataset = "RFSI"
)

rfsi_results$MeanValue <- (rfsi_results$MeanValue) * -1
rfsi_results$MeanValue <- scale(rfsi_results$MeanValue)

grace_results <- data.framgrace_results <- data.framgrace_results <- data.frame(
  Year = as.numeric(names(grace_yearly_means)),
  MeanValue = unlist(grace_yearly_means),
  Dataset = "GRACE"
)
grace_results$MeanValue <- (grace_results$MeanValue)/1000
grace_results$MeanValue <- scale(grace_results$MeanValue)

# Combine results into one data frame
combined_results <- bind_rows(rfsi_results, grace_results)
combined_results <- combined_results %>%
  filter(!(Year %in% c(2001, 2021, 2022)))


# Plot the results with the legend inside the plot area
d <- ggplot(combined_results, aes(x = Year, y = MeanValue, color = Dataset, shape = Dataset)) +
  geom_line(size = 1) +
  labs(title = "",
       x = "Date",
       y = "Standarized values") +
  scale_x_continuous(breaks = seq(min(combined_results$Year), max(combined_results$Year), by = 2)) +
  scale_color_manual(values = c("RFSI" = "deepskyblue3", "GRACE" = "sienna2")) + 
  scale_shape_manual(values = c("RFSI" = 16, "GRACE" = 17)) +  # Use different shapes for clarity
  theme_bw(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = c(0.1, 0.2),  # Position the legend inside the plot area
    legend.background = element_rect(fill = "white", color = "black"),  # Optional: add background and border to legend
    legend.box.background = element_rect(color = "black")  # Optional: add border to legend box
  )
d

#ggsave("C:\\Users\\raza\\Downloads\\RFSI-master\\RFSI-master\\brandenburg\\temp_data\\images\\cor_rfsi_grace.jpeg", d, dpi = 100, width = 9, height = 5)


####################################### heat map

latitude <- 51.96220392
longitude <- 13.94207903
#Lon == 14.30399953 & Lat == 52.70247148  2
#Lon == 13.94207903 & Lat == 51.96220392) 3
#Lon == 13.37630435 & Lat == 52.1047938    4
#Lon == 12.42743346 & Lat == 53.09865374  1

raster_dir <- paste0(base_path, "\\plot\\del")
# List all raster files
raster_files <- list.files(path = raster_dir, pattern = "monthly_rfsi_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)

# Function to extract data from a raster file at a specific point
extract_point_data <- function(file_path, lat, lon) {
  raster_data <- raster(file_path)
  
  # Convert the point coordinates to a SpatialPoints object
  point <- SpatialPoints(cbind(lon, lat), proj4string = crs(raster_data))
  
  # Extract the value at the point
  value <- extract(raster_data, point)
  
  # Extract date from filename (assuming filename format is 'monthly_rfsi_YYYY-MM.tif')
  date_part <- gsub("monthly_rfsi_|\\.tif", "", basename(file_path))
  date <- as.Date(paste0(date_part, "-01"), format = "%Y-%m-%d")
  
  # Return a data frame with the date and the extracted value
  data.frame(Date = date, Value = value)
}

# Extract data from all raster files
point_data_list <- lapply(raster_files, extract_point_data, lat = latitude, lon = longitude)
all_point_data <- do.call(rbind, point_data_list)

# Add Year and Month columns
all_point_data <- all_point_data %>%
  mutate(Year = year(Date),
         Month = month(Date, label = TRUE, abbr = TRUE))

# Reshape data for heatmap
heatmap_data <- all_point_data %>%
  select(Year, Month, Value) %>%
  spread(key = Month, value = Value)

# Convert Year and Month to factors with appropriate levels
heatmap_data$Year <- factor(heatmap_data$Year, levels = unique(heatmap_data$Year))

# Melt the data for ggplot2
heatmap_long <- melt(heatmap_data, id.vars = "Year")

heatmap_long$Year <- as.factor(heatmap_long$Year)

heatmap_long <- heatmap_long %>%
  mutate(Year = factor(Year, levels = unique(Year))) %>%
  arrange(Year)

# Create the heatmap with reduced tile size
heatmap_plot <- ggplot(heatmap_long, aes(x = variable, y = Year, fill = value)) +
  geom_tile(color = "white", size = 0.1) +  # Add white borders around tiles
  scale_fill_viridis_c(
    option = "C",
    name = "GWL (m)",  # Update legend title here
    labels = scales::comma_format()  # Format legend labels with commas
  ) +
  labs(
    #title = "13.37630435, 52.1047938",
    x = "",
    y = ""
  ) +
  theme_minimal(base_size = 15) +  # Increase base font size
  theme(
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size = 22),
    axis.text.y = element_text(size = 22),
    #axis.text.y = element_blank(),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_text(size = 20),
    legend.position = "right",  # Remove the legend
    legend.title = element_blank(),
    legend.text = element_blank()
  ) +
  #coord_fixed(ratio = 0.9) +  # Adjust the aspect ratio to make tiles narrower
  scale_y_discrete(breaks = seq(2001, 2022, by = 4))  # Set y-axis breaks to every 4 years

heatmap_plot

#ggsave(filename = "C:\\Users\\raza\\Downloads\\RFSI-master\\RFSI-master\\brandenburg\\temp_data\\images\\heatmap_plot_33.jpeg", plot = heatmap_plot, width = 15, height = 3, units = "in", dpi = 100)

