# Install and load required packages
#install.packages(c("raster", "sf", "tidyverse", "zoo"))
library(raster)
library(sf)
library(tidyverse)
library(zoo)

###################
##################      (USE this code)
##################################

base_path <- "C:\\Users\\raza\\Desktop\\rfsi_workshop\\downscaling"

# # List of CSV files
csv_files <- list.files(paste0(base_path, "\\all_data"), pattern = "\\.csv$", full.names = TRUE)

# Read and extract ST_ID and Date from each CSV file
data_list <- map(csv_files, ~ read.csv(.x, stringsAsFactors = FALSE) %>%
                   select(ST_ID, Date, Lat, Lon, GWD))

# Combine data from all CSV files
all_data1 <- bind_rows(data_list)

all_data <- all_data1

all_data$Date <- parse_date_time(all_data$Date, orders = c("%y-%b", "%b-%y", "%Y-%m"))
# Convert Date to Date type
all_data$Date <- as.Date(all_data$Date, format = "%Y-%m")

# Remove rows with missing dates
all_data <- all_data[complete.cases(all_data$Date), ]

# Group data by ST_ID
grouped_data <- all_data %>%
  group_by(ST_ID) %>%
  dplyr::summarise(Start_Date = min(Date),
                   End_Date = max(Date))

# Display the summary data
print(grouped_data)

filtered_data <- grouped_data %>%
  dplyr::filter(Start_Date == as.Date("2015-01-01") & End_Date == as.Date("2020-12-01"))

selected_ST_ID <- filtered_data$ST_ID

# Filter all_data based on selected ST_ID values
extracted_data <- all_data %>%
  filter(ST_ID %in% selected_ST_ID)



########################################################################################

# Load required libraries
library(raster)
library(dplyr)
library(ggplot2)

# Define working directories
downscaled_output_dir <- paste0(base_path, "\\downscaled_output")
rf_dir <- paste0(base_path, "\\rf")
gwsa_dir <- paste0(base_path, "\\GWSA")

# Function to extract values from raster at specified coordinates
extract_raster_values <- function(raster_file, points) {
  raster_data <- raster(raster_file)
  values <- extract(raster_data, points)
  return(values)
}

# Function to process raster files based on coordinates and date
process_raster_files <- function(coordinates_df, raster_folder, file_pattern) {
  unique_dates <- unique(coordinates_df$Date)
  result_df <- data.frame()
  
  for (date in unique_dates) {
    date_data <- coordinates_df %>% filter(Date == date)
    
    if (nrow(date_data) > 0) {
      lon <- date_data$Lon
      lat <- date_data$Lat
      ST_ID <- date_data$ST_ID
      GWD <- date_data$GWD
      
      # Construct the raster file name based on the date and pattern
      raster_file <- file.path(raster_folder, sprintf(file_pattern, format(as.Date(date), "%Y-%m")))
      
      if (file.exists(raster_file)) {
        values <- extract_raster_values(raster_file, cbind(lon, lat))
        new_rows <- data.frame(ST_ID = ST_ID, Date = as.Date(date), Lon = lon, Lat = lat, GWD = GWD, values = values)
        result_df <- bind_rows(result_df, new_rows)
      } else {
        cat("Raster file not found for date:", format(as.Date(date), "%Y-%m-%d"), "\n")
      }
    } else {
      cat("No data found for date:", format(as.Date(date), "%Y-%m-%d"), "\n")
    }
  }
  return(na.omit(result_df))
}

# Example: Process raster files with appropriate file patterns
result_rf <- process_raster_files(extracted_data, rf_dir, "GWR_%s.tiff")
result_mgwr <- process_raster_files(extracted_data, downscaled_output_dir, "GWR_%s.tiff")
result_g3p <- process_raster_files(extracted_data, gwsa_dir, "GWSA_monthly_%s.tif")

# Combine the results into a single dataframe
combine <- reduce(list(result_rf, result_mgwr, result_g3p), merge, by = c("ST_ID", "Date", "Lon", "Lat", "GWD"))

combined_data <- combine %>%
  dplyr::rename(
    ST_ID = ST_ID,
    Lon = Lon,
    Lat = Lat,
    GWD = GWD,
    rf = values.x,
    mgwr = values.y,
    g3p = values
  )

#########################################################   Performance analysis for specific well

combined_data1 <- combined_data
# Filter data for a specific station and date range
combined_data1 <- combined_data1 %>%
  filter(ST_ID == 34495032, Date >= as.Date("2015-01-01") & Date <= as.Date("2020-12-31"))

combined_data1$GWD <- combined_data1$GWD * -1

# Standardize values and apply transformations
combined_data1 <- combined_data1 %>%
  dplyr::mutate(across(c(rf, mgwr, g3p, GWD), ~ scale(.) * 1))

# Function to calculate R-squared
calculate_r2 <- function(observed, predicted) {
  cor(observed, predicted)^2
}

# Calculate R-squared for each model
r2_rf <- calculate_r2(combined_data1$GWD, combined_data1$rf)
r2_mgwr <- calculate_r2(combined_data1$GWD, combined_data1$mgwr)
r2_g3p <- calculate_r2(combined_data1$GWD, combined_data1$g3p)


# Improved Plot
d <- ggplot(combined_data1, aes(x = Date)) +
  geom_line(aes(y = GWD, color = "GWD"), size = 0.8, linetype = "dashed", na.rm = TRUE) +
  geom_line(aes(y = rf, color = "RF"), size = 0.8, linetype = "solid", na.rm = TRUE) +
  geom_line(aes(y = mgwr, color = "MGWR"), size = 0.8, linetype = "solid", na.rm = TRUE) +
  geom_line(aes(y = g3p, color = "G3P"), size = 0.8, linetype = "solid", na.rm = TRUE) +
  labs(
    title = "Comparison of Groundwater Depth and Predicted Values",
    subtitle = paste(
      "Station ID: 34495032 (2015-2020)\n",
      "R² (RF) =", round(r2_rf, 3), 
      ", R² (MGWR) =", round(r2_mgwr, 3),
      ", R² (G3P) =", round(r2_g3p, 3)
    ),
    x = "Date",
    y = "Standardized Values",
    color = ""
  ) +
  scale_color_manual(values = c("GWD" = "blue", 
                                "RF" = "red", 
                                "MGWR" = "green", 
                                "G3P" = "purple")) +
  theme_bw() +
  theme(
    text = element_text(size = 16),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 16, face = "italic"),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom",
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80")
  )

d

ggsave(filename = paste0(base_path, "\\validation_GRACE_downscaling.jpeg"), d, dpi = 100, width = 9, height = 5)

######################################################################################## Performance analysis for all wells
library(dplyr)
result_df1 <-combined_data

result_df1 <- na.omit(result_df1)

# Helper functions for RMSE, MSE, and NSE
calculate_rmse <- function(actual, predicted) {
  sqrt(mean((actual - predicted)^2, na.rm = TRUE))
}

calculate_mse <- function(actual, predicted) {
  mean((actual - predicted)^2, na.rm = TRUE)
}

calculate_nse <- function(actual, predicted) {
  1 - sum((actual - predicted)^2) / sum((actual - mean(actual))^2)
}

# Process data
result_df1 <- combined_data %>%
  na.omit() %>%
  dplyr::group_by(ST_ID) %>%
  dplyr::mutate(
    rf = scale(rf),
    mgwr = scale(mgwr),
    g3p = scale(g3p),
    GWD = scale(GWD * -1) 
  )

# Calculate correlation, RMSE, MSE, and NSE for each model
metrics_result <- result_df1 %>%
  dplyr::group_by(ST_ID) %>%
  dplyr::summarize(
    correlation_rf = cor(GWD, rf),
    correlation_mgwr = cor(GWD, mgwr),
    correlation_g3p = cor(GWD, g3p),
    
    rmse_rf = calculate_rmse(GWD, rf),
    rmse_mgwr = calculate_rmse(GWD, mgwr),
    rmse_g3p = calculate_rmse(GWD, g3p),
    
    mse_rf = calculate_mse(GWD, rf),
    mse_mgwr = calculate_mse(GWD, mgwr),
    mse_g3p = calculate_mse(GWD, g3p),
    
    nse_rf = calculate_nse(GWD, rf),
    nse_mgwr = calculate_nse(GWD, mgwr),
    nse_g3p = calculate_nse(GWD, g3p)
  )

metrics_result <- merge(metrics_result,result_df1, by = "ST_ID")
metrics_result <- distinct(metrics_result, ST_ID, .keep_all = TRUE)

# Write the results to CSV files
write.csv(metrics_result, "test_metrics.csv", row.names = FALSE)

#######################################
########################################## cc
######################################

shp <- readOGR(paste0(base_path, "\\Germany_shapefile\\Brandenburg.shp"))
shp <- spTransform(shp, CRS("+proj=longlat +datum=WGS84"))

data1 <- metrics_result

data1$class <- cut(data1$correlation_g3p, breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c("0.0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"))
data1 <- na.omit(data1)
# Plot a bar chart
ggplot(data1, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Data Classification by cc", x = "Class", y = "Count")

ggplot1 <- ggplot(data1, aes(x = Lon, y = Lat, color = class)) +
  geom_point(size = 4, shape = 16) +
  labs(title = "",
       color = "Class") +
  scale_color_manual(values = c("#FF0000", "#00FF00",  "#87CEEB","#0000FF")) +
  theme_void()+ theme(legend.title = element_text(size = 26)
                      ,legend.text=element_text(size=26))

ggplot1 <- ggplot1 + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)



data1$class <- cut(data1$correlation_mgwr, breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c("0.0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"))
data1 <- na.omit(data1)
data1 <- data1[data1$ST_ID %in% data1$ST_ID, ]
# Plot a bar chart
ggplot(data1, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Data Classification by cc", x = "Class", y = "Count")
count(data1$class)

ggplot2 <- ggplot(data1, aes(x = Lon, y = Lat, color = class)) +
  geom_point(size = 4, shape = 16) +
  labs(title = "",
       color = "Class") +
  scale_color_manual(values = c("#FF0000", "#00FF00",  "#87CEEB","#0000FF")) +
  theme_void() + theme(legend.title = element_text(size = 26)
                       ,legend.text=element_text(size=26))

ggplot2 <- ggplot2 + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

data1$class <- cut(data1$correlation_rf, breaks = c(0.0, 0.25, 0.5, 0.75, 1.0), labels = c("0.0-0.25", "0.25-0.5", "0.5-0.75", "0.75-1.0"))
data1 <- na.omit(data1)
data1 <- data1[data1$ST_ID %in% data1$ST_ID, ]
# Plot a bar chart
ggplot(data1, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Data Classification by cc", x = "Class", y = "Count")

# Now, create the plot
ggplot3 <- ggplot(data1, aes(x = Lon, y = Lat, color = class)) +
  geom_point(size = 4, shape = 16) +
  labs(title = "",
       color = "Class", size = 16) +
  scale_color_manual(values = c("#FF0000", "#00FF00",  "#87CEEB","#0000FF")) +
  theme_void()+ theme(legend.title = element_text(size = 26),
                      legend.text=element_text(size=26))

ggplot3 <- ggplot3 + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

# Arrange the plots using ggarrange
combined_plot <- ggarrange(ggplot1, ggplot2, ggplot3, ncol = 3, labels = c("A", "B", "C"), nrow = 1, common.legend = TRUE, legend = "bottom",
                           font.label = list(size = 36))

ggsave("combined_plot1_cc.png", combined_plot, units = "cm", width = 60, height = 25, dpi = 250)


a <- mean(data1$correlation)
b <-mean(data3$correlation)
((b-a)/a) *100

# Calculate percentage difference
percentage_difference <- data.frame(
  count_data1 = table(data1$class),
  count_data3 = table(data3$class),
  percentage_difference = ((table(data3$class) - table(data1$class)) / table(data1$class)) * 100
)

percentage_difference

# Extract relevant columns
change_data <- percentage_difference[, c("count_data1.Freq", "count_data2.Freq")]
# Calculate the overall change
change_data$OverallChange <- change_data$count_data2.Freq - change_data$count_data1.Freq
# Print the result
print(change_data)
mean

###################  RMSE
library(classInt)

data1$class <- cut(data1$rmse_g3p, breaks = c(0.532753, 0.711822, 1.090892, 1.369961, 1.649031), labels = c("0.53-0.71", "0.71-1.09", "1.09-1.36", "1.36-1.64"))
data1 <- na.omit(data1)
# Plot a bar chart
ggplot(data1, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Data Classification by cc", x = "Class", y = "Count")

ggplot1 <- ggplot(data1, aes(x = Lon, y = Lat, color = class)) +
  geom_point(size = 4, shape = 16) +
  labs(title = "",
       color = "Class") +
  scale_color_manual(values = c("#FF0000", "#00FF00",  "#87CEEB","#0000FF")) +
  theme_void()+ theme(legend.title = element_text(size = 26)
                      ,legend.text=element_text(size=26))

ggplot1 <- ggplot1 + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)


data1$class <- cut(data1$rmse_mgwr, breaks = c(0.532753, 0.711822, 1.090892, 1.369961, 1.649031), labels = c("0.53-0.71", "0.71-1.09", "1.09-1.36", "1.36-1.64"))
data1 <- na.omit(data1)
data1 <- data2[data1$ST_ID %in% data1$ST_ID, ]
# Plot a bar chart
ggplot(data1, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Data Classification by cc", x = "Class", y = "Count")

ggplot2 <- ggplot(data1, aes(x = Lon, y = Lat, color = class)) +
  geom_point(size = 4, shape = 16) +
  labs(title = "",
       color = "Class") +
  scale_color_manual(values = c("#FF0000", "#00FF00",  "#87CEEB","#0000FF")) +
  theme_void() + theme(legend.title = element_text(size = 26)
                       ,legend.text=element_text(size=26))

ggplot2 <- ggplot2 + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

data1$class <- cut(data1$rmse_rf, breaks = c(0.532753, 0.711822, 1.090892, 1.369961, 1.649031), labels = c("0.53-0.71", "0.71-1.09", "1.09-1.36", "1.36-1.64"))
data1 <- na.omit(data1)
data1 <- data3[data1$ST_ID %in% data1$ST_ID, ]
# Plot a bar chart
ggplot(data1, aes(x = class, fill = class)) +
  geom_bar() +
  labs(title = "Data Classification by cc", x = "Class", y = "Count")

# Now, create the plot
ggplot3 <- ggplot(data1, aes(x = Lon, y = Lat, color = class)) +
  geom_point(size = 4, shape = 16) +
  labs(title = "",
       color = "Class", size = 16) +
  scale_color_manual(values = c("#FF0000", "#00FF00",  "#87CEEB","#0000FF")) +
  theme_void()+ theme(legend.title = element_text(size = 26),
                      legend.text=element_text(size=26))

ggplot3 <- ggplot3 + geom_polygon(data = shp, aes(x = long, y = lat, group = group), colour = "black", fill = NA)

combined_plot <- ggarrange(ggplot1, ggplot2, ggplot3, ncol = 3, labels = c("A", "B", "C"), nrow = 1, common.legend = TRUE, legend = "bottom",
                           font.label = list(size = 36))

ggsave("combined_plot_rmse1.png", combined_plot, units = "cm", width = 60, height = 25, dpi = 250)


c <- mean(data1$rmse)
d <-mean(data2$rmse)
((d-c)/c) *100

# Calculate percentage difference
percentage_difference <- data.frame(
  count_data1 = table(data1$class),
  count_data3 = table(data3$class),
  percentage_difference = ((table(data2$class) - table(data3$class)) / table(data1$class)) * 100
)

percentage_difference

mean(percentage_difference$percentage_difference.Freq)

# Extract relevant columns
change_data <- percentage_difference[, c("count_data1.Freq", "count_data2.Freq")]
# Calculate the overall change
change_data$OverallChange <- change_data$count_data2.Freq - change_data$count_data1.Freq
# Print the result
print(change_data)
mean

#####################  Time Series and Spatial Accuracy 
#####################  

#####################  Mean monthly GWSA for germany

# Load necessary libraries
library(raster)
library(rgdal)
library(dplyr)
library(ggplot2)
library(tidyr)  # For pivoting data

# Define the paths
raster_path_1 <- paste0(base_path,"\\downscaled_output")
raster_path_2 <- paste0(base_path,"\\rf")
raster_path_3 <- paste0(base_path,"\\GWSA\\resampled_1km")
shapefile_path <- paste0(base_path,"\\Germany_shapefile\\DEU_adm0.shp")

# Read the shapefile
germany_shape <- readOGR(shapefile_path)

# List the raster files
raster_files_1 <- list.files(raster_path_1, pattern = "GWR_\\d{4}-\\d{2}\\.tiff$", full.names = TRUE)
raster_files_2 <- list.files(raster_path_2, pattern = "GWR_\\d{4}-\\d{2}\\.tiff$", full.names = TRUE)
raster_files_3 <- list.files(raster_path_3, pattern = "GWSA_monthly_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)

# Function to calculate mean raster value masked by the shapefile
calculate_mean_raster <- function(raster_file, mask) {
  r <- raster(raster_file)
  masked_r <- mask(r, mask)
  mean_value <- cellStats(masked_r, stat = 'mean', na.rm = TRUE)
  return(mean_value)
}

# Calculate mean values for each raster in each directory
mean_values_1 <- sapply(raster_files_1, calculate_mean_raster, mask = germany_shape)
mean_values_2 <- sapply(raster_files_2, calculate_mean_raster, mask = germany_shape)
mean_values_3 <- sapply(raster_files_3, calculate_mean_raster, mask = germany_shape)

# Extract dates from raster file names
dates_1 <- sub(".*_(\\d{4}-\\d{2})\\.tiff$", "\\1-01", raster_files_1)  # Appending '-01' to represent the first of the month
dates_2 <- sub(".*_(\\d{4}-\\d{2})\\.tiff$", "\\1-01", raster_files_2)
dates_3 <- sub(".*_(\\d{4}-\\d{2})\\.tif$", "\\1-01", raster_files_3)

# Combine into a data frame
mean_data <- data.frame(
  Date = as.Date(c(dates_1, dates_2, dates_3), format = "%Y-%m-%d"),  # Ensure the format matches
  Mean_Value = c(mean_values_1, mean_values_2, mean_values_3),
  Source = c(rep("GWSA (MGWR)", length(mean_values_1)), 
             rep("GWSA (RF)", length(mean_values_2)), 
             rep("GWSA (Origional)", length(mean_values_3)))
)

# Reshape mean_data to wide format
wide_data <- mean_data %>%
  pivot_wider(names_from = Source, values_from = Mean_Value)

# Calculate R² for MODIS 1 vs GWSA
r2_modis1_gwsa <- cor(wide_data$`GWSA (Origional)`, wide_data$`GWSA (MGWR)`, use = "complete.obs")^2

# Calculate R² for MODIS 2 vs GWSA
r2_modis2_gwsa <- cor(wide_data$`GWSA (RF)`, wide_data$`GWSA (Origional)`, use = "complete.obs")^2

a <- ggplot(mean_data, aes(x = Date, y = Mean_Value, color = Source)) +
  geom_line(size = 0.8) +  # Increase line size for better visibility
  labs(title = "", 
       x = "Date", 
       y = "GWSA (mm)") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.15),  # Position legend inside the plot
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.key = element_blank()) +  # Remove legend key background
  scale_color_manual(values = c("GWSA (MGWR)" = "blue", 
                                "GWSA (RF)" = "red", 
                                "GWSA (Origional)" = "green")) +
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  annotate("text", x = as.Date("2010-01-01"), y = max(mean_data$Mean_Value, na.rm = TRUE) - 0.1, 
           label = paste("R² (MGWR, Origional) =", round(r2_modis1_gwsa, 3)), 
           color = "blue", 
           hjust = -0.60,
           size = 5) +  # Align R² text to the left
  annotate("text", x = as.Date("2010-01-01"), y = max(mean_data$Mean_Value, na.rm = TRUE) - 15, 
           label = paste("R² (RF,Origional) =", round(r2_modis2_gwsa, 3)), 
           color = "red", 
           hjust = -0.72,
           size = 5)  # Align R² text to the left

a

# Save the plot
ggsave(filename = paste0(base_path, "\\Mean_GWSA_and_MODIS_Values.jpeg"), plot = a, dpi = 100, width = 10, height = 6)


#####################  Mean monthly GWSA for Beandenburg

# Load necessary libraries
library(raster)
library(rgdal)
library(dplyr)
library(ggplot2)
library(tidyr)  # For pivoting data

# Define the paths
raster_path_1 <- paste0(base_path, "\\downscaled_output")
raster_path_2 <-  paste0(base_path, "\\rf")
raster_path_3 <-  paste0(base_path, "\\GWSA\\resampled_1km")
shapefile_path <-  paste0(base_path, "\\Germany_shapefile\\Brandenburg.shp")

# Read the shapefile
germany_shape_longlat <- readOGR(shapefile_path)
new_crs <- CRS("+proj=longlat +datum=WGS84 +no_defs")
germany_shape <- spTransform(germany_shape_longlat, new_crs)

# List the raster files
raster_files_1 <- list.files(raster_path_1, pattern = "GWR_\\d{4}-\\d{2}\\.tiff$", full.names = TRUE)
raster_files_2 <- list.files(raster_path_2, pattern = "GWR_\\d{4}-\\d{2}\\.tiff$", full.names = TRUE)
raster_files_3 <- list.files(raster_path_3, pattern = "GWSA_monthly_\\d{4}-\\d{2}\\.tif$", full.names = TRUE)

# Function to calculate mean raster value masked by the shapefile
calculate_mean_raster <- function(raster_file, mask) {
  r <- raster(raster_file)
  masked_r <- mask(r, mask)
  mean_value <- cellStats(masked_r, stat = 'mean', na.rm = TRUE)
  return(mean_value)
}

# Calculate mean values for each raster in each directory
mean_values_1 <- sapply(raster_files_1, calculate_mean_raster, mask = germany_shape)
mean_values_2 <- sapply(raster_files_2, calculate_mean_raster, mask = germany_shape)
mean_values_3 <- sapply(raster_files_3, calculate_mean_raster, mask = germany_shape)

# Extract dates from raster file names
dates_1 <- sub(".*_(\\d{4}-\\d{2})\\.tiff$", "\\1-01", raster_files_1)  # Appending '-01' to represent the first of the month
dates_2 <- sub(".*_(\\d{4}-\\d{2})\\.tiff$", "\\1-01", raster_files_2)
dates_3 <- sub(".*_(\\d{4}-\\d{2})\\.tif$", "\\1-01", raster_files_3)

# Combine into a data frame
mean_data <- data.frame(
  Date = as.Date(c(dates_1, dates_2, dates_3), format = "%Y-%m-%d"),  # Ensure the format matches
  Mean_Value = c(mean_values_1, mean_values_2, mean_values_3),
  Source = c(rep("GWSA (MGWR)", length(mean_values_1)), 
             rep("GWSA (RF)", length(mean_values_2)), 
             rep("GWSA (Origional)", length(mean_values_3)))
)

head(mean_data)

# Reshape mean_data to wide format
wide_data <- mean_data %>%
  pivot_wider(names_from = Source, values_from = Mean_Value)

# Calculate R² for MODIS 1 vs GWSA
r2_modis1_gwsa <- cor(wide_data$`GWSA (Origional)`, wide_data$`GWSA (MGWR)`, use = "complete.obs")^2

# Calculate R² for MODIS 2 vs GWSA
r2_modis2_gwsa <- cor(wide_data$`GWSA (RF)`, wide_data$`GWSA (Origional)`, use = "complete.obs")^2

a <- ggplot(mean_data, aes(x = Date, y = Mean_Value, color = Source)) +
  geom_line(size = 0.8) +  # Increase line size for better visibility
  labs(title = "", 
       x = "Date", 
       y = "GWSA (mm)") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.15),  # Position legend inside the plot
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.key = element_blank()) +  # Remove legend key background
  scale_color_manual(values = c("GWSA (MGWR)" = "blue", 
                                "GWSA (RF)" = "red", 
                                "GWSA (Origional)" = "green")) +
  scale_x_date(date_breaks = "1 years", date_labels = "%Y") +
  annotate("text", x = as.Date("2015-01-01"), y = max(mean_data$Mean_Value, na.rm = TRUE) - 0.1, 
           label = paste("R² (MGWR, Origional) =", round(r2_modis1_gwsa, 3)), 
           color = "blue", 
           hjust = -0.60,
           size = 5) +  # Align R² text to the left
  annotate("text", x = as.Date("2015-01-01"), y = max(mean_data$Mean_Value, na.rm = TRUE) - 15, 
           label = paste("R² (RF,Origional) =", round(r2_modis2_gwsa, 3)), 
           color = "red", 
           hjust = -0.72,
           size = 5)  # Align R² text to the left

a

# Save the plot
ggsave(filename = paste0(base_path, "\\Mean_GWSA_and_brandenburg_Values.jpeg"), 
       plot = a, dpi = 100, width = 10, height = 6)


#####################  Mean monthly GWSA vs GWD for Beandenburg 

library(ggplot2)
library(dplyr)

data <- read.csv(paste0(base_path, "\\clean_2016-22.csv"))

formats <- c("ym", "my", "my", "my")
data$Date <- parse_date_time(data$Date, orders = formats)

data <- data %>%
  dplyr::mutate(YearMonth = format(Date, "%Y-%m"))

# Step 3: Group by Year-Month, then calculate mean GWD across all points
monthly_mean_gwd <- data %>%
  dplyr::group_by(YearMonth) %>%
  dplyr::summarise(mean_GWD = mean(GWD, na.rm = TRUE))

monthly_mean_gwd$YearMonth <- as.Date(paste0(monthly_mean_gwd$YearMonth, "-01"), format = "%Y-%m-%d")


#### scaling GWD, MGWR, RF and Original GW data

monthly_mean_gwd$mean_GWD <-  (monthly_mean_gwd$mean_GWD - min(monthly_mean_gwd$mean_GWD)) / 
  (max(monthly_mean_gwd$mean_GWD) - min(monthly_mean_gwd$mean_GWD)) * 2 - 1

monthly_mean_gwd$mean_GWD <- monthly_mean_gwd$mean_GWD * -1

mean_data$Mean_Value <- (mean_data$Mean_Value - min(mean_data$Mean_Value)) / 
  (max(mean_data$Mean_Value) - min(mean_data$Mean_Value)) * 2 - 1

head(mean_data)

#### combining the GWD and MGWR, RG, Original GW data

monthly_mean_gwd <- monthly_mean_gwd %>%
  dplyr::rename(Date = YearMonth)

combined_data <- merge(monthly_mean_gwd, mean_data, by = "Date", all = TRUE)


combined_data_wide <- combined_data %>%
  pivot_wider(names_from = Source, values_from = Mean_Value, 
              names_prefix = "Mean_Value_") 

combined_data_wide_renamed <- combined_data_wide %>%
  dplyr::rename(
    Date = Date,
    GWD_Mean = mean_GWD,
    NA_Mean_Value = Mean_Value_NA,
    MGWR_Mean_Value = `Mean_Value_GWSA (MGWR)`,
    RF_Mean_Value = `Mean_Value_GWSA (RF)`,
    Original_Mean_Value = `Mean_Value_GWSA (Origional)`
  )

combined_data_wide_unique <- combined_data_wide_renamed %>%
  dplyr::group_by(Date) %>%
  dplyr::summarise(
    GWD_Mean = coalesce(first(GWD_Mean), last(GWD_Mean)),
    NA_Mean_Value = coalesce(first(NA_Mean_Value), last(NA_Mean_Value)),
    MGWR_Mean_Value = coalesce(first(MGWR_Mean_Value), last(MGWR_Mean_Value)),
    RF_Mean_Value = coalesce(first(RF_Mean_Value), last(RF_Mean_Value)),
    Original_Mean_Value = coalesce(first(Original_Mean_Value), last(Original_Mean_Value)),
    .groups = 'drop'
  )

combined_data_wide_final <- combined_data_wide_unique %>%
  select(-NA_Mean_Value)

final_data <- combined_data_wide_final %>%
  drop_na()


#### Create the ggplot

# Reshape the data for plotting
final_data_long <- final_data %>%
  pivot_longer(cols = c(MGWR_Mean_Value, RF_Mean_Value, Original_Mean_Value),
               names_to = "Variable",
               values_to = "Value")

# Function to calculate R² for the specified variable
calculate_r_squared <- function(variable, data) {
  model <- lm(GWD_Mean ~ data[[variable]], data = data)
  r_squared <- summary(model)$r.squared
  return(r_squared)
}

# Calculate R² values for each variable
r_squared_values <- final_data_long %>%
  dplyr::group_by(Variable) %>%
  summarise(R_squared = calculate_r_squared(unique(Variable), final_data), .groups = 'drop')

# Calculate R² values for plotting
r_squared_values <- tibble(
  Variable = c("MGWR_Mean_Value", "Original_Mean_Value", "RF_Mean_Value"),
  R_squared = c(0.603, 0.615, 0.573)
)

a <- # Create the time series plot
  ggplot(final_data_long, aes(x = Date)) +
  geom_line(aes(y = GWD_Mean, color = "GWD Mean"), size = 0.6) +
  geom_line(aes(y = Value, color = Variable), size = 0.6) +
  labs(title = "", 
       x = "Date", 
       y = "Standarized values") +
  theme_bw() +
  theme(legend.position = c(0.12, 0.15),  # Position legend inside the plot
        legend.title = element_blank(),  # Remove legend title
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 22),
        axis.title.y = element_text(size = 22),
        axis.text.y = element_text(size = 22),
        legend.key = element_blank()) +  # Remove legend key background
  scale_color_manual(values = c("GWD Mean" = "blue", 
                                "MGWR_Mean_Value" = "green", 
                                "RF_Mean_Value" = "red", 
                                "Original_Mean_Value" = "orange"), 
                     labels = c("GWD Mean" = "GWD (Observations)", 
                                "MGWR_Mean_Value" = "MGWR", 
                                "RF_Mean_Value" = "RF", 
                                "Original_Mean_Value" = "Original")) +  # Update colors for all variables
  scale_x_date(date_breaks = "2 years", date_labels = "%Y") +
  ylim(-1,1) +
  annotate("text", x = as.Date("2017-10-01"), y = max(final_data$GWD_Mean, na.rm = TRUE) + 0.3, 
           label = paste("R² (GWD, RF) =", round(r_squared_values$R_squared[1], 3)), 
           color = "red", 
           hjust = -1.2,
           size = 5) +  # Align R² text to the left
  annotate("text", x = as.Date("2017-10-01"), y = max(final_data$GWD_Mean, na.rm = TRUE) + 0.20, 
           label = paste("R² (GWD, Original) =", round(r_squared_values$R_squared[2], 3)), 
           color = "orange", 
           hjust = -0.99,
           size = 5) +  # Align R² text to the left
  annotate("text", x = as.Date("2017-10-01"), y = max(final_data$GWD_Mean, na.rm = TRUE) + 0.1, 
           label = paste("R² (GWD, MGWR) =", round(r_squared_values$R_squared[3], 3)), 
           color = "green", 
           hjust = -1,
           size = 5)

a

# Save the plot
ggsave(filename = paste0(base_path, "\\Mean_GWSA_GWD_and_brandenburg_Values.jpeg"), plot = a, dpi = 100, width = 10, height = 6)


####### create mean yearly raster files for "RF"

# Load required libraries
library(raster)
library(terra)
library(dplyr)

# Define the path to the raster files
input_dir <- paste0(base_path, "\\rf\\")
output_dir <- paste0(base_path, "\\mean_yearly\\")

# List all .tiff files with "GWR" in their name
raster_files <- list.files(input_dir, pattern = "GWR_\\d{4}-\\d{2}.tiff$", full.names = TRUE)

# Extract years from file names
raster_years <- gsub(".*GWR_(\\d{4})-\\d{2}.tiff$", "\\1", raster_files)

# Loop through each year, calculate the mean, and save the result
unique_years <- unique(raster_years)

for (year in unique_years) {
  # Get all files for the specific year
  year_files <- raster_files[grepl(paste0("GWR_", year), raster_files)]
  
  # Load the raster files for that year
  rasters <- stack(year_files)
  
  # Calculate the mean raster
  mean_raster <- calc(rasters, fun = mean, na.rm = TRUE)
  
  # Define the output file name
  output_filename <- paste0(output_dir, "GWR_mean_", year, ".tiff")
  
  # Save the mean raster
  writeRaster(mean_raster, filename = output_filename, format = "GTiff", overwrite = TRUE)
  
  print(paste("Saved mean raster for year", year, "as", output_filename))
}

