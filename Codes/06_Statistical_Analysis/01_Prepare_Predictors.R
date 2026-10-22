data_dir <- "data/"

library(dplyr)
library(data.table)
library(terra)

LIS_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif")
LIS <- rast(LIS_path)
Pre_path <- paste0(data_dir, "Image/Global_Factors/10km/Pre_mean_10km_Global.tif")
Tem_path <- paste0(data_dir, "Image/Global_Factors/10km/Tem_mean_10km_Global.tif")
NDVI_path <- paste0(data_dir, "Image/Global_Factors/10km/NDVI_mean_10km_Global.tif")
Wind_path <- paste0(data_dir, "Image/Global_Factors/10km/WindSpeed_mean_10km_Global.tif")
Human_path <- paste0(data_dir, "Image/Global_Factors/10km/Human_10km_Global.tif")
DEM_path <- paste0(data_dir, "Image/Global_Factors/10km/DEM_mean_10km_Global.tif")
DEM_sd_path <- paste0(data_dir, "Image/Global_Factors/10km/DEM_sd_10km_Global.tif")
sand_path <- paste0(data_dir, "Image/Global_Factors/10km/sand_mean_10km_Global.tif")
sand_sd_path <- paste0(data_dir, "Image/Global_Factors/10km/sand_sd_10km_Global.tif")
water_path <- paste0(data_dir, "Image/Global_Factors/10km/water_mean_10km_Global.tif")
SOC_path <- paste0(data_dir, "Image/Global_Factors/10km/SOC_mean_10km_Global.tif")
Pre <- rast(Pre_path)

Tem <- rast(Tem_path)
NDVI <- rast(NDVI_path)
Wind <- rast(Wind_path)
Human <- rast(Human_path)
DEM <- rast(DEM_path)
DEM_sd <- rast(DEM_sd_path)
sand <- rast(sand_path)
sand_sd <- rast(sand_sd_path)
water <- rast(water_path)
SOC <- rast(SOC_path)
out_path <- paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
mean_img_list <- list(Pre, Tem, NDVI, Wind, Human, DEM, sand, water, SOC)

img_names <- c("Pre", "Tem", "NDVI", "Wind", "Human", "DEM", "sand", "water", "SOC")

for (i in seq_along(mean_img_list)) {
    image <- crop(mean_img_list[[i]], LIS) %>% extend(., LIS) %>% terra::project(., LIS, method = "average") %>% 
        mask(., LIS)
    writeRaster(image, paste0(out_path, img_names[i], "_10km.tif"), overwrite = TRUE)
}

near_img_list <- list(DEM_sd, sand_sd)
img_names <- c("DEM_sd", "sand_sd")

for (i in seq_along(near_img_list)) {
    image <- crop(near_img_list[[i]], LIS) %>% extend(., LIS) %>% terra::project(., LIS, method = "near") %>% 
        mask(., LIS)
    writeRaster(image, paste0(out_path, img_names[i], "_10km.tif"), overwrite = TRUE)
}

GroundWater <- vect(paste0(data_dir, "Image/Global_Factors/UndergroundWater_Global.shp"))
GroundWater_rast <- rasterize(GroundWater, LIS, "GmtcD50")
writeRaster(GroundWater_rast, paste0(out_path, "GroundWater.tif"), overwrite = TRUE)

Ecoregion <- vect(paste0(data_dir, "Image/Global_Factors/10km/Ecoregion_2017.shp"))
Biome2017_rast <- rasterize(Ecoregion, LIS, "BIOME_NAME")
writeRaster(Biome2017_rast, paste0(out_path, "Biome2017.tif"), overwrite = TRUE)

# MSI inputs follow the original nearest-neighbor alignment.
msi_files <- c("MSI_Precipitation_Mean.tif", "MSI_Temperature_Mean.tif", "MSI_WindSpeed_Mean.tif")
msi_names <- c("MSI_Pre", "MSI_Tem", "MSI_Wind")
for (i in seq_along(msi_files)) {
    image <- rast(paste0(data_dir, "Image/Global_Factors/10km/", msi_files[i]))
    image <- terra::project(image, LIS, method = "near") %>% mask(., LIS)
    writeRaster(image, paste0(out_path, msi_names[i], "_10km.tif"), overwrite = TRUE)
}

# Annual climate and 19 five-year windows.
Annual_Tem <- rast(paste0(data_dir, "Image/Global_Factors/10km/Annual_Tem_10km.tif"))
Annual_Pre <- rast(paste0(data_dir, "Image/Global_Factors/10km/Annual_Pre_10km.tif"))
Annual_Wind <- rast(paste0(data_dir, "Image/Global_Factors/10km/Annual_WindSpeed_10km.tif"))
mean_img_list <- list(Annual_Tem, Annual_Pre, Annual_Wind)
img_names <- c("Annual_Tem", "Annual_Pre", "Annual_Wind")

for (i in 1:length(mean_img_list)) {
    namei <- img_names[i]
    imgi <- mean_img_list[[i]]
    Projected_imgi <- crop(imgi, LIS) %>% extend(., LIS) %>% terra::project(., LIS, method = "average") %>% 
        mask(., LIS)
    writeRaster(Projected_imgi, paste0(out_path, namei, "_10km.tif"), overwrite = TRUE)
}

{
    Annual_Tem <- rast(paste0(out_path, "Annual_Tem", "_10km.tif"))
    Annual_Pre <- rast(paste0(out_path, "Annual_Pre", "_10km.tif"))
    Annual_Wind <- rast(paste0(out_path, "Annual_Wind", "_10km.tif"))
    mean_img_list <- list(Annual_Tem, Annual_Pre, Annual_Wind)
    img_names <- c("Annual_Tem", "Annual_Pre", "Annual_Wind")
    for (i in 1:length(mean_img_list)) {
        namei <- img_names[i]
        imgi <- mean_img_list[[i]]
        mw_img <- lapply(1:19, function(year) {
            imgs_5year <- imgi[[year:(year + 4)]]
            mean_img <- mean(imgs_5year, na.rm = TRUE)
            return(mean_img)
        }) %>% rast()
        writeRaster(mw_img, paste0(out_path, namei, "_MovingWindows_10km.tif"), 
            overwrite = TRUE)
    }
}
