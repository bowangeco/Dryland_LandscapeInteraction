data_dir <- "data/"

library(multimode)
library(diptest)
library(dplyr)
library(data.table)

Imgpath <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/1_Preprocessing/Rdata/")
FileList <- list.files(Imgpath, pattern = paste0("^NDVIValuesListSplit_ID.*\\.rds$"), full.names = TRUE)
Nrow_list <- sub("^.*NDVIValuesListSplit_ID.*Nrow(.*)\\.rds$", "\\1", FileList) %>% as.numeric()
FileList2 <- FileList[order(Nrow_list, decreasing = TRUE)]
result_path <- paste0(data_dir, "Coexistence_GlobalAridCutHyper_10km/Results/2_Bistable_test_Temp/")
j <- 0

for (pathi in FileList2) {
    ID <- sub("^.*NDVIValuesListSplit_ID(.*)_.*Nrow.*\\.rds$", "\\1", pathi) %>% as.numeric()
    id <- sub("^.*NDVIValuesListSplit_ID.*_(.*)Nrow.*\\.rds$", "\\1", pathi) %>% as.numeric()
    bimodal_test <- function(listi) {
        lon <- listi[1]
        lat <- listi[2]
        LIS <- listi[3]
        NDVI_values <- listi[-1:-3]
        NDVI_values_narm <- NDVI_values[!is.na(NDVI_values)]
        NDVI_values_inRange <- NDVI_values_narm[which(NDVI_values_narm > 0 & NDVI_values_narm <= 0.4)]
        row <- data.frame(x = lon, y = lat, LIS = LIS)
        if (length(NDVI_values_inRange) < 10 | length(NDVI_values_narm) < 30) {
            row[c("numb", "dip_p", "ACR_p", "mod1mean", "mod1height", "mod2mean1", "mod2mean2", "mod2height1", 
                "mod2height2", "mod2valley", "mod2valleyheight")] <- NA
        }
        else {
            numb <- length(NDVI_values_narm)
            dip_p <- dip.test(NDVI_values_narm)$p
            if (length(NDVI_values_narm) > 10000) {
                ACR_p <- multimode::modetest(sample(NDVI_values_narm, 10000), method = "ACR")$p.value
            }
            else {
                ACR_p <- multimode::modetest(NDVI_values_narm, method = "ACR")$p.value
            }
            mod1 <- locmodes(data = NDVI_values_narm, mod0 = 1, lowsup = 0, uppsup = 0.4)
            mod2 <- locmodes(data = NDVI_values_narm, mod0 = 2, lowsup = 0, uppsup = 0.4)
            mod1mean <- mod1$locations
            mod1height <- mod1$fvalue
            mod2mean1 <- mod2$locations[1]
            mod2mean2 <- mod2$locations[3]
            mod2height1 <- mod2$fvalue[1]
            mod2height2 <- mod2$fvalue[3]
            mod2valley <- mod2$locations[2]
            mod2valleyheight <- mod2$fvalue[2]
            result <- c(numb, dip_p, ACR_p, mod1mean, mod1height, mod2mean1, mod2mean2, mod2height1, 
                mod2height2, mod2valley, mod2valleyheight)
            row[c("numb", "dip_p", "ACR_p", "mod1mean", "mod1height", "mod2mean1", "mod2mean2", "mod2height1", 
                "mod2height2", "mod2valley", "mod2valleyheight")] <- result
        }
        return(row)
    }
    df_list <- readRDS(pathi)
    Total_result <- lapply(df_list, FUN = bimodal_test) %>% do.call(rbind, .)
    fwrite(Total_result, paste0(result_path, "Coexistence_ID", ID, "_id", id, "_10km.csv"), row.names = FALSE)
}
