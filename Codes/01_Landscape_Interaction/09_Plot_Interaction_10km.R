data_dir <- "data/"

library(dplyr)
library(ggplot2)
library(terra)
library(ggspatial)
library(sf)

Total_rast_10km <- rast(paste0(data_dir, "Summary_AridCutHyper/Results/1_Data/10kmProject/LIS_10km.tif"))
roi <- read_sf(paste0(data_dir, "Image/Boundary/Global_land.shp"))
roi_fixed <- st_wrap_dateline(roi, options = c("WRAPDATELINE=YES", "DATELINEOFFSET=180"))
clip_roi <- vect(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion_cutHyperArid.shp"))
dryland_roi <- read_sf(paste0(data_dir, "Image/DrylandRegion/GlobalAridRegion.shp"))
plot_img <- Total_rast_10km %>% crop(clip_roi) %>% mask(., clip_roi)

{
    plot_df <- as.data.frame(plot_img, xy = TRUE)
    colnames(plot_df) <- c("longitude", "latitude", "LIS")
    meanLIS <- mean(plot_df$LIS, na.rm = TRUE)
    p <- ggplot(plot_df[which(plot_df$latitude < 58), ], aes(latitude, LIS)) + geom_hline(yintercept = meanLIS, 
        linetype = "dashed", linewidth = 0.5, color = "black") + geom_smooth(method = "gam", size = 1, 
        color = "red", fill = "red", alpha = 0.3) + scale_color_gradient(low = "gray85", high = "grey20") + 
        guides(color = "none") + scale_x_continuous(position = "top", breaks = seq(-40, 60, 20),
        labels = function(x) ifelse(x == 0, "0°", paste0(abs(x), "°", ifelse(x > 0, "N", "S"))), limits = c(-54, 
        80), expand = c(0.02, 0.02)) + labs(y = "Landscape interaction strength", x = "Latitude") + coord_flip() + 
        theme_bw() + theme(plot.title = element_text(size = 10, hjust = 0.5), legend.key.size = unit(0.3, 
        "cm"), legend.text = element_text(size = 8, hjust = 0), legend.title = element_text(size = 8, 
        face = "bold", hjust = 0), legend.background = element_rect(colour = "NA", fill = "NA", size = 0), 
        strip.text.x = element_text(size = 10, face = "bold", vjust = 0.5), strip.text.y = element_text(size = 10, 
            face = "plain"), axis.title.x = element_text(size = 6, vjust = 0.5, color = "black"), axis.title.y = element_text(size = 6, 
            hjust = 0.5, vjust = 0.5, color = "black"), axis.text.x = element_text(size = 6, vjust = 0.5, 
            hjust = 0.5, color = "black"), axis.text.y = element_text(size = 6, color = "black", hjust = 1), 
        axis.ticks = element_line(colour = "black", size = 0.25), axis.ticks.length = unit(0.05, "cm"), 
        panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(), panel.spacing.x = unit(0.1, "mm"), panel.spacing.y = unit(0.1, 
            "mm"), legend.position = "top", legend.spacing.y = unit(0.002, "cm"), strip.background = element_blank())
    ggsave(plot = p, paste0(data_dir, "GlobalLISMap1km/Results/Merged_LIS_map/", "Landscape_Interaction_Strength_by_Latitude.jpeg"), 
        width = 3.2, height = 9, dpi = 900, units = "cm")
}

mycolors <- c("#7A0000", "#D7301F", "#FDAE61", "#F7F6B5", "#B8D9BE", "#93C6B6", "#6BB0AB", "#4A9DA2", 
    "#288B96", "#137689", "#0A617A", "#044F6D", "#013E5E", "#002E4E") %>% rev()

{
    p <- ggplot() + geom_sf(data = roi_fixed, fill = "grey92", color = NA) + geom_sf(data = dryland_roi, 
        fill = "grey45", color = NA) + layer_spatial(data = plot_img) + scale_fill_gradientn(name = "Landscape interaction strength",
        limits = c(0.4, 1), oob = scales::squish, breaks = seq(0.4, 0.9, 0.1),
        labels = c("≤0.4", "0.5", "0.6", "0.7", "0.8", "0.9"), 
        colors = mycolors, na.value = NA, guide = guide_colorbar(ticks = TRUE, ticks.colour = "black", 
            ticks.linewidth = 0.4, frame.colour = "black", frame.linewidth = 0.2, nbin = 50, draw.ulim = TRUE, 
            draw.llim = TRUE)) + labs(x = NULL, y = NULL) + scale_x_continuous(expand = c(0.03, 0.03), 
        limits = c(-154, 168)) + scale_y_continuous(expand = expansion(mult = c(0.02, 0.05)), limits = c(-56, 
        75.5)) + theme_bw() + theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"), 
        text = element_text(family = "sans"), legend.text = element_text(size = 10, hjust = 0.5, vjust = 3), 
        legend.title = element_text(size = 12, face = "plain", hjust = 0.5), legend.key.width = unit(1.2, 
            "cm"), legend.key.height = unit(0.25, "cm"), legend.background = element_rect(colour = "NA", 
            fill = "NA", linewidth = 0), legend.direction = "horizontal", legend.margin = margin(t = -12, 
            unit = "pt"), strip.background = element_rect(colour = "NA", fill = "NA", linewidth = 0), 
        strip.text.x = element_text(size = 10, face = "bold"), strip.text.y = element_text(size = 10, 
            face = "bold"), axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(), 
        plot.margin = unit(c(0.03, 0.03, 0.03, 0.03), "cm"), panel.border = element_rect(fill = NA, linewidth = 0.9), 
        panel.spacing.x = unit(2, "mm"), panel.spacing.y = unit(1, "mm"), legend.position = c(0.572, 
            0.01), )
    ggsave(plot = p, paste0(data_dir, "GlobalLISMap1km/Results/Merged_LIS_map/", "Global_Landscape_Interaction_Strength.jpeg"), 
        width = 18, height = 8, dpi = 900, units = "cm")
}
