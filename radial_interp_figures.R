## Title: Radial Interpolation Methods — FIGURES, ANIMATIONS & STATS
## Created by: Asha Kannan
## Created on: 4-8-2026
## Split: all plotting/animation/stats code. Run 01_model_building.R first
## (or source it) so that `results`, `transects`, `cover_fractions`,
## `build_veg_idw_raster()`, `master_template`, `weight_max`, and the
## tower/transect constants are available in the environment.

# Uncomment to rerun the model-building pipeline first:
# source("01_model_building.R")

# =============================================================
# Section 8: Footprint Video
# =============================================================
# Purpose: generate an animated .mp4 showing how the flux footprint
# source area changes through time (2014-2018), with transect lines
# overlaid to show sampling locations relative to the footprint.
#
#
# Output: footprint_animation_final.mp4
# =============================================================
# create fresh output directory for animation frames
frame_dir <- "frames_fpt"
unlink(frame_dir, recursive = TRUE)
dir.create(frame_dir, showWarnings = FALSE)

# tighter plot bounds - zoom in around the footprint
x_min_plot <- -106.2400
x_max_plot <- -106.2370
y_min_plot <-  34.4375
y_max_plot <-  34.4395

cat("Generating", nrow(results), "animation frames...\n")

# ... run just the single frame plot code with new bounds
for (i in 1:nrow(results)) {
  cat("frame", i, "of", nrow(results), "\n"); flush.console()
  
  tryCatch({
    # load and resample footprint to master grid
    # resampling ensures all frames have identical spatial extent
    fpt <- rast(results$filepath[i])
    fpt <- resample(fpt, master_template, method = "near")
    
    # convert raster to data frame for ggplot
    df <- as.data.frame(fpt, xy = TRUE)
    colnames(df) <- c("x", "y", "weight")
    
    # frame label for plot title
    label <- paste0(results$year[i], "-",
                    sprintf("%02d", results$month[i]), "-",
                    results$period[i])
    
    # build frame plot
    p <- ggplot(df, aes(x = x, y = y, z = weight)) +
      
      # white background within plot bounds
      annotate("rect",
               xmin = x_min_plot, xmax = x_max_plot,
               ymin = y_min_plot, ymax = y_max_plot,
               fill = "white") +
      # footprint contours colored by weight level
      geom_contour(breaks = seq(0, weight_max, length.out = 10),
                   aes(color = after_stat(level)), linewidth = 0.8) +
      
      # transect lines showing field sampling locations
      # dashed black lines radiating from tower at 180, 240, 315 degrees
      geom_segment(data = transect_lines,
                   aes(x = x_start, y = y_start,
                       xend = x_end, yend = y_end),
                   inherit.aes = FALSE,
                   color = "black", linewidth = 1, linetype = "dashed") +
      # flux tower location
      geom_point(data = tower_df,
                 aes(x = x, y = y),
                 inherit.aes = FALSE,
                 color = "black", size = 6, shape = 17) +
      # fix aspect ratio and clip to consistent bounds across frames
      coord_equal(xlim = c(x_min_plot, x_max_plot),
                  ylim = c(y_min_plot, y_max_plot),
                  expand = FALSE) +
      # color scale for contour levels
      scale_color_distiller(palette = "YlGn") +
      
      theme_classic(base_size = 14) +
      theme(
        plot.margin  = margin(t = 20, r = 10, b = 10, l = 10),
        legend.position = "none"
      ) +
      labs(title = label, x = "Longitude", y = "Latitude")
    
    # save frame
    fname <- file.path(frame_dir, sprintf("frame_%04d.png", i))
    ggsave(fname, p, width = 8, height = 6, dpi = 150)
    
    rm(fpt, df, p); gc()
    
  }, error = function(e) {
    cat("Error on frame", i, ":", e$message, "\n")
  })
}

# encode all frames into .mp4
local_output <- file.path(tempdir(), "footprint_animation_final.mp4")
av_encode_video(png_files,
                output = local_output,
                framerate = 6)

png_files <- sort(list.files(frame_dir,pattern = "frame_.*\\.png",full.names = TRUE))
av_encode_video(png_files,output = "footprint_animation_final.mp4",framerate = 6)
cat("Animation saved to footprint_animation_4_24_26.mp4\n")

# =============================================================
# Section 9: Manuscript Figure 2
# =============================================================
# Purpose: overlay contour lines on vegetation raster
#
#
# Output: Heatmap figure showing plant group relative abundance/distribution
# =============================================================
library(ggplot2)
library(patchwork)
library(terra)
library(dplyr)
library(tidyr)

# 1. compute mean annual day footprint (DAY not night)
mean_annual_day <- mean_footprint_annual <- {
  files <- results %>% 
    filter(period == "DAY") %>% 
    pull(filepath)
  
  rast_list <- lapply(files, rast)
  ref <- rast_list[[1]]
  
  resampled <- lapply(rast_list, function(r) {
    if (!all(ext(r) == ext(ref)) | !all(res(r) == res(ref))) {
      resample(r, ref, method="bilinear")
    } else {
      r
    }
  })
  stack <- rast(resampled)
  mean_fp <- mean(stack, na.rm=TRUE)
  mean_fp / global(mean_fp, "sum", na.rm=TRUE)[[1]]
}

# 2. build veg raster on mean day footprint
# expand the mean_annual_day extent before building veg raster
e <- ext(mean_annual_day)
e_expanded <- ext(
  e[1] - 0.2 * (e[2] - e[1]),
  e[2] + 0.2 * (e[2] - e[1]),
  e[3] - 0.2 * (e[4] - e[3]),
  e[4] + 0.2 * (e[4] - e[3])
)

fpt_expanded <- rast(
  ext        = e_expanded,
  resolution = res(mean_annual_day),
  crs        = crs(mean_annual_day)
)

# resample mean_annual_day onto the expanded grid so pixel centers match veg_f2
fpt_for_plot <- resample(mean_annual_day, fpt_expanded, method = "bilinear")
# renormalize after resampling so weights still sum to 1
fpt_for_plot <- fpt_for_plot / global(fpt_for_plot, "sum", na.rm = TRUE)[[1]]

# build vegetation raster
veg_f2 <- build_veg_idw_raster(fpt_expanded)

## collapsed biocrust layer
biocrust <- max(veg_f2[["cover_cyano"]],
                veg_f2[["cover_lich"]],
                veg_f2[["cover_moss"]],
                na.rm=TRUE)
names(biocrust) <- "cover_biocrust"

# 3. compute 50% and 80% contour thresholds
fpt_vals <- sort(values(fpt_for_plot), decreasing = TRUE)
fpt_vals <- fpt_vals[!is.na(fpt_vals)]
fpt_cum  <- cumsum(fpt_vals) / sum(fpt_vals)
level_50 <- fpt_vals[which.min(abs(fpt_cum - 0.50))]
level_80 <- fpt_vals[which.min(abs(fpt_cum - 0.80))]
level_90 <- fpt_vals[[which.min(abs(fpt_cum - 0.90))]]

# 4. keep only the cells inside the footprint area (to rid of radial sampling artifact)
mask_level <- level_90
fp_mask <- ifel(fpt_for_plot >= mask_level, 1, NA)

# mask all veg layers 
veg_f2 <- mask(veg_f2, fp_mask)

# crop everything to the footprint's bounded boxed area so that we don't plot empty white
fp_box <- trim(fp_mask) # drops all NA outcomes
veg_f2 <- crop(veg_f2, fp_box)
fpt_for_plot <- crop(fpt_for_plot, fp_box)

# collapsed bsc layer built from masked layers
biocrust <- max(veg_f2[["cover_cyano"]],
                veg_f2[["cover_lich"]],
                veg_f2[["cover_moss"]],
                na.rm=TRUE)
names(biocrust) <- "cover_biocrust"

# 5. compute aspect ratio from the expanded extent to use in theme and replace coord_fixed()
e2 <- ext(fpt_for_plot)
aspect <- as.numeric((e2[4] - e2[3]) / (e2[2] - e2[1]))

# 6. convert raster layer and footprint to plot df
make_plot_df <- function(layer, fpt) {
  combined <- c(layer, fpt)
  df <- as.data.frame(combined, xy=T)
  colnames(df) <- c("lon", "lat", "cover", "weight")
  return(df)
}

# 7. make one panel
make_panel <- function(layer, fpt, title_str,
                       show_x = FALSE, show_y = FALSE,
                       show_legend = FALSE,
                       fill_limits = c(0, 0.35),
                       xlab=NULL, ylab=NULL) {
  df <- make_plot_df(layer, fpt)
  
  p <- ggplot(df, aes(x=lon, y=lat)) +
    geom_raster(aes(fill=cover)) +
    geom_segment(data = transect_lines,
             aes(x = x_start, y = y_start,
                 xend = x_end, yend = y_end),
             inherit.aes = FALSE,
             color = "black", linewidth = 0.75, linetype = "dashed") +
    geom_point(data=data.frame(lon = tower_lon, lat=tower_lat),
               aes(x=lon, y=lat),
               color="black", size=3, shape=17,
               inherit.aes=FALSE) +
    geom_contour(aes(z=weight),
                 breaks=c(level_50, level_80),
                 color="black", linewidth=0.3) +
    scale_fill_distiller(palette = "Greens", direction=1,
                         limits=fill_limits,
                         na.value = "white",
                         name="Cover\nfraction") +
    theme_classic(base_size=12) +
    labs(title = title_str, x = xlab, y = ylab) +
    theme(
      aspect.ratio = aspect,
      plot.title = element_text(size=15, face="bold", hjust=0.5),
      legend.position = if (show_legend) "right" else "none",
      axis.text = element_text(size=9),
      axis.ticks = element_line(),
      legend.title = element_text(size=9), 
      plot.margin = margin(0,0,0,0)
    )
  
  if (!show_x) p <- p + theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                              axis.ticks.y = element_blank())
  return(p)
}

# 8. build all 6 panels
## consistent fill limits across all panels
fill_lims <- c(0, 0.35)

## top row
# shared scale
fill_lims_top <- c(0, 1.0)
p_biocrust <- make_panel(biocrust, fpt_for_plot, "Biocrust", show_y = TRUE, ylab = "Latitude (°)",fill_limits = fill_lims_top)
p_midveg <- make_panel(veg_f2[["cover_midveg"]], fpt_for_plot, "Midstory", fill_limits=fill_lims_top)
p_canopy <- make_panel(veg_f2[["cover_canopy"]], fpt_for_plot, "Canopy", show_legend = TRUE, fill_limits = fill_lims_top)

## bottom row
# shared scale 
fill_lims_bot <- c(0, 0.18)
p_cyano = make_panel(veg_f2[["cover_cyano"]], fpt_for_plot, "Cyanobacteria-dominated", show_x = T, show_y = T, ylab = "Latitude (°)", fill_limits = fill_lims_bot)
p_lich = make_panel(veg_f2[["cover_lich"]], fpt_for_plot, "Lichen-dominated", show_x = TRUE, xlab = "Longitude (°)", fill_limits = fill_lims_bot)
p_moss = make_panel(veg_f2[["cover_moss"]], fpt_for_plot, "Moss-dominated", show_x = TRUE, show_legend = TRUE, fill_limits = fill_lims_bot)

# 9. assemble w patchwork
# first, assemble each row separately with own legend
row1 <- (p_biocrust | p_midveg | p_canopy) +
  plot_layout(guides="collect") &
  theme(legend.position = "right", plot.margin = margin(0,0,0,0))

row2 <- (p_cyano | p_lich | p_moss) +
  plot_layout(guides="collect") &
  theme(legend.position = "right", plot.margin = margin(0,0,0,0))

# stack the two rows
library(ggh4x)

combined <- row1 / row2 &
  theme(aspect.ratio = aspect,
        legend.position = "right",
        axis.text = element_text(size=11),
        legend.title = element_text(size=11),
        legend.text = element_text(size=10),
        legend.key.size = unit(1.2, "cm"),
        plot.margin = margin(0, 0, 0, 0))



print(combined)
setwd("/Users/ashakannan/Library/CloudStorage/OneDrive-Personal/Desktop/Kannan_et_al(2026)")
ggsave("KannanFig2_edited.tif", combined, height=6, width=12, dpi=600, compression="lzw")

# Figure 1: Collapsed BSC and other PFTs ----------------------------------
# BSC = max(cyano, lich, moss) at each pixel
# use first footprint file as grid template
# the veg raster is static - LPI was sampled once
# we just need any footprint grid to define the pixel structure
fpt_template <- rast(results$filepath[1])

# expand the template extent by 20% on all sides to avoid clipping (since the 
# raster doesn't expand far enough in the western side). We fix this by rebuilding
# the veg raster using a larger template
e <- ext(fpt_template)
e_expanded <- ext(
  e[1] - 0.2 * (e[2] - e[1]),
  e[2] + 0.2 * (e[2] - e[1]),
  e[3] - 0.2 * (e[4] - e[3]),
  e[4] + 0.2 * (e[4] - e[3])
)

fpt_template_expanded <- rast(
  ext = e_expanded,
  resolution = res(fpt_template),
  crs = crs(fpt_template)
)

veg_idw <- build_veg_idw_raster(fpt_template_expanded)

# collapse BSC types into single BSC layer
# take max of cyano, lich, moss at each pixel
# represents total BSC cover regardless of type
biocrust <- max(veg_idw[["cover_cyano"]],
                veg_idw[["cover_lich"]],
                veg_idw[["cover_moss"]],
                na.rm=T)
names(biocrust) <- "cover_biocrust"

# stack w other functional groups
veg_collapsed <- c(
  biocrust,
  veg_idw[["cover_midveg"]], 
  veg_idw[["cover_canopy"]],
  veg_idw[["cover_abiotic"]]
)

# convert to long df
collapsed_df <- as.data.frame(veg_collapsed, xy=T) %>% 
  pivot_longer(cols = -c(x, y),
               names_to = "group",
               values_to = "cover") %>%
  mutate(group = recode(group,
                        cover_biocrust = "Biocrust",
                        cover_midveg   = "Midstory",
                        cover_canopy   = "Canopy",
                        cover_abiotic  = "Abiotic")) %>%
  mutate(group = factor(group,
                        levels = c("Biocrust", "Midstory", 
                                   "Canopy", "Abiotic")))


fig1 <- ggplot(collapsed_df, aes(x = x, y = y, fill = cover)) +
  geom_raster() +
  geom_point(data = tower_df, aes(x = x, y = y),
             inherit.aes = FALSE, color = "black", size = 5, shape = 17) +
  facet_wrap(~group, ncol = 2) +
  scale_fill_distiller(palette = "Greens", direction = 1,
                       na.value = "white",
                       name = "Cover fraction") +
  coord_fixed() +
  theme_classic(base_size = 14) +
  labs(x = "Longitude", y = "Latitude",
       title = "Footprint-weighted vegetation cover") +
  theme(
    axis.text.y = element_text(size = 30),
    axis.text.x =  element_text(size=20), 
    axis.title   = element_text(size = 40),
    strip.text   = element_text(size = 30, face = "bold"),
    legend.position = "right",
    legend.title=element_text(size=30),
    legend.text = element_text(size=30),
    plot.title = element_text(size=40)
  )

fig1
ggsave("fpt_bsc_grouped.png", fig1, height=20, width=40)


# Figure 2: BSC ONLY static veg raster plot -------------------------------
bsc_stack <- c(
  veg_idw[["cover_cyano"]], 
  veg_idw[["cover_lich"]],
  veg_idw[["cover_moss"]]
)

bsc_df <- as.data.frame(bsc_stack, xy=TRUE) %>% 
  pivot_longer(cols=-c(x,y),
               names_to="group",
               values_to="cover") %>% 
  mutate(group = recode(group,
                        cover_cyano = "Cyanobacteria-dominated",
                        cover_lich = "Lichen-dominated",
                        cover_moss = "Moss-dominated")) %>% 
  mutate(group = factor(group, levels=c("Cyanobacteria-dominated", 
                                        "Lichen-dominated",
                                        "Moss-dominated")))

fig_bsc <- ggplot(bsc_df, aes(x=x, y=y, fill=cover)) +
  geom_raster() +
  geom_point(data=tower_df, aes(x=x, y=y),
             inherit.aes = FALSE, color="black", size=5, shape=17) +
  facet_wrap(~group, ncol=3) +
  scale_fill_distiller(palette="Greens", direction = 1, 
                       na.value="white",
                       name="Cover fraction") +
  coord_fixed() +
  theme_classic(base_size = 14) +
  labs(x = "Longitude", y = "Latitude",
       title = "Footprint-weighted vegetation cover") +
  theme(
    axis.text.y = element_text(size = 30),
    axis.text.x =  element_text(size=20), 
    axis.title   = element_text(size = 40),
    strip.text   = element_text(size = 30, face = "bold"),
    legend.position = "right",
    legend.title=element_text(size=30),
    legend.text = element_text(size=30),
    plot.title = element_text(size=40)
  )
  
fig_bsc 
#setwd("/Users/ashakannan/Library/CloudStorage/OneDrive-Personal/Desktop/thesis/final_figures")
ggsave("/Users/ashakannan/Desktop/bsc_weighted_fractions.png", 
       fig_bsc, height = 30, width = 30, dpi = 100)
ggsave("bsc_weighted_fractions.png", fig_bsc, height = 10, width = 30, 
       dpi = 150, device = "png")
# Convert into nice figure ------------------------------------------------
# keep layers I want
veg_idw <- build_veg_idw_raster(fpt)

veg_idw <- veg_idw[[c("cover_cyano", "cover_lich", "cover_moss", "cover_midveg",
                             "cover_canopy", "cover_abiotic")]]

# convert from raster to df
library(terra)
library(ggplot2)
library(dplyr)
library(tidyr)
idw_circles <- as.data.frame(veg_idw, xy=TRUE)
idw_long <- idw_circles %>% 
  pivot_longer(cols=-c(x,y),
               names_to="group",
               values_to="cover")
# clean labels
idw_long <- idw_long %>% 
  mutate(group=recode(group,
                      cover_cyano = "Cyanobacteria-dominated",
                      cover_lich = "Lichen-dominated",
                      cover_moss = "Moss-dominated",
                      cover_midveg = "Midstory",
                      cover_canopy = "Canopy",
                      cover_abiotic = "Abiotic")) %>% 
  mutate(group=factor(group, 
                      levels=c("Moss-dominated", "Lichen-dominated", "Cyanobacteria-dominated", "Midstory", "Canopy", "Abiotic")))



library(ggplot2)

idw_weighting <- ggplot(idw_long, aes(x = x, y = y, fill = cover)) +
  geom_raster() +
  geom_point(data=tower_df, aes(x=x, y=y), inherit.aes=F, color="black", size=5, shape=17) +
  facet_wrap(~group, ncol = 3) +
  scale_fill_distiller(palette = "Greens", direction=1, na.value = "white") +
  #scale_shape_manual(values=c("Flux tower" = 17)) +
  coord_equal() +
  theme_classic(base_size=14) +
  labs(x="Distance E-W (m)", y="Distance N-S (m)") +
  theme(
    axis.text=element_text(size=16),
    axis.title = element_text(size=16),
    panel.grid = element_blank(),
    strip.text = element_text(size=20, face="bold")
  )
idw_weighting
ggsave("idw_weighting.png", idw_weighting, height=10, width=20, dpi=600)



# MAKING VIDEO of footprint changes ---------------------------------------
library(ggplot2)
library(terra)
library(dplyr)
library(tidyr)
library(av)

template <- rast(results$filepath[1])


# the raster is centered on the tower, so use raster center for the tower marker
# instead of true lat/lon
e <- ext(template)
raster_tower_lon <- (e[1] + e[2]) / 2  # -106.2384
raster_tower_lat <- (e[3] + e[4]) / 2  # 34.43837

tower_df <- data.frame(
  x = raster_tower_lon,
  y = raster_tower_lat
)
# calculate transect end points in lon/lat
# tower is the origin, transects go 100m in each azimuth
transect_lines <- data.frame(
  transect=c("Transect 1", "Transect 2", "Transect 3"),
  azimuth=c(180, 240, 315)
) %>% 
  mutate(
    # mirrow azimuth around N-S axis: 180 stuays 180, 240 becomes 120, etc
    az_mirrored = (360-azimuth) %% 360, 
    # endpoint in m relative to tower
    end_x_m = 100 * sin(azimuth*pi/180),
    end_y_m = 100 * cos(azimuth *pi/180),
    
    # convert to lon/lat
    x_start = tower_lon,
    y_start = tower_lat,
    
    x_end = raster_tower_lon + end_x_m / m_per_deg_lon,
    y_end = raster_tower_lat + end_y_m / m_per_deg_lat
  )

# rebuild master_template
all_extents <- lapply(results$filepath, function(f) ext(rast(f)))
master_ext <- all_extents[[1]]
for (e in all_extents[-1]) {
  master_ext <- union(master_ext, e)
}

master_template <- rast(
  ext = master_ext,
  resolution = res(template),
  crs = crs(template)
)

x_min_raw <- ext(master_template)[1]
x_max_raw <- ext(master_template)[2]
y_min <- ext(master_template)[3]
y_max <- ext(master_template)[4]

# after mirroring x, the bounds flip
x_min <- 2 * raster_tower_lon - x_max_raw
x_max <- 2 * raster_tower_lon - x_min_raw

frame_dir <- "frames_fpt"
unlink("frames_fpt", recursive=TRUE)
dir.create(frame_dir, showWarnings=FALSE)

# get consistent color scale across all frames
# sample max weight value across files
weight_max <- max(sapply(results$filepath, function(f) {
  global(rast(f), "max", na.rm=TRUE)[[1]]
}))

results <- results %>% 
  arrange(year, month, period)

for (i in 1:nrow(results)) {
  cat("frame", i, "of", nrow(results), "\n"); flush.console()
  
  tryCatch({
    fpt <- rast(results$filepath[i])
    fpt <- resample(fpt, master_template, method="near")
    
    df <- as.data.frame(fpt, xy=T)
    colnames(df) <- c("x", "y", "weight")
    
    
    # mirror x around the tower to flip source/upwind side
    df$x <- 2*tower_lon - df$x
    
    label <- paste0(results$year[i], "-",
                    sprintf("%02d", results$month[i]), "-",
                    results$period[i])
    
    # start contours at a small threshold to cut off the near-zero noise
    #low_threshold <- weight_max * 0.02 # only show top 98% of weight range
    
    p <- ggplot(df, aes(x=x, y=y, z=weight)) +
      annotate("rect", xmin=x_min, xmax=x_max, ymin=y_min, ymax=y_max, fill="white") +
      geom_contour(breaks=seq(0, weight_max, length.out=10), 
                    aes(color=after_stat(level)), linewidth=0.8) +
      geom_segment(data=transect_lines,aes(x=x_start, y=y_start, xend=x_end, yend=y_end),
                   inherit.aes=FALSE, color="black", linewidth=1, linetype="dashed") +
      geom_point(data=tower_df, aes(x=x, y=y), inherit.aes = FALSE,
                 color="black", size=6, shape=17) +
      coord_equal(xlim=c(x_min, x_max), ylim=c(y_min, y_max), expand=FALSE) +
      scale_color_distiller(palette="YlGn") +
      theme_classic(base_size=14) +
      theme(plot.margin=margin(t=20, r=10, b=10, l=10),
            legend.position="none") +
      labs(title=label, x="Longitude", y="Latitude")
    
    fname <- file.path(frame_dir, sprintf("frame_%04d.png", i))
    ggsave(fname, p, width=8, height=6, dpi=150)
    
    rm(fpt, df, p)
    gc()
  }, error=function(e) {
    cat("error on frame", i, ":", e$message, "\n")
  })
}

png_files <- sort(list.files(frame_dir, pattern="frame_.*\\.png", full.names=T))
av_encode_video(png_files, output="footprint_animation_mirror.mp4", framerate=6)
cat("done\n")

fpt_summary <- read_csv("/Users/ashakannan/Downloads/All_site_fpt_summary.csv")
# find US-Mpj row
fpt_summary %>% filter(grepl("US-Mpj", site_id))


# Wind rose ---------------------------------------------------------------
library(openair)

# load all gapfilled files
flux_dir <- "/Users/ashakannan/Library/CloudStorage/OneDrive-Personal/Desktop/thesis/MPJ_fluxes/NMEG"

flux_files <- list.files(flux_dir, pattern = "US-Mpj_.*_gapfilled.txt", full.names = TRUE)

flux_all <- bind_rows(lapply(flux_files, function(f) {
  read.table(f, header = TRUE, sep = ",", na.strings = c("-9999", "NA"))
}))

# filter to 2018 and valid wind data
flux_2018 <- flux_all %>%
  filter(YEAR == 2018, !is.na(WD), !is.na(WS), WS > 0)

# wind rose
windRose(flux_2018, 
         ws = "WS", 
         wd = "WD",
         main = "US-Mpj 2018 Wind Rose")

# MAKING VIDEO of crust weight fluctuations (ESA) ------------------------------------------------------------
library(dplyr)
library(tidyr)
library(gganimate)
library(av)
library(ggplot2)
frame_dir <- "frames_ts"
dir.create(frame_dir, showWarnings = FALSE)

# prepare time series data
ts_df <- results %>% 
  arrange(year, month, period) %>% 
  mutate(
    frame_idx = row_number(),
    label=paste0(year, "-", sprintf("%02d", month), "-", period)
  ) %>% 
  select(frame_idx, label, year, month, period, w_cyano, w_lich, w_moss) %>% 
  pivot_longer(cols=c(w_cyano, w_lich, w_moss),
               names_to = "crust",
               values_to = "cover") %>% 
  mutate(crust=recode(crust,
                      w_cyano="Cyanobacteria-dominated",
                      w_lich = "Lichen-dominated",
                      w_moss = "Moss-dominated"))

ymax <- max(ts_df$cover, na.rm=T) * 1.05
n_frames <- max(ts_df$frame_idx)

# add year boundary x positions for vertical lines
year_starts <- results %>% 
  arrange(year, month, period) %>% 
  mutate(frame_idx = row_number()) %>% 
  group_by(year) %>% 
  summarise(frame_idx = min(frame_idx))

for (i in 1:n_frames) {
  cat("Frame", i, "of", n_frames, "\n"); flush.console()
  
  # data up to current frame (trail)
  past <- ts_df %>% filter(frame_idx <= i)
  
  # current frame highlight points
  current <- ts_df %>% filter(frame_idx == i)
  
  p <- ggplot(past, aes(x=frame_idx, y=cover, color=crust)) +
    # year boundary lines
    geom_vline(data=year_starts, aes(xintercept=frame_idx), 
               color="grey80", linetype="dashed", linewidth=0.5, inherit.aes = F) +
    # year labels at top
    geom_text(data=year_starts, aes(x=frame_idx, label=year), y=y_max, hjust=-0.1, 
              color="grey50", size=3.5, inherit.aes = F) +
    geom_line(linewidth=1, alpha=0.8) +
    geom_point(data=current, aes(x=frame_idx, y=cover, color=crust),
               size=4, inherit.aes = F) +
    scale_color_manual(values=c(
      "Cyanobacteria-dominated" = "#8b7355",
      "Lichen-dominated" = "coral",
      "Moss-dominated" = "darkgreen"
    )) +
    scale_x_continuous(
      breaks = year_starts$frame_idx,
      labels = year_starts$year,
      limits = c(1, n_frames)
    ) +
    scale_y_continuous(limits=c(0, 0.045)) +
    theme_classic(base_size=14) +
    theme(legend.position="right",
          axis.text.x = element_text(size=14)) +
    labs(title=paste0("Footprint-weighted BSC cover: ", current$label[1]),
         x=NULL,
         y="Weighted cover fraction")
  
  fname <- file.path(frame_dir, sprintf("frame_%04d.png", i))
  ggsave(fname, p, width=12, height=6, dpi=150)
  
  rm(past, current, p)
  gc()
}
png_files <- sort(list.files(frame_dir, pattern="frame_.*\\.png", full.names=T))
av_encode_video(png_files, output="biocrust_timeseries.mp4", framerate=6)
cat("done\n")
#anim_save("biocrust_animation.gif", animation=anim_rendered)



## what percent of the footprint is biocrust? and which types of crust?
transects %>%
  filter(has_cyano | has_lich | has_moss) %>%
  summarise(
    total_crust  = n(),
    pct_cyano    = sum(has_cyano)  / n() * 100,
    pct_lich   = sum(has_lich) / n() * 100,
    pct_moss     = sum(has_moss)   / n() * 100
  )

# then apply those proportions to the mean weighted biocrust cover
results %>%
  summarise(mean_w_biocrust = mean(w_biocrust, na.rm = TRUE)) %>%
  mutate(
    pct_cyano  = mean_w_biocrust * (43/76),
    pct_lich = mean_w_biocrust * (19/76),
    pct_moss   = mean_w_biocrust * (14/76)
  )

## let's plot the variation of these percentages
library(ggplot2)
library(tidyr)
library(dplyr)

# results long form for easy plotting
results_long <- results %>% 
  filter(!is.na(w_cyano)) %>% 
  select(year, month, period, w_cyano, w_lich, w_moss) %>% 
  pivot_longer(
    cols = c(w_cyano, w_lich, w_moss),
    names_to ="crust_type",
    values_to = "weighted_cover"
  ) %>% 
  mutate(
    crust_type = case_when(
      crust_type == "w_cyano" ~ "cyanobacteria",
      crust_type == "w_lich" ~ "lichen",
      crust_type == "w_moss" ~ "moss"
    ),
    month_name = factor(month.abb[month], levels=month.abb)
  )

# plot monthly variation by crust type and day/night
ggplot(results_long, aes(x=month_name, y=weighted_cover, color=crust_type, group = crust_type)) +
  geom_line(aes(group = interaction(crust_type, year)),
            alpha=0.2, linewidth=0.4) +
  # mean across years
  stat_summary(fun=mean, geom="line", linewidth=1) +
  stat_summary(fun=mean, geom="point", size=2) +
  facet_wrap(~period, ncol=1) +
  scale_color_manual(
    values=c("cyanobacteria" = "#8b7355",
             "lichen" = "#4a7c59",
             "moss" = "darkgreen"),
    name="crust type"
  ) +
  labs(title="footprint-weighted bsc cover by month for US-Mpj",
       x="month", y="weighted cover fraction") +
  theme_classic()

# Is this significant? ----------------------------------------------------
## MANOVA: does the overall vegetation weighting vary significantly by year, month, period
## tests all 5 cover fractions simulatenously as a multivariate response
results_model <- results %>% 
  mutate(year=factor(year),
         month = factor(month),
         period = factor(period))

# fit multivariate model with all 5 cover fractions as a response
manova_fit <- manova(
  cbind(w_cyano, w_lich, w_moss, w_midveg, w_subcanopy, w_canopy, w_abiotic) ~ year + month + period,
  data=results_model
)

## Pillai's trace: test statistic
summary(manova_fit, test="Pillai")

# follow up with univariate ANOVAs for each cover type
# to see which ones are driving any significant effects
summary.aov(manova_fit)

## Year: small but significant; substaintailly smaller than seasonal variation
## yearly significance suggests 2014-2018 climatologies are slightly different but 
## generally broadly representatives (if we observe SumSq for year vs month)

## Month: strong predictor for every biotic code
## Period: largest effect sizes

## abiotic: year not significant, month only marginally. this stuff doesn't change
## since it's not biological material

## if we just want to look at BSC summed
results_model <- results_model %>% 
  mutate(w_bsc_total = w_cyano + w_lich + w_moss)

aov_bsc <- aov(w_bsc_total ~ year + month + period, data=results_model)
summary(aov_bsc)

library(kableExtra)

tidy(aov_bsc) %>%
  rename(Term = term, Df = df, `Sum Sq` = sumsq, `Mean Sq` = meansq, 
         `F value` = statistic, `P value` = p.value) %>%
  mutate(across(c(`Sum Sq`, `Mean Sq`, `F value`, `P value`), 
                ~ifelse(Term == "year", 
                        formatC(., format="f", digits=3),
                        formatC(., format="e", digits=3)))) %>%
  kable(caption="ANOVA: Biocrust cover fractions by year, month, and period",
        format="html") %>%
  kable_styling(bootstrap_options=c("striped","hover"),
                full_width=FALSE) %>%
  row_spec(which(tidy(aov_bsc)$p.value < 0.05), bold=TRUE) %>%
  save_kable("anova_biocrust.png", zoom=5, vheight=300)
# Biocrust weighting map --------------------------------------------------
library(patchwork)

# 1. compute mean footprint per month andperiod
## average all years for each month/period combination
## months: March (3), June (6), August (8), November (11)
target_months <- c(3, 6, 8, 11)
month_names <- c("March", "June", "August", "November")
target_periods <- c("DAY", "NIGHT")

# function to compute mean footprint raster for one month/period
mean_footprint <- function(month_num, period_str) {
  # get fps for this month and period across all years
  files <- results %>% 
    filter(month == month_num, period == period_str) %>% 
    pull(filepath)
  
  # load all rasters individually (b/c indvudal yearly footprint
  # rasters for the same month don't have identical extents )
  rast_list <- lapply(files, rast)
  
  # use first raster as reference grid
  # resample all others to match it
  ref <- rast_list[[1]]
  
  resampled <- lapply(rast_list, function(r) {
    if (!all(ext(r) == ext(ref)) | !all(res(r) == res(ref))) {
      resample(r, ref, method="bilinear")
    } else {
      r
    }
  })
  
  # stack and average
  stack <- rast(resampled)
  mean_fp <- mean(stack, na.rm=TRUE)
  
  # normalize so weights sum to 1
  mean_fp <- mean_fp / global(mean_fp, "sum", na.rm=TRUE)[[1]]
  
  return(mean_fp)
}
# 2. rerun fp_list
fp_list <- list()
for (m in target_months) {
  for (p in target_periods) {
    key <- paste0("m", m, "_", p)
    fp_list[[key]] <- mean_footprint(m, p)
  }
}
cat("Computed", length(fp_list), "mean footprints\n")

# 3. build vegetation raster for each footprint grid
## each mean footprint may have a slightly different extent so we can build a fresh
## veg raster for each one
## NOTE: this calls build_veg_raster() (singular "veg"), which is not defined
## anywhere in the original script — it's presumably an earlier/renamed version
## of build_veg_idw_raster() from Section 5. Left as-is per the original file;
## you'll likely want to point this at build_veg_idw_raster() instead.
veg_list <- list()
for (key in names(fp_list)) {
  veg_list[[key]] <- build_veg_raster(fp_list[[key]])
}
cat("built", length(veg_list), "vegetation rasters\n")


# Plotting function for 1 panel -------------------------------------------
# shows biocrust cover as a background color with footprint contour lines overlaid
make_panel <- function(month_num, period_str, month_name) {
  key <- paste0("m", month_num, "_", period_str)
  fpt <- fp_list[[key]]
  veg <- veg_list[[key]]
  
  # extract biocrust layer
  biocrust <- veg[["cover_biocrust"]]
  
  # convert to data frames for ggplot
  bio_df <- as.data.frame(biocrust, xy=TRUE)
  colnames(bio_df) <- c("lon", "lat", "biocrust")
  fpt_df <- as.data.frame(fpt, xy=TRUE)
  colnames(fpt_df) <- c("lon", "lat", "weight")
  
  # assign weight cols directly by position; both rasters have identical extents and resolutions
  # so rows are always in the same order
  plot_df <- bio_df
  plot_df$weight <- fpt_df$weight
  
  # compute contour levels for 50% and 80% of footprint
  # weight thresholds where cumulative sum from peak reaches 50% and 80% of total
  fpt_vals <- sort(values(fpt), decreasing=TRUE)
  fpt_vals <- fpt_vals[!is.na(fpt_vals)]
  fpt_cum <- cumsum(fpt_vals) / sum(fpt_vals) # normalize to 0-1
  
  # find weight threshold for 50% and 80% contours
  level_50 <- fpt_vals[which.min(abs(fpt_cum - 0.50))]
  level_80 <- fpt_vals[which.min(abs(fpt_cum - 0.80))]
  
  cat(key, "-level_50:", level_50, "level_80:", level_80, "\n")
  # plot
  ggplot(plot_df, aes(x=lon, y=lat)) +
    # biocrust cover as background fill
    geom_raster(aes(fill=biocrust)) +
    # footprint contour lines drawn from same df
    geom_contour(aes(z=weight),
                 breaks=c(level_50, level_80),
                 color="black", linewidth=0.6) +
    # tower location
    geom_point(data=data.frame(lon=tower_lon, lat=tower_lat),
               aes(x=lon, y=lat), color="black", size=2, shape=17) +
    # biocrust cover scale
    scale_fill_gradientn(colors=c("#f5f0e8", "#c8b89a", "#8b7355", "#4a7c59", "#2d5a1b"),
                         limits=c(0, 0.35),
                         name="Biocrust\ncover",
                         na.value="white") +
    labs(title=paste(month_name, "-", period_str),
         x=NULL, y=NULL) +
    theme_classic() +
    theme(plot.title=element_text(size=9, face="bold", hjust=0.5),
          legend.position = "right",
          axis.text=element_text(size=7),
          panel.grid=element_line(color="white", linewidth = 0.3)) +
    coord_equal()
}

# Build all 8 panels and combine ------------------------------------------
# layout: 4 rows (months) x 2 columns (day, night)
panels <- list()
for (i in seq_along(target_months)) {
  m <- target_months[i]
  name <- month_names[i]
  panels[[paste0(name, "_DAY")]] <- make_panel(m, "DAY", name)
  panels[[paste0(name, "_NIGHT")]] <- make_panel(m, "NIGHT", name)
}

# combine into 4 x 2 grid using patchwork
combined <- (panels[["March_DAY"]]    | panels[["March_NIGHT"]])   /
  (panels[["June_DAY"]]     | panels[["June_NIGHT"]])    /
  (panels[["August_DAY"]]   | panels[["August_NIGHT"]])  /
  (panels[["November_DAY"]] | panels[["November_NIGHT"]]) +
  plot_annotation(
    title    = "Footprint-weighted biocrust cover — US-Mpj",
    subtitle = "black contours show 50% and 80% footprint extent. Triangle = flux tower.",
    theme    = theme(
      plot.title    = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 9,  color = "gray40")
    )
  )

print(combined)


# check actual extent of March DAY mean footprint
fpt <- fp_list[["m3_DAY"]]
fpt_df <- as.data.frame(fpt, xy = TRUE)
colnames(fpt_df) <- c("lon", "lat", "weight")

# convert extent to meters from tower
fpt_df <- fpt_df %>%
  mutate(
    x_m = (lon - tower_lon) * m_per_deg_lon,
    y_m = (lat - tower_lat) * m_per_deg_lat,
    dist_m = sqrt(x_m^2 + y_m^2)
  )

cat("Max distance from tower in raster:", max(fpt_df$dist_m), "m\n")
cat("Raster extent in meters:\n")
cat("  x:", range(fpt_df$x_m), "\n")
cat("  y:", range(fpt_df$y_m), "\n")

# what fraction of total weight falls beyond 100m?
total_weight <- sum(fpt_df$weight, na.rm = TRUE)
beyond_100m  <- sum(fpt_df$weight[fpt_df$dist_m > 100], na.rm = TRUE)
cat("Fraction of footprint weight beyond 100m:", 
    round(beyond_100m / total_weight * 100, 1), "%\n")
