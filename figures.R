## Dryland Biocrust (BSC) Carbon Flux Model — Figure Scripts
## Asha Kannan
## Requires model output saved by DrylandBSCModel_June2026.R


# Setup -------------------------------------------------------------------
library(maptiles)
library(tidyterra)
library(terra)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(sf)
library(tidyverse)
library(lubridate)
library(patchwork)
library(viridis)
library(cowplot)

setwd("/Users/ashakannan/Documents/GitHub/DrylandBSCModel")

model_output  <- readRDS("output/model_output.rds")
model_daily   <- readRDS("output/model_daily.rds")
model_monthly <- readRDS("output/model_monthly.rds")
model_annual  <- readRDS("output/model_annual.rds")
model_input   <- readRDS("output/model_input.rds")


# Source model functions and parameters needed for ecophys curves
source("DrylandBSCModel_June2026.R")  # runs the model and generates/saves outputs
                     # (comment out saveRDS lines if re-sourcing to avoid overwrite)

# Color palette
bsc_colors <- c(
  "Moss crust"          = "#228833",
  "Cyanobacteria crust" = "goldenrod3",
  "Lichen crust"        = "#4477AA"
)

# Seasonal shading helpers
shade_monsoon <- data.frame(
  xmin = as.Date("2016-06-15"), xmax = as.Date("2016-09-30"),
  ymin = -Inf, ymax = Inf, season = "Monsoon season")
shade_shoulder <- data.frame(
  xmin = c(as.Date("2016-03-01"), as.Date("2016-10-01")),
  xmax = c(as.Date("2016-05-31"), as.Date("2016-11-30")),
  ymin = -Inf, ymax = Inf, season = "Shoulder season")


# Figure 1: Map of Site ---------------------------------------------------
#install.packages("maptiles")
#install.packages("tidyterra")
#install.packages("ggspatial")

# Site constants
tower_lat <- 34.4385
tower_lon <- -106.2377

# meters per degree at this latitude
m_per_deg_lat <- 111320
m_per_deg_lon <- 111320 * cos(tower_lat * pi / 180)

# function to convert offset in meters from tower to lon/lat
offset_lonlat <- function(dx_m, dy_m) {
  list(
    lon = tower_lon + dx_m / m_per_deg_lon,
    lat = tower_lat + dy_m / m_per_deg_lat
  )
}

# tower point
tower_sf <- st_as_sf(
  data.frame(lon=tower_lon, lat=tower_lat, label= "Flux Tower"),
  coords = c("lon", "lat"), crs=4326
)

# tower label
tower_legend_df <- data.frame(
  lon   = tower_lon,
  lat   = tower_lat,
  label = "Flux Tower"
)
tower_legend_sf <- st_as_sf(tower_legend_df, coords = c("lon", "lat"), crs = 4326)


# transect lines defined from radial interpolation code
transect_sf <- st_as_sf(
  transect_lines,
  coords = c("x_start", "y_start"), crs=4326
)

# rebuild as linestrings since st_as_sf doesn't handle line endpoints directly
transect_lines_sf <- transect_lines %>% 
  rowwise() %>% 
  mutate(geometry = st_sfc(
    st_linestring(matrix(
      c(x_start, y_start, x_end, y_end), ncol=2, byrow=TRUE 
    )),
    crs=4326
  )) %>% 
  ungroup() %>% 
  st_as_sf()

# Soil Water Sensors -- approximated from reference map
# Open Sensors: rougly N and NW of tower (20-30m)
# Juniper sensors: roughly W and SW (25-35m)
# Pinyon pine sensors: roughly E (20-30m)
# Adjust dx_m/dy_m manually once you print the figure!
sensor_coords <- bind_rows(
  # Open sensors
  data.frame(type = "Open",    dx =  -5, dy =  -8),
  data.frame(type = "Open",    dx =  -3, dy =  -3),
  data.frame(type = "Open",    dx =   3, dy =  -6),
  # Juniper sensors
  data.frame(type = "Juniper", dx = -15, dy =   3),
  data.frame(type = "Juniper", dx = -10, dy =  -8),
  data.frame(type = "Juniper", dx =  -8, dy =   8),
  # Pinyon sensors
  # Pinyon sensors
  data.frame(type = "Pinyon", dx = -11, dy = -11),  # SW, 15m (P13 area)
  data.frame(type = "Pinyon", dx = -16, dy =   8),  # WNW, ~18m (P12 area)
  data.frame(type = "Pinyon", dx = -18, dy =   4),  # W, ~18m (P2 area)
  # Flux tower
  data.frame(type = "Flux Tower", dx =   0, dy =   0)
) %>%
  rowwise() %>% 
  mutate(
    lon = tower_lon + dx / m_per_deg_lon,
    lat = tower_lat + dy / m_per_deg_lat
  ) %>% 
  ungroup()
sensors_sf <- st_as_sf(sensor_coords, coords=c("lon", "lat"), crs=4326)

# Soil pits (approx from reference map)
pit_coords <- bind_rows(
  data.frame(type = "Open",    dx =  -8, dy =  18),
  data.frame(type = "Juniper", dx = -28, dy =   5),
  data.frame(type = "Pinyon",   dx =  22, dy = -20)
) %>%
  rowwise() %>%
  mutate(
    lon = tower_lon + dx / m_per_deg_lon,
    lat = tower_lat + dy / m_per_deg_lat
  ) %>%
  ungroup()
pits_sf <- st_as_sf(pit_coords, coords=c("lon", "lat"), crs=4326)

# footprint contours from mean annual day footprint (assumes mean_annual_day alr computed)
fpt_df <- as.data.frame(mean_annual_day, xy=T)
colnames(fpt_df) <- c("lon", "lat", "weight")
fpt_vals <- sort(fpt_df$weight, decreasing=TRUE, na.last=NA)
fpt_cum <- cumsum(fpt_vals) / sum(fpt_vals)
level_50 <- fpt_vals[which.min(abs(fpt_cum - 0.50))]
level_80 <- fpt_vals[which.min(abs(fpt_cum - 0.80))]

# build biocrust cover df from veg_f2
biocrust_rast <- max(veg_f2[["cover_cyano"]],
                     veg_f2[["cover_lich"]],
                     veg_f2[["cover_moss"]],
                     na.rm = TRUE)
names(biocrust_rast) <- "biocrust"

bsc_df <- as.data.frame(biocrust_rast, xy = TRUE)
colnames(bsc_df) <- c("lon", "lat", "biocrust")

# bounding box for map tile download & zoom in around tower with 150 m buffer on all sides
buf_deg_lon <- 150 / m_per_deg_lon
buf_deg_lat <- 150 / m_per_deg_lat

bbox <- st_bbox(c(
  xmin = tower_lon - buf_deg_lon,
  xmax = tower_lon + buf_deg_lon,
  ymin = tower_lat - buf_deg_lat,
  ymax = tower_lat + buf_deg_lat
), crs=st_crs(4326))

# download satellite tiles
tile_cache_dir <- file.path(path_project_dir, "data", "basemap_tiles")
dir.create(tile_cache_dir, showWarnings = FALSE, recursive=TRUE)

tiles <- get_tiles(bbox, provider = "Esri.WorldImagery", zoom=18, crop=TRUE,
                   cachedir = tile_cache_dir, verbose = TRUE)



# Build the figure!
library(ggnewscale)
# split sensors from tower
sensors_only <- sensors_sf %>% filter(type != "Flux Tower")
tower_only   <- sensors_sf %>% filter(type == "Flux Tower")

fig_map <- ggplot() +
  # satellite basemap
  geom_spatraster_rgb(data=tiles) +
  
  # biocrust cover overlay (alpha to show veg structure
  geom_raster(data=bsc_df %>% filter(!is.na(biocrust)),
              aes(x=lon, y=lat, fill=biocrust), alpha=0.6) +
  scale_fill_gradient(low="transparent", high="mediumpurple4",
                      name="Biocrust\ncover fraction",
                      na.value="transparent") +
  new_scale_fill() +
  
  # footprint contours -- 50% and 80% as separate layers so each gets a legend label
  geom_contour(data=fpt_df %>% filter(!is.na(weight)),
               aes(x=lon, y=lat, z=weight, linetype="50% footprint"),
               breaks=level_50,
               color="black", linewidth=1.2, alpha=0.6) +
  geom_contour(data=fpt_df %>% filter(!is.na(weight)),
               aes(x=lon, y=lat, z=weight, linetype="80% footprint"),
               breaks=level_80,
               color="black", linewidth=1.2, alpha=0.6) +
  
  # transect lines
  geom_sf(data=transect_lines_sf,
          aes(linetype="Transect"),
          color="black", linewidth=1.2) +
  
  # combined linetype legend for contours + transects
  scale_linetype_manual(
    name="Footprint / Transects",
    values=c("50% footprint"="solid",
             "80% footprint"="dashed",
             "Transect"="dotted"),
    guide=guide_legend(override.aes=list(linewidth=1.6, color="black"),
                       keywidth=unit(1.2, "cm"))
  ) +
  
  # flux tower and sensors
  geom_sf(data=sensors_only,
          aes(shape=type, fill=type), color="black",
          size=3.2, stroke=0.6) +
  geom_sf(data=tower_only,
          aes(shape=type, fill=type), color="black",
          size=4.8, stroke=1.2) +
  scale_fill_manual(
    values=c(
      "Open" = "yellow",
      "Juniper" = "maroon3",
      "Pinyon" = "darkblue",
      "Flux Tower" = "white"
    ),
    name="Site sensors"
  ) +
  scale_shape_manual(values=c(
    "Open" = 22,      # filled square (was 0)
    "Juniper" = 21,   # filled circle (was 1)
    "Pinyon" = 25,    # filled triangle (was 2)
    "Flux Tower" = 24 # keep as filled triangle w/ outline — distinct from Pinyon
  ),
  name="Site sensors") +
  # north arrow and scale bar via annotation
  annotation_scale(location = "br", width_hint = 0.3,
                   text_col = "white", line_col = "white",
                   bar_cols = c("white", "grey50"),
                   text_cex=1) +
  annotation_north_arrow(location = "tl",
                         style = north_arrow_minimal(text_col = "white",
                                                     line_col = "white",
                                                     fill = "white",
                                                     text_size=8)) +
  
  coord_sf(
    xlim = c(tower_lon - buf_deg_lon, tower_lon + buf_deg_lon),
    ylim = c(tower_lat - buf_deg_lat, tower_lat + buf_deg_lat),
    expand = FALSE
  ) +
  
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(
    axis.text = element_text(size=8),
    axis.title = element_text(size=12),
    legend.text = element_text(size=12),
    legend.title = element_text(size=12),
    legend.background = element_rect(fill = "white", color=NA),
    legend.position = "right",
    legend.key.height=unit(0.8, "cm")
  )

ggsave("KannanFig1.tif", fig_map, width = 8, height = 8, units="in", dpi = 600, compression="lzw")

# Figure 3 — Ecophysiological response curves + site environment ----------
cy_Pg_T <- function(t) pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=t, T_k=T_k_cyanolich, Pgmax=Pgmax_cyanolich) *
                              f_H(Tv=t, delta_S=delta_S_cyanolich, delta_Hd=delta_Hd_cyanolich))
ch_Pg_T <- function(t) pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=t, T_k=T_k_chlorolich, Pgmax=Pgmax_chlorolich) *
                              f_H(Tv=t, delta_S=delta_S_chlorolich, delta_Hd=delta_Hd_chlorolich))


temp_seq <- seq(0, 45, 1)
temp_curve_all <- tibble(temp = temp_seq) %>%
  rowwise() %>%
  mutate(
    Moss  = pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=temp, T_k=T_k_moss,  Pgmax=Pgmax_moss) *
                   f_H(Tv=temp, delta_S=delta_S_moss,  delta_Hd=delta_Hd_moss)),
    Cyano = pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=temp, T_k=T_k_cyano, Pgmax=Pgmax_cyano) *
                   f_H(Tv=temp, delta_S=delta_S_cyano, delta_Hd=delta_Hd_cyano)),
    Cyanolichen  = cy_Pg_T(temp),
    Chlorolichen = ch_Pg_T(temp),
    Lich = w_cyanolich * Cyanolichen + w_chlorolich * Chlorolichen
  ) %>%
  ungroup() %>%
  pivot_longer(-temp, names_to="type", values_to="Pg_T") %>%
  mutate(type = recode(type,
                       Moss  = "Moss crust",
                       Cyano = "Cyanobacteria crust",
                       Lich  = "Lichen crust"))

rd_curve_all <- tibble(temp = temp_seq) %>%
  rowwise() %>%
  mutate(
    Moss  = calc_RD_fxn_T_rowwise(temp=temp, r_r=r_r_moss,  R_o=RD_min_moss),
    Cyano = calc_RD_fxn_T_rowwise(temp=temp, r_r=r_r_cyano, R_o=RD_min_cyano),
    Cyanolichen  = calc_RD_fxn_T_rowwise(temp=temp, r_r=r_r_cyanolich,  R_o=R0_cyanolich),
    Chlorolichen = calc_RD_fxn_T_rowwise(temp=temp, r_r=r_r_chlorolich, R_o=R0_chlorolich),
    Lich = w_cyanolich * Cyanolichen + w_chlorolich * Chlorolichen
  ) %>%
  ungroup() %>%
  pivot_longer(-temp, names_to="type", values_to="RD_T") %>%
  mutate(type = recode(type,
                       Moss  = "Moss crust",
                       Cyano = "Cyanobacteria crust",
                       Lich  = "Lichen crust"))

par_seq <- seq(0, 3000, 1)
light_curve_all <- tibble(PAR = par_seq) %>%
  mutate(
    `Moss crust`          = calc_PN_fxn_I_rowwise(PAR=PAR, f_Io_Icomp=f_Io_Icomp_moss,
                                                  b=b_moss,  g=g_moss,  Icomp=Icomp_moss)  * moss_light_scale,
    `Cyanobacteria crust` = calc_PN_fxn_I_rowwise(PAR=PAR, f_Io_Icomp=f_Io_Icomp_cyano,
                                                  b=b_cyano, g=g_cyano, Icomp=Icomp_cyano) * cyano_light_scale,
    Cyanolichen  = calc_PN_fxn_I_rowwise(PAR=PAR, f_Io_Icomp=f_cyanolich, b=b_cyanolich, g=g_cyanolich, Icomp=Icomp_cyanolich) * cyanolich_light_scale,
    Chlorolichen = calc_PN_fxn_I_rowwise(PAR=PAR, f_Io_Icomp=f_chlorolich, b=b_chlorolich, g=g_chlorolich, Icomp=Icomp_chlorolich) * chlorolich_light_scale,
    `Lichen crust` = w_cyanolich * Cyanolichen + w_chlorolich * Chlorolichen
  ) %>%
  pivot_longer(-PAR, names_to="type", values_to="Pg_I")

wc_seq <- seq(0, 1, 0.01)
water_curve_all <- tibble(WC = wc_seq) %>%
  mutate(
    `Moss crust` = ifelse(WC < 0.05 / max(model_output$SWC_moss, na.rm=TRUE), 0,
                          (low_water_response_fxn(WC=WC, WC_k=WC_k_moss, Pgmax=1) *
                             high_water_response_fxn(WC=WC, slope=slope_moss,
                                                     intercept=intercept_moss, Pgmax=1)) / WC_scale_moss),
    `Cyanobacteria crust` = {
      raw  <- wc_cyano_coef[1] + wc_cyano_coef[2]*WC + wc_cyano_coef[3]*WC^2
      peak <- max(wc_cyano_coef[1] + wc_cyano_coef[2]*seq(0,1,0.001) +
                    wc_cyano_coef[3]*seq(0,1,0.001)^2)
      threshold <- 0.05 / max(model_output$SWC_cyano, na.rm=TRUE)
      ifelse(WC < threshold, 0, pmax(0, raw / peak))
    },
    `Lichen crust` = {
      threshold <- 0.05 / max(model_output$SWC_lich, na.rm=TRUE)
      ifelse(WC < threshold, 0,
             w_psora * calc_WC_frac_psora(WC) + w_collema * calc_WC_frac_collema(WC))
    }
  ) %>%
  pivot_longer(-WC, names_to="type", values_to="WC_frac")

keep <- c("Moss crust", "Cyanobacteria crust", "Lichen crust")

temp_curve_all  <- temp_curve_all  %>% filter(type %in% keep)
rd_curve_all    <- rd_curve_all    %>% filter(type %in% keep)
light_curve_all <- light_curve_all %>% filter(type %in% keep)
# day of year climatology
env_doy <- model_input %>%
  filter(PPFD_IN > 0) %>%
  mutate(doy = yday(as.Date(datetime))) %>%
  group_by(doy) %>%
  summarise(
    mean_T   = mean(SOILT_O_5, na.rm=TRUE),
    min_T    = min(SOILT_O_5,  na.rm=TRUE),
    max_T    = max(SOILT_O_5,  na.rm=TRUE),
    mean_PAR = mean(PPFD_IN,   na.rm=TRUE),
    min_PAR  = min(PPFD_IN,    na.rm=TRUE),   # actual instantaneous min
    max_PAR  = max(PPFD_IN,    na.rm=TRUE),   # actual instantaneous max
    mean_SWC = mean(SWC_O_5,   na.rm=TRUE),
    min_SWC  = min(SWC_O_5,    na.rm=TRUE),
    max_SWC  = max(SWC_O_5,    na.rm=TRUE),
    .groups = "drop"
  )

p_temp <- ggplot(temp_curve_all, aes(x=temp, y=Pg_T, color=type, linetype="Carbon assimilation")) +
  geom_line(data=rd_curve_all, aes(x=temp, y=RD_T, color=type, linetype="Dark respiration"),
            linewidth=1.5) +
  geom_line(linewidth=2) +
  scale_color_manual(values=bsc_colors) +
  scale_linetype_manual(
    values=c("Carbon assimilation" = "solid", "Dark respiration" = "dashed"),
    breaks=c("Carbon assimilation", "Dark respiration"),
    name=NULL)+
  guides(color = "none",linetype=guide_legend(override.aes = list(color="black", linewidth=1.2))) +
  labs(x="Temperature (°C)", y=expression(CO[2]~exchange~(µmol~CO[2]~m^-2~s^-1)),
       color=NULL, tag="a)") +
  theme_classic() +
  theme(legend.position=c(0.75,1.1),
        legend.key.width=unit(1.2, "cm"),
        legend.direction="horizontal",
        legend.text = element_text(size=16),
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=25),
        axis.text    = element_text(size=18),
        plot.tag     = element_text(size=14, face="bold"))

p_light <- ggplot(light_curve_all, aes(x=PAR, y=Pg_I, color=type)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=bsc_colors) +
  labs(x=expression("PPFD"~("µmol photons"~m^-2~s^-1)),
       y=expression(P[net]~(µmol~CO[2]~m^-2~s^-1)),
       color=NULL, tag="g)") +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text    = element_text(size=18),
        plot.tag     = element_text(size=14, face="bold"))

p_water <- ggplot(water_curve_all, aes(x=WC, y=WC_frac, color=type)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=bsc_colors) +
  labs(x="Relative Water Content", y="Proportion of Pgmax",
       color=NULL, tag="d)") +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size=25),
        axis.title.y = element_text(size=30),
        axis.text    = element_text(size=18),
        plot.tag     = element_text(size=14, face="bold"))

p4 <- ggplot(env_doy, aes(x=doy)) +
  geom_ribbon(aes(ymin=min_T, ymax=max_T), fill="#F44336", alpha=0.2) +
  geom_line(aes(y=mean_T), color="#F44336", linewidth=2) +
  scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                     labels=month.abb) +
  labs(x="Month", y="Mean Soil Temperature (°C)", tag="b)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=25),
        axis.title = element_text(size=20),
        axis.text  = element_text(size=18),
        plot.tag   = element_text(size=14, face="bold"))

p5 <- ggplot(env_doy, aes(x=doy)) +
  geom_ribbon(aes(ymin=min_PAR, ymax=max_PAR), fill="#FFC107", alpha=0.2) +
  geom_line(aes(y=mean_PAR), color="#FFC107", linewidth=2) +
  scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                     labels=month.abb) +
  labs(x="Month", y=expression("Mean Daytime PPFD"~("µmol photons"~m^-2~s^-1)), tag="h)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=25),
        axis.title = element_text(size=20),
        axis.text  = element_text(size=18),
        plot.tag   = element_text(size=14, face="bold"))

p6 <- ggplot(env_doy, aes(x=doy)) +
  geom_ribbon(aes(ymin=min_SWC, ymax=max_SWC), fill="#2196F3", alpha=0.2) +
  geom_line(aes(y=mean_SWC), color="#2196F3", linewidth=2) +
  scale_x_continuous(breaks=c(1,32,60,91,121,152,182,213,244,274,305,335),
                     labels=month.abb) +
  labs(x="Month", y="Mean Soil Water Content (m³/m³)", tag="e)") +
  theme_classic() +
  theme(axis.title.x = element_text(size=25),
        axis.title = element_text(size=20),
        axis.text  = element_text(size=18),
        plot.tag   = element_text(size=14, face="bold"))

legend_plot <- ggplot(temp_curve_all, aes(x=temp, y=Pg_T, color=type)) +
  geom_line(linewidth=10) +
  scale_color_manual(values=bsc_colors, name=NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal", 
        legend.text = element_text(size=50),
        legend.justification = "left",
        legend.box.just="left")
shared_legend <- cowplot::get_legend(legend_plot)

main_grid <- (p_temp | p4) / (p_water | p6) / (p_light | p5) +
  plot_annotation(theme=theme(plot.title=element_text(size=16, face="bold"))) &
  theme(plot.tag=element_text(size=30, face="bold"),
        plot.tag.position=c(0.18, 0.95))

fig3 <- cowplot::plot_grid(shared_legend, main_grid, ncol=1, rel_heights=c(0.15, 1))
fig3
setwd("/Users/ashakannan/Library/CloudStorage/OneDrive-Personal/Desktop/Kannan_et_al(2026)")
ggsave("KannanFig3_edited.pdf", fig3, width=20, height=20, dpi=600, bg="white")


# Figure 4 — BSC carbon flux timeseries -----------------------------------
ecosystem_doy <- model_daily %>%
  mutate(doy = yday(date)) %>%
  filter(GPP_daily > -1) %>%
  group_by(doy) %>%
  summarise(GPP_mean = mean(GPP_daily, na.rm=TRUE), 
            GPP_sd = sd(GPP_daily, na.rm=TRUE),
            .groups="drop")

bsc_doy <- model_daily %>%
  mutate(doy = yday(date)) %>%
  group_by(doy) %>%
  summarise(
    GPP_mean = mean(Pg_BSC_daily,   na.rm=TRUE),
    GPP_sd = sd(Pg_BSC_daily, na.rm=TRUE),
    R_mean   = mean(RD_BSC_daily,   na.rm=TRUE),
    R_sd = sd(RD_BSC_daily, na.rm=TRUE),
    NEE_mean = mean(Pnet_BSC_daily, na.rm=TRUE),
    NEE_sd = sd(Pnet_BSC_daily, na.rm=TRUE),
    .groups="drop"
  )

bsc_org_doy <- model_daily %>%
  mutate(doy = yday(date)) %>%
  group_by(doy) %>%
  summarise(
    Moss_mean = mean(Pg_moss_daily, na.rm = TRUE),
    Moss_sd = sd(Pg_moss_daily, na.rm = TRUE),
    Cyano_mean = mean(Pg_cyano_daily, na.rm=TRUE),
    Cyano_sd = sd(Pg_cyano_daily, na.rm=TRUE),
    Lich_mean  = mean(Pg_lich_daily,  na.rm=TRUE),
    Lich_sd = sd(Pg_lich_daily, na.rm=TRUE),
    .groups="drop"
  )

doy_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)

p_a <- ggplot(ecosystem_doy, aes(x=doy, y=GPP_mean)) +
  geom_ribbon(aes(ymin=GPP_mean - GPP_sd, ymax=GPP_mean + GPP_sd),
              fill="burlywood4", alpha=0.2) +
  geom_line(aes(color="Ecosystem GPP"), linewidth=1.2) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_continuous(breaks=doy_breaks, labels=month.abb) +
  scale_color_manual(values=c("Ecosystem GPP"="burlywood4")) +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  labs(x=NULL, y=expression("GPP"~(gC~m^-2~d^-1)), tag="a)") +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        legend.text = element_text(size=30),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        plot.tag=element_text(size=14, face="bold"),
        plot.margin=margin(4, 4, 0, 4))

p_b <- ggplot(bsc_doy, aes(x=doy)) +
  geom_ribbon(aes(ymin=GPP_mean - GPP_sd, ymax=GPP_mean + GPP_sd), fill="#069494", alpha=0.2) +
  geom_ribbon(aes(ymin=-(R_mean - R_sd), ymax=-(R_mean + R_sd)), fill="tomato3", alpha=0.2) +
  geom_ribbon(aes(ymin=NEE_mean - NEE_sd, ymax=NEE_mean + NEE_sd), fill="black", alpha=0.2) +
  geom_line(aes(y=GPP_mean,  color="Biocrust GPP"), linewidth=1.0, alpha=0.8) +
  geom_line(aes(y=-R_mean,   color="Biocrust Respiration"),   linewidth=1.0, alpha=0.8) +
  geom_line(aes(y=NEE_mean,  color="Biocrust Net C Exchange"), linewidth=1.2) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_continuous(breaks=doy_breaks, labels=month.abb) +
  scale_color_manual(values=c("Biocrust GPP"="#069494", "Biocrust Respiration"="tomato3", "Biocrust Net C Exchange"="black")) +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  labs(x=NULL, y=expression("Carbon exchange"~(gC~m^-2~d^-1)), color=NULL, tag="b)") +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="right",
        legend.text = element_text(size=30),
        axis.text = element_text(size=20),
        axis.title = element_text(size=20),
        plot.tag=element_text(size=14, face="bold"),
        plot.margin=margin(0, 4, 0, 4))

p_c <- ggplot(bsc_org_doy, aes(x=doy)) +
  geom_ribbon(aes(ymin=Moss_mean - Moss_sd, ymax=Moss_mean + Moss_sd),
              fill=bsc_colors["Moss crust"], alpha = 0.2) +
  geom_ribbon(aes(ymin=Cyano_mean - Cyano_sd, ymax=Cyano_mean + Cyano_sd),
              fill=bsc_colors["Cyanobacteria crust"], alpha=0.2) +
  geom_ribbon(aes(ymin=Lich_mean - Lich_sd, ymax=Lich_mean + Lich_sd),
              fill=bsc_colors["Lichen crust"], alpha=0.2) +
  geom_line(aes(y=Moss_mean, color="Moss crust"), linewidth=1.0) +
  geom_line(aes(y=Cyano_mean, color="Cyanobacteria crust"), linewidth=1.0) +
  geom_line(aes(y=Lich_mean, color="Lichen crust"), linewidth=1.0) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_continuous(breaks=doy_breaks, labels=month.abb) +
  scale_color_manual(values=bsc_colors) +
  labs(x="Day of year", y=expression("Biocrust GPP"~(gC~m^-2~d^-1)), color=NULL, tag="c)") +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  theme_classic(base_size=13) +
  theme(legend.position="right",
        axis.title=element_text(size=20),
        axis.text = element_text(size=20),
        legend.text = element_text(size=30),
        plot.tag=element_text(size=14, face="bold"),
        plot.margin=margin(0, 4, 4, 4))

fig4 <- p_a / p_b / p_c + plot_layout(heights=c(1, 1, 1)) &
  theme(plot.tag=element_text(size=30, face="bold"),
        plot.tag.position=c(0.12, 0.95))
fig4
ggsave("KannanFig4.pdf", fig4, width=15, height=20, dpi=600)


# Figure 5 — BSC percent of ecosystem GPP --------------------------------
# Interannual sensitivity ribbon
full_years <- c(2010, 2011, 2013, 2015, 2017, 2018, 2019, 2020, 2021)
gpp_threshold <- 0.5 # unified denominator-stability filter for both panels
min_years <- 4 # min contributiing years required to plot panel b's mean and ribbon
doy_breaks <- c(1,32,60,91,121,152,182,213,244,274,305,335)

pct_doy <- model_output %>%
  mutate(date      = as.Date(datetime),
         year      = lubridate::year(date),
         month_day = format(date, "%m-%d")) %>%
  filter(year %in% full_years) %>% 
  group_by(month_day) %>%
  summarise(
    n_years        = n_distinct(year),
    Pg_daily       = mean(Pg_BSC   * 0.021618 * 48, na.rm=TRUE),
    Pg_moss_daily  = mean(Pg_moss  * 0.021618 * 48, na.rm=TRUE),
    Pg_lich_daily  = mean(Pg_lich  * 0.021618 * 48, na.rm=TRUE),
    Pg_cyano_daily = mean(Pg_cyano * 0.021618 * 48, na.rm=TRUE),
    GPP_daily      = mean(GPP_DT   * 0.021618 * 48, na.rm=TRUE),
    .groups="drop") %>%
  filter(GPP_daily > gpp_threshold) %>%
  mutate(
    date       = as.Date(paste0("2016-", month_day)),
    BSC_frac   = Pg_daily       / GPP_daily * 100,
    moss_frac  = Pg_moss_daily  / GPP_daily * 100,
    lich_frac  = Pg_lich_daily  / GPP_daily * 100,
    cyano_frac = Pg_cyano_daily / GPP_daily * 100
  ) %>%
  filter(n_years >= min_years)   # drop sparse days entirely instead of NA-ing them

p_pct <- ggplot(pct_doy, aes(x=date)) +
  geom_rect(data=shade_monsoon,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_rect(data=shade_shoulder,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_line(aes(y=BSC_frac,   color="Total Biocrust"),           linewidth=2) +
  geom_line(aes(y=moss_frac,  color="Moss crust"),          linewidth=1.5) +
  geom_line(aes(y=lich_frac,  color="Lichen crust"),        linewidth=1.5) +
  geom_line(aes(y=cyano_frac, color="Cyanobacteria crust"), linewidth=1.5) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_date(date_labels="%b", date_breaks="1 month") +
  scale_fill_manual(values=c("Monsoon season"="lightskyblue2", "Shoulder season"="darkgrey"),
                    name="Season") +
  scale_color_manual(
    values=c("Total Biocrust"="black", "Moss crust"="#228833",
             "Cyanobacteria crust"="goldenrod3", "Lichen crust"="#4477AA"),
    breaks=c("Total Biocrust", "Moss crust", "Lichen crust", "Cyanobacteria crust")) +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  labs(x=NULL, y="Percent of Ecosystem GPP (%)", color=NULL, tag="a)") +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_blank(),
        legend.position="right",
        plot.tag=element_text(size=20, face="bold"),
        plot.margin=margin(4, 4, 0, 4))

bsc_daily_yearly <- model_output %>%
  mutate(date      = as.Date(datetime),
         year      = lubridate::year(date),
         month_day = format(date, "%m-%d")) %>%
  group_by(year, month_day) %>%
  summarise(
    Pg_daily  = mean(Pg_BSC * 0.021618 * 48, na.rm=TRUE),
    GPP_daily = mean(GPP_DT * 0.021618 * 48, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    date     = as.Date(paste0("2016-", month_day)),
    BSC_frac = Pg_daily / GPP_daily * 100,
    is_full  = year %in% full_years
  ) %>%
  filter(is.finite(BSC_frac), GPP_daily > gpp_threshold) %>%
  arrange(year, month_day) %>%
  group_by(year) %>% 
  mutate(doy=as.numeric(format(date, "%j")),
         gap=c(0,diff(doy)) > 3,
         segment=cumsum(gap)) %>% 
  ungroup()
  
bsc_daily_mean <- bsc_daily_yearly %>%
  filter(is_full) %>%
  group_by(month_day, date) %>%
  summarise(
    n_years     = sum(!is.na(BSC_frac)),
    mean_frac_ratio  = mean(BSC_frac, na.rm=TRUE),        # mean-of-ratios (for ribbon spread)
    sd_frac_ratio    = sd(BSC_frac, na.rm=TRUE),
    Pg_pooled   = mean(Pg_daily,  na.rm=TRUE),             # ratio-of-means (for central line)
    GPP_pooled  = mean(GPP_daily, na.rm=TRUE),
    .groups="drop"
  ) %>%
  mutate(
    central_frac = Pg_pooled / GPP_pooled * 100,           # matches panel a's method
    ribbon_low   = pmax(0, mean_frac_ratio - sd_frac_ratio),
    ribbon_high  = mean_frac_ratio + sd_frac_ratio
  ) %>%
  filter(n_years >= min_years)

yrs <- as.character(sort(unique(bsc_daily_yearly$year)))
# make a color palette with 14 very different hues
year_cols <- setNames(unname(palette.colors(length(yrs), palette = "Polychrome 36")), yrs)
year_ltys <- setNames(rep(c("solid", "dashed", "dotdash"), length.out = length(yrs)), yrs)

p_sens <- ggplot() +
  geom_rect(data=shade_monsoon,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_rect(data=shade_shoulder,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_line(data = bsc_daily_yearly,
            aes(x = date, y = BSC_frac, 
                group = interaction(year, segment),
                color = as.factor(year)),
            alpha = 0.7, linewidth = 0.8) +
  geom_ribbon(data=bsc_daily_mean,
              aes(x=date, ymin=ribbon_low, ymax=ribbon_high),
              fill="darkgreen", alpha=0.6) +
  geom_line(data=bsc_daily_mean,
            aes(x=date, y=central_frac),
            color="black", linewidth=1.5) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_date(date_labels="%b", date_breaks="1 month") +
  scale_color_manual(
    values=setNames(
      viridis::turbo(length(unique(bsc_daily_yearly$year)), begin=0.2, end=0.9),
      as.character(sort(unique(bsc_daily_yearly$year)))),
    name="Year") +
  guides(color = guide_legend(override.aes = list(linewidth = 4))) +
  scale_fill_manual(values=c("Monsoon season"="lightskyblue2", "Shoulder season"="darkgrey"), guide="none") +
  coord_cartesian(ylim=c(0, 30)) +
  labs(x="Month", y="Percent of Ecosystem GPP (%)", tag="b)") +
  theme_classic(base_size=20) +
  theme(legend.position="right",
        plot.tag=element_text(size=20, face="bold"),
        plot.margin=margin(0, 4, 4, 4))

fig5 <- p_pct / p_sens + plot_layout(heights=c(1, 1)) &
  theme(plot.tag=element_text(size=30, face="bold"),
        plot.tag.position=c(0.12, 0.95))
fig5

ggsave("KannanFig5.pdf", fig5, width=12, height=10, dpi=600)



# Supplemental Figure -----------------------------------------------------
## limiting factor analysis
crusts <- c(moss = "Moss", lich = "Lichen", cyano = "Cyanobacteria")

model_output <- model_output %>%
  ungroup() %>%
  mutate(
    # weighted temperature ceiling (absolute Pg units, normalized later by the 99th percentile)
    Pg_T_lich  = w_cyanolich * Pg_T_cyanolich + w_chlorolich * Pg_T_chlorolich,
    # weighted light fraction (already 0-1, so Pgmax_lich is 1 for get_scalars)
    Pg_I_lich  = w_cyanolich * LF_cyanolich + w_chlorolich * LF_chlorolich,
    Pgmax_lich = 1
  )

model_output <- model_output %>%
  ungroup() %>%
  mutate(
    Pg_T_lich = w_cyanolich * Pg_T_cyanolich + w_chlorolich * Pg_T_chlorolich,
    LF_lich   = w_cyanolich * LF_cyanolich   + w_chlorolich * LF_chlorolich
  )

get_scalars <- function(d, x) {
  d %>%
    transmute(datetime, year = year(datetime), month = month(datetime),
              crust = crusts[[x]],
              T = .data[[paste0("Pg_T_", x)]],
              W = .data[[paste0("WC_frac_", x)]],
              I = .data[[paste0("LF_", x)]])
}

day_data <- model_output %>% filter(period == "DAY", !is.na(PPFD_IN), PPFD_IN > 0)

lim <- bind_rows(lapply(names(crusts), function(x) get_scalars(day_data, x))) %>%
  filter(!is.na(T), !is.na(W), !is.na(I)) %>%
  group_by(crust) %>%
  mutate(across(c(T, W, I), ~ pmin(1, .x / quantile(.x, 0.99, na.rm = TRUE)))) %>%   # normalize
  ungroup() %>%
  mutate(f_min = pmin(T, W, I),
         limiting = c("Water", "Temperature", "Light")[max.col(-cbind(W, T, I), ties.method = "first")],
         md = month * 100 + day(datetime),
         season = case_when(md >= 615 & md <= 930 ~ "Monsoon",
                            month %in% c(3:5, 10:11) ~ "Shoulder",
                            TRUE ~ "Other"))

# overall: should be close to your 52.9% / 69.7% / 72.9% (they will move after the Pg fix)
lim %>% count(crust, limiting) %>% group_by(crust) %>% mutate(pct = round(100 * n / sum(n), 1))

# by season
lim %>% filter(season != "Other") %>%
  count(crust, season, limiting) %>%
  group_by(crust, season) %>% mutate(pct = round(100 * n / sum(n), 1))

lim_long <- lim %>% 
  pivot_longer(
    cols = c(W, T, I),
    names_to = "factor",
    values_to = "scalar"
  ) %>% 
  mutate(factor = recode(
    factor,
    W = "Water",
    T = "Temperature",
    I = "Light"),
    is_limiting = case_when(factor == limiting ~ TRUE,
                            TRUE ~ FALSE
    )
  )

pal <- c(Water = "#2196F3", Temperature = "#F44336", Light = "#FFC107")

# Plot A: share of daytime timesteps limited by each factor, by month
lim_pct <- lim %>% 
  count(year, crust, month, limiting, name="n") %>% 
  complete(year, crust, month, limiting, fill=list(n=0)) %>% 
  group_by(year, crust, month) %>% 
  filter(sum(n) > 0) %>% 
  mutate(pct = 100 * n / sum(n)) %>% 
  ungroup() %>% 
  group_by(crust, month, limiting) %>% 
  summarise(pct=mean(pct), .groups="drop")

make_lim_panel <- function(crust_name, title_str, tag_str) {
  ggplot(filter(lim_pct, crust == crust_name),
         aes(x=month, y=pct, color=limiting, group = limiting)) +
    geom_line(linewidth=1.5) +
    geom_point(size=3) +
    scale_x_continuous(breaks=1:12, labels=month.abb) +
    scale_y_continuous(limits=c(0,100), breaks=seq(0, 100, 25)) +
    scale_color_manual(values = pal) +
    labs(title=title_str, x="Month", y="Average daytime timesteps (%)", tag=tag_str) +
    theme_classic()+
    theme(legend.position = "none",
          plot.title   = element_text(face = "bold", size = 25, hjust = 0.5),
          axis.title.x = element_text(size = 25),
          axis.title.y = element_text(size = 20),
          axis.text    = element_text(size = 18),
          plot.tag     = element_text(size = 14, face = "bold"))
}

pa_moss <- make_lim_panel("Moss", "Moss-dominated", "c)")
pa_lich <- make_lim_panel("Lichen", "Lichen-dominated", "f)")
pa_cyano <- make_lim_panel("Cyanobacteria", "Cyanobacteria-dominated", "i)")

lim_order <- c("Temperature", "Water", "Light")

lim_legend_plot <- ggplot(lim_pct, aes(x = month, y = pct, color = limiting)) +
  geom_line(linewidth = 10) +
  scale_color_manual(values = pal, breaks = lim_order, name = NULL) +
  theme_classic() +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.text = element_text(size = 50), legend.justification = "left",
        legend.box.just = "left")
lim_legend <- cowplot::get_legend(lim_legend_plot)

legend_col <- cowplot::plot_grid(shared_legend, lim_legend, ncol = 1)


main_grid <- (p_temp | p4 | pa_moss) / (p_water | p6 | pa_lich) / (p_light | p5 | pa_cyano) +
  plot_annotation(theme = theme(plot.title = element_text(size = 16, face = "bold"))) &
  theme(plot.tag = element_text(size = 30, face = "bold"),
        plot.tag.position = c(0.18, 0.95))

legend_col <- cowplot::plot_grid(shared_legend, lim_legend, ncol = 1)
fig3 <- cowplot::plot_grid(legend_col, main_grid, ncol = 1, rel_heights = c(0.25, 1))
fig3
ggsave("KannanFig3.pdf", fig3, width = 28, height = 20, dpi = 600, bg = "white")



