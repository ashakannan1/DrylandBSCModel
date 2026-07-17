## Dryland Biocrust (BSC) Carbon Flux Model — Figure Scripts
## Asha Kannan
## Requires model output saved by model.R


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

setwd("/Users/ashakannan/Documents/GitHub/Dryland-Bryo-BSC-Model")

model_output  <- readRDS("data/model_output.rds")
model_daily   <- readRDS("data/model_daily.rds")
model_monthly <- readRDS("data/model_monthly.rds")
model_annual  <- readRDS("data/model_annual.rds")
model_input   <- readRDS("data/model_input.rds")


# Source model functions and parameters needed for ecophys curves
source("R/model.R")  # loads all parameters and functions without re-running model
                     # (comment out saveRDS lines if re-sourcing to avoid overwrite)

# Color palette
bsc_colors <- c(
  "Moss crust"          = "#228833",
  "Cyanobacteria crust" = "#EE6677",
  "Lichen crust"        = "#4477AA"
)

# Seasonal shading helpers
shade_monsoon <- data.frame(
  xmin = as.Date("2016-07-01"), xmax = as.Date("2016-08-15"),
  ymin = -Inf, ymax = Inf, season = "Monsoon season")
shade_shoulder <- data.frame(
  xmin = c(as.Date("2016-03-01"), as.Date("2016-08-23")),
  xmax = c(as.Date("2016-04-01"), as.Date("2016-10-15")),
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
  data.frame(type = "Open",    dx =  -5, dy =  25),
  data.frame(type = "Open",    dx = -15, dy =  15),
  data.frame(type = "Open",    dx =  10, dy =  20),
  # Juniper sensors
  data.frame(type = "Juniper", dx = -30, dy =  10),
  data.frame(type = "Juniper", dx = -25, dy = -10),
  data.frame(type = "Juniper", dx = -20, dy =  25),
  # Pinyon sensors
  data.frame(type = "Pinyon",   dx =  25, dy = -15),
  data.frame(type = "Pinyon",   dx =  30, dy =  10),
  data.frame(type = "Pinyon",   dx =  20, dy = -25),
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
tiles <- get_tiles(bbox, provider = "Esri.WorldImagery", zoom=18, crop=TRUE)



# Build the figure!
fig_map <- ggplot() +
  # satellite basemap
  geom_spatraster_rgb(data=tiles) +
  
  # biocrust cover overlay (alpha to show veg structure
  geom_raster(data=bsc_df %>% filter(!is.na(biocrust)),
              aes(x=lon, y=lat, fill=biocrust), alpha=0.6) +
  scale_fill_gradient(low="transparent", high="mediumpurple4",
                      name="Biocrust\ncover fraction",
                      na.value="transparent") +
  
  # footprint contours -- 50% and 80% as separate layers so each gets a legend label
  geom_contour(data=fpt_df %>% filter(!is.na(weight)),
               aes(x=lon, y=lat, z=weight, linetype="50% footprint"),
               breaks=level_50,
               color="black", linewidth=3, alpha=0.6) +
  geom_contour(data=fpt_df %>% filter(!is.na(weight)),
               aes(x=lon, y=lat, z=weight, linetype="80% footprint"),
               breaks=level_80,
               color="black", linewidth=3, alpha=0.6) +
  
  # transect lines
  geom_sf(data=transect_lines_sf,
          aes(linetype="Transect"),
          color="black", linewidth=3) +
  
  # combined linetype legend for contours + transects
  scale_linetype_manual(
    name="Footprint / Transects",
    values=c("50% footprint"="solid",
             "80% footprint"="dashed",
             "Transect"="dotted"),
    guide=guide_legend(override.aes=list(linewidth=1.2, color="black"))
  ) +
  
  # flux tower and sensors
  geom_sf(data=sensors_sf,
          aes(shape=type, color=type), fill="white",
          size=14, stroke=3) +
  scale_color_manual(
    values=c(
      "Open" = "salmon4",
      "Juniper" = "maroon3",
      "Pinyon" = "darkblue",
      "Flux Tower" = "black"
    ),
    name="Site sensors"
  ) +
  scale_shape_manual(values=c(
    "Open" = 0,
    "Juniper" = 1,
    "Pinyon" = 2,
    "Flux Tower" = 24
  ),
  name="Site sensors") +
  # north arrow and scale bar via annotation
  annotation_scale(location = "br", width_hint = 0.3,
                   text_col = "white", line_col = "white",
                   bar_cols = c("white", "grey50"),
                   text_cex=1.5) +
  annotation_north_arrow(location = "tl",
                         style = north_arrow_minimal(text_col = "white",
                                                     line_col = "white",
                                                     fill = "white",
                                                     text_size=20)) +
  
  coord_sf(
    xlim = c(tower_lon - buf_deg_lon, tower_lon + buf_deg_lon),
    ylim = c(tower_lat - buf_deg_lat, tower_lat + buf_deg_lat),
    expand = FALSE
  ) +
  
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(
    axis.text = element_text(size=20),
    axis.title = element_text(size=30),
    legend.text = element_text(size=30),
    legend.title = element_text(size=30),
    legend.background = element_rect(fill = "white", color=NA),
    legend.position = "right",
    legend.key.height=unit(2, "cm")
  )

ggsave("fig_site_map.png", fig_map, width = 20, height = 20, dpi = 300)

  
  

# Figure 3 — Ecophysiological response curves + site environment ----------

temp_seq <- seq(0, 45, 1)
temp_curve_all <- tibble(temp = temp_seq) %>%
  rowwise() %>%
  mutate(
    Moss  = pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=temp, T_k=T_k_moss,  Pgmax=Pgmax_moss) *
                   f_H(Tv=temp, delta_S=delta_S_moss,  delta_Hd=delta_Hd_moss)),
    Cyano = pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=temp, T_k=T_k_cyano, Pgmax=Pgmax_cyano) *
                   f_H(Tv=temp, delta_S=delta_S_cyano, delta_Hd=delta_Hd_cyano)),
    Lich  = pmax(0, calc_Pg_fxn_T_Tk_rowwise(temp=temp, T_k=T_k_lich,  Pgmax=Pgmax_lich) *
                   f_H(Tv=temp, delta_S=delta_S_lich,  delta_Hd=delta_Hd_lich))
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
    Lich  = calc_RD_fxn_T_rowwise(temp=temp, r_r=r_r_lich,  R_o=RD_min_lich)
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
    `Lichen crust`        = calc_PN_fxn_I_rowwise(PAR=PAR, f_Io_Icomp=f_Io_Icomp_lich,
                                                  b=b_lich,  g=g_lich,  Icomp=Icomp_lich)  * lich_light_scale
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

env_monthly <- model_input %>%
  filter(PPFD_IN > 0) %>%
  group_by(month) %>%
  summarise(
    mean_T   = mean(SOILT_O_5, na.rm=TRUE),
    min_T    = min(SOILT_O_5,  na.rm=TRUE),
    max_T    = max(SOILT_O_5,  na.rm=TRUE),
    mean_PAR = mean(PPFD_IN,   na.rm=TRUE),
    min_PAR  = min(PPFD_IN,    na.rm=TRUE),
    max_PAR  = max(PPFD_IN,    na.rm=TRUE),
    mean_SWC = mean(SWC_O_5,   na.rm=TRUE),
    min_SWC  = min(SWC_O_5,    na.rm=TRUE),
    max_SWC  = max(SWC_O_5,    na.rm=TRUE),
    .groups = "drop"
  )

p_temp <- ggplot(temp_curve_all, aes(x=temp, y=Pg_T, color=type)) +
  geom_line(data=rd_curve_all, aes(x=temp, y=RD_T, color=type),
            linewidth=1.5, linetype="dashed") +
  geom_line(linewidth=2) +
  scale_color_manual(values=bsc_colors) +
  labs(x="Temperature (°C)", y=expression(CO[2]~exchange~(µmol~CO[2]~m^-2~s^-1)),
       color=NULL, tag="a)") +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=25),
        axis.text    = element_text(size=18),
        plot.tag     = element_text(size=14, face="bold"))

p_light <- ggplot(light_curve_all, aes(x=PAR, y=Pg_I, color=type)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=bsc_colors) +
  labs(x=expression("PPFD"~("µmol photons"~m^-2~s^-1)),
       y=expression(P[net]~(µmol~CO[2]~m^-2~s^-1)),
       color=NULL, tag="c)") +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text    = element_text(size=18),
        plot.tag     = element_text(size=14, face="bold"))

p_water <- ggplot(water_curve_all, aes(x=WC, y=WC_frac, color=type)) +
  geom_line(linewidth=2) +
  scale_color_manual(values=bsc_colors) +
  labs(x="Relative Water Content", y=expression(P[frac]),
       color=NULL, tag="e)") +
  theme_classic() +
  theme(legend.position="none",
        axis.title.x = element_text(size=30),
        axis.title.y = element_text(size=30),
        axis.text    = element_text(size=18),
        plot.tag     = element_text(size=14, face="bold"))

p4 <- ggplot(env_monthly, aes(x=month)) +
  geom_ribbon(aes(ymin=min_T, ymax=max_T), fill="#F44336", alpha=0.2) +
  geom_line(aes(y=mean_T), color="#F44336", linewidth=2) +
  geom_point(aes(y=mean_T), color="#F44336", size=3) +
  scale_x_continuous(breaks=1:12, labels=month.abb) +
  labs(x=NULL, y="Mean Soil Temperature (°C)", tag="b)") +
  theme_classic() +
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=18),
        plot.tag   = element_text(size=14, face="bold"))

p5 <- ggplot(env_monthly, aes(x=month)) +
  geom_ribbon(aes(ymin=min_PAR, ymax=max_PAR), fill="#FFC107", alpha=0.2) +
  geom_line(aes(y=mean_PAR), color="#FFC107", linewidth=2) +
  geom_point(aes(y=mean_PAR), color="#FFC107", size=3) +
  scale_x_continuous(breaks=1:12, labels=month.abb) +
  labs(x=NULL, y=expression("Mean Daytime PPFD"~("µmol photons"~m^-2~s^-1)), tag="d)") +
  theme_classic() +
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=18),
        plot.tag   = element_text(size=14, face="bold"))

p6 <- ggplot(env_monthly, aes(x=month)) +
  geom_ribbon(aes(ymin=min_SWC, ymax=max_SWC), fill="#2196F3", alpha=0.2) +
  geom_line(aes(y=mean_SWC), color="#2196F3", linewidth=2) +
  geom_point(aes(y=mean_SWC), color="#2196F3", size=3) +
  scale_x_continuous(breaks=1:12, labels=month.abb) +
  labs(x=NULL, y="Mean Soil Water Content (m³/m³)", tag="f)") +
  theme_classic() +
  theme(axis.title = element_text(size=20),
        axis.text  = element_text(size=18),
        plot.tag   = element_text(size=14, face="bold"))

legend_plot <- ggplot(temp_curve_all, aes(x=temp, y=Pg_T, color=type)) +
  geom_line(linewidth=10) +
  scale_color_manual(values=bsc_colors, name=NULL) +
  theme_classic() +
  theme(legend.position = "bottom",
        legend.text = element_text(size=50))
shared_legend <- cowplot::get_legend(legend_plot)

main_grid <- (p_temp | p4) / (p_light | p5) / (p_water | p6) +
  plot_annotation(theme=theme(plot.title=element_text(size=16, face="bold"))) &
  theme(plot.tag=element_text(size=30, face="bold"),
        plot.tag.position=c(0.18, 0.95))

fig3 <- cowplot::plot_grid(main_grid, shared_legend, ncol=1, rel_heights=c(1, 0.05))
fig3
ggsave("figures/fig3_ecophys.png", fig3, width=20, height=20, dpi=300, bg="white")


# Figure 4 — BSC carbon flux timeseries -----------------------------------

ecosystem_doy <- model_daily %>%
  mutate(doy = yday(date)) %>%
  filter(GPP_daily > -1) %>%
  group_by(doy) %>%
  summarise(GPP_mean = mean(GPP_daily, na.rm=TRUE), .groups="drop")

bsc_doy <- model_daily %>%
  mutate(doy = yday(date)) %>%
  group_by(doy) %>%
  summarise(
    GPP_mean = mean(Pg_BSC_daily,   na.rm=TRUE),
    R_mean   = mean(RD_BSC_daily,   na.rm=TRUE),
    NEE_mean = mean(Pnet_BSC_daily, na.rm=TRUE),
    .groups="drop"
  )

bsc_org_doy <- model_daily %>%
  mutate(doy = yday(date)) %>%
  group_by(doy) %>%
  summarise(
    Moss  = mean(Pg_moss_daily,  na.rm=TRUE),
    Cyano = mean(Pg_cyano_daily, na.rm=TRUE),
    Lich  = mean(Pg_lich_daily,  na.rm=TRUE),
    .groups="drop"
  ) %>%
  pivot_longer(-doy, names_to="type", values_to="Pg") %>%
  mutate(type = recode(type,
                       Moss  = "Moss crust",
                       Cyano = "Cyanobacteria crust",
                       Lich  = "Lichen crust"))

doy_breaks <- c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335)

p_a <- ggplot(ecosystem_doy, aes(x=doy, y=GPP_mean)) +
  geom_line(aes(color="Ecosystem GPP"), linewidth=1.2) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_continuous(breaks=doy_breaks, labels=month.abb) +
  scale_color_manual(values=c("Ecosystem GPP"="burlywood4")) +
  labs(x=NULL, y=expression("GPP"~(gC~m^-2~d^-1)), tag="a)") +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.title=element_blank(),
        plot.tag=element_text(size=14, face="bold"),
        plot.margin=margin(4, 4, 0, 4))

p_b <- ggplot(bsc_doy, aes(x=doy)) +
  geom_line(aes(y=GPP_mean,  color="BSC GPP"), linewidth=1.0, alpha=0.8) +
  geom_line(aes(y=-R_mean,   color="BSC R"),   linewidth=1.0, alpha=0.8) +
  geom_line(aes(y=NEE_mean,  color="BSC NEE"), linewidth=1.2) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_continuous(breaks=doy_breaks, labels=month.abb) +
  scale_color_manual(values=c("BSC GPP"="darkgreen", "BSC R"="tomato3", "BSC NEE"="black")) +
  labs(x=NULL, y=expression("Carbon exchange"~(gC~m^-2~d^-1)), color=NULL, tag="b)") +
  theme_classic(base_size=13) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="right",
        plot.tag=element_text(size=14, face="bold"),
        plot.margin=margin(0, 4, 0, 4))

p_c <- ggplot(bsc_org_doy, aes(x=doy, y=Pg, color=type)) +
  geom_line(linewidth=1.0) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_continuous(breaks=doy_breaks, labels=month.abb) +
  scale_color_manual(values=bsc_colors) +
  labs(x="Day of year", y=expression("BSC GPP"~(gC~m^-2~d^-1)), color=NULL, tag="c)") +
  theme_classic(base_size=13) +
  theme(legend.position="right",
        plot.tag=element_text(size=14, face="bold"),
        plot.margin=margin(0, 4, 4, 4))

fig4 <- p_a / p_b / p_c + plot_layout(heights=c(1, 1, 1))
fig4
ggsave("figures/fig4_timeseries.png", fig4, width=10, height=12, dpi=300)


# Figure 5 — BSC percent of ecosystem GPP --------------------------------

pct_doy <- model_output %>%
  mutate(date      = as.Date(datetime),
         month_day = format(date, "%m-%d")) %>%
  group_by(month_day) %>%
  summarise(
    Pg_daily       = mean(Pg_BSC   * 0.021618 * 48, na.rm=TRUE),
    Pg_moss_daily  = mean(Pg_moss  * 0.021618 * 48, na.rm=TRUE),
    Pg_lich_daily  = mean(Pg_lich  * 0.021618 * 48, na.rm=TRUE),
    Pg_cyano_daily = mean(Pg_cyano * 0.021618 * 48, na.rm=TRUE),
    GPP_daily      = mean(GPP_DT   * 0.021618 * 48, na.rm=TRUE),
    .groups="drop") %>%
  filter(GPP_daily > 0.01) %>%
  mutate(
    date       = as.Date(paste0("2016-", month_day)),
    BSC_frac   = Pg_daily       / GPP_daily * 100,
    moss_frac  = Pg_moss_daily  / GPP_daily * 100,
    lich_frac  = Pg_lich_daily  / GPP_daily * 100,
    cyano_frac = Pg_cyano_daily / GPP_daily * 100
  )

p_pct <- ggplot(pct_doy, aes(x=date)) +
  geom_rect(data=shade_monsoon,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_rect(data=shade_shoulder,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_line(aes(y=BSC_frac,   color="Total BSC"),           linewidth=2) +
  geom_line(aes(y=moss_frac,  color="Moss crust"),          linewidth=1.5) +
  geom_line(aes(y=lich_frac,  color="Lichen crust"),        linewidth=1.5) +
  geom_line(aes(y=cyano_frac, color="Cyanobacteria crust"), linewidth=1.5) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_date(date_labels="%b", date_breaks="1 month") +
  scale_fill_manual(values=c("Monsoon season"="lightskyblue2", "Shoulder season"="darkgrey")) +
  scale_color_manual(
    values=c("Total BSC"="black", "Moss crust"="#228833",
             "Cyanobacteria crust"="#EE6677", "Lichen crust"="#4477AA"),
    breaks=c("Total BSC", "Moss crust", "Lichen crust", "Cyanobacteria crust")) +
  labs(x=NULL, y="Percent of Ecosystem GPP (%)", color=NULL, tag="a)") +
  theme_classic(base_size=20) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        legend.position="right",
        plot.tag=element_text(size=20, face="bold"),
        plot.margin=margin(4, 4, 0, 4))

# Interannual sensitivity ribbon
full_years <- c(2010, 2011, 2013, 2015, 2017, 2018, 2019, 2020, 2021)

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
  filter(is.finite(BSC_frac), GPP_daily > 0.5, BSC_frac > 0) %>%
  arrange(year, month_day) %>%
  group_by(year) %>%
  mutate(
    date_numeric = as.numeric(as.Date(paste0("2016-", month_day))),
    gap          = c(0, diff(date_numeric)),
    BSC_frac     = ifelse(gap > 3, NA, BSC_frac)  # break line if gap > 3 days
  ) %>%
  ungroup()

bsc_daily_mean <- bsc_daily_yearly %>%
  filter(is_full) %>%
  group_by(month_day, date) %>%
  summarise(
    mean_frac   = mean(BSC_frac, na.rm=TRUE),
    ribbon_low  = pmax(0, mean_frac - sd(BSC_frac, na.rm=TRUE)),
    ribbon_high = mean_frac + sd(BSC_frac, na.rm=TRUE),
    .groups="drop"
  )

p_sens <- ggplot() +
  geom_rect(data=shade_monsoon,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_rect(data=shade_shoulder,
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=season),
            alpha=0.3, inherit.aes=FALSE) +
  geom_line(data=bsc_daily_yearly,
            aes(x=date, y=BSC_frac, group=year, color=as.factor(year)),
            alpha=0.4, linewidth=0.8) +
  geom_ribbon(data=bsc_daily_mean,
              aes(x=date, ymin=ribbon_low, ymax=ribbon_high),
              fill="darkgreen", alpha=0.6) +
  geom_line(data=bsc_daily_mean,
            aes(x=date, y=mean_frac),
            color="black", linewidth=1.5) +
  geom_hline(yintercept=0, linetype="dotted", color="gray50") +
  scale_x_date(date_labels="%b", date_breaks="1 month") +
  scale_color_manual(
    values=setNames(
      viridis::mako(length(unique(bsc_daily_yearly$year)), begin=0.2, end=0.9),
      as.character(sort(unique(bsc_daily_yearly$year)))),
    name="Year") +
  scale_fill_manual(values=c("Monsoon season"="lightskyblue2", "Shoulder season"="darkgrey")) +
  coord_cartesian(ylim=c(0, 30)) +
  labs(x=NULL, y="Percent of Ecosystem GPP (%)", tag="b)") +
  theme_classic(base_size=20) +
  theme(legend.position="right",
        plot.tag=element_text(size=20, face="bold"),
        plot.margin=margin(0, 4, 4, 4))

fig5 <- p_pct / p_sens + plot_layout(heights=c(1, 1))
fig5
ggsave("figures/fig5_percent_gpp.png", fig5, width=12, height=10, dpi=300)


