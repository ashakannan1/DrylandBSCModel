## Title: Angular Interpolation Methods — MODEL BUILDING
## Created by: Asha Kannan
## Created on: 4-8-2026
## Split: model-building logic only (see 02_figures.R for plots/animations/stats)

# =============================================================
# Section 1: Libraries + Packages
# =============================================================
## This script uses here::here() to set working directory. Open via its .Rproj file (or run from project root)
## and paths below should resolve automatically
## Expected folders: data/mpj_transects.xlsx (LPI data) and data/monthly_footprint_climatology_weight_map (Chu et al. 2021 .tif files)
## Run Section 1 to install/load required packages

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(terra)
library(stringr)
library(av)

# here::here() anchors all paths to the project root so this script will run on the same machine/OS without editing setwd() calls
#library(here)
# or, set to your local wd
# setwd("/path/to/your/wd/where/transectdata/stored")

# =============================================================
# Section 2: Site Constants & Coordinate References
# =============================================================
# The flux tower at US-Mpj is the spatial origin for all analyses.
# All transect coordinates, footprint raster, and vegetation grids are defined
# relative to this location.
# =============================================================
# tower coordinates
tower_lat <- 34.4385
tower_lon <- -106.2377

# tower location as df for ggplot plotting
tower_df <- data.frame(x=tower_lon, y=tower_lat)

# unit conversions at this latitude
# meters per degree stays constant for latitude
# meters per degree longitude shrinks with cos(lat) toward poles
m_per_deg_lat <- 111320 # stay constant
m_per_deg_lon <- 111320 * cos(tower_lat * pi/180)


# =============================================================
# Section 3: Transect Definitions
# =============================================================
# 3 transects radiate from flux tower at fixed compass azimuths
# each transect extends 100m from tower and was surveyed using LPI methodology
# at 0.2m intervals
#
# T1: 180 (due south)
# T2: 240 (SW)
# T3: 315 (NW)
#
# Coordinate convention:
# transect end points computed in standard compass math
# x_m = distance * sin(azimuth) --> positive = east
# y_m = distance * cos(azimuth) --> positive = north
#
# =============================================================
transect_info <- data.frame(
  sheet = c("Transect 1", "Transect 2", "Transect 3"),
  label = c("South", "Southwest", "Northwest"),
  azimuth = c(180, 240, 315),
  stringsAsFactors = FALSE
)

# Compute transect endpoints in lon/lat for plotting
transect_lines <- transect_info %>% 
  mutate(
    # endpoint in meters from tower using true compass azimuth
    end_x_m = 100 * sin(azimuth * pi / 180),
    end_y_m = 100 * cos(azimuth * pi / 180),
    
    # start at raster-center tower position
    x_start = tower_lon,
    y_start = tower_lat,
    
    # convert to lon/lat
    # compass directions for S, SW, NW azimuths
    x_end = tower_lon + end_x_m / m_per_deg_lon,
    y_end = tower_lat + end_y_m / m_per_deg_lat,
  )

# =============================================================
# Section 4: Transect Data Loading & Cover Fraction Calculation
# =============================================================
# Purpose: read veg LPI survey data from 3 transects, classify each intercept 
# point into a functional group, and compute percent cover fractions binned into 8m
# intervals matching the footprint raster pixel size
#
# Input: mpj_transects.xlsx (expected at data/mpj_transects.xlsx)
# - one sheet per transect
# - cols A:K contain dist and veg data
# - cols L + ignored
#
# Output: cover_fractions
# - one row per transect per 8m dist bin
# - cols for each functional group cover fraction (0-1)
# =============================================================
# -------------------------------------------------------------------------
# FUNCTION: read_transect
# reads one transect sheet from the excel file and converts each veg point to (x,y) coordinates
# in meters relative to the tower
#
# Arguments: sheet_name: name of excel sheet; azimuth: compass bearing of transect in degrees
# Returns: df with one row per veg point every 0.2m along transect, cols for distance, 
# all veg layers, transect ID, azimuth, and x/y coordinates

read_transect <- function(sheet_name, azimuth) {
  df <- read_excel(
    path = "mpj_transects.xlsx",
    sheet = sheet_name,
    range=cell_cols("A:K") # ignore cols L and beyond
  )
  
  # rename columns to clean names
  colnames(df) <- c("dist_m", "canopy", "canopy2", "top_veg", "veg2", "veg3", "other_veg", "soil_surf", "soil2")
  
  # drop any rows where dist_M is NA (empty rows at bottom of excel file)
  df <- df[!is.na(df$dist_m),]
  
  # add transect id and azimuth (metadata)
  df$transect <- sheet_name
  df$azimuth <- azimuth
  
  # convert polar transect (distance + azimuth) into local Cartesian coordinates centered at flux tower (0,0)
  # x = distance * sin(azimuth): positive = east
  # y = distance * cos(azimuth): positive = north
  az_rad <- azimuth * pi/180
  df$x_m <- df$dist_m * sin(az_rad)
  df$y_m <- df$dist_m * cos(az_rad)
  
  return(df)
}

# read all 3 transects & stack within the same df
# drop empty cols from Excel explor
transects <- bind_rows(
  read_transect("Transect 1", 180),
  read_transect("Transect 2", 240),
  read_transect("Transect 3", 315)
) %>% 
  select(-`...10`, -`...11`)

## check: should have 1500 rows (500 points * 3 transects) and 13 cols
head(transects)

# Canopy Code Standardization ---------------------------------------------
# some codes have multiple versions in raw data
# standardize to consistent labels before classifying

## standardize canopy codes
clean_canopy <- function(x) {
  case_when(
    x == "JUNIP canopy"          ~ "JUNIP canopy",
    x == "PINP canopy"           ~ "PINP canopy",
    x == "PJ canopy"             ~ "PJ canopy",
    x == "pinyon-j combo canopy" ~ "PJ canopy",
    TRUE                         ~ x
  )
}
transects <- transects %>%
  mutate(
    canopy  = clean_canopy(canopy),
    canopy2 = clean_canopy(canopy2)
  ) 

# Functional group codes --------------------------------------------------
# define classification groups based on IDs used in the LPI excel file
cyano_codes <- c("dkcy", "ltcy")
lich_codes <- c("LICH1", "LICH2", "PSOR", "TONI", "LICH3", "LICH4", "PELT")
moss_codes <- c("MOSS1", "MOSS3")

# define microhabitat of each crust type
transects <- transects %>% 
  mutate(
    # crust type: check soil_surf and soil2 cols
    # a point can have multiple crust types in diff layers
    # we assign priority based on percent weighting: cyano > moss > lich
    has_lich = soil_surf %in% lich_codes | soil2 %in% lich_codes,
    has_moss = soil_surf %in% moss_codes | soil2 %in% moss_codes,
    has_cyano = soil_surf %in% cyano_codes | soil2 %in% cyano_codes,
    
    # microhab: check canopy and canopy2 cols
    # open = no canopy overhead
    microhabitat = case_when(
      canopy == "PINP canopy" | canopy2 == "PINP canopy" ~ "under_pinp",
      canopy == "JUNIP canopy" | canopy2 == "JUNIP canopy" ~ "under_junip",
      canopy == "PJ canopy" | canopy2 == "PJ canopy" ~ "under_pj",
      TRUE ~ "open"
    )
  )

# define what codes correspond to diff functional groups
crust_codes <- c("dkcy", "ltcy", "LICH1", "LICH2", "PSOR", "TONI", "LICH3", "LICH4", "PELT", "MOSS1", "MOSS3")
midveg_codes <- c("grass", "OPUN", "yucca", "snake", "LEAF")
subcanopy_codes <- c("PINJ", "JUNIP", "OAK")
canopy_codes <- c("PINJ canopy", "JUNIP canopy", "PJ canopy")
exclude_codes <- c("EQUIP", "treeLICH")
abiotic_codes <- c("BARE", "EQUIP", "DUFF", "GRAV", "COBB", "BOUL", "STON")

# function to check if any column in a row matches a set of codes
any_matches <- function(df, cols, codes) {
  rowSums(sapply(df[cols], function(col) col %in% codes), na.rm=TRUE) > 0
}

# all veg cols to search across
veg_cols <- c("canopy", "canopy2", "top_veg", "veg2", "veg3", "other_veg",
              "soil_surf", "soil2")

# total points (denominator)
n_total <- nrow(transects)

# count occurrences of each code in surface layers
surface_cover <- transects %>%
  mutate(surface = ifelse(!is.na(soil_surf), soil_surf, soil2)) %>%
  count(surface) %>%
  mutate(pct_cover = n / n_total * 100) %>%
  arrange(desc(pct_cover))

sum(surface_cover$pct_cover)  # should be 100%
print(surface_cover)

# add binary presence columns for each functional group - does this sect include this veg type
transects <- transects %>% 
  mutate(
    is_cyano = any_matches(., veg_cols, cyano_codes),
    is_lich = any_matches(., veg_cols, lich_codes),
    is_moss = any_matches(., veg_cols, moss_codes),
    is_midveg = any_matches(., veg_cols, midveg_codes),
    is_subcanopy = any_matches(., veg_cols, subcanopy_codes),
    is_canopy = any_matches(., veg_cols, canopy_codes),
    is_exclude = any_matches(., veg_cols, exclude_codes),
    is_abiotic = any_matches(., veg_cols, abiotic_codes))

# how many crust points fall in each crust type × microhabitat combination?
transects %>%
  filter(has_cyano | has_lich | has_moss) %>%
  mutate(crust_type = case_when(
    has_cyano   ~ "cyano",
    has_lich ~ "lich",
    has_moss  ~ "moss"
  )) %>%
  mutate(has_midveg = is_midveg,
         has_canopy = is_canopy) %>% 
  count(crust_type, has_midveg, has_canopy)

# summarize the total number of veg types per row
transects %>%
  summarise(
    n_total      = n(),
    n_cyano      = sum(is_cyano),
    n_lich       = sum(is_lich),
    n_moss       = sum(is_moss),
    n_midveg     = sum(is_midveg),
    n_subcanopy  = sum(is_subcanopy),
    n_canopy     = sum(is_canopy),
    n_exclude    = sum(is_exclude),
    n_abiotic   = sum(is_abiotic)
  )

# Distance binning --------------------------------------------------------
# bin each transect point into 8m intervals for matching footprint pixel size
# bin_center stores midpoint distance of each bin for spatial lookup
# define bin width to match footprint pixel size 
bin_width <- 8

# create distance bins along each transects; breaks go from 0-100m in 8 m steps
# cut() assigns each point to its bin
# labels = F gives numeric bin index instead of interval string
transects <- transects %>% 
  mutate(
    dist_bin = cut(dist_m,
                   breaks=seq(0, 104, by=bin_width), ## we use 104 to capture all points thru 100
                   labels=FALSE,
                   include.lowest=TRUE),
    
    # also store the bin center distance for later spatial analysis
    bin_center = (dist_bin - 0.5) * bin_width
  )

# Cover Fraction Calculation ----------------------------------------------
# for each transect x distance bin, compute the fraction of included points 
# belonging to each functional group
#
# note: cover fractions do NOT sum to 1 because multiple vertical layers can occur
# at a single point (eg BSC under canopy) and both can contribute to ecosystem GPP
# at the same location

cover_fractions <- transects %>% 
  filter(!is_exclude) %>%  ## remove equipment points
  group_by(transect, azimuth, dist_bin, bin_center) %>% 
  summarise(
    n_points = n(), ## total valid points found in a bin
    cover_cyano = sum(is_cyano) / n(),
    cover_lich = sum(is_lich) / n(),
    cover_moss = sum(is_moss) / n(),
    cover_midveg = sum(is_midveg) / n(),
    cover_subcanopy = sum(is_subcanopy) / n(),
    cover_canopy = sum(is_canopy) / n(),
    cover_abiotic = sum(is_abiotic) / n(),
    .groups="drop"
  )

# verification: check bin counts per transect (expect ~13 bins each)
cover_fractions %>%
  group_by(transect) %>%
  summarise(n_bins = n())

# =============================================================
# Section 5: Vegetation IDW Raster Builder
# =============================================================
# Purpose: for any given footprint raster grid, build a matching veg cover raster
# using angular inverse distance weighting (IDW) across the 3 transects
# 
# For each pixel in the footprint grid, we want to 
# 1. calculate its distance and bearing from the tower
# 2. look up cover fractions at that distance bin from each of the transects
# 3. weight each transect's cover by inverse angular dist. (closer trnansect
# azimuth = higher weight)
# 4. compute weighted average cover for each functional group
# 
## Pixels beyond 100 meters use the outermost distance bin (bin 13) as a constant extrapolation
#
# Arguments: fpt_template (a SpatRaster whose grid defines the output extent, resolution,
# and coordinate reference system)
#
# Returns: a 7-layer SpatRaster with 1 layer per functional group
# Dependencies: cover_fractions, tower_lat, tower_lon, m_per_deg_lon, m_per_deg_lat
# =============================================================
# FUNCTION: Build Vegetation Raster ---------------------------------------
build_veg_idw_raster <- function(fpt_template) {
  # get pixel center coordinates from footprint grid
  # df <- as.data.frame(fpt_template, xy=TRUE)
  # colnames(df) <- c("lon", "lat", "fpt_val")
  
  # use crds() instead of as.data.frame() so it works on empty template rasters
  coords <- crds(fpt_template, na.rm=FALSE)
  df <- data.frame(lon=coords[,1], lat=coords[,2], fpt_val = NA)
  
  # convert pixel centers to meters relative to tower
  # calculate bearing from tower to each pixel
  pixel_coords <- df %>% 
    mutate(
      # convert to meters relative to tower
      x_m = (lon - tower_lon) * m_per_deg_lon, #e-w
      y_m = (lat - tower_lat) * m_per_deg_lat, #n-s
      
      # distance from tower in meters
      dist_m = sqrt(x_m^2 + y_m^2),
      
      # compass bearing from tower (0=N, clockwise)
      bearing=(atan2(x_m, y_m)*180 / pi+ 360) %% 360,
      
      # assign distance bin, capped at 13 for pixels beyond 100 meters
      dist_bin = pmin(ceiling(dist_m / 8), 13))
  
  # FUNCTION: Angular IDW weighting -----------------------------------------
  ## for each pixel, calculate the angular distance to each transect azimuth
  ## and calculate inverse distance weights

  # Angular distance is the shortest arc between 2 bearings on a 360 deg circle
  # eg: angular_diff (350, 10) = 20, not 340
  angular_diff <- function(bearing, azimuth) {
    diff <- abs(bearing - azimuth) %% 360
    pmin(diff, 360-diff) # take shorter arc
  }
  
  # for each pixel, calculate angular distance to each transect and compute normalized
  # inverse distance weights
  # enforce min angular dist of 0.1 deg to avoid division by 0 when pixel falls exactly
  # on a transect azimuth
  pixel_coords <- pixel_coords %>% 
    mutate(
      # angular dist to each transect in degs
      ang_dist_t1 = pmax(angular_diff(bearing, 180), 0.1),
      ang_dist_t2 = pmax(angular_diff(bearing, 240), 0.1),
      ang_dist_t3 = pmax(angular_diff(bearing, 315), 0.1),
      
      # inverse angular distance weights (closer = higher weight)
      w1_raw = 1/ang_dist_t1,
      w2_raw = 1/ang_dist_t2,
      w3_raw = 1/ang_dist_t3,
      
      # normalize weights so they sum to 1 for each pixel
      w_sum = w1_raw + w2_raw + w3_raw,
      w1 = w1_raw / w_sum,
      w2 = w2_raw / w_sum,
      w3 = w3_raw / w_sum
    )
  
  # Cover Fraction Lookup ---------------------------------------------------
  # for each pixel, look up cover fraction at its distance bin from each of the 3 
  # transects separately, then combine using the angular IDW weights
  
  # pull cover fractions for each transect at each distance bin
  cf_t1 <- cover_fractions %>%
    filter(transect == "Transect 1") %>%
    select(dist_bin, cover_cyano, cover_lich, cover_moss, cover_midveg, 
           cover_subcanopy, cover_canopy, cover_abiotic)
  cf_t2 <- cover_fractions %>%
    filter(transect == "Transect 2") %>%
    select(dist_bin, cover_cyano, cover_lich, cover_moss, cover_midveg, 
           cover_subcanopy, cover_canopy, cover_abiotic)
  
  cf_t3 <- cover_fractions %>%
    filter(transect == "Transect 3") %>%
    select(dist_bin, cover_cyano, cover_lich, cover_moss, cover_midveg, 
           cover_subcanopy, cover_canopy, cover_abiotic)
  
  # join each transect's cover fractions to pixel data by dist_bin
  pixel_coords <- pixel_coords %>% 
    left_join(cf_t1, by="dist_bin", suffix=c("", "_t1")) %>% 
    rename(
      cyano_t1 = cover_cyano,
      lich_t1 = cover_lich,
      moss_t1 = cover_moss,
      mid_t1 = cover_midveg,
      sub_t1 = cover_subcanopy,
      can_t1 = cover_canopy,
      abi_t1 = cover_abiotic
    ) %>% 
    left_join(cf_t2, by="dist_bin", suffix=c("", "_t2")) %>% 
    rename(
      cyano_t2 = cover_cyano,
      lich_t2 = cover_lich,
      moss_t2 = cover_moss,
      mid_t2 = cover_midveg,
      sub_t2 = cover_subcanopy,
      can_t2 = cover_canopy,
      abi_t2 = cover_abiotic
    ) %>% 
    left_join(cf_t3, by="dist_bin", suffix=c("", "_t3")) %>% 
    rename(
      cyano_t3 = cover_cyano,
      lich_t3 = cover_lich,
      moss_t3 = cover_moss,
      mid_t3 = cover_midveg,
      sub_t3 = cover_subcanopy,
      can_t3 = cover_canopy,
      abi_t3 = cover_abiotic
    )
  
  # WEIGHTED AVERAGE --------------------------------------------------------
  ## combine the three transects' cover fractions using angular IDW weights to get 
  ## final coverage fractions
  pixel_coords <- pixel_coords %>% 
    mutate(
      cover_cyano = w1*cyano_t1 + w2*cyano_t2 + w3*cyano_t3,
      cover_lich = w1*lich_t1 + w2*lich_t2 + w3*lich_t3,
      cover_moss = w1*moss_t1 + w2*moss_t2 + w3*moss_t3,
      cover_midveg = w1*mid_t1 + w2*mid_t2 + w3*mid_t3,
      cover_subcanopy = w1*sub_t1 + w2*sub_t2 + w3*sub_t3,
      cover_canopy = w1*can_t1 + w2*can_t2 + w3*can_t3,
      cover_abiotic = w1*abi_t1 + w2*abi_t2 + w3*abi_t3
    ) %>% 
    # mask pixels beyond 100 m
    mutate(across(starts_with("cover_"), ~ifelse(dist_m > 100, NA, .)))
  
  # BUILD OUTPUT RASTER STACK -----------------------------------------------
  # one layer per functional group, matched to footprint grid
  make_layer <- function(col_name) {
    r <- fpt_template
    values(r) <- pixel_coords[[col_name]]
    names(r) <- col_name
    return(r)
  }
  
  veg <- c(
    make_layer("cover_cyano"),
    make_layer("cover_lich"),
    make_layer("cover_moss"),
    make_layer("cover_midveg"),
    make_layer("cover_subcanopy"),
    make_layer("cover_canopy"),
    make_layer("cover_abiotic")
  )
  
  return(veg)
} # end build_veg_idw_raster

# ====================================
# Section 6: Footprint File Loading & Metadata
# ====================================
# Purpose: Chu et al. (2021) for all US-Mpj .tif files and parse temporal metadata
# from filename convention: US-Mpj_YYYY_MM_[DAY/NIGHT]_fpt_weight.tif
#
# Output: results
# - one row per fpt file
# - cols: filepath, year, month, period
# - empty cols for weighted cover fractions
# ====================================
# path to Chu et al. footprint raster files
fpt_dir <- here("monthly_footprint_climatology_weight_map")

# get all US-Mpj footprint files
fpt_files <- list.files(
  path     = fpt_dir,
  pattern  = "^US-Mpj",
  full.names = TRUE
)
cat("Found", length(fpt_files), "footprint files\n")

# FUNCTION: Extract seasonal/timing from file name ------------------------
# extract year, month, day/night from file name
# filename format: US-Mpj_YYYY_MM_DAY_fpt_weight.tif or US-Mpj_YYY_MM_NIGHT_fpt_weight.tif

parse_filename <- function(filepath) {
  fname <- basename(filepath)
  parts <- str_split(fname, "_")[[1]]
  
  # parts: [1] = US-MPj; [2] = YYYY; [3]=MM; [4]=DAY/NIGHT; [5]=fpt
  year <- as.integer(parts[2])
  
  # check if it's montly or annual
  if (parts[3] %in% c("DAY", "NIGHT")) {
    # if we wanted to use an annual file that doesn't include months, which we do have
    month <- NA
    period <- parts[3]
  } else {
    month <- as.integer(parts[3])
    period <- parts[4]
  }
  data.frame(filepath, year, month, period, stringsAsFactors=FALSE)
}

# parse all filenames into a metadata table
fpt_metadata <- bind_rows(lapply(fpt_files, parse_filename))

# sort by year -> month -> day then night
results <- fpt_metadata %>%
  mutate(period_order = ifelse(period == "DAY", 1, 2)) %>%
  arrange(year, month, period_order) %>%
  select(-period_order)

# initialize empty columns for weighted cover results
# these will be filled in the overlay loop in Section 7
results$w_cyano     <- NA
results$w_lich      <- NA
results$w_moss      <- NA
results$w_midveg    <- NA
results$w_subcanopy <- NA
results$w_canopy    <- NA
results$w_abiotic   <- NA
  
# build master template raster covering the union of all footprint extents
# this ensures all frames in the animation use the same spatial grid
template     <- rast(results$filepath[1])
all_extents  <- lapply(results$filepath, function(f) ext(rast(f)))
master_ext   <- all_extents[[1]]
for (ex in all_extents[-1]) master_ext <- terra::union(master_ext, ex)
  
master_template <- rast(
  ext        = master_ext,
  resolution = res(template),
  crs        = crs(template)
) 

# consistent color scale for animation: find max weight across all files
# so contour levels are comparable frame to frame
weight_max <- max(sapply(results$filepath, function(f) {
  global(rast(f), "max", na.rm = TRUE)[[1]]
}))
cat("Max footprint weight across all files:", weight_max, "\n")

x_min <- ext(master_template)[1]
x_max <- ext(master_template)[2]
y_min     <- ext(master_template)[3]
y_max     <- ext(master_template)[4]

cat("Plot bounds (mirrored):\n")
cat("  x:", x_min, "to", x_max, "\n")
cat("  y:", y_min, "to", y_max, "\n")

# ====================================
# Section 7: Weighted Cover Overlay Loop
# =============================================================
# Purpose: for each monthly footprint file, build a vegetation
# raster matched to that footprint's grid, multiply footprint
# weights by vegetation cover layers, and sum across all pixels
# to get a single weighted cover fraction per functional group.
#
# This produces a time series of footprint-weighted vegetation cover 
# for each month/year, which can then be related to tower flux measurements.
#
# Output: results data frame with w_cyano, w_lich, w_moss,
# w_midveg, w_subcanopy, w_canopy, w_abiotic columns
# filled in for each footprint file
# =============================================================
cat("Processing", nrow(results), "footprint files...\n")

for (i in 1:nrow(results)) {
  if (i %% 10 == 0) cat("Processing file", i, "of", nrow(results), "\n")
  
  tryCatch({
    # load this month's footprint weight raster
    fpt <- rast(results$filepath[i])
    
    # build vegetation raster matched to this footprint's grid
    veg <- build_veg_idw_raster(fpt)
    
    # normalize footprint weights so they sum to 1. this ensures weighted cover
    # fractions are on a 0-1 scale
    fpt_norm <- fpt / global(fpt, "sum", na.rm=TRUE)[[1]]
    
    # multiply each vegetation layer by footprint weights, then sum across pixels
    weighted <- veg * fpt_norm
    cover <- global(weighted, "sum", na.rm = TRUE)
    
    # store results
    results$w_cyano[i] <- cover$sum[1]
    results$w_lich[i] <- cover$sum[2]
    results$w_moss[i] <- cover$sum[3]
    results$w_midveg[i] <- cover$sum[4]
    results$w_subcanopy[i] <- cover$sum[5]
    results$w_canopy[i] <- cover$sum[6]
    results$w_abiotic[i] <- cover$sum[7]
  }, error=function(e) {
    # if a file fails, print which one and continue
    cat("Error on file:", basename(results$filepath[i]), "-", e$message, "\n")
  })
}

# check results
cat("Done. Summary of weighted cover fractions:\n")
summary(results[, c("w_cyano", "w_lich", "w_moss",
                    "w_midveg", "w_canopy", "w_abiotic")])

# check for any failed files
n_failed <- sum(is.na(results$w_cyano))
if (n_failed > 0) {
  cat("Warning:", n_failed, "files failed to process\n")
} else {
  cat("All files processed successfully\n")
}

## write results to project-relative outputs folder
dir.create(here("outputs"), showWarnings=FALSE)
write.csv(results, here("outputs", "results_coverfrac.csv"), row.names = FALSE)

# =============================================================
# END OF MODEL BUILDING
# =============================================================
# Downstream products (figures, animations, stats) need:
#   - transects, cover_fractions          (Sections 3-4)
#   - build_veg_idw_raster()              (Section 5)
#   - results, master_template, weight_max, x_min/x_max/y_min/y_max (Sections 6-7)
#   - tower_lat, tower_lon, tower_df, m_per_deg_lat, m_per_deg_lon, transect_lines (Section 2-3)
# These can be re-loaded from "results_coverfrac.csv" rather than rerun if the
# heavy overlay loop (Section 7) has already been completed once.
# =============================================================
