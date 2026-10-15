README
================
This repository hosts all the data, data analysis tools, model functions, and model framework for the ecophysiology-based model framework to optimize biocrusts in dryland carbon cycling models, funded by Middlebury College. <br>
***Data Folders***  
1. **LICOR** This folder includes all of the raw and standardized data
from infrared gas analysis lab experiments. This data is further
organized by crust type, including moss cultures, combined CO2/water
curves on cultures, combined CO2/water curves on cyanobacteria crust,
cyanobacteria crust single factor curves, and lichen crust curves.
Within a given folder, there are raw data files outputted directly from
the LICOR-6800. Defined below are the relevant columns:  
- elapsed (seconds)  
- A_scaled (mumol/m2/s): carbon assimilation scaled by dry weight and
mass of sample  
- A (mumol/kg/s): carbon assimilation recorded by LICOR machine  
- Qin (mumol/m2/s): light level observed by LICOR machine  
- CO2_r (ppm): CO2 level observed by LICOR machine  
- Tair (degC): temperature observed by LICOR machine  

These files are used as a foundation to develop relevant mathematical
relationships. This data has been cleaned and scaled appropriately by
dry weights and masses. These data were collected by Asha Kannan
(2025-2026), and cleaned with protocol developed by Kirsten Coe.  

2.  **NMEG** This folder includes all the flux tower data files, with
    both gaps and gapfilled, available. Our model pulls the gapfilled
    files. All files are named with the following format:
    US-Mpj\_\[year\]\_gapfilled.txt. Defined below are the relevant
    columns:  

- TIMESTAMP_START (YYYYMMDDHHMM) and TIMESTAMP_END (YYYYMMDDHHMM): used
  as timesteps in the model  
- PPFD_IN (mumol/m2/s): light used in the model  
- GPP_DT (mumol/m2/s): ecosystem GPP (gapfilled) in model  
- RECO_DT (mumol/m2/s): ecosystem respiration (gapfilled) in model  
- TA_F (degC): temperature recorded by flux tower  
- P_F_FLAG (–)  

These files are used in DrylandBSCModel_June2026.R in additional to
angular_interpolation_building.R. These data were provided by Rae DeVan
& Marcy Litvak at the University of New Mexico.  

3.  **monthly_footprint_climatology_weight_map**  
    This folder includes all the monthly footprint climatology .tif
    files downloaded from Chu et al. (2021) Supplementary Dataset. See
    more information here: <https://zenodo.org/records/4015350>. All
    files are in the following format:
    US-Mpj\_\[YYYY\]*\[MM\]*\[DAY/NIGHT\]\_fpt_weight.tif.  

These files are used in angular_interpolation_building.R.  

4.  **PJ_soil_allyears_30min.csv**<br> This file contains relevant soil
    moisture and soil temperature information from sensors at the US-Mpj
    site. Columns are defined in this format for temperature sensors:
    SOILT\_\[P,O,J\]\[1,2,3\]\_\[5,10,30\]\_AVG, where P/O/J relates to
    the relevant canopy cover (pinyon pine, open, or juniper), 1/2/3
    refers to the unit number, and 5/10/30 refers to the depth of the
    relevant sensor, either 5, 10, or 30 cm below the surface. <br>

This file was used as environmental inputs in
DrylandBSCModel_June2026.R. These data were provided by Rae DeVan &
Marcy Litvak at the University of New Mexico. <br>

5.  **mpj_transects.xslx** <br> This file contains vegetation
    identification as measured through the LPI methodology in the field
    by Kirsten Coe and Asha Kannan in October 2025. These data are used
    to build the vegetation raster in angular_interpolation_building.R.
    <br>

6.  **results_coverfrac.csv**<br> This file contains the results of the
    cover fractions of different vegetation types that is an output of
    angular_interpolation_building.R. Defined below are the relevant
    columns: <br>

- year, month, period (day or night), w_cyano (weight of cyanobacteria
  crust), w_lich (weight of lichen crust), w_moss (weight of moss
  crust), w_midveg (weight of mid-level vegetation), w_subcanopy (weight
  of pinyon/juniper but not as a canopy), w_canopy (weight of
  overhanging canopy), w_abiotic (weight of nonphotosynthetic material,
  eg hardware, fallen leaves, etc)  

***Data Analysis***  
1. **RollingSlope.R**  
Authors: Casey Schine and Lauren Marsala  
Purpose: time-series method that identifies points where the calculated
slope changes too drastically compared to other nearby data points. Used
to clean outliers in LICOR data. Arguments: x, y, width  

2.  **PhotosynthesisFunctions.R**  
    Authors: Casey Schine and Asha Kannan, adapted from Webb et
    al. (1974); Farquhar et al. (2001); and Ye et al. (2021)  
    Purpose: Functions for temperature, light, and water developed to
    describe biocrust photosynthesis response  

Temperature: calc_Pg_fxn_T_Tk_rowwise  
Arguments: temp (measured soil temperature); T_k (controls rate of
rise); Pgmax (max photosynthesis)  
Output: out_P  
Temperature: f_H  
Arguments: Tv (measured temperature); delta_S (entropy); delta_Hd
(deactivation enthalpy); R_gas (universal gas constant)  
Output: f_H (high-temperature deactivation term)  
Dark respiration: calc_RD_fxn_T_rowwise  
Arguments: temp (measured temperature); r_r (temperature sensitivity
coefficient); R_o (base respiration at 0C)  
Output: out_R  

Light: calc_PN_fxn_I_rowwise  
Arguments: PAR (measured light level); f_Io_Icomp (apparent quantum
yield); b (photoinhibition coefficient); g (curvature coefficient);
Icomp (light compensation point)  
Output: PN  

Water: low_water_response_fxn  
Arguments: WC (normalized water content); WC_k (threshold water
content)  
Output: low-water activation (out_P)  
Water: high_water_response_fxn  
Arguments: WC (normalized water content); toggled intercepts for
high-water activation  
Output: high-water activation (out_P)  

3.  **ye_et_al_light_fit.xlsx**  
    This solver spreadsheet, developed by Casey Schine, was adapted from
    Ye et al. (2007) for fitting light curves.  
    Inputs: Bin LICOR data for each light level, then average per bin
    such that there is one value for carbon assimilation per light
    level. Instructions are in Sheet 2 of the spreadsheet.  
    Outputs: Save all parameters that are outputted by the spreadsheet,
    including the sum of squares.  

4.  **angular_interpolation_building.R**  
    Script developed by Asha Kannan, to build the vegetation raster
    based on field-collected spatial vegetation cover data.  
    Packages required: “readxl”, “dplyr”, “tidyr”, “ggplot2”, “terra”,
    “stringr”, “av”, “here”  
    Inputs required: mpj_transects.xlsx;
    monthly_footprint_climatology_weight_map  
    Output: results_coverfrac.csv  

5.  **DrylandBSCModel_June2026.R**  
    Script developed by Asha Kannan, as ecophysiology-based model
    framework.  
    Packages required: “tidyverse”, “readxl”, “minpack.lm”, “lubridate”,
    “zoo”, “patchwork”  
    Inputs required: NMEG, PJ_soil_allyears_30min.csv,
    results_coverfrac.csv  
    Outputs: model_daily, model_monthly, model_annual  
