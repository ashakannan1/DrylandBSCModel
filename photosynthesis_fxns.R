## Functions for calculating photosynthesis as a function of temperature, light, water
## Written by Casey Schine, modified by Asha Kannan
## Date: 2025-2026


# Temperature/Dark Respiration --------------------------------------------
## Gross photosynthesis as a function of temperature
calc_Pg_fxn_T_Tk_rowwise <- function(temp, T_k, Pgmax){
  out_P <- Pgmax*(1-exp(-temp/T_k))
  return(out_P)
}

## High temperature deactiviation term
f_H <- function(Tv, delta_S, delta_Hd, R_gas = 8.314) {
  (1 + exp((298.15 * delta_S - delta_Hd) / (298.15 * R_gas))) /
    (1 + exp((delta_S * (Tv+273.15) - delta_Hd) / (R_gas * (Tv+273.15))))
}

## Dark respiration as a function of temperature
calc_RD_fxn_T_rowwise <- function(temp, r_r, R_o){
  if (temp<=0) {
    out_R=0
  } else {
    out_R <- R_o*exp(r_r*temp)
  }
  return(out_R)
}


# Light -------------------------------------------------------------------
## Net photosynthesis as a function of light
## adapted from Ye (2007)
calc_PN_fxn_I_rowwise <- function(PAR, f_Io_Icomp, b, g, Icomp) {
  PN <- f_Io_Icomp * ((1 - b * PAR) / (1 + g * PAR)) * (PAR - Icomp)
  return(PN)
}

## Gross light response: 0 in the dark, rises linearly to RD_T at Icomp
## Pn + RD_T at Icomp
calc_Pg_I <- function(PAR, Pn_I, RD_T, Icomp) {
  PAR <- pmax(PAR, 0)
  ifelse(PAR < Icomp, RD_T * PAR / Icomp, pmax(0, Pn_I + RD_T))
}

# Normalized gross photosynthesis response to light
## Converts net light response to gross photosynthesis used a fixed respiration rate,
## then normalizes the gross curve so that its maximum is 1. Returns a dimensionless light 
## scalar (0-1)
make_light_fraction <- function(f_Io_Icomp, b, g, Icomp, RD_ref, PAR_max=3000) {
  ## build the reference curve
  PAR_seq <- seq(0, PAR_max, by=1)
  Pn_seq <- calc_PN_fxn_I_rowwise(PAR = PAR_seq, f_Io_Icomp = f_Io_Icomp, b = b, g = g, Icomp = Icomp)
  Pg_max <- max(calc_Pg_I(PAR = PAR_seq, Pn_I = Pn_seq, RD_T = RD_ref, Icomp = Icomp), na.rm = TRUE)
  
  # return a function that only needs PAR
  function(PAR) {
    Pn <- calc_PN_fxn_I_rowwise(PAR = PAR, f_Io_Icomp = f_Io_Icomp, b = b, g = g, Icomp = Icomp)
    Pg <- calc_Pg_I(PAR=PAR, Pn_I = Pn, RD_T = RD_ref, Icomp=Icomp)
    pmin(1, pmax(0, Pg / Pg_max))
  }
}

# Water (moss) ------------------------------------------------------------
## Gross photosynthesis as a function of water
## high and low water response (used for moss only)
low_water_response_fxn <- function(WC, WC_k, Pgmax){
  out_P <- Pgmax*(1-exp(-WC/WC_k))
  return(out_P)
}

high_water_response_fxn <- function(WC, slope=-0.4, 
                                    intercept= 0.9, Pgmax){
  out_P <- Pgmax * ((slope*WC)+intercept)
  return(out_P)
}
