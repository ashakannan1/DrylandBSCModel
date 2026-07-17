## Functions for calculating photosynthesis as a function of temperature, light, water
## Written by Casey Schine, modified by Asha Kannan
## Date: 2025-2026


# Temperature/Dark Respiration --------------------------------------------
## Gross photosynthesis as a function of temperature
calc_Pg_fxn_T_Tk_rowwise <- function(temp, T_k, Pgmax){
  out_P <- Pgmax*(1-exp(-temp/T_k))
  return(out_P)
}

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
## adapted from Ye (2021)
calc_PN_fxn_I_rowwise <- function(PAR, f_Io_Icomp, b, g, Icomp) {
  PN <- f_Io_Icomp * ((1 - b * PAR) / (1 + g * PAR)) * (PAR - Icomp)
  return(pmax(PN, 0))
}


# Water (moss) ------------------------------------------------------------
## Gross photosynthesis as a function of water
## high and low water response
low_water_response_fxn <- function(WC, WC_k, Pgmax){
  out_P <- Pgmax*(1-exp(-WC/WC_k))
  return(out_P)
}

high_water_response_fxn <- function(WC, slope=-0.4, 
                                    intercept= 0.9, Pgmax){
  out_P <- Pgmax * ((slope*WC)+intercept)
  return(out_P)
}



