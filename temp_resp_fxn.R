# temperature model function
# Casey Schine & Asha Kannan (March 2026)
library(tidyverse)

calc_Pg_fxn_T_Tk_rowwise <- function(temp, T_k, Pgmax){
  out_P <- Pgmax*(1-exp(-temp/T_k))
  return(out_P)
}

f_H <- function(Tv, delta_S, delta_Hd, R_gas = 8.314) {
  (1 + exp((298.15 * delta_S - delta_Hd) / (298.15 * R_gas))) /
    (1 + exp((delta_S * (Tv+273.15) - delta_Hd) / (R_gas * (Tv+273.15))))
}

#total temp response is equal to A*calc_Pg_fxn_T_Tk_rowwise*f_H

temp <- c(seq(1,40,1))

temp_resp_tbl <- tibble(temp)

temp_resp_tbl <- temp_resp_tbl |>
  mutate(low_temp_resp_Tk5 =calc_Pg_fxn_T_Tk_rowwise(temp=temp, T_k=5, Pgmax=1), ## initially Pg
         high_temp_resp=f_H(Tv=temp,delta_S=647, delta_Hd=200000),
         both_resp=low_temp_resp_Tk5*high_temp_resp)

# deltaH = 149250
# deltaS = 485


ggplot(temp_resp_tbl) +
  geom_line(aes(x=temp, y=both_resp))
