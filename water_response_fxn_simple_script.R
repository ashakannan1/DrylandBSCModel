low_water_response_fxn <- function(WC, WC_k, Pgmax){
  out_P <- Pgmax*(1-exp(-WC/WC_k))
  return(out_P)
}

high_water_response_fxn <- function(WC, slope=-0.4, 
                                    intercept= 0.9, Pgmax){
  out_P <- Pgmax * ((slope*WC)+intercept)
  return(out_P)
}


WC<- c(seq(0,1,0.01))

WC_resp_tbl <- tibble(WC) 

WC_resp_tbl <- WC_resp_tbl|>
  mutate(schine_fit_WCk_low = low_water_response_fxn(WC=WC, WC_k=0.2, Pgmax=1),
         schine_fit_WCk_med = low_water_response_fxn(WC=WC, WC_k=0.3, Pgmax=1),
         schine_fit_WCk_high = low_water_response_fxn(WC=WC, WC_k=0.4, Pgmax=1))


WC_resp_tbl <- WC_resp_tbl |>
  mutate(high_water_resp_highslope = high_water_response_fxn(
    WC=WC, slope=-0.6, Pgmax=1),
    combo_resp = schine_fit_WCk_high*high_water_resp_highslope,
    combo_resp_scaled=combo_resp/max(combo_resp))

ggplot() +
  geom_line(data = WC_resp_tbl,
            aes(x = WC, y = schine_fit_WCk_low), color="darkblue") +
  geom_line(data = WC_resp_tbl,
            aes(x = WC, y = schine_fit_WCk_med), color="blue") +
  geom_line(data = WC_resp_tbl,
            aes(x = WC, y = schine_fit_WCk_high), color="skyblue") +
  geom_line(data = WC_resp_tbl,
            aes(x = WC, y = high_water_resp_highslope), color="red") +
  geom_line(data = WC_resp_tbl,
            aes(x = WC, y = combo_resp_scaled), color="black") +
  theme_tufte()