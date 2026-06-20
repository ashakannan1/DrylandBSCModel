 # functions for calculating P_N as a function of light, temperature, and WC

calc_Pg_fxn_T_Tk_rowwise <- function(temp, T_k, Pgmax){
  if (temp<=0) {
    out_P=0
  } else {
    if (T_k<4){ T_k=4 }
    out_P <- Pgmax*(1-exp(-temp/T_k))
  }
  return(out_P)
}

calc_RD_fxn_T_rowwise <- function(temp, r_r, R_o){
  if (temp<=0) {
    out_R=0
  } else {
    out_R <- R_o*exp(r_r*temp)
  }
  return(out_R)
}


# PN = {Pgmax[1 - exp(-f(Io) x I/Pgmax)]} - RD from Webb
calc_PN_fxn_I_rowwise <- function(PAR, P_gmax_t, R_D_t, I_k) {
  PN = P_gmax_t*(1-exp(-PAR/I_k))-R_D_t
  return(PN)
}


# PN = {Pgmax[1 - exp(-f(Io) x I/Pgmax)]} - RD from Webb
calc_Pg_fxn_I_rowwise <- function(PAR, P_gmax_t, I_k) {
  PN = P_gmax_t*(1-exp(-PAR/I_k))
  return(PN)
}


#calculate the fraction of max P or R posssible by BWC
calc_var_frac_from_BWC <- function(x, var=1, group) {
  
  if (group=="Feather") {
    var_bwc <- var * ((0.4158*(x^3))-(2.0601*(x^2))+(2.7684*x)-0.1914)
  }

  if (group=="Sphagnum") {
    var_bwc <- var*((0.8446*(x^3))-(3.7429*(x^2))+(4.9793*x)-1.0807)
  }
  
  return(var_bwc)
}

