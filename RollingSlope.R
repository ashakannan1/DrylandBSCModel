## Rolling slope
## developed by: Casey Schine, Lauren Marsala
rolling_slope <- function(x, y, width = 7) {
  rollapplyr(1:length(y), width, function(i) {
    coef(lm(y[i] ~ x[i]))[2]
  }, fill = NA)
}

