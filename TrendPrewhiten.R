# Prewhitening to eliminate the influence of serial correlation on the Sen's slope trend estimator and Mann-Kendall test 
# Reference
# Zhang X, Vincent L A, Hogg W D, et al. Temperature and precipitation trends in Canada during the 20th century[J]. Atmosphere-ocean, 2000, 38(3): 395-429.
# Wang X L, Swail V A. Changes of Extreme Wave Heights in Northern Hemisphere Oceans and Related Atmospheric Circulation Regimes[J]. Journal of Climate, 2001, 14(10): 2204-2221.

rm(list = ls())
library(wql)  # MannKen function, which is a Sen's Slope Estimator and M-K test function
TrendPrewhiten <- function(Y) {
  n <- length(Y)
  c0 <- stats::acf(Y, plot = FALSE, na.action = na.pass)$acf[2]
  if (c0 <= 0.05) {
    b <- mannKen(Y)$sen.slope  # the final slope
    p <- mannKen(Y)$p.value    # the final p-value
    break
  } else {
    W <- rep(NA, n)
    for(i in 2:n) {
      W[i-1] <- (Y[i] - c0 * Y[i-1]) / (1 - c0)
    }
    b0 <- mannKen(W)$sen.slope
    p0 <- mannKen(W)$p.value
    W  <- Y - b0 * c(1:n)
    c1 <- stats::acf(W, plot = FALSE, na.action = na.pass)$acf[2]
    if(c1 < 0.05) {
      b <- b0  # the final slope
      p <- p0  # the final p-value
      break
    } else {
      W <- rep(NA, n)
      for(i in 2:n) {
        W[i-1] <- (Y[i] - c1 * Y[i-1]) / (1 - c1)
      }
      b1 <- mannKen(W)$sen.slope
      p1 <- mannKen(W)$p.value
      if (abs(c1 - c0) < 1e-4 & abs(b1 - b0) < 1e-4) {
        b <- b1
        p <- p1
        break
      } else {
        while (abs(c1 - c0) > 1e-4 | abs(b1 - b0) > 1e-4) {
          c0 <- c1
          b0 <- b1
          W  <- Y - b0 * c(1:n)
          c1 <- stats::acf(W, plot = FALSE, na.action = na.pass)$acf[2]
          if (c1 < 0.05) {
            b <- b1
            p <- p1
            break
          } else {
            W <- rep(NA, n)
            for(i in 2:n) {
              W[i-1] <- (Y[i] - c1 * Y[i-1]) / (1 - c1)
            }
            b1 <- mannKen(W)$sen.slope
            p1 <- mannKen(W)$p.value
            if (abs(c1 - c0) < 1e-4 & abs(b1 - b0) < 1e-4) {
              b <- b1
              p <- p1
              break
            }
          }
        }
      }
    }
  }
  result <- list(slope = b, p.value = p)
  return(result)
}
