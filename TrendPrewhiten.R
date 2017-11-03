# Prewhitening to eliminate the influence of serial correlation on the Sen's 
# slope trend estimator and Mann-Kendall test.

# Reference:
# Zhang X, Vincent L A, Hogg W D, et al. Temperature and precipitation trends in Canada during the 20th century[J]. Atmosphere-ocean, 2000, 38(3): 395-429.
# Wang X L, Swail V A. Changes of Extreme Wave Heights in Northern Hemisphere Oceans and Related Atmospheric Circulation Regimes[J]. Journal of Climate, 2001, 14(10): 2204-2221.

# Note:
# The function "mannKen" of wql package is to estimate slope and to do mk test.
# the stats:acf and TSA:acf has a small difference, some people may have already
# loaded the TSA package, so in order to avoid confusion, it is need to use "::"
# noted that the function comes from which package.
# The precedure usually converges within 1 to 19 iterations. There may be a few of 
# time series that do not converge, so setting a maximum number of iterations 100

#-------------------------------------------------------------------------------
# TrendPrewhiten begin
TrendPrewhiten <- function(Y) {
  # step 1
  n <- length(Y)
  r0 <- stats::acf(Y, plot = FALSE, na.action = na.pass)$acf[2]
  k <- 1  # to calculate the iteration times
  if (r0 <= 0.05) {
    r <- r0  # the final autocorrelation coefficient
    b <- wql::mannKen(Y)$sen.slope  # the final slope
    p <- wql::mannKen(Y)$p.value    # the final p-value
  } else {
    # prewhitening
    W <- rep(NA, n-1)
    W <- (Y[1:(n-1)] - r0 * Y[2:n]) / (1 - r0)
    # estimating slope and p.vlue
    b0 <- wql::mannKen(W)$sen.slope  
    p0 <- wql::mannKen(W)$p.value    
    
    # below is step 2 & 3
    repeat {
      # detrending
      W  <- Y - b0 * c(1:n)
      r1 <- stats::acf(W, plot = FALSE, na.action = na.pass)$acf[2]
      if (r1 <= 0.05) {
        r <- r0  # the final autocorrelation coefficient
        b <- b0  # the final slope
        p <- p0  # the final p-value
        break
      } else {
        # re-prewhitening
        W <- rep(NA, n-1)
        W <- (Y[1:(n-1)] - r1 * Y[2:n]) / (1 - r1)
        # re-estimating slope and p.vlue
        b1 <- wql::mannKen(W)$sen.slope  
        p1 <- wql::mannKen(W)$p.value
        if (abs(r1 - r0) < 1e-4 & abs(b1 - b0) < 1e-4) {
          r <- r1
          b <- b1
          p <- p1
          break
        } else {
          r0 <- r1
          b0 <- b1
          p0 <- p1
        }
      }
      k <- k + 1  # calculating the iteration times
      if (k > 100) {
        r <- r1
        b <- b1
        p <- p1
        break # if the iterations more than 100, then stop the iteration
      }
    }
  }
  return(list(icor = r, slope = b, p.value = p, iteration = k))
}
# end of function TrendPrewhiten------------------------------------------------
#-------------------------------------------------------------------------------
