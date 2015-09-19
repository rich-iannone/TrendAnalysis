#' Determine trend information using Sen's methodology and Kendall tests
#' @description Using a vector of data, The Thiel-Sen method for trend determination will be carried out with a Kendall test for significance. 
#' @param y a vector of evenly-spaced data points to be taken for analysis.
#' @return a data frame.
#' @import Kendall
#' @export get_trend

get_trend <- function(y){
  
  # Define the 'slope_diff' function
  slope_diff <- function(i, xx, yy, n){
    (yy[1:(n - i)] - yy[(i + 1):n])/(xx[1:(n - i)] - xx[(i + 1):n])
  }
  
  # Define the 'get_ci' function
  get_ci <- function(object, parm, level = 0.95, ...){
    
    res <- rep(0,4)
    
    dim(res) <- c(2,2)
    
    slopes <- sort(object$slopes)
    
    intercepts <- sort(object$intercepts)
    
    uquantile <- 1 - (1 - level)/2
    
    zstat <- qnorm(uquantile)
    
    k <- Kendall(object$x, object$y)
    
    c_alpha <- sqrt(k$varS) * zstat
    
    n_prime <- length(slopes)
    
    isect_sd <- sqrt(var(intercepts))
    
    isect_mean <- mean(intercepts)
    
    idx <- c(round((n_prime - c_alpha)/2), round((n_prime + c_alpha)/2))
    
    rownames(res) <- names(object$coefficients)
    
    colnames(res) <- as.character(c((1 - level)/2, 1 - (1 - level)/2))
    
    res[2,] <- slopes[idx]
    
    res[1,] <- isect_mean + c(-isect_sd * zstat, isect_sd * zstat)
    
    return(res)
  }
  
  # Create a data frame for which the data will reside
  trend_df <- as.data.frame(mat.or.vec(nr = 1, nc = 11))
  
  colnames(trend_df) <- c("trend", "signif", "lbound", "ubound",
                          "y_int", "trend_p", "tau", "autocor",
                          "valid_frac", "linest_slope", "linest_y_int")
  
  n <- length(y)
  
  x = 1:length(y)
  
  t <- x
  
  t_prime <- t[1:(n - 1)]
  
  dmap <- which(!is.na(y))
  
  ynm <- as.numeric(y[dmap])
  
  tnm <- as.numeric(t[dmap])
  
  if (length(dmap) <= 3 | length(which(ynm != 0)) < 3 | length(dmap)/n < 0.1) {
    trend_df[1,] <- NA
    return(trend_df)
  }
  
  slopes <- unlist(lapply(1:(n - 1), slope_diff, t, y, n))
  
  sni <- which(is.finite(slopes))
  
  median_slope <- median(slopes[sni])
  
  intercepts <- y - median_slope * t
  
  median_intercept <- median(intercepts)
  
  res <- list(coefficients = c(median_intercept, median_slope),
              slopes = slopes, 
              intercepts = intercepts,
              rank = 2,
              residuals = (y - median_slope * t + median_intercept),
              x = t,
              y = y)

  class(res) = c("zyp", "lm")
  
  xt_prime <- y[1:n] - median_slope * t
  
  ac <- acf(xt_prime, lag.max = 1, plot = FALSE, na.action = na.pass)$acf[2]
  
  if (is.na(ac)){
    trend_df[1,] <- NA
    return(trend_df)
  }
  
  yt_prime <- (xt_prime[2:n] - ac * xt_prime[1:(n - 1)])/(1 - ac)
  
  yt <- yt_prime[1:(n - 1)] + median_slope * t_prime
  
  dmap_prime <- which(!is.na(yt))
  
  ytnm <- as.numeric(yt[dmap_prime])
  
  Kend <- Kendall(t_prime[dmap_prime], ytnm)
  
  tau <- Kend[1]
  
  Bsig <- Kend[2]
  
  ci <- get_ci(res)
  
  trend_df[1,1] <- as.numeric(median_slope)
  trend_df[1,2] <- as.numeric(Bsig)
  trend_df[1,3] <- as.numeric(ci[2, 1])
  trend_df[1,4] <- as.numeric(ci[2, 2])
  trend_df[1,5] <- as.numeric(res$coefficients[1])
  trend_df[1,6] <- as.numeric(median_slope) * n
  trend_df[1,7] <- as.numeric(tau)
  trend_df[1,8] <- as.numeric(ac)
  trend_df[1,9] <- as.numeric(length(dmap)/length(y))
  trend_df[1,10] <- as.numeric(lm(y ~ t)$coefficients[2])
  trend_df[1,11] <- as.numeric(lm(y ~ t)$coefficients[1])
  
  return(trend_df)
}
