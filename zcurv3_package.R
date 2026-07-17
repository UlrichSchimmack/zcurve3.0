

#rm(list = ls())
#options(scipen = 999)  # Disable scientific notation




### INSTALL PACKAGES (only once – manually run if needed)
if (1 == 2) {  # This block is ignored unless manually changed to (1 == 1)
  install.packages("KernSmooth")
  install.packages("parallel")
} # END install block

### LOAD LIBRARIES
library(parallel)
library(KernSmooth)



########################################################################
### DRAFT: new zcurve() package interface + main function
### NOT YET intEGRATED into zcurve.V3.90.R -- for review/testing only.
########################################################################



#######################################################################
### FITTING CONTROL
#######################################################################

zcurve_control <- function(
  boot_iter            = 0,       # bootstrap iterations (500+ recommended for final models)
  parallel             = TRUE,    # parallelize bootstrap/EM where available

  em_criterion         = 1e-6,    # EM convergence threshold
  em_max_iter          = 2000,    # max EM iterations
  em_max_iter_boot     = 500,     # max EM iterations per bootstrap replicate

  ci_alpha             = .05,     # CI level for reported intervals (default 95%)

  augment              = TRUE,    # bias correction at the lower bound (OF method)
  n_bars               = 512,     # density-estimation resolution (OF method)
  bw_est               = 0.05,    # kernel bandwidth used for fitting (lower = less smoothing)
  spline_width         = .5       # blending of truncated and normal KD (OF method)
) {
  list(
    boot_iter = boot_iter, parallel = parallel,
    em_criterion = em_criterion, em_max_iter = em_max_iter, em_max_iter_boot = em_max_iter_boot,
    ci_alpha = ci_alpha,
    augment = augment, n_bars = n_bars, bw_est = bw_est, spline_width = spline_width
  )
} ### EOF zcurve_control


#######################################################################
### MAIN FITTING FUNCTION
#######################################################################

zcurve <- function(
  ### DATA INPUT (supply exactly one consistent combination)
  zval         = NULL,
  tval         = NULL,
  df           = NULL,
  pval         = NULL,
  n            = NULL,
  yi           = NULL,
  sei          = NULL,
  cluster_id   = NULL,

  ### CORE MODEL CHOICES
  est_method   = "OF",           # "SQP", "OF", "EM", "ML_gamma"
  directional  = FALSE,          # FALSE: sign of test_stat doesn't matter (default meta-science use)
  folded       = TRUE,           # TRUE: input is |z| or |t| (matches directional = FALSE almost always)

  alpha        = 0.05,           # significance threshold
  crit         = qnorm(1 - alpha / 2),   # two-sided critical value, derived from alpha (z-scale default)
  int_beg      = crit,                   # selection window start, derived from crit
  int_end      = 6,              # right-truncation boundary for the modeling grid (test_stat > int_end: power = 1)
  edr_ci_adj   = .05,
  err_ci_adj   = .03,

  ### MODEL STRUCTURE (component grid)
  ncp          = 0:6,             # component locations
  k_ncp        = length(ncp),     # number of components
  z_sd         = rep(1, k_ncp),   # one SD per component
  z_sd_fixed   = TRUE,            # hold component SDs fixed (EM method)
  ncp_fixed    = TRUE,            # hold component locations fixed (EM method)
  w_fixed      = FALSE,           # hold component weights fixed (EM method)

  x_lim_min    = 0,               # lower bound of the modeling grid
  x_lim_max    = 6,               # upper bound of the modeling grid (plot x-axis)
  int_loc      = 1,               # bin width for local-power/local-ES reporting (0 disables)

  ### OUTPUT / DISPLAY
  version      = "V3.90",         # set to NULL or "" to suppress
  date         = "2026.07.16",    # set to NULL or "" to suppress
  show_plot    = TRUE,            # draw the plot automatically when zcurve() is called
  show_summary = TRUE,            # print the text summary automatically

  ### FITTING CONTROL
  control      = zcurve_control()
) {


### begin of actual function


##############################################
### Get Densities
##############################################

build_dens <- function(d_x, ncp, z_sd, df = NULL, 
                       curve_type = "z", folded = TRUE) {
  
  xmat  <- matrix(rep(d_x, each = length(ncp)), nrow = length(ncp))
  mumat <- matrix(rep(ncp, times = length(d_x)), nrow = length(ncp))
  
  if (curve_type == "z") {
    
    if (is.null(z_sd)) {
      z_sd <- rep(1, length(ncp))
    }
    
    if (length(z_sd) == 1) {
      z_sd <- rep(z_sd, length(ncp))
    }
    
    sdmat <- matrix(rep(z_sd, times = length(d_x)), nrow = length(ncp))
    

    if (folded) {
      
      # Standard z-curve: folded normal density.
      Dens <- dnorm(xmat, mean = mumat, sd = sdmat) +
              dnorm(-xmat, mean = mumat, sd = sdmat)
      
    } else {

      # mass of component k inside [int_beg, int_end]
      denom <- pnorm(int_end, mean = mumat, sd = sdmat) -
               pnorm(int_beg, mean = mumat, sd = sdmat)
      denom <- pmax(denom, .Machine$double.xmin)

      dens <- dnorm(xmat, mean = mumat, sd = sdmat) / denom

     
    }
    
  } else {
    
    if (is.null(df)) {
      stop("df must be supplied when curve type is not 'z'.")
    }
    
    if (length(df) == 1) {
      df <- rep(df, length(ncp))
    }
    
    dfmat <- matrix(rep(df, times = length(d_x)), nrow = length(ncp))
    
    if (folded) {
      
      # Standard z-curve: folded noncentral t density.
      dens <- dt(xmat, df = dfmat, ncp = mumat) +
              dt(-xmat, df = dfmat, ncp = mumat)
      
    } else {
      
      # Directional model: signed positive noncentral t density,
      # truncated at zero.
      denom <- pt(0, df = dfmat, ncp = mumat, lower.tail = FALSE)
      denom <- pmax(denom, .Machine$double.xmin)
      
      dens <- dt(xmat, df = dfmat, ncp = mumat) / denom
      
      # Only relevant if d_x accidentally contains negative values.
      dens[xmat < 0] <- 0
    }
  }

    
  dens

}

###

#ddd

get_densities = function(int, bw = 0.20, from = 1.96, to = 6, width = 1, augment = TRUE) {

  #from = 0; to = 6;bw = .05
  #
  #summary(int)

  n = length(seq(from,to,.02))
  grid <- seq(from, to, length.out = n)
  dens <- numeric(n)

  if (augment) {

    bnd.width <- width

    splice <- from + bnd.width
    bnd.idx <- which(grid < splice)
    int.idx <- which(grid >= splice)
    
    ## ---- Truncated normal kernel for boundary zone ----
    int.s <- int - from
    for (i in bnd.idx) {
      z.s <- grid[i] - from
      raw <- dnorm(int.s, z.s, bw)
      trunc.corr <- pnorm(z.s / bw)
      dens[i] <- mean(raw) / trunc.corr
    }
 
    ## ---- bkde for interior ----
    bkde.out <- bkde(int, bandwidth = bw, range.x = c(from, to), gridsize = n)
    dens[int.idx] <- bkde.out$y[int.idx]
    
    ## ---- Calibrate: scale bkde to match truncated normal at splice ----
    tn.at.splice <- dens[max(bnd.idx)]
    bkde.at.splice <- bkde.out$y[min(int.idx)]
    if (bkde.at.splice > 0) {
      scale.factor <- tn.at.splice / bkde.at.splice
      dens[int.idx] <- dens[int.idx] * scale.factor
    
    
      ## ---- Blend over transition ----
      blend.n <- min(5, length(bnd.idx))
      for (j in 1:blend.n) {
        idx <- max(bnd.idx) - blend.n + j
        w <- j / (blend.n + 1)
        dens[idx] <- (1 - w) * dens[idx] + w * bkde.out$y[idx] * scale.factor
      }
    
      dens <- pmax(0, dens)
    }

    
  } else {

    bkde.out <- bkde(int, bandwidth = bw, range.x = c(from, to), gridsize = n)
    dens <- bkde.out$y

  }
  
  D <- data.frame(ZX = grid, ZY = dens)
  D_bw <- D$ZX[2] - D$ZX[1]
  D$ZY <- D$ZY / sum(D$ZY * D_bw)

  #plot(D$ZX,D$ZY)

  return(D)

} # EOF Densities 

#######################################################
### End of Get Densities
#######################################################




#####################################################################
### Compute Power Function (New: Discrete & Continuous
#####################################################################


get_estimates = function(
    est_inp, test_stat, yi, sei, int_beg, int_end,
    x_lim_min, x_lim_max) {


  #-------------------------------------------------------------------
  # helper: estimate a sampling-error function sei(z) over the z grid
  #-------------------------------------------------------------------

  estimate_sei_function <- function(
      test_stat, sei, grid, int_beg, int_end, min_k = 10) {

    dat <- data.frame(test_stat = abs(test_stat), sei = sei)
    dat <- dat[is.finite(dat$test_stat) & is.finite(dat$sei) & dat$sei > 0, ]

    dat_full <- dat[dat$test_stat <= int_end, ]

    if (nrow(dat_full) < min_k) {
      return(data.frame(grid = grid, sei_grid = rep(NA_real_, length(grid))))
    }

    # full-range SE function (fallback / non-selecting case)
    fit_full <- lm(sei ~ test_stat, data = dat_full)
    b1       <- coef(fit_full)[["test_stat"]]
    if (b1 > 0) b1 <- 0
    sei_out  <- coef(fit_full)[["(Intercept)"]] + b1 * grid

    # selecting case: fit on the significant interval and extrapolate
    if (int_beg > 0) {
      dat_int <- dat[dat$test_stat >= int_beg & dat$test_stat <= int_end, ]
    if (nrow(dat_int) >= min_k) {
      fit_int <- lm(sei ~ test_stat, data = dat_int)
      b1      <- coef(fit_int)[["test_stat"]]
      if (b1 > 0) b1 <- 0
      sei_out <- coef(fit_int)[["(Intercept)"]] + b1 * grid
      }
    }

    sei_out <- pmax(sei_out, min(dat$sei))   # floor at smallest observed SE
    data.frame(grid = grid, sei_grid = sei_out)
  }



  #-------------------------------------------------------------------
  # helper: empirical-Bayes effect-size estimates
  #-------------------------------------------------------------------

  get_es_estimates <- function(
      grid,
      pow_grid,
      sei_grid,
      w_grid,
      ncp_ext,
      w_all_ext,
      test_stat,
      yi,
      directional,
      folded,
      crit,
      int_beg,
      int_end,
      x_lim_min,
      x_lim_max,
      int_loc = 0.2
     ) {

    xz <- grid[2] - grid[1]

    ok       <- is.finite(grid) & is.finite(w_grid) & is.finite(sei_grid)
    stopifnot(length(grid) == length(w_grid), length(grid) == length(sei_grid))
    grid     <- grid[ok]
    w_grid   <- w_grid[ok]
    sei_grid <- sei_grid[ok]

    # normalize reconstructed density to probability weights
    w_grid <- w_grid / sum(w_grid, na.rm = TRUE)

    # ---- Empirical-Bayes posterior-mean ncp given observed z ---------------
    #   ncp_hat(z) = sum_k w_k phi(z - mu_k) mu_k / sum_k w_k phi(z - mu_k)
    if (!folded) {
      Lik <- outer(grid, ncp_ext, function(x, mu) dnorm(x - mu))
    } else {
      Lik <- outer(grid, ncp_ext, function(x, mu) dnorm(x - mu) + dnorm(x + mu))
    }

    num <- as.numeric(Lik %*% (w_all_ext * ncp_ext))
    den <- as.numeric(Lik %*% w_all_ext)

    ncp_grid <- num / den
    ncp_grid[!is.finite(ncp_grid)] <- NA_real_

    if (directional) {
      es_grid <- ncp_grid * sei_grid
    } else {
      es_grid <- abs(ncp_grid) * sei_grid
    }

    # ---- effect-size heterogeneity (weighted SD of the shrunken ES) --------
    m      <- sum(w_grid * es_grid)
    es_tau <- sqrt(sum(w_grid * (es_grid - m)^2))

    # ---- Extreme-value (tail) correction -----------------------------------
    if (directional) {
      es_ext_neg <- mean(yi[test_stat < -int_end])
    } else {
      es_ext_neg <- abs(mean(yi[test_stat < -int_end]))
    }
    w_ext_neg <- mean(test_stat < -int_end)

    es_ext_pos <- mean(yi[test_stat > int_end])
    w_ext_pos  <- mean(test_stat > int_end)

    grid_share <- 1 - w_ext_neg - w_ext_pos

    grid_ext    <- grid
    es_grid_ext <- es_grid
    w_grid_ext  <- w_grid * grid_share
    if (directional) sig_ext <- grid > crit else sig_ext <- abs(grid) > crit

    if (w_ext_neg > 0) {
      grid_ext    <- c(min(grid) - xz, grid_ext)
      es_grid_ext <- c(es_ext_neg,      es_grid_ext)
      w_grid_ext  <- c(w_ext_neg,       w_grid_ext)
      sig_ext     <- c(TRUE,            sig_ext)   # extreme tail is always significant
    }
    if (w_ext_pos > 0) {
      grid_ext    <- c(grid_ext,    max(grid) + xz)
      es_grid_ext <- c(es_grid_ext, es_ext_pos)
      w_grid_ext  <- c(w_grid_ext,  w_ext_pos)
      sig_ext     <- c(sig_ext,     TRUE)
    }

    # ---- weighted MEAN over all reconstructed results ----------------------
    es_mean_all <- sum(es_grid_ext * w_grid_ext, na.rm = TRUE)

    # ---- weighted mean over reconstructed SIGNIFICANT results only ---------
    if (any(sig_ext) && sum(w_grid_ext[sig_ext], na.rm = TRUE) > 0) {
      w_sig       <- w_grid_ext[sig_ext] / sum(w_grid_ext[sig_ext], na.rm = TRUE)
      es_mean_sig <- sum(es_grid_ext[sig_ext] * w_sig, na.rm = TRUE)
    } else {
      es_mean_sig <- NA_real_
    }

    # ---- weighted MEDIAN over all reconstructed results --------------------
    o      <- order(es_grid_ext)
    es_srt <- es_grid_ext[o]
    w_srt  <- w_grid_ext[o]
    keep   <- is.finite(es_srt) & is.finite(w_srt) & w_srt > 0
    es_srt <- es_srt[keep]
    w_srt  <- w_srt[keep]

    if (length(es_srt) >= 2 && sum(w_srt) > 0) {
      cum_w     <- cumsum(w_srt) / sum(w_srt)
      es_median <- approx(cum_w, es_srt, xout = 0.5, rule = 2,
                          ties = list("ordered", mean))$y
    } else {
      es_median <- NA_real_
    }

    # ---- local effect sizes by z-bin ----------------------------------------
    z_breaks <- seq(x_lim_min, x_lim_max, by = int_loc)
    if (tail(z_breaks, 1) < int_end) z_breaks <- c(z_breaks, int_end)

    local_es <- matrix(NA_real_, nrow = 5, ncol = length(z_breaks) - 1)
    rownames(local_es) <- c("z_min", "z_max", "weight", "mean_z", "es")
    colnames(local_es) <- paste0(head(z_breaks, -1), "-", tail(z_breaks, -1))

    for (j in seq_len(length(z_breaks) - 1)) {
      in_bin <- grid >= z_breaks[j] & grid < z_breaks[j + 1]
      if (j == length(z_breaks) - 1) {
        in_bin <- grid >= z_breaks[j] & grid <= z_breaks[j + 1]
      }
      if (!any(in_bin)) next

      bin_weight <- sum(w_grid[in_bin], na.rm = TRUE)
      if (!is.finite(bin_weight) || bin_weight <= 0) next

      w_bin <- w_grid[in_bin] / bin_weight

      local_es["z_min",  j] <- z_breaks[j]
      local_es["z_max",  j] <- z_breaks[j + 1]
      local_es["mean_z", j] <- sum(grid[in_bin] * w_bin, na.rm = TRUE)
      local_es["weight", j] <- bin_weight
      local_es["es",     j] <- sum(es_grid[in_bin] * w_bin, na.rm = TRUE)
    }

    list(
      es_mean_all   = es_mean_all,
      es_median_all = es_median,
      es_mean_sig   = es_mean_sig,
      es_tau        = es_tau,
      local_es      = local_es[5, ],
      local_w       = local_es[3, ]
    )

  } ### EOF get_es_estimates


  ################################
  ### Power Computation
  ################################

  compute_power <- function(est_inp, int_beg, int_end, grid,
        test_stat, crit, directional, alpha,
        x_lim_min, x_lim_max, int_loc, folded,
        z_sd = NULL, df = NULL) {

    w_inp <- est_inp$w_inp
    ncp   <- est_inp$ncp
    z_sd  <- est_inp$z_sd

    xz      <- grid[2] - grid[1]
    k_tests <- length(test_stat)

    # extreme-value input PROPORTIONS (of total)
    ext_neg_inp <- sum(test_stat < -int_end) / k_tests
    ext_pos_inp <- sum(test_stat >  int_end) / k_tests

    # component power (two-tailed, with sign error)
    pow_neg <- pnorm(-crit, ncp)
    pow_pos <- 1 - pnorm(crit, ncp)
    pow_dir <- c(pow_neg[ncp < 0], pow_pos[ncp >= 0])
    pow     <- pow_neg + pow_pos

    # power conditional on selection interval; sel_crit floors window bound at 0
    sel_crit <- max(0, int_beg)
    if (directional)  pow_sel <- pnorm(ncp,      sel_crit) + pnorm(-sel_crit, ncp)
    if (!directional) pow_sel <- pnorm(abs(ncp), sel_crit) + pnorm(-sel_crit, abs(ncp))

    # selection correction: input -> population weights (grid only)
    w_all <- w_inp / pow_sel
    w_all <- w_all / sum(w_all)

    # ---- extend vectors with extreme components (used only for EDR/ERR) ----
    add_neg <- (ext_neg_inp > 0 & directional == FALSE)
    add_pos <- (ext_pos_inp > 0)

    # grid weight after removing whatever extreme mass is present (computed once)
    grid_share <- 1 - (if (add_neg) ext_neg_inp else 0) -
                      (if (add_pos) ext_pos_inp else 0)

    ncp_ext     <- ncp
    pow_ext     <- pow
    pow_dir_ext <- pow_dir
    pow_sel_ext <- pow_sel
    w_inp_ext   <- w_inp * grid_share

    if (add_neg) {
      ncp_ext     <- c(-(int_end + 1), ncp_ext)
      pow_ext     <- c(1, pow_ext)
      pow_dir_ext <- c(1, pow_dir_ext)
      pow_sel_ext <- c(1, pow_sel_ext)
      w_inp_ext   <- c(ext_neg_inp, w_inp_ext)
    }

    if (add_pos) {
      ncp_ext     <- c(ncp_ext, int_end + 1)
      pow_ext     <- c(pow_ext, 1)
      pow_dir_ext <- c(pow_dir_ext, 1)
      pow_sel_ext <- c(pow_sel_ext, 1)
      w_inp_ext   <- c(w_inp_ext, ext_pos_inp)
    }

    stopifnot(length(pow_ext)     == length(w_inp_ext),
              length(pow_dir_ext) == length(w_inp_ext),
              length(pow_sel_ext) == length(w_inp_ext))

    w_all_ext <- w_inp_ext / pow_sel_ext
    w_all_ext <- w_all_ext / sum(w_all_ext)

    # EDR: average power over all studies (pre-selection)
    EDR <- sum(w_all_ext * pow_ext)

    # ERR: average directional power over significant studies
    w_sig_ext <- w_all_ext * pow_ext
    w_sig_ext <- w_sig_ext / sum(w_sig_ext)
    ERR <- sum(w_sig_ext * pow_dir_ext)

    ERR <- min(max(ERR, alpha / 2), 1)
    EDR <- min(max(EDR, alpha),     1)

    # ---- mixture density at each grid point (grid-only ncp and w_all) ----
    # test_stat is folded/absolute, so the density of |Z| (or |T|) for a
    # component at mu needs the mirror term when folded = TRUE. build_dens
    # handles the fold, the component SDs (zsds), and the z-vs-t choice (df).
    Dens    <- build_dens(d_x = grid, ncp = ncp, z_sd = z_sd, df = df,
                          curve_type = if (is.null(df)) "z" else "t",
                          folded = folded)

    lik_mat <- t(Dens)                       # rows = grid, cols = components
    w_grid  <- as.vector(lik_mat %*% w_all)

    # local power: posterior-weighted average of component power
    # (numerator and denominator share the same w_all-based lik_mat, so this
    #  ratio is scale-invariant -- no normalization of w_grid needed or wanted)
    pow_grid <- as.vector((lik_mat %*% (w_all * pow)) / w_grid)

    # bin into intervals, density-weighted average within each bin
    int  <- seq(x_lim_min, x_lim_max, by = int_loc)
    bins <- cut(grid, breaks = int, include.lowest = TRUE)

    local_power <- as.vector(tapply(seq_along(grid), bins, function(idx) {
      denom <- sum(w_grid[idx])
      if (denom <= 0 || !is.finite(denom)) return(NA_real_)
      sum(pow_grid[idx] * w_grid[idx]) / denom
    }))

    midpoints          <- (int[-length(int)] + int[-1]) / 2
    names(local_power) <- paste0("lp.", midpoints)


    # ---- heterogeneity of the non-centrality distribution (z-scale) --------
    # (a) SD of the fitted mixture (grid-only components, so the placeholder
    #     extreme anchors don't inflate it)
    ncp_mean  <- sum(w_all * ncp)
    ncp_var   <- sum(w_all * ncp^2) - ncp_mean^2
    ncp_tau_a <- sqrt(max(ncp_var, 0))

    # (b) method-of-moments: Var(observed z) - 1 (remove sampling variance)
    wg        <- w_grid / sum(w_grid)          # w_grid is an unnormalized density here
    z_mean    <- sum(wg * grid)
    z_var     <- sum(wg * (grid - z_mean)^2) - 1
    ncp_tau_b <- sqrt(max(z_var, 0))

    # conservative combined estimate
    ncp_tau <- min(ncp_tau_a, ncp_tau_b)

    # shape diagnostic: predicted vs observed nonsignificant z
    upper <- qnorm(1 - alpha)
    lower <- 0
    idx   <- grid > lower & grid < upper

    d_med  <- NA_real_
    d_mean <- NA_real_
    if (sum(idx) > 1) {
      X_ns <- grid[idx]
      pred <- w_grid[idx]
      pred <- pred / (sum(pred) * xz)

      pred_cdf  <- cumsum(pred) * xz
      pred_cdf  <- pred_cdf / pred_cdf[length(pred_cdf)]
      pred_med  <- approx(pred_cdf, X_ns, xout = 0.5, rule = 2)$y
      pred_mean <- sum(X_ns * pred) * xz

      obs <- test_stat[test_stat > lower & test_stat < upper]
      if (length(obs) > 0) {
        d_med  <- median(obs) - pred_med
        d_mean <- mean(obs)   - pred_mean
      }
    }

    list(
      EDR            = EDR,
      ERR            = ERR,
      w_all          = w_all,
      w_all_ext      = w_all_ext,
      pow_ext        = pow_ext,
      pow_dir_ext    = pow_dir_ext,
      pow_sel_ext    = pow_sel_ext,
      ncp_ext        = ncp_ext,
      grid           = grid,
      w_grid         = w_grid,
      pow_grid       = pow_grid,
      local_power    = local_power,
      ncp_tau        = ncp_tau,
      ncp_tau_a      = ncp_tau_a,
      ncp_tau_b      = ncp_tau_b,
      shape_d_median = d_med,
      shape_d_mean   = d_mean
    )

  } ### EOF compute_power


  ################################################################
  ### Main Power # ppp2
  ################################################################

  xz     <- 0.01
  grid <- seq(x_lim_min, x_lim_max, by = xz)

  est_out <- compute_power(
    est_inp     = est_inp,
    int_beg     = int_beg,
    int_end     = int_end,
    grid        = grid,
    test_stat   = test_stat,
    crit        = crit,
    directional = directional,
    folded      = folded,
    alpha       = alpha,
    x_lim_min   = x_lim_min,
    x_lim_max   = x_lim_max,
    int_loc     = int_loc
   )

  est_out$local_es      <- NULL
  est_out$es_mean_all   <- NULL
  est_out$es_median_all <- NULL
  est_out$es_mean_sig   <- NULL
  est_out$es_tau        <- NULL

  if (!is.null(yi)) {

    sei_est <- estimate_sei_function(
      test_stat = test_stat,
      sei       = sei,
      grid      = grid,
      int_beg   = int_beg,
      int_end   = int_end
    )

   es_est <- get_es_estimates(
      sei_grid    = sei_est[,2],
      grid        = est_out$grid,
      w_grid      = est_out$w_grid,
      ncp_ext     = est_out$ncp_ext,
      w_all_ext   = est_out$w_all_ext,
      test_stat   = test_stat,  
      yi          = yi,         
      directional = directional,
      folded      = folded,
      crit        = crit,
      int_beg     = int_beg,
      int_end     = int_end,
      int_loc     = int_loc,
      x_lim_min   = x_lim_min,
      x_lim_max   = x_lim_max
    )
 
    est_out$local_es      <- es_est$local_es
    est_out$es_mean_all   <- es_est$es_mean_all
    est_out$es_median_all <- es_est$es_median_all   # FIX: now propagated
    est_out$es_mean_sig   <- es_est$es_mean_sig
    est_out$es_tau        <- es_est$es_tau

  } # EOF es 

  est_out
 
  return(est_out)
 
} ### EOF Compute Estimates
 



############################################################
### SLOPE DIAGNOSTICS
############################################################

get_slope <- function(int, from, to, bw, spline_width = .5) {

  d <- get_densities(int, bw = bw, from = from, to = to,
       width = spline_width, augment = TRUE)
  fit <- lm(d$ZY ~ d$ZX)
  slope <- unname(coef(fit)[2])
  se <- unname(summary(fit)$coefficients[2, 2])
  ci <- confint(fit)[2, ]
  
  out <- c(slope = slope, se = se, ci.lb = ci[1], ci.ub = ci[2])
  out
} # EOF get_slope


######################################################
### END OF SLOPE DIAGNOSTICS
######################################################



##################################################
### ROOT OF FUNCTION
##################################################


#ooo

run_OF <- function(test_stat,yi=NULL,sei=NULL,ODR,int_beg,int_end,ncp,z_sd, 
                    cola = "springgreen4",bw_est = .2,
                    dens_pre = NULL, d_x_pre = NULL, width = 1, 
                    ncp_fixed, z_sd_fixed, w_fixed,
                    x_lim_min,x_lim_max) {


  testing = FALSE
  #testing = TRUE


  int = test_stat[test_stat >= int_beg & test_stat <= int_end]
  #hist(int)

  ODR = mean(test_stat > 1.96)

  k_ncp <- length(ncp)

  #from = int_beg; to = int_end; width = 1
  ## ---- Observed density ----
  dens <- get_densities(int, bw = bw_est, from = int_beg, 
                         to = int_end, width = width, augment = augment)

  d_x       <- dens[, 1]
  o_d_y     <- dens[, 2]
  n_bars    <- length(d_x) 
  bar_width <- d_x[2] - d_x[1]

  use_pre_dens <- k_ncp > 1 &&
                     ncp_fixed == TRUE &&
                     z_sd_fixed == TRUE


  ## ---- theta = [weights, ncp, z_sd] — always ----

  startval <- c(rep(1/k_ncp, k_ncp), ncp, z_sd)

  lowlim   <- c(if (k_ncp == 1) 1 else rep(0, k_ncp),
                if (ncp_fixed) ncp else rep(0,k_ncp),
                if (z_sd_fixed) z_sd else rep(1,k_ncp)
			)
  highlim  <- c(rep(1, k_ncp),
                if (ncp_fixed) ncp else rep(max(ncp),k_ncp),
                if (z_sd_fixed) z_sd else rep(10,k_ncp)
			)

  ## ---- Positions (always the same) ----
  idx_wt  <- 1:k_ncp
  idx_ncp <- (k_ncp + 1):(2 * k_ncp)
  idx_z_sd <- (2 * k_ncp + 1):(3 * k_ncp)


  curve_fitting <- function(theta) {

    wt <- theta[idx_wt]
    wt <- wt / sum(wt)
 
    if (use_pre_dens && ncp_fixed && z_sd_fixed) {
      dens_now <- dens_pre
    } else {
      dens_now <- build_dens(d_x, theta[idx_ncp], theta[idx_z_sd], df, 
           curve_type, folded=folded)
    }
    dens_now <- dens_now / (rowSums(dens_now) * bar_width)

    e_d_y <- colSums(dens_now * wt)

    rmse = sqrt(mean((e_d_y - o_d_y)^2))

    scale = 1
    ## ---- Diagnostic plot ----
    if (testing && runif(1) > 0.6) {
      plot(d_x, o_d_y * scale, type = 'l',
           xlim = c(x_lim_min, x_lim_max), ylim = c(0, ymax),
           xlab = "", ylab = "", axes = FALSE)
      lines(d_x, e_d_y * scale, col = "red3")
    }


    rmse

    } #EOF curve_fitting

  ## ---- Optimize ----
  auto <- nlminb(startval, curve_fitting,
                 lower = lowlim, upper = highlim,
                 control = list(eval.max = 1000, abs.tol = 1e-10))

  of_fit =  list(
    w_inp   = auto$par[idx_wt]/sum(auto$par[idx_wt]),
    ncp     = auto$par[idx_ncp],
    z_sd    = auto$par[idx_z_sd],
    fit     = auto$objective
  )

  of_fit$w_inp

  est_inp = of_fit
 
  ## ---- Compute power for each bootstrap result ----
  res_est = get_estimates(
        est_inp = est_inp, 
        test_stat = test_stat,
        int_beg = int_beg,
        int_end = int_end, 
        yi = yi, 
        sei = sei,
        x_lim_min = x_lim_min,
        x_lim_max = x_lim_max)

  ret = list(
    w_inp           = est_inp$w_inp,
    ncp             = est_inp$ncp,
    z_sd            = est_inp$z_sd,
    model_fit       = est_inp$fit,
    ODR             = ODR,   
    EDR             = res_est$EDR,
    ERR             = res_est$ERR,
    w_all           = res_est$w_all, 
    local_power     = res_est$local_power,
    ncp_tau         = res_est$ncp, 
    shape_d_median  = res_est$shape_d_median,
    shape_d_mean    = res_est$shape_d_mean,
    local_es        = res_est$local_es,
    es_mean_all     = res_est$es_mean_all,
    es_median_all   = res_est$es_median_all,
    es_mean_sig     = res_est$es_mean_sig,
    es_tau          = res_est$es_tau
  )

  ret
  
  return(ret)


} ### EOF run_OF



##################################################
### ROOT EM FUNCTION
##################################################


#eee

run_EM <- function(
  test_stat,
  yi, 
  sei, 
  int_beg,
  int_end,
  ncp,
  z_sd,
  ncp_fixed,
  z_sd_fixed,
  w_fixed,
  x_lim_min,
  x_lim_max,
  max_iter   = 500,
  tol        = 1e-8) {


n_starts = 1

components = length(ncp)
w.inp = rep(1/components,components)

int = test_stat[test_stat >= int_beg & test_stat <= int_end]
ODR = mean(test_stat > 1.96);ODR

EM_EXT_SIGNED <- function(
  int,
  ncp,
  z_sd,
  int_beg = -Inf,        # signed lower window bound (e.g. -8, or -Inf for none)
  int_end =  Inf,        # signed upper window bound
  ncp_fixed,
  z_sd_fixed,
  w_fixed,
  max_iter = 1000,
  tol = 1e-8
) {

  x <- int
  x <- x[!is.na(x)]
  k.int <- length(x)
  components <- length(ncp)

  loglik_old <- -Inf
  loglik_trace <- numeric(max_iter)

  # ---- Newton update for a single component's mean (signed, truncated) ----
  update_ncp_newton_signed <- function(
    mu_init, sigma, tau_k, x, int_beg, int_end,
    max_iter = 20, tol = 1e-8, lower = -10, upper = 10
  ) {
    mu  <- mu_init
    n_k <- sum(tau_k)
    for (iter in 1:max_iter) {
      # no folding: sufficient statistic is just tau-weighted sum of x
      S1 <- sum(tau_k * x)

      # single truncation interval [int_beg, int_end] on the signed axis
      a <- (int_beg - mu) / sigma
      b <- (int_end - mu) / sigma
      Cmu <- pnorm(b) - pnorm(a)
      if (Cmu <= 0 || !is.finite(Cmu)) return(mu)

      phi_a <- dnorm(a); phi_b <- dnorm(b)
      delta <- (phi_a - phi_b) / Cmu
      kappa <- (a * phi_a - b * phi_b) / Cmu

      U <- (S1 - n_k * mu) / sigma^2 - n_k * delta / sigma
      H <- -n_k / sigma^2 - n_k / sigma^2 * (delta^2 + kappa)
      if (H >= 0) H <- -n_k / sigma^2

      step   <- U / H
      mu_new <- max(lower, min(upper, mu - step))
      if (abs(mu_new - mu) < tol) break
      mu <- mu_new
    }
    return(mu)
  }

  # ---- Newton update for a single component's sd (signed, truncated) ----
  update_z_sd_newton_signed <- function(
    sigma_init, mu, tau_k, x, int_beg, int_end,
    max_iter = 20, tol = 1e-6, lower = 0.5, upper = 10
  ) {
    sigma <- sigma_init
    n_k   <- sum(tau_k)
    for (iter in 1:max_iter) {
      # no folding: second moment about mu directly
      S2 <- sum(tau_k * (x - mu)^2)

      a <- (int_beg - mu) / sigma
      b <- (int_end - mu) / sigma
      Cmu <- pnorm(b) - pnorm(a)
      if (Cmu <= 0 || !is.finite(Cmu)) return(sigma)

      phi_a <- dnorm(a); phi_b <- dnorm(b)
      kappa <- (a * phi_a - b * phi_b) / Cmu

      U <- S2 / sigma^3 - n_k / sigma - n_k * kappa / sigma
      H <- -2 * n_k / sigma^2

      step      <- U / H
      sigma_new <- max(lower, min(upper, sigma - step))
      if (abs(sigma_new - sigma) < tol) break
      sigma <- sigma_new
    }
    return(sigma)
  }

  ######## EM start #####
  for (iter in 1:max_iter) {

    # ---- E-STEP (signed, no fold) ----
    log_g <- matrix(NA_real_, k.int, components)
    for (k in 1:components) {
      mu    <- ncp[k]
      sigma <- z_sd[k]

      # single signed density, no + dnorm(x, -mu)
      log_dens <- dnorm(x, mu, sigma, log = TRUE)

      # single truncation interval on the signed axis
      a <- (int_beg - mu) / sigma
      b <- (int_end - mu) / sigma
      norm_const <- pnorm(b) - pnorm(a)
      if (norm_const <= 0 || !is.finite(norm_const)) return(NULL)

      log_g[, k] <- log_dens - log(norm_const)
    }

    log_num <- sweep(log_g, 2, log(w.inp), "+")
    log_denom <- apply(log_num, 1, function(v) {
      m <- max(v); m + log(sum(exp(v - m)))
    })
    tau <- exp(log_num - log_denom)

    # ---- M-STEP ----
    n_k <- colSums(tau)
    if (!W_FIXED) {
      w.inp <- n_k / k.int
      w.inp <- pmax(w.inp, 1e-12)
      w.inp <- w.inp / sum(w.inp)
    }
    if (!NCP_FIXED) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        ncp[k] <- update_ncp_newton_signed(
          mu_init = ncp[k], sigma = z_sd[k], tau_k = tau[, k],
          x = x, int_beg = int_beg, int_end = int_end
        )
      }
    }
    if (!z_sd_FIXED) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        z_sd[k] <- max(1, update_z_sd_newton_signed(
          sigma_init = z_sd[k], mu = ncp[k], tau_k = tau[, k],
          x = x, int_beg = int_beg, int_end = int_end
        ))
      }
    }

    # ---- LOG-LIKELIHOOD ----
    loglik <- sum(log_denom)
    loglik_trace[iter] <- loglik
    if (!is.finite(loglik)) return(NULL)
    if (abs(loglik - loglik_old) < tol) break
    loglik_old <- loglik
  }

  list(w.inp = w.inp, ncp = ncp, z_sd = z_sd, fit = loglik)
}




EM_EXT_FOLDED <- function(
  int,
  ncp,
  z_sd,
  int_beg = 1.96,
  int_end = max(int),
  ncp_fixed,
  z_sd_fixed,
  w_fixed,
  max_iter = 1000,
  tol = 1e-8
) {


x <- int
x <- x[!is.na(x)]
k.int <- length(x)
components <- length(ncp)

loglik_old <- -Inf
loglik_trace <- numeric(max_iter)



update_ncp_newton_folded <- function(
  mu_init,
  sigma,
  tau_k,
  x,
  int_beg,
  int_end,
  max_iter = 20,
  tol      = 1e-8,
  lower    = -10,
  upper    = 10
) {
  mu  <- mu_init
  n_k <- sum(tau_k)

  for (iter in 1:max_iter) {
    # --- posterior fold probabilities (depend on current mu) ---
    f_pos <- dnorm(x,  mu, sigma)
    f_neg <- dnorm(x, -mu, sigma)
    p_pos <- f_pos / (f_pos + f_neg)

    # signed sufficient statistic
    r         <- 2 * p_pos - 1
    S1_signed <- sum(tau_k * x * r)

    # --- truncation correction ---
    a1 <- ( int_beg - mu) / sigma
    b1 <- ( int_end - mu) / sigma
    a2 <- (-int_end - mu) / sigma
    b2 <- (-int_beg - mu) / sigma

    Cmu <- (pnorm(b1) - pnorm(a1)) +
           (pnorm(b2) - pnorm(a2))
    if (Cmu <= 0 || !is.finite(Cmu)) return(mu)

    phi_a1 <- dnorm(a1); phi_b1 <- dnorm(b1)
    phi_a2 <- dnorm(a2); phi_b2 <- dnorm(b2)

    delta <- (phi_a1 - phi_b1 + phi_a2 - phi_b2) / Cmu
    kappa <- (a1*phi_a1 - b1*phi_b1 +
              a2*phi_a2 - b2*phi_b2) / Cmu

    # --- score ---
    U <- (S1_signed - n_k * mu) / sigma^2 -
         n_k * delta / sigma

    # --- Hessian (with folding term) ---
    p_neg  <- 1 - p_pos
    H_fold <- 4 / sigma^4 * sum(tau_k * x^2 * p_pos * p_neg)
    H      <- H_fold - n_k / sigma^2 -
              n_k / sigma^2 * (delta^2 + kappa)

    if (H >= 0) H <- -n_k / sigma^2   # safeguard

    step   <- U / H
    mu_new <- max(lower, min(upper, mu - step))

    if (abs(mu_new - mu) < tol) break
    mu <- mu_new
  }
  return(mu)
}


update_z_sd_newton_folded <- function(
  sigma_init,
  mu,
  tau_k,
  x,
  int_beg,
  int_end,
  max_iter = 20,
  tol      = 1e-6,
  lower    = 0.5,
  upper    = 10
) {
  sigma <- sigma_init
  n_k   <- sum(tau_k)

  for (iter in 1:max_iter) {
    # --- posterior fold probabilities (depend on current sigma) ---
    f_pos <- dnorm(x,  mu, sigma)
    f_neg <- dnorm(x, -mu, sigma)
    p_pos <- f_pos / (f_pos + f_neg)
    p_neg <- 1 - p_pos

    # signed second moment
    S2_signed <- sum(tau_k * (p_pos * (x - mu)^2 +
                              p_neg * (x + mu)^2))

    # --- truncation correction ---
    a1 <- ( int_beg - mu) / sigma
    b1 <- ( int_end - mu) / sigma
    a2 <- (-int_end - mu) / sigma
    b2 <- (-int_beg - mu) / sigma

    Cmu <- (pnorm(b1) - pnorm(a1)) +
           (pnorm(b2) - pnorm(a2))
    if (Cmu <= 0 || !is.finite(Cmu)) return(sigma)

    phi_a1 <- dnorm(a1); phi_b1 <- dnorm(b1)
    phi_a2 <- dnorm(a2); phi_b2 <- dnorm(b2)

    kappa <- (a1*phi_a1 - b1*phi_b1 +
              a2*phi_a2 - b2*phi_b2) / Cmu

    # --- score ---
    U <- S2_signed / sigma^3 - n_k / sigma -
         n_k * kappa / sigma

    # --- Hessian (conservative approximation) ---
    H <- -2 * n_k / sigma^2

    step      <- U / H
    sigma_new <- max(lower, min(upper, sigma - step))

    if (abs(sigma_new - sigma) < tol) break
    sigma <- sigma_new
  }
  return(sigma)
}

    ######## EM start ##### 

    for (iter in 1:max_iter) {

    # ============================================================
    # E-STEP
    # ============================================================
    log_g <- matrix(NA_real_, k.int, components)
    for (k in 1:components) {
      mu    <- ncp[k]
      sigma <- z_sd[k]
      l1 <- dnorm(x,  mu, sigma, log = TRUE)
      l2 <- dnorm(x, -mu, sigma, log = TRUE)
      m  <- pmax(l1, l2)
      log_fold <- m + log(exp(l1 - m) + exp(l2 - m))
      a1 <- ( int_beg - mu) / sigma
      b1 <- ( int_end - mu) / sigma
      a2 <- (-int_end - mu) / sigma
      b2 <- (-int_beg - mu) / sigma
      norm_const <-
        (pnorm(b1) - pnorm(a1)) +
        (pnorm(b2) - pnorm(a2))
      if (norm_const <= 0 || !is.finite(norm_const))
        return(NULL)
      log_g[, k] <- log_fold - log(norm_const)
    }
    log_num <- sweep(log_g, 2, log(w.inp), "+")
    log_denom <- apply(log_num, 1, function(v) {
      m <- max(v)
      m + log(sum(exp(v - m)))
    })
    tau <- exp(log_num - log_denom)
    # ============================================================
    # M-STEP
    # ============================================================
    n_k <- colSums(tau)
    if (!W_FIXED) {
      w.inp <- n_k / k.int
      w.inp <- pmax(w.inp, 1e-12)
      w.inp <- w.inp / sum(w.inp)
    }
    if (!NCP_FIXED) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        ncp[k] <- update_ncp_newton_folded(
          mu_init = ncp[k],
          sigma   = z_sd[k],
          tau_k   = tau[, k],
          x       = x,
          int_beg = int_beg,
          int_end = int_end
        )
      }
    }
    if (!z_sd_fixed) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        z_sd[k] <- max(1, update_z_sd_newton_folded(
          sigma_init = z_sd[k],
          mu         = ncp[k],
          tau_k      = tau[, k],
          x          = x,
          int_beg    = int_beg,
          int_end    = int_end
        ))
      }
    }

    # ============================================================
    # LOG-LIKELIHOOD
    # ============================================================

    loglik <- sum(log_denom)
    loglik_trace[iter] <- loglik

    if (!is.finite(loglik))
      return(NULL)

    if (abs(loglik - loglik_old) < tol)
      break

    loglik_old <- loglik
  }

  list(
    w.inp        = w.inp,
    ncp          = ncp,
    z_sd         = z_sd,
    fit          = loglik
  )

  
}



### Main
### eee2

  components <- length(ncp)
  ncp_init   <- ncp
  z_sd_init  <- z_sd

  best_fit <- NULL
  best_ll  <- -Inf

  print("Running EM")
  print(paste0("Folded: ",Folded))

  for (s in 1:n_starts) {

      raw <- rgamma(components, shape = 1)   # random Dirichlet
      w_s <- raw / sum(raw)

      if(Folded) {

        fit <- EM_EXT_FOLDED(
          int        = int,
          ncp        = ncp,
          z_sd       = z_sd,
          int_beg    = int_beg,
          int_end    = int_end,
          ncp_fixed  = ncp_fixed,
          zsds_fixed = zsds_fixed,
          w_fixed    = w_fised,
          max_iter   = max_iter,
          tol        = tol
        ) 

      } else {

        fit <- EM_EXT_SIGNED(
          int        = int,
          ncp        = ncp,
          z_sd       = z_sd,
          int_beg    = int_beg,
          int_end    = int_end,
          ncp_fixed  = ncp_fixed,
          z_sd_fixed = z_sd_fixed,
          w_fixed    = w_fixed,
          max_iter   = max_iter,
          tol        = tol
         ) 

      }

    if (!is.null(fit) && is.finite(fit$fit) && fit$fit > best_ll) {
      best_ll  <- fit$fit
      best_fit <- fit
      best_fit$start <- s   # track which start won
    }
  }

  em.est = best_fit

  #print(em.est)

  cp.inp = em.est

  ## ---- Compute power
  cp.res = get_estimates(
        cp.inp = cp.inp, 
        test_stat = test_stat,
        int_beg = int_beg,
        int_end = int_end, 
        yi = yi, 
        sei = sei,
        x_lim_min = x_lim_min,
        x_lim_max = x_lim_max)

  ret = list(
    w.inp           = em.est$w.inp,
    ncp             = em.est$ncp,
    z_sd            = em.est$z_sd,
    model.fit       = em.est$fit,
    ODR             = ODR,   
    EDR             = cp.res$EDR,
    ERR             = cp.res$ERR,
    w_all           = cp.res$w_all, 
    local_power     = cp.res$local_power,
    shape_d_median  = cp.res$shape_d_median,
    shape_d_mean    = cp.res$shape_d_mean,
    local_es        = cp.res$local_es,
    es_mean_all     = cp.res$es_mean_all,
    es_median_all   = cp.res$es_median_all,
    es_mean_sig     = cp.res$es_mean_sig,
    es_tau          = cp.res$es_tau
  )

  ret

  return(ret)

}




##################################################
### Bootstrap
##################################################

### BOOTSTRAP FUNCTION FOR CI

#xxx

run_bootstrap <- function(
                   test_stat,
                   yi,
                   sei,
                   crit,   
                   ncp,
                   z_sd,
                   int_beg,
                   int_end,
                   directional,
                   folded,
                   dens_pre,
                   d_x_pre,
                   width = 1,
                   boot_iter,
                   ncp_fixed,
                   z_sd_fixed,
                   w_fixed,
                   cluster_id,
                   est_method,
                   x_lim_min,
                   x_lim_max,
                   augment
                  ) {

    

      k_ncp = length(ncp)

      if (!is.null(cluster_id)) {
   
          ## ---- Cluster Bootstrap ----

          stopifnot(length(test_stat) == length(cluster_id))	
          z_cluster_rows  <- split(seq_along(test_stat), cluster_id)
          z_cluster.names <- names(z_cluster_rows)
          n_z_clusters    <- length(z_cluster_names)

      } else {

          z_cluster_rows = NULL; z_cluster_names = NULL; n_z_clusters = NULL

      }

 
      var_list = c("est_method","test_stat","crit",
           "cluster_id","z_cluster_rows","z_cluster_names","n_z_clusters",
           "ncp", "z_sd","int_beg", "int_end", 
           "curve_type","directional","folded","x_lim_min", "x_lim_max", 
           "ncp_fixed", "z_sd_fixed","w_fixed",
           "get_estimates","alpha","int_loc"
           )


      if (!is.null(yi)) var_list = c(var_list,"yi","sei")

      boot_res = NULL

      if (est_method == "OF") {

         var_list = c(var_list,"run_OF","build_dens", "get_densities",
                  "bw_est","width","dens_pre", "d_x_pre","augment")

         ## ---- Set up cluster ----
         ncores <- max(1, parallel::detectCores() * .8)

         cl <- parallel::makeCluster(ncores)
         on.exit(parallel::stopCluster(cl), add = TRUE)

         parallel::clusterExport(cl, varlist = var.list, 
           envir = environment() )
         invisible(parallel::clusterEvalQ(cl, library(KernSmooth)))
         invisible(parallel::clusterEvalQ(cl, testing <- FALSE))
         parallel::clusterSetRNGStream(cl, iseed = 42)

         boot_res <- parallel::parLapply(cl, 1:boot_iter, function(i) {

         if (!is.null(cluster_id)) {
              sampled  <- sample(z_cluster_names, n_z_clusters, replace = TRUE)
              rows     <- unlist(z_cluster_rows[sampled], use_names = FALSE)
              boot_all <- test_stat[rows]
          } else {
              rows     <- sample(seq_along(test_stat), length(test_stat), replace = TRUE)
              boot_all <- test_stat[rows]
          }

          boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]

          if (!is.null(yi)) { yi_boot <- yi[rows]; sei_boot <- sei[rows] 
           } else { yi_boot = NULL; sei_boot = NULL }

          fit_res <- #tryCatch(
                 run_OF(
                   test_stat  = boot_all,
                   yi         = yi, # yi_boot,
                   sei        = sei, #, _boot,
                   int_beg    = int_beg,
                   int_end    = int_end,
                   ncp        = ncp,
                   z_sd       = z_sd,
                   dens_pre   = dens_pre,
                   d_x_pre    = d_x_pre,
                   bw_est     = bw_est, 
                   width      = width,
                   ncp_fixed  = ncp_fixed,
                   z_sd_fixed = z_sd_fixed,
                   w_fixed    = w_fixed,
                   x_lim_min  = x_lim_min,
                   x_lim_max  = x_lim_max,
                   augment    = augment
                 )
	        #  ,error = function(e) NULL )

          fit_res
  
      }) # EOF boot_res

      if(is.null(boot_res)) stop("boostrap failure.")

    } # EOF OF

     # boot_res

      if (est_method == "EM") {

         var.list = c(var.list,"run_EM")

         ## ---- Set up cluster ----
         ncores <- max(1, parallel::detectCores() * .8)
         cl <- parallel::makeCluster(ncores)
         on.exit(parallel::stopCluster(cl), add = TRUE)

         parallel::clusterExport(cl, varlist = var.list, 
           envir = environment() )

          parallel::clusterEvalQ(cl, library(stats))


         boot_res <- parallel::parLapply(cl, 1:boot.iter, function(i) {


         if (!is.null(cluster.id)) {
              sampled  <- sample(z.cluster.names, n.z.clusters, replace = TRUE)
              rows     <- unlist(z.cluster.rows[sampled], use.names = FALSE)
              boot_all <- test_stat[rows]
          } else {
             rows     <- sample(seq_along(test_stat), length(test_stat), replace = TRUE)
             boot_all <- test_stat[rows]
          }

          boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]

          if (!is.null(yi)) { yi_boot <- yi[rows]; sei_boot <- sei[rows] 
           } else { yi_boot = NULL; sei_boot = NULL }

          fit.EM <- #tryCatch(
                 run_EM(
                 test_stat  = boot_all,
                 yi         = yi, # yi_boot,
                 sei        = sei, #, _boot,
                 int_beg    = int_beg,
                 int_end    = int_end,
                 ncp        = ncp,
                 z_sd       = z_sd,
                 ncp_fixed  = ncp_fixed,
                 z_sd_fixed = z_sd_fixed,
                 w_fixed    = w_fixed,
                 x_lim_min  = x_lim_min,
                 x_lim_max  = x_lim_max
                 )
	         #,error = function(e) NULL )

          fit.EM 

         }) # close boot_res

    } # EOF EM 

    #boot_res 

    #est_method = "ML_gamma"; boot.iter = 500
    if (est_method == "ML_gamma") {

         n_starts = 2
         if (est_method == "ML_gamma") var.list = c(var.list,"run_ML_gamma","n_starts")

         ## ---- Set up cluster ----
         ncores <- max(1, parallel::detectCores() * .8)
         cl <- parallel::makeCluster(ncores)
         on.exit(parallel::stopCluster(cl), add = TRUE)

         parallel::clusterExport(cl, varlist = var.list, 
           envir = environment() )
         invisible(parallel::clusterEvalQ(cl, library(KernSmooth)))
         invisible(parallel::clusterEvalQ(cl, TESTING <- FALSE))
         parallel::clusterSetRNGStream(cl, iseed = 42)

         boot_res <- parallel::parLapply(cl, 1:boot.iter, function(i) {

           if (!is.null(cluster.id)) {
               sampled <- sample(z.cluster.names, n.z.clusters, replace = TRUE)
               boot_all <- unlist(z.clusters[sampled], use.names = FALSE)
               odr_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]
           } else {
               boot_all <- sample(test_stat, length(test_stat), replace = TRUE)
               odr_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= int_beg & boot_all <= int_end]
           }


           gamma_fit <- tryCatch(
                 run_ML_gamma(boot_all, int_beg = int_beg, int_end = int_end, crit = crit,
                   n_starts = n_starts)
               ,error = function(e) NULL )

           fit.res = list(
             ODR    = odr_boot,
             gamma_edr = if (!is.null(gamma_fit)) gamma_fit$gamma_edr[1] else NA,
             gamma_err = if (!is.null(gamma_fit)) gamma_fit$gamma_err[1] else NA
             #fit = gamma_fit$negloglik
             )

          fit.res
  
      }) # EOF boot_res

   } # EOF ML_gamma 


   boot_res <- boot_res[!sapply(boot_res, is.null)]

   #summary(sapply(boot_res, function(x) x$gamma_edr))
   #summary(sapply(boot_res, function(x) x$gamma_err))

   #length(boot_res)
   if (length(boot_res) == 0)
      stop("All bootstrap runs failed.")

  return(boot_res)

}


##################################################
### End of New Extended Zcurve
##################################################

#zzz

run_zcurve = function(est_method,test_stat,crit,int_beg,int_end,
           ncp,z_sd,ncp_fixed,z_sd_fixed,w_fixed,
           cluster_id,boot.iter=0,edr_ci_adj=0,err_ci_adj=0,
           yi=NULL,sei=NULL,folded,directional,alpha,int_loc,df,
           x_lim_min,x_lim_max,spline_width,augment) {


  zcurve_time = system.time({

  int = test_stat[test_stat >= int_beg & test_stat <= int_end]

  ODR = mean(test_stat > crit)

  width = spline_width

  d_x_pre  <- NULL
  dens_pre <- NULL

 if(est_method == "OF") { 

 ## ---- Precompute Dens if both fixed ----
 if (ncp_fixed && z_sd_fixed) {
      densy     <- get_densities(int, bw = bw_est, from = int_beg,
                             to = int_end, width = width, augment = augment)
      d_x_pre   <- densy[, 1]
      dens_pre  <- build_dens(d_x_pre, ncp, z_sd, df, 
           curve_type = curve_type, folded=folded)
      bar_width <- d_x_pre[2] - d_x_pre[1]
      dens_pre  <- dens_pre / (rowSums(dens_pre) * bar_width)
 }

         res_of <- run_OF(
                 test_stat  = test_stat, 
                 yi         = yi,
                 sei        = sei,
                 int_beg    = int_beg,
                 int_end    = int_end,   
                 ncp        = ncp,
                 z_sd       = z_sd,
                 dens_pre   = dens_pre,
                 d_x_pre    = d_x_pre,
                 bw_est     = bw_est,
                 width      = width,
                 ncp_fixed  = ncp_fixed,
                 z_sd_fixed = z_sd_fixed,
                 w_fixed    = w_fixed,
                 x_lim_min  = x_lim_min,
                 x_lim_max  = x_lim_max,
                 augment    = augment
           )

        res_of
        res_run = res_of 
  }  


   if (est_method == "EM") { 

         d_x_pre = NULL
         dens    = NULL

         res.em <- run_EM(
                 test_stat  = test_stat, 
                 yi         = yi,
                 sei        = sei,
                 int_beg    = int_beg,
                 int_end    = int_end,   
                 ncp        = ncp,
                 z_sd       = z_sd,
                 ncp_fixed  = ncp_fixed,
                 z_sd_fixed = z_sd_fixed,
                 w_fixed    = w_fixed,
                 x_lim_min  = x_lim_min,
                 x_lim_max  = x_lim_max
        )


        res.em
        res.run = res.em
     
   }

    res_pe <- list(
      ODR              = ODR,
      EDR              = res_run$EDR,
      ERR              = res_run$ERR,
      ncp              = res_run$ncp,
      z_sd             = res_run$z_sd,
      w_inp            = res_run$w_inp,
      w_all            = res_run$w_all,
      local_power      = res_run$local_power,
      ncp_tau          = res_run$ncp_tau,
      model_fit        = res_run$model_fit,
      ODR_EDR_D        = ODR - res_run$EDR,
      shape_d_median   = res_run$shape_d_median,
      shape_d_mean     = res_run$shape_d_mean,
      local_es_pe      = res_run$local_es,
      es_mean_all_pe   = res_run$es_mean_all,
      es_median_all_pe = res_run$es_median_all,
      es_mean_sig_pe   = res_run$es_mean_sig,
      es_tau_pe        = res_run$es_tau
     )


  #################################

  ### finsh # res.pe

  #################################


  ### start boostrap

  if (boot_iter > 0) {  

       boot_res <- run_bootstrap(
          test_stat   = test_stat, 
          yi          = yi,
          sei         = sei,
          crit        = crit,
          ncp         = res_pe$ncp,
          z_sd        = res_pe$z_sd,
          int_beg     = int_beg,
          int_end     = int_end,
          boot_iter   = boot_iter,
          directional = directional,
          folded      = folded,
          dens_pre    = dens_pre,
          d_x_pre     = d_x_pre,
          ncp_fixed   = ncp_fixed,
          z_sd_fixed  = z_sd_fixed,
          w_fixed     = w_fixed,
          cluster_id  = cluster_id,
          est_method  = est_method,
          x_lim_min   = x_lim_min,
          x_lim_max   = x_lim_max,
          augment     = augment)

          ODR_boot   <- sapply(boot_res, function(x) x$ODR)
          ODR_boot[is.nan(ODR_boot)] = NA
          ODR_ci = c(
             quantile(ODR_boot, probs = .025, na.rm = TRUE),
             quantile(ODR_boot, probs = .975, na.rm = TRUE)
          )


          EDR_boot   <- sapply(boot_res, function(x) x$EDR)

          EDR_boot_lb <- pmax(alpha, EDR_boot - edr_ci_adj)
          EDR_boot_ub <- pmin(1,     EDR_boot + edr_ci_adj)

          EDR_ci = c(
             quantile(EDR_boot_lb, probs = .025, na.rm = TRUE),
             quantile(EDR_boot_ub, probs = .975, na.rm = TRUE)
          )


          ERR_boot   <- sapply(boot_res, function(x) x$ERR)

          ERR_boot_lb <- pmax(alpha, ERR_boot - err_ci_adj)
          ERR_boot_ub <- pmin(1,     ERR_boot + err_ci_adj)

          ERR_ci = c(
             quantile(ERR_boot_lb, probs = .025, na.rm = TRUE),
             quantile(ERR_boot_ub, probs = .975, na.rm = TRUE)
          )


          ODR_EDR_D_ci  <- c(
               quantile(ODR_boot - EDR_boot_lb, .025, na.rm = TRUE),
               quantile(ODR_boot - EDR_boot_ub, .975, na.rm = TRUE)
           )
          ODR_EDR_D = c(res_pe$ODR_EDR_D,ODR_EDR_D_ci)

          if(length(res_pe$ncp) == 1) {
             ncp_ci = quantile(sapply(boot_res, function(x) x$ncp), c(.025, .975))   # CI straight from the percentiles
             ncp = c(res_pe$ncp,ncp_ci)
          } else {
             ncp_ci = apply(sapply(boot_res, function(x) x$ncp),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             ncp  = rbind(res_pe$ncp, ncp_ci) 
          }

          if(length(res_pe$z_sd) == 1) {
             z_sd_ci = quantile(sapply(boot_res, function(x) x$z_sd), c(.025, .975))   # CI straight from the percentiles
             z_sd = c(res_pe$z_sd,z_sd_ci)
          } else {
             z_sd_ci = apply(sapply(boot_res, function(x) x$z_sd),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             z_sd  = rbind(res_pe$z_sd, z_sd_ci) 
          }

          if(length(res_pe$w_all) == 1) {
             w_all_ci = quantile(sapply(boot_res, function(x) x$w_all), c(.025, .975))   # CI straight from the percentiles
             w_all = c(res_pe$w_all,w_all_ci)
          } else {
             w_all_ci = apply(sapply(boot_res, function(x) x$w_all),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             w_all  = rbind(res_pe$w_all, w_all_ci) 
          }

          if(length(res_pe$w_inp) == 1) {
             w_inp_ci = quantile(sapply(boot_res, function(x) x$w_inp), c(.025, .975))   # CI straight from the percentiles
             w_inp = c(res_pe$w_inp,w_inp_ci)
          } else {
             w_inp_ci = apply(sapply(boot_res, function(x) x$w_inp),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             w_inp  = rbind(res_pe$w_inp, w_inp_ci) 
          }

          if(length(res_pe$local_power) == 1) {
             local_power_ci = quantile(sapply(boot_res, function(x) x$local_power), c(.025, .975))   # CI straight from the percentiles
             local_power = c(res_pe$local_power,local_power_ci)
          } else {
             local_power_ci = apply(sapply(boot_res, function(x) x$local_power),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             local_power  = rbind(res_pe$local_power, local_power_ci) 
          }


#          if(!is.na(res_pe$shape_d_median)) {
#            shape_d_median = c(res_pe$shape_d_median,quantile(sapply(boot_res, function(x) x$shape_d_median), c(.025, .975)) )   # CI straight from the percentiles
#            shape_d_mean   = c(res_pe$shape_d_mean,quantile(sapply(boot_res, function(x) x$shape_d_mean), c(.025, .975)) )   # CI straight from the percentiles
#          } else {
#            shape_d_median = c(NA,NA,NA)
#            shape_d_mean   = c(NA,NA,NA)
#          }


          ODR = c(res_pe$ODR,ODR_ci)
          EDR = c(res_pe$EDR,EDR_ci)
          ERR = c(res_pe$ERR,ERR_ci)

          EDR[2] = max(alpha,EDR[2])
          EDR[3] = min(1,EDR[3])

          ERR[2] = max(alpha/2,ERR[2])
          ERR[3] = min(1,ERR[3])

          ###

          if(!is.null(yi)) {

            es_mean_all_boot = sapply(boot_res, function(x) x$es_mean_all)

            es_mean_all_ci = quantile(es_mean_all_boot, probs = c(.025,.975), na.rm=TRUE)
            es_mean_all = c(res_pe$es_mean_all_pe,es_mean_all_ci)

            es_median_all_boot = sapply(boot_res, function(x) x$es_median_all)
            es_median_all_ci = quantile(es_median_all_boot, probs = c(.025,.975), na.rm=TRUE)
            es_median_all = c(res_pe$es_median_all_pe,es_median_all_ci)

            es_tau_boot = sapply(boot_res, function(x) x$es_tau)
            es_tau_ci = quantile(es_tau_boot, probs = c(.025,.975), na.rm=TRUE)
            es_tau = c(res_pe$es_tau_pe,es_tau_ci)


            es_mean_sig_boot = sapply(boot_res, function(x) x$es_mean_sig)

            es_mean_sig_ci = quantile(es_mean_sig_boot, probs = c(.025,.975), na.rm=TRUE)
            es_mean_sig = c(res_pe$es_mean_sig_pe,es_mean_sig_ci)

            es.loc = sapply(boot_res, function(x) x$local_es)

            if(length(res_pe$local_es_pe) == 1) {
               local_es_ci = quantile(es.loc, c(.025, .975),na.rm=TRUE)   # CI straight from the percentiles
               local_es = c(res_pe$local_es_pe,local.es_ci)
            } else {
               local_es_ci = apply(es.loc,1,function(x) quantile(x,c(.025, .975),na.rm=TRUE) )   # CI straight from the percentiles
               local_es  = rbind(res_pe$local_es, local_es_ci) 
            }


          }# EOF if yi 


   } else {
        
          ODR = c(res_pe$ODR,ODR_ci)
          EDR = c(res_pe$EDR,EDR_ci)
          ERR = c(res_pe$ERR,ERR_ci)

          EDR[2] = max(alpha,EDR[2])
          EDR[3] = min(1,EDR[3])

          ERR[2] = ERR[2] - CI.ERR.MIN.ADJ
          ERR[3] = ERR[3] + CI.ERR.MIN.ADJ

          ERR[2] = max(alpha/2,ERR[2])
          ERR[3] = min(1,ERR[3])

        ODR = t(c(res_pe$ODR,rep(NA,2)))
        EDR = t(c(res_pe$EDR,rep(NA,2)))
        ERR = t(c(res_pe$ERR,rep(NA,2)))
        ncp = t(cbind(res_pe$ncp, matrix(NA,length(res_pe$w_all),2)))
        z_sd = t(cbind(res_pe$z_sd, matrix(NA,length(res_pe$w_inp),2)))
        w_all = t(cbind(res_pe$w_all, matrix(NA,length(res_pe$w_all),2)))
        w_inp = t(cbind(res_pe$w_inp, matrix(NA,length(res_pe$w_inp),2)))
        local_power = t(cbind(res_pe$local_power,matrix(NA,length(res_pe$local_power),2)))
        local_es = t(cbind(res_pe$local_es_pe,matrix(NA,length(res_pe$local_es_pe),2)))
        ODR_EDR_D = c(res_pe$ODR_EDR_D,NA,NA)
        shape_d_median  = c(res_pe$shape_d_median,NA,NA)
        shape_d_mean    = c(res_pe$shape_d_mean,NA,NA)
        es_mean_all     = c(res_pe$es_mean_all_pe,NA,NA)
        es_median_all   = c(res_pe$es_median_all_pe,NA,NA)
        es_mean_sig     = c(res_pe$es_mean_sig_pe,NA,NA)
        es_tau          = c(res_pe$es_tau_pe,NA,NA)

  } # EOF no bootstrap

}) # EOF Timer

  zcurve.res = list(
            ODR                = ODR,
            EDR                = EDR,
            ERR                = ERR,
            ncp                = ncp,
            z_sd               = z_sd,
            w_inp              = w_inp,
            w_all              = w_all,
            local_power        = local_power,
            local_es           = local_es,
            ODR_EDR_D          = ODR_EDR_D,
            shape_d_median     = shape_d_median,
            shape_d_mean       = shape_d_mean,
            es_mean_all        = es_mean_all,
            es_median_all      = es_median_all,  
            es_mean_sig        = es_mean_sig,
            es_tau             = es_tau,
            time               = zcurve_time
          )



#zcurve.res

#print(zcurve.res)

return(zcurve.res)

}


############################################

Write_Local_Power = function(local_power) {

  if (int.loc == 0) return(invisible(NULL))

  # Midpoints of each bin — these are the correct label positions
  midpoints = seq(x_lim_min, x_lim_max - int_loc, by = int_loc) + int_loc / 2

  # Format labels
  lab = paste0(round(local_power * 100), "%")

  mtext(lab, side = 1, line = 1.5, at = midpoints, cex = 1.0, las = 1)

} ### EOF Write.Local.Power


# l_p = local_power[1,];l_es = local_es[1,]

Write_Local_Power_ES = function(l_p, l_es) {

  if (int.loc == 0) return(invisible(NULL))

  # Midpoints of each bin — these are the correct label positions
  midpoints = seq(x_lim_min, x_lim_max - int_loc, by = int_loc) + int_loc / 2

  #l_es = zres.tmz$local_es[1,]
   l_es

  # Format labels
  lab1 <- paste0(round(l_p * 100), "%")
  lab2 <- sprintf("%.2f", l_es)

  lab1 = c("local power",lab1)
  lab2 = c("local es",lab2)

  location1 = c(-1,midpoints)
  location2 = c(-1,midpoints-.05)

  mtext(lab1, side = 1, line = -.5, at = location1, cex = 1.0, las = 1)
  mtext(lab2, side = 1, line =  .5, at = location2, cex = 1.0, las = 1)

} ### EOF Write.Local.Power



###################################################
#### Begin Draw Histogram
###################################################

Draw.Histogram.3.8 = function(w,cola="blue3",
	Write.CI = FALSE) {

	#w = w.all;cola = "blue3"; 

	z.hist = test_stat[test_stat > x.lim.min & test_stat < x.lim.max -.04] + .04
   

	int.start = round(int_beg,1)
	if (round(int_beg,2) == 1.96) {
         Int.start = 2
	}

	n.breaks = seq(x.lim.min,x.lim.max,hist.bar_width);n.breaks

	par(cex.axis = 1)
	par(cex.lab = 1)
	par(family = "sans") 
	par(font.axis = 1) #"Microsoft Himalaya")
	par(font.lab = 1) # Microsoft Himalaya")

	col1 = adjustcolor(cola, alpha.f = 0.2)
	col2 = adjustcolor(cola, alpha.f = 0.3)

    if(Show.Y.Label) y.label = "Density" else y.label = ""

	hist(z.hist,breaks=n.breaks,freq=FALSE,
		col=col1,border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(0,ymax),
		ylab=y.label,xlab="",axes=FALSE,main=Title,lwd=1)

    if(Show.Y.Label) {

  	  axis(2, ylim = c(0,ymax))

	  axis(1, line = 2)  # draw x-axis lower

    }

	par(new=TRUE)

	bars1 = table(cut(z.hist,seq(x.lim.min,x.lim.max,.1)))
	bars1 = bars1 / sum(bars1)
	sum(bars1)
	names(bars1) = as.character(seq(x.lim.min+.1,x.lim.max,.1)-.1)

	bar.names <- as.numeric(sub("^\\(([^,]+),.*$", "\\1", names(bars1)))
	bar.names

	bars2 = bars1[which(bar.names == int.start):length(bars1)]
	bars2
 	sum(bars2)
	
	scale = sum(bars1)/sum(bars2)

	ymax.scale = ymax*scale

	hist(z.hist[z.hist > round(int_beg,1) & z.hist < int_end],
		breaks=n.breaks,freq=FALSE,
		col=col2,border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(0,ymax.scale),
		ylab="",xlab="",main=Title,lwd=1,axes=FALSE)


	abline(h=0)

######################################### 
######################################### 

if (Show.Text) {

#par(family = fam[2])


par(new=TRUE)
hist(c(0),main="",ylim=c(0,ymax),ylab="",xlab="",xlim=c(x.lim.min,x.lim.max),
	density=0,border="white",axes=FALSE)

	min.z = min(test_stat)
	max.z = max(test_stat)
	n.z = length(test_stat)
	n.z.sig = length(test_stat[test_stat > crit])
	n.not.shown = length(test_stat[test_stat > x.lim.max])

	### set location parameters for writing text
	y.line.factor 
	y.text = ymax
	y.line = ymax/100*y.line.factor

	### write copyright 
	i = 0
	text(x.lim.min,ymax-y.line*i,"Schimmack, Bartos, Brunner",pos=4,cex =letter.size)
	i = i + 1.5
	text(x.lim.min,ymax-y.line*i,version,pos=4,cex =letter.size)
	i = i + 1.5
	text(x.lim.min,ymax-y.line*i,date,pos=4,cex =letter.size)

	########

	#TESTING

	i = 0
	results.x = x.lim.max 

	if (!is.null(n.z.clusters)) {
		text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z.clusters,big.mark=",")),
			" clusters "),pos=2,cex = letter.size)
		i = i + 1.5
	}

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z,big.mark=",")),
		" tests "),pos=2,cex = letter.size)
	i = i + 1.5

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z.sig,big.mark=","))," significant "),
	  pos=2,cex = letter.size)
	i = i + 1.5
	text(results.x,y.text-y.line*i,paste0(n.not.shown," z > ",x.lim.max),
		pos=2,cex = letter.size)
	

#	par(family = fam[1])

	#############################################
	#############################################

	ODR = round(results$ODR[1]*100)
	ODR.low = round(results$ODR[2]*100)
	ODR.high = round(results$ODR[3]*100)
	EDR = round(results$EDR[1]*100)
	EDR.low = round(results$EDR[2]*100)
	EDR.high = round(results$EDR[3]*100)
	ERR = round(results$ERR[1]*100)
	ERR.low = round(results$ERR[2]*100)
	ERR.high = round(results$ERR[3]*100)
	FDR = results$FDR[1]*100
	FDR.low = results$FDR[2]*100
	FDR.high = results$FDR[3]*100

	ODR.a <- sprintf("%3d", round(ODR))
	ODR.b <- sprintf("%3d", round(ODR.low))
	ODR.c <- sprintf("%3d", round(ODR.high))
	EDR.a <- sprintf("%3d", round(EDR))
	EDR.b <- sprintf("%3d", round(EDR.low))
	EDR.c <- sprintf("%3d", round(EDR.high))
	ERR.a <- sprintf("%3d", round(ERR))
	ERR.b <- sprintf("%3d", round(ERR.low))
	ERR.c <- sprintf("%3d", round(ERR.high))
	FDR.a <- sprintf("%3d", round(FDR))
	FDR.b <- sprintf("%3d", round(FDR.low))
	FDR.c <- sprintf("%3d", round(FDR.high))




	if (Write.CI) {

	i = i + 3
	text(results.x,y.text-y.line*i,
		paste0(round((1-CI.ALPHA)*100),"% CI     "),
		pos=2,cex=letter.size.1)

	x.loc = results.x + c(-1.5,-1.1,-1.0,-.6,-.5,-.1,0)



	res4text = c()
	res4text = rbind(res4text,c("ODR: ","EDR: ","ERR :","FDR :"))
	res4text = rbind(res4text,c(ODR.a,EDR.a,ERR.a,FDR.a))
	res4text = rbind(res4text,rep("[",4))
	res4text = rbind(res4text,c(ODR.b,EDR.b,ERR.b,FDR.b))
	res4text = rbind(res4text,rep(",",4))
	res4text = rbind(res4text,c(ODR.c,EDR.c,ERR.c,FDR.c))
	res4text = rbind(res4text,rep("]",4))

	dim(res4text)
	res4text

	x.i = 1
	for (x.i in 1:length(x.loc)) {
		for (ii in 1:4) {
			text(x.loc[x.i],y.text-y.line*(i+2*ii),res4text[x.i,ii]
			,pos=2,cex=letter.size.1)
		}
	}

	##############################################
	# End of IF CI show results with CI
	} else {  
	# End of IF CI show results without CI

 	
	i = i + 3
	text(results.x,y.text-y.line*i,
		paste("ODR:",ODR.a,"%"),
		pos=2,cex=letter.size.1)

	i = i + 2
	text(results.x,y.text-y.line*i,
		paste("EDR:",EDR.a,"%"),
		pos=2,cex=letter.size.1)

	i = i + 2
	text(results.x,y.text-y.line*i,
		paste("ERR:",ERR.a,"%"),
		pos=2,cex=letter.size.1) 

	i = i + 2
	text(results.x,y.text-y.line*i,
		paste("FDR:",FDR.a,"%"),
		pos=2,cex=letter.size.1) 


	} # End of CI

} # End of Show.Text



} # End of Histogram


 
######################################### 
######################################### 
######################################### 


Draw.KD = function(z.draw,w,Write.CI=FALSE,cola="blue",Lwidth=5) {

	#Draw.Histogram(w.all,results=res,cola=col.hist)

	#z.draw = test_stat
	#cola = "blue"

	x.adj.max = x.lim.max + 3 * bw.draw

	### densitiy line

	summary(z.draw)

	d.all = Get.Densities(z.draw[z.draw >= x.lim.min & z.draw < x.lim.max],
		bw=bw.draw,from=x.lim.min,to=x.lim.max,augment=augment)
	summary(d.all)

	bar_width = d.all[2,1] - d.all[1,1]
	d.all = d.all[d.all[,1] > x.lim.min & d.all[,1] <= x.lim.max,]
	dim(d.all)
	sum(d.all[,2])*bar_width
	d.all.X = d.all[,1]
	d.all.Y = d.all[,2]
	summary(d.all.X)


	summary(z.draw)

	d.sig = Get.Densities(z.draw[z.draw >= int_beg & z.draw < x.lim.max],
		bw=bw.draw,from=int_beg,to=x.lim.max,augment=augment)
	
	bar_width = d.sig[2,1] - d.sig[1,1]
	dim(d.sig)
	d.sig = d.sig[d.sig[,1] > int_beg & d.sig[,1] <= x.lim.max,]

	#dim(d.sig)
	#sum(d.sig[,2])*bar_width

	d.sig.X = d.sig[,1]
	d.sig.Y = d.sig[,2]
	#summary(d.sig.X)

    #plot(d.sig.X,d.sig.Y)



	d.fit = length(z.draw[z.draw > int_beg & z.draw < x.lim.max]) /
		length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max])

	d.fit

	d.sig.Y = d.sig.Y * d.fit

	sum(d.sig[,2])*bar_width
	

	### draw the density of the observed values in the selected region
	par(new=TRUE)
	plot(d.sig.X[d.sig.X > int_beg +.05 & d.sig.X < x.lim.max - .05],
		d.sig.Y[d.sig.X > int_beg + .05 & d.sig.X < x.lim.max - .05],
		type="l",col=adjustcolor(cola, alpha.f = 0.8),lty=2,lwd=Lwidth,
		xlim =c(x.lim.min,x.lim.max),ylim=c(0,ymax),xlab="",ylab="",
		axes=FALSE)
	### draw vertical line for Beginning

	int_beg
	#height = max(d.sig.Y[which(round(d.sig.X,1) == round(int_beg+.02,1))]);height
	#segments(int_beg+.02,0,int_beg+.02,height,lty=1,lwd=3,col=cola)


} # End of Draw.KD

###################################################
#### End Draw.KD
###################################################



###################################################
#### Begin Draw Curve
###################################################

Draw.Curve.All = function(results,cola="black",Lwidth=4,
	Ltype=1,x.start=x.lim.min,x.end=x.lim.max,Show.Gamma=FALSE) {

#x.start = x.lim.min; x.end = x.lim.max;cola = "black";Lwidth=4;Ltype=1

if(length(results$ncp) > length(ncp)) {
  ncp = results$ncp[1,]
  zsds = results$zsds[1,]
  w = results$w.all[1,]
} else {
  ncp = results$ncp
  zsds = results$zsds
  w = results$w.all
}

width = 1

bar_width = .01
D.X = seq(x.lim.min,x.lim.max,bar_width) 

D.Y = build_dens(D.X, ncp, zsds, df, CURVE.TYPE,Folded=Folded) 
D.Y = D.Y / (rowSums(D.Y) * bar_width)
dim(D.Y)
D.Y = colSums(D.Y * w)


d.dense.sel = sum(D.Y[D.X > x.lim.min & D.X >= int_beg]*bar_width)
d.dense.sel

d.hist.sel = mean(as.numeric(test_stat > x.lim.min & test_stat > int_beg & test_stat < int_end)) /
	mean(as.numeric(test_stat > x.lim.min & test_stat < int_end))
d.hist.sel

check = D.Y[D.X > crit & D.X < int_beg]
length(check)
EPJS = sum(D.Y[D.X > crit & D.X < int_beg])/sum(D.Y[D.X > crit])
EPJS

OFJS = sum(as.numeric(test_stat > crit & test_stat < int_beg))
OFJS

OFS = sum(as.numeric(test_stat > crit & test_stat < int_end))
OFS

OPJS = OFJS/OFS
OPJS

scale = d.hist.sel/d.dense.sel;scale

lines(D.X[which(D.X >= x.start & D.X < x.end)],
	D.Y[which(D.X >= x.start & D.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max),ylim=c(0,ymax))


if (Show.Gamma & !is.na(results$gamma_shape) & !is.na(results$gamma_rate)) {

     

   D.Y.g = Get.Gamma.Density(D.X,shape = results$gamma_shape,rate = results$gamma_rate)
   D.Y.g <- D.Y.g / sum(D.Y.g * bar_width)

   d.dense.gamma.sel = sum(D.Y.g[D.X > x.lim.min & D.X >= int_beg]*bar_width)
   d.dense.gamma.sel

   scale.gamma = d.hist.sel/d.dense.gamma.sel;scale.gamma

   #plot(D.X,D.Y.g)

   lines(D.X[which(D.X >= x.start & D.X < x.end)],
	  D.Y.g[which(D.X >= x.start & D.X < x.end)]*scale.gamma,
	  lty=Ltype,col=col.gamma,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max),ylim=c(0,ymax))


}

boundary = int_beg
if (round(int_beg,2) == 1.96) {boundary = 2}

#segments(boundary,0,boundary,
#	D.Y[which(round(D.X,2) == round(boundary,2))[1]]*scale,
#	col=cola,lwd=3)


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		sig.level = crit
		if(two.sided == FALSE) sig.level = c(-sig.level,sig.level)

		sig.lty = rep(2,length(sig.level));sig.lty

		i = 1
		for (i in 1:length(sig.level)) {
             sig.col = "red"
             sig.lty = 1
			height = max(scale*D.Y[which(round(D.X,1) == round(sig.level[i],1))]);height
#			segments(sig.level[i]-.002,0,sig.level[i]-.002,height,lty=1,lwd=3,col="white")
			segments(sig.level[i]-.002,0,sig.level[i]-.002,height,lty=sig.lty,lwd=4,col=sig.col)
		}

} # End of Show.Significance


} # End of Draw.Curve.All


##################################################


# ============================================================
# Gamma ML Model for Estimating ERR
# Based on Brunner & Schimmack (2020) Appendix
# Simplified to z-score space (no sample sizes needed)
# ============================================================
#
# Model: NCPs ~ Gamma(shape, rate)
#        z_observed ~ N(ncp, 1), selected for z > 1.96
#
# This is the z-space equivalent of Jerry's gamma ML model.
# Jerry's version conditions on known sample sizes and estimates
# the effect size distribution. This version works directly
# with NCPs in z-space, which is what z-curve does too.
# ============================================================

#z = test_stat

# --- Main analysis function ---
run_ML_gamma <- function(
  z,
  int_beg = 1.96,
  int_end = Inf,
  crit = 1.96,
  n_starts = 1,
  start_shape = NULL,
  start_rate = NULL,
  two.sided = TRUE
) {

  
  p_extreme_obs <- mean(z >= int_end)
  p_sig_trunc_obs <- mean(z >= crit & z < int_end)
  p_all_trunc_obs <- mean(z < int_end)

  # keep only selected range
  z_obs <- z[z >= int_beg & z < int_end]

  # probability of observed absolute z falling in selected interval
  prob_select_ncp <- function(ncp, int_beg, int_end) {
    if (two.sided) {
      # P(int_beg <= |Z| < int_end | ncp)
      upper <- if (is.infinite(int_end)) 1 else pnorm(int_end, mean = ncp, sd = 1)
      lower <- pnorm(int_beg, mean = ncp, sd = 1)

      pos <- upper - lower

      upper_neg <- pnorm(-int_beg, mean = ncp, sd = 1)
      lower_neg <- if (is.infinite(int_end)) 0 else pnorm(-int_end, mean = ncp, sd = 1)

      neg <- upper_neg - lower_neg

      pos + neg
    } else {
      # one-sided positive z selection
      upper <- if (is.infinite(int_end)) 1 else pnorm(int_end, mean = ncp, sd = 1)
      lower <- pnorm(int_beg, mean = ncp, sd = 1)
      upper - lower
    }
  }

  get_ncp_upper <- function(shape, rate, cap = 50) {
    ncp_upper <- qgamma(.999999, shape = shape, rate = rate)

    if (!is.finite(ncp_upper) || ncp_upper <= 0) {
      ncp_upper <- cap
    }

    min(ncp_upper, cap)
  }

  # Probability of observing a z-value in any interval [lo, hi),
  # integrated over the gamma distribution of ncp
  prob_select_range <- function(shape, rate, lo, hi) {

    if (lo <= 0 && is.infinite(hi)) {
      return(1)
    }

    ncp_upper <- get_ncp_upper(shape, rate)

    integrate(
      function(ncp) {
        prob_select_ncp(
          ncp = ncp,
          int_beg = lo,
          int_end = hi
        ) *
          dgamma(ncp, shape = shape, rate = rate)
      },
      lower = 0,
      upper = ncp_upper,
      rel.tol = 1e-8,
      subdivisions = 1000
    )$value
  }


  # Probability of observing a z-value in the fitted interval
  # [int_beg, int_end). This is the likelihood denominator.
  prob_select <- function(shape, rate) {
    prob_select_range(
      shape = shape,
      rate = rate,
      lo = int_beg,
      hi = int_end
    )
  }

  compute_edr_gamma_trunc <- function(shape, rate, zmax = int_end) {
    p_sig <- prob_interval_gamma(shape, rate, crit, zmax)
    p_all <- prob_interval_gamma(shape, rate, 0, zmax)
    p_sig / p_all
  }

  get_ncp_upper <- function(shape, rate, prob = .999999, cap = 50) {
    upper <- qgamma(prob, shape = shape, rate = rate)

    if (!is.finite(upper) || upper <= 0) {
      upper <- cap
    }

    min(upper, cap)
  }

  marginal_density <- function(zj, shape, rate) {

    ncp_upper <- get_ncp_upper(shape, rate)

    integrate(
      function(ncp) {
        if (two.sided) {
          dens <- dnorm(zj, mean = ncp, sd = 1) +
            dnorm(-zj, mean = ncp, sd = 1)
        } else {
          dens <- dnorm(zj, mean = ncp, sd = 1)
        }

        dens * dgamma(ncp, shape = shape, rate = rate)
      },
      lower = 0,
      upper = ncp_upper,
      rel.tol = 1e-8,
      subdivisions = 1000
    )$value
  }

  negloglik_gamma <- function(params) {
    shape <- exp(params[1])
    rate  <- exp(params[2])

    # Optional but useful to prevent degenerate optimizer escapes
    if (!is.finite(shape) || !is.finite(rate) ||
        shape <= 0 || rate <= 0 ||
        shape > 100 || rate > 1000) {
      return(1e10)
    }

    p_sel <- tryCatch(
      prob_select(shape, rate),
      error = function(e) NA_real_
    )

    if (!is.finite(p_sel) || p_sel <= 1e-15) {
      return(1e10)
    }

    log_densities <- vapply(
      z_obs,
      function(zj) {
        md <- tryCatch(
          marginal_density(zj, shape, rate),
          error = function(e) NA_real_
        )

        if (!is.finite(md) || md <= 0) {
          return(-1e10)
        }

        log(md)
      },
      numeric(1)
    )

    nll <- -(sum(log_densities) - length(z_obs) * log(p_sel))

    if (!is.finite(nll)) 1e10 else nll
  }


  best_fit <- NULL
  best_nll <- Inf

  start_list <- list()

  if (!is.null(start_shape) && !is.null(start_rate)) {
    start_list[[1]] <- log(c(start_shape, start_rate))
  }

  for (i in seq_len(n_starts)) {
    start_list[[length(start_list) + 1]] <- log(c(
      runif(1, 0.2, 10),
      runif(1, 0.05, 10)
    ))
  }

  for (start_params in start_list) {
    fit <- tryCatch(
      optim(
        par = start_params,
        fn = negloglik_gamma,
        method = "Nelder-Mead",
        hessian = TRUE,
        control = list(maxit = 5000)
      ),
      error = function(e) NULL
    )

    if (!is.null(fit) && fit$value < best_nll) {
      best_fit <- fit
      best_nll <- fit$value
    }
  }

  if (!is.null(best_fit)) {

    shape <- exp(best_fit$par[1])
    rate  <- exp(best_fit$par[2])

    # --- EDR and ERR on the same truncated z-range as the fitted density ---

    compute_edr_gamma <- function(shape, rate) {

      # Numerator: P(crit <= |Z| < int_end)
      p_sig <- prob_select_range(
        shape = shape,
        rate = rate,
        lo = crit,
        hi = int_end
      )

      # Denominator: P(0 <= |Z| < int_end)
      p_all <- prob_select_range(
        shape = shape,
        rate = rate,
        lo = 0,
        hi = int_end
      )

      p_sig / p_all
    }


    compute_err_gamma <- function(shape, rate) {

      ncp_upper <- get_ncp_upper(shape, rate)

      numerator <- integrate(
        function(ncp) {
          p_sig_ncp <- prob_select_ncp(
            ncp = ncp,
            int_beg = crit,
            int_end = int_end
          )

          p_sig_ncp^2 * dgamma(ncp, shape = shape, rate = rate)
        },
        lower = 0,
        upper = ncp_upper,
        rel.tol = 1e-8,
      subdivisions = 1000
      )$value

      denominator <- prob_select_range(
        shape = shape,
        rate = rate,
        lo = crit,
        hi = int_end
      )

      numerator / denominator
    }

    edr_trunc <- compute_edr_gamma(shape, rate)
    err_trunc <- compute_err_gamma(shape, rate)

    edr <- p_all_trunc_obs * edr_trunc + p_extreme_obs
    err <- (p_sig_trunc_obs * err_trunc + p_extreme_obs) /
            (p_sig_trunc_obs + p_extreme_obs)

  # --- Delta-method SEs for ERR and EDR ---
  # Parameters are optimized on the log scale:
  # theta = (log(shape), log(rate)).
  # The Hessian returned by optim is therefore the observed information matrix
  # for theta. SEs for functions of shape and rate are obtained by the delta method.

    eps_shape <- max(1e-5, abs(shape) * 1e-4)
    eps_rate  <- max(1e-5, abs(rate)  * 1e-4)

    H_inv <- solve(best_fit$hessian)

    # ERR
    err0 <- compute_err_gamma(shape, rate)

    derr_dshape <- (compute_err_gamma(shape + eps_shape, rate) - err0) / eps_shape
    derr_drate  <- (compute_err_gamma(shape, rate + eps_rate) - err0) / eps_rate

    # Chain rule because the Hessian is for log(shape), log(rate):
    # d f / d log(shape) = d f / d shape * shape
    # d f / d log(rate)  = d f / d rate  * rate
    grad_err <- c(
      derr_dshape * shape,
      derr_drate  * rate
    )

    se_err <- sqrt(as.numeric(t(grad_err) %*% H_inv %*% grad_err))

    # EDR
    edr0 <- compute_edr_gamma(shape, rate)

    dedr_dshape <- (compute_edr_gamma(shape + eps_shape, rate) - edr0) / eps_shape
    dedr_drate  <- (compute_edr_gamma(shape, rate + eps_rate) - edr0) / eps_rate

    grad_edr <- c(
      dedr_dshape * shape,
      dedr_drate  * rate
    )

    se_edr <- sqrt(as.numeric(t(grad_edr) %*% H_inv %*% grad_edr))

    gamma_mean_ncp = shape/rate
    gamma_sd_ncp = sqrt(shape)/ rate

    se_edr[is.nan(se_edr)] <- NA
    se_err[is.nan(se_err)] <- NA

   } else { # EOF shape and rate were available 

     gamma_mean_ncp = NA
     gamma_sd_ncp = NA
     gamma_edr = c(NA,NA)
     gamma_err = c(NA,NA)

   }
    
   gamma_res = list(
      gamma_shape = shape,
      gamma_rate = rate,
      gamma_mean_ncp = gamma_mean_ncp,
      gamma_sd_ncp = gamma_sd_ncp,
      gamma_edr = c(edr,se_edr),
      gamma_err = c(err,se_err),
      negloglik = best_nll,
      convergence = best_fit$convergence
    )

    gamma_res

}

#################

# Gamma-implied density of absolute z-values
#z = seq(0,6,.01); shape = results$gamma_shape;rate = results$gamma_rate
#delta_max= NULL; scale=NULL;x.end = x.lim.max;x.start = x.lim.min

Get.Gamma.Density <- function(z, shape, rate = NULL, 
  scale = NULL, delta_max = NULL) {

  if (is.null(rate) && is.null(scale)) {
    stop("Provide either rate or scale.")
  }
  if (is.null(rate)) {
    rate <- 1 / scale
  }

  if (is.null(delta_max)) {
    delta_max <- qgamma(.99999, shape = shape, rate = rate)
  }

  if (shape < 0.1) {
    # H0 for all studies: predicted density is truncated N(0,1)
    gamma.density <- dnorm(z) * 2
  } else {
    gamma.density = sapply(z, function(zz) {
      integrate(
        function(delta) {
          # density of |Z| at zz given delta
          folded_density <- dnorm(zz, mean = delta, sd = 1) +
            dnorm(-zz, mean = delta, sd = 1)

          folded_density * dgamma(delta, shape = shape, rate = rate)
        },
        lower = 0,
        upper = delta_max,
        rel.tol = 1e-8
      )$value
    })

  }

    return(gamma.density)

}


#################################################
#################################################
#################################################

#For quick testing
#test_stat = abs(c(rnorm(500,0,1),rnorm(100,2,1)));hist(test_stat)
#cluster.id = NULL;yi=NULL;sei=NULL


#################################################
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
#################################################


  ### --- curve type / crit derivation ---
  curve_type <- "z"
  # crit already set via the formal-parameter default (qnorm(1 - alpha/2)) or user override

  ### --- determine test_stat from whichever input was supplied ---
  n_inputs_given <- sum(!is.null(zval), !is.null(tval), !is.null(pval), !is.null(yi))
  if (n_inputs_given == 0) stop("Supply one of zval, tval (+ df), pval, or yi (+ sei).")
  if (n_inputs_given > 1)  stop("Supply exactly one of zval, tval (+ df), pval, or yi (+ sei) -- not more than one.")

  if (!is.null(zval)) test_stat <- zval
  if (!is.null(pval)) test_stat <- qnorm(1 - pval / 2)
  if (!is.null(yi))   test_stat <- yi / sei

  if (!is.null(tval)) {
    if (is.null(df)) stop("df must be supplied when tval is used.")
    curve_type <- "t"
    if (missing(crit)) crit <- qt(alpha / 2, df, lower.tail = FALSE)
    test_stat <- tval
  }

  if (curve_type == "t" && est_method != "OF") {
    stop("Genuine t-curve fitting (tval + df) is only supported for est_method = \"OF\".")
  }

  ### --- remove NA (one mask, built from test_stat, not from zval/tval specifically) ---
  sel <- !is.na(test_stat)
  if (!is.null(cluster_id)) sel <- sel & !is.na(cluster_id)
  if (!is.null(yi))         sel <- sel & !is.na(yi)
  if (!is.null(sei))        sel <- sel & !is.na(sei)

  test_stat  <- test_stat[sel]
  if (!is.null(cluster_id)) cluster_id <- cluster_id[sel]
  if (!is.null(yi))         yi         <- yi[sel]
  if (!is.null(sei))        sei        <- sei[sel]

  ### --- fold to |test_stat| when sign doesn't matter ---
  if (!directional) {
    test_stat <- abs(test_stat)
    if (!is.null(yi)) yi <- abs(yi)
  }

  n_tests <- length(test_stat)
  n_sig   <- sum(test_stat > crit)
  if (n_sig < 20) stop("Insufficient data: need at least 20 significant test statistics.")

  ### --- extreme-value diagnostic ---
  extreme <- c(
    Ext.Neg = sum(test_stat < -int_end) / max(1, sum(test_stat < -crit)),
    Ext.Pos = sum(test_stat >  int_end) / max(1, sum(test_stat >  crit))
  )

  n_clusters <- NULL
  if (!is.null(cluster_id)) n_clusters <- length(unique(cluster_id))

  ### --- density-fit diagnostic (informational only -- not fed into any adjustment) ---
  slope_width <- 1
  slope_int   <- test_stat[test_stat >= int_beg + 2 * control$bw_est &
                            test_stat <  int_beg + 2 * control$bw_est + slope_width]
  slope_k <- length(slope_int)
  slope   <- rep(NA, 6)
  if (slope_k > 10) {
    s <- get_slope(int = slope_int, bw = control$bw_est, 
         from = int_beg, to = int_beg + slope_width)
    reliability    <- s["slope"]^2 / (s["slope"]^2 + s["se"]^2)
    adjusted_slope <- reliability * s["slope"]
    slope <- c(slope_k, adjusted_slope, s)
  }
  names(slope) <- c("slope k", "slope adj", "slope", "slope se", "slope lb", "slope ub")

  ### --- fit ---
  res_main <- run_zcurve(
    est_method   = est_method,
    test_stat    = test_stat,
    crit         = crit,
    int_beg      = int_beg,
    int_end      = int_end,
    cluster_id   = cluster_id,
    boot.iter    = control$boot_iter,
    ncp          = ncp,
    z_sd         = z_sd,
    ncp_fixed    = ncp_fixed,
    z_sd_fixed   = z_sd_fixed,
    w_fixed      = w_fixed,
    augment      = augment,
    directional  = directional,   # NOTE: run_zcurve() doesn't accept this yet -- see file header
    folded       = folded,
    alpha        = alpha,         # NOTE: same
    int_loc      = int_loc,       # NOTE: same
    df           = df,            # NOTE: same
    yi           = yi,
    sei          = sei,
    x_lim_min    = x_lim_min,
    x_lim_max    = x_lim_max,
    spline_width = control$spline_width
  )

  ### --- assemble output ---
  FDR <- as.numeric((1 / res_main$EDR - 1) * (alpha / (1 - alpha)))
  FDR <- FDR[c(1, 3, 2)]
  names(FDR) <- c("FDR_pe", "FDR_lb", "FDR_ub")

  POW_h1_sig <- as.numeric((res_main$ERR - FDR * alpha) / (1 - FDR))
  names(POW_h1_sig) <- c("POW_H1_SIG_pe", "POW_H1_SIG_lb", "POW_H1_SIG_ub")

  ODR_EDR_D <- as.numeric(res_main$ODR_EDR_D)
  names(ODR_EDR_D) <- c("ODR_EDR_D_pe", "ODR_EDR_D_lb", "ODR_EDR_D_ub")

  res_main$ODR <- as.numeric(res_main$ODR); names(res_main$ODR) <- c("ODR_pe", "ODR_lb", "ODR_ub")
  res_main$EDR <- as.numeric(res_main$EDR); names(res_main$EDR) <- c("EDR_pe", "EDR_lb", "EDR_ub")
  res_main$ERR <- as.numeric(res_main$ERR); names(res_main$ERR) <- c("ERR_pe", "ERR_lb", "ERR_ub")

#  names(res_main$shape_d_mean)   <- c("shape_mean_pe", "shape_mean_lb", "shape_mean_ub")
#  names(res_main$shape_d_median) <- c("shape_median_pe", "shape_median_lb", "shape_median_ub")

  names(res_main$es_mean_all)   <- c("es_mean_all_pe", "es_mean_all_lb", "es_mean_all_ub")
  names(res_main$es_median_all) <- c("es_median_all_pe", "es_median_all_lb", "es_median_all_ub")
  names(res_main$es_mean_sig)   <- c("es_mean_sig_pe", "es_mean_sig_lb", "es_mean_sig_ub")
  names(res_main$es_tau)        <- c("es_tau_pe", "es_tau_lb", "es_tau_ub")

  if (length(res.main$local_power) > k_ncp) {
    colnames(res_main$local_power) <- paste0("local_power", seq_len(ncol(res_main$local_power)))
    rownames(res_main$local_power) <- c("pe", "lb", "ub")
  }

  if (length(res_main$local_es) > k_ncp) {
    colnames(res_main$local_es) <- paste0("local_es", seq_len(ncol(res_main$local_es)))
    rownames(res_main$local_es) <- c("pe", "lb", "ub")
  }

  if (k_ncp > 1) {
    colnames(res_main$ncp)   <- paste0("ncp",   seq_len(k_ncp)); rownames(res_main$ncp)   <- c("pe", "lb", "ub")
    colnames(res_main$zsds)  <- paste0("zsds",  seq_len(k_ncp)); rownames(res_main$zsds)  <- c("pe", "lb", "ub")
    colnames(res_main$w.inp) <- paste0("w.inp", seq_len(k_ncp)); rownames(res_main$w.inp) <- c("pe", "lb", "ub")
    colnames(res_main$w.all) <- paste0("w.all", seq_len(k_ncp)); rownames(res_main$w.all) <- c("pe", "lb", "ub")
  }

  results <- list(
    clusters       = n_clusters,
    slope          = slope,
    extreme        = extreme,
    ODR            = res_main$ODR,
    EDR            = res_main$EDR,
    ERR            = res_main$ERR,
    FDR            = FDR,
    POW_h1_sig     = POW_h1_sig,
    ncp            = res_main$ncp,
    zsds           = res_main$zsds,
    w.inp          = res_main$w.inp,
    w.all          = res_main$w.all,
    local_power    = res_main$local_power,
    local_es       = res_main$local_es,
    ODR_EDR_D      = ODR_EDR_D,
    shape_d_median = res_main$shape_d_median,
    shape_d_mean   = res_main$shape_d_mean,
    fit            = res_main$fit,
    gamma_shape    = res_main$gamma_shape,
    gamma_rate     = res_main$gamma_rate,
    gamma_edr      = res_main$gamma_edr,
    gamma_err      = res_main$gamma_err,
    es_mean_all    = res_main$es_mean_all,
    es_median_all  = res_main$es_median_all,
    es_mean_sig    = res_main$es_mean_sig,
    es_tau         = res_main$es_tau

#   version        = if (is.null(version) || version == "") NULL else version,
#   date           = if (is.null(date)    || date    == "") NULL else date
  )

  class(results) <- "zcurve"

  show_summary = TRUE
  show_plot = TRUE

  if (isTRUE(show_summary)) print(results)
# if (isTRUE(show_plot))    plot(results)

  invisible(results)

} ### EOF zcurve


######################################################################
######################################################################
######################################################################

yi = rnorm(1000,.4,.2)
sei = rep(.2,length(yi))
test_stat = yi / sei
zval = test_stat
summary(zval)
sd(zval)

fit = zcurve(zval)






boot.iter = 0,
parallel = TRUE,
em_criterion = 1e-6,
em_max_iter  = 2000,
em_max_iter_boot = 500,
   ci_alpha = .05,
   augment = TRUE,
   n_bars = 512,
   bw_est = .05,
   spline_width = .5
)
   

tval= NULL
yi = NULL
sei = NULL
cluster_id = NULL
pval = NULL

directional = TRUE
folded = FALSE
test_stat 
est_method = "OF"
int_beg = 0
int_end = 6
alpha = .05
crit = qnorm(1-alpha/2)
ncp_fixed = TRUE
z_sd_fixed = TRUE
w_fixed = FALSE
augment = TRUE
bw_est = .05
ncp = 0:6
k_ncp = length(ncp)
z_sd = rep(1,k_ncp)
curve_type = "z"
folded = FALSE
x_lim_min = 0
x_lim_max = 6
ymax = 0.6


from = int_beg
to = int_beg + 1
bw = control$bw_est
int_loc = 1
grid = seq(0,6,.01)

boot_iter = 100