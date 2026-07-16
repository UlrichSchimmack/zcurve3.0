

#rm(list = ls())

########################################################################
### SETTING PARAMETERS FOR Z-CURVE MODEL
########################################################################


version <- "zcurve3(3.90)"
date <- "2026.07.15"  # Version label to appear on plots
note <- "Effect Size Estimation with Bayesian Shrinkage"

print(version)
print(date)
print(note)

# Optional cleanup 
# rm(list = ls())
# options(scipen = 999)  # Disable scientific notation


### INSTALL PACKAGES (only once – manually run if needed)
if (1 == 2) {  # This block is ignored unless manually changed to (1 == 1)
  install.packages("KernSmooth")
  install.packages("parallel")
} # END install block

### LOAD LIBRARIES
library(parallel)
library(KernSmooth)

########################################################################
### GLOBAL PARAMETERS
########################################################################

TESTING <- FALSE             # Toggle for development/debugging mode

CURVE.TYPE = "z"             # Set to "t" for t-distribution
df = c()

Directional = FALSE          # Assumes sign does not matter, set to TRUE for meta-analysis with directional z-values
Folded      = TRUE           # Assumes absolute z-values as input

### SPEED CONTROL

parallel <- TRUE             # Default - serial option not available
max_iter <- 1e6              # Max iterations for model estimation
max_iter_boot <- 1e5         # Max iterations for bootstrapped estimates

EM.criterion <- 1e-3         # Convergence threshold for EM algorithm
EM.max.iter <- 1000          # Max iterations for EM

Plot.Fitting <- FALSE        # Plot fitting curve (only for Est.Method = "OF" or "EXT")

### PLOT SETTINGS

Title <- ""                  # Optional plot title

letter.size <- 1             # Text size in plots
letter.size.1 <- letter.size # Used for version labels in plot
y.line.factor <- 3           # Controls spacing of plot text

x.lim.min <- 0               # X-axis lower bound
x.lim.max <- 6               # X-axis upper bound
ymax <- 0.6                  # Y-axis upper bound

Show.Histogram <- TRUE       # Toggle histogram in plot
Show.Text <- TRUE            # Toggle model results in plot
Show.Curve.All <- TRUE       # Show predicted z-curve
Show.Curve.Sig <- FALSE      # Option: show z-curve only for significant values
Show.Significance <- TRUE    # Show z = critical value line
Show.Gamma = FALSE           # Show gamma distribtuion based on ML_gamma estimates
Show.KD <- FALSE             # Toggle kernel density overlay (density method only)
Show.Y.Label <- TRUE         # Option to suppress to add manually

sig.levels <- c()            # Optional: mark additional p-value thresholds on plot

int.loc <- 1                 # Plot local power intervals below x-axis (set 0 to disable)
hist.bar.width <- 0.2        # Width of histogram bars

Heterogeneity.Test = FALSE   # run z-curve with single normal component to get SD estimate
min.heterogeneity = .5       # criterion for use of homogenous model 

### MODEL PARAMETERS

alpha <- 0.05                        # Significance level
crit <- qnorm(1 - alpha / 2)         # Corresponding two-sided critical z

two.sided <- TRUE                   # Assume two-sided z-values (use abs(z)); not yet compatible with signed z-values

# Color scheme
col.curve <- "red"
col.hist  <- "blue3"
col.kd    <- "black"
col.gamma <- "chartreuse3"

Est.Method <- "OF"                  # Estimation method: "OF", "EM", "ML_gamma" 
Run.Gamma  <- FALSE

Int.Beg <- 1.96                     # Default: critical value for alpha = .05
Int.End <- 6                        # End of modeling interval (z > 6 = power = 1)

ncp <- 0:6                          # Component locations (z-values at which densities are centered)
components = length(ncp)
zsds = rep(1,components)            # one SD for each component

just <- 0.8                         # Cutoff for "just significant" z-values (used in optional bias test)

ZSDS_FIXED <- TRUE                  # Fix SD values for EXT method 
NCP_FIXED <- TRUE                   # Fix non-central parameter(NCP) means values for EXT method
W_FIXED   <- FALSE                  # Fix weights for EXT method

### DENSITY-BASED SETTINGS (Only used with Est.Method = "OF")

n.bars <- 512                       # Number of bars in histogram

Augment <- TRUE                     # Apply correction for bias at lower bound

bw.est  <- 0.05                     # Bandwidth for kernel density (lower = less smoothing, higher = more smoothing)
bw.draw <- 0.05   		        	 # Bandwith of kernel density in plot
spline.width <- .5                  # blending of truncated and normal KD

cluster.id = NULL
### INPUT RESTRICTIONS

MAX.INP.Z <- 100                    # Optionally restrict very large z-values (set Inf to disable)

### CONFIDENCE INTERVALS / BOOTSTRAPS

boot.iter <- 0                      # Number of bootstrap iterations (suggest 500+ for final models)

CI.ALPHA <- 0.05                    # CI level (default = 95%)
CI.EDR.MIN.ADJ = .05                # Minimum Adjustmnet for EDR CI
CI.ERR.MIN.ADJ = .03                # Minimum Adjustment for ERR CI


##################################################################
### END OF SETTING DEFAULT PARAMETERs
##################################################################


Zing = function(val.input,cluster.id=c(),df=c(),p=c(), lp=c(),
  yi=NULL, sei=NULL )  {

##################################################################
##################################################################
### AAA Start of Zing Code
##################################################################
##################################################################


run_SQP <- function(
  val.input,       
  ncp,                       # grid of component means (z-space)
  zsds,                      # grid of component SDs (z-space)
  Int.Beg,
  Int.End,
  Folded  = TRUE
) {

  INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
  #hist(INT)
  ODR = mean(val.input > 1.96);ODR

  ncp_k <- length(ncp)

  comp_sd <- sqrt(zsds^2 + 1)          # component sd convolved with unit sampling sd

  # ---- density of each z under each component, in the fitted window -----
  dens <- function(zz, mu, sdv) {
    if (Folded) dnorm(zz, mu, sdv) + dnorm(zz, -mu, sdv)
    else        dnorm(zz, mu, sdv)
  }

  L <- outer(seq_along(INT), seq_along(ncp),
             function(i, ncp_k) dens(INT[i], ncp[ncp_k], comp_sd[ncp_k]))

  # ---- selection prob = mass of component k INSIDE [Int.Beg, Int.End] ----
  #      (double truncation, not just significance)
  psel <- vapply(seq_along(ncp), function(k) {
    if (Folded) {
      # Folded model: fitted values are |z|, so the window [Int.Beg, Int.End]
      # applies to magnitude. Mass = both signed tails mapping to that |z| range.
      hi <- pnorm( Int.End, ncp[k], comp_sd[k]) - pnorm( Int.Beg, ncp[k], comp_sd[k])
      lo <- pnorm(-Int.Beg, ncp[k], comp_sd[k]) - pnorm(-Int.End, ncp[k], comp_sd[k])
      hi + lo
    } else {
      # Signed model: fitted values are signed z in the window [Int.Beg, Int.End].
      # Mass is simply the probability of landing in that signed interval.
      # For unselected data, set Int.Beg = -Int.End so this is ~1 (no truncation).
      pnorm(Int.End, ncp[k], comp_sd[k]) - pnorm(Int.Beg, ncp[k], comp_sd[k])
    }
  }, numeric(1))
  psel <- pmax(psel, .Machine$double.xmin)

  # renormalize columns so recovered weights are UNSELECTED population g
  L <- sweep(L, 2, pmax(psel, 1e-12), "/")

  # ---- convex solve: global optimum, no multistart ----------------------
  w.inp <- mixsqp::mixsqp(L, log = FALSE, control = list(verbose = FALSE))$x
  w.inp / sum(w.inp)
  names(w.inp) = ncp

  hist(INT)

  sqp.est = list(
    w.inp  = w.inp,
    ncp    = ncp,
    zsds   = zsds
  )

 
  cp.inp = sqp.est
  cp.inp$val.input = val.input

  ## ---- Compute power for each bootstrap result ----
  cp.res = Compute.Power.Z.General(
        cp.inp = cp.inp, 
        Int.Beg = Int.Beg,
        Int.End = Int.End, 
        yi = yi, 
        sei = sei,
        x.lim.min = x.lim.min,
        x.lim.max = x.lim.max)
  cp.res

  ret = list(
    w.inp           = sqp.est$w.inp,
    ncp             = sqp.est$ncp,
    zsds            = sqp.est$zsds,
    model.fit       = sqp.est$fit,
    ODR             = ODR,   
    EDR             = cp.res$EDR,
    ERR             = cp.res$ERR,
    w.all           = cp.res$w.all, 
    local_power     = cp.res$local_power,
    shape_d_median  = cp.res$shape_d_median,
    shape_d_mean    = cp.res$shape_d_mean,
    local_es        = cp.res$local_es,
    es_mean_all     = cp.res$es_mean_all,
    es_median_all   = cp.res$es_median_all,
    es_mean_sig     = cp.res$es_mean_sig,
    es_tau          = cp.res$es_tau
  )

  
  #print(cp.res$mean_all)
  #print(ret$es_mean_all)

  return(ret)
}



##############################################
### Get Densities
##############################################

#D.X = seq(Int.Beg,Int.End,.01)
build.dens <- function(D.X, ncp, zsds = NULL, df = NULL, 
                       CURVE.TYPE = "z", Folded = TRUE) {
  
  xmat  <- matrix(rep(D.X, each = length(ncp)), nrow = length(ncp))
  mumat <- matrix(rep(ncp, times = length(D.X)), nrow = length(ncp))
  
  if (CURVE.TYPE == "z") {
    
    if (is.null(zsds)) {
      zsds <- rep(1, length(ncp))
    }
    
    if (length(zsds) == 1) {
      zsds <- rep(zsds, length(ncp))
    }
    
    sdmat <- matrix(rep(zsds, times = length(D.X)), nrow = length(ncp))
    

    if (Folded) {
      
      # Standard z-curve: folded normal density.
      Dens <- dnorm(xmat, mean = mumat, sd = sdmat) +
              dnorm(-xmat, mean = mumat, sd = sdmat)
      
    } else {

      # mass of component k inside [Int.Beg, Int.End]
      denom <- pnorm(Int.End, mean = mumat, sd = sdmat) -
               pnorm(Int.Beg, mean = mumat, sd = sdmat)
      denom <- pmax(denom, .Machine$double.xmin)

      Dens <- dnorm(xmat, mean = mumat, sd = sdmat) / denom

     
    }
    
  } else {
    
    if (is.null(df)) {
      stop("df must be supplied when CURVE.TYPE is not 'z'.")
    }
    
    if (length(df) == 1) {
      df <- rep(df, length(ncp))
    }
    
    dfmat <- matrix(rep(df, times = length(D.X)), nrow = length(ncp))
    
    if (Folded) {
      
      # Standard z-curve: folded noncentral t density.
      Dens <- dt(xmat, df = dfmat, ncp = mumat) +
              dt(-xmat, df = dfmat, ncp = mumat)
      
    } else {
      
      # Directional model: signed positive noncentral t density,
      # truncated at zero.
      denom <- pt(0, df = dfmat, ncp = mumat, lower.tail = FALSE)
      denom <- pmax(denom, .Machine$double.xmin)
      
      Dens <- dt(xmat, df = dfmat, ncp = mumat) / denom
      
      # Only relevant if D.X accidentally contains negative values.
      Dens[xmat < 0] <- 0
    }
  }

  #plot(D.X,Dens[4,])  
    
  Dens

}

###

#ddd

Get.Densities = function(INT, bw = 0.20, from = 1.96, to = 6, width = 1, Augment = TRUE) {

  #from = -8; to = 8;bw = .05
  #
  #summary(INT)

  n = length(seq(from,to,.02))
  grid <- seq(from, to, length.out = n)
  dens <- numeric(n)

  if (Augment) {

    bnd.width <- width

    splice <- from + bnd.width
    bnd.idx <- which(grid < splice)
    int.idx <- which(grid >= splice)
    
    ## ---- Truncated normal kernel for boundary zone ----
    INT.s <- INT - from
    for (i in bnd.idx) {
      z.s <- grid[i] - from
      raw <- dnorm(INT.s, z.s, bw)
      trunc.corr <- pnorm(z.s / bw)
      dens[i] <- mean(raw) / trunc.corr
    }
 
    ## ---- bkde for interior ----
    bkde.out <- bkde(INT, bandwidth = bw, range.x = c(from, to), gridsize = n)
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

    bkde.out <- bkde(INT, bandwidth = bw, range.x = c(from, to), gridsize = n)
    dens <- bkde.out$y

  }
  
  D <- data.frame(ZX = grid, ZY = dens)
  bar.width <- D$ZX[2] - D$ZX[1]
  D$ZY <- D$ZY / (sum(D$ZY * bar.width))

  #plot(D$ZX,D$ZY)

  return(D)

} # EOF Densities 

#######################################################
### End of Get Densities
#######################################################




#####################################################################
### Compute Power Function (New: Discrete & Continuous
#####################################################################

#ppp

#####################################################################
### Compute Power Function (Discrete & Continuous)
#####################################################################
#
# Fixes applied vs. the prior version:
#   * median_all now returns the weighted MEDIAN (was mistakenly the mean).
#   * ncp passed as cp.inp$ncp in the es call (was an undefined `ncp`).
#   * int.loc passed by name (was positional -> landed in wrong slot).
#   * out$median_all is now propagated from the es result (was dropped).
#   * All enclosing-scope dependencies made EXPLICIT arguments:
#       - get_zcurve_es_estimates now receives:
#           val.input, yi, Directional, xz
#       - estimate_sei_function default z_min = Int.Beg (was x.lim.max, a bug).
#   * Directional/Folded no longer mixed: the es function uses Directional
#     consistently for the abs() folding of extreme means and es.X.
#   * Extremes block: length invariant checked against z_grid, es.X, wX.
#   * Continuous branch left in but self-contained args noted; still gated off.
#####################################################################


Compute.Power.Z.General = function(
    cp.inp, val.input, yi, sei, Int.Beg, Int.End,
    x.lim.min, x.lim.max) {

   # crit, Directional, Folded, alpha,
   #  int.loc) {


  #-------------------------------------------------------------------
  # helper: estimate a sampling-error function sei(z) over the z grid
  #-------------------------------------------------------------------
  estimate_sei_function <- function(
      z,
      sei,
      z_grid,
      crit    = 1.96,
      z_max   = max(z_grid, na.rm = TRUE),
      z_min   = Int.Beg,          # FIX: was x.lim.max (wrong threshold)
      Int.Beg = 0,    
      lambda  = 0.5,
      min_k   = 10) {

    dat <- data.frame(z = z, sei = sei)
    dat <- dat[is.finite(dat$z) & is.finite(dat$sei) & dat$sei > 0, ]

    # Fit over the full available positive z range
    dat_full <- dat[dat$z >= 0 & dat$z <= z_max, ]

    if (nrow(dat_full) < min_k) {
      return(data.frame(z_grid = z_grid,
                        sei_out = rep(NA_real_, length(z_grid))))
    }

    fit_full <- lm(log(sei) ~ z, data = dat_full)
    sei_full <- exp(predict(fit_full, newdata = data.frame(z = z_grid)))

    # Default: full-range SE function
    sei_out <- sei_full

    # Optional shrinkage extrapolation below z_min when selecting (Int.Beg > 0)
    if (Int.Beg > 0) {
      dat_sig <- dat[dat$z >= Int.Beg & dat$z <= z_max, ]
      if (nrow(dat_sig) >= min_k) {
        fit_sig  <- lm(log(sei) ~ z, data = dat_sig)
        b1       <- coef(fit_sig)[["z"]]
        sei_crit <- exp(predict(fit_sig, newdata = data.frame(z = crit)))
        sei_low  <- sei_crit * exp(lambda * b1 * (z_grid - crit))

        sei_out[z_grid <  z_min] <- sei_low[z_grid <  z_min]
        sei_out[z_grid >= z_min] <- sei_full[z_grid >= z_min]
      }
    }

    data.frame(z_grid = z_grid, sei_out = sei_out)
  }



  #-------------------------------------------------------------------
  # helper: empirical-Bayes effect-size estimates
  #-------------------------------------------------------------------

  get_zcurve_es_estimates <- function(
      loc.p,
      sei_est,
      z_grid,
      wd.X,
      ncp,
      w.all,
      val.input,
      yi,
      Directional,
      crit    = 1.96,
      Int.Beg = 0,
      Int.End = 6,
      int.loc = 0.2,
      Folded  = FALSE,
      xz      = 0.01) {

    comp_ncp <- ncp
    comp_w   <- w.all

    z_grid  <- as.numeric(sei_est[, 1])
    wX      <- as.numeric(wd.X)
    sei_out <- as.numeric(sei_est[, 2])

    ok      <- is.finite(z_grid) & is.finite(wX) & is.finite(sei_out)
    stopifnot(length(z_grid) == length(wX), length(z_grid) == length(sei_out))
    z_grid  <- z_grid[ok]
    wX      <- wX[ok]
    sei_out <- sei_out[ok]

    # normalize reconstructed density to probability weights
    wX <- wX / sum(wX, na.rm = TRUE)

    # ---- Empirical-Bayes posterior-mean ncp given observed z ---------------
    #   ncp_hat(z) = sum_k w_k phi(z - mu_k) mu_k / sum_k w_k phi(z - mu_k)
    if (!Folded) {
      Lik <- outer(z_grid, comp_ncp, function(x, mu) dnorm(x - mu))
    } else {
      Lik <- outer(z_grid, comp_ncp, function(x, mu) dnorm(x - mu) + dnorm(x + mu))
    }

    num <- as.numeric(Lik %*% (comp_w * comp_ncp))
    den <- as.numeric(Lik %*% comp_w)

    local_ncp.X <- num / den
    local_ncp.X[!is.finite(local_ncp.X)] <- NA_real_

    if (Directional) {
      es.X <- local_ncp.X * sei_out
    } else {
      es.X <- abs(local_ncp.X) * sei_out
    }

    # ---- heterogeneity (tau) of the fitted mixture -------------------------
    w_g      <- w.all / sum(w.all)
    ncp_mean <- sum(w_g * ncp)
    ncp_tau2 <- sum(w_g * ncp^2) - ncp_mean^2
    ncp_tau  <- sqrt(ncp_tau2)
    es_tau   <- ncp_tau * mean(sei_out)

    # ---- Extreme-value (tail) correction -----------------------------------
    if (Directional) {
      es.ext.neg <- mean(yi[val.input < -Int.End])
    } else {
      es.ext.neg <- abs(mean(yi[val.input < -Int.End]))
    }
    w.ext.neg <- mean(val.input < -Int.End)

    es.ext.pos <- mean(yi[val.input > Int.End])
    w.ext.pos  <- mean(val.input > Int.End)

    grid_share <- 1 - w.ext.neg - w.ext.pos

    z_grid_ext <- z_grid
    wX_ext     <- wX * grid_share
    es_ext     <- es.X
    sig_ext    <- z_grid >= crit

    if (w.ext.neg > 0) {
      z_grid_ext <- c(min(z_grid) - xz, z_grid_ext)
      wX_ext     <- c(w.ext.neg, wX_ext)
      es_ext     <- c(es.ext.neg, es_ext)
      sig_ext    <- c(TRUE, sig_ext)
    }

    if (w.ext.pos > 0) {
      z_grid_ext <- c(z_grid_ext, max(z_grid) + xz)
      wX_ext     <- c(wX_ext, w.ext.pos)
      es_ext     <- c(es_ext, es.ext.pos)
      sig_ext    <- c(sig_ext, TRUE)
    }

    # ---- weighted MEAN over all reconstructed results ----------------------
    es_mean_all <- sum(es_ext * wX_ext, na.rm = TRUE)

    # ---- weighted mean over reconstructed SIGNIFICANT results only ---------
    if (any(sig_ext) && sum(wX_ext[sig_ext], na.rm = TRUE) > 0) {
      wX_sig           <- wX_ext
      wX_sig[!sig_ext] <- 0
      wX_sig           <- wX_sig / sum(wX_sig, na.rm = TRUE)
      es_mean_sig      <- sum(es_ext * wX_sig, na.rm = TRUE)
    } else {
      es_mean_sig <- NA_real_
    }

    # ---- weighted MEDIAN over all reconstructed results --------------------
    o      <- order(es_ext)
    es.srt <- es_ext[o]
    w.srt  <- wX_ext[o]
    keep   <- is.finite(es.srt) & is.finite(w.srt) & w.srt > 0
    es.srt <- es.srt[keep]
    w.srt  <- w.srt[keep]

    if (length(es.srt) >= 2 && sum(w.srt) > 0) {
      cum.w     <- cumsum(w.srt) / sum(w.srt)
      es_median <- approx(cum.w, es.srt, xout = 0.5, rule = 2,
                          ties = list("ordered", mean))$y
    } else {
      es_median <- NA_real_
    }

    # ---- local effect sizes by z-bin ----------------------------------------
    z_breaks <- seq(x.lim.min, Int.End, by = int.loc)
    if (tail(z_breaks, 1) < Int.End) z_breaks <- c(z_breaks, Int.End)

    local_es <- matrix(NA_real_, nrow = 5, ncol = length(z_breaks) - 1)
    rownames(local_es) <- c("z_min", "z_max", "weight", "mean_z", "es")
    colnames(local_es) <- paste0(head(z_breaks, -1), "-", tail(z_breaks, -1))

    for (j in seq_len(length(z_breaks) - 1)) {
      in_bin <- z_grid >= z_breaks[j] & z_grid < z_breaks[j + 1]
      if (j == length(z_breaks) - 1) {
        in_bin <- z_grid >= z_breaks[j] & z_grid <= z_breaks[j + 1]
      }
      if (!any(in_bin)) next

      bin_weight <- sum(wX[in_bin], na.rm = TRUE)
      if (!is.finite(bin_weight) || bin_weight <= 0) next

      w_bin <- wX[in_bin] / bin_weight

      local_es["z_min",  j] <- z_breaks[j]
      local_es["z_max",  j] <- z_breaks[j + 1]
      local_es["mean_z", j] <- sum(z_grid[in_bin] * w_bin, na.rm = TRUE)
      local_es["weight", j] <- bin_weight
      local_es["es",     j] <- sum(es.X[in_bin] * w_bin, na.rm = TRUE)
    }

    list(
      es_mean_all   = es_mean_all,
      es_median_all = es_median,
      es_mean_sig   = es_mean_sig,
      es_tau        = es_tau,
      local_es      = local_es[5, ],
      local_w       = local_es[3, ]
    )
  } ### EOF get_zcurve_es_estimates




  ################################
  ### Power: Discrete Components
  ################################
  Compute.Power.Z.Discrete <- function(cp.inp, Int.Beg, Int.End, z_grid,
                                       val.input, crit, Directional, alpha,
                                       x.lim.min, x.lim.max, int.loc) {


  Compute.Power.Z.Discrete <- function(cp.inp, Int.Beg, Int.End, z_grid,
                                       val.input, crit, Directional, alpha,
                                       x.lim.min, x.lim.max, int.loc,
                                       Folded = FALSE, xz = 0.01) {

    w.inp      <- cp.inp$w.inp
    ncp        <- cp.inp$ncp
    components <- length(ncp)

    n_total <- length(val.input)

    # extreme-value input PROPORTIONS (of total)
    ext.neg.inp <- sum(val.input < -Int.End) / n_total
    ext.pos.inp <- sum(val.input > Int.End) / n_total

    # component power (two-tailed, with sign error)
    pow.neg <- pnorm(-crit, ncp)
    pow.pos <- 1 - pnorm(crit, ncp)
    pow.dir <- c(pow.neg[ncp < 0], pow.pos[ncp >= 0])
    pow     <- pow.neg + pow.pos

    # power conditional on selection interval; sel.crit floors window bound at 0
    sel.crit <- max(0, Int.Beg)
    if (Directional)  pow.sel <- pnorm(ncp,      sel.crit) + pnorm(-sel.crit, ncp)
    if (!Directional) pow.sel <- pnorm(abs(ncp), sel.crit) + pnorm(-sel.crit, abs(ncp))

    # selection correction: input -> population weights (grid only)
    w.all <- w.inp / pow.sel
    w.all <- w.all / sum(w.all)

    # significant composition (grid only)
    w.sig <- w.all * pow
    w.sig <- w.sig / sum(w.sig)

    # extend vectors with extreme components (used only for EDR/ERR below)
    if (ext.neg.inp > 0 & Directional == FALSE) {
      ncp.ext = c(-(Int.End+1),ncp)
      pow.ext     <- c(1, pow)
      pow.dir.ext <- c(1, pow.dir)
      pow.sel.ext <- c(1, pow.sel)
      w.inp.ext   <- c(ext.neg.inp,w.inp * (1 - ext.neg.inp))
    } else {
      ncp.ext = ncp
      pow.ext = pow
      pow.dir.ext = pow.dir
      pow.sel.ext  = pow.sel
      w.inp.ext   <- w.inp
    } 

    if (ext.pos.inp > 0)  {
      ncp.ext = c(ncp,Int.End + 1)
      pow.ext     <- c(pow.ext, 1)
      pow.dir.ext <- c(pow.dir, 1)
      pow.sel.ext <- c(pow.sel, 1)
      w.inp.ext   <- c(w.inp * (1 - ext.pos.inp),ext.pos.inp)
    } else {
      ncp.ext = ncp.ext
      pow.ext = pow.ext
      pow.dir.ext = pow.dir.ext
      pow.sel.ext  = pow.sel.ext
      w.inp.ext   <- w.inp.ext
    } 

    stopifnot(length(pow.ext)     == length(w.inp.ext),
              length(pow.dir.ext) == length(w.inp.ext),
              length(pow.sel.ext) == length(w.inp.ext))

    w.all.ext <- w.inp.ext / pow.sel.ext
    w.all.ext <- w.all.ext / sum(w.all.ext)

    # EDR: average power over all studies (pre-selection)
    EDR <- sum(w.all.ext * pow.ext)

    # ERR: average directional power over significant studies
    w.sig.ext <- w.all.ext * pow.ext
    w.sig.ext <- w.sig.ext / sum(w.sig.ext)
    ERR <- sum(w.sig.ext * pow.dir.ext)

    ERR <- min(max(ERR, alpha / 2), 1)
    EDR <- min(max(EDR, alpha),     1)

    # mixture density at each grid point (grid-only ncp and w.all)
    # val.input is folded/absolute (two.sided <- TRUE), so the density of
    # |Z| for a component at mu needs the mirror term when Folded = TRUE.
    if (Folded) {
      lik.mat <- outer(z_grid, ncp, function(x, mu) dnorm(x - mu) + dnorm(x + mu))
    } else {
      lik.mat <- outer(z_grid, ncp, function(x, mu) dnorm(x - mu))
    }
    wd.X <- as.vector(lik.mat %*% w.all)

    # local power: posterior-weighted average of component power
    # (numerator and denominator share the same w.all-based lik.mat, so this
    #  ratio is scale-invariant -- no normalization of wd.X needed or wanted)
    loc.p <- as.vector((lik.mat %*% (w.all * pow)) / wd.X)

    # bin into intervals, density-weighted average within each bin
    int  <- seq(x.lim.min, x.lim.max, by = int.loc)
    bins <- cut(z_grid, breaks = int, include.lowest = TRUE)

    local_power <- as.vector(tapply(seq_along(z_grid), bins, function(idx) {
      denom <- sum(wd.X[idx])
      if (denom <= 0 || !is.finite(denom)) return(NA_real_)
      sum(loc.p[idx] * wd.X[idx]) / denom
    }))

    midpoints          <- (int[-length(int)] + int[-1]) / 2
    names(local_power) <- paste0("lp.", midpoints)

    # shape diagnostic: predicted vs observed nonsignificant z
    dx    <- z_grid[2] - z_grid[1]
    upper <- qnorm(1 - alpha)
    lower <- 0
    idx   <- z_grid > lower & z_grid < upper

    d_med  <- NA_real_
    d_mean <- NA_real_
    if (sum(idx) > 1) {
      X_ns <- z_grid[idx]
      pred <- wd.X[idx]
      pred <- pred / (sum(pred) * dx)

      pred_cdf  <- cumsum(pred) * dx
      pred_cdf  <- pred_cdf / pred_cdf[length(pred_cdf)]
      pred_med  <- approx(pred_cdf, X_ns, xout = 0.5, rule = 2)$y
      pred_mean <- sum(X_ns * pred) * dx

      obs <- val.input[val.input > lower & val.input < upper]
      if (length(obs) > 0) {
        d_med  <- median(obs) - pred_med
        d_mean <- mean(obs)   - pred_mean
      }
    }

    list(
      EDR            = EDR,
      ERR            = ERR,
      w.all          = w.all,
      w.all.ext      = w.all.ext,
      pow.ext        = pow.ext,
      pow.dir.ext    = pow.dir.ext,
      pow.sel.ext    = pow.sel.ext,
      ncp.ext        = ncp.ext,
      z_grid         = z_grid,
      wd.X           = wd.X,
      loc.p          = loc.p,
      local_power    = local_power,
      shape_d_median = d_med,
      shape_d_mean   = d_mean
    )
  } ### EOF Compute.Power.Z.Discrete




  ########################
  ### Main Power # ppp2
  ########################
  xz     <- 0.01
  z_grid <- seq(x.lim.min, Int.End, by = 0.01)

  out <- Compute.Power.Z.Discrete(
    cp.inp      = cp.inp,
    Int.Beg     = Int.Beg,
    Int.End     = Int.End,
    z_grid      = z_grid,
    val.input   = val.input,
    crit        = crit,
    Directional = Directional,
    alpha       = alpha,
    x.lim.min   = x.lim.min,
    x.lim.max   = x.lim.max,
    int.loc     = int.loc
   )

  #out$wd.X
  #sum(out$wd.X)*xz
  #w.all = out$w.all
  #out$local_power
  #se.bin = tapply(sei,cut(val.input,ncp),mean)
  #ncp.bin = qnorm(out$local_power,1.96)
  #ncp.bin
  #ncp.bin*se.bin
  #zval.bin = tapply(val.input,cut(val.input,ncp),mean)
  #zval.bin
  #yi.bin = tapply(yi,cut(val.input,ncp),mean)
  #yi.bin
 
  out$local_es      <- NULL
  out$es_mean_all   <- NULL
  out$es_median_all <- NULL
  out$es_mean_sig   <- NULL
  out$es_tau        <- NULL

  if (!is.null(yi)) {

    sei_est <- estimate_sei_function(
      z       = val.input,
      sei     = sei,
      z_grid  = out$z_grid,
      crit    = crit,
      z_max   = Int.End,
      z_min   = Int.Beg,
      Int.Beg = Int.Beg
    )

   es_est <- get_zcurve_es_estimates(
      loc.p       = out$loc.p,
      sei_est     = sei_est,
      z_grid      = out$z_grid,
      wd.X        = out$wd.X,
      ncp         = out$ncp.ext,
      w.all       = out$w.all.ext,
      pow         = out$pow.ext,
      pow.dir     = out$pow.dir.ext,
      pow.sel     = out$pow.sel.ext,
      val.input   = val.input,        # explicit
      yi          = yi,               # explicit
      Directional = Directional,      # explicit
      crit        = crit,
      Int.Beg     = Int.Beg,
      Int.End     = Int.End,
      int.loc     = int.loc,          # FIX: passed by name
      Folded      = Folded,
      xz          = xz
    )
 
    out$local_es      <- es_est$local_es
    out$es_mean_all   <- es_est$es_mean_all
    out$es_median_all <- es_est$es_median_all   # FIX: now propagated
    out$es_mean_sig   <- es_est$es_mean_sig
    out$es_tau        <- es_est$es_tau

    out$local_power 
    out$local_es

  } # EOF es 

  #round(cbind(w.all[1:6],zval.bin,ncp.bin,se.bin,yi.bin,out$local_es),2)

  out
 
  return(out)
 
} ### EOF Compute.Power.Z.General
 
### EOF Compute Power
 




############################################################
### SLOPE DIAGNOSTICS
############################################################

get_slope <- function(INT, from, to, bw) {

  d <- Get.Densities(INT, bw = bw, from = from, to = to,
       width = spline.width, Augment = TRUE)
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

#ncp = 2.5; zsds = 2; NCP_FIXED = FALSE; ZSDS_FIXED = FALSE; W_FIXED = TRUE; w.inp = 1

run_OF <- function(val.input,yi=NULL,sei=NULL,ODR,Int.Beg,Int.End,ncp,zsds, cola = "springgreen4",bw.est = .2,
                    Dens.pre = NULL, D.X.pre = NULL, width = 1, 
                    NCP_FIXED = TRUE, ZSDS_FIXED = TRUE, W_FIXED = FALSE,
                    x.lim.min,x.lim.max) {

  #width = 1
  #Folded = FALSE
  #Int.Beg = 1.96

  INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
  #summary(INT)
  ODR = mean(val.input > 1.96);ODR

  components <- length(ncp)

  #from = Int.Beg; to = Int.End; width = 1
  ## ---- Observed density ----
  densy <- Get.Densities(INT, bw = bw.est, from = Int.Beg, 
                         to = Int.End, width = width, Augment = Augment)

  D.X   <- densy[, 1]
  O.D.Y <- densy[, 2]
  n.bars    <- length(D.X)
  bar.width <- D.X[2] - D.X[1]

  #hist(INT,freq=FALSE,xlim=c(Int.Beg,Int.End),ylim=c(0,.5));par(new=TRUE);plot(D.X,O.D.Y,xlim=c(Int.Beg,Int.End),ylim=c(0,.5));abline(v=0)
  
  use.precomputed <- components > 1 &&
                     NCP_FIXED == TRUE &&
                     ZSDS_FIXED == TRUE


  if(use.precomputed) Dens.pre = build.dens(D.X, ncp, zsds, 
     Folded = Folded)
  #plot(D.X,Dens.pre[5,])

  ## ---- theta = [weights, ncp, zsds] — always ----
  n.wt <- components

  startval <- c(rep(1/n.wt, n.wt), ncp, zsds)

  lowlim   <- c(if (n.wt == 1) 1 else rep(0, n.wt),
                if (NCP_FIXED) ncp else rep(0,components),
                if (ZSDS_FIXED) zsds else rep(1,components)
			)
  highlim  <- c(rep(1, n.wt),
                if (NCP_FIXED) ncp else rep(6,components),
                if (ZSDS_FIXED) zsds else rep(6,components)
			)

  ## ---- Positions (always the same) ----
  idx.wt  <- 1:n.wt
  idx.ncp <- (n.wt + 1):(n.wt + components)
  idx.zsd <- (n.wt + components + 1):(n.wt + 2 * components)

  #theta = WT
  #startval = theta

  ## ---- Fitting function ----
  curve.fitting <- function(theta) {

    wt <- theta[idx.wt]
    wt <- wt / sum(wt)
 
    if (use.precomputed && NCP_FIXED && ZSDS_FIXED) {
      Dens.now <- Dens.pre
    } else {
      Dens.now <- build.dens(D.X, theta[idx.ncp], theta[idx.zsd], df, 
           CURVE.TYPE, Folded=Folded)
    }
    Dens.now <- Dens.now / (rowSums(Dens.now) * bar.width)

    E.D.Y <- colSums(Dens.now * wt)
    #plot(D.X,E.D.Y)
    rmse = sqrt(mean((E.D.Y - O.D.Y)^2))

    scale = 1
    ## ---- Diagnostic plot ----
    if (TESTING && runif(1) > 0.6) {
      plot(D.X, O.D.Y * scale, type = 'l',
           xlim = c(x.lim.min, x.lim.max), ylim = c(0, ymax),
           xlab = "", ylab = "", axes = FALSE)
      lines(D.X, E.D.Y * scale, col = "red1")
    }


    rmse

  }

  #TESTING = FALSE
  #TESTING = TRUE

  ## ---- Optimize ----
  auto <- nlminb(startval, curve.fitting,
                 lower = lowlim, upper = highlim,
                 control = list(eval.max = 1000, abs.tol = 1e-20))

  auto$par

  #abline(v=0)

  ## ---- Extract — positions are guaranteed ----
  WT <- auto$par[idx.wt]
  WT <- WT / sum(WT)

  of.est =  list(
    w.inp   = WT,
    ncp     = auto$par[idx.ncp],
    zsds    = auto$par[idx.zsd],
    fit     = auto$objective
  )

  #print(of.est)

  cp.inp = of.est
 
  #print(cp.inp)
  #print("OF power")


  ## ---- Compute power for each bootstrap result ----
  cp.res = Compute.Power.Z.General(
        cp.inp = cp.inp, 
        val.input = val.input,
        Int.Beg = Int.Beg,
        Int.End = Int.End, 
        yi = yi, 
        sei = sei,
        x.lim.min = x.lim.min,
        x.lim.max = x.lim.max)

  #print(cp.res$w.all)
  #sd(cp.res$w.all)
  #print(cp.res$es_mean_all)
  #print(cp.res$es_median_all)
  #print(cp.res$es_tau)
  #print(sim_run["true_tau"])

  ret = list(
    w.inp           = of.est$w.inp,
    ncp             = of.est$ncp,
    zsds            = of.est$zsds,
    model.fit       = of.est$fit,
    ODR             = ODR,   
    EDR             = cp.res$EDR,
    ERR             = cp.res$ERR,
    w.all           = cp.res$w.all, 
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


} ### EOF run_OF



##################################################
### ROOT EM FUNCTION
##################################################


#eee

run_EM <- function(
  val.input,
  yi, 
  sei, 
  Int.Beg,
  Int.End,
  ncp,
  zsds,
  NCP_FIXED  = TRUE,
  ZSDS_FIXED = TRUE,
  W_FIXED    = FALSE,
  x.lim.min,
  x.lim.max,
  max_iter   = 500,
  tol        = 1e-8) {


n_starts = 1

components = length(ncp)
w.inp = rep(1/components,components)

INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
ODR = mean(val.input > 1.96);ODR

EM_EXT_SIGNED <- function(
  INT,
  ncp,
  zsds,
  Int.Beg = -Inf,        # signed lower window bound (e.g. -8, or -Inf for none)
  Int.End =  Inf,        # signed upper window bound
  NCP_FIXED  = TRUE,
  ZSDS_FIXED = TRUE,
  W_FIXED    = FALSE,
  max_iter = 1000,
  tol = 1e-8
) {

  x <- INT
  x <- x[!is.na(x)]
  k.int <- length(x)
  components <- length(ncp)

  loglik_old <- -Inf
  loglik_trace <- numeric(max_iter)

  # ---- Newton update for a single component's mean (signed, truncated) ----
  update_ncp_newton_signed <- function(
    mu_init, sigma, tau_k, x, Int.Beg, Int.End,
    max_iter = 20, tol = 1e-8, lower = -10, upper = 10
  ) {
    mu  <- mu_init
    n_k <- sum(tau_k)
    for (iter in 1:max_iter) {
      # no folding: sufficient statistic is just tau-weighted sum of x
      S1 <- sum(tau_k * x)

      # single truncation interval [Int.Beg, Int.End] on the signed axis
      a <- (Int.Beg - mu) / sigma
      b <- (Int.End - mu) / sigma
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
  update_zsds_newton_signed <- function(
    sigma_init, mu, tau_k, x, Int.Beg, Int.End,
    max_iter = 20, tol = 1e-6, lower = 0.5, upper = 10
  ) {
    sigma <- sigma_init
    n_k   <- sum(tau_k)
    for (iter in 1:max_iter) {
      # no folding: second moment about mu directly
      S2 <- sum(tau_k * (x - mu)^2)

      a <- (Int.Beg - mu) / sigma
      b <- (Int.End - mu) / sigma
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
      sigma <- zsds[k]

      # single signed density, no + dnorm(x, -mu)
      log_dens <- dnorm(x, mu, sigma, log = TRUE)

      # single truncation interval on the signed axis
      a <- (Int.Beg - mu) / sigma
      b <- (Int.End - mu) / sigma
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
          mu_init = ncp[k], sigma = zsds[k], tau_k = tau[, k],
          x = x, Int.Beg = Int.Beg, Int.End = Int.End
        )
      }
    }
    if (!ZSDS_FIXED) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        zsds[k] <- max(1, update_zsds_newton_signed(
          sigma_init = zsds[k], mu = ncp[k], tau_k = tau[, k],
          x = x, Int.Beg = Int.Beg, Int.End = Int.End
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

  list(w.inp = w.inp, ncp = ncp, zsds = zsds, fit = loglik)
}




EM_EXT_FOLDED <- function(
  INT,
  ncp,
  zsds,
  Int.Beg = 1.96,
  Int.End = max(INT),
  NCP_FIXED  = TRUE,
  ZSDS_FIXED = TRUE,
  W_FIXED    = FALSE,
  max_iter = 1000,
  tol = 1e-8
) {


x <- INT
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
  Int.Beg,
  Int.End,
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
    a1 <- ( Int.Beg - mu) / sigma
    b1 <- ( Int.End - mu) / sigma
    a2 <- (-Int.End - mu) / sigma
    b2 <- (-Int.Beg - mu) / sigma

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


update_zsds_newton_folded <- function(
  sigma_init,
  mu,
  tau_k,
  x,
  Int.Beg,
  Int.End,
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
    a1 <- ( Int.Beg - mu) / sigma
    b1 <- ( Int.End - mu) / sigma
    a2 <- (-Int.End - mu) / sigma
    b2 <- (-Int.Beg - mu) / sigma

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
      sigma <- zsds[k]
      l1 <- dnorm(x,  mu, sigma, log = TRUE)
      l2 <- dnorm(x, -mu, sigma, log = TRUE)
      m  <- pmax(l1, l2)
      log_fold <- m + log(exp(l1 - m) + exp(l2 - m))
      a1 <- ( Int.Beg - mu) / sigma
      b1 <- ( Int.End - mu) / sigma
      a2 <- (-Int.End - mu) / sigma
      b2 <- (-Int.Beg - mu) / sigma
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
          sigma   = zsds[k],
          tau_k   = tau[, k],
          x       = x,
          Int.Beg = Int.Beg,
          Int.End = Int.End
        )
      }
    }
    if (!ZSDS_FIXED) {
      for (k in 1:components) {
        if (n_k[k] < 1e-10) next
        zsds[k] <- max(1, update_zsds_newton_folded(
          sigma_init = zsds[k],
          mu         = ncp[k],
          tau_k      = tau[, k],
          x          = x,
          Int.Beg    = Int.Beg,
          Int.End    = Int.End
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
    zsds         = zsds,
    fit          = loglik
  )

  
}



### Main
### eee2

  components <- length(ncp)
  ncp_init   <- ncp
  zsds_init  <- zsds

  best_fit <- NULL
  best_ll  <- -Inf

  print("Running EM")
  print(paste0("Folded: ",Folded))

  for (s in 1:n_starts) {

      raw <- rgamma(components, shape = 1)   # random Dirichlet
      w_s <- raw / sum(raw)

      if(Folded) {

        fit <- EM_EXT_FOLDED(
          INT        = INT,
          ncp        = ncp,
          zsds       = zsds,
          Int.Beg    = Int.Beg,
          Int.End    = Int.End,
          NCP_FIXED  = NCP_FIXED,
          ZSDS_FIXED = ZSDS_FIXED,
          W_FIXED    = W_FIXED,
          max_iter   = max_iter,
          tol        = tol
        ) 

      } else {

        fit <- EM_EXT_SIGNED(
          INT        = INT,
          ncp        = ncp,
          zsds       = zsds,
          Int.Beg    = Int.Beg,
          Int.End    = Int.End,
          NCP_FIXED  = NCP_FIXED,
          ZSDS_FIXED = ZSDS_FIXED,
          W_FIXED    = W_FIXED,
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
  cp.res = Compute.Power.Z.General(
        cp.inp = cp.inp, 
        val.input = val.input,
        Int.Beg = Int.Beg,
        Int.End = Int.End, 
        yi = yi, 
        sei = sei,
        x.lim.min = x.lim.min,
        x.lim.max = x.lim.max)

  ret = list(
    w.inp           = em.est$w.inp,
    ncp             = em.est$ncp,
    zsds            = em.est$zsds,
    model.fit       = em.est$fit,
    ODR             = ODR,   
    EDR             = cp.res$EDR,
    ERR             = cp.res$ERR,
    w.all           = cp.res$w.all, 
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
                   val.input,
                   yi,
                   sei,
                   crit,   
                   ncp,
                   zsds,
                   Int.Beg,
                   Int.End,
                   Directional = FALSE,
                   Folded = TRUE,
                   Dens.pre = NULL,
                   D.X.pre = NULL,
                   width = 1,
                   boot.iter  = 500,
                   NCP_FIXED  = TRUE,
                   ZSDS_FIXED = TRUE,
                   W_FIXED    = FALSE,
                   cluster.id = NULL,
                   Est.Method = "OF",
                   x.lim.min,
                   x.lim.max
                  ) {


    
      #NCP_FIXED = FALSE; ZSDS_FIXED = FALSE; W_FIXED = TRUE; ncp = 2.5; zsds = 2
      #ncp = 0:6;zsds = rep(1,length(ncp))
      #cluster.id = NULL
      #boot.iter = 100
    
      #INT = val.input[val.input >= Int.Beg & val.input <= Int.End]

      width = 1
      z = NULL

      components = length(ncp)

      if (!is.null(cluster.id)) {

          print("USING CLUSTERS")

          ## ---- Cluster Bootstrap ----

          stopifnot(length(val.input) == length(cluster.id))	
          z.cluster.rows  <- split(seq_along(val.input), cluster.id)
          z.cluster.names <- names(z.cluster.rows)
          n.z.clusters    <- length(z.cluster.names)

      } else {

          z.cluster.rows = NULL; z.cluster.names = NULL; n.z.clusters = NULL

      }

 
      #Est.Method = "OF";width = 1
      #Est.Method = "EM"

      var.list = c("Est.Method","val.input","crit",
           "cluster.id","z.cluster.rows","z.cluster.names","n.z.clusters",
           "ncp", "zsds","Int.Beg", "Int.End", "Augment",
           "CURVE.TYPE","Directional","Folded","x.lim.min", "x.lim.max", 
           "NCP_FIXED", "ZSDS_FIXED","W_FIXED",
           "Compute.Power.Z.General","alpha","int.loc",
           "x.lim.min","x.lim.max")


      if (!is.null(yi)) var.list = c(var.list,"yi","sei")

      boot_res = NULL

      if (Est.Method == "OF") {

         var.list = c(var.list,"run_OF","build.dens", "Get.Densities",
                  "bw.est","width","Dens.pre", "D.X.pre")

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
              sampled  <- sample(z.cluster.names, n.z.clusters, replace = TRUE)
              rows     <- unlist(z.cluster.rows[sampled], use.names = FALSE)
              boot_all <- val.input[rows]
          } else {
             rows     <- sample(seq_along(val.input), length(val.input), replace = TRUE)
             boot_all <- val.input[rows]
          }

          boot_sample <- boot_all[boot_all >= Int.Beg & boot_all <= Int.End]

          if (!is.null(yi)) { yi_boot <- yi[rows]; sei_boot <- sei[rows] 
           } else { yi_boot = NULL; sei_boot = NULL }

          fit <- #tryCatch(
                 run_OF(
                 val.input  = boot_all,
                 yi         = yi, # yi_boot,
                 sei        = sei, #, _boot,
                 Int.Beg    = Int.Beg,
                 Int.End    = Int.End,
                 ncp        = ncp,
                 zsds       = zsds,
                 Dens.pre   = Dens.pre,
                 D.X.pre    = D.X.pre,
                 bw.est     = bw.est, 
                 width      = width,
                 NCP_FIXED  = NCP_FIXED,
                 ZSDS_FIXED = ZSDS_FIXED,
                 W_FIXED    = W_FIXED,
                 x.lim.min  = x.lim.min,
                 x.lim.max  = x.lim.max
                 )
	         #,error = function(e) NULL )

          fit.res = fit
  
      }) # EOF boot_res

      if(is.null(boot_res)) stop("boostrap failure.")

    } # EOF OF

     # boot_res

      if (Est.Method == "EM") {

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
              boot_all <- val.input[rows]
          } else {
             rows     <- sample(seq_along(val.input), length(val.input), replace = TRUE)
             boot_all <- val.input[rows]
          }

          boot_sample <- boot_all[boot_all >= Int.Beg & boot_all <= Int.End]

          if (!is.null(yi)) { yi_boot <- yi[rows]; sei_boot <- sei[rows] 
           } else { yi_boot = NULL; sei_boot = NULL }

          fit.EM <- #tryCatch(
                 run_EM(
                 val.input  = boot_all,
                 yi         = yi, # yi_boot,
                 sei        = sei, #, _boot,
                 Int.Beg    = Int.Beg,
                 Int.End    = Int.End,
                 ncp        = ncp,
                 zsds       = zsds,
                 NCP_FIXED  = NCP_FIXED,
                 ZSDS_FIXED = ZSDS_FIXED,
                 W_FIXED    = W_FIXED,
                 x.lim.min  = x.lim.min,
                 x.lim.max  = x.lim.max
                 )
	         #,error = function(e) NULL )

          fit.EM 

         }) # close boot_res

    } # EOF EM 

    #boot_res 

    #Est.Method = "ML_gamma"; boot.iter = 500
    if (Est.Method == "ML_gamma") {

         n_starts = 2
         if (Est.Method == "ML_gamma") var.list = c(var.list,"run_ML_gamma","n_starts")

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
               boot_sample <- boot_all[boot_all >= Int.Beg & boot_all <= Int.End]
           } else {
               boot_all <- sample(val.input, length(val.input), replace = TRUE)
               odr_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= Int.Beg & boot_all <= Int.End]
           }


           gamma_fit <- tryCatch(
                 run_ML_gamma(boot_all, Int.Beg = Int.Beg, Int.End = Int.End, crit = crit,
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

#yi=NULL;sei=NULL;
#width = 1
#ncp = 2.5; zsds = 3;NCP_FIXED=FALSE;ZSDS_FIXED=FALSE;W_FIXED=FALSE
#ncp = 0:6;zsds=rep(1,length(ncp))
#NCP_FIXED = FALSE; ZSDS_FIXED = FALSE; W_FIXED = TRUE; ncp = 2.5; zsds = 2; Run.Gamma=FALSE
#yi = NULL;sei = NULL
#cluster.id = NULL
#boot.iter = 0
#hist(val.input)
#val.input = c(rnorm(1500,0),rnorm(500,2),rnorm(100,4))


run_zcurve = function(Est.Method,val.input,crit,Int.Beg,Int.End,
           ncp,zsds,NCP_FIXED=TRUE,ZSDS_FIXED = TRUE,W_FIXED = FALSE,
           cluster.id,boot.iter=0,edr_adj_ci=0,Run.Gamma=FALSE,
           yi=NULL,sei=NULL,Folded=Folded,
           x.lim.min,x.lim.max) {

  zcurve.time = system.time({

  ### start
  #Int.Beg = 0

  INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
  #summary(val.input)
  #summary(INT) 

  ODR = mean(val.input > crit)

  width = spline.width

  D.X.pre  <- NULL
  Dens.pre <- NULL

 ## ---- Precompute Dens if both fixed ----
 if (Est.Method == "OF" && NCP_FIXED && ZSDS_FIXED) {
      densy     <- Get.Densities(INT, bw = bw.est, from = Int.Beg,
                             to = Int.End, width = width, Augment = Augment)
      D.X.pre   <- densy[, 1]
      Dens.pre  <- build.dens(D.X.pre, ncp, zsds, df, 
           CURVE.TYPE, Folded=Folded)
      bar.width <- D.X.pre[2] - D.X.pre[1]
      Dens.pre  <- Dens.pre / (rowSums(Dens.pre) * bar.width)
 }

 if(Est.Method == "OF") { 

         print("Using OF")

         res.of <- run_OF(
                 val.input = val.input, 
                 yi = yi,
                 sei = sei,
                 Int.Beg = Int.Beg,
                 Int.End = Int.End,   
                 ncp = ncp,
                 zsds = zsds,
                 Dens.pre = Dens.pre,
                 D.X.pre  = D.X.pre,
                 bw.est = bw.est,
                 width = width,
                 NCP_FIXED  = NCP_FIXED,
                 ZSDS_FIXED = ZSDS_FIXED,
                 W_FIXED = W_FIXED,
                 x.lim.min  = x.lim.min,
                 x.lim.max  = x.lim.max
        )

        res.of
        res.run = res.of 
  }  


   if (Est.Method == "EM") { 

      print("Using EM")

      #summary(val.input)
         res.em <- run_EM(
                 val.input = val.input, 
                 yi = yi,
                 sei = sei,
                 Int.Beg = Int.Beg,
                 Int.End = Int.End,   
                 ncp = ncp,
                 zsds = zsds,
                 NCP_FIXED  = NCP_FIXED,
                 ZSDS_FIXED = ZSDS_FIXED,
                 W_FIXED = W_FIXED,
                 x.lim.min  = x.lim.min,
                 x.lim.max  = x.lim.max
        )


     #cp.inp = res.em
     res.em
     res.run = res.em
     
   }

  if (Est.Method == "SQP") { 

      print("Using SQP")

      res.sqp <- run_SQP(
        val.input = val.input, 
        ncp = ncp,
        zsds = zsds,
        Int.Beg = Int.Beg,
        Int.End = Int.End,
        Folded = FALSE
      )

     res.sqp
     res.run = res.sqp

  }

  #res.run

  # comp = round(cbind(ncp,res.of$w.inp,res.em$w.inp,res.sqp),2) 
  # colnames(comp) = c("ncp","OF","EM","SQP")
  # comp


  #res.run


  #plot(es_pe$w.est,es_pe$w.obs,xlim=c(0,.25),ylim=c(0,.25))

    res.pe <- list(
      ODR              = ODR,
      EDR              = res.run$EDR,
      ERR              = res.run$ERR,
      ncp              = res.run$ncp,
      zsds             = res.run$zsds,
      w.inp            = res.run$w.inp,
      w.all            = res.run$w.all,
      local_power      = res.run$local_power,
      model.fit        = res.run$model.fit,
      gamma_shape      = NA,
      gamma_rate       = NA,
      gamma_edr        = c(NA,NA),
      gamma_er         = c(NA,NA),
      ODR_EDR_D        = ODR - res.run$EDR,
      shape_d_median   = res.run$shape_d_median,
      shape_d_mean     = res.run$shape_d_mean,
      local_es_pe      = res.run$local_es,
      es_mean_all_pe   = res.run$es_mean_all,
      es_median_all_pe = res.run$es_median_all,
      es_mean_sig_pe   = res.run$es_mean_sig,
      es_tau_pe        = res.run$es_tau
     )

     print(res.pe$es_tau_pe)

     #res.pe 
     #Run.Gamma = FALSE

     if(Run.Gamma)  {

       gamma_fit <- tryCatch(
             run_ML_gamma(val.input, Int.Beg = Int.Beg, Int.End = Int.End, crit = crit,
                n_starts = n_starts)
             ,error = function(e) NULL )

       if(!is.null(gamma_fit)) {

          res.pe$gamma_shape   = gamma_fit$gamma_shape
          res.pe$gamma_rate    = gamma_fit$gamma_rate
          res.pe$gamma_edr     = gamma_fit$gamma_edr
          res.pe$gamma_err     = gamma_fit$gamma_err

       }

     }

  #################################

  ### finsh # res.pe

  #################################


  ### start boostrap

  #Est.Method = "OF"; boot.iter = 50
  #cluster.id = NULL
  if (boot.iter > 0) {  

       print("Begin of Bootstrap")
   
       boot_res = NULL

       boot_res <- run_bootstrap(
          val.input   = val.input, 
          yi          = yi,
          sei         = sei,
          crit        = crit,
          ncp         = res.pe$ncp,
          zsds        = res.pe$zsds,
          Int.Beg     = Int.Beg,
          Int.End     = Int.End,
          boot.iter   = boot.iter,
          Directional = Directional,
          Folded      = Folded,
          Dens.pre    = Dens.pre,
          D.X.pre     = D.X.pre,
          NCP_FIXED   = NCP_FIXED,
          ZSDS_FIXED  = ZSDS_FIXED,
          W_FIXED     = W_FIXED,
          cluster.id  = cluster.id,
          Est.Method  = Est.Method,
          x.lim.min   = x.lim.min,
          x.lim.max   = x.lim.max)	

        ###

        if(Est.Method == "ML_gamma") { 

          boot_mat <- do.call(rbind, lapply(boot_list, function(x) {
              unlist(x) }))

          ci_lo <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
          ci_hi <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

          ci <- cbind(ci_lo, ci_hi)

          gamma_edr = t(c(res.pe$gamma_edr,ci["gamma_edr",]))
          gamma_err = t(c(res.pe$gamma_err,ci["gamma_err",]))

          ODR_EDR_D <- ODR_pe - EDR_pe

          zcurve.res = list(
                ODR             = ODR,
                EDR             = res.pe$EDR,
                ERR             = res.pe$ERR,
                ncp             = rbind(res.pe$ncp,matrix(NA,2,length(res.pe$ncp))),
                zsds            = rbind(res.pe$zsds,matrix(NA,2,length(res.pe$zsds))),
                w.inp           = rbind(res.pe$w.inp,matrix(NA,2,length(res.pe$w.inp))),
                w.all           = rbind(res.pe$w.all,matrix(NA,2,length(res.pe$w.all))),
                local_power     = rbind(res.pe$local_power,matrix(NA,2,length(res.pe$local_power))), 
                local_es        = rbind(res.pe$local_es_pe,matrix(NA,2,length(res.pe$local_es_pe))), 
                fit             = res.pe$fit,
#                gamma_shape     = res.pe$gamma_shape,
#                gamma_rate      = res.pe$gamma_rate,      
#                gamma_edr       = gamma_edr,
#                gamma_err       = gamma_err,
                ODR_EDR_D       = c(res.pe$ODR_EDR_D,NA,NA),
                es_mean_all     = c(res.pe$es_mean_all_pe,NA,NA),
                es_median_all   = c(res.pe$es_median_all_pe,NA,NA),
                es_mean_sig     = c(res.pe$es_mean_sig_pe,NA,NA),
                es_tau          = c(res.pe$es_tau_pe,NA,NA)
            )


         } else {

          print("End of Bootstrap")

          ODR_boot   <- sapply(boot_res, function(x) x$ODR)
          ODR_boot[is.nan(ODR_boot)] = NA
          ODR_ci = c(
             quantile(ODR_boot, probs = .025, na.rm = TRUE),
             quantile(ODR_boot, probs = .975, na.rm = TRUE)
          )


          EDR_boot   <- sapply(boot_res, function(x) x$EDR)

          ERR_boot   <- sapply(boot_res, function(x) x$ERR)

          EDR_boot_lb <- pmax(alpha, EDR_boot - CI.EDR.MIN.ADJ)
          EDR_boot_ub <- pmin(1,     EDR_boot + CI.EDR.MIN.ADJ)

          EDR_ci = c(
             quantile(EDR_boot_lb, probs = .025, na.rm = TRUE),
             quantile(EDR_boot_ub, probs = .975, na.rm = TRUE)
          )

 
          ERR_ci = c(
             quantile(ERR_boot, probs = .025, na.rm = TRUE),
             quantile(ERR_boot, probs = .975, na.rm = TRUE)
          )


          ODR_EDR_D_ci  <- c(
               quantile(ODR_boot - EDR_boot_lb, .025, na.rm = TRUE),
               quantile(ODR_boot - EDR_boot_ub, .975, na.rm = TRUE)
           )
          ODR_EDR_D = c(res.pe$ODR_EDR_D,ODR_EDR_D_ci)

          if(!is.null(yi)) {

            es_mean_all_boot = sapply(boot_res, function(x) x$es_mean_all)

            if(max(es_mean_all_boot) + .05 < res.pe$es_mean_all_pe) {
               print("check CI")    
               print("check CI")    
               print("check CI")    
               print(summary(es_mean_all_boot))
               print(res.pe$es_mean_all_pe)
               print("check CI")    
               print("check CI")    
               print("check CI")    
               #stop("check ES CI")
            } 

            es_mean_all_ci = quantile(es_mean_all_boot, probs = c(.025,.975), na.rm=TRUE)
            es_mean_all = c(res.pe$es_mean_all_pe,es_mean_all_ci)


            es_median_all_boot = sapply(boot_res, function(x) x$es_median_all)
            es_median_all_ci = quantile(es_median_all_boot, probs = c(.025,.975), na.rm=TRUE)
            es_median_all = c(res.pe$es_median_all_pe,es_median_all_ci)

            es_tau_boot = sapply(boot_res, function(x) x$es_tau)
            es_tau_ci = quantile(es_tau_boot, probs = c(.025,.975), na.rm=TRUE)
            es_tau = c(res.pe$es_tau_pe,es_tau_ci)


            es_mean_sig_boot = sapply(boot_res, function(x) x$es_mean_sig)

            es_mean_sig_ci = quantile(es_mean_sig_boot, probs = c(.025,.975), na.rm=TRUE)
            es_mean_sig = c(res.pe$es_mean_sig_pe,es_mean_sig_ci)

            es.loc = sapply(boot_res, function(x) x$local_es)

            if(length(res.pe$local_es_pe) == 1) {
               local_es_ci = quantile(es.loc, c(.025, .975),na.rm=TRUE)   # CI straight from the percentiles
               local_es = c(res.pe$local_es_pe,local.es_ci)
            } else {
               local_es_ci = apply(es.loc,1,function(x) quantile(x,c(.025, .975),na.rm=TRUE) )   # CI straight from the percentiles
               local_es  = rbind(res.pe$local_es, local_es_ci) 
            }


          } # EOF if yi 

          #print(es_mean_all)
          #print(es_mean_sig)
          #print(local_es)

          if(length(res.pe$ncp) == 1) {
             ncp_ci = quantile(sapply(boot_res, function(x) x$ncp), c(.025, .975))   # CI straight from the percentiles
             ncp = c(res.pe$ncp,ncp_ci)
          } else {
             ncp_ci = apply(sapply(boot_res, function(x) x$ncp),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             ncp  = rbind(res.pe$ncp, ncp_ci) 
          }

          if(length(res.pe$zsds) == 1) {
             zsds_ci = quantile(sapply(boot_res, function(x) x$zsds), c(.025, .975))   # CI straight from the percentiles
             zsds = c(res.pe$zsds,zsds_ci)
          } else {
             zsds_ci = apply(sapply(boot_res, function(x) x$zsds),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             zsds  = rbind(res.pe$zsds, zsds_ci) 
          }

          if(length(res.pe$w.all) == 1) {
             w.all_ci = quantile(sapply(boot_res, function(x) x$w.all), c(.025, .975))   # CI straight from the percentiles
             w.all = c(res.pe$w.all,w.all_ci)
          } else {
             w.all_ci = apply(sapply(boot_res, function(x) x$w.all),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             w.all  = rbind(res.pe$w.all, w.all_ci) 
          }


          if(length(res.pe$w.inp) == 1) {
             w.inp_ci = quantile(sapply(boot_res, function(x) x$w.inp), c(.025, .975))   # CI straight from the percentiles
             w.inp = c(res.pe$w.inp,w.inp_ci)
          } else {
             w.inp_ci = apply(sapply(boot_res, function(x) x$w.inp),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             w.inp  = rbind(res.pe$w.inp, w.inp_ci) 
          }

          if(length(res.pe$local_power) == 1) {
             local_power_ci = quantile(sapply(boot_res, function(x) x$local_power), c(.025, .975))   # CI straight from the percentiles
             local_power = c(res.pe$local_power,local_power_ci)
          } else {
             local_power_ci = apply(sapply(boot_res, function(x) x$local_power),1,function(x) quantile(x,c(.025, .975)) )   # CI straight from the percentiles
             local_power  = rbind(res.pe$local_power, local_power_ci) 
          }



          if(!is.na(res.pe$shape_d_median)) {
            shape_d_median = c(res.pe$shape_d_median,quantile(sapply(boot_res, function(x) x$shape_d_median), c(.025, .975)) )   # CI straight from the percentiles
            shape_d_mean   = c(res.pe$shape_d_mean,quantile(sapply(boot_res, function(x) x$shape_d_mean), c(.025, .975)) )   # CI straight from the percentiles
          } else {
            shape_d_median = c(NA,NA,NA)
            shape_d_mean   = c(NA,NA,NA)
          }

          ODR = c(res.pe$ODR,ODR_ci)
          EDR = c(res.pe$EDR,EDR_ci)
          ERR = c(res.pe$ERR,ERR_ci)

          EDR[2] = max(alpha,EDR[2])
          EDR[3] = min(1,EDR[3])

          ERR[2] = ERR[2] - CI.ERR.MIN.ADJ
          ERR[3] = ERR[3] + CI.ERR.MIN.ADJ

          ERR[2] = max(alpha/2,ERR[2])
          ERR[3] = min(1,ERR[3])

       } # End of OF or EM

  } else { # End of if Bootstrap


        print(res.pe$es_median_all_pe)

        ODR = t(c(res.pe$ODR,rep(NA,2)))
        EDR = t(c(res.pe$EDR,rep(NA,2)))
        ERR = t(c(res.pe$ERR,rep(NA,2)))
        ncp = t(cbind(res.pe$ncp, matrix(NA,length(res.pe$w.all),2)))
        zsds = t(cbind(res.pe$zsds, matrix(NA,length(res.pe$w.inp),2)))
        w.all = t(cbind(res.pe$w.all, matrix(NA,length(res.pe$w.all),2)))
        w.inp = t(cbind(res.pe$w.all, matrix(NA,length(res.pe$w.inp),2)))
        local_power = t(cbind(res.pe$local_power,matrix(NA,length(res.pe$local_power),2)))
        local_es = t(cbind(res.pe$local_es_pe,matrix(NA,length(res.pe$local_es_pe),2)))
        ODR_EDR_D = c(res.pe$ODR_EDR_D,NA,NA)
        shape_d_median  = c(res.pe$shape_d_median,NA,NA)
        shape_d_mean    = c(res.pe$shape_d_mean,NA,NA)
        es_mean_all     = c(res.pe$es_mean_all_pe,NA,NA)
        es_median_all   = c(res.pe$es_median_all_pe,NA,NA)
        es_mean_sig     = c(res.pe$es_mean_sig_pe,NA,NA)
        es_tau          = c(res.pe$es_tau_pe,NA,NA)

  } # EOF no bootstrap

  zcurve.res = list(
            ODR                = ODR,
            EDR                = EDR,
            ERR                = ERR,
            ncp                = ncp,
            zsds               = zsds,
            w.inp              = w.inp,
            w.all              = w.all,
            local_power        = local_power,
            local_es           = local_es,
            gamma_shape        = res.pe$gamma_shape,
            gamma_rate         = res.pe$gamma_rate,
            gamma_edr          = res.pe$gamma_edr,
            gamma_err          = res.pe$gamma_err,
            ODR_EDR_D          = ODR_EDR_D,
            shape_d_median     = shape_d_median,
            shape_d_mean       = shape_d_mean,
            es_mean_all        = es_mean_all,
            es_median_all      = es_median_all,  
            es_mean_sig        = es_mean_sig,
            es_tau             = es_tau
          )


});print(zcurve.time) ### End of Timer


#print(zcurve.res$es_tau)


#zcurve.res

#print(zcurve.res)

return(zcurve.res)

}


############################################

Write_Local_Power = function(local_power) {

  if (int.loc == 0) return(invisible(NULL))

  # Midpoints of each bin — these are the correct label positions
  midpoints = seq(x.lim.min, x.lim.max - int.loc, by = int.loc) + int.loc / 2

  # Format labels
  lab = paste0(round(local_power * 100), "%")

  mtext(lab, side = 1, line = 1.5, at = midpoints, cex = 1.0, las = 1)

} ### EOF Write.Local.Power


# l_p = local_power[1,];l_es = local_es[1,]

Write_Local_Power_ES = function(l_p, l_es) {

  if (int.loc == 0) return(invisible(NULL))

  # Midpoints of each bin — these are the correct label positions
  midpoints = seq(x.lim.min, x.lim.max - int.loc, by = int.loc) + int.loc / 2

  #l_es = zres.tmz$local_es[1,]
   l_es

  # Format labels
  lab1 <- paste0(round(l_p * 100), "%")
  lab2 <- sprintf("%.2f", l_es)
#  lab2 <- paste0(paste(x[-length(x)], collapse = ", "), ", and ", x[length(x)])

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

	z.hist = val.input[val.input > x.lim.min & val.input < x.lim.max -.04] + .04
   

	int.start = round(Int.Beg,1)
	if (round(Int.Beg,2) == 1.96) {
         Int.start = 2
	}

	n.breaks = seq(x.lim.min,x.lim.max,hist.bar.width);n.breaks

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

	hist(z.hist[z.hist > round(Int.Beg,1) & z.hist < Int.End],
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

	min.z = min(val.input)
	max.z = max(val.input)
	n.z = length(val.input)
	n.z.sig = length(val.input[val.input > crit])
	n.not.shown = length(val.input[val.input > x.lim.max])

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

	#z.draw = val.input
	#cola = "blue"

	x.adj.max = x.lim.max + 3 * bw.draw

	### densitiy line

	summary(z.draw)

	d.all = Get.Densities(z.draw[z.draw >= x.lim.min & z.draw < x.lim.max],
		bw=bw.draw,from=x.lim.min,to=x.lim.max,Augment=Augment)
	summary(d.all)

	bar.width = d.all[2,1] - d.all[1,1];bar.width
	d.all = d.all[d.all[,1] > x.lim.min & d.all[,1] <= x.lim.max,]
	dim(d.all)
	sum(d.all[,2])*bar.width
	d.all.X = d.all[,1]
	d.all.Y = d.all[,2]
	summary(d.all.X)


	summary(z.draw)

	d.sig = Get.Densities(z.draw[z.draw >= Int.Beg & z.draw < x.lim.max],
		bw=bw.draw,from=Int.Beg,to=x.lim.max,Augment=Augment)
	
	bar.width = d.sig[2,1] - d.sig[1,1];bar.width

	dim(d.sig)
	d.sig = d.sig[d.sig[,1] > Int.Beg & d.sig[,1] <= x.lim.max,]

	#dim(d.sig)
	#sum(d.sig[,2])*bar.width

	d.sig.X = d.sig[,1]
	d.sig.Y = d.sig[,2]
	#summary(d.sig.X)

    #plot(d.sig.X,d.sig.Y)



	d.fit = length(z.draw[z.draw > Int.Beg & z.draw < x.lim.max]) /
		length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max])

	d.fit

	d.sig.Y = d.sig.Y * d.fit

	sum(d.sig[,2])*bar.width
	

	### draw the density of the observed values in the selected region
	par(new=TRUE)
	plot(d.sig.X[d.sig.X > Int.Beg +.05 & d.sig.X < x.lim.max - .05],
		d.sig.Y[d.sig.X > Int.Beg + .05 & d.sig.X < x.lim.max - .05],
		type="l",col=adjustcolor(cola, alpha.f = 0.8),lty=2,lwd=Lwidth,
		xlim =c(x.lim.min,x.lim.max),ylim=c(0,ymax),xlab="",ylab="",
		axes=FALSE)
	### draw vertical line for Beginning

	Int.Beg
	#height = max(d.sig.Y[which(round(d.sig.X,1) == round(Int.Beg+.02,1))]);height
	#segments(Int.Beg+.02,0,Int.Beg+.02,height,lty=1,lwd=3,col=cola)


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

#print("DRAW")
#print("DRAW")
#print("DRAW")
#print(ncp)
#print(zsds)
#print(w.all)

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

bar.width = .01
D.X = seq(x.lim.min,x.lim.max,bar.width) 

D.Y = build.dens(D.X, ncp, zsds, df, CURVE.TYPE,Folded=Folded) 
D.Y = D.Y / (rowSums(D.Y) * bar.width)
dim(D.Y)
D.Y = colSums(D.Y * w)


d.dense.sel = sum(D.Y[D.X > x.lim.min & D.X >= Int.Beg]*bar.width)
d.dense.sel

d.hist.sel = mean(as.numeric(val.input > x.lim.min & val.input > Int.Beg & val.input < Int.End)) /
	mean(as.numeric(val.input > x.lim.min & val.input < Int.End))
d.hist.sel

check = D.Y[D.X > crit & D.X < Int.Beg]
length(check)
EPJS = sum(D.Y[D.X > crit & D.X < Int.Beg])/sum(D.Y[D.X > crit])
EPJS

OFJS = sum(as.numeric(val.input > crit & val.input < Int.Beg))
OFJS

OFS = sum(as.numeric(val.input > crit & val.input < Int.End))
OFS

OPJS = OFJS/OFS
OPJS

scale = d.hist.sel/d.dense.sel;scale

lines(D.X[which(D.X >= x.start & D.X < x.end)],
	D.Y[which(D.X >= x.start & D.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max),ylim=c(0,ymax))


if (Show.Gamma & !is.na(results$gamma_shape) & !is.na(results$gamma_rate)) {

     

   D.Y.g = Get.Gamma.Density(D.X,shape = results$gamma_shape,rate = results$gamma_rate)
   D.Y.g <- D.Y.g / sum(D.Y.g * bar.width)

   d.dense.gamma.sel = sum(D.Y.g[D.X > x.lim.min & D.X >= Int.Beg]*bar.width)
   d.dense.gamma.sel

   scale.gamma = d.hist.sel/d.dense.gamma.sel;scale.gamma

   #plot(D.X,D.Y.g)

   lines(D.X[which(D.X >= x.start & D.X < x.end)],
	  D.Y.g[which(D.X >= x.start & D.X < x.end)]*scale.gamma,
	  lty=Ltype,col=col.gamma,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max),ylim=c(0,ymax))


}

boundary = Int.Beg
if (round(Int.Beg,2) == 1.96) {boundary = 2}

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

#z = val.input

# --- Main analysis function ---
run_ML_gamma <- function(
  z,
  Int.Beg = 1.96,
  Int.End = Inf,
  crit = 1.96,
  n_starts = 1,
  start_shape = NULL,
  start_rate = NULL,
  two.sided = TRUE
) {

  
  p_extreme_obs <- mean(z >= Int.End)
  p_sig_trunc_obs <- mean(z >= crit & z < Int.End)
  p_all_trunc_obs <- mean(z < Int.End)

  # keep only selected range
  z_obs <- z[z >= Int.Beg & z < Int.End]

  # probability of observed absolute z falling in selected interval
  prob_select_ncp <- function(ncp, Int.Beg, Int.End) {
    if (two.sided) {
      # P(Int.Beg <= |Z| < Int.End | ncp)
      upper <- if (is.infinite(Int.End)) 1 else pnorm(Int.End, mean = ncp, sd = 1)
      lower <- pnorm(Int.Beg, mean = ncp, sd = 1)

      pos <- upper - lower

      upper_neg <- pnorm(-Int.Beg, mean = ncp, sd = 1)
      lower_neg <- if (is.infinite(Int.End)) 0 else pnorm(-Int.End, mean = ncp, sd = 1)

      neg <- upper_neg - lower_neg

      pos + neg
    } else {
      # one-sided positive z selection
      upper <- if (is.infinite(Int.End)) 1 else pnorm(Int.End, mean = ncp, sd = 1)
      lower <- pnorm(Int.Beg, mean = ncp, sd = 1)
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
          Int.Beg = lo,
          Int.End = hi
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
  # [Int.Beg, Int.End). This is the likelihood denominator.
  prob_select <- function(shape, rate) {
    prob_select_range(
      shape = shape,
      rate = rate,
      lo = Int.Beg,
      hi = Int.End
    )
  }

  compute_edr_gamma_trunc <- function(shape, rate, zmax = Int.End) {
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

      # Numerator: P(crit <= |Z| < Int.End)
      p_sig <- prob_select_range(
        shape = shape,
        rate = rate,
        lo = crit,
        hi = Int.End
      )

      # Denominator: P(0 <= |Z| < Int.End)
      p_all <- prob_select_range(
        shape = shape,
        rate = rate,
        lo = 0,
        hi = Int.End
      )

      p_sig / p_all
    }


    compute_err_gamma <- function(shape, rate) {

      ncp_upper <- get_ncp_upper(shape, rate)

      numerator <- integrate(
        function(ncp) {
          p_sig_ncp <- prob_select_ncp(
            ncp = ncp,
            Int.Beg = crit,
            Int.End = Int.End
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
        hi = Int.End
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
#val.input = abs(c(rnorm(500,0,1),rnorm(100,2,1)));hist(val.input)
#cluster.id = NULL;yi=NULL;sei=NULL
#################################################
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
#################################################

#print(paste0("Title: ",Title))

#cluster.id = NULL


skip_zcurve = FALSE

if (length(p) > 0) val.input = qnorm(1-p/2)

if (!is.null(yi) > 0) val.input = yi/sei


if (CURVE.TYPE == "z") { 
  crit = qnorm(alpha/2,lower.tail=FALSE)
 } else {
  crit = qt(alpha/2,df,lower.tail=FALSE)
}


### remove NA
missing = is.na(val.input)
val.input = val.input[!missing]
if (!is.null(cluster.id)) cluster.id = cluster.id[!missing]
if (!is.null(yi)) yi = yi[!missing]
if (!is.null(sei)) sei = sei[!missing]

check = c(length(val.input),length(cluster.id),length(yi),length(sei))
check = check[check > 0]
check = check == sum(check)/length(check)
if(mean(check) < 1) stop("Check data: unequal number of cases")

#if (two.sided) val.input = abs(val.input)

### create set with z-scores in the interval used for model fitting

if (length(val.input[val.input > Int.Beg & val.input < Int.End]) > 10) {
  print("!START!") 
} else {
  print("Insufficient Data")
  skip_zcurve = TRUE
}    

n.z = length(val.input);
#print(n.z);return(n.z)

n.z.sig = length(val.input[abs(val.input) > crit])

pow = pnorm(val.input,crit) + pnorm(-crit, val.input)
MOP = mean(pow)

###

extreme = c(
	length(val.input[val.input < -Int.End])/length(val.input[val.input < -crit]),
	length(val.input[val.input > Int.End])/length(val.input[val.input > crit])
)
names(extreme) = c("Ext.Neg","Ext.Pos")
(extreme*100)

components = length(ncp)

### create set with z-scores in the interval used for model fitting
#INT = val.input[val.input >= Int.Beg & val.input <= Int.End]


slope.width = 1
slope.int = val.input[val.input >= Int.Beg + 2 * bw.est & val.input < Int.Beg + 2 * bw.est + slope.width]
slope.k = length(slope.int)
slope = rep(NA,6)
if (slope.k > 10) { 
     s = get_slope(INT = slope.int, bw=bw.est, from=Int.Beg, to = Int.Beg + slope.width)
     reliability <- s["slope"]^2 / (s["slope"]^2 + s["se"]^2)
     adjusted.slope <- reliability * s["slope"]
     slope = c(slope.k,adjusted.slope,s)
}
names(slope) = c("slope.k","slope.adj","slope","slope.se","slope.lb","slope.ub")
#print("Slope check")
#print(slope)

g_ncp = ncp
g_zsds = zsds
g_ncp_fixed = TRUE
g_zsds_fixed = TRUE
#ncp = g_ncp;zsds=g_zsds;NCP_FIXED = g_ncp_fixed;ZSDS_FIXED = g_zsds_fixed

heterogeneity = rep(NA,3)

if(!skip_zcurve) {

### Heterogeneity Estimation
if(Heterogeneity.Test) {
  print("Heterogeneity Test")
  heterogeneity.test = run_zcurve(Est.Method="OF",val.input=val.input,crit=crit,
     Int.Beg = Int.Beg,Int.End=Int.End,ncp = 2.5, zsds = 2,
     NCP_FIXED=FALSE,ZSDS_FIXED=FALSE,W_FIXED=TRUE,
     cluster.id = cluster.id, boot.iter = 500,Run.Gamma=FALSE,
     yi=NULL,sei=NULL,Folded=Folded)
  heterogeneity = sqrt(heterogeneity.test$zsds^2 - 1)
  #print("Hetero OK")
  #print(heterogeneity)
  if(length(heterogeneity) < 3) break

  if (mean(heterogeneity) < min.heterogeneity) { 
    print("QuickTest")
    quick.test = run_zcurve(Est.Method="OF",val.input=val.input,crit=crit,Int.Beg=Int.Beg,Int.End=Int.End,
      cluster.id=cluster.id,boot.iter = 500,ncp=2.5,zsds=1,
      NCP_FIXED=FALSE,ZSDS_FIXED=TRUE,Run.Gamma=FALSE, yi=NULL, sei=NULL,
      Folded=Folded)
    quick.test$ncp #quick.test$ncp[1] = 1.5
    if (all(abs(ncp - quick.test$ncp[1]) > .1)) {
      ncp_val = round(quick.test$ncp[1], 2)
      ncp = sort(unique(c(0, seq(ncp_val, 0, -1), seq(ncp_val, 6, 1))))
      components = length(ncp)
      zsds = rep(1, components)
    }
}

if (heterogeneity[3] < .2) { # high confidence that data are homo
  ncp = quick.test$ncp[1]
  zsds = 1
  NCP_FIXED=FALSE
  ZSDS_FIXED=TRUE
}


} # EOF Heterogeneity Test

ncp
zsds
NCP_FIXED
ZSDS_FIXED

sel.width = .5
k.sel = length(val.input[val.input > Int.Beg & val.input < Int.Beg + sel.width])
#k.sel = 1000
edr.factor <- max(.01,1/sqrt(k.sel)) 
edr.factor
edr.adj.ci = 1 * edr.factor 
edr.adj.ci = max(CI.EDR.MIN.ADJ,edr.adj.ci)
edr.adj.ci

#}
### BBB 

#res = zcurve.res
#boot.iter = 500

print("MAIN")
print(Folded)
#Est.Method = "OF"
print(Est.Method)
res.main = run_zcurve(Est.Method = Est.Method,val.input = val.input,
                 crit=crit,Int.Beg=Int.Beg,Int.End=Int.End,
                 cluster.id=cluster.id,boot.iter = boot.iter,ncp = ncp, zsds = zsds, 
                 NCP_FIXED=NCP_FIXED,ZSDS_FIXED=ZSDS_FIXED,edr_adj_ci=edr.adj.ci,
                 Run.Gamma=Run.Gamma,yi=yi,sei=sei,Folded=Folded,
                 x.lim.min=x.lim.min,x.lim.max=x.lim.max  )

#res.main$w.all
#print(res.main$local_es)
#print(res.main$es_mean_all)
#print(res.main$es_median_all)
#print(res.main$es_tau)

### adjustments for difficult cases 

adjustments = rep(NA,5)

if (Est.Method != "ML_gamma" & boot.iter > 0) {

  #res$EDR;res$ERR;res$ncp;res$w.inp

  #edr.adj.ci;res$EDR

  #res$EDR[2] <- res$EDR[2] - edr.adj.ci
  #res$EDR[3] <- res$EDR[3] + edr.adj.ci

  #print(slope)

  if(!is.na(slope[1])) { 
    slope.sign = sign(slope[2]) 
    slope.adj = slope["slope.adj"]
    if(slope.sign < 0) {
       if(slope.adj < -1) slope.adj = -1
       slope.adj = -1 - slope.adj
       slope.adj = slope.adj * edr.factor
       #res$EDR[3] = res$EDR[3] + abs(slope.adj)
     } else {
       if(slope.adj > 3) slope.adj = 1
       slope.adj = slope.adj * edr.factor  
       #res$EDR[2] = res$EDR[2] - abs(slope.adj)
     }
  } else {
    slope.sign = NA
    slope.adj = NA
    #res$EDR[2] = res$EDR[2] - .15
    #res$EDR[3] = res$EDR[3] + .15
  }

  #if(res$EDR[2] < alpha) res$EDR[2] = alpha
  #if(res$EDR[3] > 1) res$EDR[3] = 1

  ###

  err.factor = .6
  sel.width = 1
  k.sel = length(val.input[val.input > Int.Beg & val.input < Int.Beg + sel.width])
  err.adj.ci <- 1/sqrt(k.sel) * err.factor

  err.adj.ci = max(CI.ERR.MIN.ADJ,err.adj.ci)

  #print(paste0("ERR adj: ",err.adj.ci))

  #print(res$ERR)

  #res$ERR[2] <- res$ERR[2] - err.adj.ci
  #res$ERR[3] <- res$ERR[3] + err.adj.ci

  #if(res$ERR[2] < alpha/2) res$ERR[2] = alpha / 2
  #if(res$ERR[3] > 1) res$ERR[3] = 1

  #print(res$ERR)

  ###

  #res$EDR;res$ERR

  adjustments = c(edr.factor,edr.adj.ci,slope.adj*edr.factor,err.factor,err.adj.ci)

  }


FDR = as.numeric((1/res.main$EDR - 1)*(alpha/(1-alpha)))
FDR = FDR[c(1,3,2)]
names(FDR) = c("FDR_pe","FDR_lb","FDR_ub")

POW_h1_sig <- as.numeric((res.main$ERR - FDR * alpha) / (1 - FDR))
names(POW_h1_sig) = c("POW_H1_SIG_pe","POW_H1_SIG_lb","POW_H1_SIG_ub")

ODR_EDR_D = as.numeric(res.main$ODR_EDR_D)
names(ODR_EDR_D) = c("ODR_EDR_D_pe","ODR_EDR_D_lb","ODR_EDR_D_ub")



} # skip zcurve

#skip_zcurve

#res.main

###

#skip_zcurve

if(skip_zcurve) {

  res.main = list(
   ODR = rep(NA,3),
   EDR = rep(NA,3),
   ERR = rep(NA,3),
   ncp = matrix(NA,3,7),
   zsds = matrix(NA,3,7),
   w.inp = matrix(NA,3,7),
   w.all = matrix(NA,3,7),
   local_power = matrix(NA,3,12),
   local_es = matrix(NA,3,12),
   fit = rep(NA,3),
   gamma_shape = rep(NA,3),
   gamma_rate = rep(NA,3),
   gamma_edr = rep(NA,3),
   gamma_err = rep(NA,3),
   shape_d_median = rep(NA,3),
   shape_d_mean   = rep(NA,3),
   es_mean_all    = rep(NA,3),
   es_median_all  = rep(NA,3),
   es_mean_sig    = rep(NA,3),
   es_tau         = rep(NA,3)
   )

POW_h1_sig = rep(NA,3)
FDR = rep(NA,3)
ODR_EDR_D = rep(NA,3)
if(Heterogeneity.Test) heterogeneity = rep(NA,3)
adjustments = rep(NA,5)

} # EOF skip

names(ODR_EDR_D) = c("ODR_EDR_D_pe","ODR_EDR_D_lb","ODR_EDR_D_ub")

res.main$ODR = as.numeric(res.main$ODR)
res.main$EDR = as.numeric(res.main$EDR)
res.main$ERR = as.numeric(res.main$ERR)

names(res.main$ODR) = c("ODR_pe","ODR_lb","ODR_ub")
names(res.main$EDR) = c("EDR_pe","EDR_lb","EDR_ub")
names(res.main$ERR) = c("ERR_pe","ERR_lb","ERR_ub")

names(res.main$shape_d_mean) = c("shape_mean_pe","shape_mean_lb","shape_mean_ub")
names(res.main$shape_d_median) = c("shape_median_pe","shape_median_lb","shape_median_ub")

names(res.main$es_mean_all)   = c("es_mean_all_pe","es_mean_all_lb","es_mean_all_ub")
names(res.main$es_median_all) = c("es_median_all_pe","es_median_all_lb","es_median_all_ub")
names(res.main$es_mean_sig)   = c("es_mean_sig_pe","es_mean_sig_lb","es_mean_sig_ub")
names(res.main$es_tau)        = c("es_tau_pe","es_tau_lb","es_tau_ub")



if(Heterogeneity.Test) {
  heterogeneity = as.numeric(heterogeneity)
  names(heterogeneity) = c("ncp_sd_pe","ncp_sd_lb","ncp_sd_ub")
}

names(adjustments) = c("edr.factor","edr.adj.ci","slope.adj","err.factor","err.adj")


if (length(res.main$local_power) > length(ncp)) {
   colnames(res.main$local_power) = paste0("local_power",1:ncol(res.main$local_power))
   rownames(res.main$local_power) = c("pe","lb","ub")
}

if (length(res.main$local_es) > components) {
   colnames(res.main$local_es) = paste0("local_es",1:ncol(res.main$local_es))
   rownames(res.main$local_es) = c("pe","lb","ub")
}


if(length(ncp) > 1) {
  colnames(res.main$ncp) = paste0("ncp",1:length(ncp))
  rownames(res.main$ncp) = c("pe","lb","ub")

  colnames(res.main$zsds) = paste0("zsds",1:length(ncp))
  rownames(res.main$zsds) = c("pe","lb","ub")

  colnames(res.main$w.inp) = paste0("w.inp",1:length(ncp))
  rownames(res.main$w.inp) = c("pe","lb","ub")

  colnames(res.main$w.all) = paste0("w.all",1:length(ncp))
  rownames(res.main$w.all) = c("pe","lb","ub")
}


n.z.clusters = NULL
if (!is.null(cluster.id)) n.z.clusters <- length(unique(cluster.id))


results = list(
		clusters       = n.z.clusters,
		slope          = slope,
		ODR            = res.main$ODR,
		EDR            = res.main$EDR,
		ERR            = res.main$ERR,
		FDR            = FDR,
		POW_h1_sig     = POW_h1_sig,
		ncp            = res.main$ncp,
		zsds           = res.main$zsds,
		heterogeneity  = heterogeneity,
		w.inp          = res.main$w.inp,
		w.all          = res.main$w.all,	
		local_power    = res.main$local_power,
		local_es       = res.main$local_es,
		ODR_EDR_D      = ODR_EDR_D,
		shape_d_median = res.main$shape_d_median,
		shape_d_mean   = res.main$shape_d_mean,
		adjustments    = adjustments,
		fit            = res.main$fit,         
		gamma_shape    = res.main$gamma_shape,
		gamma_rate     = res.main$gamma_rate,
		gamma_edr      = res.main$gamma_edr,
		gamma_err      = res.main$gamma_err,
         es_mean_all    = res.main$es_mean_all,
         es_median_all  = res.main$es_median_all,
         es_mean_sig    = res.main$es_mean_sig,
         es_tau         = res.main$es_tau
      )

#fff

#results
#results$w.all
#summary(val.input)

##########################################
### This Code is Used to Create Graphic
##########################################

###sss

if (Show.Histogram & !skip_zcurve) {

     par("mar") 
     old_mar <- par("mar")
     if(length(results$local_power) > 0) par(mar = old_mar[1] + c(2.5,0,0,0) )
     if(length(results$local_es) > 0)    par(mar = old_mar[1] + c(3.5,0,0,0) )
     par("mar")
  
	if (boot.iter == 0) {
		#print("Draw Histogram NO CI")
		Draw.Histogram.3.8(results$w.all[1,],cola=col.hist,Write.CI = FALSE)
	} else {
		#print("Draw Histogram WITH CI")
		Draw.Histogram.3.8(results$w.all[1,],cola=col.hist,Write.CI = TRUE)
	}

	if (Show.Curve.All) { 

		Draw.Curve.All(results,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max,Show.Gamma=Show.Gamma)
		Draw.Curve.All(results,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End,Show.Gamma=Show.Gamma)

	} # EOF Show.Curve.All


	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)


	if (length(results$local_power) > 0) { 
       if(!is.null(yi)) {
          Write_Local_Power_ES(results$local_power[1,],results$local_es[1,])	
       } else {   
          Write_Local_Power(results$local_power[1,])
       }	
    } 

    par(mar = old_mar)



} # End of Show.Histogram	for bootstrap


### reset to global specifications
ncp  = g_ncp
zsds = g_zsds
NCP_FIXED = g_ncp_fixed
ZSDS_FIXED = g_zsds_fixed

 
return(results)

} #End of Zing


######################################################################
######################################################################
######################################################################


