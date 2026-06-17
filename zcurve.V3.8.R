

#rm(list = ls())

########################################################################
### SETTING PARAMETERS FOR Z-CURVE MODEL
########################################################################

version <- "zcurve3(3.8)"
date <- "2026.06.13"  # Version label to appear on plots
note <- "Now with Effect Size Estimation !!!"

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

int.loc <- 0.5               # Plot local power intervals below x-axis (set 0 to disable)
hist.bar.width <- 0.2        # Width of histogram bars

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

ZSDS.FIXED <- TRUE                  # Fix SD values for EXT method 
NCP.FIXED <- TRUE                   # Fix non-central parameter(NCP) means values for EXT method
W.FIXED   <- FALSE                  # Fix weights for EXT method

### DENSITY-BASED SETTINGS (Only used with Est.Method = "OF")

n.bars <- 512                       # Number of bars in histogram

Augment <- TRUE                     # Apply correction for bias at lower bound

bw.est  <- 0.05                     # Bandwidth for kernel density (lower = less smoothing, higher = more smoothing)
bw.draw <- 0.05   		        	 # Bandwith of kernel density in plot
spline.width <- .5                  # blending of truncated and normal KD

### INPUT RESTRICTIONS

MAX.INP.Z <- 100                    # Optionally restrict very large z-values (set Inf to disable)

### CONFIDENCE INTERVALS / BOOTSTRAPS

boot.iter <- 0                      # Number of bootstrap iterations (suggest 500+ for final models)

CI.ALPHA <- 0.05                    # CI level (default = 95%)
CI.EDR.MIN.ADJ = .05                # Minimum Adjustmnet for EDR CI
CI.ERR.MIN.ADJ = .03                # Minimum Adjustment for ERR CI
ODR_EDR_D_adj = .05                 # Adjustment for ODR > EDR bias test  	


TEST4BIAS <- TRUE                   # Enable optional bias test

### DISPLAY FINAL STATUS

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

#odr = res$ODR[1];edr = res$EDR
#init_mu = NULL;init_tau = NULL
#selection = "one.sided"
fit_fixed_weight_selection <- function(yi, sei, odr, edr,
                                       crit = 1.96,
                                       selection = "one.sided",
                                       init_mu = NULL,
                                       init_tau = NULL) {


  w_from_odr_edr <- function(odr, edr, eps = 1e-6) {
    odr <- pmin(pmax(odr, eps), 1 - eps)
    edr <- pmin(pmax(edr, eps), 1 - eps)
  
    # One-sided publication bias only: if EDR >= ODR, set w = 1.
    edr <- pmin(edr, odr)
  
   w <- edr * (1 - odr) / (odr * (1 - edr))
    pmin(pmax(w, 0), 1)
  }

  w = w_from_odr_edr(odr,edr)

  
  stopifnot(length(yi) == length(sei))
  
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi <- yi[keep]
  sei <- sei[keep]
  
  eps <- 1e-10
  w <- pmin(pmax(w, eps), 1)
  
  if (is.null(init_mu)) {
    init_mu <- sum(yi / sei^2) / sum(1 / sei^2)
  }
  
  if (is.null(init_tau)) {
    init_tau <- sd(yi)
    if (!is.finite(init_tau) || init_tau <= 0) init_tau <- 0.10
  }

  ###
  get_estimates = function(w) {

    negloglik <- function(par) {
      mu <- par[1]
      tau <- exp(par[2])
      sd_i <- sqrt(sei^2 + tau^2)
  
      if (selection == "one.sided") {
        # Assumes yi is sign-coded so positive values are in the predicted direction.
        sig <- yi / sei > crit
      
        p_sig <- 1 - pnorm(crit * sei, mean = mu, sd = sd_i)
        p_nsig <- 1 - p_sig
      
      } else {
        # Two-sided selection based on |z| > crit.
        sig <- abs(yi / sei) > crit
      
        p_lower <- pnorm(-crit * sei, mean = mu, sd = sd_i)
        p_upper <- 1 - pnorm(crit * sei, mean = mu, sd = sd_i)
        p_sig <- p_lower + p_upper
        p_nsig <- 1 - p_sig
      }
    
      A <- p_sig + w * p_nsig
      A <- pmax(A, eps)
    
      log_weight <- ifelse(sig, 0, log(w))
    
      ll <- dnorm(yi, mean = mu, sd = sd_i, log = TRUE) +
        log_weight -
        log(A)
    
      -sum(ll)
    }

    fit <- optim(
      par = c(init_mu, log(init_tau)),
      fn = negloglik,
      method = "BFGS",
      hessian = TRUE,
      control = list(maxit = 1000)
    )
  
    mu_hat <- fit$par[1]
    tau_hat <- exp(fit$par[2])

    ret = c(w,mu_hat,tau_hat)

  }

  ret = get_estimates(w)
  
}





run_bias_tests <- function()  {

ODR = length(val.input[val.input > crit])/length(val.input)

Credibility = c(round(ODR*100,2),round(MOP*100,2),
    round(pbinom(n.z-n.z*ODR,n.z,1-MOP),6))
names(Credibility) = c("ODR","MOP","Credibility_p_value")
#Credibility
#ODR;MOP
Rindex = c(round(ODR*100,2),round(min(1,(2*MOP - ODR))*100,2),
    round(pbinom(n.z-n.z*ODR,n.z,1-min(1,(2*MOP - ODR))),6))
names(Rindex) = c("ODR","Rindex","Rindex_p_value")

delta  <- 0.2 *crit   # 
k.above <- sum(val.input > crit & val.input <= crit + delta)  # just significant
k.below <- sum(val.input >= crit - delta & val.input < crit)   # just nonsignificant
# Under H0 (no bias), these should be ~equal
# Binomial test with p = 0.5
Caliper.p = NA
if (k.above > 0 & k.below > 0) Caliper.p <- round(binom.test(k.above, k.above + k.below, 
               p = 0.5, alternative = "greater")$p.value,6)
Caliper = c(k.above,k.below,Caliper.p)
names(Caliper) = c("k(above)","k(below)","Caliper_p_value")

ZBIAS = bias.test(val.input,caliper.width = .8)

  ## rdd bias test

  z = val.input
  bw = 2 # slope.int 
  bin.width = .20

  
  # --- Step 1: Bin the z-values ---
  # Create bins that don't straddle the cutoff
  breaks.below <- seq(crit - bw, crit, by = bin.width)
  breaks.above <- seq(crit, crit + bw, by = bin.width)
  breaks <- unique(c(breaks.below, breaks.above))
  
  h <- hist(z[z >= (crit - bw) & z <= (crit + bw)], 
            breaks = breaks, plot = FALSE)
  
  mids   <- h$mids
  counts <- h$counts
  
  # --- Step 2: Set up RDD regression ---
  # Running variable: distance from cutoff
  x <- mids - crit
  
  # Treatment: above the cutoff
  D <- as.numeric(mids >= crit)
  
  # Interaction model: allows different slopes on each side
  # counts ~ x + D + x:D
  # Coefficient on D = the jump at the cutoff
  fit <- lm(counts ~ x * D)
  
  # --- Step 3: Extract the discontinuity estimate ---
  coefs <- summary(fit)$coefficients
  
  jump   <- coefs["D", "Estimate"]
  se     <- coefs["D", "Std. Error"]
  t.stat <- coefs["D", "t value"]
  
  # One-sided test: bias means jump > 0 (more above than below)
  p.value <- pt(t.stat, df = fit$df.residual, lower.tail = FALSE)
  
  rdd_res <- c(jump = jump, se = se, p.RDD = p.value)
  names(rdd_res) <- c("jump", "se", "RDD_p_value")


  ### hybrid_bias_test <- function(z, crit = qnorm(0.975), band = 0.5, bw = 0.20) {
  # z:    all z-values (significant and non-significant)
  # crit: significance threshold
  # band: width of band on each side of crit
  # bw:   kernel bandwidth for density estimation on right side

  band = 1
  
  z.sig <- z[z >= crit]
  z.ns  <- z[z < crit]
  
  # --- Right side: density estimation with truncated kernel ---
  z.shifted <- z.sig - crit  # shift so boundary is at 0
  
  # Evaluate density at grid points in [0, band]
  grid <- seq(0, band, length.out = 50)
  step <- grid[2] - grid[1]
  d.grid <- numeric(length(grid))
  
  for (i in seq_along(grid)) {
    raw <- dnorm(z.shifted, grid[i], bw)
    tc  <- pnorm(grid[i] / bw)
    if (tc > 0.001) {
      d.grid[i] <- mean(raw) / tc
    }
  }
  
  # Density at crit from right side
  d.at.crit <- d.grid[1]
  
  # Average density in right band (captures slope)
  d.avg.right <- mean(d.grid)
  
  # --- Left side: uniform assumption ---
  # Under uniform, density is constant everywhere below crit
  # Under H0 (continuity), that constant = d.at.crit
  # So average density in left band = d.at.crit
  d.avg.left <- d.at.crit
  
  # --- Observed counts in bands ---
  k.above <- sum(z.sig >= crit & z.sig < crit + band)
  k.below <- sum(z.ns >= (crit - band) & z.ns < crit)
  k.total <- k.above + k.below
  
  # --- Binomial test ---
  # Under H0: probability of being in left band
  # = avg density left / (avg density left + avg density right)
  # This adjusts for slope on the significant side
  p.below.exp <- d.avg.left / (d.avg.left + d.avg.right)
  
  # One-sided: fewer below than expected = publication bias
  p.value <- pbinom(k.below, k.total, p.below.exp)
  
  hybrid <- c(k.below/(k.below + k.above), p.exp = round(p.below.exp, 4), 
              p.hybrid = p.value)

  names(hybrid) = c("OBS", "EXP","hybrid_p_value")

  hybrid


BIAS.RESULTS = list(
   Credibility = Credibility,
   Rindex      = Rindex,
   ZBIAS       = ZBIAS,
   Caliper     = Caliper,
   RDD         = rdd_res,
   Hybrid      = hybrid
)

#print("EOF BIAS TEST")
#print(BIAS.RESULTS)

BIAS.RESULTS

return(BIAS.RESULTS)

}




#####################################################################
### Z-curve 3.0 - NEW EM Estimation for Extended (free w, ncp, zsds)
#####################################################################


##############################################
### Get Densities
##############################################

build.dense <- function(D.X, ncp, zsds = NULL, df = NULL, CURVE.TYPE = "z") {
  
  xmat  <- matrix(rep(D.X, each = length(ncp)), nrow = length(ncp))
  mumat <- matrix(rep(ncp, times = length(D.X)), nrow = length(ncp))
  
  if (CURVE.TYPE == "z") {
    sdmat <- matrix(rep(zsds, times = length(D.X)), nrow = length(ncp))
    
    Dens <- dnorm(xmat, mean = mumat, sd = sdmat) +
            dnorm(-xmat, mean = mumat, sd = sdmat)
    
  } else {
    Dens <- dt(xmat, df = df, ncp = mumat) +
            dt(-xmat, df = df, ncp = mumat)
  }
  
  Dens
}


###

Get.Densities = function(INT, bw = 0.20, from = 1.96, to = 6, width = 1, Augment = TRUE) {

  #from = 1.96; to = 6;bw = .05

  n = 401
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

} # EOF `nsities 

#######################################################
### End of Get Densities
#######################################################




#####################################################################
### Compute Power Function (New: Discrete & Continuous
#####################################################################

Compute.Power.Z.General = function(cp.input, Int.Beg, Int.End) {

  Compute.Power.Z.Discrete = function(cp.input, Int.Beg, Int.End) {

  w.inp      = cp.input$w.inp
  ncp        = cp.input$ncp
  components = length(ncp)

  # Extreme value correction
  ext.inp = length(val.input[val.input > Int.End]) /
            length(val.input[val.input > Int.Beg])

  # Component power (two-tailed, with sign error)
  pow.dir = pnorm(abs(ncp), crit)
  pow     = pow.dir + (1 - pnorm(ncp, -crit))

  # Power conditional on selection interval
  pow.sel = pnorm(ncp, Int.Beg) + pnorm(-ncp, Int.Beg)

  # Extend vectors with extreme component
  pow.ext     = c(pow,     1)
  pow.dir.ext = c(pow.dir, 1)
  pow.sel.ext = c(pow.sel, 1)
  w.inp.ext   = c(w.inp * (1 - ext.inp), ext.inp)

  # Correct for selection: w.all = w.inp / pow.sel, normalized
  w.inp.ext   = c(w.inp * (1 - ext.inp), ext.inp)
  pow.sel.ext = c(pow.sel, 1)

  w.all.ext   = (w.inp.ext / pow.sel.ext)
  w.all.ext   = w.all.ext / sum(w.all.ext)

  w.all = w.all.ext[1:components]/(1-ext.inp)

  # EDR: average power over all studies (pre-selection)
  EDR = sum(w.all.ext * pow.ext)

  # ERR: average directional power over significant studies
  w.sig.ext = w.all.ext * pow.ext
  w.sig.ext = w.sig.ext / sum(w.sig.ext)

  w.sig = w.sig.ext[1:components]/(1-w.sig.ext[components+1])

  ERR = sum(w.sig.ext * pow.dir.ext)

  ERR = min(max(ERR, alpha / 2), 1)
  EDR = min(max(EDR, alpha),     1)


  if (1 == 2) { ### Don't run

  if (components > 4 ) {
     if (ncp[2] < 2 & ncp[1] == 0) { 

     drop_idx = 2  
     lo <- drop_idx - 1
     hi <- drop_idx + 1
    
     w.new = w.all.ext
     w.new[lo] <- w.new[lo] + w.new[drop_idx] * (pow[hi] - pow[drop_idx]) / (pow[hi] - pow[lo])
     w.new[hi] <- w.new[hi] + w.new[drop_idx] * (pow[drop_idx] - pow[lo]) / (pow[hi] - pow[lo])
     w.new[drop_idx] <- 0
  
     FDR2 <- (w.new[1] * alpha) / sum(w.new * pow.ext)

     h1_idx <- which(ncp != 0)
     POW_h1_sig <- sum(w.new[h1_idx] * pow[h1_idx]^2) / sum(w.new[h1_idx] * pow.ext[h1_idx]) 
     
   }}

   } ### End of Don't Run

#####

  local.power = NULL

  if (int.loc > 0) {

  X     = seq(x.lim.min, Int.End, by = 0.01)

  # Mixture density at each grid point
  lik.mat = outer(X, ncp, function(x, mu) dnorm(x, mu))
  wd.X    = lik.mat %*% w.all  # column vector, length(X)

  # Local power: posterior-weighted average of component power
  loc.p   = (lik.mat %*% (w.all * pow)) / wd.X

  # Bin into intervals, density-weighted average within each bin
  int  = seq(x.lim.min, x.lim.max, by = int.loc)
  bins = cut(X, breaks = int, include.lowest = TRUE)

  local.power = as.vector(tapply(seq_along(X), bins, function(idx) {
    sum(loc.p[idx] * wd.X[idx]) / sum(wd.X[idx])
  }))

  midpoints           = (int[-length(int)] + int[-1]) / 2
  names(local.power)  = paste0("lp.", midpoints)

} ### end local power

####

  out = list(
    EDR         = EDR,
    ERR         = ERR,
    w.all       = w.all,
    local.power = local.power
  )

out

} ### EOF Compute.Power.Z.Discrete


Compute.Power.Z.Continuous = function(cp.input, Int.Beg, Int.End) {

  deci   = 2
  zx.bw  = 1 / (10^deci)
  zx     = seq(0, Int.End, zx.bw)

  components = length(cp.input$ncp)
  ncp        = round(cp.input$ncp, deci)
  ncp.sd     = sqrt(cp.input$zsds^2 - 1)

  # Floor ncp.sd at grid step to prevent underflow
  ncp.sd = max(zx.bw, ncp.sd)

  pow.sel = pnorm(Int.Beg, ncp, lower.tail = FALSE) + pnorm(-Int.Beg, ncp)
  pow.sel[pow.sel < .05] = .05

  w.all = cp.input$w.inp / pow.sel
  w.all = w.all / sum(w.all)

  # Extreme values correction
  ext.inp = length(val.input[val.input > Int.End]) /
            length(val.input[val.input > Int.Beg])

  # Power at each grid point
  pow.zx = pnorm(abs(zx), crit) + pnorm(-crit, abs(zx))

  # Mixture density over grid
  lik.mat = outer(zx, ncp, function(z, mu) dnorm(z, mu, ncp.sd))
  wd.mat  = lik.mat * matrix(w.all, nrow = length(zx), ncol = components, byrow = TRUE)
  wd.all  = rowSums(wd.mat)
  wd.all  = wd.all / sum(wd.all)

  wd.all.ext = wd.all * (1 - ext.inp)

  EDR       = sum(pow.zx * wd.all.ext) + ext.inp
  ERR.num   = sum(pow.zx * wd.all.ext * pow.zx) + ext.inp
  ERR.denom = sum(wd.all.ext * pnorm(abs(zx), crit)) + ext.inp
  ERR       = ERR.num / ERR.denom

  ERR = min(max(ERR, alpha / 2), 1)
  EDR = min(max(EDR, alpha),     1)

  ### local power

  if (int.loc == 0) {
      local.power = NULL
  } else {

  # Posterior weights: rows = grid points, cols = components
  row.sums            = rowSums(wd.mat)
  zero.rows           = row.sums == 0
  wd.mat[zero.rows, ] = matrix(w.all, nrow = sum(zero.rows),
                               ncol = components, byrow = TRUE)
  row.sums[zero.rows] = 1
  post.mat            = sweep(wd.mat, 1, row.sums, "/")

  # Posterior mean NCP for each (grid point, component) combination
  # Bayesian update: prior N(ncp_k, ncp.sd), likelihood z ~ N(mu, 1)
  mu_post_mat  = outer(zx, ncp, function(z, mu_k)
                   (z * ncp.sd^2 + mu_k) / (ncp.sd^2 + 1))

  # Power at each posterior mean NCP (two-tailed with sign error)
  pow_post_mat = pnorm(mu_post_mat - crit) + pnorm(-mu_post_mat - crit)

  # Local power: posterior-weighted average of posterior-mean power
  lp           = rowSums(post.mat * pow_post_mat)

  # Bin into intervals, simple mean within each bin
  int         = seq(0, Int.End, by = int.loc)
  bins        = cut(zx, breaks = int, include.lowest = TRUE)
  local.power = as.vector(tapply(lp, bins, mean))
  midpoints   = (int[-length(int)] + int[-1]) / 2
  names(local.power) = paste0("lp.", midpoints)

} ### end local power

  out = list(
    EDR         = EDR,
    ERR         = ERR,
    w.all       = w.all,
    local.power = local.power
  )

out

} ### EOF Compute.Power.Z.Continuous


  ### Main

  if (max(as.numeric(cp.input$zsds)) < 1.01) {
      out = Compute.Power.Z.Discrete(cp.input, Int.Beg, Int.End)
  } else {
      out = Compute.Power.Z.Continuous(cp.input, Int.Beg, Int.End)
  }


return(out)

} ### EOF Main


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


run_OF <- function(INT,Int.Beg,Int.End,ncp,zsds, cola = "springgreen4",bw.est = .2,
                    Dense.pre = NULL, D.X.pre = NULL, width = 1, 
                    NCP.FIXED = TRUE, ZSDS.FIXED = TRUE, W.FIXED = FALSE) {

  #Dense.pre = NULL; D.X.pre = NULL;width = spline.width

  components <- length(ncp)

  #from = Int.Beg; to = Int.End; width = 1
  ## ---- Observed density ----
  densy <- Get.Densities(INT, bw = bw.est, from = Int.Beg, 
                         to = Int.End, width = width, Augment = Augment)

  D.X   <- densy[, 1]
  O.D.Y <- densy[, 2]
  n.bars    <- length(D.X)
  bar.width <- D.X[2] - D.X[1]

  #plot(D.X,O.D.Y)

  use.precomputed <- !is.null(Dense.pre) &&
                     !is.null(D.X.pre) &&
                     identical(D.X, D.X.pre)

  ## ---- theta = [weights, ncp, zsds] — always ----
  n.wt <- components

  startval <- c(rep(1/n.wt, n.wt), ncp, zsds)

  lowlim   <- c(if (n.wt == 1) 1 else rep(0, n.wt),
                if (NCP.FIXED) ncp else rep(0,components),
                if (ZSDS.FIXED) zsds else rep(1,components)
			)
  highlim  <- c(rep(1, n.wt),
                if (NCP.FIXED) ncp else rep(6,components),
                if (ZSDS.FIXED) zsds else rep(6,components)
			)

  ## ---- Positions (always the same) ----
  idx.wt  <- 1:n.wt
  idx.ncp <- (n.wt + 1):(n.wt + components)
  idx.zsd <- (n.wt + components + 1):(n.wt + 2 * components)

  #theta = startval

  ## ---- Fitting function ----
  curve.fitting <- function(theta) {
    wt <- theta[idx.wt]
    wt <- wt / sum(wt)
 
    if (use.precomputed && NCP.FIXED && ZSDS.FIXED) {
      Dens.now <- Dense.pre
    } else {
      Dens.now <- build.dense(D.X, theta[idx.ncp], theta[idx.zsd], df, CURVE.TYPE)
      Dens.now <- Dens.now / (rowSums(Dens.now) * bar.width)
    }

    E.D.Y <- colSums(Dens.now * wt)
    #plot(D.X,E.D.Y)
    rmse = sqrt(mean((E.D.Y - O.D.Y)^2))

    scale = 1
    ## ---- Diagnostic plot ----
    if (TESTING && runif(1) > 0.9) {
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

  ## ---- Extract — positions are guaranteed ----
  WT <- auto$par[idx.wt]
  WT <- WT / sum(WT)

  out =  list(
    w.inp   = WT,
    ncp     = auto$par[idx.ncp],
    zsds    = auto$par[idx.zsd],
    fit     = auto$objective
  )

out

} ### EOF run_OF



##################################################
### ROOT EM FUNCTION
##################################################

run_EM <- function(
  INT,
  ncp,
  zsds,
  Int.Beg    = 1.96,
  Int.End    = max(INT),
  NCP.FIXED  = TRUE,
  ZSDS.FIXED = TRUE,
  W.FIXED    = FALSE,
  max_iter   = 200,
  tol        = 1e-6
) {


n_starts = 10

EM_EXT_FOLDED <- function(
  INT,
  ncp,
  zsds,
  w.inp,
  Int.Beg = 1.96,
  Int.End = max(INT),
  NCP.FIXED  = TRUE,
  ZSDS.FIXED = TRUE,
  W.FIXED    = FALSE,
  max_iter = 200,
  tol = 1e-6
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
  tol      = 1e-8,
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
    if (!W.FIXED) {
      w.inp <- n_k / k.int
      w.inp <- pmax(w.inp, 1e-12)
      w.inp <- w.inp / sum(w.inp)
    }
    if (!NCP.FIXED) {
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
    if (!ZSDS.FIXED) {
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

  components <- length(ncp)
  ncp_init   <- ncp
  zsds_init  <- zsds

  best_fit <- NULL
  best_ll  <- -Inf

  for (s in 1:n_starts) {

      raw <- rgamma(components, shape = 1)   # random Dirichlet
      w_s <- raw / sum(raw)

      fit <- EM_EXT_FOLDED(
        INT        = INT,
        ncp        = ncp,
        zsds       = zsds,
        w.inp      = w_s,
        Int.Beg    = Int.Beg,
        Int.End    = Int.End,
        NCP.FIXED  = NCP.FIXED,
        ZSDS.FIXED = ZSDS.FIXED,
        W.FIXED    = W.FIXED,
        max_iter   = max_iter,
        tol        = tol
       ) 

    if (!is.null(fit) && is.finite(fit$fit) && fit$fit > best_ll) {
      best_ll  <- fit$fit
      best_fit <- fit
      best_fit$start <- s   # track which start won
    }
  }

  return(best_fit)
}




##################################################
### Bootstrap
##################################################

### BOOTSTRAP FUNCTION FOR CI

run_bootstrap <- function(
                   val.input,
                   crit,   
                   ncp,
                   zsds,
                   Int.Beg,
                   Int.End,
                   Dense.pre = NULL,
                   D.X.pre = NULL,
                   width = 1,
                   boot.iter  = 500,
                   NCP.FIXED  = TRUE,
                   ZSDS.FIXED = TRUE,
                   W.FIXED    = FALSE,
                   cluster.id = NULL,
                   Est.Method = "OF"
                  ) {

    
      #NCP.FIXED = FALSE; ZSDS.FIXED = FALSE; W.FIXED = TRUE; ncp = 2.5; zsds = 2
      #cluster.id = NULL
      #boot.iter = 100
    
       
      INT = val.input[val.input >= Int.Beg & val.input <= Int.End]

      components = length(ncp)

      if (!is.null(cluster.id)) {

          print("USING CLUSTERS")

          ## ---- Cluster Bootstrap ----

          stopifnot(length(val.input) == length(cluster.id))	
          z.clusters <- split(val.input, cluster.id)
          z.cluster.names <- names(z.clusters)
          n.z.clusters <- length(z.cluster.names)

      } else {

         z.clusters = NULL; z.cluster.names = NULL; n.z.clusters = NULL

      }

      #Est.Method = "OF";width = 1
      #Est.Method = "EM"

      var.list = c("Est.Method","val.input","crit",
           "INT","cluster.id","z.clusters","z.cluster.names","n.z.clusters",
           "ncp", "zsds","Int.Beg", "Int.End", "Augment",
           "CURVE.TYPE","x.lim.min", "x.lim.max", 
           "NCP.FIXED", "ZSDS.FIXED","W.FIXED")


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
               sampled <- sample(z.cluster.names, n.z.clusters, replace = TRUE)
               boot_all <- unlist(z.clusters[sampled], use.names = FALSE)
               odr_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= Int.Beg & boot_all <= Int.End]
          } else {
               boot_all <- sample(val.input, length(val.input), replace = TRUE)
               odr_boot <- mean(boot_all >= crit)
               boot_sample <- boot_all[boot_all >= Int.Beg & boot_all <= Int.End]
          }

          fit <- run_EM(
              INT        = boot_sample,
              ncp        = ncp,
              zsds       = zsds,
              Int.Beg    = Int.Beg,
              Int.End    = Int.End,
              NCP.FIXED  = NCP.FIXED,
              ZSDS.FIXED = ZSDS.FIXED,
              W.FIXED    = W.FIXED)

          fit.res = list(
             ODR    = odr_boot,
             w.inp  = fit$w.inp,
             ncp    = fit$ncp,
             zsds   = fit$zsds,
             fit = fit$fit
             )


          fit.res

         }) # close boot_res

    } # EOF EM 

    if (Est.Method == "OF") {

         if (Est.Method == "OF") { 
           var.list = c(var.list,"run_OF","build.dense", "Get.Densities",
                  "bw.est","width","Dense.pre", "D.X.pre")
         } 

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

           fit <- tryCatch(
                 run_OF(
                 INT        = boot_sample,
                 Int.Beg    = Int.Beg,
                 Int.End    = Int.End,
                 ncp        = ncp,
                 zsds       = zsds,
                 Dense.pre  = Dense.pre,
                 D.X.pre    = D.X.pre,
                 bw.est     = bw.est, 
                 width      = width,
                 NCP.FIXED  = NCP.FIXED,
                 ZSDS.FIXED = ZSDS.FIXED,
                 W.FIXED    = W.FIXED
                 )
	         ,error = function(e) NULL )

          fit.res = list(
             ODR    = odr_boot,
             w.inp  = fit$w.inp,
             ncp    = fit$ncp,
             zsds   = fit$zsds,
             fit = fit$fit
             )

          fit.res
  
      }) # EOF boot_res

    } # EOF OF


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

#ncp = 2.5; zsds = 3;NCP.FIXED=FALSE;ZSDS.FIXED=FALSE;W.FIXED=FALSE
run_zcurve = function(Est.Method,val.input,crit,Int.Beg,Int.End,
           ncp,zsds,NCP.FIXED=TRUE,ZSDS.FIXED = TRUE,W.FIXED = FALSE,
           cluster.id,boot.iter=0,edr_adj_ci=0,Run.Gamma=FALSE,
           yi,sei) {

#NCP.FIXED = FALSE; ZSDS.FIXED = FALSE; W.FIXED = TRUE; ncp = 2.5; zsds = 2; Run.Gamma=FALSE
#yi = NULL;sei = NULL
#cluster.id = NULL

zcurve.time = system.time({

  INT = val.input[val.input >= Int.Beg & val.input <= Int.End]

  ODR = mean(val.input > crit)

  if(Est.Method != "EM") { 

      width = spline.width

      ## ---- Precompute Dens if both fixed ----
      if (NCP.FIXED && ZSDS.FIXED) {
        densy     <- Get.Densities(INT, bw = bw.est, from = Int.Beg,
                                 to = Int.End, width = width, Augment = Augment)
        D.X.pre   <- densy[, 1]
        Dense.pre  <- build.dense(D.X.pre, ncp, zsds, df, CURVE.TYPE)
        bar.width <- D.X.pre[2] - D.X.pre[1]
        Dense.pre  <- Dense.pre / (rowSums(Dense.pre) * bar.width)
      } else {
        D.X.pre  <- NULL
        Dense.pre <- NULL
      }

      ## ---- Point estimate ----
      res.run <- run_OF(
                 INT = INT, 
                 Int.Beg = Int.Beg,
                 Int.End = Int.End,   
                 ncp = ncp,
                 zsds = zsds,
                 Dense.pre = Dense.pre,
                 D.X.pre  = D.X.pre,
                 bw.est = bw.est,
                 width = width,
                 NCP.FIXED  = NCP.FIXED,
                 ZSDS.FIXED = ZSDS.FIXED,
                 W.FIXED = W.FIXED)

  }  else { 

      res.run <- run_EM(INT = INT, 
        ncp = ncp,
        zsds = zsds,
        Int.Beg = Int.Beg,
        Int.End = Int.End,
        NCP.FIXED = NCP.FIXED,
        ZSDS.FIXED = ZSDS.FIXED,
        W.FIXED = W.FIXED
      )

  }

    ncp = res.run$ncp
    zsds = res.run$zsds
    
    ## Compute Power 
    cp.input <- list(
      w.inp    = res.run$w.inp,
      ncp      = res.run$ncp,
      zsds     = res.run$zsds
    )

    cp.res <- Compute.Power.Z.General(cp.input, Int.Beg = Int.Beg, Int.End = Int.End)

    use.ncp = res.run$ncp
    use.zsds = res.run$zsds

    if(!is.null(yi)) { mu_tau = fit_fixed_weight_selection (
      yi, sei, crit, odr = ODR, edr = cp.res$EDR)
    } else { mu_tau = c(NA,NA,NA) }

    #mu_tau

    res.pe <- list(
      ODR         = ODR,
      EDR         = cp.res$EDR,
      ERR         = cp.res$ERR,
      ncp         = res.run$ncp,
      zsds        = res.run$zsds,
      w.inp       = res.run$w.inp,
      w.all       = cp.res$w.all,
      local.power = cp.res$local.power,
      fit         = res.run$fit,
      gamma_shape = NA,
      gamma_rate  = NA,
      gamma_edr   = c(NA,NA),
      gamma_er    = c(NA,NA),
      sel_ns_w    = mu_tau[1],
      mu          = mu_tau[2],
      tau         = mu_tau[3],
      ODR_EDR_D   = ODR - cp.res$EDR
     )

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

  ### res.pe

  ### boostrap

  #Est.Method = "OF"; boot.iter = 100
  if (boot.iter > 0) {  

       boot_res <- run_bootstrap(val.input = val.input, crit=crit,
          ncp        = use.ncp,
          zsds       = use.zsds,
          Int.Beg    = Int.Beg,
          Int.End    = Int.End,	
          boot.iter  = boot.iter,
          Dense.pre  = Dense.pre,
          D.X.pre    = D.X.pre,
          NCP.FIXED  = NCP.FIXED,
          ZSDS.FIXED = ZSDS.FIXED,
          W.FIXED    = W.FIXED,
          cluster.id = cluster.id,
          Est.Method = Est.Method)	

          valid_boot <- vapply(boot_res, function(x) {
             if (is.null(x$w.inp) || is.null(x$ncp) || is.null(x$zsds)) return(FALSE)
             w <- as.numeric(x$w.inp)
             n <- as.numeric(x$ncp)
             s <- as.numeric(x$zsds)

             length(w) > 0 &&
             length(n) > 0 &&
             length(s) > 0 &&
             any(is.finite(w)) &&
             any(is.finite(n)) &&
             any(is.finite(s)) &&
             sum(w, na.rm = TRUE) > 0
             }, logical(1))

             #mean(valid_boot)
             if(mean(valid_boot) < .90) stop

          boot_res <- boot_res[valid_boot]

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
                local.power     = rbind(res.pe$local.power,matrix(NA,2,length(res.pe$local.power))), 
                fit             = res.pe$fit,
                gamma_shape     = res.pe$gamma_shape,
                gamma_rate      = res.pe$gamma_rate,      
                gamma_edr       = gamma_edr,
                gamma_err       = gamma_err,
                ODR_EDR_D       = c(res.pe$ODR_EDR_D,NA,NA)
            )


       } else {

          ## ---- Compute power for each bootstrap result ----
           boot_power <- lapply(boot_res, function(x)
             Compute.Power.Z.General(x, Int.Beg = Int.Beg, Int.End = Int.End)
           )

          ODR_boot   <- sapply(boot_res, function(x) x$ODR)

          EDR_boot   <- sapply(boot_power, function(x) x$EDR)

          EDR_boot_lb <- pmax(alpha, EDR_boot - CI.EDR.MIN.ADJ)
          EDR_boot_ub <- pmin(1,     EDR_boot + CI.EDR.MIN.ADJ)

          ODR_ci = c(
             quantile(ODR_boot, probs = .025, na.rm = TRUE),
             quantile(ODR_boot, probs = .975, na.rm = TRUE)
          )

          EDR_ci = c(
             quantile(EDR_boot_lb, probs = .025, na.rm = TRUE),
             quantile(EDR_boot_ub, probs = .975, na.rm = TRUE)
          )


          ODR_EDR_D_ci  <- c(
               quantile(ODR_boot - EDR_boot_lb, .025),
               quantile(ODR_boot - EDR_boot_ub, .975)
           )

          

          if(!is.null(yi)) {

            mu_tau_lb <- apply(
              cbind(ODR_boot, EDR_boot_lb),1,
              function(x) fit_fixed_weight_selection(yi, sei, x[1], x[2])
            )

           mu_tau_ub <- apply(
             cbind(ODR_boot, EDR_boot_ub),1,
             function(x) fit_fixed_weight_selection(yi, sei, x[1], x[2])
           )

            # If rows are: 1 = sel_ns_w, 2 = mu, 3 = tau
            sel_ns_w_ci <- c(
               quantile(mu_tau_lb[1, ], probs = .025, na.rm = TRUE),
               quantile(mu_tau_ub[1, ], probs = .975, na.rm = TRUE)
            )

            # Envelope approach for mu and tau
            mu_lower_draws <- pmin(mu_tau_lb[2, ], mu_tau_ub[2, ], na.rm = TRUE)
            mu_upper_draws <- pmax(mu_tau_lb[2, ], mu_tau_ub[2, ], na.rm = TRUE)

            tau_lower_draws <- pmin(mu_tau_lb[3, ], mu_tau_ub[3, ], na.rm = TRUE)
            tau_upper_draws <- pmax(mu_tau_lb[3, ], mu_tau_ub[3, ], na.rm = TRUE)

            mu_ci <- c(
             quantile(mu_lower_draws, probs = .025, na.rm = TRUE),
             quantile(mu_upper_draws, probs = .975, na.rm = TRUE)
            )

            tau_ci <- c(
             quantile(tau_lower_draws, probs = .025, na.rm = TRUE),
             quantile(tau_upper_draws, probs = .975, na.rm = TRUE)
            )

          } else {

            sel_ns_w_ci = c(NA,NA)
            mu_ci = c(NA,NA)
            tau_ci = c(NA,NA)
        
          } # EOF mu and tau

          lcbind <- function(list1, list2) {
            stopifnot(length(list1) == length(list2))
            lapply(seq_along(list1), function(i) c(list1[[i]], list2[[i]]))
          }

          boot_list <- lcbind(boot_res, boot_power)

          boot_mat <- do.call(rbind, lapply(boot_list, function(x) {
            unlist(x) }))

          ci_lo <- apply(boot_mat, 2, quantile, probs = 0.025, na.rm = TRUE)
          ci_hi <- apply(boot_mat, 2, quantile, probs = 0.975, na.rm = TRUE)

          ci <- cbind(ci_lo, ci_hi)

          ODR = c(res.pe$ODR,ODR_ci)
          EDR = c(res.pe$EDR,EDR_ci)
          ERR = c(res.pe$ERR,ci["ERR",])

          EDR[2] = max(alpha,EDR[2])
          EDR[3] = min(1,EDR[3])

          ERR[2] = ERR[2] - CI.ERR.MIN.ADJ
          ERR[3] = ERR[3] + CI.ERR.MIN.ADJ

          ERR[2] = max(alpha/2,ERR[2])
          ERR[3] = min(1,ERR[3])


          sel_ns_w  = c(res.pe$sel_ns_w,sel_ns_w_ci)
          mu  = c(res.pe$mu[1],mu_ci)
          tau = c(res.pe$tau[1],tau_ci)


          ODR_EDR_D = c(res.pe$ODR_EDR_D,ODR_EDR_D_ci)

          if(length(res.pe$ncp) == 1) {
             ncp = t(c(res.pe$ncp,ci["ncp",]))
          } else {
             ncp  = t(cbind(res.pe$ncp, ci[grepl("ncp",rownames(ci)),]) )
          }


          if(length(res.pe$zsds) == 1) {
             zsds = t(c(res.pe$zsds,ci["zsds",]))
          } else { 
             zsds  = t(cbind(res.pe$zsds, ci[grepl("zsds",rownames(ci)),]))
          }


          if(length(res.pe$w.all) == 1) {
             w.all = t(c(res.pe$w.all,ci["w.all",]))
          } else { 
             w.all  = t(cbind(res.pe$w.all, ci[grepl("w.all",rownames(ci)),]))
          }


          if(length(res.pe$w.inp) == 1) {
             w.inp = t(c(res.pe$w.inp,ci["w.inp",]))
          } else { 
             w.inp  = t(cbind(res.pe$w.inp, ci[grepl("w.inp",rownames(ci)),]))
          }

          local.power = t(cbind(res.pe$local.power, ci[grepl("local",rownames(ci)),]))

          zcurve.res = list(
            ODR                = ODR,
            EDR                = EDR,
            ERR                = ERR,
            ncp                = ncp,
            zsds               = zsds,
            w.inp              = w.inp,
            w.all              = w.all,
            local.power        = local.power,
            gamma_shape        = res.pe$gamma_shape,
            gamma_rate         = res.pe$gamma_rate,
            gamma_edr          = res.pe$gamma_edr,
            gamma_err          = res.pe$gamma_err,
            ODR_EDR_D          = ODR_EDR_D,
            sel_ns_w           = sel_ns_w,
            mu                 = mu,
            tau                = tau
          )

          if(length(zcurve.res$zsds) < 3) break

       } # End of OF or EM

  } else { # End of if Bootstrap

        ODR = t(c(res.pe$ODR,rep(NA,2)))
        EDR = t(c(res.pe$EDR,rep(NA,2)))
        ERR = t(c(res.pe$ERR,rep(NA,2)))
        ncp = t(cbind(res.pe$ncp, matrix(NA,length(res.pe$w.all),2)))
        zsds = t(cbind(res.pe$zsds, matrix(NA,length(res.pe$w.inp),2)))
        w.all = t(cbind(res.pe$w.all, matrix(NA,length(res.pe$w.all),2)))
        w.inp = t(cbind(res.pe$w.all, matrix(NA,length(res.pe$w.inp),2)))
        local.power = t(cbind(res.pe$local.power,matrix(NA,length(res.pe$local.power),2)))

        zcurve.res = list(
            ODR                = ODR,
            EDR                = EDR,
            ERR                = ERR,
            ncp                = ncp,
            zsds               = zsds,
            w.inp              = w.inp,
            w.all              = w.all,
            local.power        = local.power,
            gamma_shape        = res.pe$gamma_shape,
            gamma_rate         = res.pe$gamma_rate,
            gamma_edr          = res.pe$gamma_edr,
            gamma_err          = res.pe$gamma_err,
            ODR_EDR_D          = rep(NA,3)
        )

    } # EOF no bootstrap

})

print(zcurve.time) ### End of Timer

#print(zcurve.res)

return(zcurve.res)

}

###################

bias.test = function(val.input, caliper.width = .8) {

  #clu.id = NULL

  just = caliper.width

  ncp.bias = 0:6
  zsds.bias = rep(1,7)

  bias.res = run_zcurve(Est.Method="OF",val.input=val.input,crit=crit,Int.Beg = 0,Int.End = Int.End,
   ncp=ncp.bias,zsds=zsds.bias,cluster.id=NULL,boot.iter = 0,yi=NULL,sei=NULL)

  bar.width = .001
  D.X = seq(0,Int.End,bar.width)
  dense = build.dense(D.X,ncp.bias,zsds.bias,CURVE.TYPE=CURVE.TYPE)
  E.D.Y <- as.numeric(crossprod(bias.res$w.all[1,], dense))

  prob1 = sum(E.D.Y[D.X > crit & D.X < crit + just]*bar.width);prob1
  prob2 = sum(E.D.Y[D.X >= crit + just & D.X < Int.End]*bar.width);prob2
  prob3 = sum(E.D.Y[D.X < crit & D.X > crit - just]*bar.width);prob3
  prob4 = sum(E.D.Y[D.X <= crit - just]*bar.width);prob4
  sum(prob1,prob2,prob3,prob4)
  prob = prob1 / (prob1 + prob2)
  
  bias.int = val.input[val.input <= Int.End]
  k = length(bias.int)
  sig.k <- sum(bias.int > crit)
  just.sig.k = sum(bias.int > crit & bias.int < crit + just)
  just.sig.k
  ojs = just.sig.k/sig.k

  ### binomial
  p.bias.binomial = NA
  p.bias.binomial = pbinom(sig.k - just.sig.k, sig.k, prob = 1-prob)

  EJST = c(round(ojs*100,2),round(prob*100,2),
    p.bias.binomial)
  names(EJST) = c("EJST OJS","EJST EJS","EJS p-value")

  # (a) deficit of just-non-significant
  jns.k <- sum(bias.int > crit - just & bias.int < crit)
  #jns.k;k;prob3;pbinom(ns.k, k, prob = prob3)
  p.missing.just.ns <- pbinom(jns.k, k, prob = prob3)
  MJNST = c(round(jns.k/k*100,2),round(prob3*100,2),p.missing.just.ns)
  names(MJNST) = c("MJNS OBS","MJNS EXP","MJNS p-value")

  #combined EJST MJNST: ratio of just-sig to just-ns vs model prediction
  obs.prop <- just.sig.k / (just.sig.k + jns.k)
  exp.prop <- prob1 / (prob1 + prob3)  # from the fitted density
  p.ZCT = NA
  if(jns.k > 0) p.ZCT <- binom.test(just.sig.k, just.sig.k + jns.k, p = exp.prop,
                         alternative = "greater")$p.value
  ZCT = c(round(obs.prop*100,2),round(exp.prop*100,2),p.ZCT)
  names(ZCT) = c("ZCT OBS","ZCT EXP","ZCT p-value") 


ZBIAS = list(
   EJST = EJST,
   MJNST = MJNST,
   ZCT = ZCT
)

ZBIAS

return(ZBIAS)

} #EOF bias test


############################################

Write.Local.Power = function(local.power) {

  if (int.loc == 0) return(invisible(NULL))

  # Midpoints of each bin — these are the correct label positions
  midpoints = seq(x.lim.min, x.lim.max - int.loc, by = int.loc) + int.loc / 2

  # Format labels
  lab = paste0(round(local.power * 100), "%")

  # Add bottom margin for the extra label row
  old_mar = par("mar")
  par(mar = old_mar + c(1.5, 0, 0, 0))

  mtext(lab, side = 1, line = 1.5, at = midpoints, cex = 1.0, las = 1)

  par(mar = old_mar)

} ### EOF Write.Local.Power



###################################################
#### Begin Draw Histogram
###################################################

Draw.Histogram = function(w,cola="blue3",
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

ncp = results$ncp[1,]
zsds = results$zsds[1,]
w = results$w.all[1,]

width = 1

bar.width = .01
D.X = seq(x.lim.min,x.lim.max,bar.width) 

D.Y = build.dense(D.X, ncp, zsds, df, CURVE.TYPE) 
D.Y = D.Y / (rowSums(D.Y) * bar.width)
D.Y = colSums(D.Y * w)

#plot(D.X,D.Y,type="l",ylim=c(0,1));par(new=TRUE);plot(D.X,D.Y.g,type="l",lty=2,col="chartreuse3",ylim=c(0,1))

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

skip_zcurve = FALSE

if (length(p) > 0) {
  val.input = qnorm(1-p/2)
}


if (CURVE.TYPE == "z") { 
  crit = qnorm(alpha/2,lower.tail=FALSE)
 } else {
  crit = qt(alpha/2,df,lower.tail=FALSE)
}


### remove NA
missing = is.na(val.input)
val.input = val.input[!missing]
if (!is.null(cluster.id)) cluster.id = cluster.id[!missing]

#if (two.sided) val.input = abs(val.input)

### create set with z-scores in the interval used for model fitting

if (length(val.input[val.input > Int.Beg & val.input < Int.End]) > 10) {
  print("!START!") 
} else {
  print("Insufficient Data")
  skip_zcurve = TRUE
}    

n.z = length(val.input);n.z
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


if(!skip_zcurve) {

g_ncp = ncp
g_zsds = zsds

#cluster.id = NULL
### Heterogeneity Estimation
heterogeneity.test = run_zcurve(Est.Method="OF",val.input=val.input,crit=crit,
     Int.Beg = Int.Beg,Int.End=Int.End,ncp = 2.5, zsds = 2,
     NCP.FIXED=FALSE,ZSDS.FIXED=FALSE,W.FIXED=TRUE,
     cluster.id = cluster.id, boot.iter = 500,Run.Gamma=FALSE,
     yi=NULL,sei=NULL)

heterogeneity = sqrt(heterogeneity.test$zsds^2 - 1)
#heterogeneity

if(length(heterogeneity) < 3) break

if (mean(heterogeneity) < min.heterogeneity) { 
  quick.test = run_zcurve(Est.Method="OF",val.input=val.input,crit=crit,Int.Beg=Int.Beg,Int.End=Int.End,
      cluster.id=cluster.id,boot.iter = 0,ncp=2.5,zsds=1,
      NCP.FIXED=FALSE,ZSDS.FIXED=TRUE,Run.Gamma=FALSE, yi=NULL, sei=NULL)
  quick.test$ncp #quick.test$ncp[1] = 1.5
  if (all(abs(ncp - quick.test$ncp[1]) > .1)) {
    ncp_val = round(quick.test$ncp[1], 2)
    ncp = sort(unique(c(0, seq(ncp_val, 0, -1), seq(ncp_val, 6, 1))))
    components = length(ncp)
    zsds = rep(1, components)
  }
}


if (heterogeneity[1,3] < .2) { # high confidence that data are homo
  ncp = quick.test$ncp[1]
  zsds = 3
  NCP.FIXED=FALSE
  ZSDS.FIXED=FALSE
}

ncp
zsds
NCP.FIXED
ZSDS.FIXED

sel.width = .5
k.sel = length(val.input[val.input > Int.Beg & val.input < Int.Beg + sel.width])
#k.sel = 1000
edr.factor <- max(.01,1/sqrt(k.sel)) 
edr.factor
edr.adj.ci = 1 * edr.factor 
edr.adj.ci = max(CI.EDR.MIN.ADJ,edr.adj.ci)
edr.adj.ci

BIAS = NULL
if (TEST4BIAS) BIAS = run_bias_tests()

#}
### BBB 

#boot.iter = 500
print("MAIN")
res = run_zcurve(Est.Method = Est.Method,val.input = val.input,crit=crit,Int.Beg=Int.Beg,Int.End=Int.End,
                 cluster.id=cluster.id,boot.iter = boot.iter,ncp = ncp, zsds = zsds, 
                 NCP.FIXED=NCP.FIXED,ZSDS.FIXED=ZSDS.FIXED,edr_adj_ci=edr.adj.ci,
                 Run.Gamma=Run.Gamma,yi=yi,sei=sei)

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


FDR = as.numeric((1/res$EDR - 1)*(alpha/(1-alpha)))
FDR = FDR[c(1,3,2)]
names(FDR) = c("FDR_pe","FDR_lb","FDR_ub")

POW_h1_sig <- as.numeric((res$ERR - FDR * alpha) / (1 - FDR))
names(POW_h1_sig) = c("POW_H1_SIG_pe","POW_H1_SIG_lb","POW_H1_SIG_ub")

ODR_EDR_D = as.numeric(res$ODR_EDR_D)
names(ODR_EDR_D) = c("ODR_EDR_D_pe","ODR_EDR_D_lb","ODR_EDR_D_ub")


if(boot.iter > 0 & is.na(res$ODR_EDR_D[1])) stop


} # skip zcurve

#skip_zcurve

#res

###############################################

###############################################
###############################################

#skip_zcurve

if(skip_zcurve) {

  res = list(
   ODR = rep(NA,3),
   EDR = rep(NA,3),
   ERR = rep(NA,3),
   ncp = matrix(NA,3,7),
   zsds = matrix(NA,3,7),
   w.inp = matrix(NA,3,7),
   w.all = matrix(NA,3,7),
   local.power = matrix(NA,3,12),
   fit = rep(NA,3),
   gamma_shape = rep(NA,3),
   gamma_rate = rep(NA,3),
   gamma_edr = rep(NA,3),
   gamma_err = rep(NA,3)
   )

BIAS = list(
  Credibility = rep(NA,3),
  Rindex      = rep(NA,3),
  Caliper = rep(NA,3))

BIAS$ZBIAS = list(
  EJST    = rep(NA,3),
  MJNST   = rep(NA,3),    
  ZCT     = rep(NA,3),
  RDD     = rep(NA,3),
  Hybrid  = rep(NA,3)
)

BIAS

POW_h1_sig = rep(NA,3)
FDR = rep(NA,3)
ODR_EDR_D = rep(NA,3)
heterogeneity = rep(NA,3)
adjustments = rep(NA,5)

}

names(ODR_EDR_D) = c("ODR_EDR_D_pe","ODR_EDR_D_lb","ODR_EDR_D_ub")

res$ODR = as.numeric(res$ODR)
res$EDR = as.numeric(res$EDR)
res$ERR = as.numeric(res$ERR)

names(res$ODR) = c("ODR_pe","ODR_lb","ODR_ub")
names(res$EDR) = c("EDR_pe","EDR_lb","EDR_ub")
names(res$ERR) = c("ERR_pe","ERR_lb","ERR_ub")

heterogeneity = as.numeric(heterogeneity)
names(heterogeneity) = c("ncp_sd_pe","ncp_sd_lb","ncp_sd_ub")

names(adjustments) = c("edr.factor","edr.adj.ci","slope.adj","err.factor","err.adj")


if (length(res$local.power) > 3) {
   colnames(res$local.power) = paste0("local.power",1:ncol(res$local.power))
   rownames(res$local.power) = c("pe","lb","ub")
}

if(length(ncp) > 1) {
  colnames(res$ncp) = paste0("ncp",1:length(ncp))
  rownames(res$ncp) = c("pe","lb","ub")

  colnames(res$zsds) = paste0("zsds",1:length(ncp))
  rownames(res$zsds) = c("pe","lb","ub")

  colnames(res$w.inp) = paste0("w.inp",1:length(ncp))
  rownames(res$w.inp) = c("pe","lb","ub")

  colnames(res$w.all) = paste0("w.all",1:length(ncp))
  rownames(res$w.all) = c("pe","lb","ub")
}


n.z.clusters = NULL
if (!is.null(cluster.id)) n.z.clusters <- length(unique(cluster.id))

results = list(
		clusters       = n.z.clusters,
		slope          = slope,
		ODR            = res$ODR,
		EDR            = res$EDR,
		ERR            = res$ERR,
		FDR            = FDR,
		POW_h1_sig     = POW_h1_sig,
		ncp            = res$ncp,
		zsds           = res$zsds,
		heterogeneity  = heterogeneity,
		w.inp          = res$w.inp,
		w.all          = res$w.all,	
		local.power    = res$local.power,
		BIAS           = BIAS,
		ODR_EDR_D      = ODR_EDR_D,
         adjustments    = adjustments,
		fit            = res$fit,         
		gamma_shape    = res$gamma_shape,
		gamma_rate     = res$gamma_rate,
		gamma_edr      = res$gamma_edr,
		gamma_err      = res$gamma_err,
         sel_ns_w       = res$sel_ns_w,
         mu             = res$mu,
         tau            = res$tau
      )


#print(results)


##########################################
### This Code is Used to Create Graphic
##########################################


if (Show.Histogram & !skip_zcurve) {

	if (boot.iter == 0) {
		#print("Draw Histogram NO CI")
		Draw.Histogram(results$w.all[1,],cola=col.hist,Write.CI = FALSE)
	} else {
		#print("Draw Histogram WITH CI")
		Draw.Histogram(results$w.all[1,],cola=col.hist,Write.CI = TRUE)
	}

	if (Show.Curve.All) { 


		Draw.Curve.All(results,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max,Show.Gamma=Show.Gamma)
		Draw.Curve.All(results,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End,Show.Gamma=Show.Gamma)

	} # EOF Show.Curve.All


	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)


	if (length(results$local.power) > 0 && !is.na(results$local.power[1])) Write.Local.Power(results$local.power[1,])	


} # End of Show.Histogram	for bootstrap


### reset to global specifications
ncp  = g_ncp
zsds = g_zsds

return(results)

} #End of Zing


######################################################################
######################################################################
######################################################################


zcurve.es.estimation = function(N,pow) {

ZCURVE.EM.HOM = c()

pow.i = 1
for (pow.i in 1:length(pow)) {

	
	if (pow[pow.i] <= .05) { ZCURVE.EM.HOM = c(ZCURVE.EM.HOM,0) } else {

		n = round(N/2)

		avg.n = mean(N)/2;avg.n
		d.1 = pwr.t.test(
		n = avg.n, d = NULL, sig.level = 0.05, power = pow[pow.i],
		type = "two.sample", alternative = "two.sided")$d 
		d.1

		d.2 = unlist(lapply(n,function(x) pwr.t.test(
		n = x, d = NULL, sig.level = 0.05, power = pow[pow.i],
		type = "two.sample", alternative = "two.sided")$d ))
		d.2 = mean(d.2)
		d.2

		### do not use d.1
		ZCURVE.EM.HOM = c(ZCURVE.EM.HOM,d.2)
	} # end else

} # End pow.i

#print(ZCURVE.EM.HOM)

return(ZCURVE.EM.HOM)

} # End of Function



