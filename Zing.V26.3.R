
#rm(list = ls())

########################################################################
### SETTING PARAMETERS FOR Z-CURVE MODEL
########################################################################

version <- "Version.3.4 (2026.03.18)"   # Version label to appear on plots

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
ymin <- 0                    # OUTDATED Y-axis lower bound (for label space)

Show.Histogram <- TRUE       # Toggle histogram in plot
Show.Text <- TRUE            # Toggle model results in plot
Show.Curve.All <- TRUE       # Show predicted z-curve
Show.Curve.Sig <- FALSE      # Option: show z-curve only for significant values
Show.Significance <- TRUE    # Show z = critical value line
Show.KD <- FALSE             # Toggle kernel density overlay (density method only)

sig.levels <- c()            # Optional: mark additional p-value thresholds on plot

int.loc <- 0.5               # Plot local power intervals below x-axis (set 0 to disable)
hist.bar.width <- 0.2        # Width of histogram bars

### CONSOLE OUTPUT

Show.Iterations <- FALSE      # Do not use with parallel Show iterations for slow procedures (e.g., EXT, TEST4HETEROGENEITY)

### MODEL PARAMETERS

alpha <- 0.05                        # Significance level
crit <- qnorm(1 - alpha / 2)       # Corresponding two-sided critical z

two.sided <- TRUE                   # Assume two-sided z-values (use abs(z)); not yet compatible with signed z-values

# Color scheme
col.curve <- "violetred3"
col.hist <- "blue3"
col.kd <- "black"

Est.Method <- "OF"                  # Estimation method: "OF", "EM", or "EXT"
                                    # Clustered Data: "CLU-W" (weighted),"CLU-B" (bootstrap)   
Int.Beg <- 1.96                     # Default: critical value for alpha = .05
Int.End <- 6                        # End of modeling interval (z > 6 = power = 1)

ncp <- 0:6                          # Component locations (z-values at which densities are centered)
components <- length(ncp)           # Number of components
zsd <- 1                            # SD of standard normal z-distribution
zsds = rep(zsd,components)          # one SD for each component

just <- 0.8                         # Cutoff for "just significant" z-values (used in optional bias test)

ZSDS.FIXED <- TRUE                  # Fix SD values for EXT method 
NCP.FIXED <- TRUE                   # Fix non-central parameter(NCP) means values for EXT method
W.FIXED   <- FALSE                  # Fix weights for EXT method

fixed.false.positives <- 0          # If > 0, constrains proportion of false positives (e.g., weight for z = 0 component)

### DENSITY-BASED SETTINGS (Only used with Est.Method = "OF")

n.bars <- 512                       # Number of bars in histogram

Augment <- TRUE                     # Apply correction for bias at lower bound

bw.est  <- 0.05                     # Bandwidth for kernel density (lower = less smoothing, higher = more smoothing)
bw.draw <-  .20   		        	 # Bandwith of kernel density in plot

### INPUT RESTRICTIONS

MAX.INP.Z <- Inf                    # Optionally restrict very large z-values (set Inf to disable)

### CONFIDENCE INTERVALS / BOOTSTRAPS

boot.iter <- 0                      # Number of bootstrap iterations (suggest 500+ for final models)
ERR.CI.adjust <- 0.03               # Conservative widening of confidence intervals for ERR
EDR.CI.adjust <- 0.05               # Conservative widening for EDR

CI.ALPHA <- 0.05                    # CI level (default = 95%)

### CI levels for Heterogeneity Test
fit.ci <- c(.01, .025, .05, .10, .17, .20, .50, .80, .83, .90, .95, .975, .99)  # CI levels for model fit test

TEST4BIAS <- FALSE                  # Enable optional bias test
TEST4HETEROGENEITY <- 0             # Optional heterogeneity test (slow) — set number of bootstrap iterations


### DISPLAY FINAL STATUS

print(version)
print("Parameter OK")

##################################################################
### END OF SETTING DEFAULT PARAMETERs
##################################################################


Zing = function(val.input,df=c(),cluster.id=c(),p=c(), lp=c()) {

##################################################################
##################################################################
### AAA Start of Zing Code
##################################################################
##################################################################


#####################################################################
### Z-curve 3.0 - NEW EM Estimation for Extended (free w, ncp, zsds)
#####################################################################


round_list <- function(x, digits = 3) {
  if (is.numeric(x)) {
    round(x, digits)
  } else if (is.list(x)) {
    lapply(x, round_list, digits = digits)
  } else {
    x
  }
}

##################################################
### Begin of New Extended Zcurve
##################################################

EM_EXT_FOLDED <- function(
  INT,
  ncp,
  zsds = rep(1, length(ncp)),
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
    w.inp.est    = w.inp,
    ncp          = ncp,
    zsds         = zsds,
    loglik       = loglik,
    loglik_trace = loglik_trace[1:iter],
    iter         = iter
  )
}



##################################################
### End of Core EM
##################################################


EM_EXT_FOLDED_MULTISTART <- function(
  INT,
  ncp,
  zsds       = rep(1, length(ncp)),
  Int.Beg    = 1.96,
  Int.End    = max(INT),
  NCP.FIXED  = TRUE,
  ZSDS.FIXED = TRUE,
  W.FIXED    = FALSE,
  max_iter   = 200,
  tol        = 1e-6,
  n_starts   = 10,
  w_start    = NULL   # optional external starting weights
) {

  components <- length(ncp)
  ncp_init   <- ncp
  zsds_init  <- zsds

  best_fit <- NULL
  best_ll  <- -Inf

  for (s in 1:n_starts) {

    # --- starting weights ---
    if (s == 1 && !is.null(w_start)) {
      # use supplied weights for first start
      w_s <- w_start
      w_s <- pmax(w_s, 1e-3)       # small floor to avoid zeros
      w_s <- w_s / sum(w_s)
    } else if (s == 1) {
      w_s <- rep(1/components, components)   # uniform
    } else {
      raw <- rgamma(components, shape = 1)   # random Dirichlet
      w_s <- raw / sum(raw)
    }


fit <- tryCatch(
      EM_EXT_FOLDED(
        INT        = INT,
        ncp        = ncp_init,
        zsds       = zsds_init,
        w.inp      = w_s,
        Int.Beg    = Int.Beg,
        Int.End    = Int.End,
        NCP.FIXED  = NCP.FIXED,
        ZSDS.FIXED = ZSDS.FIXED,
        W.FIXED    = W.FIXED,
        max_iter   = max_iter,
        tol        = tol
      ),

      error = function(e) { message("Start ", s, ": ", e$message); NULL }
#      error = function(e) NULL
    )

    if (!is.null(fit) && is.finite(fit$loglik) && fit$loglik > best_ll) {
      best_ll  <- fit$loglik
      best_fit <- fit
      best_fit$start <- s   # track which start won
    }
  }

  return(best_fit)
}


##################################################
### Bootstrap
##################################################

run_bootstrap_list <- function(INT,
                               ncp.start,
                               zsds.start,
                               Int.Beg,
                               Int.End,
                               boot.iter  = 1000,
                               NCP.FIXED  = TRUE,
                               ZSDS.FIXED = TRUE,
                               W.FIXED    = FALSE,
                               n_starts   = 10,
                               cores      = round(parallel::detectCores() * .8)) {

  ncores <- max(1, cores)

  cl <- parallel::makeCluster(cores)

  parallel::clusterExport(
    cl,
    varlist = c("EM_EXT_FOLDED",
                "EM_EXT_FOLDED_MULTISTART",
                "INT",
                "ncp.start",
                "zsds.start",
                "Int.Beg",
                "Int.End",
                "NCP.FIXED",
                "ZSDS.FIXED",
                "W.FIXED",
                "n_starts"),
    envir = environment()
  )

  parallel::clusterEvalQ(cl, library(stats))

  boot_res <- parallel::parLapply(cl, 1:boot.iter, function(i) {

    boot_sample <- sample(INT, length(INT), replace = TRUE)

    fit <- tryCatch(
      EM_EXT_FOLDED_MULTISTART(
        INT        = boot_sample,
        ncp        = ncp.start,
        zsds       = zsds.start,
        Int.Beg    = Int.Beg,
        Int.End    = Int.End,
        NCP.FIXED  = NCP.FIXED,
        ZSDS.FIXED = ZSDS.FIXED,
        W.FIXED    = W.FIXED,
        n_starts   = n_starts
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) return(NULL)

    list(
      w.inp  = fit$w.inp.est,
      ncp    = fit$ncp,
      zsds   = fit$zsds,
      fit = fit$loglik
    )
  })

  parallel::stopCluster(cl)

  boot_res <- boot_res[!sapply(boot_res, is.null)]

  if (length(boot_res) == 0)
    stop("All bootstrap runs failed.")

  return(boot_res)
}


##################################################
### End of New Extended Zcurve
##################################################


run.new.EM = function(val.input, w_start, boot.iter=0, NCP.FIXED=TRUE,ZSDS.FIXED = TRUE) {

INT = val.input[val.input > Int.Beg & val.input < Int.End]

res.run <- EM_EXT_FOLDED_MULTISTART(INT = INT,
              ncp = ncp,
              zsds = zsds,
			 w_start = w_start,
              Int.Beg = Int.Beg,
              Int.End = Int.End,
              NCP.FIXED = NCP.FIXED,
              ZSDS.FIXED = ZSDS.FIXED
		)


res.run

cp.input = list(
	w.inp = res.run$w.inp.est,
	ncp.est = res.run$ncp,
	zsds.est = res.run$zsds    
)

cp.input

cp.res = Compute.Power.Z.General(cp.input,Int.Beg=Int.Beg)

res.pe = list(
   EDR = cp.res$EDR,
   ERR = cp.res$ERR,
   ncp = res.run$ncp,
   zsds = res.run$zsds,
   w.inp = res.run$w.inp,
   w.sig = cp.res$w.sig,
   w.all = cp.res$w.all,
   local.power = cp.res$local.power,
   fit = res.run$loglik
)

######################

if (boot.iter > 0) {

boot_res <- run_bootstrap_list(INT = INT,
                          ncp.start = ncp,
                          zsds.start = zsds,
                          Int.Beg = Int.Beg,
                          Int.End = Int.End,		
                          boot.iter = boot.iter,
                          NCP.FIXED = NCP.FIXED,
                          ZSDS.FIXED = ZSDS.FIXED,
                          W.FIXED = FALSE)


boot_power <- lapply(boot_res, function(x)
  Compute.Power.Z.General(x, Int.Beg = Int.Beg)
)

lcbind <- function(list1, list2) {
  stopifnot(length(list1) == length(list2))
  lapply(seq_along(list1), function(i) c(list1[[i]], list2[[i]]))
}


boot_list = lcbind(boot_res,boot_power)


  lquantile <- function(boot_list, probs = c(0.025, 0.975)) {
    fields <- names(boot_list[[1]])
    
    result <- lapply(fields, function(f) {
      vals <- lapply(boot_list, function(x) x[[f]])
      if (all(sapply(vals, is.null))) return(NULL)
      mat <- do.call(rbind, vals)
      apply(mat, 2, quantile, probs = probs, na.rm = TRUE)
    })
  
    names(result) <- fields
    result[!sapply(result, is.null)]
  }


  get.ci <- function(boot_list, probs = c(0.025, 0.975)) {
    ci <- lquantile(boot_list, probs = probs)
  }

  res.ci = get.ci(boot_list)

} else {  ### if no bootstrap ->

  res.ci = list(
    w.inp = matrix(NA,2,components), 
    w.all = matrix(NA,2,components),
    w.sig = matrix(NA,2,components),
    ncp = matrix(NA,2,components),
    zsds = matrix(NA,2,components),
    fit = matrix(NA,2,1),
    EDR = matrix(NA,2,1),
    ERR = matrix(NA,2,1),
    local.power = matrix(NA,2,length(res.pe$local.power))
  )

}


combine.pe.ci = function(res.pe,res.ci) {

  fields <- names(res.pe)
  result <- lapply(fields, function(f) {
    pe <- res.pe[[f]]
    ci_f <- res.ci[[f]]
    
    # If scalar field, quantile returns a vector not a matrix
    if (is.null(dim(ci_f))) {
      rbind(estimate = pe, ci_lb = ci_f[1], ci_ub = ci_f[2])
    } else {
      rbind(estimate = pe, ci_lb = ci_f[1, ], ci_ub = ci_f[2, ])
    }
  })
  
  names(result) <- fields
  result
}

zcurve.res <- combine.pe.ci(res.pe,res.ci)

zcurve.res$EDR[2] = zcurve.res$EDR[2] - EDR.CI.adjust
zcurve.res$EDR[3] = zcurve.res$EDR[3] + EDR.CI.adjust

zcurve.res$ERR[2] = zcurve.res$ERR[2] - ERR.CI.adjust
zcurve.res$ERR[3] = zcurve.res$ERR[3] + ERR.CI.adjust
 
zcurve.res

return(zcurve.res)

} # End of New Zcurve


#zcurve.res

##########################################

### TESTING NEW ZCURVE

if (1 == 2) {

source(zcurve3)

boot.iter = 100

zcurve.res = run.new.EM(val.input,boot.iter,NCP.FIXED,ZSDS.FIXED)

zcurve.res.3d = round_list(zcurve.res)
zcurve.res.3d

}


############################################################
### Functions used in Zing
############################################################


############################################################
### SLOPE DIAGNOSTICS
############################################################

get_slope <- function(INT, Int.Beg = 1.96, width = 1, min_in_window = 10, bw = 0.1) {

from <- Int.Beg + 2 * bw
to <- Int.Beg + 2 * bw + width

k.slope = length(INT[INT >= from & INT <= to])
  
if (k.slope >= min_in_window) {

  d <- Get.Densities(INT, bw = bw, from = from, to = to,
       width = width, Augment = FALSE)
  fit <- lm(d$ZY ~ d$ZX)
  slope <- unname(coef(fit)[2])

  out = c(k.slope = k.slope, slope = slope)

} else { 
  out = c(k.slope = k.slope, slope = NA)
}

out

} # EOF get_slope

############################################################
### END OF SLOPE DIAGNOSTICS


build.Dens <- function(D.X, ncp, zsds, df, CURVE.TYPE) {

  n.bars <- length(D.X)
  components <- length(ncp)
  Dens <- matrix(0, nrow = components, ncol = n.bars)
  
  if (CURVE.TYPE == "z") {
    for (i in 1:n.bars) {
      for (j in 1:components) {
        Dens[j, i] <- dnorm(D.X[i], ncp[j], zsds[j]) +
                       dnorm(-D.X[i], ncp[j], zsds[j])
      }
    }
  } else {
    for (i in 1:n.bars) {
      for (j in 1:components) {
        Dens[j, i] <- dt(D.X[i], df, ncp[j]) +
                       dt(-D.X[i], df, ncp[j])
      }
    }
  }

  Dens
  
  return(Dens)
}


##############################################
### Get Densities
##############################################

Get.Densities = function(INT, bw = 0.20, from = 1.96, to = 6, width = 1, Augment = TRUE) {

  #from = 1.96; to = 6;bw = .2

  n = 401
  grid <- seq(from, to, length.out = n)
  dens <- numeric(n)

  width = 0
  
  if (Augment) {

    bnd.width <- 2 * bw

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
    bkde.out <- bkde(INT, bandwidth = .05, range.x = c(from, to), gridsize = n)
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

} # EOF Get.Densities 

#######################################################
### End of Get Densities
#######################################################


#######################################################
### Old Fashioned Density Method
#######################################################

old.fashioned <- function(val.input, cola = "springgreen4",bw.est = .2,
                          Dens.precomputed = NULL, D.X.precomputed = NULL,
                          width = 1, 
                          NCP.FIXED = TRUE, ZSDS.FIXED = TRUE) {


  #Dens.precomputed = NULL; D.X.precomputed = NULL

  components <- length(ncp)

  ## ---- Observed density ----
  INT   <- val.input[val.input >= Int.Beg & val.input <= Int.End]
  densy <- Get.Densities(INT, bw = bw.est, from = Int.Beg,
                         to = Int.End, Augment = Augment)

  D.X   <- densy[, 1]
  O.D.Y <- densy[, 2]
  n.bars    <- length(D.X)
  bar.width <- D.X[2] - D.X[1]

  #plot(D.X,O.D.Y)

  use.precomputed <- !is.null(Dens.precomputed) &&
                     !is.null(D.X.precomputed) &&
                     identical(D.X, D.X.precomputed)

  ## ---- theta = [weights, ncp, zsds] — always ----
  n.wt <- components
  if (fixed.false.positives > 0 & 0 %in% ncp) n.wt <- components - 1


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
    if (fixed.false.positives > 0)
      wt <- c(fixed.false.positives, wt * (1 - fixed.false.positives))

    if (use.precomputed && NCP.FIXED && ZSDS.FIXED) {
      Dens.now <- Dens.precomputed
    } else {
      Dens.now <- build.Dens(D.X, theta[idx.ncp], theta[idx.zsd], df, CURVE.TYPE)
      Dens.now <- Dens.now / (rowSums(Dens.now) * bar.width)
    }

    E.D.Y <- colSums(Dens.now * wt)
    #plot(D.X,E.D.Y)
    rmse = sqrt(mean((E.D.Y - O.D.Y)^2))

    scale = 1
    ## ---- Diagnostic plot ----
    if (TESTING && runif(1) > 0.9) {
      plot(D.X, O.D.Y * scale, type = 'l',
           xlim = c(x.lim.min, x.lim.max), ylim = c(ymin, ymax),
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
  if (fixed.false.positives > 0)
    WT <- c(fixed.false.positives, WT * (1 - fixed.false.positives))

out =  list(
    weights = WT,
    ncp     = auto$par[idx.ncp],
    zsds    = auto$par[idx.zsd],
    fit     = auto$objective
  )

out

}

############################


boot.old.fashioned <- function(val.input, boot.iter = 10000, ncores = NULL, seed = 42,
                                NCP.FIXED = TRUE, ZSDS.FIXED = TRUE) {


  #seed = 42

  if (is.null(ncores)) ncores <- max(1, parallel::detectCores() *.8)

  ncores <- max(1, cores)

  ## ---- Pre-filter to fitting range ----
  INT <- val.input[val.input >= Int.Beg & val.input <= Int.End]

  ## ---- Precompute Dens only when both fixed ----
  if (NCP.FIXED && ZSDS.FIXED) {
    densy     <- Get.Densities(INT, bw = bw.est, from = Int.Beg,
                               to = Int.End, Augment = Augment)
    D.X.pre   <- densy[, 1]
    Dens.pre  <- build.Dens(D.X.pre, ncp, zsds, df, CURVE.TYPE)
    bar.width <- D.X.pre[2] - D.X.pre[1]
    Dens.pre  <- Dens.pre / (rowSums(Dens.pre) * bar.width)
  } else {
    D.X.pre  <- NULL
    Dens.pre <- NULL
  }

  ## ---- Point estimate ----
  point.est <- old.fashioned(INT,
                             Dens.precomputed = Dens.pre,
                             D.X.precomputed  = D.X.pre,
                             bw.est=.1,width = 1,
                             NCP.FIXED = NCP.FIXED,
                             ZSDS.FIXED = ZSDS.FIXED)
 
  ## ---- Parallel bootstrap ----
  cl <- parallel::makeCluster(ncores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(cl,
    varlist = c(
      # Functions
      "old.fashioned", 
      "build.Dens", 
      "Get.Densities",
      # Data
      "INT",
      # Parameters
      "ncp", 
      "zsds", 
      "df", 
      "bw.est",
      "Int.Beg", 
      "Int.End", 
      "Augment",
      "CURVE.TYPE", 
      "fixed.false.positives",
      # Plotting params
      "x.lim.min", 
      "x.lim.max", 
      "ymin", 
      "ymax",
      "TESTING",
      # Bootstrap-specific
      "NCP.FIXED", 
      "ZSDS.FIXED",
      "Dens.pre", 
      "D.X.pre"
    ),
    envir = environment()
  )

  invisible(parallel::clusterEvalQ(cl, library(KernSmooth)))
  invisible(parallel::clusterEvalQ(cl, TESTING <- FALSE))
  parallel::clusterSetRNGStream(cl, iseed = seed)

  boot_res <- parallel::parLapply(cl, 1:boot.iter, function(i) {
    boot_sample <- sample(INT, length(INT), replace = TRUE)
    fit <- tryCatch(
      old.fashioned(boot_sample,
                    Dens.precomputed = Dens.pre,
                    D.X.precomputed  = D.X.pre,
                    bw.est = bw.est, width = width,
                    NCP.FIXED  = NCP.FIXED,
                    ZSDS.FIXED = ZSDS.FIXED)
	,error = function(e) NULL )
    fit
  })

  ## ---- Drop failed iterations ----
  boot_res <- boot_res[!sapply(boot_res, is.null)]
  if (length(boot_res) == 0) stop("All bootstrap iterations failed")

  ## ---- Collect ----
  boot.weights <- do.call(rbind, lapply(boot_res, `[[`, "weights"))
  boot.ncp     <- do.call(rbind, lapply(boot_res, `[[`, "ncp"))
  boot.zsds    <- do.call(rbind, lapply(boot_res, `[[`, "zsds"))
  boot.fit     <- sapply(boot_res, `[[`, "fit")

  ## ---- CIs ----
  ci.fun <- function(mat, point) {
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
    data.frame(
      estimate = point,
      lo95     = apply(mat, 2, quantile, 0.025),
      hi95     = apply(mat, 2, quantile, 0.975))
  }

  list(
    point  = point.est,
    ci     = list(weights = ci.fun(boot.weights, point.est$weights),
                  ncp     = ci.fun(boot.ncp,     point.est$ncp),
                  zsds    = ci.fun(boot.zsds,    point.est$zsds)),
    boot   = list(weights = boot.weights,
                  ncp     = boot.ncp,
                  zsds    = boot.zsds,
                  fit     = boot.fit),
    B.ok   = length(boot_res),
    B      = B,
    ncores = ncores
  )
}


### NEW OLD FASHIONED

run.new.OF <- function(val.input, boot.iter = 0, NCP.FIXED = TRUE, ZSDS.FIXED = TRUE,bw.est=.2) {

    INT <- val.input[val.input > Int.Beg & val.input < Int.End]

    ## ---- Precompute Dens if both fixed ----
    if (NCP.FIXED && ZSDS.FIXED) {
      densy     <- Get.Densities(INT, bw = bw.est, from = Int.Beg,
                                 to = Int.End, Augment = Augment)
      D.X.pre   <- densy[, 1]
      Dens.pre  <- build.Dens(D.X.pre, ncp, zsds, df, CURVE.TYPE)
      bar.width <- D.X.pre[2] - D.X.pre[1]
      Dens.pre  <- Dens.pre / (rowSums(Dens.pre) * bar.width)
    } else {
      D.X.pre  <- NULL
      Dens.pre <- NULL
    }


  ## ---- Point estimate ----
  res.run <- old.fashioned(INT,
                      Dens.precomputed = Dens.pre,
                      D.X.precomputed  = D.X.pre,
					bw.est = bw.est, width = width,
                      NCP.FIXED  = NCP.FIXED,
                      ZSDS.FIXED = ZSDS.FIXED)

  cp.input <- list(
    w.inp    = res.run$weights,
    ncp.est  = res.run$ncp,
    zsds.est = res.run$zsds
  )

  cp.res <- Compute.Power.Z.General(cp.input, Int.Beg = Int.Beg)

  res.pe <- list(
    EDR     = cp.res$EDR,
    ERR     = cp.res$ERR,
    ncp     = res.run$ncp,
    zsds    = res.run$zsds,
    w.inp   = res.run$weights,
    w.sig   = cp.res$w.sig,
    w.all   = cp.res$w.all,
    local.power = cp.res$local.power,
    fit     = res.run$fit
  )

  #res.pe

  ######################
  if (boot.iter > 0) {   #boot.iter = 50

  #print("Do New OF bootstrap")

    ## ---- Precompute Dens if both fixed ----
    if (NCP.FIXED && ZSDS.FIXED) {
      densy     <- Get.Densities(INT, bw = bw.est, from = Int.Beg,
                                 to = Int.End, Augment = Augment)
      D.X.pre   <- densy[, 1]
      Dens.pre  <- build.Dens(D.X.pre, ncp, zsds, df, CURVE.TYPE)
      bar.width <- D.X.pre[2] - D.X.pre[1]
      Dens.pre  <- Dens.pre / (rowSums(Dens.pre) * bar.width)
    } else {
      D.X.pre  <- NULL
      Dens.pre <- NULL
    }

    ## ---- Set up cluster ----
    ncores <- max(1, parallel::detectCores() * .8)
    cl <- parallel::makeCluster(ncores)
    on.exit(parallel::stopCluster(cl), add = TRUE)

    parallel::clusterExport(cl,
      varlist = c("old.fashioned", "build.Dens", "Get.Densities",
                  "INT",
                  "ncp", "zsds", "df", "bw.est",
                  "Int.Beg", "Int.End", "Augment",
                  "CURVE.TYPE", "fixed.false.positives",
                  "x.lim.min", "x.lim.max", "ymin", "ymax",
                  "bw.est",
                  "NCP.FIXED", "ZSDS.FIXED",
                  "Dens.pre", "D.X.pre"),
      envir = environment())

    invisible(parallel::clusterEvalQ(cl, library(KernSmooth)))
    invisible(parallel::clusterEvalQ(cl, TESTING <- FALSE))
    parallel::clusterSetRNGStream(cl, iseed = 42)

    ## ---- Bootstrap ----
    boot_res <- parallel::parLapply(cl, 1:boot.iter, function(i) {
      boot_sample <- sample(INT, length(INT), replace = TRUE)
      fit <- tryCatch(
        old.fashioned(boot_sample,
                     Dens.precomputed = Dens.pre,
                     D.X.precomputed  = D.X.pre,
                     bw.est = bw.est, width = width,
                     NCP.FIXED  = NCP.FIXED,
                     ZSDS.FIXED = ZSDS.FIXED)
		,error = function(e) NULL
      )
      fit
    })

    #boot_res

    ## ---- Drop failures ----
    boot_res <- boot_res[!sapply(boot_res, is.null)]

    ## ---- Reformat to match EM structure for Compute.Power.Z.General ----
    boot_res <- lapply(boot_res, function(x) {
      list(w.inp    = x$weights,
           ncp      = x$ncp,
           zsds     = x$zsds,
           fit      = x$fit)
    })

	boot_res

    ## ---- Compute power for each bootstrap result ----
    boot_power <- lapply(boot_res, function(x)
      Compute.Power.Z.General(x, Int.Beg = Int.Beg)
    )

    lcbind <- function(list1, list2) {
      stopifnot(length(list1) == length(list2))
      lapply(seq_along(list1), function(i) c(list1[[i]], list2[[i]]))
    }
    boot_list <- lcbind(boot_res, boot_power)

    lquantile <- function(boot_list, probs = c(0.025, 0.975)) {
      fields <- names(boot_list[[1]])
    
     result <- lapply(fields, function(f) {
      vals <- lapply(boot_list, function(x) x[[f]])
      if (all(sapply(vals, is.null))) return(NULL)
      mat <- do.call(rbind, vals)
      apply(mat, 2, quantile, probs = probs, na.rm = TRUE)
    })
  
    names(result) <- fields
    result[!sapply(result, is.null)]
  }

    get.ci <- function(boot_list, probs = c(0.025, 0.975)) {
      ci <- lquantile(boot_list, probs = probs)
    }

    res.ci <- get.ci(boot_list)

  } else {
    components <- length(ncp)
    res.ci <- list(
      w.inp   = matrix(NA, 2, components),
      w.all   = matrix(NA, 2, components),
      w.sig   = matrix(NA, 2, components),
      ncp     = matrix(NA, 2, components),
      zsds    = matrix(NA, 2, components),
      fit     = matrix(NA, 2, 1),
      EDR     = matrix(NA, 2, 1),
      ERR     = matrix(NA, 2, 1),
      local.power = matrix(NA, 2, length(res.pe$local.power))
    )
  }

  combine.pe.ci <- function(res.pe, res.ci) {
    fields <- names(res.pe)
    result <- lapply(fields, function(f) {
      pe   <- res.pe[[f]]
      ci_f <- res.ci[[f]]
      if (is.null(dim(ci_f))) {
        rbind(estimate = pe, ci_lb = ci_f[1], ci_ub = ci_f[2])
      } else {
        rbind(estimate = pe, ci_lb = ci_f[1, ], ci_ub = ci_f[2, ])
      }
    })
    names(result) <- fields
    result
  }

  zcurve.res <- combine.pe.ci(res.pe, res.ci)
  zcurve.res$EDR[2] <- zcurve.res$EDR[2] - EDR.CI.adjust
  zcurve.res$EDR[3] <- zcurve.res$EDR[3] + EDR.CI.adjust
  zcurve.res$ERR[2] <- zcurve.res$ERR[2] - ERR.CI.adjust
  zcurve.res$ERR[3] <- zcurve.res$ERR[3] + ERR.CI.adjust

  return(zcurve.res)
}

##############################################################
##############################################################
##############################################################


### function to detect bias

test.bias = function(w.all) {

#Int.Beg = 2 + .4

D.X = seq(x.lim.min,x.lim.max,.01)
n.bars = length(D.X)
bar.width = D.X[2]-D.X[1]

components = length(ncp)
means = ncp
sds = zsds
w = w.all


### get the densities for each interval and each non-centrality parameter
Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:components) {
		if (CURVE.TYPE == "z") {
			Dens = c(Dens,
			dnorm(D.X[i],means[j],sds[j]) +
			dnorm(-D.X[i],means[j],sds[j]) 
			)
		} else {
			sds = sds - 1
			sds[sds <= 0] = 0
			sds
			if (max(sds) > 0) {
				Dens = c(Dens,
				dt(D.X[i],df,means[j]) +
				dt(-D.X[i],df,means[j])
				)
			} else { 
				Dens = c(Dens,
				dt(D.X[i],df,means[j]) +
				dt(-D.X[i],df,means[j])
				)
			}
		}
	}
}

Dens = matrix(Dens,length(ncp),byrow=FALSE)

### rescale the densities for the range of z-values to 1
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)
dim(Dens)

### compute the new estimated density distribution
E.D.Y = colSums(Dens*w)
length(E.D.Y)

E.D.Y = E.D.Y / sum(E.D.Y * bar.width)
sum(E.D.Y * bar.width)

prob1 = sum(E.D.Y[D.X > crit & D.X < crit + just]*bar.width);prob1
prob2 = sum(E.D.Y[D.X > crit + just & D.X < Int.End]*bar.width);prob2
prob3 = sum(E.D.Y[D.X < crit]*bar.width);prob3
sum(prob1,prob2,prob3)
prob = prob1 / (prob1 + prob2)
prob
adj = .00
prob.adj = prob + adj 

k = length(val.input[val.input < Int.End])

sig.k = length(val.input[val.input > crit & val.input < Int.End])
sig.k

#Int.Beg = 0

just.sig.k = length(val.input[val.input > crit & val.input < crit + just])
just.sig.k

just.sig.k/sig.k

### binomial
p.bias.binomial = 1 - pbinom(just.sig.k - 1, sig.k, prob = prob.adj);p.bias.binomial


### for all
#prob = prob1 / (prob1 + prob2 + prob3)
#just.sig.k/k
#p.bias.binomial = 1 - pbinom(just.sig.k - 1, k, prob = prob.adj);p.bias.binomial


### chi2 approximation
#O = just.sig.k
#E = sig.k*prob;E
#chi2.val <- (O - E)^2 / (E * (1 - prob));chi2.val
#p.bias.chi2 <- 1 - pchisq(chi2.val, df = 1)

bias.res = c(just.sig.k/sig.k,prob.adj,p.bias.binomial)

return(bias.res)

} #EOF bias test

#Test Bias
#test.bias(w.all)


### function to run the zcurve package to get EM weights and fit

run.zcurve = function(Est.Method = "OF",kd.model="KD2",K=6,
	alpha = .05,Int.Beg = 1.96,Int.End = 6, boot.iter = 0,parallel=TRUE) {

	meth = Est.Method
	if (Est.Method == "OF") meth = "density"

	biter = boot.iter
	if (boot.iter == 0) biter = FALSE

	z.res = zcurve(val.input,bootstrap=biter,method=meth,
		control=list(parallel = parallel,
		sig_level=alpha,a = Int.Beg,b = Int.End,mu=ncp))

	return(z.res)

} ### End run.zcurve

if (TESTING) {
	meth = "density"
	z.res = run.zcurve(val.input,Est.Method=meth,boot.iter=0)
	z.res
	#str(z.res)
	round(z.res$fit$weights,3)
	w.inp = z.res$fit$weights
	round(w.inp,3)
}

############################################
############################################
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

  mtext(lab, side = 1, line = 1.8, at = midpoints, cex = 1.0, las = 1)
  #mtext(lab, side = 1, line = 2, at = midpoints, cex = 1.0, las = 1)

  par(mar = old_mar)

} ### EOF Write.Local.Power



###################################################
#### Begin Draw Histogram
###################################################

Draw.Histogram = function(w,cola="blue3",
	Write.CI = FALSE) {

	#w = w.all;cola = "blue3"; 

	int.start = round(Int.Beg,1)
	if (round(Int.Beg,2) == 1.96) int.start = 2

	z.hist = val.input[val.input > x.lim.min & val.input < x.lim.max]
	if (round(Int.Beg,2) == 1.96) {
		z.hist = val.input[val.input > x.lim.min & val.input < x.lim.max -.04] + .04
	}
	table(z.hist > 1.96)

	n.breaks = seq(x.lim.min,x.lim.max,hist.bar.width);n.breaks

	par(cex.axis = 1)
	par(cex.lab = 1)
	par(family = "sans") 
	par(font.axis = 1) #"Microsoft Himalaya")
	par(font.lab = 1) # Microsoft Himalaya")

	col1 = adjustcolor(cola, alpha.f = 0.2)
	col2 = adjustcolor(cola, alpha.f = 0.3)

	hist(z.hist,breaks=n.breaks,freq=FALSE,
		col=col1,border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		ylab="Density",xlab="",axes=FALSE,main=Title,lwd=1)

	axis(2, ylim = c(0,ymax))

	axis(1, line = 2)  # draw x-axis lower

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

	ymax.scale = (ymax+ymin)*scale

	hist(z.hist[z.hist > round(Int.Beg,1) & z.hist < Int.End],
		breaks=n.breaks,freq=FALSE,
		col=col2,border="white",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax.scale),
		ylab=,xlab="",main=Title,lwd=1,axes=FALSE)



######################################### 
######################################### 

if (Show.Text) {

#par(family = fam[2])


par(new=TRUE)
hist(c(0),main="",ylim=c(ymin,ymax),ylab="",xlab="",xlim=c(x.lim.min,x.lim.max),
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
	text(x.lim.min,ymax-y.line*i,"z-curve 3.0 (2025)",pos=4,cex=letter.size)
	i = i + 1.5
	text(x.lim.min,ymax-y.line*i,version,pos=4,cex =letter.size)

	########

	#TESTING

	i = 0
	results.x = x.lim.max 

#	text(results.x,y.text-y.line*i,paste0("Range: ",
#		round(min.z,2)," to ",round(max.z,2)),
#		pos=2,cex=letter.size)
#	i = i + 2

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z,big.mark=",")),
		" tests "),pos=2,cex = letter.size)
	i = i + 1.5

	if (substring(Est.Method,1,3) == "CLU") {
		text(results.x+0.1,y.text-y.line*i,paste0(toString(format(cluster.k,big.mark=",")),
			" articles "),pos=2,cex = letter.size)
		i = i + 1.5
	}

	text(results.x+0.1,y.text-y.line*i,paste0(toString(format(n.z.sig,big.mark=","))," significant "),
	  pos=2,cex = letter.size)
	i = i + 1.5
	text(results.x,y.text-y.line*i,paste0(n.not.shown," z > ",x.lim.max),
		pos=2,cex = letter.size)
	

#	par(family = fam[1])

	#############################################
	#############################################


	if (Write.CI) {

	if(TEST4BIAS & !is.null(bias) ) {
		if (bias[3] < .00005) { bias.res = "EJS, p < .0001" } else {  
			bias.res = paste0("EJS, p = ",sub("^0","",formatC(bias[3],format="f",digits=4))) }
	bias.res

	}	

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

	if (ODR.high > 100) ODR.high = 100
	if(ERR.low < alpha/2) ERR.low = alpha/2
	if(EDR.low < alpha) EDR.low = alpha
	
	i = i + 3
	text(results.x,y.text-y.line*i,
		paste0(round((1-CI.ALPHA)*100),"% CI     "),
		pos=2,cex=letter.size.1)

	x.loc = results.x + c(-1.5,-1.1,-1.0,-.6,-.5,-.1,0)

	ODR.a = str_pad(round(ODR), width = 3, pad = " ");ODR.a
	ODR.b = str_pad(round(ODR.low), width = 3, pad = " ");ODR.b
	ODR.c = str_pad(round(ODR.high), width = 3, pad = " ");ODR.c

	EDR.a = str_pad(round(EDR), width = 3, pad = " ");EDR.a
	EDR.b = str_pad(round(EDR.low), width = 3, pad = " ");EDR.b
	EDR.c = str_pad(round(EDR.high), width = 3, pad = " ");EDR.c

	ERR.a = str_pad(round(ERR), width = 3, pad = " ");ERR.a
	ERR.b = str_pad(round(ERR.low), width = 3, pad = " ");ERR.b
	ERR.c = str_pad(round(ERR.high), width = 3, pad = " ");ERR.c

	FDR.a = str_pad(round(FDR), width = 3, pad = " ");FDR.a
	FDR.b = str_pad(round(FDR.low), width = 3, pad = " ");FDR.b
	FDR.c = str_pad(round(FDR.high), width = 3, pad = " ");FDR.c

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

	if(TEST4BIAS) {
		i = i + 10
		text(results.x,y.text-y.line*i,results$bias[3],
			pos=2,cex=letter.size.1)
		}

	##############################################
	# End of IF CI show results with CI
	} else {  
	# End of IF CI show results without CI

	
	ODR = results$ODR[1]*100;ODR
	EDR = results$EDR[1]*100;ERR
	ERR = results$ERR[1]*100;EDR
	FDR = results$FDR[1]*100;FDR

	if(TEST4BIAS & !is.null(bias) ) { 
		if (results$bias[3] < .00005) { bias.res = "EJS, p < .0001" } else { 
			bias.res = paste0("EJS, p = ",sub("^0","",formatC(bias[3],format="f",digits=4))) }
		bias.res
	}	

	i = i + 3
	ODR.string = str_pad(round(ODR), width = 3, pad = " ");ODR.string
	text(results.x,y.text-y.line*i,
		paste("ODR:",ODR.string,"%"),
		pos=2,cex=letter.size.1)

	i = i + 2
	EDR.string = str_pad(round(EDR), width = 3, pad = " ");EDR.string
	text(results.x,y.text-y.line*i,
		paste("EDR:",EDR.string,"%"),
		pos=2,cex=letter.size.1)

	i = i + 2
	ERR.string = str_pad(round(ERR), width = 3, pad = " ");ERR.string
	text(results.x,y.text-y.line*i,
		paste("ERR:",ERR.string,"%"),
		pos=2,cex=letter.size.1) 

	i = i + 2
	FDR.string = str_pad(round(FDR), width = 3, pad = " ");FDR.string
	text(results.x,y.text-y.line*i,
		paste("FDR:",FDR.string,"%"),
		pos=2,cex=letter.size.1) 

	if(TEST4BIAS) {
		i = i + 2
		text(results.x,y.text-y.line*i,bias.res,
			pos=2,cex=letter.size.1)
	}


#	i = i + 2
#	text(results.x,y.text-y.line*i,
#		paste0("OSR:   ",OSR,"%"),
#		pos=2,cex=letter.size.1)

	i = i + 2

	if (Est.Method == "OF" & boot.iter > 0 & Write.CI == FALSE) text(results.x,ymax-y.line*i,
		"WAIT FOR BOOSTRAPPED CIs",
		pos=2,cex=letter.size,col="red")

	i = i + 2
	if (TEST4HETEROGENEITY) text(results.x,ymax-y.line*i,
		"WAIT FOR HETEROGENEITY BOOSTRAPPED CIs",
		pos=2,cex=letter.size,col="red")



	} # End of CI

} # End of Show.Text


abline(h=0)


} # End of Histogram

 
######################################### 
######################################### 
######################################### 


Draw.KD = function(z.draw,w,Write.CI=FALSE,cola="blue",Lwidth=5) {

	#ymin=-.015;ymax = .6
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

	d.fit = length(z.draw[z.draw > Int.Beg & z.draw < x.lim.max]) /
		length(z.draw[z.draw > x.lim.min & z.draw < x.lim.max])

	d.fit

	d.sig.Y = d.sig.Y * d.fit

	sum(d.sig[,2])*bar.width
	
#	### draw the density of the observed values in the selected region
#	par(new=TRUE)
#	plot(d.all.X[d.all.X > x.lim.min+.05 & d.all.X < x.lim.max - .05],
#		d.all.Y[d.all.X > x.lim.min+.05 & d.all.X < x.lim.max - .05],
#		type="l",col=adjustcolor(cola, alpha.f = 0.3),lty=1,lwd=Lwidth,
#		xlim =c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),xlab="",ylab="")


	### draw the density of the observed values in the selected region
	par(new=TRUE)
	plot(d.sig.X[d.sig.X > Int.Beg +.05 & d.sig.X < x.lim.max - .05],
		d.sig.Y[d.sig.X > Int.Beg + .05 & d.sig.X < x.lim.max - .05],
		type="l",col=adjustcolor(cola, alpha.f = 0.8),lty=2,lwd=Lwidth,
		xlim =c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),xlab="",ylab="",
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
	Ltype=1,x.start=x.lim.min,x.end=x.lim.max) {


#x.start = x.lim.min; x.end = x.lim.max;cola = "black";Lwidth=4;Ltype=1

ncp = results$ncp[1,]
zsds = results$zsds[1,]
w = results$w.all[1,]

width = 1

bar.width = .01
D.X = seq(x.lim.min,x.lim.max,bar.width) 

D.Y = build.Dens(D.X, ncp, zsds, df, CURVE.TYPE) 
D.Y = D.Y / (rowSums(D.Y) * bar.width)
D.Y = colSums(D.Y * w)

#plot(D.X,D.Y)

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
	lty=Ltype,col=cola,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))

boundary = Int.Beg
if (round(Int.Beg,2) == 1.96) {boundary = 2}

segments(boundary,0,boundary,
	D.Y[which(round(D.X,2) == round(boundary,2))[1]]*scale,
	col=cola,lwd=3)


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		sig.level = crit
		if(two.sided == FALSE) sig.level = c(-sig.level,sig.level)

		sig.lty = rep(2,length(sig.level));sig.lty

		i = 1
		for (i in 1:length(sig.level)) {
			height = max(scale*D.Y[which(round(D.X,1) == round(sig.level[i],1))]);height
			segments(sig.level[i],0,sig.level[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance


} # End of Draw.Curve.All


###################################################
###################################################
###################################################



#########################################################################
### This Function Computes Power from Weights and Non-Centrality Parameters
#########################################################################

### Claude.26.03.14

Compute.Power.Z.Discrete = function(cp.input, Int.Beg = 1.96) {

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

  w.all = w.all.ext[1:components]

  # EDR: average power over all studies (pre-selection)
  EDR = sum(w.all.ext * pow.ext)

  # ERR: average directional power over significant studies
  w.sig.ext = w.all.ext * pow.ext
  w.sig.ext = w.sig.ext / sum(w.sig.ext)
  ERR = sum(w.sig.ext * pow.dir.ext)

  ERR = min(max(ERR, alpha / 2), 1)
  EDR = min(max(EDR, alpha),     1)

  EDR
  ERR


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

out = (list(
    EDR     = EDR,
    ERR     = ERR,
    w.all   = w.all.ext[1:components],
    w.sig   = w.sig.ext[1:components],
    local.power = local.power
  ))


out


} ### EOF Compute.Power.Z.Discrete


#####

Compute.Power.Z.Continuous = function(cp.input, Int.Beg) {

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

  EDR
  ERR


  local.power = NULL

if (int.loc > 0) {

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
    EDR     = EDR,
    ERR     = ERR,
    w.all   = w.all,
    local.power = local.power
  )

out

} ### EOF Compute.Power.Z.Continuous

#####

Compute.Power.Z.General = function(cp.input, Int.Beg = 1.96) {
  if (max(as.numeric(cp.input$zsds)) < 1.01) {
    Compute.Power.Z.Discrete(cp.input, Int.Beg)
  } else {
    Compute.Power.Z.Continuous(cp.input, Int.Beg)
  }

}


#####################################
#####################################
#####################################

#For quick testing
#val.input = rnorm(1000,2,2)

#####################################
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
### BB ZingStart #START #Begin of Main Program 
#####################################

if (length(p) > 0) {
  val.input = qnorm(1-p/2)
}

if (substring(Est.Method,1,3) == "CLU") {
	if (length(cluster.id > 0)) {
	   clu.inp = data.frame(val.input,cluster.id)
	   clu.inp
	   summary(clu.inp)
	   clu.inp = clu.inp[complete.cases(clu.inp),]
	   z.clu.input = zcurve_data(paste("z = ", clu.inp[,1]), id = clu.inp[,2])
	   tab = table(cluster.id);tab
	   cluster.k = length(tab);cluster.k
	} else{
       print("Cluster IDs are missing!")
	}
} 

if (CURVE.TYPE == "z") { 
  crit = qnorm(alpha/2,lower.tail=FALSE)
 } else {
  crit = qt(alpha/2,df,lower.tail=FALSE)
}

if (length(Int.Beg) == 0) Int.Beg = crit


components = length(ncp)


### remove NA
val.input = val.input[!is.na(val.input)];length(val.input)


if (two.sided) val.input = abs(val.input)

### limit range
val.input[val.input > MAX.INP.Z] = MAX.INP.Z
length(val.input)

if (length(val.input) > 10) print("!START!") else print("Insufficient Data")

print(paste0("Title: ",Title))

### create set with z-scores in the interval used for model fitting
INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
summary(INT)

n.z = length(val.input);n.z
n.z.sig = length(val.input[abs(val.input) > crit])

ODR = n.z.sig/n.z;ODR

ODR.se = sqrt(ODR * (1 - ODR) / n.z)
ODR.se

ODR.low = round((ODR-1.96*ODR.se),3);ODR.low
ODR.high = round((ODR+1.96*ODR.se),3);ODR.high


ODR.res = c(ODR,ODR.low,ODR.high)
ODR.res


extreme = c(
	length(val.input[val.input < -Int.End])/length(val.input[val.input < -crit]),
	length(val.input[val.input > Int.End])/length(val.input[val.input > crit])
)
names(extreme) = c("Ext.Neg","Ext.Pos")
(extreme*100)


slope = get_slope(INT)
slope

w_start = NULL


### End of Preparation
###############################################

### BBB 

###############################################
###############################################

### NEW OLD FASHIONED Z-CURVE 

if (Est.Method == "OF") {
	res = run.new.OF(val.input,boot.iter = boot.iter,bw=bw.est,NCP.FIXED,ZSDS.FIXED)
} # EOF Est.Method OF

### NEW EM Z-Curve 
if(Est.Method == "EM") {

	cp.input = run.new.OF(val.input,boot.iter = 0)

	floor = .1
	w_from_density <- cp.input$w.inp[1,]
	w_start <- w_from_density + 0.1  # small floor
	w_start <- w_start / sum(w_start)  # renormalize
	w_start

	components = length(ncp)	

	res = run.new.EM(val.input,w_start,boot.iter,NCP.FIXED,ZSDS.FIXED)

}


res$FDR = (1/res$EDR - 1)*(alpha/(1-alpha))
res$FDR = res$FDR[c(1,3,2)]

bias = NULL
if(TEST4BIAS) { 
  bias = test.bias(res$w.all[1,]) 
  names(bias) = c("OBS.JS","EXP.JS","EJS.p")
}

results = list(
		slope = slope,
		ODR = ODR.res,
		EDR = res$EDR,
		ERR = res$ERR,
		FDR = res$FDR,
		ncp = res$ncp,
		zsds = res$zsds,
		w.inp = res$w.inp,
		w.all = res$w.all,	
		local.power = res$local.power,
		bias = bias,
		fit = res$fit
	  )


#results

####################################################################

if (Est.Method %in% c("CLU", "CLU-W","CLU-B") & boot.iter >= 0) {

	bb = Int.End
    if (bb > 6) bb = 6
	
	print("USING CLUSTER OPTION BOOTSTRAP")

	method = "b"

    if (Est.Method == "CLU-W") print("CLU-W NOT SUPPORTED") #method = "w"

    z.clu = zcurve_clustered(z.clu.input,method = method, bootstrap=boot.iter,
            control=list(method="EM",a = Int.Beg,b = bb,
            sig_level=alpha,max_iter = max_iter,max_iter_boot = max_iter_boot))

	summary(z.clu)     

	fit = z.clu$fit$Q

	w.inp =	 summary(z.clu, type="parameters")$coefficients[(components+1):(2*components)]
	round(w.inp,3)

	cp.res = Compute.Power.Z.General(Int.Beg = Int.Beg,c(w.inp,ncp,zsds))
	round(cp.res,3)

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	round(w.all,3)
	sum(w.all)

	local.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(local.power,3)

	#str(summary(z.clu))

	EDR = summary(z.clu)$coefficients[2,];EDR
	ERR = summary(z.clu)$coefficients[1,];ERR
	FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2)[c(1,3,2)];FDR

} # EOF Cluster Method 



##########################################
### This Code is Used to Create Graphic
##########################################

if (Show.Histogram) {

	if (boot.iter == 0) {
		#print("Draw Histogram NO CI")
		Draw.Histogram(results$w.all[1,],cola=col.hist,Write.CI = FALSE)
	} else {
		#print("Draw Histogram WITH CI")
		Draw.Histogram(results$w.all[1,],cola=col.hist,Write.CI = TRUE)
	}

	if (Show.Curve.All) { 

		Draw.Curve.All(results,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(results,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)

	} # EOF Show.Curve.All


	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)


	if (length(results$local.power) > 0 && !is.na(results$local.power[1])) Write.Local.Power(results$local.power[1,])	


} # End of Show.Histogram	for bootstrap


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

