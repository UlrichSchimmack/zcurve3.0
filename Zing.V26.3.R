
#rm(list = ls())

########################################################################
### SETTING PARAMETERS FOR Z-CURVE MODEL
########################################################################

version <- "Version 2026.03.15"   # Version label to appear on plots

# Optional cleanup 
# rm(list = ls())
# options(scipen = 999)  # Disable scientific notation

### INSTALL PACKAGES (only once – manually run if needed)
if (1 == 2) {  # This block is ignored unless manually changed to (1 == 1)
  install.packages("pwr")
  install.packages("zcurve")
  install.packages("KernSmooth")
  install.packages("parallel")
  install.packages("stringr")
} # END install block

### LOAD LIBRARIES
library(parallel)
library(KernSmooth)
library(zcurve)
library(stringr)
library(pwr)


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
bw.draw <- 0.10              # Smoothing for kernel density display

### CONSOLE OUTPUT

Show.Iterations <- FALSE      # Do not use with parallel Show iterations for slow procedures (e.g., EXT, TEST4HETEROGENEITY)

### MODEL PARAMETERS

alpha <- 0.05                        # Significance level
crit <- qnorm(1 - alpha / 2)       # Corresponding two-sided critical z

two.sided <- TRUE                   # Assume two-sided z-values (use abs(z)); not yet compatible with signed z-values

# Color scheme
col.curve <- "violetred3"
col.hist <- "blue3"
col.kd <- "green3"

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
Augment.Regression <- FALSE         # Use Slope for Augmentation
Augment.Factor <- 1                 # Amount of augmentation

bw.est <- 0.05                      # Bandwidth for kernel density (lower = less smoothing, higher = more smoothing)
bw.aug <-  .20					 # Width of Augmentation interval

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
                               B          = 1000,
                               NCP.FIXED  = TRUE,
                               ZSDS.FIXED = TRUE,
                               W.FIXED    = FALSE,
                               n_starts   = 10,
                               cores      = round(parallel::detectCores() * .7)) {

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

  boot_res <- parallel::parLapply(cl, 1:B, function(i) {

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
      loglik = fit$loglik
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


run.new.zcurve = function(val.input, w_start, NCP.FIXED=TRUE,ZSDS.FIXED = TRUE) {


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
   loc.pow = cp.res$loc.pow,
   loglik = res.run$loglik
)

######################

if (boot.iter > 0) {

boot_res <- run_bootstrap_list(INT = INT,
                          ncp.start = ncp,
                          zsds.start = zsds,
                          Int.Beg = Int.Beg,
                          Int.End = Int.End,		
                          B = boot.iter,
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
    mat <- do.call(rbind, lapply(boot_list, function(x) x[[f]]))
    apply(mat, 2, quantile, probs = probs, na.rm = TRUE)
  })
  
  names(result) <- fields
  result
}

get.ci <- function(boot_list, probs = c(0.025, 0.975)) {
  ci <- lquantile(boot_list, probs = probs)
}

res.ci = get.ci(boot_list)

} else {

res.ci = list(
 w.inp = matrix(NA,2,components), 
 w.all = matrix(NA,2,components),
 w.sig = matrix(NA,2,components),
 ncp = matrix(NA,2,components),
 zsds = matrix(NA,2,components),
 loglik = matrix(NA,2,1),
 EDR = matrix(NA,2,1),
 ERR = matrix(NA,2,1),
 loc.pow = matrix(NA,2,length(res.pe$loc.pow))
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

boot.iter = 0

zcurve.res = run.new.zcurve(val.input,NCP.FIXED,ZSDS.FIXED)

zcurve.res.3d = round_list(zcurve.res)
zcurve.res.3d

}



##################################################################

### Functions used in Zing

#######################################################################


############################################################
############################################################
############################################################
############################################################
### SLOPE DIAGNOSTICS

adj_density <- function(zval,
                        Int.Beg = 1.96,
                        z.crit = 1.96,
                        bw = 0.1,
                        aug.bw = 0.2,
                        n.bars = 512,
                        max.z = NULL,
                        d.x.max = Int.Beg + 1) {
  
  zval <- abs(zval[is.finite(zval)])
  if (length(zval) == 0) return(NULL)
  
  if (is.null(max.z)) max.z <- max(zval)
  
  if (Int.Beg < max.z) {
    
    AUG <- c()
    
    if (Int.Beg >= z.crit + 2*bw) {
      AUG <- zval[zval > (Int.Beg - 2*bw) & zval < Int.Beg]
    }
    
    if (Int.Beg < z.crit + 2*bw) {
      n.AUG <- round(length(zval[zval > Int.Beg & zval < (Int.Beg + aug.bw)]))
      if (n.AUG > 0) {
        AUG <- seq(Int.Beg - aug.bw, Int.Beg - 0.01, aug.bw / n.AUG)
      }
    }
    
    Z.INT.USE <- c(zval, AUG)
    
  } else {
    # If Int.Beg is beyond the data, there is nothing to estimate
    return(NULL)
  }
  
  density(Z.INT.USE, n = n.bars, bw = bw, from = Int.Beg, to = d.x.max)
}


slope_adj <- function(zval, Int.Beg = 1.96, width = 1, min_in_window = 10, ...) {
  
  zval <- abs(zval[is.finite(zval)])
  k.slope <- sum(zval >= Int.Beg & zval <= (Int.Beg + width))
  
  slope <- NA_real_
  
  if (k.slope >= min_in_window) {
    
    d <- adj_density(zval, Int.Beg = Int.Beg, d.x.max = Int.Beg + width, ...)
    if (!is.null(d)) {
      
      idx <- d$x < (Int.Beg + width) & is.finite(d$y)
      
      if (sum(idx) >= 5) {
        fit <- lm(d$y[idx] ~ d$x[idx])
        slope <- unname(coef(fit)[2])
      }
    }
  }
  
  c(k.slope = k.slope, slope = slope)
}

############################################################
### END OF SLOPE DIAGNOSTICS




get.t.density = function(tx) {     

n.bars = length(tx);n.bars
bar.width = tx[2]-tx[1]

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncp)) {
		Dens = c(Dens,
		dt(tx[i],df,ncp[j]) +
		dt(-tx[i],df,ncp[j])
		)
	}
}
Dens = matrix(Dens,length(ncp),byrow=FALSE)
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)

return(Dens)

} # EOF get.t.density 

####################################

get.z.density = function(zx) {     

n.bars = length(zx);n.bars
bar.width = zx[2] - zx[1]

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncp)) {
		Dens = c(Dens,
		dnorm(zx[i],ncp[j],zsds[j]) +
		dnorm(-zx[i],ncp[j],zsds[j])
		)
	}
}
Dens = matrix(Dens,length(ncp),byrow=FALSE)
row.sum.dens = rowSums(Dens)
Dens = Dens/(row.sum.dens * bar.width)
dim(Dens)

return(Dens)

} #EOF get.z.density


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

get.ci.info = function(Est.Method = "EM") {

#val.input = abs(c(rnorm(1000,1.1),rnorm(1000,2.8)))
#boot.iter = 0

if (Est.Method %in% c("OF","EM")) {

#Int.Beg = 2
zres = run.zcurve(val.input, Est.Method=Est.Method, alpha = alpha,boot.iter = boot.iter,
	Int.Beg = Int.Beg, Int.End = Int.End,parallel=parallel)

ODR = ODR.res
ERR = summary(zres)$coefficients[1,];ERR
EDR = summary(zres)$coefficients[2,];EDR
FDR = (1/EDR - 1)*(alpha/(1-alpha));FDR
FDR = FDR[c(1,3,2)];FDR

# bootstrap selected weights (B x 7)
Wboot <- zres$boot$weights

pow_sel <- function(ncp, Int.Beg) {
  pnorm(ncp, mean = Int.Beg) + pnorm(-ncp, mean = Int.Beg)  # Φ(ncp-c)+Φ(-ncp-c)
}

w_all_from_w_inp <- function(w_inp, ncp, Int.Beg) {
  ps <- pow_sel(ncp, Int.Beg)
  x  <- w_inp / ps
  x / sum(x)
}

# bootstrap distribution of w.all
w_all_boot <- t(apply(Wboot, 1, function(wrow) {
  w_all_from_w_inp(wrow, ncp, Int.Beg)
}))

# point estimate from the fitted selected weights (replace with your object’s slot)
w_inp_hat <- zres$fit$weights
w_all_hat <- w_all_from_w_inp(w_inp_hat, ncp, Int.Beg)

# percentile 95% CIs
w_all_ci <- t(apply(w_all_boot, 2, quantile, probs = c(.025, .975), na.rm = TRUE))
colnames(w_all_ci) <- c("low", "high")

w.all = round(cbind(estimate = w_all_hat, w_all_ci), 3)
w.all

fit.val = summary(zres)$model$fit_index
fit.val

res.ci = list(
	slope = slope,
    ODR = ODR,
    EDR = EDR,
    ERR = ERR,
    FDR = FDR,
    ncp.est = ncp,
    zsds.est = zsds,
    w.all = w.all,
	bias = bias,
    fit = fit.val
)

#res.ci

}
 

if (Est.Method == "EXT") {

res.ext = EXT.boot(ZSDS.FIXED=ZSDS.FIXED)
res.ext

w.all.with.ci = cbind(results$w.all,res.ext[6:(5+components),1],
		res.ext[6:(5+components),2])
w.all.with.ci

res.ci = list(
	slope = slope,
    ODR = ODR.res,
    EDR = c(results$EDR,res.ext[1,]),
    ERR = c(results$ERR,res.ext[2,]),
    FDR = c(results$FDR,res.ext[3,]),
    ncp.est = c(results$ncp,res.ext[4,]),
    zsds.est = c(results$zsds,res.ext[5,]),
    w.all = w.all.with.ci,
	bias = bias,
    fit = c(results$fit,res.ext[nrow(res.ext),])
)

#res.ci

}

return(res.ci)

}

#######################################################################


Write.Local.Power = function(loc.power) {

  if (int.loc == 0) return(invisible(NULL))

  # Midpoints of each bin — these are the correct label positions
  midpoints = seq(x.lim.min, x.lim.max - int.loc, by = int.loc) + int.loc / 2

  # Format labels
  lab = paste0(round(loc.power * 100), "%")
#  lab = paste0(format(round(loc.power * 100), nsmall = 0), "%")

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


	if (is.na(bias[3])) TEST4BIAS = FALSE

	if(TEST4BIAS) { 
		if (bias[3] < .00005) { bias.res = "EJS, p < .0001" } else {  
			bias.res = paste0("EJS, p = ",sub("^0","",formatC(bias[3],format="f",digits=4))) }
		bias.res
	}	

	ODR = results$ODR[1]*100
	ODR.low = results$ODR[2]*100
	ODR.high = results$ODR[3]*100
	EDR = results$EDR[1]*100
	EDR.low = results$EDR[2]*100
	EDR.high = results$EDR[3]*100
	ERR = results$ERR[1]*100
	ERR.low = results$ERR[2]*100
	ERR.high = results$ERR[3]*100
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

	if (is.na(results$bias[3])) TEST4BIAS = FALSE

	if(TEST4BIAS) { 
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

	if (is.na(bias[3])) TEST4BIAS = FALSE

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
		bw=bw.draw,d.x.min=x.lim.min,d.x.max=x.lim.max,Augment=Augment)
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
		bw=bw.draw,d.x.min=Int.Beg,d.x.max=x.lim.max,Augment=Augment)
	
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
		type="l",col=adjustcolor(cola, alpha.f = 0.4),lty=1,lwd=Lwidth,
		xlim =c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),xlab="",ylab="",
		axes=FALSE)
	### draw vertical line for Beginning

	Int.Beg
	height = max(d.sig.Y[which(round(d.sig.X,1) == round(Int.Beg+.02,1))]);height
	segments(Int.Beg+.02,0,Int.Beg+.02,height,lty=1,lwd=3,col=cola)




} # End of Draw.KD

###################################################
#### End Draw.KD
###################################################



###################################################
#### Begin Draw Zcurve
###################################################


Draw.Curve.All = function(w,cola="black",Lwidth=4,
	Ltype=1,x.start=x.lim.min,x.end=x.lim.max) {


#Draw.Histogram()

bar.width = .01
D.X = seq(x.lim.min,x.lim.max,bar.width)

if (CURVE.TYPE == "z") {
  Dens = get.z.density(D.X)
} else {
  Dens = get.t.density(D.X)
}

D.Y = colSums(Dens*w)

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
#### Begin Draw Zcurve SDG1
###################################################

Draw.Curve.All.SDG1 = function(w.inp,ncp=ncp,zsds=zsds,cola=col.curve,
	Ltype=1,Lwidth=4,x.start=x.lim.min,x.end=x.lim.max) {


#results$w.inp
if(length(w.inp) > components) w.inp = w.inp[1,]
if(length(ncp) > components) ncp = ncp[1,]
if(length(zsds) > components) zsds = zsds[1,]

components = length(ncp)

if (components == 1) { 
	D.X = seq(x.lim.min,x.lim.max,.01)
	if (CURVE.TYPE == "z") {
		D.Y = dnorm(D.X,ncp,zsds)
	} else {
		D.Y = dt(D.Y,df,ncp)
	}
	d.sim = cbind(D.X,D.Y)
} else {
    print("WARNING NOT SUPPOSED TO DO THIS")
    print(components)
    print(w.inp)
    print(k.sim)
    print(ncp)
    print(zsds)
}
dim(d.sim)

bar.width = d.sim[2,1] - d.sim[1,1];bar.width
d.sim = d.sim[d.sim[,1] > x.lim.min & d.sim[,1] <= x.lim.max,]
dim(d.sim)
d.sim[,2] = d.sim[,2]/(sum(d.sim[,2])*bar.width)
sum(d.sim[,2])*bar.width

d.sim.X = d.sim[,1]
d.sim.Y = d.sim[,2]

#plot(d.sim.X,d.sim.Y,ylim=c(ymin,ymax))

summary(d.sim.X)
sum(d.sim.Y*bar.width)

d.hist.sel = mean(as.numeric(val.input > x.lim.min & val.input > Int.Beg & val.input < Int.End)) /
	mean(as.numeric(val.input > x.lim.min & val.input < Int.End))
d.hist.sel

d.dense.sel = sum(d.sim.Y[d.sim.X > x.lim.min & d.sim.X >= Int.Beg]*bar.width)
d.dense.sel

scale = d.hist.sel/d.dense.sel;scale

lines(d.sim.X[which(d.sim.X >= x.start & d.sim.X < x.end)],
	d.sim.Y[which(d.sim.X >= x.start & d.sim.X < x.end)]*scale,
	lty=Ltype,col=cola,lwd=Lwidth,xlim=c(x.lim.min,x.lim.max)
	,ylim=c(ymin,ymax))


if (Show.Significance) {

	if (length(sig.levels) == 0) sig.levels = alpha

		cz.alpha = -qnorm(sig.levels/2);cz.alpha
		if(two.sided == FALSE) cz.alpha = c(-cz.alpha,cz.alpha);cz.alpha

		sig.lty = rep(2,length(cz.alpha));sig.lty

		i = 1
		for (i in 1:length(cz.alpha)) {
			cz.alpha[i]
			height = max(d.sim.Y[which(round(d.sim.X,1) == round(cz.alpha[i],1))]);height
			segments(cz.alpha[i],0,cz.alpha[i],height,lty=sig.lty,lwd=2,col="firebrick3")
		}

} # End of Show.Significance



} # End of Draw.Curve.SDG1

###################################################




##############################################
### Get Densities
##############################################

Get.Densities = function(INT,bw="nrd0",d.x.min=0,d.x.max=6,Augment=TRUE) {

### find the maximum z-score. This is only needed if the maximum z-score is below Int.End
max.z = Int.End
if (max(INT) < d.x.max) max.z = max(INT)
max.z

if(Augment.Regression & slope[1] > 20 & Int.Beg > 1) { 
	Augment.Factor = 1 - d.reg;Augment.Factor
} 

### Augment z-scores on the left side of Interval to avoid downward trend 
### of kernal density function (assumes values go to 0)

if (Augment == TRUE) { 

	AUG = c()
	n.AUG = round(Augment.Factor*length(INT[INT > d.x.min & INT < d.x.min+bw.aug]));n.AUG
	if (n.AUG > 0) AUG = seq(d.x.min-bw.aug,d.x.min,bw.aug/n.AUG)

	Z.INT.USE = c(INT,AUG)

} else {
	Z.INT.USE = INT[INT > d.x.min & INT <= max.z]
}

Z.Density = bkde(Z.INT.USE,bandwidth=bw,range=c(d.x.min-bw.aug,d.x.max)) 
val.max = d.x.max
D = data.frame(Z.Density$x,Z.Density$y)
colnames(D) = c("ZX","ZY")
#plot(D$ZX,D$ZY);abline(v=d.x.min)


D = D[D$ZX > d.x.min & D$ZX < val.max,]
dim(D)
#plot(D$ZX,D$ZY)

bar.width = D$ZX[2] - D$ZX[1]
D$ZY = D$ZY/(sum(D$ZY*bar.width)) 
sum(D$ZY*bar.width)

return(D)


}  ### End of Get Densities 


#######################################################
### End of Get Densities
#######################################################



#######################################################
### Old Fashioned Density Method
#######################################################

#weight = startval

old.fashioned = function(val.input,cola = "springgreen2") {

curve.fitting = function(theta,RetEst=FALSE)    {

### get the weights and rescale 
weight = theta
weight = weight/sum(weight)

if (fixed.false.positives > 0) weight = c(fixed.false.positives,weight*(1-fixed.false.positives))
sum(weight)

### compute the new estimated density distribution
E.D.Y = colSums(Dens*weight)

### compare to observed density distribution
rmse = sqrt(mean((E.D.Y-O.D.Y)^2))

### return either fit if continue or estimates if finished
value = rmse
if(RetEst) value = E.D.Y


### showing the fitting of the function in a plot
if(Plot.Fitting) {

	rval = runif(1)
	if (rval > .9) {

	tit = ""
	xL = ""
	yL = ""
	plot(D.X,O.D.Y*scale,type='l',
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		main=tit,xlab=xL,ylab=yL,axes=FALSE)
	lines(D.X,E.D.Y*scale,lty=1,col="red1",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),axes=FALSE)
	#points(D.X,z.est,pch=20,col="red1",ylim=c(0,ymax),)

	}

} ### End of Plot Fitting

### return value to optimization function
return(value)

} ### End of curve.fitting

############################
#### End of Fitting Function
############################

components = length(ncp)

### get the densities for each interval and each non-centrality parameter

INT = val.input[val.input >= Int.Beg & val.input <= Int.End]

densy = Get.Densities(INT,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment)

D.X = densy[,1]
O.D.Y = densy[,2]

#plot(D.X,O.D.Y,type="l")

n.bars = length(D.X)

bar.width = D.X[2] - D.X[1]
bar.width

### Finish getting observed densities 

Dens	= c()
for(i in 1:n.bars) {
	for (j in 1:length(ncp)) {
		if (CURVE.TYPE == "z") {
			Dens = c(Dens,
			dnorm(D.X[i],ncp[j],zsds[j]) +
			dnorm(-D.X[i],ncp[j],zsds[j])
			)
		} else {
			Dens = c(Dens,
			dt(D.X[i],df,ncp[j]) +
			dt(-D.X[i],df,ncp[j])
			)
		}
	}
}
summary(Dens)
Dens = matrix(Dens,length(ncp),byrow=FALSE)
sum.dens = rowSums(Dens)
Dens = Dens/(sum.dens * bar.width)
dim(Dens)

startval = rep(1/(components),components)
#startval[1] = 1
startval = startval/sum(startval)
startval

lowlim = c(rep(0,components),ncp,zsds);lowlim
highlim = c(rep(1,components),ncp,zsds);highlim

if (fixed.false.positives > 0 & 0 %in% ncp) {
  startval = rep(1/(components-1),components-1)
  lowlim = c(rep(0,components-1),ncp,zsds);lowlim
  highlim = c(rep(1,components-1),ncp,zsds);highlim
}


d.hist.sel = mean(as.numeric(val.input > x.lim.min & val.input > Int.Beg & val.input < Int.End)) /
	mean(as.numeric(val.input > x.lim.min & val.input < Int.End))
d.hist.sel

scale = d.hist.sel

#TESTING = TRUE
if (TESTING) Plot.Fitting = TRUE else Plot.Fitting = FALSE


auto = nlminb(startval,curve.fitting,lower=lowlim,upper=highlim,
	control=list(eval.max=1000,abs.tol = 1e-20))

fit = auto$objective;fit

### get the estimated weights 
WT = auto$par
WT = WT/sum(WT)
if (fixed.false.positives > 0) WT = c(fixed.false.positives,WT*(1-fixed.false.positives))
sum(WT)

res = c(WT,fit);res

return(res)

} ### End of old.fashioned

######################################################
### End of Old Fashioned Density Method
#######################################################



#########################################################################
### This Function Computes Power from Weights and Non-Centrality Parameters
#########################################################################

Compute.Power.Z.Old = function(cp.input,Int.Beg=crit) {

crit = qnorm(alpha/2,lower.tail=FALSE)

ext.all = length(val.input[val.input > Int.End]) / 
	length(val.input)
ext.all

ext.sig = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > crit])
ext.sig

ext.inp = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > Int.Beg])
ext.inp

w.inp = cp.input$w.inp

ncp = cp.input$ncp
zsds = cp.input$zsds

components = length(ncp)

### get power values for the components (ncp)
pow.dir = pnorm(abs(ncp),crit);pow.dir 

### get the opposite sign probability
sign.error = 1-pnorm(ncp,-crit);round(sign.error,3)

pow = pow.dir + sign.error
pow.ext = c(pow,1)

pow.sel = pnorm(ncp,Int.Beg) + pnorm(-ncp,Int.Beg);pow.sel
pow.sel.ext = c(pow.sel,1)

w.inp.ext = c(w.inp*(1-ext.inp), ext.inp)          # estimated weights

w.all.ext   <- (w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext)) / sum(w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext))

EDR = sum(w.all.ext*pow.ext);EDR

w.sig.ext = w.all.ext * pow.ext
w.sig.ext = w.sig.ext / EDR
sum(w.sig.ext)

pow.dir.ext = c(pow.dir,1)

w.all = w.all.ext[1:components]
w.sig = w.sig.ext[1:components]


ERR = sum(w.sig.ext*pow.dir.ext);ERR

if (ERR > 1) ERR = 1
if (EDR > 1) EDR = 1

if (ERR < alpha/2) ERR = alpha/2
if (EDR < alpha) EDR = alpha 


ERR.pos = NA
ERR.neg = NA

EDR.pos = NA
EDR.neg = NA


### now we compute mean power as a function of z-scores continuously
### this is only performed if local power is requested (int.loc > 0)

local.power = NA

if (int.loc > 0) {

	bar.width = .01 # how fine should be the resolution
	X = seq(x.lim.min,Int.End,bar.width);summary(X)	 # set of z-scores 

	W.D.Sum = unlist(lapply(X,function(x) sum(dnorm(x,ncp)*w.all) ))
	#plot(Z.X,Z.W.D.Sum)

	loc.p = unlist(lapply(X,function(x) 
		sum(dnorm(x,ncp)*w.all*(pow.dir+sign.error)) /	sum(dnorm(x,ncp)*w.all)	))

	int = seq(x.lim.min,x.lim.max,int.loc)
	local.power = c()
	i = 1
	for (i in 1:(length(int)-1)) local.power = c(local.power,
		sum(loc.p[X > int[i] & X < int[i+1]]*
			W.D.Sum[X > int[i] & X < int[i+1]])/
		sum(W.D.Sum[X > int[i] & X < int[i+1]])	 )		

	local.power
	names(local.power) = paste0("lp.",1:length(local.power))
}

local.power

res = list(
   EDR = EDR,
   ERR = ERR,
   w.sig = w.sig,
   w.all = w.all,
   loc.pow = local.power
)

### to be past back to the main program from this function
return(res)

}

### Finish Compute.Power.Z

#####################################
#####################################
#####################################

Compute.Power.T = function(para,Int.Beg=crit,BOOT=FALSE) {

#para = para.est.OF

t.crit = qt(1-alpha/2,df);t.crit

t.ext.all = length(val.input[val.input > Int.End]) / 
	length(val.input)
t.ext.all

t.ext.sig = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > crit])
t.ext.sig

t.ext.inp = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > Int.Beg])
t.ext.inp

w.inp = para[1:components]
nct = para[(1+components):(2*components)]

### get power values for the components (ncp)
pow.dir = pt(t.crit,df,nct,lower.tail=FALSE);pow.dir 

### get the opposite sign probability
sign.error = pt(-t.crit,df,nct,lower.tail=TRUE);round(sign.error,3)

pow = pow.dir + sign.error
pow.ext = c(pow,1)

pow.sel = pt(Int.Beg,df,nct,lower.tail=FALSE) + pt(-Int.Beg,df,nct);pow.sel
pow.sel.ext = c(pow.sel,1)
pow.sel

w.inp.ext = c(w.inp*(1-t.ext.inp), t.ext.inp)          # estimated weights
w.inp.ext

w.all.ext   <- (w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext)) / sum(w.inp.ext + w.inp.ext*(1-pow.sel.ext)
   /(pow.sel.ext))

EDR = sum(w.all.ext*pow.ext);EDR

w.sig.ext = w.all.ext * pow.ext
w.sig.ext = w.sig.ext / EDR
sum(w.sig.ext)

pow.dir.ext = c(pow.dir,1)

w.all = w.all.ext[1:components]
w.sig = w.sig.ext[1:components]

ERR = sum(w.sig.ext*pow.dir.ext);ERR

if (ERR > 1) ERR = 1
if (EDR > 1) EDR = 1

if (ERR < alpha/2) ERR = alpha/2
if (EDR < alpha) EDR = alpha 


ERR.pos = NA
ERR.neg = NA

EDR.pos = NA
EDR.neg = NA

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,w.all,w.sig)
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncp[1:components]),paste0("w.sig",ncp[1:components]) )
res


local.power <- c()
for (i in 1:(length(int)-1)) {
  segment <- (ncp.grid >= int[i]) & (ncp.grid < int[i+1])
  denom <- sum(dd[segment])
  local.power <- c(local.power, if (denom == 0) NA else sum(power[segment] * dd[segment]) / denom)
}

local.power = NA

if (int.loc > 0) {

	bar.width = .01 # how fine should be the resolution
	Z.X = seq(x.lim.min,Int.End,bar.width);summary(Z.X)	 # set of z-scores 

	Z.W.D.Sum = unlist(lapply(Z.X,function(x) sum(dt(x,df,ncp)*w.all) ))
	#plot(Z.X,Z.W.D.Sum)

	loc.p = unlist(lapply(Z.X,function(x) 
		sum(dt(x,df,ncp)*w.all*(pow.dir+sign.error)) /	
			sum(dt(x,df,ncp)*w.all)	))

	loc.p
	int = seq(x.lim.min,x.lim.max,int.loc)
	local.power = c()
	i = 1
	for (i in 1:(length(int)-1)) local.power = c(local.power,
		sum(loc.p[Z.X > int[i] & Z.X < int[i+1]]*
			Z.W.D.Sum[Z.X > int[i] & Z.X < int[i+1]])/
		sum(Z.W.D.Sum[Z.X > int[i] & Z.X < int[i+1]])	 )		

	local.power
	names(local.power) = paste0("lp.",1:length(local.power))
	res = c(res,local.power)
}

return(res)

} # EOF Compute.Power.T


#####################################
#####################################
#####################################

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

local.power

####

  return(list(
    EDR     = EDR,
    ERR     = ERR,
    w.all   = w.all.ext[1:components],
    w.sig   = w.sig.ext[1:components],
    loc.pow = local.power
  ))

} ### EOF Compute.Power.Z.Discrete


#####

Compute.Power.Z.Continuous = function(cp.input, Int.Beg) {

  deci   = 2
  zx.bw  = 1 / (10^deci)
  zx     = seq(0, Int.End, zx.bw)

  components = length(cp.input$ncp)
  ncp        = round(cp.input$ncp, deci)
  ncp.sd     = max(1, sqrt(cp.input$zsds^2 - 1))

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


  # Local power
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

  local.power

  return(list(
    EDR     = EDR,
    ERR     = ERR,
    w.all   = w.all,
    loc.pow = local.power
  ))

} ### EOF Compute.Power.Z.Continuous

#####

Compute.Power.Z.General = function(cp.input, Int.Beg = 1.96) {
  if (max(cp.input$zsds) < 1.01) {
    Compute.Power.Z.Discrete(cp.input, Int.Beg)
  } else {
    Compute.Power.Z.Continuous(cp.input, Int.Beg)
  }
}


#####################################
#####################################
#####################################



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


slope = slope_adj(val.input)
slope

w_start = NULL


### End of Preparation
###############################################

### BBB 

###############################################
###############################################


### OLD FASHIONED Z-CURVE 

### run always to get starting values, unless one component model

if(components > 1 & ZSDS.FIXED == TRUE & NCP.FIXED == TRUE) {

#summary(val.input)
para.est.OF = old.fashioned(val.input)

fit = para.est.OF[components+1];fit

cp.input  = list(
	w.inp = para.est.OF[1:components],
	ncp.est = ncp,
	zsds.est = zsds
)


if (CURVE.TYPE == "z") {
  cp.res = Compute.Power.Z.General(cp.input,Int.Beg=Int.Beg)
} else {
  cp.res = Compute.Power.T(cp.input,Int.Beg=Int.Beg)
}

EDR = cp.res$EDR;EDR
ERR = cp.res$ERR;ERR

FDR = (1/EDR - 1)*(alpha/(1-alpha))

w.all = cp.res$w.all

w.sig = cp.res$w.sig

loc.power = cp.res$loc.pow

bias = c(NA,NA,NA)
if(TEST4BIAS) { 
	bias = test.bias(w.all) 
}
names(bias) = c("OBS.JS","EXP.JS","EJS.p")

results = list(
	slope = slope,
	ODR = ODR,
	EDR = EDR,
	ERR = ERR,
	FDR = FDR,
    ncp = cp.input$ncp.est,
    zsds = cp.input$zsds.est,
    w.all = w.all,
    bias = bias,
    fit = fit
  )

floor = .1
w_from_density <- cp.input$w.inp
w_start <- w_from_density + 0.1  # small floor
w_start <- w_start / sum(w_start)  # renormalize
w_start


#print(results)

} # EOF Est.Method OF
 


##########################################################

 
### if requested, run slower EM
### Est.Method = "EM"
if(Est.Method == "EM") {

	components = length(ncp)	
	#summary(val.input)
	#Int.Beg = 2.4; alpha = .005

	res.em = run.zcurve(val.input,Est.Method="EM",boot.iter=boot.iter,alpha=alpha,
		Int.Beg=Int.Beg,Int.End=Int.End,parallel=parallel)

	res.em

	fit = res.em$fit$Q
	fit

	w.inp = res.em$fit$weights

	cp.input  = list(
		w.inp = res.em$fit$weights,
		ncp.est = ncp,
		zsds.est = zsds
	)

	cp.input

	cp.res = Compute.Power.Z.General(cp.input,Int.Beg=Int.Beg)
	cp.res


	summary(res.em)     

	EDR = cp.res$EDR;EDR
	ERR = cp.res$ERR;ERR

	w.all = cp.res$w.all

	w.sig = cp.res$w.sig

	loc.power = cp.res$loc.pow

	FDR = (1/EDR - 1)*(alpha/(1-alpha));FDR

bias = c(NA,NA,NA)
if(TEST4BIAS) { 
	bias = test.bias(w.all) 
}
names(bias) = c("OBS.JS","EXP.JS","EJS.p")


results = list(
	slope = slope,
    ODR = ODR,
	EDR = EDR,
	ERR = ERR,
	FDR = FDR,
    ncp = cp.input$ncp.est,
    zsds = cp.input$zsds.est,
    w.all = w.all,
    bias = bias,
    fit = fit
  )


#results


} 


#####################################################
#####################################################
#####################################################

### Est.Method = "NEW" ### New EM algorithm
if(Est.Method == "NEW") {

	components = length(ncp)	

	res.new = run.new.zcurve(val.input,w_start,NCP.FIXED,ZSDS.FIXED)

	res.new

	w.inp = res.new$w.inp[1,]

	EDR = res.new$EDR[1];EDR
	ERR = res.new$ERR[1];ERR

	FDR = (1/res.new$EDR - 1)*(alpha/(1-alpha))
	FDR = FDR[c(1,3,2)]

	w.all = if (is.null(dim(res.new$w.all))) res.new$w.all else res.new$w.all[1,]

	loc.power = if (is.null(dim(res.new$loc.pow))) res.new$loc.power else res.new$loc.pow[1,]

	bias = c(NA,NA,NA)
	if(TEST4BIAS) { 
		bias = test.bias(w.all) 
	}


	names(bias) = c("OBS.JS","EXP.JS","EJS.p")

	results.new = list(
		slope = slope,
		ODR = ODR.res,
		EDR = res.new$EDR,
		ERR = res.new$ERR,
		FDR = FDR,
		ncp = res.new$ncp,
		zsds = res.new$zsds,
		w.inp = res.new$w.inp,
		w.all = res.new$w.all,	
		bias = bias,
		fit = res.new$loglik
	  )

	results = list(
		slope = slope,
		ODR = ODR.res,
		EDR = res.new$EDR[1],
		ERR = res.new$ERR[1],
		FDR = FDR[1],
		ncp = res.new$ncp[1,1],
		zsds = res.new$zsds[1,1],
		w.all = res.new$w.all[1,1],	
		bias = bias,
		fit = res.new$loglik[1,1]
	  )

	results = results.new


} 

#results
#results.new

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

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(loc.power,3)

	#str(summary(z.clu))

	EDR = summary(z.clu)$coefficients[2,];EDR
	ERR = summary(z.clu)$coefficients[1,];ERR
	FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2)[c(1,3,2)];FDR

} # EOF Cluster Method 


########################################################## 



##########################################
### This Code is Used to Create Graphic
##########################################

#CCC

if (Show.Histogram & sum(extreme,na.rm=TRUE) < .95) { 

	#print("Show Histogram")

	Draw.Histogram(w.all,cola=col.hist)

	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)

	zsds.check = max(results$zsds[1,])

	if (Show.Curve.All & zsds.check < 1.05) {

		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)

		}

	if (Show.Curve.All & zsds.check > 1.05) {

		print("Showing SDG1 PLOT")	

		Draw.Curve.All.SDG1(w=w.all,ncp=results$ncp,zsds=results$zsds,cola=col.curve,
			Ltype = 3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All.SDG1(w=w.all,ncp=results$ncp,zsds=results$zsds,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (length(loc.power) > 0 && !is.na(loc.power[1])) Write.Local.Power(loc.power)	

	} # End of Show.Histogram



########################################################################
### code for confidence intervals 
########################################################################

### If Confidence Intervals are requested, compute CI (boot.iter > 0)
if (boot.iter > 0 & Est.Method %in% c("OF","EM","EXT") )  {

	#boot.iter = 50
	#Est.Method = "EM"

	res.with.ci = get.ci.info(Est.Method = Est.Method)

	results = res.with.ci

}


if (boot.iter > 0 & substring(Est.Method,1,3) == "CLU") {

	results = rbind(c(ODR,ODR.low,ODR.high),EDR,ERR,FDR)
	rownames(results) = c("ODR","EDR","ERR","FDR")
	round(results,3)

}


if (boot.iter > 0 & Est.Method == "NEW") results = results.new

zsds.check = max(results$zsds[1,])

if (boot.iter > 0 & Show.Histogram) {

	Draw.Histogram(w.all,cola=col.hist,Write.CI = TRUE)

	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)

	if (Show.Curve.All & zsds.check < 1.05) {
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (Show.Curve.All & zsds.check > 1.05) {
		Draw.Curve.All.SDG1(w=w.all,ncp=results$ncp,zsds=results$zsds,cola=col.curve,
			Ltype = 3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All.SDG1(w=w.all,ncp=results$ncp,zsds=results$zsds,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}


	if (length(loc.power) > 0 && !is.na(loc.power[1])) Write.Local.Power(loc.power)	


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

