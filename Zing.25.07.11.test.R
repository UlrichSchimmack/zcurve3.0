
########################################################################
### SETTING PARAMETERS FOR Z-CURVE MODEL
########################################################################

version <- "Version 2026.02.20"   # Version label to appear on plots

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

ZSDS.FIXED <- FALSE                 # Fix SD values for EXT method 
NCP.FIXED <- FALSE                  # Fix non-central parameter(NCP) means values for EXT method
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


Compute.Power.General <- function(cp.input, Int.Beg = 1.96) {

  w.inp    <- cp.input$w.inp
  ncp.est  <- cp.input$ncp.est
  zsds.est <- cp.input$zsds.est

  components <- length(ncp.est)

  # --------------------------------
  # Expected power per component
  # --------------------------------

  beta <- numeric(components)

  for (k in 1:components) {

    a1 <- ( Int.Beg - ncp.est[k]) / zsds.est[k]
    a2 <- (-Int.Beg - ncp.est[k]) / zsds.est[k]

    beta[k] <- (1 - pnorm(a1)) + pnorm(a2)
  }

  beta[beta < 1e-12] <- 1e-12

  # --------------------------------
  # Pre-selection weights
  # --------------------------------

  w.all <- w.inp / beta
  w.all <- w.all / sum(w.all)

  # --------------------------------
  # Expected discovery rate
  # --------------------------------

  EDR <- sum(w.all * beta)

  # --------------------------------
  # Weights among significant studies
  # --------------------------------

  w.sig <- w.all * beta / EDR

  # --------------------------------
  # Expected replication rate
  # --------------------------------

  ERR <- sum(w.sig * beta)

  # --------------------------------
  # Local power
  # --------------------------------

  loc.pow <- beta
  names(loc.pow) <- paste0("lp.", 1:components)

  list(
    EDR     = EDR,
    ERR     = ERR,
    w.sig   = w.sig,
    w.all   = w.all,
    loc.pow = loc.pow
  )
}




round_list <- function(x, digits = 3) {
  if (is.numeric(x)) {
    round(x, digits)
  } else if (is.list(x)) {
    lapply(x, round_list, digits = digits)
  } else {
    x
  }
}



########################################################
######## NEW EM ALGORITHM

update_ncp_newton_folded <- function(
  mu_init,
  sigma,
  n_k,
  S1_k,
  Int.Beg,
  Int.End,
  max_iter = 20,
  tol = 1e-8,
  lower = -10,
  upper = 10
) {

  mu <- mu_init

  for (iter in 1:max_iter) {

    # truncation limits
    a1 <- (Int.Beg - mu) / sigma
    b1 <- (Int.End - mu) / sigma
    a2 <- (-Int.End - mu) / sigma
    b2 <- (-Int.Beg - mu) / sigma

    Phi_a1 <- pnorm(a1)
    Phi_b1 <- pnorm(b1)
    Phi_a2 <- pnorm(a2)
    Phi_b2 <- pnorm(b2)

    Cmu <- (Phi_b1 - Phi_a1) +
           (Phi_b2 - Phi_a2)

    if (Cmu <= 0 || !is.finite(Cmu))
      return(mu)   # fail safely

    phi_a1 <- dnorm(a1)
    phi_b1 <- dnorm(b1)
    phi_a2 <- dnorm(a2)
    phi_b2 <- dnorm(b2)

    delta <- (phi_a1 - phi_b1 +
              phi_a2 - phi_b2) / Cmu

    kappa <- (a1 * phi_a1 - b1 * phi_b1 +
              a2 * phi_a2 - b2 * phi_b2) / Cmu

    # score
    U <- (S1_k - n_k * mu) / sigma^2 -
         n_k * delta / sigma

    # hessian
    H <- - n_k / sigma^2 -
         n_k / sigma^2 * (delta^2 + kappa)

    step <- U / H
    mu_new <- mu - step

    mu_new <- max(lower, min(upper, mu_new))

    if (abs(mu_new - mu) < tol)
      break

    mu <- mu_new
  }

  return(mu)
}


update_zsds_newton_folded <- function(
  sigma_init,
  mu,
  n_k,
  S1_k,
  S2_k,
  Int.Beg,
  Int.End,
  max_iter = 20,
  tol = 1e-8
) {

  eta <- log(sigma_init)

  for (iter in 1:max_iter) {

    sigma <- exp(eta)

    a1 <- (Int.Beg - mu) / sigma
    b1 <- (Int.End - mu) / sigma
    a2 <- (-Int.End - mu) / sigma
    b2 <- (-Int.Beg - mu) / sigma

    Phi_a1 <- pnorm(a1)
    Phi_b1 <- pnorm(b1)
    Phi_a2 <- pnorm(a2)
    Phi_b2 <- pnorm(b2)

    C <- (Phi_b1 - Phi_a1) +
         (Phi_b2 - Phi_a2)

    if (C <= 0 || !is.finite(C))
      return(sigma)

    phi_a1 <- dnorm(a1)
    phi_b1 <- dnorm(b1)
    phi_a2 <- dnorm(a2)
    phi_b2 <- dnorm(b2)

    delta <- (phi_a1 - phi_b1 +
              phi_a2 - phi_b2) / C

    kappa <- (a1 * phi_a1 - b1 * phi_b1 +
              a2 * phi_a2 - b2 * phi_b2) / C

    Q2 <- S2_k -
          2 * mu * S1_k +
          n_k * mu^2

    # score in sigma
    U_sigma <- - n_k / sigma +
               Q2 / sigma^3 -
               n_k * ( (Int.Beg - mu)/sigma * delta +
                       kappa ) / sigma

    # convert to eta scale
    U_eta <- U_sigma * sigma

    # simple stable Hessian approximation
    H_eta <- -2 * n_k -
             3 * Q2 / sigma^2

    step <- U_eta / H_eta
    eta_new <- eta - step

    if (abs(eta_new - eta) < tol)
      break

    eta <- eta_new
  }

  return(max(0.05, exp(eta)))
}


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

  for (iter in 1:max_iter) {

    # ============================================================
    # E-STEP
    # ============================================================

    log_g <- matrix(NA_real_, k.int, components)

    for (k in 1:components) {

      mu    <- ncp[k]
      sigma <- zsds[k]

      # ----- folded density -----
      l1 <- dnorm(x,  mu, sigma, log = TRUE)
      l2 <- dnorm(x, -mu, sigma, log = TRUE)

      m  <- pmax(l1, l2)
      log_fold <- m + log(exp(l1 - m) + exp(l2 - m))

      # ----- two-sided truncation -----
      a1 <- ( Int.Beg - mu) / sigma
      b1 <- ( Int.End - mu) / sigma
      a2 <- (-Int.End - mu) / sigma
      b2 <- (-Int.Beg - mu) / sigma

      norm_const <-
        (pnorm(b1) - pnorm(a1)) +
        (pnorm(b2) - pnorm(a2))

	    print(c(mu, sigma))
	    print(norm_const)

 	  if (norm_const <= 0 || !is.finite(norm_const)) {
	    print("BAD NORM CONST")
	    stop("norm_const failure")
		}

      log_g[, k] <- log_fold - log(norm_const)

    }

    log_num <- sweep(log_g, 2, log(w.inp), "+")
    log_denom <- apply(log_num, 1, function(v) {
      m <- max(v)
      m + log(sum(exp(v - m)))
    })

    tau <- exp(log_num - log_denom)

    # ============================================================
    # M-STEP (weights only by default)
    # ============================================================


# ============================================================
# M-STEP
# ============================================================

# sufficient statistics
n_k  <- colSums(tau)
S1_k <- colSums(tau * x)
S2_k <- colSums(tau * x^2)

# ---- update weights ----
if (!W.FIXED) {
  w.inp <- n_k / k.int
  w.inp <- pmax(w.inp, 1e-12)
  w.inp <- w.inp / sum(w.inp)
}

# ---- update ncp ----
if (!NCP.FIXED) {

  for (k in 1:components) {

    if (n_k[k] < 1e-10) next  # guard zero-mass components

    ncp[k] <- update_ncp_newton_folded(
      mu_init = ncp[k],
      sigma   = zsds[k],
      n_k     = n_k[k],
      S1_k    = S1_k[k],
      Int.Beg = Int.Beg,
      Int.End = Int.End
    )
  }
}

# ---- update zsds ----
if (!ZSDS.FIXED) {

  for (k in 1:components) {

    if (n_k[k] < 1e-10) next

    zsds[k] <- update_zsds_newton_folded(
      sigma_init = zsds[k],
      mu         = ncp[k],
      n_k        = n_k[k],
      S1_k       = S1_k[k],
      S2_k       = S2_k[k],
      Int.Beg    = Int.Beg,
      Int.End    = Int.End
    )
  }
}



    # ============================================================
    # LOG-LIKELIHOOD
    # ============================================================

    loglik <- sum(log_denom)
    loglik_trace[iter] <- loglik

    if (!is.finite(loglik))
      return(NULL)


	if (max(w.inp) > 1 - 1e-8) {
	    if (abs(loglik - loglik_old) < tol) {
	        break
	    }
	}

    loglik_old <- loglik
  }

  list(
    w.inp.est    = w.inp,
    ncp.est      = ncp,
    zsds.est     = zsds,
    loglik       = loglik,
    loglik_trace = loglik_trace[1:iter],
    iter         = iter
  )
}



############################# NOT FOLDED BELOW




update_ncp_newton <- function(mu_init, sigma, n_k, S1_k, c,
                              max_iter = 20, tol = 1e-8,
                              lower = -10, upper = 10) {

  mu <- mu_init

  for (iter in 1:max_iter) {

    a <- (c - mu) / sigma
    denom <- 1 - pnorm(a)
    lambda <- dnorm(a) / denom

    # score
    U <- (S1_k - n_k * mu) / sigma^2 -
         n_k * lambda / sigma

    # hessian
    H <- - n_k / sigma^2 -
         n_k / sigma^2 * lambda * (lambda - a)

    step <- U / H
    mu_new <- mu - step

    # clamp to bounds
    mu_new <- max(lower, min(upper, mu_new))

    if (abs(mu_new - mu) < tol)
      break

    mu <- mu_new
  }

  return(mu)
}





EM_EXT <- function(INT,
                   ncp,
                   zsds = rep(1, length(ncp)),
                   w.inp,
                   Int.Beg = 1.96,
                   Int.End = 6,
                   components = length(ncp),
                   NCP.FIXED  = TRUE,
                   ZSDS.FIXED = TRUE,
                   W.FIXED    = FALSE,
                   max_iter = 200,
                   tol = 1e-6) {

  k.int <- length(INT)
  loglik_old <- -Inf
  loglik_trace <- numeric(max_iter)

  for (iter in 1:max_iter) {

    # ============================================================
    # E-STEP
    # ============================================================

	log_f <- {
	  l1 <- dnorm(x, mu, sigma, log = TRUE)
	  l2 <- dnorm(-x, mu, sigma, log = TRUE)
	  m  <- pmax(l1, l2)
	  m + log(exp(l1 - m) + exp(l2 - m))
	}

    log_g <- matrix(NA_real_, k.int, components)

    for (k in 1:components) {

		a <- (Int.Beg - ncp[k]) / zsds[k]
		b <- (Int.End - ncp[k]) / zsds[k]

		log_norm <- log(pnorm(b) - pnorm(a) )

		log_g[, k] <-
		  dnorm(INT, ncp[k], zsds[k], log = TRUE) -
		  log_norm

    }

    log_num <- sweep(log_g, 2, log(w.inp), "+")
    
    log_denom <- apply(log_num, 1, function(x) {
      m <- max(x)
      m + log(sum(exp(x - m)))
    })

    tau <- exp(log_num - log_denom)

    # ============================================================
    # SUFFICIENT STATISTICS
    # ============================================================

    n_k  <- colSums(tau)
    S1_k <- colSums(tau * INT)
    S2_k <- colSums(tau * INT^2)

    # ============================================================
    # M-STEP
    # ============================================================

    # ---- Update weights ----
    if (!W.FIXED) {
      w.inp <- n_k / k.int
    }

    # ---- Update ncp ----
    if (!NCP.FIXED) {

      for (k in 1:components) {

        mu <- ncp[k]
        sigma <- zsds[k]

        for (inner in 1:20) {


		a <- (Int.Beg - mu) / sigma
		b <- (Int.End - mu) / sigma

		phi_a <- dnorm(a)
		phi_b <- dnorm(b)

		Phi_a <- pnorm(a)
		Phi_b <- pnorm(b)

		denom <- Phi_b - Phi_a

		delta <- (phi_a - phi_b) / denom

          # Score
          U <- (S1_k[k] - n_k[k] * mu) / sigma^2 -
               n_k[k] * delta / sigma


		delta <- (phi_a - phi_b) / denom
		kappa <- (a * phi_a - b * phi_b) / denom

         # Hessian
		H <- - n_k[k] / sigma^2 -
		     n_k[k] / sigma^2 * (delta^2 + kappa)


          step <- U / H
          mu_new <- mu - step

          if (abs(mu_new - mu) < 1e-8)
            break

          mu <- mu_new
        }

        ncp[k] <- mu
      }
    }

    # ---- Update zsds ----
    if (!ZSDS.FIXED) {

      for (k in 1:components) {

        mu <- ncp[k]
        eta <- log(zsds[k])

        for (inner in 1:20) {

          sigma <- exp(eta)

          a <- (Int.Beg - mu) / sigma
          denom <- 1 - pnorm(a)
          lambda <- dnorm(a) / denom

          Q2 <- S2_k[k] -
                2 * mu * S1_k[k] +
                n_k[k] * mu^2

          # Score in sigma
          U_sigma <- - n_k[k] / sigma +
                     Q2 / sigma^3 -
                     n_k[k] * lambda *
                     (Int.Beg - mu) / sigma^2

          U_eta <- U_sigma * sigma

          # Hessian in eta
          H_eta <- -2 * n_k[k] -
                   3 * Q2 / sigma^2 +
                   n_k[k] * lambda *
                   (Int.Beg - mu) / sigma

          step <- U_eta / H_eta
          eta_new <- eta - step

          if (abs(eta_new - eta) < 1e-8)
            break

          eta <- eta_new
        }

        zsds[k] <- max(0.05, exp(eta))
      }
    }

    # ============================================================
    # LOG-LIKELIHOOD
    # ============================================================

    loglik <- sum(log_denom)

    loglik_trace[iter] <- loglik

    if (abs(loglik - loglik_old) < tol)
      break

    loglik_old <- loglik
  }


	if (!is.finite(loglik)) {
	  print(iter)
 	 stop("loglik became non-finite")
	}



  list(
    ncp.est      = ncp,
    zsds.est     = zsds,
    w.inp.est    = w.inp,
    loglik       = loglik,
    loglik_trace = loglik_trace,
    iter         = iter
   )


}


################################################

bootstrap_one <- function(vals,
                          ncp.start,
                          zsds.start,
                          Int.Beg,
						Int.End,
						w.start = rep(1/length(ncp.start), length(ncp.start)),
                          NCP.FIXED,
                          ZSDS.FIXED,
                          W.FIXED) {

  n <- length(vals)
  boot_sample <- sample(vals, n, replace = TRUE)


  fit <- EM_EXT_FOLDED(INT = boot_sample,
           ncp = ncp.start,
           zsds = zsds.start,
           w.inp = w.start,
           Int.Beg = Int.Beg,
		  Int.End = Int.End,
           NCP.FIXED = NCP.FIXED,
           ZSDS.FIXED = ZSDS.FIXED,
           W.FIXED = W.FIXED)

  if (inherits(fit, "try-error"))
    return(rep(NA,
               length(ncp.start) +
               length(ncp.start) +
               length(zsds.start)))

  c(fit$w.inp,
    fit$ncp,
    fit$zsds)
}


#####################

run_bootstrap_list <- function(INT,
                               ncp.start,
                               zsds.start,
                               Int.Beg,
                               Int.End,
                               B = 1000,
                               NCP.FIXED = TRUE,
                               ZSDS.FIXED = TRUE,
                               W.FIXED = FALSE,
                               cores = round(parallel::detectCores() * .7)) {

  cl <- parallel::makeCluster(cores)

  parallel::clusterExport(
    cl,
    varlist = c("EM_EXT_FOLDED",
                "INT",
                "ncp.start",
                "zsds.start",
                "Int.Beg",
                "Int.End",
                "NCP.FIXED",
                "ZSDS.FIXED",
                "W.FIXED"),
    envir = environment()
  )

  parallel::clusterEvalQ(cl, library(stats))

boot_res <- parallel::parLapply(cl, 1:B, function(i) {

  boot_sample <- sample(INT, length(INT), replace = TRUE)

  fit <- try(
    EM_EXT_FOLDED(
      INT = boot_sample,
      ncp = ncp.start,
      zsds = zsds.start,
      w.inp = rep(1/length(ncp.start), length(ncp.start)),
      Int.Beg = Int.Beg,
      Int.End = Int.End,
      NCP.FIXED = NCP.FIXED,
      ZSDS.FIXED = ZSDS.FIXED,
      W.FIXED = W.FIXED
    ),
    silent = TRUE
  )

  if (inherits(fit, "try-error"))
    return(NULL)

  list(
    w.inp = fit$w.inp.est,
    ncp.est   = fit$ncp.est,
    zsds.est  = fit$zsds.est,
    loglik = fit$loglik
  )
})

parallel::stopCluster(cl)

  # remove failed runs
  boot_res <- boot_res[!sapply(boot_res, is.null)]

  if (length(boot_res) == 0)
    stop("All bootstrap runs failed.")

  return(boot_res)
}



###

run_bootstrap <- function(INT,
                          ncp.start,
                          zsds.start,
                          Int.Beg,
                          Int.End,					
                          B = 1000,
                          NCP.FIXED = TRUE,
                          ZSDS.FIXED = TRUE,
                          W.FIXED = FALSE,
                          cores = round(detectCores() *.7) ) {

  cl <- makeCluster(cores)

  clusterExport(cl,
             varlist = c("EM_EXT_FOLDED",
                          "INT",
                          "bootstrap_one",
                          "ncp.start",
                          "zsds.start",
                          "Int.Beg",
						"Int.End",
                          "NCP.FIXED",
                          "ZSDS.FIXED",
                          "W.FIXED"),
              envir = environment())



  clusterEvalQ(cl, library(stats))

  results <- parLapply(cl, 1:B, function(i) {
    bootstrap_one(
      vals = INT,
      ncp.start = ncp.start,
      zsds.start = zsds.start,
      Int.Beg = Int.Beg, 
      Int.End = Int.End,
      NCP.FIXED = NCP.FIXED,
      ZSDS.FIXED = ZSDS.FIXED,
      W.FIXED = W.FIXED
    )
  })

  stopCluster(cl)

  boot_mat <- do.call(rbind, results)

  colnames(boot_mat) <-
    c(paste0("w", 1:length(ncp.start)),
      paste0("ncp", 1:length(ncp.start)),
      paste0("zsds", 1:length(zsds.start)))

  boot_mat
}


#############

run.new.zcurve = function() {


fit <- EM_EXT_FOLDED(INT = INT,
              ncp = ncp,
              w.inp = rep(1/length(ncp),length(ncp)),
              Int.Beg = Int.Beg
		)
fit

cp.input = list(
	w.inp = fit$w.inp.est,
	ncp.est = fit$ncp.est,
	zsds.est = fit$zsds.est
)

cp.input

cp.res = Compute.Power.General(cp.input,Int.Beg=Int.Beg)

zcurve.res = cp.res

######################

if (boot.iter > 0) {

boot_res <- run_bootstrap_list(INT = INT,
                          ncp.start = ncp,
                          zsds.start = rep(1, length(ncp)),
                          Int.Beg = Int.Beg,
                          Int.End = Int.End,		
                          B = boot.iter,
                          NCP.FIXED = TRUE,
                          ZSDS.FIXED = TRUE,
                          W.FIXED = FALSE)



boot_power <- lapply(boot_res, function(x)
  Compute.Power.Z(x, Int.Beg = Int.Beg)
)


param.mat <- do.call(
  rbind,
  lapply(boot_power, function(x) {
    c(
      EDR = x$EDR,
      ERR = x$ERR,
      setNames(x$w.sig, paste0("w.sig", seq_along(x$w.sig))),
      setNames(x$w.all, paste0("w.all", seq_along(x$w.all))),
      x$loc.pow
    )
  })
)

CI.all <- apply(
  param.mat,
  2,
  quantile,
  probs = c(CI.ALPHA/2, 1-CI.ALPHA/2),
  na.rm = TRUE
)

combine_with_ci <- function(cp.res, CI.all) {

  list(
    EDR = c(
      est   = cp.res$EDR,
      ci.lb = CI.all[paste0(as.character(round(CI.ALPHA/2*100,1)),"%"), "EDR"],
      ci.ub = CI.all[paste0(as.character(round(100-CI.ALPHA/2*100,1)),"%"), "EDR"]
    ),

    ERR = c(
      est   = cp.res$ERR,
      ci.lb = CI.all[paste0(as.character(round(CI.ALPHA/2*100,1)),"%"), "ERR"],
      ci.ub = CI.all[paste0(as.character(round(100-CI.ALPHA/2*100,1)),"%"), "ERR"]
    ),

    w.sig = cbind(
      est   = cp.res$w.sig,
      ci.lb = CI.all[paste0(as.character(round(CI.ALPHA/2*100,1)),"%"), grep("^w.sig", colnames(CI.all))],
      ci.ub = CI.all[paste0(as.character(round(100-CI.ALPHA/2*100,1)),"%"), grep("^w.sig", colnames(CI.all))]
    ),

    w.all = cbind(
      est   = cp.res$w.all,
      ci.lb = CI.all[paste0(as.character(round(CI.ALPHA/2*100,1)),"%"), grep("^w.all", colnames(CI.all))],
      ci.ub = CI.all[paste0(as.character(round(100-CI.ALPHA/2*100,1)),"%"), grep("^w.all", colnames(CI.all))]
    ),

    loc.pow = cbind(
      est   = cp.res$loc.pow,
      ci.lb = CI.all[paste0(as.character(round(CI.ALPHA/2*100,1)),"%"), grep("^lp\\.", colnames(CI.all))],
      ci.ub = CI.all[paste0(as.character(round(100-CI.ALPHA/2*100,1)),"%"), grep("^lp\\.", colnames(CI.all))]
    )
  )
}

zcurve.res <- combine_with_ci(cp.res,CI.all)

}

return(zcurve.res)

} # End of New Zcurve




##########################################

if (1 == 2) {

source(zcurve3)

val.input = val.input[!is.na(val.input)]
INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
summary(INT)

boot.iter = 500

zcurve.res = run.new.zcurve()

zcurve.res.3d = round_list(zcurve.res)
zcurve.res.3d

zcurve.res.3d$w.all
palliation_zcurve3_default$w.all

zcurve.res.3d$EDR
palliation_zcurve3_default$EDR

zcurve.res.3d$ERR
palliation_zcurve3_default$ERR

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


###########################################################


### FUN.2 - HETEROGENEITY TEST

run.heterogeneity.test = function(val.input,boot.run = 500, 
	fit.ci = c(.01,.025,.05,.10,.17,.20,.50,.80,.83,.90,.95,.975,.99)
	) {


res.boot = c()

boot.i = 1
for (boot.i in 1:boot.run ) {
	print(boot.i)

	zboot = sample(val.input,replace=TRUE)
	summary(zboot)

	ncp = 2
	zsds = 1
	components = length(ncp)
	para.est.EXT = extended.curve(zboot,ncp,zsds);para.est.EXT
	fit.hom = para.est.EXT[(3*components+1)];fit.hom

	NCP.FIXED = FALSE
	ncp = c(0,3)
	zsds = c(1,1)
	components = length(ncp)
	para.est.EXT = extended.curve(zboot,ncp,zsds);para.est.EXT
	fit.het.1 = para.est.EXT[length(para.est.EXT)];fit.het.1

	ncp = 2
	zsds = 5
	components = length(ncp)
	para.est.EXT = extended.curve(zboot,ncp,zsds);para.est.EXT
	fit.het.2 = para.est.EXT[4];fit.het.2
	sd.het.2 = para.est.EXT[3];sd.het.2

	delta.fit.1 = fit.hom - fit.het.1
	delta.fit.2 = fit.hom - fit.het.2
	delta.fit.3 = fit.het.2 - fit.het.1

	res.boot = rbind(res.boot,c(fit.hom,fit.het.1,fit.het.2,sd.het.2,
		delta.fit.1,delta.fit.2,delta.fit.3))

} # EOF boot iterations

fit.hom.ci = quantile(res.boot[,1],fit.ci);fit.hom.ci
fit.het.1.ci = quantile(res.boot[,2],fit.ci);fit.het.1.ci
fit.het.2.ci = quantile(res.boot[,3],fit.ci);fit.het.2.ci

sd.het.2.ci = quantile(res.boot[,4],fit.ci);sd.het.2.ci

delta.fit.1.ci = quantile(res.boot[,5],fit.ci);delta.fit.1.ci
delta.fit.2.ci = quantile(res.boot[,6],fit.ci);delta.fit.2.ci
delta.fit.3.ci = quantile(res.boot[,7],fit.ci);delta.fit.3.ci

return.res = cbind(
	fit.hom.ci,fit.het.1.ci,fit.het.2.ci,
	sd.het.2.ci,
	delta.fit.1.ci,delta.fit.2.ci,delta.fit.3.ci)

dim(return.res)

colnames(return.res) = c(
	"fit.hom","fit.het.1","fit.het.2",
	"sd.het.2",	
	"delta.fit.2-1","delta.fit.3-1","delta.fit.2-3"
)

return(return.res)

} # EOF run.heterogeneity.test



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
p.bias.binomial
#print(p.bias.binomial)


### for all
#prob = prob1 / (prob1 + prob2 + prob3)
#just.sig.k/k
#p.bias.binomial = 1 - pbinom(just.sig.k - 1, k, prob = prob.adj);p.bias.binomial


### chi2 approximation
#O = just.sig.k
#E = sig.k*prob;E
#chi2.val <- (O - E)^2 / (E * (1 - prob));chi2.val
#p.bias.chi2 <- 1 - pchisq(chi2.val, df = 1)
#print(p.bias.chi2)

bias.res = c(just.sig.k/sig.k,prob.adj,p.bias.binomial)
#print(bias.res)

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

EXT.boot = function(ZSDS.FIXED = TRUE)	{

	#boot.iter = 100

	components = length(ncp)

	INT = val.input[val.input >= Int.Beg & val.input <= Int.End]

	print("Running parallel bootstrap")

	if (components > 4) {
		print("More than 4 components. Really intended?")
	} else {

	ncores <- floor(parallel::detectCores() * 0.7)
	cl <- parallel::makeCluster(ncores)
	on.exit(parallel::stopCluster(cl), add = TRUE)


	# Make RNG reproducible across workers
	parallel::clusterSetRNGStream(cl, 123)

	# Export everything the workers need
	parallel::clusterExport(cl, varlist = c(
	 "INT",
	 "Get.Densities",
	  "extended.curve","slope",
	  "Augment.Regression","Augment","Augment.Factor","bw.aug","bkde",
	  "W.FIXED","NCP.FIXED","ZSDS.FIXED","TESTING","CURVE.TYPE",
	  "startval","components","ncp","zsds","crit","Int.Beg","Int.End","bw.est","alpha",
	  "Plot.Fitting"
	), envir = environment())

	# If you use packages inside the iteration, load them on workers
	parallel::clusterEvalQ(cl, { NULL })  # e.g., library(stats)

	run.time = system.time({

	results_cp <- parallel::parLapply(cl, 1:boot.iter, function(boot) {

	  val.sample <- sample(INT, size = length(INT), replace = TRUE)

	  para.est.EXT <- extended.curve(vals=val.sample, startval, ncp, zsds)

	list(
	  w.inp    = para.est.EXT$w.inp,
	  ncp.est  = para.est.EXT$ncp.est,
	  zsds.est = para.est.EXT$zsds.est,
	  fit  = para.est.EXT$fit
	)


	}) # End of System Time


	})  # End of Parallel Bootstrap

	} # End of "too many components"


	print(run.time)

	print("End of bootstrap")


	boot_power <- lapply(results_cp, function(cp.input)
  		Compute.Power.Free(cp.input,Int.Beg = 1.96)
	)

	boot.res1 <- do.call(rbind, boot_power)

	CIs = c()
	CIs = rbind(CIs,quantile(boot.res1[,1],c(CI.ALPHA/2,1-CI.ALPHA/2)) )
	CIs[1,1] = CIs[1,1]-ERR.CI.adjust
	CIs[1,2] = CIs[1,2]+ERR.CI.adjust
	if (CIs[1,1] < alpha/2) CIs[1,1] = alpha/2
	if (CIs[1,2] > 1) CIs[1,2] = 1

	CIs = rbind(CIs,quantile(boot.res1[,4],c(CI.ALPHA/2,1-CI.ALPHA/2)) )
	CIs[2,1] = CIs[2,1]-EDR.CI.adjust
	CIs[2,2] = CIs[2,2]+EDR.CI.adjust
	if (CIs[2,1] < alpha) CIs[2,1] = alpha
	if (CIs[2,2] > 1) CIs[2,2] = 1

	FDR = (1/CIs[1,] - 1)*(alpha/(1-alpha))
	CIs = rbind(CIs,FDR[2:1])

	CIs

	# extract bootstrap draws into matrices (rows = boot iterations)
	w.all.mat   <- do.call(rbind, lapply(results_cp, `[[`, "w.inp"))
	NCP.mat <- do.call(rbind, lapply(results_cp, `[[`, "ncp.est"))
	ZSD.mat <- do.call(rbind, lapply(results_cp, `[[`, "zsds.est"))
	fit.mat <- do.call(rbind, lapply(results_cp, `[[`, "fit"))


	# CIs by component
	CI.w.all <- apply(w.all.mat,   2, quantile, probs = c(CI.ALPHA/2, 1 - CI.ALPHA/2), na.rm = TRUE)
	CI.NCP <- apply(NCP.mat, 2, quantile, probs = c(CI.ALPHA/2, 1 - CI.ALPHA/2), na.rm = TRUE)
	CI.ZSD <- apply(ZSD.mat, 2, quantile, probs = c(CI.ALPHA/2, 1 - CI.ALPHA/2), na.rm = TRUE)
	CI.fit <- apply(fit.mat, 2, quantile, probs = c(CI.ALPHA/2, 1 - CI.ALPHA/2), na.rm = TRUE)

	# assemble the parameter CI block in the order you want
	CIs <- rbind(CIs, t(CI.NCP), t(CI.ZSD), t(CI.w.all), t(CI.fit))

	CIs

	rownames(CIs) <- c(
      "EDR","ERR","FDR",
	  rep("NCP", ncol(NCP.mat)),
	  rep("ZSD", ncol(ZSD.mat)),
	  rep("WALL", ncol(w.all.mat)),
	  "fit"	
	)


	CIs = CIs[c(1,2,3:nrow(CIs)),]

	round(CIs,2)

	return(CIs)


} # End function EXT.boot

#####################################
#####################################
#####################################

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
	names(loc.power) = paste0("LP",seq(1,length(loc.power)))
	int = seq(x.lim.min,x.lim.max-int.loc,int.loc)+int.loc/2

	old_mar <- par("mar")
	par(mar = old_mar + c(1, 0, 0, 0))  # increase bottom margin by 2 lines

	new.y = 0

#for (i in 1:length(int)) text(int[i],new.y,
#	paste0(format(round(loc.power[i]*100),nsmall=0),"%"),srt=90,pos=1
#                              ,cex=letter.size)
#par(xpd = FALSE)

lab = c()
for (i in 1:length(int)) lab[i] = paste0("   ",format(round(loc.power[i]*100),nsmall=0),"%")
lab

x_pos  <- seq(x.lim.min, x.lim.max, length.out = length(lab))
x_pos

#axis(1)

mtext(lab, side=1, line=.8, at=x_pos, cex=1.0, las=2)

par(mar = old_mar) 

}


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

	######

	#print("Check write CI")
	#print(Write.CI)

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
	#print("Check ERR.low")
	#print(ERR.low)
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

	# Witout CI
	#print("Writing without CI")

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

#x.start = x.lim.min; x.end = x.lim.max;Ltype = 3;Lwidth = 4;cola = col.hist

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
	k.sim = 500000
	sim.val = c()
	for (i in 1:components) {
		if (CURVE.TYPE == "z") { 
			sim.val = c(sim.val,rnorm(w.inp[i]*k.sim,ncp[i],zsds[i]))
		} else {
			sim.val = c(sim.val,rt(w.inp[i]*k.sim,df,ncp[i]))
		}
	}
	sim.val = abs(sim.val)
	d.sim = Get.Densities(sim.val[sim.val >= x.lim.min & sim.val < x.lim.max],
		bw=bw.draw,d.x.min=x.lim.min,d.x.max=x.lim.max,Augment=Augment)
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
	#print("Not Augmenting")
	Z.INT.USE = INT[INT > d.x.min & INT <= max.z]
}

#print(summary(Z.INT.USE))

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

#print(summary(D))

#print("End Augment check")

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

#print("Create Dens")
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



################################################################
### USE EXTENDED CURVE (Est.Method = "EXT"  (slower than OF)
################################################################

extended.curve = function(vals,startval=NULL,ncp=ncp,zsds=zsds) {

### this is the actual fitting function
extended.curve.fitting = function(para,RetEst=FALSE,Fixed.Null=FALSE)    {

### get the weights
weights = para[1:components]
means = para[(components+1):(2*components)]
sds = para[(components*2+1):(3*components)]

i = 1
j = 1
k = 1
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
summary(Dens)


### compute the new estimated density distribution
E.D.Y = colSums(Dens*weights)

### compare to observed density distribution
rmse = sqrt(mean((E.D.Y-O.D.Y)^2));rmse


### return either fit if continue or estimates if finished
value = rmse
if(RetEst) value = est


### showing the fitting of the function in a plot
if(Plot.Fitting) {

	rval = runif(1)
	if (rval > .4) {

	tit = ""
	xyy = cbind(D.X,O.D.Y,E.D.Y)
	plot(xyy[,1],xyy[,2],type='l',
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax),
		main=tit,xlab='Z')
	lines(xyy[,1],xyy[,3],lty=1,col="red1",
		xlim=c(x.lim.min,x.lim.max),ylim=c(ymin,ymax))
	}

}


### return value to optimization function
return(value)

} 

#vals = val.sample
densy = Get.Densities(vals,bw=bw.est,d.x.min=Int.Beg,d.x.max=Int.End,Augment=Augment)

components = length(ncp);components

D.X = densy[,1]
O.D.Y = densy[,2]

#plot(D.X,O.D.Y)

n.bars = length(D.X)
n.bars

bar.width = D.X[2] - D.X[1]
bar.width

#startval = NULL

if (is.null(startval)) {

if (W.FIXED) {
	startval = w.fix
} else {
	startval = rep(1/components, components)
}

startval = c(startval,ncp)

startval = c(startval,zsds)

} # End of startval

###

if (W.FIXED) {
	lowlim = w.fix
	highlim = w.fix
} else {
	lowlim = rep(0,components)
	highlim = rep(1,components)
}

if (NCP.FIXED) {
	lowlim = c(lowlim,ncp)
	highlim = c(highlim,ncp)
} else {
	lowlim = c(lowlim,rep(0,components))
	highlim = c(highlim,rep(Int.End,components))
}

if (ZSDS.FIXED) { 
	lowlim = c(lowlim,zsds)
	highlim = c(highlim,zsds)
} else {
	lowlim = c(lowlim,rep(1,components))
	highlim = c(highlim,zsds)
}

if(components == 1) lowlim[1] = 1

#print(startval);print(lowlim);print(highlim)

#TESTING = TRUE
if (TESTING == TRUE) Plot.Fitting = TRUE

auto = nlminb(startval,extended.curve.fitting,lower=lowlim,upper=highlim,control=list(eval.max=1000))

para = auto$par

w.inp = para[1:components]
w.inp = w.inp / sum(w.inp)
ncp.est = para[(components+1):(2*components)]
zsds.est = para[(components*2+1):(3*components)]

fit = auto$objective

res = list(
    ncp.est = ncp.est,
	zsds.est = zsds.est,
	w.inp = w.inp,
    fit = fit
)

res

return(res)

} 

######################################################
### End of Extended Zcurve
#######################################################


#######################################################
### Begin Old Fashioned Zcurve (Est.Method = "OF" 
#######################################################


#########################################################################
### This Function Computes Power from Weights and Non-Centrality Parameters
#########################################################################

Compute.Power.Z = function(cp.input,Int.Beg=crit) {

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
est.ncp = cp.input$ncp.est
est.zsds = cp.input$zsds.est

components = length(est.ncp)

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


res = list(
   EDR = EDR,
   ERR = ERR,
   w.sig = w.sig,
   w.all = w.all,
   loc.pow = local.power
)

#print(res)

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

print(res)

return(res)

} # EOF Compute.Power.T


#####################################
#####################################
#####################################




###############################################################

Compute.Power.General = function(cp.input,Int.Beg=1.96) {

ext.all = length(val.input[val.input > Int.End]) / 
	length(val.input)
ext.all

ext.sig = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > crit])
ext.sig

ext.inp = length(val.input[val.input > Int.End]) / 
	length(val.input[val.input > Int.Beg])
ext.inp


components = length(cp.input$ncp.est)

w.inp = cp.input$w.inp
est.ncp = cp.input$ncp.est
est.zsds = cp.input$zsds.est
est.zsds[est.zsds < 1] = 1  
est.tau = sqrt(est.zsds^2 - 1);est.tau

###

zx.bw = .01
zx = seq(0,Int.End,zx.bw)
pow.zx.dir = pnorm(abs(zx),1.96) 
pow.zx.sign.error = + pnorm(-1.96,abs(zx))
pow.zx = pow.zx.dir + pow.zx.sign.error

i = 1
est.wd.all = c()
for (i in 1:components) {
	if (est.tau[i] < .01) {
		wd = rep(0,length(zx))
		wd[which(round(zx,2) == round(est.ncp[i],2))] = w.inp[i]
	} else {
		wd = dnorm(zx,est.ncp[i],est.tau[i])*w.inp[i]
	}
	wd = wd / pow.zx
	wd = wd / sum(wd)
	est.wd.all = rbind(est.wd.all,wd)
}

dim(est.wd.all)

est.wd.all = colMeans(est.wd.all)
sum(est.wd.all)



est.wd.all.ext = c(est.wd.all*(1-ext.all),ext.all)
pow.zx.dir.ext = c(pow.zx.dir,1)

#cbind(est.wd.all.ext,pow.zx.dir.ext)

pow.zx.ext = c(pow.zx,1)
EDR = sum(pow.zx.ext*est.wd.all.ext);EDR

est.wd.sig.ext = est.wd.all.ext*pow.zx.ext
est.wd.sig.ext = est.wd.sig.ext/sum(est.wd.sig.ext)

ERR = sum(pow.zx.dir.ext*est.wd.sig.ext)/sum(est.wd.sig.ext);ERR

if (ERR > 1) ERR = 1
if (EDR > 1) EDR = 1

if (ERR < alpha/2) ERR = alpha/2
if (EDR < alpha) EDR = alpha 

ERR
EDR

EDR.pos = NA
EDR.neg = NA

ERR.pos = NA
ERR.neg = NA

res.est = c(EDR,EDR.pos,EDR.neg,ERR,ERR.pos,ERR.neg)
res.est
res = c(res.est,rep(NA,components),rep(NA,components))
names(res) = c("EDR","EDR.pos","EDR.neg","ERR","ERR.pos","ERR.neg",
paste0("w.all.",ncp[1:components]),paste0("w.sig.",ncp[1:components]))
res

if (int.loc > 0) {


# ---- local true-power by ncp segments (0..Int.End) ----
bar.width <- 0.01
ncp.grid <- seq(0, Int.End, by = bar.width)

# build density over TRUE ncp (sampling variance removed): mixture of normals on ncp axis
dd <- 0
for (k in 1:components) {
  if (est.tau[k] < 0.1) {
    # nearly discrete component
    dd_k <- rep(0, length(ncp.grid))
    dd_k[which(round(ncp.grid, 2) == round(est.ncp[k], 2))] <- 1
  } else {
    dd_k <- dnorm(ncp.grid, mean = est.ncp[k], sd = est.tau[k])
  }
  dd <- dd + est.cw.all[k] * dd_k
}
dd <- dd / sum(dd)

power <- pnorm(ncp.grid,crit) + pnorm(-crit, ncp.grid) 

# segment edges
int <- seq(0, Int.End, by = int.loc)

local.power <- c()
for (i in 1:(length(int)-1)) {
  segment <- (ncp.grid >= int[i]) & (ncp.grid < int[i+1])
  denom <- sum(dd[segment])
  local.power <- c(local.power, if (denom == 0) NA else sum(power[segment] * dd[segment]) / denom)
}
names(local.power) <- paste0("lp.", seq_along(local.power))

# append to output
res <- c(res, local.power)

#print("Compute Extended Power End")

} ### end of local power


### to be past back to the main program from this function
return(res)

} ### EOF Compute.Power for SDG1


#####################################
#####################################
#####################################



#####################################
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
### BBB ZingStart #START #Begin of Main Program 
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
print(crit)

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


### End of Preparation
###############################################
###############################################
###############################################


### OLD FASHIONED Z-CURVE 

#zsds = rep(1,components)

if (Est.Method %in% c("OF","EM","CLU")) {

#summary(val.input)
para.est.OF = old.fashioned(val.input)

fit = para.est.OF[components+1];fit

cp.input  = list(
	w.inp = para.est.OF[1:components],
	ncp.est = ncp,
	zsds.est = zsds
)


if (CURVE.TYPE == "z") {
  cp.res = Compute.Power.Z(cp.input,Int.Beg=Int.Beg)
} else {
  cp.res = Compute.Power.T(cp.input,Int.Beg=Int.Beg)
}

EDR = cp.res$EDR
ERR = cp.res$ERR

w.all = cp.res$w.all

w.sig = cp.res$w.sig

loc.power = cp.res$loc.pow

#print("Finish Old Fashioned")

} # EOF if Est.Method OF
 
#########################################
### EXTENDED Z-CURVE

if (Est.Method == "EXT") {

	components = length(ncp)

	INT = val.input[val.input >= Int.Beg & val.input <= Int.End]
	startval = NULL

	para.est.EXT = extended.curve(INT,startval,ncp,zsds);para.est.EXT

	fit = para.est.EXT$fit

	cp.input = list(
		w.inp = para.est.EXT$w.inp,
		ncp.est = para.est.EXT$ncp.est,
		zsds.est = para.est.EXT$zsds.est
	)

	startval = c(cp.input$w.inp,cp.input$ncp.est,cp.input$zsds.est)
	startval

	if (max(cp.input$zsds.est) < 1.05) {
		if (CURVE.TYPE == "z") {
			cp.res = Compute.Power.Z(cp.input,Int.Beg=Int.Beg)
		} else {
			cp.res = Compute.Power.T(cp.input,Int.Beg=Int.Beg)
		}
	} else {
		print("SD > 1")
		cp.res = Compute.Power.SDG1(cp.input,Int.Beg=Int.Beg)
	}			

	round(cp.res,3)	

	EDR = cp.res[1];EDR
	ERR = cp.res[4];ERR

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	if (components == 1) w.all = 1

	if (Int.Beg == 0) w.all = w.inp

	w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
	if (components == 1) w.sig = 1
	round(w.sig,3)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	loc.power

#	print("Extended Version Completed")

} # EOF Extended Version 


##########################################################

 
### if requested, run slower EM
### Est.Method = "EM"
if(Est.Method == "EM"& boot.iter == 0) {

	components = length(ncp)	
	#summary(val.input)
	#print("Fitting EM")

	#Int.Beg = 2.4; alpha = .005

	res.em = run.zcurve(val.input,Est.Method="EM",boot.iter=boot.iter,alpha=alpha,
		Int.Beg=Int.Beg,Int.End=Int.End,parallel=parallel)
	summary(res.em)     

	fit = res.em$fit$Q
	fit

	para.est.EM = res.em$fit$weights
	para.est.EM

	para.est.EM = c(para.est.EM,ncp,zsds)

	cp.res = Compute.Power.Z(para.est.EM,Int.Beg=Int.Beg)
	cp.res

	if (1 == 2) {

		para.est.OF = old.fashioned(val.input)
		fit = para.est.OF[components+1];fit
		para.est.OF = c(para.est.OF[1:components],ncp,zsds)
		cp.res.of = Compute.Power.Z(para.est.OF,Int.Beg=Int.Beg)

		summary(res.em)$coefficients[2:1]
		cp.res[c(1,4)]
		cp.res.of[c(1,4)]

	}


	EDR = cp.res[1];EDR
	ERR = cp.res[4];ERR

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	w.all
	sum(w.all)

	w.sig = cp.res[which(substring(names(cp.res),1,5) == "w.sig")]
	round(w.sig,3)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(loc.power,3)


	#print("Finished EM")
} 



### use only for testing, replaced by OF
if(Est.Method == "density" & boot.iter == 0) {

	print(paste("Est.Method = ",Est.Method))

	components = length(ncp)	
	#summary(val.input)
	#print("Fitting EM")
	z.res = run.zcurve(val.input,Est.Method="density",boot.iter=boot.iter,
		Int.Beg=Int.Beg,Int.End=Int.End,parallel=parallel)
	summary(z.res)     
	para.est.EM = c(summary(z.res, type="parameters")$coefficients)
	para.est.EM = para.est.EM[c((components+1):(2*components),1:components)]
	para.est.EM = c(para.est.EM,ncp)
	para.est.EM
	
	fit = z.res$fit$objective

	w.inp =	 summary(z.res, type="parameters")$coefficients[(components+1):(2*components)]
	round(w.inp,3)

	#w.all = w.sig/(pow.dir+sign.error)
	#w.all = w.all/sum(w.all)
	#round(w.all,3)

	cp.res = Compute.Power.Z(Int.Beg = Int.Beg,c(w.inp,ncp,zsds))
	round(cp.res,3)

	w.all = cp.res[which(substring(names(cp.res),1,5) == "w.all")]
	round(w.all,3)
	sum(w.all)

	loc.power = cp.res[which(substring(names(cp.res),1,2) == "lp")]
	round(loc.power,3)

	EDR = z.res$coefficients[2];EDR
	ERR = z.res$coefficients[1];ERR

	print("Finished density")

} #EOF density


### ADD RESULTS

FDR = round((1/EDR - 1)*(alpha/(1-alpha)),2);
names(FDR) = "FDR"
FDR

bias = c(NA,NA,NA)
if(TEST4BIAS) { 
	bias = test.bias(w.all) 
}
names(bias) = c("OBS.JS","EXP.JS","EJS.p")

### END OF ADD RESULTS


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

	cp.res = Compute.Power.Z(Int.Beg = Int.Beg,c(w.inp,ncp,zsds))
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



##########################################
### This Code is Used to Create Graphic
##########################################

print("Show Histogram")

if (Show.Histogram & sum(extreme,na.rm=TRUE) < .95) { 

	Draw.Histogram(w.all,cola=col.hist)

	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)

	if (Show.Curve.All & max(zsds) < 1.05) {

		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)

		}

	if (Show.Curve.All & max(zsds) > 1.05) {
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype = 3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (length(loc.power > 0)) Write.Local.Power(loc.power)

	} # End of Show.Histogram


########################################################################
########################################################################
########################################################################

res.het = NA
if (TEST4HETEROGENEITY > 0) {
	boot.run = TEST4HETEROGENEITY
	if (!is.numeric(TEST4HETEROGENEITY)) print("!!! Need to provide number of iterations!!!")
	else res.het = run.heterogeneity.test(val.input,boot.run=boot.run,fit.ci)
	}



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

	#print("MAKE RESTEXT 4 CLU")

	results = rbind(c(ODR,ODR.low,ODR.high),EDR,ERR,FDR)
	rownames(results) = c("ODR","EDR","ERR","FDR")
	round(results,3)

}

if (boot.iter > 0 & Show.Histogram) {

##############################
	print(results$EDR) 
#################################################

	Draw.Histogram(w.all,cola=col.hist,Write.CI = TRUE)

	if (Show.KD) Draw.KD(val.input,w.all,cola=col.kd)

	if (Show.Curve.All & Est.Method != "DF" & max(zsds) < 1.05) {
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All(w=w.all,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (Show.Curve.All & Est.Method != "DF" & max(zsds) > 1.05) {
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype = 3,Lwidth = 4,x.start=x.lim.min,x.end=x.lim.max)
		Draw.Curve.All.SDG1(w=w.all,ncp=ncp,zsds=zsds,cola=col.curve,
			Ltype=1,Lwidth = 4,x.start=Int.Beg,x.end=Int.End)
		}

	if (Show.Curve.Sig) {
		Draw.Curve.Sig(val.input,ncp=ncp,zsds=zsds,w.sig,cola=col.curve,Ltype=3,x.start=x.lim.min)
		Draw.Curve.Sig(val.input,ncp=ncp,zsds=zsds,w.sig,cola=col.curve,Ltype=1,x.start=Int.Beg)
		}
		#	par(family = fam[1])	
		if (length(loc.power > 0)) Write.Local.Power(loc.power)

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

