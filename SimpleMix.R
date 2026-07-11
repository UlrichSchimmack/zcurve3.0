

### New simpler simulation 
print("SimpleMix")
print("2006-07-09")



select_positive_only <- function(
  sim,                    # output of sim_mixture (must return `ref`, not just dat.fit)
  k_fit   = 1000,
  crit    = 0
) {

  ref <- sim$ref          # <-- need sim_mixture to return the full ref pool

  # ---- selection: positive AND significant ------------------------------
  crit <- ref$crit
  keep <- ref$obs_d > crit    # positive-significant only (t > +crit)
  sel  <- ref[keep, ]

  if (nrow(sel) < k_fit)
    stop(sprintf("Only %d positive-significant in pool, need %d; raise k_ref",
                 nrow(sel), k_fit))

  # ---- true values for the SELECTED population --------------------------
  # (means among the selected set, from the large pool = near-exact)
  df       <- sel$Ns - 2
  crit.sel <- sel$crit
  pow.dir  <- 1 - pt( crit.sel, df = df, ncp = sel$ncp)   # directional power
  # replication = significant again, same (positive) direction:
  true_ERR_sel  <- mean(pow.dir)                          # among selected
  true_mean_sel <- mean(sel$pop_d)                        # mean true effect, selected
  true_FDR_sel  <- mean(sel$h == 0)                       # fraction null among selected

  # EDR is a PRE-selection quantity: mean two-sided power over the WHOLE pool,
  # NOT the selected set. Compute from ref, not sel.
  df.ref  <- ref$Ns - 2
  pow.two.ref <- (1 - pt(ref$crit, df=df.ref, ncp=ref$ncp)) +
                       pt(-ref$crit, df=df.ref, ncp=ref$ncp)
  true_EDR <- mean(pow.two.ref)

  # ---- draw the small fitting sample from the SELECTED set --------------
  fit.idx <- sample(seq_len(nrow(sel)), size = k_fit, replace = FALSE)
  dat.fit <- sel[fit.idx, ]

  list(
    dat           = dat.fit,        # positive-significant only -> feed z-curve
    true_EDR      = true_EDR,       # pre-selection (from full pool)
    true_ERR_sel  = true_ERR_sel,   # replication rate of selected studies
    true_mean_sel = true_mean_sel,  # mean effect among selected (z-curve es target)
    true_FDR_sel  = true_FDR_sel    # FDR among selected
  )
}


#######################

SimpleMix <- function(
  k         = 1000,       # small sample actually fed to z-curve
  p_h0      = 0.40,
  d         = 0.40,
  SD        = 0.20,
  N         = 100,
  N_low     = NULL,
  N_high    = NULL,
  min_d     = -1,
  max_d     = 1,
  min_z     = -8,
  max_z     = 8,
  pop_dist  = "normal",
  pop_dist_beta_shape1 = 1,
  pop_dist_beta_shape2 = 10,
  alpha     = 0.05,
  only.pos  = FALSE,  
  sel.crit  = 0
) {

  #helper: get beta distribution of effect sizes
  r_true_beta_shape <- function(k, mean_d, sd_d, shape1 = 2, shape2 = 5,
                              exact_moments = TRUE) {
  
    if (shape1 <= 0 || shape2 <= 0) {
      stop("shape1 and shape2 must be positive.")
    }
  
    # Draw beta values on [0, 1]
    x <- rbeta(k, shape1 = shape1, shape2 = shape2)
  
    if (exact_moments) {
      # Force exact realized mean and SD in this simulated population
      pop_d <- as.numeric(scale(x)) * sd_d + mean_d
    } else {
      # Match mean and SD in expectation using beta moments
      beta_mean <- shape1 / (shape1 + shape2)
      beta_sd <- sqrt(shape1 * shape2 /
                      ((shape1 + shape2)^2 * (shape1 + shape2 + 1)))
    
      pop_d <- ((x - beta_mean) / beta_sd) * sd_d + mean_d
    }
    attr(pop_d, "shape1") <- shape1
    attr(pop_d, "shape2") <- shape2
    attr(pop_d, "target_mean") <- mean_d
    attr(pop_d, "target_sd") <- sd_d
    attr(pop_d, "realized_mean") <- mean(pop_d)
    attr(pop_d, "realized_sd") <- sd(pop_d)
    attr(pop_d, "realized_min") <- min(pop_d)
    attr(pop_d, "realized_max") <- max(pop_d)
  
    pop_d
  }

  ########
 
  # ---- helper: generate one pool of size k_gen --------------------------
  gen_pool <- function(k_gen) {
    if (is.null(N_low) || is.null(N_high)) {
      Ns <- rep(N, k_gen)
    } else {
      Ns <- round(runif(k_gen, N_low, N_high))
    }
    se <- 2 / sqrt(Ns)

    n_h0 <- round(k_gen * p_h0)
    n_h1 <- k_gen - n_h0
    h    <- c(rep(0, n_h0), rep(1, n_h1))

    ### Effect Sizes
    pop_d <- numeric(k_gen)
    pop_d[h == 0] <- 0
    if (pop_dist == "normal") {
      pop_d[h == 1] <- rnorm(n_h1, d, SD)
    } else if (pop_dist == "beta") {
      pop_d[h == 1] <- r_true_beta_shape(k = n_h1, mean_d = d, sd_d = SD)
    } else if (pop_dist == "t.dist") {
      shape_df = 5
      t_std <- rt(n_h1, shape_df) / sqrt(shape_df / (shape_df - 2))    # now mean 0, SD 1
      pop_d[h == 1] <- d + SD * t_std
    } else stop("bad pop_dist")

    length(pop_d)
    length(se)
    ncp   <- pop_d / se
    df    <- Ns - 2
    t     <- rt(k_gen, df = df, ncp = ncp)
    obs_d <- t * se
    p     <- 2 * pt(abs(t), df = df, lower.tail = FALSE)
    z     <- qnorm(1 - p/2) * sign(t)

    crit   <- qt(1 - alpha/2, df = df)
    is.sig <- abs(t) > crit

    data.frame(h, Ns, se, pop_d, ncp, t, p, obs_d, z, is.sig, crit)
  }

  # START
  # ---- LARGE reference pool: compute true values from this --------------

  k_gen = k*100
  k_fit = k

  ref <- gen_pool(k_gen)

  sel1 = ref$z > min_z & ref$z < max_z
  sel2 = ref$pop_d > min_d & ref$pop_d < max_d
  table(sel1,sel2)
  sel = sel1 & sel2
  table(sel)
  ref = ref[sel,]


  # --- pre-selection truths from the FULL ref (EDR especially) ---
  df.ref  <- ref$Ns - 2
  pow.two <- (1 - pt(ref$crit, df=df.ref, ncp=ref$ncp)) +
                   pt(-ref$crit, df=df.ref, ncp=ref$ncp)
  pow.dir <- 1 - pt(ref$crit, df=df.ref, ncp=ref$ncp)
  true_EDR  <- mean(pow.two)                       # pre-selection, from full ref
  true_mean <- mean(ref$pop_d)                     # all studies


  # --- SELECT on the large ref pool ---
  if (only.pos) {
    keep <- ref$t > sel.crit          # sel.crit=0 -> positive; pass ref$crit for pos-sig
  } else {
    keep <- rep(TRUE, nrow(ref))
  }
  sel <- ref[keep, ]
  if (nrow(sel) < k_fit)
    stop(sprintf("Only %d selected, need %d; raise k_gen", nrow(sel), k_fit))

  # --- post-selection truths from sel ---
  true_mean_sig <- mean(sel$pop_d)
  true_FDR      <- mean(sel$h == 0)
  # --- pre-selection truths from ref ---
  true_mean_h1  <- mean(ref$pop_d[ref$h == 1])
  true_ERR      <- sum(pow.two * pow.dir) / sum(pow.two)

  # --- draw k_fit from the SELECTED pool (fixed index range) ---
  fit.idx <- sample(seq_len(nrow(sel)), size = k_fit, replace = FALSE)
  dat.fit <- sel[fit.idx, ]
  
  list(
    dat          = dat.fit,
    true_EDR     = true_EDR,
    true_ERR     = true_ERR,
    true_mean    = true_mean,
    true_mean_h1 = true_mean_h1,
    true_mean_sig= true_mean_sig,
    true_FDR     = true_FDR
  )

}



