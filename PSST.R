

#rm(list = ls())

#yi = rnorm(1000,0,.4)
#sei = rnorm(1000,.2,0)
#sig_rule = "two.sided"
#mu = NULL
#tau2 = NULL
#use_abs_mu = TRUE

ess_tests <- function(yi, sei, mu = NULL, tau2 = NULL, odr = NULL,
                      alpha = .05,
                      crit = qnorm(.975),
                      sig_rule = c("positive", "two.sided"),
                      use_abs_mu = TRUE,
                      exact_tess = FALSE) {

  sig_rule <- match.arg(sig_rule)

  #sig_rule <- "two.sided"

  stopifnot(length(yi) == length(sei))
  keep <- is.finite(yi) & is.finite(sei) & sei > 0
  yi  <- yi[keep]
  sei <- sei[keep]

  k <- length(yi)
  if (k == 0) stop("No valid studies.")

  # UWLS / fixed-effect estimate if not supplied
  if (is.null(mu)) {
    w <- 1 / sei^2
    mu <- if(use_abs_mu) sum(w * abs(yi))/sum(w) else sum(w * yi) / sum(w) 
  }

  # DerSimonian-Laird tau^2 if not supplied
  if (is.null(tau2)) {
    w <- 1 / sei^2
    Q <- sum(w * (yi - mu)^2)
    C <- sum(w) - sum(w^2) / sum(w)
    tau2 <- max(0, (Q - (k - 1)) / C)
  }

  # Stanley et al. equation (3)
  Zi <- (crit * sei - mu) / sqrt(sei^2 + tau2)
  

  # Stanley et al. equation (4)
  p_sig_expected_i <- 1 - pnorm(Zi)
  Esig <- sum(p_sig_expected_i)
  pi_exp <- Esig / k

  # Observed number of significant results
  if(is.null(odr)) {
    if (sig_rule == "positive") {
      sig_obs <- yi / sei > crit
    } else {
      sig_obs <- abs(yi / sei) > crit
    }
    SS <- sum(sig_obs)
    p_obs <- SS / k
  } else {
    p_obs = odr
  }

  # Excess statistical significance
  ESS <- (SS - Esig) / k

  # PSST: Stanley et al. equation (5)
  z_psst <- (p_obs - pi_exp) / sqrt(pi_exp * (1 - pi_exp) / k)
  p_psst <- 1 - pnorm(z_psst)

  # TESS: Stanley et al. equation (6)
  z_tess <- (ESS - .05) / sqrt(.05 * (1 - .05) / k)
  p_tess <- 1 - pnorm(z_tess)

  # Optional exact binomial version for TESS.
  # Null: ESS <= .05, so expected significant count under null is Esig + .05*k.
  # This follows the logic of testing whether observed SS exceeds Esig by more than 5%.
  if (exact_tess) {
    p0_tess <- min(1, max(0, pi_exp + .05))
    p_tess_exact <- binom.test(SS, k, p = p0_tess, alternative = "greater")$p.value
  } else {
    p_tess_exact <- NA_real_
  }

  list(
    k = k,
    mu = mu,
    tau2 = tau2,
    tau = sqrt(tau2),
    SS = SS,
    p_obs = p_obs,
    Esig = Esig,
    p_exp = pi_exp,
    ESS = ESS,
    PSST_Z = z_psst,
    PSST_p = p_psst,
    TESS_Z = z_tess,
    TESS_p = p_tess,
    TESS_p_exact = p_tess_exact,
    TESS_sig = p_tess < alpha,
    PSST_sig = p_psst < alpha,
    TESSPSST_sig = (p_tess < alpha) | (p_psst < alpha)
  )
}


TES.het <- function(ess, ses) {

  k <- length(ess)
  
  # UWLS estimate of mean effect
  w <- 1 / ses^2
  uwls <- sum(w * ess) / sum(w)
  
  # DerSimonian-Laird tau^2
  Q <- sum(w * (ess - uwls)^2)
  tau2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  
  # Expected significance: Equation (3)
  z_i <- (1.96 * ses - abs(uwls)) / sqrt(ses^2 + tau2)
  power_i <- 1 - pnorm(z_i)
  
  # Expected and observed proportions
  prop.exp.sig <- mean(power_i)
  prop.obs.sig <- mean(abs(ess / ses) > 1.96)
  
  # Binomial test
  p_tess <- 1 - pbinom(k * prop.obs.sig - 1, k, prop.exp.sig)
  
  out <- c(prop.obs.sig, prop.exp.sig, p_tess)
  names(out) <- c("ODR", "TES.het EDR", "TES.het p-value")

  return(out)
}


