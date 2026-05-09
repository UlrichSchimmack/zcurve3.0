




### Load weightr in main program. 

if (!"weightr" %in% loadedNamespaces()) {
  library(weightr)
}


### Clustered Step-Function Selection Model (Weight-R with parallel bootstrap)

#nboot = 500;steps=c(.5,.05,.025,.005);ncores=16;seed = 2026
#yi = dat$d; vi = dat$mn

run_boot_cluster_weightr <- function(yi, vi, cluster, 
                                      nboot = 500, 
                                      steps = c(.05),
                                      ncores = parallel::detectCores() - 1,
                                      seed = 2026) {
  
  # --- Original estimate ---
  orig <- weightfunct(yi, vi, steps = steps)
  n_par <- length(orig[[2]]$par)  # mean + tau + n_weights
  n_weights <- n_par - 2

  # Extract estimates
  tau_orig  <- sqrt(orig[[2]]$par[1])   # heterogeneity
  est_mean_orig  <- orig[[2]]$par[2]   # adjusted mean
  weights_orig <- orig[[2]]$par[-(1:2)]  # weight parameters

  # Build interval labels from sorted steps
  steps_sorted <- sort(steps)
  boundaries <- c(0, steps_sorted, 1)
  interval_labels <- paste0("p in (", boundaries[1:length(boundaries)-1], 
                            ", ", boundaries[2:length(boundaries)], "]")
  
#  cat("Original adjusted mean:  ", round(est_mean_orig, 4), "\n")
#  cat("Tau (heterogeneity):     ", round(tau_orig, 4), "\n")
#  cat("Weights:\n")
#  cat(sprintf("  Reference (=1):  %s\n", interval_labels[1]))
#  for (j in seq_len(n_weights)) {
#    cat(sprintf("  Weight %d = %.4f:  %s\n", j, weights_orig[j], interval_labels[j+1]))
#  }
#  cat("\n")
  
  # --- Cluster bootstrap ---
  clusters <- unique(cluster)
  k_clust  <- length(clusters)
  
  cat("Clusters:", k_clust, " | Effects:", length(yi), 
      " | Bootstrap iterations:", nboot, "\n")
  
  set.seed(seed)
  boot_samples <- lapply(1:nboot, function(i) {
    sample(clusters, k_clust, replace = TRUE)
  })
  
  expand_clusters <- function(resampled_ids) {
    idx <- unlist(lapply(resampled_ids, function(cid) which(cluster == cid)))
    return(idx)
  }
  
  # Single bootstrap iteration
  boot_one <- function(b) {
    idx <- expand_clusters(boot_samples[[b]])
    yi_b <- yi[idx]
    vi_b <- vi[idx]
  
    res <- tryCatch({
      fit <- weightfunct(yi_b, vi_b, steps = steps)
      c(adj   = fit[[2]]$par[1],
        unadj = fit[[1]]$par[1],
        fit[[2]]$par[-1])          # tau + all weight params
    }, error = function(e) {
      rep(NA, 1 + n_par)
    })

    return(res)
  }

  # Run in parallel
  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = c("yi", "vi", "cluster", "clusters", "k_clust",
                               "boot_samples", "expand_clusters", "steps", "n_par"),
                envir = environment())
  clusterEvalQ(cl, library(weightr))
  
  cat("Running", nboot, "bootstrap iterations on", ncores, "cores...\n")
  t0 <- proc.time()
  
  boot_results <- parLapply(cl, 1:nboot, boot_one)
  
  stopCluster(cl)
  elapsed <- (proc.time() - t0)[3]
  cat("Done in", round(elapsed, 1), "seconds\n\n")
  
  # --- Combine results ---
  boot_mat <- do.call(rbind, boot_results)

  boot_mat
 
  colnames(boot_mat) <- c("adj_tau2", "unadj_tau2", "adj_mean",
                         paste0("weight_", seq_len(n_weights)))
  
  n_fail <- sum(is.na(boot_mat[,1]))
  cat("Convergence failures:", n_fail, "of", nboot, 
      "(", round(100 * n_fail / nboot, 1), "%)\n")
  
  # Drop failures for CIs
  boot_clean <- boot_mat[!is.na(boot_mat[,1]), , drop = FALSE]

  pred_draws <- unlist(lapply(1:nrow(boot_clean), function(b) {
    rnorm(100, 
        mean = boot_clean[b, "adj_mean"], 
        sd   = sqrt(boot_clean[b, "adj_tau2"]))
   }))
  pred.interval = quantile(pred_draws, c(.025, .975))
 
  # Percentile CIs and bootstrap medians
  ci_mean_adj   <- quantile(boot_clean[,"adj_mean"], c(.025, .975))
  ci_tau_adj   <- quantile(sqrt(boot_clean[,"adj_tau2"]), c(.025, .975))
  
  med_adj   <- median(boot_clean[,"adj_mean"])
  med_tau   <- median(sqrt(boot_clean[,"adj_tau2"]))
  
  cat("\n--- Results ---\n")
  cat(sprintf("                   Original  Median   [95%% CI]\n"))
  cat(sprintf("Adjusted mean:     %.4f    %.4f   [%.4f, %.4f]\n", 
              est_mean_orig, med_adj, ci_mean_adj[1], ci_mean_adj[2]))
  cat(sprintf("Tau:               %.4f    %.4f   [%.4f, %.4f]\n", 
              tau_orig, med_tau, ci_tau_adj[1], ci_tau_adj[2]))

  for (j in seq_len(n_weights)) {
    wname <- paste0("weight_", j)
    ci_w <- quantile(boot_clean[, wname], c(.025, .975))
    med_w <- median(boot_clean[, wname])
    cat(sprintf("Weight %d (%s):\n                   %.4f    %.4f   [%.4f, %.4f]\n",
                j, interval_labels[j+1], weights_orig[j], med_w, ci_w[1], ci_w[2]))
  }

  cat("\n")
 
  cat(sprintf("Prediction Interval ranges from %.2f to %.2f\n",
      pred.interval[1], pred.interval[2]))

  
  # Return everything
  invisible(list(
    original = orig,
    est_orig = est_mean_orig,
    tau_orig = tau_orig,
    weights_orig = weights_orig,
    ci_mean_adj = ci_mean_adj,
    ci_tau_adj = ci_tau_adj,
    med_adj = med_adj,
    med_tau = med_tau,
    boot_matrix = boot_mat,
    n_fail = n_fail,
    interval_labels = interval_labels,
    pred.interval = pred.interval
  ))
}
