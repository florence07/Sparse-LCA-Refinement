clip01 <- function(x, eps = 1e-8) {
  pmin(pmax(x, eps), 1 - eps)
}

loglik_grouped <- function(s, w, z, theta) {
  sum(s * log(theta[z]) + (w - s) * log1p(-theta[z]))
}

compute_theta <- function(s, w, z, k, eps = 1e-8) {
  theta <- numeric(k)
  for (g in seq_len(k)) {
    idx <- which(z == g)
    sw <- sum(w[idx])
    theta[g] <- if (sw > 0) sum(s[idx]) / sw else 0.5
  }
  clip01(theta, eps)
}

repair_empty_groups <- function(z, p, w, k) {
  present <- sort(unique(z))
  missing <- setdiff(seq_len(k), present)
  if (!length(missing)) return(z)

  for (g_new in missing) {
    group_size <- tabulate(z, nbins = k)
    group_mass <- sapply(seq_len(k), function(g) sum(w[z == g]))

    donor <- which.max(ifelse(group_size > 1, group_mass, -Inf))
    if (!is.finite(group_mass[donor])) {
      donor <- which.max(group_size)
    }

    donor_idx <- which(z == donor)
    donor_center <- if (sum(w[donor_idx]) > 0) {
      weighted.mean(p[donor_idx], w[donor_idx])
    } else {
      mean(p[donor_idx])
    }

    take <- donor_idx[which.max(abs(p[donor_idx] - donor_center))]
    z[take] <- g_new
  }

  z
}

init_partition_grid <- function(p, w, k) {
  K <- length(p)
  
  if (length(w) != K) stop("length(w) must equal length(p)")
  if (k < 1 || k > K) stop("k must be between 1 and K")
  
  # sort classes by empirical probabilities
  ord <- order(p)
  
  # special case: one group
  if (k == 1) {
    z <- rep(1L, K)
    return(list(z))
  }
  
  # choose k-1 cut points among K-1 adjacent gaps
  cut_list <- combn(K - 1, k - 1, simplify = FALSE)
  
  z_list <- vector("list", length(cut_list))
  
  for (m in seq_along(cut_list)) {
    cuts <- cut_list[[m]]
    
    starts <- c(1L, cuts + 1L)
    ends   <- c(cuts, K)
    
    z_sorted <- integer(K)
    for (g in seq_len(k)) {
      z_sorted[starts[g]:ends[g]] <- g
    }
    
    # map back to original class order
    z <- integer(K)
    z[ord] <- z_sorted
    
    z_list[[m]] <- z
  }
  
  z_list
}

normalize_partition_labels <- function(z) {
  z <- as.integer(z)
  as.integer(match(z, sort(unique(z))))
}

unique_partition_list <- function(z_list) {
  if (!length(z_list)) {
    return(list())
  }

  z_norm <- lapply(z_list, normalize_partition_labels)
  key <- vapply(z_norm, function(v) paste(v, collapse = ","), character(1))
  z_norm[!duplicated(key)]
}

adjacent_merge_inits <- function(z) {
  z <- normalize_partition_labels(z)
  k <- max(z)

  if (k <= 1L) {
    return(list(rep(1L, length(z))))
  }

  out <- vector("list", k - 1L)
  for (g in seq_len(k - 1L)) {
    z_new <- z
    z_new[z == (g + 1L)] <- g
    z_new[z_new > (g + 1L)] <- z_new[z_new > (g + 1L)] - 1L
    out[[g]] <- as.integer(z_new)
  }

  key <- vapply(out, function(v) paste(v, collapse = ","), character(1))
  out[!duplicated(key)]
}

adjacent_merge_candidates_by_theta <- function(z, theta) {
  z <- normalize_partition_labels(z)
  k <- max(z)
  theta <- as.numeric(theta)

  if (length(theta) != k) {
    stop("length(theta) must equal the number of groups in z.")
  }

  if (k <= 1L) {
    return(list(rep(1L, length(z))))
  }

  ord <- order(theta, seq_len(k))
  relabel <- integer(k)
  relabel[ord] <- seq_len(k)
  z_ord <- relabel[z]

  out <- vector("list", k - 1L)
  for (g in seq_len(k - 1L)) {
    z_new <- z_ord
    z_new[z_ord == (g + 1L)] <- g
    z_new[z_new > (g + 1L)] <- z_new[z_new > (g + 1L)] - 1L
    out[[g]] <- as.integer(z_new)
  }

  unique_partition_list(out)
}

score_partition <- function(s, w, z, eps = 1e-8) {
  z <- normalize_partition_labels(z)
  k <- max(z)
  theta <- compute_theta(s, w, z, k, eps)
  Q <- loglik_grouped(s, w, z, theta)

  list(
    k = k,
    z = z,
    theta = theta,
    beta_class = theta[z],
    Q = Q,
    converged = TRUE,
    iter = 0L
  )
}

score_partition_candidates <- function(s, w, p, N_eff, z_list, eps = 1e-8,
                                       select_by = c("BIC", "Q")) {
  select_by <- match.arg(select_by)
  z_list <- unique_partition_list(z_list)
  if (!length(z_list)) {
    stop("z_list is empty.")
  }

  candidates <- lapply(z_list, function(z) {
    cand <- score_partition(s, w, z, eps = eps)
    c(
      cand,
      list(
        BIC = -2 * cand$Q + log(N_eff) * cand$k,
        s = s,
        w = w,
        p = p
      )
    )
  })

  if (select_by == "Q") {
    scores <- vapply(candidates, function(obj) obj$Q, numeric(1))
    best_idx <- which.max(scores)
  } else {
    scores <- vapply(candidates, function(obj) obj$BIC, numeric(1))
    best_idx <- which.min(scores)
  }

  list(
    best = candidates[[best_idx]],
    loglik = vapply(candidates, function(obj) obj$Q, numeric(1)),
    candidates = candidates
  )
}

run_alt_from_init <- function(s, w, p, k, z_init,
                              maxit = 200,
                              tol = 1e-10,
                              eps = 1e-8) {
  z <- normalize_partition_labels(z_init)

  q_old <- -Inf
  converged <- FALSE
  last_it <- maxit

  for (it in seq_len(maxit)) {
    z <- repair_empty_groups(z, p, w, k)
    theta <- compute_theta(s, w, z, k, eps)

    ll_mat <- outer(s, log(theta)) + outer(w - s, log1p(-theta))
    z_new <- max.col(ll_mat, ties.method = "first")
    z_new <- repair_empty_groups(z_new, p, w, k)

    theta_new <- compute_theta(s, w, z_new, k, eps)
    q_new <- loglik_grouped(s, w, z_new, theta_new)

    if (abs(q_new - q_old) < tol) {
      z <- z_new
      theta <- theta_new
      q_old <- q_new
      converged <- TRUE
      last_it <- it
      break
    }

    z <- z_new
    theta <- theta_new
    q_old <- q_new
    last_it <- it
  }

  Q <- loglik_grouped(s, w, z, theta)
  list(
    k = k,
    z = z,
    theta = theta,
    beta_class = theta[z],
    Q = Q,
    converged = converged,
    iter = last_it
  )
}

# Refinement: item and model-level

refine_item_alt <- function(xj, gamma, k,
                            init_z_list = NULL,
                            maxit = 200,
                            tol = 1e-10,
                            eps = 1e-8) {
  obs <- !is.na(xj)
  x_obs <- xj[obs]
  G_obs <- gamma[obs, , drop = FALSE]
  N_eff <- sum(obs)

  if (N_eff == 0) stop("This item is completely missing.")

  K <- ncol(G_obs)
  stopifnot(k >= 1, k <= K)

  # sufficient statistics
  w <- colSums(G_obs)
  s <- colSums(G_obs * x_obs)
  p <- ifelse(w > 0, s / w, 0.5)
  p <- clip01(p, eps)

  if (k == 1L) {
    theta <- clip01(sum(s) / sum(w), eps)
    z <- rep(1L, K)
    Q <- loglik_grouped(s, w, z, theta)
    best <- list(
      k = k,
      z = z,
      theta = theta,
      beta_class = rep(theta, K),
      Q = Q,
      BIC = -2 * Q + log(N_eff) * k,
      s = s,
      w = w,
      p = p,
      converged = TRUE,
      iter = 1L
    )
    return(list(best = best, loglik = Q))
  }

  if (k == K) {
    z <- seq_len(K)
    theta <- p
    Q <- loglik_grouped(s, w, z, theta)
    best <- list(
      k = k,
      z = z,
      theta = theta,
      beta_class = theta,
      Q = Q,
      BIC = -2 * Q + log(N_eff) * k,
      s = s,
      w = w,
      p = p,
      converged = TRUE,
      iter = 1L
    )
    return(list(best = best, loglik = Q))
  }

  best <- NULL
  loglik <- NULL
  z_list <- if (is.null(init_z_list)) {
    init_partition_grid(p, w, k)
  } else {
    z_custom <- lapply(init_z_list, function(z0) {
      z_norm <- normalize_partition_labels(z0)
      if (length(z_norm) != K) {
        stop("Each init partition in init_z_list must have length equal to ncol(gamma).")
      }
      if (max(z_norm) != k) {
        stop("Each init partition in init_z_list must contain exactly k groups.")
      }
      z_norm
    })
    if (!length(z_custom)) {
      stop("init_z_list is empty.")
    }
    z_custom
  }
  
  for (z in z_list) {
    cand <- run_alt_from_init(
      s = s, w = w, p = p, k = k, z_init = z,
      maxit = maxit, tol = tol, eps = eps
    )

    cand <- c(
      cand,
      list(
        BIC = -2 * cand$Q + log(N_eff) * k,
        s = s,
        w = w,
        p = p
      )
    )
    
    loglik <- c(loglik, cand$Q)

    if (is.null(best) || cand$Q > best$Q) {
      best <- cand
    }
  }
  
  list(best = best, loglik = loglik)
}

refine_lca_binary <- function(X, gamma,
                              k_grid = seq_len(ncol(gamma)),
                              c_pen_grid = c(1,5,10),
                              c_pen = NULL,
                              partition_search = c("full_grid", "stepwise_merge"),
                              maxit = 200,
                              tol = 1e-10,
                              eps = 1e-8,
                              verbose = TRUE) {
  X <- as.matrix(X)
  gamma <- as.matrix(gamma)

  storage.mode(X) <- "double"
  storage.mode(gamma) <- "double"
  partition_search <- match.arg(partition_search)
  
  if (nrow(X) != nrow(gamma)) {
    stop("X and gamma must have the same number of rows.")
  }

  # normalize posterior just in case
  rs <- rowSums(gamma)
  gamma <- gamma / rs

  J <- ncol(X)
  K <- ncol(gamma)

  k_grid <- sort(unique(k_grid))
  k_grid <- k_grid[k_grid >= 1 & k_grid <= K]

  if (!is.null(c_pen)) {
    c_pen_grid <- c_pen
  }
  c_pen_grid <- as.numeric(c_pen_grid)
  c_pen_grid <- c_pen_grid[is.finite(c_pen_grid) & c_pen_grid > 0]
  if (!length(c_pen_grid)) {
    stop("c_pen_grid must contain at least one positive value.")
  }

  item_names <- colnames(X)
  if (is.null(item_names)) item_names <- paste0("item", seq_len(J))

  class_names <- colnames(gamma)
  if (is.null(class_names)) class_names <- paste0("class", seq_len(K))

  c_pen_names <- paste0("c_pen=", c_pen_grid)
  
  selected_k <- matrix(
    NA_integer_, nrow = J, ncol = length(c_pen_grid),
    dimnames = list(item_names, c_pen_names)
  )
  
  beta_refined <- array(
    NA_real_, dim = c(J, K, length(c_pen_grid)), 
    dimnames = list(item_names, class_names, c_pen_names)
  )
  
  group_id <- array(
    NA_integer_, dim = c(J, K, length(c_pen_grid)),
    dimnames = list(item_names, class_names, c_pen_names)
  )
  
  bic_table <- array(
    NA_real_, dim = c(J, length(k_grid), length(c_pen_grid)),
    dimnames = list(item_names, paste0("k=", k_grid), c_pen_names)
  )
  
  loglik_table <- matrix(
    NA_real_, nrow = J, ncol = length(k_grid),
    dimnames = list(item_names, paste0("k=", k_grid))
  )
  
  fits_all <- vector("list", J)
  names(fits_all) <- item_names
  
  for (j in seq_len(J)) {
    if (verbose) {
      message(sprintf("Refining %s (%d/%d)", item_names[j], j, J))
    }

    if (partition_search == "full_grid") {
      fits_j <- lapply(k_grid, function(k) {
        refine_item_alt(
          xj = X[, j],
          gamma = gamma,
          k = k,
          maxit = maxit,
          tol = tol,
          eps = eps
        )
      })
    } else {
      obs <- !is.na(X[, j])
      x_obs <- X[obs, j]
      G_obs <- gamma[obs, , drop = FALSE]
      N_eff <- sum(obs)

      if (N_eff == 0) {
        stop(sprintf("Item %s is completely missing.", item_names[j]))
      }

      w <- colSums(G_obs)
      s <- colSums(G_obs * x_obs)
      p <- ifelse(w > 0, s / w, 0.5)
      p <- clip01(p, eps)

      k_desc <- sort(seq.int(1L, K), decreasing = TRUE)
      fits_desc <- vector("list", length(k_desc))
      prev_best <- NULL

      for (idx in seq_along(k_desc)) {
        k <- k_desc[idx]

        if (k == K) {
          z_list <- list(seq_len(K))
        } else {
          z_list <- adjacent_merge_candidates_by_theta(
            z = prev_best$z,
            theta = prev_best$theta
          )
        }

        fits_desc[[idx]] <- score_partition_candidates(
          s = s,
          w = w,
          p = p,
          N_eff = N_eff,
          z_list = z_list,
          eps = eps,
          select_by = "Q"
        )

        prev_best <- fits_desc[[idx]]$best
      }
      fits_map <- setNames(fits_desc, paste0("k=", k_desc))
      fits_j <- lapply(k_grid, function(k) fits_map[[paste0("k=", k)]])
    }
    names(fits_j) <- paste0("k=", k_grid)
    
    bics <- vapply(fits_j, function(obj) obj$best$BIC, numeric(1))
    
    for (c in seq_along(c_pen_grid)) {
      bics_c <- bics + 2*log(c_pen_grid[c]) * k_grid
      bic_table[j, , c] <- bics_c
      
      best <- fits_j[[which.min(bics_c)]]$best
      
      selected_k[j,c] <- best$k
      beta_refined[j, ,c] <- best$beta_class
      group_id[j, ,c] <- best$z
    }
    
    loglik <- vapply(fits_j, function(obj) obj$best$Q, numeric(1))
    loglik_table[j, ] <- loglik
    fits_all[[j]] <- fits_j
  }

  list(
    selected_k = selected_k,
    beta_refined = beta_refined,
    group_id = group_id,
    bic_table = bic_table,
    loglik_table = loglik_table,
    fits = fits_all,
    partition_search = partition_search
  )
}


logsumexp_vec <- function(x) {
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Sparse LCA re-estimation

compute_sparse_beta_from_group <- function(X, gamma, group_id, eps = 1e-8) {
  X <- as.matrix(X)
  gamma <- as.matrix(gamma)
  group_id <- as.matrix(group_id)

  storage.mode(X) <- "double"
  storage.mode(gamma) <- "double"
  storage.mode(group_id) <- "integer"

  N <- nrow(X)
  J <- ncol(X)
  K <- ncol(gamma)

  if (nrow(gamma) != N) stop("X and gamma must have the same number of rows.")
  if (nrow(group_id) != J || ncol(group_id) != K) {
    stop("group_id must be a J x K matrix.")
  }

  item_names <- colnames(X)
  if (is.null(item_names)) item_names <- paste0("item", seq_len(J))
  class_names <- colnames(gamma)
  if (is.null(class_names)) class_names <- paste0("class", seq_len(K))

  beta <- matrix(NA_real_, nrow = J, ncol = K,
                 dimnames = list(item_names, class_names))
  theta_list <- vector("list", J)
  names(theta_list) <- item_names

  for (j in seq_len(J)) {
    obs <- !is.na(X[, j])
    xj <- X[obs, j]
    G  <- gamma[obs, , drop = FALSE]

    s <- as.numeric(crossprod(G, xj))
    w <- colSums(G)

    z <- as.integer(group_id[j, ])
    groups <- sort(unique(z))
    theta_j <- numeric(length(groups))
    names(theta_j) <- paste0("g", groups)

    for (ii in seq_along(groups)) {
      g <- groups[ii]
      idx <- which(z == g)
      sw <- sum(w[idx])
      theta_g <- if (sw > 0) sum(s[idx]) / sw else 0.5
      theta_g <- min(1 - eps, max(eps, theta_g))
      theta_j[ii] <- theta_g
      beta[j, idx] <- theta_g
    }

    theta_list[[j]] <- theta_j
  }

  list(beta = beta, theta = theta_list)
}

compute_sparse_posterior <- function(X, eta, beta, eps = 1e-8) {
  X <- as.matrix(X)
  beta <- as.matrix(beta)

  storage.mode(X) <- "double"
  storage.mode(beta) <- "double"

  N <- nrow(X)
  J <- ncol(X)
  K <- length(eta)

  if (nrow(beta) != J || ncol(beta) != K) {
    stop("beta must be a J x K matrix.")
  }

  eta <- as.numeric(eta)
  eta <- eta / sum(eta)
  eta <- pmax(eta, eps)
  eta <- eta / sum(eta)

  log_eta <- log(eta)
  log_beta <- log(pmin(pmax(beta, eps), 1 - eps))
  log_1mb  <- log1p(-pmin(pmax(beta, eps), 1 - eps))

  log_post <- matrix(0, nrow = N, ncol = K)
  loglik_i <- numeric(N)

  for (i in seq_len(N)) {
    obs <- !is.na(X[i, ])
    xi <- X[i, obs]
    if (!length(xi)) {
      lp <- log_eta
    } else {
      lp <- log_eta +
        colSums(matrix(xi, nrow = sum(obs), ncol = K) * log_beta[obs, , drop = FALSE] +
                  matrix(1 - xi, nrow = sum(obs), ncol = K) * log_1mb[obs, , drop = FALSE])
    }
    lse <- logsumexp_vec(lp)
    log_post[i, ] <- lp - lse
    loglik_i[i] <- lse
  }

  posterior <- exp(log_post)
  colnames(posterior) <- colnames(beta)
  list(
    posterior = posterior,
    loglik = sum(loglik_i),
    loglik_i = loglik_i
  )
}

reestimate_sparse_lca <- function(X, group_id,
                                  eta_init = NULL,
                                  beta_init = NULL,
                                  gamma_init = NULL,
                                  maxit = 300,
                                  tol = 1e-8,
                                  eps = 1e-8,
                                  verbose = TRUE) {
  X <- as.matrix(X)
  group_id <- as.matrix(group_id)

  storage.mode(X) <- "double"
  storage.mode(group_id) <- "integer"

  N <- nrow(X)
  J <- ncol(X)
  K <- ncol(group_id)

  if (nrow(group_id) != J) {
    stop("group_id must have nrow(group_id) = ncol(X).")
  }

  item_names <- colnames(X)
  if (is.null(item_names)) item_names <- paste0("item", seq_len(J))
  class_names <- colnames(group_id)
  if (is.null(class_names)) class_names <- paste0("class", seq_len(K))

  rownames(group_id) <- item_names
  colnames(group_id) <- class_names

  if (!is.null(gamma_init)) {
    gamma <- as.matrix(gamma_init)
    if (nrow(gamma) != N || ncol(gamma) != K) {
      stop("gamma_init must be an N x K matrix.")
    }
    rs <- rowSums(gamma)
    gamma <- gamma / rs
    eta <- colMeans(gamma)
    beta_obj <- compute_sparse_beta_from_group(X, gamma, group_id, eps = eps)
    beta <- beta_obj$beta
  } else {
    if (is.null(eta_init)) eta_init <- rep(1 / K, K)
    eta <- as.numeric(eta_init)
    eta <- eta / sum(eta)

    if (is.null(beta_init)) {
      # crude start from marginal means
      p0 <- colMeans(X, na.rm = TRUE)
      p0[is.na(p0)] <- 0.5
      beta <- matrix(rep(p0, K), nrow = J, ncol = K)
      rownames(beta) <- item_names
      colnames(beta) <- class_names
      # project to sparse structure
      fake_gamma <- matrix(rep(eta, each = N), nrow = N, ncol = K)
      beta <- compute_sparse_beta_from_group(X, fake_gamma, group_id, eps = eps)$beta
    } else {
      beta <- as.matrix(beta_init)
      if (nrow(beta) != J || ncol(beta) != K) {
        stop("beta_init must be a J x K matrix.")
      }
      rownames(beta) <- item_names
      colnames(beta) <- class_names
      # enforce the supplied sparse structure at the start
      post0 <- compute_sparse_posterior(X, eta, beta, eps = eps)$posterior
      beta <- compute_sparse_beta_from_group(X, post0, group_id, eps = eps)$beta
    }
  }

  ll_trace <- numeric(maxit)
  converged <- FALSE
  posterior <- NULL

  for (it in seq_len(maxit)) {
    # E-step
    estep <- compute_sparse_posterior(X, eta, beta, eps = eps)
    posterior <- estep$posterior
    ll_new <- estep$loglik
    ll_trace[it] <- ll_new

    # M-step: eta
    eta <- colMeans(posterior)
    eta <- pmax(eta, eps)
    eta <- eta / sum(eta)

    # M-step: beta under fixed sparse structure
    beta_obj <- compute_sparse_beta_from_group(X, posterior, group_id, eps = eps)
    beta <- beta_obj$beta

    if (verbose) {
      message(sprintf("sparse-EM iter %d: logLik = %.8f", it, ll_new))
    }

    if (it > 1 && abs(ll_trace[it] - ll_trace[it - 1]) < tol) {
      converged <- TRUE
      ll_trace <- ll_trace[seq_len(it)]
      break
    }
  }

  if (!converged) {
    ll_trace <- ll_trace[seq_len(maxit)]
  }

  list(
    eta = eta,
    beta = beta,
    posterior = posterior,
    group_id = group_id,
    theta = beta_obj$theta,
    logLik = tail(ll_trace, 1),
    logLik_trace = ll_trace,
    converged = converged,
    iter = length(ll_trace)
  )
}

# Simulation helper functions

sample_dirichlet <- function(alpha, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  alpha <- as.numeric(alpha)
  if (any(alpha <= 0)) {
    stop("All Dirichlet concentration parameters must be positive.")
  }

  x <- rgamma(length(alpha), shape = alpha, rate = 1)
  x / sum(x)
}


generate_lca_data_from_true_parameters <- function(
    N,
    Pi_true,
    eta,
    m_true = NULL,
    group_true = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)

  Pi_true <- as.matrix(Pi_true)
  J <- nrow(Pi_true)
  K <- ncol(Pi_true)

  eta <- as.numeric(eta)
  if (length(eta) != K) {
    stop("Length of eta must equal ncol(Pi_true).")
  }
  eta <- eta / sum(eta)

  class_id <- sample(seq_len(K), size = N, replace = TRUE, prob = eta)

  X <- matrix(NA_integer_, nrow = N, ncol = J)
  for (i in seq_len(N)) {
    X[i, ] <- rbinom(J, size = 1, prob = Pi_true[, class_id[i]])
  }

  colnames(X) <- rownames(Pi_true)

  list(
    X = X,
    class_id = class_id,
    eta = eta,
    Pi_true = Pi_true,
    m_true = m_true,
    group_true = group_true
  )
}



generate_true_item_parameters <- function(
    J = 64,
    K = 8,
    seed = NULL,
    low_range = c(0.18, 0.28),
    gap2_range = c(0.48, 0.58),
    gap3_range = c(0.24, 0.31),
    max_prob = 0.88
) {
  if (!is.null(seed)) set.seed(seed)

  if (K == 4) {
    class_profiles <- matrix(c(
      0, 0,
      0, 1,
      1, 0,
      1, 1
    ), nrow = 4, byrow = TRUE)

  } else if (K == 8) {
    class_profiles <- matrix(c(
      0, 0, 0,
      0, 0, 1,
      0, 1, 0,
      1, 0, 0,
      0, 1, 1,
      1, 0, 1,
      1, 1, 0,
      1, 1, 1
    ), nrow = 8, byrow = TRUE)

  } else {
    stop("Only K = 4 or K = 8 is supported.")
  }

  D <- ncol(class_profiles)

  colnames(class_profiles) <- paste0("alpha", seq_len(D))
  rownames(class_profiles) <- paste0("class", seq_len(K))

  beta_true <- matrix(
    0,
    nrow = J,
    ncol = D + 1,
    dimnames = list(
      paste0("item", seq_len(J)),
      c("beta0", paste0("beta", seq_len(D)))
    )
  )

  Pi_true <- matrix(
    NA_real_,
    nrow = J,
    ncol = K,
    dimnames = list(
      paste0("item", seq_len(J)),
      paste0("class", seq_len(K))
    )
  )

  group_true <- matrix(
    NA_integer_,
    nrow = J,
    ncol = K,
    dimnames = dimnames(Pi_true)
  )

  level_probs <- vector("list", J)

  if (K == 4) {
    m_true <- rep(2L, J)
  } else {
    m_true <- c(rep(3L, floor(J / 4 )), rep(2L, J - floor(J / 4)))
    m_true <- sample(m_true)
  }

  active_dim <- rep(NA_integer_, J)
  zero_dim <- rep(NA_integer_, J)
  partition_type <- rep(NA_character_, J)
  singleton_class <- rep(NA_integer_, J)

  for (j in seq_len(J)) {

    p_low <- runif(1, low_range[1], low_range[2])

    if (m_true[j] == 2L) {

      gap_max <- min(gap2_range[2], max_prob - p_low)
      gap <- runif(1, gap2_range[1], gap_max)

      if (K == 4) {
        item_idx_2 <- sum(m_true[seq_len(j)] == 2L)
        partition_type[j] <- rep(c("2-2", "2-2", "1-3"), length.out = J)[item_idx_2]

        if (partition_type[j] == "2-2") {
          r <- ((sum(partition_type[seq_len(j)] == "2-2", na.rm = TRUE) - 1L) %% D) + 1L
          active_dim[j] <- r
          score <- class_profiles[, r]
          beta_true[j, "beta0"] <- p_low
          beta_true[j, paste0("beta", r)] <- gap
        } else {
          one_three_idx <- sum(partition_type[seq_len(j)] == "1-3", na.rm = TRUE)
          singled <- c(1L, 4L)[((one_three_idx - 1L) %% 2L) + 1L]
          singleton_class[j] <- singled
          score <- if (singled == 1L) {
            c(0L, 1L, 1L, 1L)
          } else {
            c(0L, 0L, 0L, 1L)
          }
          beta_true[j, "beta0"] <- p_low
        }
      } else {
        partition_type[j] <- "4-4"
        r <- ((sum(m_true[seq_len(j)] == 2L) - 1) %% D) + 1
        active_dim[j] <- r
        score <- class_profiles[, r]
        beta_true[j, "beta0"] <- p_low
        beta_true[j, paste0("beta", r)] <- gap
      }

      Pi_true[j, ] <- p_low + gap * score
      group_true[j, ] <- 1L + score
      level_probs[[j]] <- c(low = p_low, high = p_low + gap)

    } else {

      r <- ((sum(m_true[seq_len(j)] == 3L) - 1) %% D) + 1
      zero_dim[j] <- r
      others <- setdiff(seq_len(D), r)

      gap_max <- min(gap3_range[2], (max_prob - p_low) / 2)
      gap <- runif(1, gap3_range[1], gap_max)

      score <- rowSums(class_profiles[, others, drop = FALSE])

      beta_true[j, "beta0"] <- p_low
      beta_true[j, paste0("beta", others)] <- gap

      Pi_true[j, ] <- p_low + gap * score
      group_true[j, ] <- 1L + score
      level_probs[[j]] <- c(
        low = p_low,
        mid = p_low + gap,
        high = p_low + 2 * gap
      )
    }
  }

  names(level_probs) <- paste0("item", seq_len(J))

  out <- list(
    Pi_true = Pi_true,
    beta_true = beta_true,
    m_true = m_true,
    group_true = group_true,
    active_dim = active_dim,
    partition_type = partition_type,
    singleton_class = singleton_class,
    level_probs = level_probs
  )

  if (K == 8) {
    out$zero_dim <- zero_dim
  }

  out
}
