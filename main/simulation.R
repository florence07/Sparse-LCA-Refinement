
# Utilities: extraction & permutations


extract_polca_prob <- function(fit) {
  J <- length(fit$probs)
  K <- length(fit$P)
  out <- matrix(NA_real_, nrow = J, ncol = K)
  for (j in seq_len(J)) {
    out[j, ] <- fit$probs[[j]][, 2]
  }
  rownames(out) <- names(fit$probs)
  colnames(out) <- paste0("class", seq_len(K))
  out
}

all_permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))

  out <- lapply(seq_along(x), function(i) {
    rest <- x[-i]
    subp <- all_permutations(rest)
    cbind(x[i], subp)
  })
  do.call(rbind, out)
}

best_label_switch <- function(true_vec, est_vec) {
  stopifnot(length(true_vec) == length(est_vec))

  labels <- sort(unique(c(true_vec, est_vec)))
  g <- length(labels)

  true_idx <- match(true_vec, labels)
  est_idx <- match(est_vec, labels)

  if (g == 1L) {
    return(list(
      switched = est_vec,
      accuracy = 1,
      permutation = labels
    ))
  }

  perms <- all_permutations(seq_len(g))
  best_acc <- -Inf
  best_perm <- perms[1, ]

  for (r in seq_len(nrow(perms))) {
    perm <- perms[r, ]
    switched_idx <- perm[est_idx]
    acc <- mean(switched_idx == true_idx)
    if (acc > best_acc) {
      best_acc <- acc
      best_perm <- perm
    }
  }

  switched <- labels[best_perm[est_idx]]
  names(best_perm) <- as.character(labels)

  list(
    switched = switched,
    accuracy = best_acc,
    permutation = best_perm
  )
}

permute_prob_columns <- function(prob_mat, class_perm) {
  prob_mat <- as.matrix(prob_mat)
  if (length(class_perm) != ncol(prob_mat)) return(prob_mat)

  out <- prob_mat[, class_perm, drop = FALSE]
  colnames(out) <- colnames(prob_mat)
  out
}

best_prob_permutation <- function(true_prob, est_prob) {
  true_prob <- as.matrix(true_prob)
  est_prob <- as.matrix(est_prob)

  if (!all(dim(true_prob) == dim(est_prob))) {
    return(list(
      permutation = seq_len(ncol(est_prob)),
      est_permuted = est_prob,
      mse = NA_real_
    ))
  }

  K <- ncol(true_prob)
  perms <- all_permutations(seq_len(K))
  best_perm <- perms[1, ]
  best_mse <- Inf

  for (r in seq_len(nrow(perms))) {
    perm <- perms[r, ]
    mse_r <- mean((true_prob - est_prob[, perm, drop = FALSE])^2)
    if (mse_r < best_mse) {
      best_mse <- mse_r
      best_perm <- perm
    }
  }

  list(
    permutation = best_perm,
    est_permuted = est_prob[, best_perm, drop = FALSE],
    mse = best_mse
  )
}


# Utilities: alignment & metrics

make_empty_align_result <- function(n_item, n_class) {
  list(
    switched = matrix(NA_integer_, nrow = n_item, ncol = n_class),
    item_accuracy = rep(NA_real_, n_item),
    overall_accuracy = NA_real_
  )
}

align_group_with_fixed_classes <- function(true_group, est_group_aligned) {
  true_group <- as.matrix(true_group)
  est_group_aligned <- as.matrix(est_group_aligned)

  if (!all(dim(true_group) == dim(est_group_aligned))) {
    return(make_empty_align_result(nrow(true_group), ncol(true_group)))
  }

  J <- nrow(true_group)
  switched <- matrix(NA_integer_, nrow = J, ncol = ncol(true_group))
  item_accuracy <- numeric(J)

  for (j in seq_len(J)) {
    obj <- best_label_switch(true_group[j, ], est_group_aligned[j, ])
    switched[j, ] <- as.integer(obj$switched)
    item_accuracy[j] <- obj$accuracy
  }

  rownames(switched) <- rownames(true_group)
  colnames(switched) <- colnames(true_group)

  list(
    switched = switched,
    item_accuracy = item_accuracy,
    overall_accuracy = mean(switched == true_group)
  )
}

adjusted_rand_index <- function(true_labels, est_labels) {
  true_labels <- as.integer(true_labels)
  est_labels <- as.integer(est_labels)

  ok <- !(is.na(true_labels) | is.na(est_labels))
  true_labels <- true_labels[ok]
  est_labels <- est_labels[ok]

  n <- length(true_labels)
  if (n < 2L) return(NA_real_)

  tab <- table(true_labels, est_labels)
  nij2 <- sum(tab * (tab - 1) / 2)
  ai <- rowSums(tab)
  bj <- colSums(tab)
  ai2 <- sum(ai * (ai - 1) / 2)
  bj2 <- sum(bj * (bj - 1) / 2)
  n2 <- n * (n - 1) / 2

  expected <- (ai2 * bj2) / n2
  max_index <- 0.5 * (ai2 + bj2)
  denom <- max_index - expected

  if (denom == 0) return(1)
  (nij2 - expected) / denom
}

safe_mean <- function(x) {
  if (all(is.na(x))) return(NA_real_)
  mean(x, na.rm = TRUE)
}

classify_mj_selection <- function(m_true, selected_mj) {
  m_true <- as.integer(m_true)
  selected_mj <- as.integer(selected_mj)
  out <- rep(NA_character_, length(selected_mj))
  ok <- !(is.na(m_true) | is.na(selected_mj))
  diff <- selected_mj - m_true

  out[ok & diff == 0L] <- "correct"
  out[ok & diff > 0L] <- "overselect"
  out[ok & diff < 0L] <- "underselect"
  out
}

summarize_mj_selection_by_item <- function(item_mj_selection) {
  if (is.null(item_mj_selection) || nrow(item_mj_selection) == 0L) {
    return(data.frame())
  }

  key_cols <- c("N", "J", "k", "c_pen", "item", "m_true")
  missing_cols <- setdiff(key_cols, names(item_mj_selection))
  if (length(missing_cols) > 0L) {
    stop(sprintf(
      "item_mj_selection is missing columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  row_key <- do.call(
    paste,
    c(item_mj_selection[key_cols], list(sep = "\r"))
  )

  out <- lapply(split(item_mj_selection, row_key, drop = TRUE), function(df) {
    selected_mj <- as.integer(df$selected_mj)
    m_true <- as.integer(df$m_true)
    mj_difference <- selected_mj - m_true
    n_valid <- sum(!is.na(selected_mj) & !is.na(m_true))
    n_correct <- sum(mj_difference == 0L, na.rm = TRUE)
    n_overselect <- sum(mj_difference > 0L, na.rm = TRUE)
    n_underselect <- sum(mj_difference < 0L, na.rm = TRUE)
    n_mj_error <- n_overselect + n_underselect

    data.frame(
      N = df$N[1],
      J = df$J[1],
      k = df$k[1],
      c_pen = df$c_pen[1],
      item = df$item[1],
      m_true = df$m_true[1],
      n_reps = nrow(df),
      n_valid = n_valid,
      n_correct = n_correct,
      n_mj_error = n_mj_error,
      n_overselect = n_overselect,
      n_underselect = n_underselect,
      correct_rate = if (n_valid > 0L) n_correct / n_valid else NA_real_,
      overselect_rate = if (n_valid > 0L) n_overselect / n_valid else NA_real_,
      underselect_rate = if (n_valid > 0L) n_underselect / n_valid else NA_real_,
      overselect_share_among_errors = if (n_mj_error > 0L) {
        n_overselect / n_mj_error
      } else {
        NA_real_
      },
      underselect_share_among_errors = if (n_mj_error > 0L) {
        n_underselect / n_mj_error
      } else {
        NA_real_
      },
      mean_signed_mj_error = if (n_valid > 0L) {
        mean(mj_difference, na.rm = TRUE)
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out[order(out$N, out$J, out$k, out$c_pen, out$item), , drop = FALSE]
}


# Utilities: refine-output normalization

normalize_refine_output <- function(ref_sim, c_pen_grid, K_expected) {
  selected_k <- ref_sim$selected_k
  if (is.null(dim(selected_k))) {
    selected_k <- matrix(as.integer(selected_k), ncol = 1)
    rownames(selected_k) <- names(ref_sim$selected_k)
  } else {
    selected_k <- as.matrix(selected_k)
    storage.mode(selected_k) <- "integer"
  }

  c_pen_names <- colnames(selected_k)
  if (is.null(c_pen_names)) {
    c_pen_names <- paste0("c_pen=", c_pen_grid[seq_len(ncol(selected_k))])
    colnames(selected_k) <- c_pen_names
  }

  beta_refined <- ref_sim$beta_refined
  if (length(dim(beta_refined)) == 2L) {
    beta_refined <- array(
      beta_refined,
      dim = c(nrow(beta_refined), ncol(beta_refined), 1L),
      dimnames = list(rownames(beta_refined), colnames(beta_refined), c_pen_names[1])
    )
  }

  group_id <- ref_sim$group_id
  if (length(dim(group_id)) == 2L) {
    group_id <- array(
      group_id,
      dim = c(nrow(group_id), ncol(group_id), 1L),
      dimnames = list(rownames(group_id), colnames(group_id), c_pen_names[1])
    )
  }

  bic_table <- ref_sim$bic_table
  if (length(dim(bic_table)) == 2L) {
    bic_table <- array(
      bic_table,
      dim = c(nrow(bic_table), ncol(bic_table), 1L),
      dimnames = list(rownames(bic_table), colnames(bic_table), c_pen_names[1])
    )
  }

  loglik_table <- as.matrix(ref_sim$loglik_table)

  stopifnot(ncol(beta_refined) == K_expected)

  list(
    selected_k = selected_k,
    beta_refined = beta_refined,
    group_id = group_id,
    bic_table = bic_table,
    loglik_table = loglik_table,
    c_pen_names = c_pen_names
  )
}

fit_polca_unrestricted_grid <- function(f_sim, X12, class_grid, nrep, maxiter,
                                        verbose = TRUE, N = NA_integer_,
                                        seed = NA_integer_) {
  fits <- vector("list", length(class_grid))
  names(fits) <- paste0("K=", class_grid)

  for (idx in seq_along(class_grid)) {
    k <- class_grid[idx]
    if (verbose) {
      message(sprintf("[N %d][seed %d] Fitting unrestricted poLCA with k = %d", N, seed, k))
    }

    fits[[idx]] <- poLCA::poLCA(
      formula = f_sim,
      data = X12,
      nclass = k,
      nrep = nrep,
      maxiter = maxiter,
      verbose = FALSE
    )
  }

  fit_grid <- data.frame(
    candidate_K = class_grid,
    BIC = vapply(fits, function(fit) fit$bic, numeric(1)),
    logLik = vapply(fits, function(fit) fit$llik, numeric(1)),
    npar = vapply(fits, function(fit) {
      if (is.null(fit$npar)) NA_real_ else fit$npar
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  fit_grid$selected_by_BIC <- fit_grid$candidate_K == fit_grid$candidate_K[which.min(fit_grid$BIC)]

  list(
    fits = fits,
    fit_grid = fit_grid,
    selected_K = fit_grid$candidate_K[which.min(fit_grid$BIC)]
  )
}

make_unrestricted_parameter_tables <- function(fit, seed, N, J, K, true_eta, true_beta) {
  beta_hat <- extract_polca_prob(fit)
  align <- best_prob_permutation(true_beta, beta_hat)
  class_perm <- align$permutation
  beta_aligned <- permute_prob_columns(beta_hat, class_perm)

  nu_hat <- as.numeric(fit$P)
  nu_aligned <- if (length(class_perm) == length(nu_hat)) {
    nu_hat[class_perm]
  } else {
    rep(NA_real_, length(true_eta))
  }

  beta_long <- do.call(
    rbind,
    lapply(seq_len(nrow(beta_aligned)), function(j) {
      data.frame(
        seed = seed,
        N = N,
        J = J,
        K = K,
        item = rownames(beta_aligned)[j],
        class = colnames(beta_aligned),
        beta_true = as.numeric(true_beta[j, ]),
        beta_unrestricted = as.numeric(beta_aligned[j, ]),
        squared_error = as.numeric((true_beta[j, ] - beta_aligned[j, ])^2),
        stringsAsFactors = FALSE
      )
    })
  )
  rownames(beta_long) <- NULL

  nu_df <- data.frame(
    seed = seed,
    N = N,
    J = J,
    K = K,
    class = paste0("class", seq_along(true_eta)),
    nu_true = as.numeric(true_eta),
    nu_unrestricted = nu_aligned,
    squared_error = (as.numeric(true_eta) - nu_aligned)^2,
    stringsAsFactors = FALSE
  )

  list(
    beta = beta_long,
    nu = nu_df,
    beta_aligned = beta_aligned,
    nu_aligned = nu_aligned,
    class_perm = class_perm,
    MSE_beta = if (is.finite(align$mse)) align$mse else NA_real_,
    MSE_nu = mean(nu_df$squared_error, na.rm = TRUE)
  )
}

stack_unrestricted_tables <- function(unrestricted_params) {
  nu_tbl <- data.frame(
    seed = unrestricted_params$nu$seed,
    N = unrestricted_params$nu$N,
    J = unrestricted_params$nu$J,
    K = unrestricted_params$nu$K,
    param_type = "nu",
    item = NA_character_,
    class = unrestricted_params$nu$class,
    true_value = unrestricted_params$nu$nu_true,
    estimate = unrestricted_params$nu$nu_unrestricted,
    squared_error = unrestricted_params$nu$squared_error,
    stringsAsFactors = FALSE
  )

  beta_tbl <- data.frame(
    seed = unrestricted_params$beta$seed,
    N = unrestricted_params$beta$N,
    J = unrestricted_params$beta$J,
    K = unrestricted_params$beta$K,
    param_type = "beta",
    item = unrestricted_params$beta$item,
    class = unrestricted_params$beta$class,
    true_value = unrestricted_params$beta$beta_true,
    estimate = unrestricted_params$beta$beta_unrestricted,
    squared_error = unrestricted_params$beta$squared_error,
    stringsAsFactors = FALSE
  )

  rbind(nu_tbl, beta_tbl)
}

# Main: single simulation

run_one_simulation <- function(
    seed = 500,
    N = 500,
    J = 12,
    K = 4,
    eta = rep(1 / K, K),
    true_params = NULL,
    class_grid = NULL,
    c_pen_grid = c(1,5,10,15,20),
    polca_nrep = 100,
    polca_maxiter = 5000,
    partition_search = "full_grid",
    refine_maxit = 300,
    refine_tol = 1e-10,
    sparse_maxit = 300,
    sparse_tol = 1e-8,
    dir = "output/simulations",
    verbose = TRUE
) {
  if (is.null(class_grid)) {
    class_grid <- seq.int(max(2L, K - 1L), K + 1L)
  }

  stopifnot(!is.null(true_params))
  
  sim <- generate_lca_data_from_true_parameters(
    N = N,
    Pi_true = true_params$Pi_true,
    eta = eta,
    m_true = true_params$m_true,
    group_true = true_params$group_true,
    seed = seed
  )

  X <- as.matrix(sim$X)
  stopifnot(all(X %in% c(0, 1)))

  # poLCA expects categories 1,2 for binary data
  X12 <- as.data.frame(X + 1L)
  names(X12) <- colnames(X)

  f_sim <- as.formula(
    paste0("cbind(", paste(names(X12), collapse = ", "), ") ~ 1")
  )

  out_dir <- file.path(dir, sprintf("N_%d", N))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  summary_file <- sprintf("%s/summary_seed%d.csv", out_dir, seed)
  unrestricted_grid_file <- sprintf("%s/unrestricted_fit_grid_seed%d.csv", out_dir, seed)
  unrestricted_params_file <- sprintf("%s/unrestricted_params_seed%d.csv", out_dir, seed)
  item_selection_file <- sprintf("%s/item_mj_selection_seed%d.csv", out_dir, seed)

  class_grid <- sort(unique(as.integer(class_grid)))
  class_grid <- class_grid[class_grid >= 2L]
  if (!K %in% class_grid) {
    class_grid <- sort(unique(c(class_grid, K)))
  }

  unrestricted <- fit_polca_unrestricted_grid(
    f_sim = f_sim,
    X12 = X12,
    class_grid = class_grid,
    nrep = polca_nrep,
    maxiter = polca_maxiter,
    verbose = verbose,
    N = N,
    seed = seed
  )
  best_k <- unrestricted$selected_K
  best_fit <- unrestricted$fits[[paste0("K=", K)]]

  sim_prob <- sim$Pi_true
  unrestricted_params <- make_unrestricted_parameter_tables(
    fit = best_fit,
    seed = seed,
    N = N,
    J = J,
    K = K,
    true_eta = eta,
    true_beta = sim_prob
  )
  mse_beta_step1 <- unrestricted_params$MSE_beta
  mse_nu_step1 <- unrestricted_params$MSE_nu
  polca_prob_step1 <- unrestricted_params$beta_aligned

  unrestricted_fit_grid <- unrestricted$fit_grid
  unrestricted_fit_grid$seed <- seed
  unrestricted_fit_grid$N <- N
  unrestricted_fit_grid$J <- J
  unrestricted_fit_grid$true_K <- K
  unrestricted_fit_grid$selected_true_K <- unrestricted_fit_grid$candidate_K == K & unrestricted_fit_grid$selected_by_BIC

  unrestricted_param_table <- stack_unrestricted_tables(unrestricted_params)

  ref_sim <- refine_lca_binary(
    X = X,
    gamma = best_fit$posterior,
    k_grid = seq_len(ncol(best_fit$posterior)),
    c_pen_grid = c_pen_grid,
    partition_search = partition_search,
    maxit = refine_maxit,
    tol = refine_tol,
    verbose = verbose
  )

  ref_norm <- normalize_refine_output(ref_sim, c_pen_grid, ncol(best_fit$posterior))
  selected_k_mat <- ref_norm$selected_k
  beta_refined_arr <- ref_norm$beta_refined
  group_id_arr <- ref_norm$group_id
  bic_table_arr <- ref_norm$bic_table
  loglik_table <- ref_norm$loglik_table
  c_pen_names <- ref_norm$c_pen_names
  c_pen_vals <- suppressWarnings(as.numeric(sub("^.*=", "", c_pen_names)))
  c_pen_vals[is.na(c_pen_vals)] <- c_pen_grid[seq_len(sum(is.na(c_pen_vals)))]

  polca_prob_base <- extract_polca_prob(best_fit)
  item_names <- rownames(selected_k_mat)
  if (is.null(item_names)) item_names <- paste0("item", seq_len(nrow(selected_k_mat)))
  k_seq <- seq_len(ncol(best_fit$posterior))

  summary_rows <- vector("list", ncol(selected_k_mat))
  item_mj_selection_list <- vector("list", ncol(selected_k_mat))

  refine_candidates_df <- do.call(
    rbind,
    lapply(seq_len(ncol(selected_k_mat)), function(c_idx) {
      c_pen_val <- c_pen_vals[c_idx]

      do.call(
        rbind,
        lapply(seq_len(nrow(selected_k_mat)), function(j) {
          selected_mj <- as.integer(selected_k_mat[j, c_idx])
          m_true_j <- as.integer(sim$m_true[j])
          data.frame(
            c_pen = c_pen_val,
            item = item_names[j],
            m_true = m_true_j,
            candidate_mj = as.integer(k_seq),
            logLik = as.numeric(loglik_table[j, k_seq]),
            BIC = as.numeric(bic_table_arr[j, k_seq, c_idx]),
            selected_mj = selected_mj,
            selected_mj_difference = selected_mj - m_true_j,
            selected_mj_type = classify_mj_selection(m_true_j, selected_mj),
            is_selected = as.integer(k_seq == selected_mj),
            stringsAsFactors = FALSE
          )
        })
      )
    })
  )

  for (c_idx in seq_len(ncol(selected_k_mat))) {
    c_pen_val <- c_pen_vals[c_idx]

    selected_k_vec <- as.integer(selected_k_mat[, c_idx])
    beta_refined_mat <- beta_refined_arr[, , c_idx]
    group_id_mat <- group_id_arr[, , c_idx]

    sfit_sim <- reestimate_sparse_lca(
      X = X,
      group_id = group_id_mat,
      eta_init = best_fit$P,
      beta_init = beta_refined_mat,
      gamma_init = best_fit$posterior,
      maxit = sparse_maxit,
      tol = sparse_tol,
      verbose = verbose
    )

    prob_align <- best_prob_permutation(sim_prob, sfit_sim$beta)
    class_perm_prob <- prob_align$permutation
    group_id_perm <- if (length(class_perm_prob) == ncol(group_id_mat)) {
      group_id_mat[, class_perm_prob, drop = FALSE]
    } else {
      group_id_mat
    }

    align_obj <- align_group_with_fixed_classes(
      true_group = sim$group_true,
      est_group_aligned = group_id_perm
    )

    sfit_prob <- permute_prob_columns(sfit_sim$beta, class_perm_prob)

    mse_per_item <- if (all(dim(sim_prob) == dim(sfit_prob))) {
      rowMeans((sim_prob - sfit_prob)^2)
    } else {
      rep(NA_real_, nrow(sim_prob))
    }

    mse_beta <- if (all(dim(sim_prob) == dim(sfit_prob))) {
      if (is.finite(prob_align$mse)) prob_align$mse else mean((sim_prob - sfit_prob)^2)
    } else {
      NA_real_
    }

    mse_nu <- if (length(class_perm_prob) == length(eta)) {
      mean((as.numeric(eta) - as.numeric(sfit_sim$eta)[class_perm_prob])^2)
    } else {
      NA_real_
    }

    ari_per_item <- if (all(dim(sim$group_true) == dim(align_obj$switched))) {
      vapply(
        seq_len(nrow(sim$group_true)),
        function(j) adjusted_rand_index(sim$group_true[j, ], align_obj$switched[j, ]),
        numeric(1)
      )
    } else {
      rep(NA_real_, nrow(sim$group_true))
    }

    compare_item <- data.frame(
      seed = seed,
      N = N,
      J = J,
      k = K,
      c_pen = c_pen_val,
      item = item_names,
      m_true = sim$m_true,
      selected_mj = selected_k_vec,
      selected_k = selected_k_vec,
      mj_difference = selected_k_vec - sim$m_true,
      mj_selection_type = classify_mj_selection(sim$m_true, selected_k_vec),
      is_mj_correct = as.integer(selected_k_vec == sim$m_true),
      is_mj_error = as.integer(selected_k_vec != sim$m_true),
      is_mj_overselect = as.integer(selected_k_vec > sim$m_true),
      is_mj_underselect = as.integer(selected_k_vec < sim$m_true),
      mse_beta_per_item = mse_per_item,
      mse_beta_step1_per_item = if (all(dim(sim_prob) == dim(polca_prob_step1))) {
        rowMeans((sim_prob - polca_prob_step1)^2)
      } else {
        rep(NA_real_, nrow(sim_prob))
      },
      group_recovery_per_item = align_obj$item_accuracy,
      partition_ARI_per_item = ari_per_item,
      true_group_pattern = apply(sim$group_true, 1, function(v) paste(v, collapse = "-")),
      selected_group_pattern_switched = apply(align_obj$switched, 1, function(v) paste(v, collapse = "-")),
      stringsAsFactors = FALSE
    )

    item_mj_selection_list[[c_idx]] <- compare_item[, c(
      "seed",
      "N",
      "J",
      "k",
      "c_pen",
      "item",
      "m_true",
      "selected_mj",
      "mj_difference",
      "mj_selection_type",
      "is_mj_correct",
      "is_mj_error",
      "is_mj_overselect",
      "is_mj_underselect"
    )]

    summary_rows[[c_idx]] <- data.frame(
      seed = seed,
      N = N,
      J = J,
      k = K,
      c_pen = c_pen_val,
      unrestricted_selected_k = best_k,
      unrestricted_selected_true_k = best_k == K,
      exact_mj_match_rate = mean(compare_item$m_true == compare_item$selected_k),
      mj_correct_count = sum(compare_item$is_mj_correct, na.rm = TRUE),
      mj_error_count = sum(compare_item$is_mj_error, na.rm = TRUE),
      mj_overselect_count = sum(compare_item$is_mj_overselect, na.rm = TRUE),
      mj_underselect_count = sum(compare_item$is_mj_underselect, na.rm = TRUE),
      MSE_beta_step1 = mse_beta_step1,
      MSE_nu_step1 = mse_nu_step1,
      MSE_beta = mse_beta,
      MSE_nu = mse_nu,
      Recovery_rate_switched = align_obj$overall_accuracy,
      ARI_partition_mean = safe_mean(compare_item$partition_ARI_per_item),
      sparse_logLik = sfit_sim$logLik,
      polca_logLik = best_fit$llik
    )

  }

  summary_row <- do.call(rbind, summary_rows)
  rownames(summary_row) <- NULL
  item_mj_selection <- do.call(rbind, item_mj_selection_list)
  rownames(item_mj_selection) <- NULL

  write.csv(
    summary_row,
    file = summary_file,
    row.names = FALSE
  )
  write.csv(
    unrestricted_fit_grid,
    file = unrestricted_grid_file,
    row.names = FALSE
  )
  write.csv(
    unrestricted_param_table,
    file = unrestricted_params_file,
    row.names = FALSE
  )
  write.csv(
    item_mj_selection,
    file = item_selection_file,
    row.names = FALSE
  )

  list(
    summary = summary_row,
    unrestricted_fit_grid = unrestricted_fit_grid,
    item_mj_selection = item_mj_selection
  )
}

# Main: parallel wrapper

run_simulation_reps_parallel <- function(
    N,
    n_reps = 100,
    dir = "output/simulations",
    workers = NULL,
    base_seed = 10000,
    ...
) {
  if (is.null(workers)) {
    workers <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
  }
  workers <- max(1L, min(as.integer(workers), n_reps))
  extra_args <- list(...)
  
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  seeds <- base_seed + seq_len(n_reps)

  cl <- parallel::makeCluster(workers, type = "PSOCK")
  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterEvalQ(cl, {
    source("script/functions.R")
    NULL
  })

  parallel::clusterExport(
    cl,
    varlist = c(
      "extract_polca_prob",
      "all_permutations",
      "best_label_switch",
      "align_group_with_fixed_classes",
      "best_prob_permutation",
      "normalize_refine_output",
      "permute_prob_columns",
      "adjusted_rand_index",
      "safe_mean",
      "classify_mj_selection",
      "summarize_mj_selection_by_item",
      "fit_polca_unrestricted_grid",
      "make_unrestricted_parameter_tables",
      "stack_unrestricted_tables",
      "run_one_simulation"
    ),
    envir = environment()
  )

  reps <- parallel::parLapply(cl, seeds, function(sd, n_val, extra_args) {
    do.call(
      run_one_simulation,
      c(list(seed = sd, N = n_val, dir = dir), extra_args)
    )
  }, n_val = N, extra_args = extra_args)

  summary_df <- do.call(
    rbind,
    lapply(seq_along(reps), function(i) {
      out_i <- reps[[i]]$summary
      out_i$rep_id <- i
      out_i
    })
  )
  rownames(summary_df) <- NULL
  summary_df$N <- N
  
  item_mj_selection <- do.call(
    rbind,
    lapply(seq_along(reps), function(i) {
      out_i <- reps[[i]]$item_mj_selection
      if (is.null(out_i) || nrow(out_i) == 0L) return(NULL)
      out_i$rep_id <- i
      out_i
    })
  )

  if (!is.null(item_mj_selection) && nrow(item_mj_selection) > 0L) {
    rownames(item_mj_selection) <- NULL
    item_mj_selection$N <- N
    item_mj_error_by_item <- summarize_mj_selection_by_item(item_mj_selection)
  } else {
    item_mj_selection <- data.frame()
    item_mj_error_by_item <- data.frame()
  }

  unrestricted_fit_grid <- do.call(
    rbind,
    lapply(seq_along(reps), function(i) {
      out_i <- reps[[i]]$unrestricted_fit_grid
      if (is.null(out_i) || nrow(out_i) == 0L) return(NULL)
      out_i$rep_id <- i
      out_i
    })
  )
  if (is.null(unrestricted_fit_grid)) {
    unrestricted_fit_grid <- data.frame()
  } else {
    rownames(unrestricted_fit_grid) <- NULL
  }

  list(
    N = N,
    seeds = seeds,
    summary = summary_df,
    unrestricted_fit_grid = unrestricted_fit_grid,
    item_mj_selection = item_mj_selection,
    item_mj_error_by_item = item_mj_error_by_item
  )
}

# Main: simulation driver

source("functions.R")
# ---- Requested simulation design ----
design <- data.frame(
  K = c(4),
  J = c(32)
)

N_grid <- c(500, 750, 1000,1500,2000)
c_pen_grid <- c(1, 5, 10, 20,40,80,160,320)
n_reps <- 100
workers <- 32
base_seed <- 1000L


dir <- "output2/resimulation_K4"
dir.create(dir, recursive = TRUE, showWarnings = FALSE)

all_settings_summary <- list()
all_settings_mj_error <- list()
all_settings_unrestricted <- list()

for (s in seq_len(nrow(design))) {
  K_val <- design$K[s]
  J_val <- design$J[s]
  
  # Draw one eta per (J, K) with concentration > 1 to avoid tiny class weights.
  eta_val <- sample_dirichlet(rep(8, K_val), seed = 20260402 + K_val)

  # Generate true item parameters once per (J, K), reused across all N and reps.
  true_params_val <- generate_true_item_parameters(
    J = J_val,
    K = K_val,
    low_range = c(0.18, 0.38),
    gap2_range = c(0.40, 0.55),
    seed = 20260402 + K_val + J_val
  )

  true_pi_df <- as.data.frame(true_params_val$Pi_true)
  true_pi_df$item <- rownames(true_params_val$Pi_true)
  true_pi_df$m_true <- true_params_val$m_true
  true_pi_df$group_pattern <- apply(true_params_val$group_true, 1, function(v) paste(v, collapse = "-"))
  true_pi_df <- true_pi_df[, c("item", "m_true", "group_pattern", paste0("class", seq_len(K_val)))]

  write.csv(
    true_pi_df,
    file = file.path(dir,sprintf("true_item_parameters_K%d_J%d.csv", K_val, J_val)),
    row.names = FALSE
  )

  eta_df <- data.frame(
    class = paste0("class", seq_len(K_val)),
    eta = eta_val,
    K = K_val,
    J = J_val
  )
  write.csv(
    eta_df,
    file = file.path(dir,sprintf("true_eta_K%d_J%d.csv", K_val, J_val)),
    row.names = FALSE
  )
  
  summary_by_N <- vector("list", length(N_grid))
  names(summary_by_N) <- paste0("N_", N_grid)
  mj_error_by_N <- vector("list", length(N_grid))
  names(mj_error_by_N) <- paste0("N_", N_grid)
  unrestricted_by_N <- vector("list", length(N_grid))
  names(unrestricted_by_N) <- paste0("N_", N_grid)

  for (i in seq_along(N_grid)) {
    N_val <- N_grid[i]

    run_out <- run_simulation_reps_parallel(
        N = N_val,
        n_reps = n_reps,
        workers = workers,
        dir = dir,
        base_seed = base_seed + 100000L * s + 1000L * i,
        J = J_val,
        K = K_val,
        eta = eta_val,
        true_params = true_params_val,
        class_grid = seq.int(max(2L, K_val - 1L), K_val + 1L),
        c_pen_grid = c_pen_grid,
        partition_search = "stepwise_merge",
        verbose = TRUE
      )

    summary_N <- run_out$summary
    summary_N$setting_id <- s
    rownames(summary_N) <- NULL
    summary_by_N[[i]] <- summary_N

    unrestricted_N <- run_out$unrestricted_fit_grid
    if (!is.null(unrestricted_N) && nrow(unrestricted_N) > 0L) {
      unrestricted_N$setting_id <- s
      rownames(unrestricted_N) <- NULL
      unrestricted_by_N[[i]] <- unrestricted_N
    } else {
      unrestricted_by_N[[i]] <- data.frame()
    }

    item_mj_error_N <- run_out$item_mj_error_by_item
    if (!is.null(item_mj_error_N) && nrow(item_mj_error_N) > 0L) {
      item_mj_error_N$setting_id <- s
      rownames(item_mj_error_N) <- NULL
      mj_error_by_N[[i]] <- item_mj_error_N
    } else {
      mj_error_by_N[[i]] <- data.frame()
    }
  }

  setting_summary_all <- do.call(rbind, summary_by_N)
  rownames(setting_summary_all) <- NULL
  all_settings_summary[[s]] <- setting_summary_all

  setting_mj_error_all <- do.call(rbind, mj_error_by_N)
  if (!is.null(setting_mj_error_all) && nrow(setting_mj_error_all) > 0L) {
    rownames(setting_mj_error_all) <- NULL
    all_settings_mj_error[[s]] <- setting_mj_error_all
  } else {
    all_settings_mj_error[[s]] <- data.frame()
  }

  setting_unrestricted_all <- do.call(rbind, unrestricted_by_N)
  if (!is.null(setting_unrestricted_all) && nrow(setting_unrestricted_all) > 0L) {
    rownames(setting_unrestricted_all) <- NULL
    all_settings_unrestricted[[s]] <- setting_unrestricted_all
  } else {
    all_settings_unrestricted[[s]] <- data.frame()
  }

}

summary_all_settings <- do.call(rbind, all_settings_summary)
rownames(summary_all_settings) <- NULL
write.csv(
  summary_all_settings,
  file = file.path(dir, "summary_all_settings.csv"),
  row.names = FALSE
)

summary_all_mj_errors <- do.call(rbind, all_settings_mj_error)
if (!is.null(summary_all_mj_errors) && nrow(summary_all_mj_errors) > 0L) {
  rownames(summary_all_mj_errors) <- NULL
  write.csv(
    summary_all_mj_errors,
    file = file.path(dir, "mj_selection_errors_all_settings.csv"),
    row.names = FALSE
  )
}

summary_all_unrestricted <- do.call(rbind, all_settings_unrestricted)
if (!is.null(summary_all_unrestricted) && nrow(summary_all_unrestricted) > 0L) {
  rownames(summary_all_unrestricted) <- NULL
  write.csv(
    summary_all_unrestricted,
    file = file.path(dir, "unrestricted_fit_grid_all_settings.csv"),
    row.names = FALSE
  )
}
