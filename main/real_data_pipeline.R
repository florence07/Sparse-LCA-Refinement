source("script/functions.R")

if (!requireNamespace("poLCA", quietly = TRUE)) {
  stop("The poLCA package is required for script/real_data_pipeline.R.")
}

extract_polca_prob <- function(fit) {
  J <- length(fit$probs)
  K <- length(fit$P)
  out <- matrix(NA_real_, nrow = J, ncol = K)

  for (j in seq_len(J)) {
    prob_j <- as.matrix(fit$probs[[j]])
    if (ncol(prob_j) < 2L) {
      stop("This helper expects binary items coded as 1/2 inside poLCA.")
    }
    out[j, ] <- prob_j[, 2]
  }

  rownames(out) <- names(fit$probs)
  colnames(out) <- paste0("class", seq_len(K))
  out
}

data_path <- "promis/PROMIS 1 Wave 1.Rdata"
output_dir <- "output/real_data"
class_grid <- 7:9
penalties <- c(10, 20, 40)
default_penalty <- 20

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load(data_path)
if (!exists("table")) {
  stop("Expected object `table` after loading ", data_path)
}

item_names <- grep("^SRPPER", names(table), value = TRUE)
dat_raw <- table[, c("caseid", item_names)]
dat_raw <- dat_raw[stats::complete.cases(dat_raw), , drop = FALSE]

labels_df <- data.frame(
  item = item_names,
  label = vapply(item_names, function(x) attr(dat_raw[[x]], "label"), character(1)),
  stringsAsFactors = FALSE
)
write.csv(labels_df, file.path(output_dir, "labels_list.csv"), row.names = FALSE)

X <- dat_raw[item_names]
X[] <- lapply(X, function(x) ifelse(x < 4, 0L, 1L))
X <- as.data.frame(X, check.names = FALSE)
rm(table, dat_raw)

formula_lca <- as.formula(
  paste0("cbind(", paste(names(X), collapse = ", "), ") ~ 1")
)
X_polca <- X + 1L

fit_grid <- data.frame(
  k = class_grid,
  AIC = NA_real_,
  BIC = NA_real_,
  stringsAsFactors = FALSE
)

fits <- vector("list", length(class_grid))
for (i in seq_along(class_grid)) {
  k <- class_grid[[i]]
  fit <- poLCA::poLCA(
    formula = formula_lca,
    data = X_polca,
    nclass = k,
    nrep = 100,
    maxiter = 2000,
    calc.se = FALSE,
    verbose = FALSE
  )
  fits[[i]] <- fit
  fit_grid$AIC[i] <- fit$aic
  fit_grid$BIC[i] <- fit$bic
}

write.csv(fit_grid, file.path(output_dir, "unrestricted_fit_grid.csv"), row.names = FALSE)

best_idx <- which.min(fit_grid$BIC)
best_fit <- fits[[best_idx]]
K <- class_grid[[best_idx]]

gamma_hat <- best_fit$posterior
initial_beta_raw <- extract_polca_prob(best_fit)
initial_order <- order(colMeans(initial_beta_raw))
initial_beta <- initial_beta_raw[, initial_order, drop = FALSE]
initial_eta <- best_fit$P[initial_order]
initial_class_avg <- colMeans(initial_beta)

write.csv(initial_eta, file.path(output_dir, "initial_eta.csv"), row.names = FALSE)
write.csv(initial_beta, file.path(output_dir, "initial_beta.csv"), row.names = TRUE)
write.csv(initial_class_avg, file.path(output_dir, "initial_class_avg.csv"), row.names = FALSE)

refinement <- refine_lca_binary(
  X = X,
  gamma = gamma_hat,
  k_grid = seq_len(K),
  c_pen_grid = penalties,
  tol = 1e-10,
  verbose = TRUE
)

write.csv(refinement$selected_k, file.path(output_dir, "selected_k.csv"), row.names = TRUE)

for (penalty in penalties) {
  penalty_key <- sprintf("c_pen=%s", penalty)
  group_id <- refinement$group_id[, , penalty_key]
  beta_refined <- refinement$beta_refined[, , penalty_key]

  sparse_fit <- reestimate_sparse_lca(
    X = X,
    group_id = group_id,
    eta_init = best_fit$P,
    beta_init = beta_refined,
    gamma_init = gamma_hat,
    maxit = 300,
    tol = 1e-8,
    verbose = TRUE
  )

  class_order <- order(colMeans(sparse_fit$beta))
  beta_out <- sparse_fit$beta[, class_order, drop = FALSE]
  group_out <- group_id[, class_order, drop = FALSE]
  eta_out <- sparse_fit$eta[class_order]
  class_avg_out <- colMeans(beta_out)

  write.csv(eta_out, file.path(output_dir, sprintf("eta_%s.csv", penalty)), row.names = FALSE)
  write.csv(beta_out, file.path(output_dir, sprintf("fin_beta_%s.csv", penalty)), row.names = TRUE)
  write.csv(group_out, file.path(output_dir, sprintf("fin_group_%s.csv", penalty)), row.names = TRUE)
  write.csv(
    class_avg_out,
    file.path(output_dir, sprintf("final_class_avg_%s.csv", penalty)),
    row.names = FALSE
  )

  if (penalty == default_penalty) {
    write.csv(eta_out, file.path(output_dir, "eta.csv"), row.names = FALSE)
    write.csv(beta_out, file.path(output_dir, "fin_beta.csv"), row.names = TRUE)
    write.csv(group_out, file.path(output_dir, "fin_group.csv"), row.names = TRUE)
  }
}

message(
  "Saved real-data outputs for K = ", K,
  " and penalties ", paste(penalties, collapse = ", "),
  "."
)
