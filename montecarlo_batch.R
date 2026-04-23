# MONTE CARLO SIMULATION - DPD (batch grid)
# Adapted to run all combinations:
# sparse in {0.0, 0.1, ..., 1.0}
# N in {20, 100}
# T0 in {5, 10, 20, 30}
# gamma in {0.2, 0.8}
# num_iterations fixed at 1000

rm(list = ls())
options(scipen = 999)

suppressPackageStartupMessages({
  library(plm)
  library(dplyr)
  library(glmnet)
  library(fastDummies)
  library(FEShR)
})

# -----------------------------
# Global configuration
# -----------------------------
num_iterations <- 1000
rho <- 0.5
output_dir <- "outputs"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sparse_grid <- seq(0, 1, by = 0.1)
N_grid <- c(20, 100)
T0_grid <- c(5, 10, 20, 30)
gamma_grid <- c(0.2, 0.8)

# -----------------------------
# Helper functions
# -----------------------------
safe_colmean <- function(M) {
  if (is.null(dim(M))) return(mean(M, na.rm = TRUE))
  apply(M, 2, function(x) mean(x, na.rm = TRUE))
}

make_year_blocks <- function(year, K = 5) {
  yrs <- sort(unique(year))
  K_eff <- min(K, length(yrs))
  split(yrs, cut(seq_along(yrs), breaks = K_eff, labels = FALSE))
}

make_lambda_seq <- function(X, y, alpha, penalty_factor) {
  fit0 <- glmnet(
    X, y,
    alpha = alpha,
    penalty.factor = penalty_factor,
    intercept = FALSE,
    standardize = TRUE
  )
  fit0$lambda
}

rolling_origin_select_lambda <- function(X, y, years, year_blocks, alpha,
                                         penalty_factor, lambda_seq, use_1se = TRUE) {
  K <- length(year_blocks)
  if (K < 2) return(lambda_seq[which.min(lambda_seq)])

  val_blocks <- 2:K
  L <- matrix(NA_real_, nrow = length(lambda_seq), ncol = length(val_blocks))

  for (jj in seq_along(val_blocks)) {
    k <- val_blocks[jj]
    tr <- years %in% unlist(year_blocks[1:(k - 1)])
    va <- years %in% unlist(year_blocks[k])

    fit <- glmnet(
      X[tr, , drop = FALSE], y[tr],
      alpha = alpha,
      lambda = lambda_seq,
      penalty.factor = penalty_factor,
      intercept = FALSE,
      standardize = TRUE
    )

    pred <- predict(fit, newx = X[va, , drop = FALSE], s = lambda_seq)
    L[, jj] <- colMeans((y[va] - pred)^2)
  }

  mean_loss <- rowMeans(L, na.rm = TRUE)
  se_loss <- apply(L, 1, sd, na.rm = TRUE) / sqrt(ncol(L))
  j_min <- which.min(mean_loss)

  if (!use_1se) return(lambda_seq[j_min])

  thresh <- mean_loss[j_min] + se_loss[j_min]
  max(lambda_seq[mean_loss <= thresh])
}

eval_metrics <- function(y, yhat) {
  y <- as.numeric(y)
  yhat <- as.numeric(yhat)
  ok <- is.finite(y) & is.finite(yhat)
  if (sum(ok) == 0) return(c(bias = NA_real_, var = NA_real_, mse = NA_real_))

  err <- y[ok] - yhat[ok]
  c(
    bias = mean(err),
    var = var(yhat[ok]),
    mse = mean(err^2)
  )
}

compute_shrunk_lsdv <- function(in.sample.df, shrink_type = c("EBMLE", "URE"),
                                centering = c("gen", "0")) {
  shrink_type <- match.arg(shrink_type)
  centering <- match.arg(centering)

  in.sample.df$id <- as.factor(in.sample.df$id)

  fit <- lm(Y ~ lagY + X + id, data = in.sample.df)

  ids <- levels(in.sample.df$id)
  coefs <- coef(fit)
  V <- vcov(fit)

  alpha_hat <- numeric(length(ids))
  names(alpha_hat) <- ids
  alpha_hat[ids[1]] <- coefs["(Intercept)"]

  if (length(ids) >= 2) {
    for (j in 2:length(ids)) {
      cname <- paste0("id", ids[j])
      alpha_hat[ids[j]] <- coefs["(Intercept)"] + coefs[cname]
    }
  }

  M_list <- vector("list", length(ids))
  names(M_list) <- ids
  M_list[[ids[1]]] <- matrix(V["(Intercept)", "(Intercept)"], 1, 1)

  if (length(ids) >= 2) {
    for (j in 2:length(ids)) {
      cname <- paste0("id", ids[j])
      var_alpha_j <- V["(Intercept)", "(Intercept)"] + V[cname, cname] +
        2 * V["(Intercept)", cname]
      M_list[[ids[j]]] <- matrix(var_alpha_j, 1, 1)
    }
  }

  y_mat <- matrix(alpha_hat, nrow = 1)

  shrink_fit <- FEShR::fe_shrink(
    y = y_mat,
    M = M_list,
    centering = centering,
    type = shrink_type,
    n_init_vals = 1,
    all_init_vals = FALSE,
    optim_control = list(maxit = 500)
  )

  alpha_shrunk <- as.numeric(shrink_fit$thetahat)
  names(alpha_shrunk) <- ids

  list(
    fit_lm = fit,
    alpha_fe = alpha_hat,
    alpha_shrunk = alpha_shrunk,
    shrink_fit = shrink_fit
  )
}

predict_gmm_levels_1step <- function(gmm_obj, df) {
  g <- as.numeric(coef(gmm_obj)["lag(Y, 1)"])
  b <- as.numeric(coef(gmm_obj)["X"])

  dYhat <- g * as.numeric(df$LDY) + b * as.numeric(df$DX)
  as.numeric(df$lagY) + dYhat
}

format_num <- function(x, digits = 5) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

df_to_latex <- function(df, caption, label, digits = 5) {
  cols <- colnames(df)
  align <- paste0("l", paste(rep("r", ncol(df) - 1), collapse = ""))

  body <- apply(df, 1, function(row) {
    vals <- mapply(function(v, j) {
      if (j <= 4) as.character(v) else format_num(as.numeric(v), digits = digits)
    }, row, seq_along(row), SIMPLIFY = TRUE)
    paste(vals, collapse = " & ")
  })

  header <- paste(cols, collapse = " & ")

  paste0(
    "\\begin{table}[!htbp]\\centering\n",
    "\\caption{", caption, "}\n",
    "\\label{", label, "}\n",
    "\\begin{tabular}{", align, "}\n",
    "\\hline\n",
    header, " \\\\\n",
    "\\hline\n",
    paste(body, collapse = " \\\\\n"), " \\\\\n",
    "\\hline\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
}

simulate_one_config <- function(sparse, N, T0, gamma, num_iterations, rho = 0.5) {
  T <- T0 + 10
  beta <- 1 - gamma
  Ttr <- floor(T0 * 0.8)

  panel <- data.frame(
    id = rep(1:N, each = T),
    time = rep(1:T, times = N),
    X = NA_real_
  )

  set.seed(123)
  panel$X[panel$time == 1] <- 0
  for (i in 1:N) {
    for (t in 2:T) {
      panel$X[(panel$id == i) & (panel$time == t)] <-
        rho * panel$X[(panel$id == i) & (panel$time == t - 1)] +
        rnorm(1, mean = 0, sd = 2)
    }
  }

  panel <- panel %>% arrange(id, time)
  panel$FE <- rep(rnorm(N, mean = 0, sd = 1), each = T)

  num_of_zeros <- round(sparse * N)
  if (num_of_zeros > 0) {
    indices_to_replace <- sample(N, num_of_zeros)
    panel$FE[panel$id %in% indices_to_replace] <- 0
  }

  models <- c("OLS", "LSDV", "EBMLE", "URE", "AH", "GMM1", "GMM2", "LASSO", "RIDGE", "ELASTIC_NET")

  gamma_hat <- matrix(NA_real_, nrow = num_iterations, ncol = length(models), dimnames = list(NULL, models))
  in_bias <- in_var <- in_mse <- matrix(NA_real_, nrow = num_iterations, ncol = length(models), dimnames = list(NULL, models))
  out_bias <- out_var <- out_mse <- matrix(NA_real_, nrow = num_iterations, ncol = length(models), dimnames = list(NULL, models))

  for (it in 1:num_iterations) {
    set.seed(123 + it)

    panel_new <- panel
    panel_new$Y <- NA_real_
    panel_new$Y[panel_new$time == 1] <- 0
    panel_new$e <- rnorm(N * T, mean = 0, sd = 1)

    for (j in 1:N) {
      for (t in 2:T) {
        idx <- (panel_new$id == j) & (panel_new$time == t)
        idx_lag <- (panel_new$id == j) & (panel_new$time == t - 1)
        panel_new$Y[idx] <-
          gamma * panel_new$Y[idx_lag] +
          beta * panel_new$X[idx] +
          panel_new$FE[idx] +
          panel_new$e[idx]
      }
    }

    panel_new <- panel_new %>%
      group_by(id) %>%
      mutate(
        lagY = dplyr::lag(Y, 1),
        DY = Y - dplyr::lag(Y, 1),
        L2Y = dplyr::lag(lagY, 1),
        LDY = dplyr::lag(DY, 1),
        lagX = dplyr::lag(X, 1),
        DX = X - lagX
      ) %>%
      ungroup()

    panel_new <- as.data.frame(panel_new) %>% filter(time > 10)
    panel_new$time <- panel_new$time - 10

    panel_new <- fastDummies::dummy_cols(
      panel_new,
      select_columns = "id",
      remove_first_dummy = FALSE,
      remove_selected_columns = FALSE
    )

    in.sample.df <- panel_new %>% filter(time <= Ttr)
    out.sample.df <- panel_new %>% filter(time > Ttr)

    in.sample <- pdata.frame(in.sample.df, index = c("id", "time"))
    out.of.sample <- pdata.frame(out.sample.df, index = c("id", "time"))

    y_in <- as.numeric(in.sample$Y)
    y_out <- as.numeric(out.of.sample$Y)

    X_in <- as.matrix(in.sample %>% select(lagY, X, starts_with("id_")))
    X_out <- as.matrix(out.of.sample %>% select(lagY, X, starts_with("id_")))

    penalty_factors <- c(0, 0, rep(1, N))

    OLS <- LSDV <- EBMLE <- URE <- AH <- GMM1 <- GMM2 <- NULL

    try({
      OLS <- plm(Y ~ lagY + X, data = in.sample, model = "pooling")
      LSDV <- plm(Y ~ lagY + X, data = in.sample, model = "within", effect = "individual")
      EBMLE <- compute_shrunk_lsdv(in.sample.df, shrink_type = "EBMLE", centering = "gen")
      URE <- compute_shrunk_lsdv(in.sample.df, shrink_type = "URE", centering = "gen")
      AH <- plm(DY ~ LDY - 1 + DX | L2Y + DX, data = in.sample, model = "pooling")
    }, silent = TRUE)

    try({
      GMM1 <- pgmm(
        Y ~ lag(Y, 1) + X | lag(Y, 2:(Ttr - 2)),
        data = in.sample,
        effect = "individual",
        model = "onestep",
        transformation = "d",
        collapse = FALSE
      )
    }, silent = TRUE)

    try({
      GMM2 <- pgmm(
        Y ~ lag(Y, 1) + X | lag(Y, 2:(Ttr - 2)),
        data = in.sample,
        effect = "individual",
        model = "twosteps",
        transformation = "d",
        collapse = FALSE
      )
    }, silent = TRUE)

    years_vec <- as.numeric(attr(in.sample, "index")$time)
    blocks <- make_year_blocks(years_vec, K = 5)

    l_seq_lasso <- make_lambda_seq(X_in, y_in, alpha = 1, penalty_factor = penalty_factors)
    best_lambda_lasso <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
      alpha = 1, penalty_factor = penalty_factors, lambda_seq = l_seq_lasso, use_1se = TRUE
    )

    l_seq_ridge <- make_lambda_seq(X_in, y_in, alpha = 0, penalty_factor = penalty_factors)
    best_lambda_ridge <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
      alpha = 0, penalty_factor = penalty_factors, lambda_seq = l_seq_ridge, use_1se = TRUE
    )

    alpha_en <- 0.5
    l_seq_en <- make_lambda_seq(X_in, y_in, alpha = alpha_en, penalty_factor = penalty_factors)
    best_lambda_en <- rolling_origin_select_lambda(X_in, y_in, years_vec, blocks,
      alpha = alpha_en, penalty_factor = penalty_factors, lambda_seq = l_seq_en, use_1se = TRUE
    )

    LASSO <- glmnet(X_in, y_in, alpha = 1, lambda = best_lambda_lasso, penalty.factor = penalty_factors, intercept = FALSE)
    RIDGE <- glmnet(X_in, y_in, alpha = 0, lambda = best_lambda_ridge, penalty.factor = penalty_factors, intercept = FALSE)
    ELASTIC_NET <- glmnet(X_in, y_in, alpha = alpha_en, lambda = best_lambda_en, penalty.factor = penalty_factors, intercept = FALSE)

    if (!is.null(OLS)) gamma_hat[it, "OLS"] <- coef(OLS)["lagY"]
    if (!is.null(LSDV)) gamma_hat[it, "LSDV"] <- coef(LSDV)["lagY"]
    if (!is.null(AH)) gamma_hat[it, "AH"] <- coef(AH)["LDY"]
    if (!is.null(GMM1)) gamma_hat[it, "GMM1"] <- coef(GMM1)["lag(Y, 1)"]
    if (!is.null(GMM2)) gamma_hat[it, "GMM2"] <- coef(GMM2)["lag(Y, 1)"]
    if (!is.null(EBMLE)) gamma_hat[it, "EBMLE"] <- coef(EBMLE$fit_lm)["lagY"]
    if (!is.null(URE)) gamma_hat[it, "URE"] <- coef(URE$fit_lm)["lagY"]
    gamma_hat[it, "LASSO"] <- as.numeric(coef(LASSO)["lagY", ])
    gamma_hat[it, "RIDGE"] <- as.numeric(coef(RIDGE)["lagY", ])
    gamma_hat[it, "ELASTIC_NET"] <- as.numeric(coef(ELASTIC_NET)["lagY", ])

    metric_update <- function(name, y, yhat, container_bias, container_var, container_mse, idx) {
      m <- eval_metrics(y, yhat)
      container_bias[idx, name] <- m["bias"]
      container_var[idx, name] <- m["var"]
      container_mse[idx, name] <- m["mse"]
      list(bias = container_bias, var = container_var, mse = container_mse)
    }

    if (!is.null(OLS)) {
      yhat <- as.numeric(coef(OLS)["(Intercept)"] + coef(OLS)["lagY"] * in.sample$lagY + coef(OLS)["X"] * in.sample$X)
      tmp <- metric_update("OLS", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    if (!is.null(LSDV)) {
      fe_i <- fixef(LSDV, effect = "individual")
      yhat <- as.numeric(fe_i[as.character(in.sample$id)] + coef(LSDV)["lagY"] * in.sample$lagY + coef(LSDV)["X"] * in.sample$X)
      tmp <- metric_update("LSDV", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    if (!is.null(EBMLE)) {
      yhat <- as.numeric(EBMLE$alpha_shrunk[as.character(in.sample.df$id)] + coef(EBMLE$fit_lm)["lagY"] * in.sample.df$lagY + coef(EBMLE$fit_lm)["X"] * in.sample.df$X)
      tmp <- metric_update("EBMLE", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    if (!is.null(URE)) {
      yhat <- as.numeric(URE$alpha_shrunk[as.character(in.sample.df$id)] + coef(URE$fit_lm)["lagY"] * in.sample.df$lagY + coef(URE$fit_lm)["X"] * in.sample.df$X)
      tmp <- metric_update("URE", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    if (!is.null(AH)) {
      yhat <- as.numeric(fitted(AH)) + as.numeric(in.sample$lagY)
      tmp <- metric_update("AH", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    if (!is.null(GMM1)) {
      yhat <- predict_gmm_levels_1step(GMM1, in.sample)
      tmp <- metric_update("GMM1", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    if (!is.null(GMM2)) {
      yhat <- predict_gmm_levels_1step(GMM2, in.sample)
      tmp <- metric_update("GMM2", y_in, yhat, in_bias, in_var, in_mse, it)
      in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    }
    tmp <- metric_update("LASSO", y_in, as.numeric(predict(LASSO, newx = X_in, s = best_lambda_lasso)), in_bias, in_var, in_mse, it)
    in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    tmp <- metric_update("RIDGE", y_in, as.numeric(predict(RIDGE, newx = X_in, s = best_lambda_ridge)), in_bias, in_var, in_mse, it)
    in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse
    tmp <- metric_update("ELASTIC_NET", y_in, as.numeric(predict(ELASTIC_NET, newx = X_in, s = best_lambda_en)), in_bias, in_var, in_mse, it)
    in_bias <- tmp$bias; in_var <- tmp$var; in_mse <- tmp$mse

    if (!is.null(OLS)) {
      yhat <- as.numeric(coef(OLS)["(Intercept)"] + coef(OLS)["lagY"] * out.of.sample$lagY + coef(OLS)["X"] * out.of.sample$X)
      tmp <- metric_update("OLS", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    if (!is.null(LSDV)) {
      fe_i <- fixef(LSDV, effect = "individual")
      yhat <- as.numeric(fe_i[as.character(out.of.sample$id)] + coef(LSDV)["lagY"] * out.of.sample$lagY + coef(LSDV)["X"] * out.of.sample$X)
      tmp <- metric_update("LSDV", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    if (!is.null(EBMLE)) {
      yhat <- as.numeric(EBMLE$alpha_shrunk[as.character(out.sample.df$id)] + coef(EBMLE$fit_lm)["lagY"] * out.sample.df$lagY + coef(EBMLE$fit_lm)["X"] * out.sample.df$X)
      tmp <- metric_update("EBMLE", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    if (!is.null(URE)) {
      yhat <- as.numeric(URE$alpha_shrunk[as.character(out.sample.df$id)] + coef(URE$fit_lm)["lagY"] * out.sample.df$lagY + coef(URE$fit_lm)["X"] * out.sample.df$X)
      tmp <- metric_update("URE", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    if (!is.null(AH)) {
      yhat <- as.numeric(coef(AH)["LDY"] * out.of.sample$LDY + coef(AH)["DX"] * out.of.sample$DX + out.of.sample$lagY)
      tmp <- metric_update("AH", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    if (!is.null(GMM1)) {
      yhat <- predict_gmm_levels_1step(GMM1, out.of.sample)
      tmp <- metric_update("GMM1", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    if (!is.null(GMM2)) {
      yhat <- predict_gmm_levels_1step(GMM2, out.of.sample)
      tmp <- metric_update("GMM2", y_out, yhat, out_bias, out_var, out_mse, it)
      out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    }
    tmp <- metric_update("LASSO", y_out, as.numeric(predict(LASSO, newx = X_out, s = best_lambda_lasso)), out_bias, out_var, out_mse, it)
    out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    tmp <- metric_update("RIDGE", y_out, as.numeric(predict(RIDGE, newx = X_out, s = best_lambda_ridge)), out_bias, out_var, out_mse, it)
    out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
    tmp <- metric_update("ELASTIC_NET", y_out, as.numeric(predict(ELASTIC_NET, newx = X_out, s = best_lambda_en)), out_bias, out_var, out_mse, it)
    out_bias <- tmp$bias; out_var <- tmp$var; out_mse <- tmp$mse
  }

  gamma_bias <- gamma - safe_colmean(gamma_hat)
  gamma_var <- apply(gamma_hat, 2, function(x) mean((x - mean(x, na.rm = TRUE))^2, na.rm = TRUE))
  gamma_mse <- apply(gamma_hat, 2, function(x) mean((x - gamma)^2, na.rm = TRUE))

  table1 <- rbind(Bias = gamma_bias, Var = gamma_var, MSE = gamma_mse)
  table2 <- rbind(Bias = safe_colmean(in_bias), Var = safe_colmean(in_var), MSE = safe_colmean(in_mse))
  table3 <- rbind(Bias = safe_colmean(out_bias), Var = safe_colmean(out_var), MSE = safe_colmean(out_mse))

  list(table1 = table1, table2 = table2, table3 = table3)
}

bind_by_scenario <- function(result_map, table_name) {
  rows <- list()
  idx <- 1
  for (scenario_name in names(result_map)) {
    tab <- result_map[[scenario_name]][[table_name]]
    df <- as.data.frame(round(tab, 5))
    df$Measure <- rownames(df)
    df$Scenario <- scenario_name
    rownames(df) <- NULL
    df <- df[, c("Scenario", "Measure", colnames(tab))]
    rows[[idx]] <- df
    idx <- idx + 1
  }
  do.call(rbind, rows)
}

# -----------------------------
# Batch run and output
# -----------------------------
for (sparse in sparse_grid) {
  message("Running sparse = ", sparse)
  scenario_results <- list()

  for (N in N_grid) {
    for (T0 in T0_grid) {
      for (gamma in gamma_grid) {
        scenario_name <- paste0("N", N, "_T", T0, "_g", gamma)
        message("  Scenario: ", scenario_name)
        scenario_results[[scenario_name]] <- simulate_one_config(
          sparse = sparse,
          N = N,
          T0 = T0,
          gamma = gamma,
          num_iterations = num_iterations,
          rho = rho
        )
      }
    }
  }

  table1_df <- bind_by_scenario(scenario_results, "table1")
  table2_df <- bind_by_scenario(scenario_results, "table2")
  table3_df <- bind_by_scenario(scenario_results, "table3")

  sparse_tag <- gsub("\\.", "p", formatC(sparse, digits = 1, format = "f"))

  write.csv(table1_df, file.path(output_dir, paste0("table1_sparse_", sparse_tag, ".csv")), row.names = FALSE)
  write.csv(table2_df, file.path(output_dir, paste0("table2_sparse_", sparse_tag, ".csv")), row.names = FALSE)
  write.csv(table3_df, file.path(output_dir, paste0("table3_sparse_", sparse_tag, ".csv")), row.names = FALSE)

  tex1 <- df_to_latex(table1_df,
    caption = paste0("Gamma estimation metrics (sparse = ", sparse, ")"),
    label = paste0("tab:table1_sparse_", sparse_tag)
  )
  tex2 <- df_to_latex(table2_df,
    caption = paste0("In-sample prediction metrics (sparse = ", sparse, ")"),
    label = paste0("tab:table2_sparse_", sparse_tag)
  )
  tex3 <- df_to_latex(table3_df,
    caption = paste0("Out-of-sample prediction metrics (sparse = ", sparse, ")"),
    label = paste0("tab:table3_sparse_", sparse_tag)
  )

  writeLines(tex1, file.path(output_dir, paste0("table1_sparse_", sparse_tag, ".tex")))
  writeLines(tex2, file.path(output_dir, paste0("table2_sparse_", sparse_tag, ".tex")))
  writeLines(tex3, file.path(output_dir, paste0("table3_sparse_", sparse_tag, ".tex")))
}

message("Finished. Outputs saved in: ", normalizePath(output_dir))
