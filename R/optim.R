# =============================================================================
# optim.R
# Parameter optimisation wrappers around UKF_blend.
#
# Functions:
#   iterative_param_optim()  -- Repeated UKF passes until convergence
#   optim_params()           -- optim()-based (L-BFGS-B or SANN) outer optimiser
#
# Bugs fixed vs. original:
#   - Convergence norm changed from abs(sum(diff)) to true L2 norm
#     to prevent positive/negative cancellation hiding non-convergence.
#   - Non-trace path returns best parameters (lowest chi-sq), not last ones.
#   - Verbose log uses a temp file to avoid clobbering previous runs.
# =============================================================================

source("R/constants.R")
source("R/ukf_engine.R")


# -----------------------------------------------------------------------------
#' iterative_param_optim
#'
#' Repeatedly runs UKF_blend from the previous run's parameter estimates until
#' the L2 parameter change falls below param_tol or MAXSTEPS is reached.
#'
#' @param param_guess  Initial parameter vector (length N_p).
#' @param t_dummy      Scalar dummy time variable.
#' @param ts_data      Time-series matrix (T x (1+N_y)).
#' @param ode_model    ODE model function(t, x, p).
#' @param N_p          Number of model parameters.
#' @param N_y          Number of observed state variables.
#' @param dt           RK4 sub-step size.
#' @param dT           Data time step.
#' @param param_tol    L2 convergence tolerance (default 1e-3).
#' @param MAXSTEPS     Maximum iterations (default 1000).
#' @param R_scale      Observation noise scale.
#' @param Q_scale      Process noise scale.
#' @param forcePositive  Clamp parameters positive if TRUE.
#' @param seeded       Reproducible noise injection if TRUE.
#' @param verbose      Write per-iteration CSV log if TRUE.
#' @param log_file     Path for verbose CSV log (default: temp file).
#' @param param_lower  Optional numeric vector (length N_p) of per-parameter lower
#'                     bounds applied via pmax() after each UKF pass.  Use this to
#'                     keep a, b within the physiological BOLD frequency range.
#'                     NULL (default) applies no lower clipping beyond forcePositive.
#' @param param_upper  Optional numeric vector (length N_p) of per-parameter upper
#'                     bounds applied via pmin() after each UKF pass.
#'                     NULL (default) applies no upper clipping.
#' @param chisq_plateau_tol  Chi-square change below which iteration is declared a
#'                     plateau and stopped early.  Defaults to the static fallback in
#'                     UKF_CONSTANTS$CHISQ_PLATEAU_TOL (1e-5).  Pass the data-driven
#'                     value CHISQ_PLATEAU_TOL_DYN (= N_y * R_scale * 1e-3, computed
#'                     in §3.3) for a threshold calibrated to the actual noise floor.
#' @return List: par, value, param_norm, steps, param_est, xhat,
#'               chisq (vector over iterations), stop_reason (string).
#' @export
iterative_param_optim <- function(param_guess,
                                   t_dummy, ts_data, ode_model,
                                   N_p, N_y, dt, dT,
                                   param_tol    = UKF_CONSTANTS$PARAM_TOL_DEFAULT,
                                   MAXSTEPS     = UKF_CONSTANTS$MAXSTEPS_DEFAULT,
                                   R_scale      = 0.3,
                                   Q_scale      = 0.015,
                                   forcePositive = FALSE,
                                   seeded        = FALSE,
                                   verbose       = FALSE,
                                   log_file      = tempfile(fileext = ".csv"),
                                   param_lower   = NULL,
                                   param_upper   = NULL,
                                   chisq_plateau_tol = UKF_CONSTANTS$CHISQ_PLATEAU_TOL) {

  stopifnot(length(param_guess) == N_p, param_tol > 0, MAXSTEPS >= 1)
  if (!is.null(param_lower)) stopifnot(length(param_lower) == N_p)
  if (!is.null(param_upper)) stopifnot(length(param_upper) == N_p)

  best_param  <- param_guess
  best_chisq  <- Inf
  best_xhat   <- NULL
  chisq_trace <- numeric(0)
  param_norm  <- NA_real_
  steps       <- 0L
  done        <- FALSE

  if (verbose) {
    header <- c("step", paste0("param", seq_len(N_p)), "chisq", "param_norm")
    write.csv(as.data.frame(t(header)), log_file,
              row.names = FALSE, quote = FALSE)
  }

  while (!done) {
    ukf_run <- UKF_blend(t_dummy, ts_data, ode_model,
                          N_p, N_y, param_guess, dt, dT,
                          R_scale, Q_scale,
                          forcePositive = forcePositive,
                          seeded        = seeded)

    param_new    <- as.numeric(ukf_run$param_est)

    # ── Physiological bounds clipping ──────────────────────────────────────
    # Applied BEFORE storing as best_param and BEFORE using as next param_guess.
    # This keeps a, b in the BOLD frequency band across all iterations.
    if (!is.null(param_lower)) param_new <- pmax(param_new, param_lower)
    if (!is.null(param_upper)) param_new <- pmin(param_new, param_upper)

    chisq_trace  <- c(chisq_trace, ukf_run$chisq)

    # FIX: true L2 norm — prevents +/- cancellation masking non-convergence
    param_norm   <- sqrt(sum((param_new - param_guess)^2))

    if (ukf_run$chisq < best_chisq) {
      best_chisq <- ukf_run$chisq
      best_param <- param_new
      best_xhat  <- ukf_run$xhat
    }

    steps <- steps + 1L

    if (verbose) {
      row <- c(steps, param_new, ukf_run$chisq, param_norm)
      write.table(t(row), log_file,
                  sep = ",", row.names = FALSE, col.names = FALSE,
                  quote = FALSE, append = TRUE)
    }

    chisq_plateau <- length(chisq_trace) > 1 &&
      abs(diff(tail(chisq_trace, 2))) < chisq_plateau_tol

    done <- (param_norm < param_tol) || (steps >= MAXSTEPS) || chisq_plateau
    param_guess <- param_new
  }

  stop_reason <- if (param_norm < param_tol)          "param_tol"
                 else if (steps >= MAXSTEPS)           "maxsteps"
                 else                                  "chisq_plateau"

  list(
    par         = best_param,
    value       = best_chisq,
    param_norm  = param_norm,
    steps       = steps,
    param_est   = best_param,
    xhat        = best_xhat,
    chisq       = chisq_trace,    # chi-square at every iteration
    stop_reason = stop_reason
  )
}


# -----------------------------------------------------------------------------
#' optim_params
#'
#' Wraps optim() (L-BFGS-B or SANN) around UKF_blend to minimise chi-square
#' as an outer objective function.
#'
#' L-BFGS-B is the recommended method for the coupled-oscillator UKF because:
#'   - The chi-square landscape is smooth and locally convex near the optimum.
#'   - Three parameters (a, b, k) are few enough for quasi-Newton methods.
#'   - Hard box bounds on a, b, k are enforced exactly (unlike the soft clipping
#'     used in iterative_param_optim).
#'   - Convergence is by gradient norm — a well-defined, reproducible criterion.
#'
#' @param param_guess  Initial parameter vector (length N_p).  Use the single-pass
#'                     UKF estimate as a warm start for each pair.
#' @param method       "L-BFGS-B" (default) or "SANN".
#' @param lower_lim    Lower bound vector length N_p (L-BFGS-B only).
#' @param upper_lim    Upper bound vector length N_p (L-BFGS-B only).
#' @param maxit        Max optimizer iterations (default 200).
#' @param factr        L-BFGS-B relative function-decrease tolerance (default 1e3;
#'                     lower = tighter convergence; 1e7 is optim() default).
#' @param temp         Annealing temperature (SANN only, default 20).
#' @param t_dummy, ts_data, ode_model, N_p, N_y, dt, dT  See UKF_blend.
#' @param R_scale, Q_scale, forcePositive, seeded  See UKF_blend.
#' @return List: par, value, param_est, xhat, convergence (0=ok), message,
#'               counts (number of UKF evaluations).
#' @export
optim_params <- function(param_guess,
                          method    = "L-BFGS-B",
                          lower_lim = NULL, upper_lim = NULL,
                          maxit     = 200L,
                          factr     = 1e3,
                          temp      = 20,
                          t_dummy, ts_data, ode_model,
                          N_p, N_y, dt, dT,
                          R_scale = 0.3, Q_scale = 0.015,
                          forcePositive = FALSE, seeded = FALSE) {

  ukf_obj   <- NULL
  CHISQ_MAX <- 1e6   # penalty returned when UKF produces NaN/Inf

  chisq_objective <- function(par_vec) {
    result <- tryCatch(
      UKF_blend(t_dummy, ts_data, ode_model,
                N_p, N_y, par_vec, dt, dT,
                R_scale, Q_scale,
                forcePositive = forcePositive,
                seeded        = seeded),
      error = function(e) NULL
    )
    if (is.null(result) || !is.finite(result$chisq)) return(CHISQ_MAX)
    ukf_obj <<- result
    result$chisq
  }

  if (method == "SANN") {
    opt <- optim(param_guess, chisq_objective,
                 method  = "SANN",
                 control = list(maxit = maxit, temp = temp))
  } else {
    if (is.null(lower_lim) || is.null(upper_lim))
      stop("lower_lim and upper_lim are required for L-BFGS-B.")
    opt <- optim(param_guess, chisq_objective,
                 method  = "L-BFGS-B",
                 lower   = lower_lim,
                 upper   = upper_lim,
                 control = list(maxit = maxit, factr = factr))
  }

  list(
    par         = opt$par,
    value       = opt$value,
    param_est   = opt$par,
    xhat        = if (!is.null(ukf_obj)) ukf_obj$xhat else NULL,
    convergence = opt$convergence,    # 0 = converged; 1 = maxit hit; 52/53 = error
    message     = if (!is.null(opt$message)) opt$message else "",
    counts      = unname(opt$counts["function"])  # number of UKF evaluations
  )
}


# -----------------------------------------------------------------------------
#' print_optim_results
#'
#' Pretty-prints the output of iterative_param_optim.
#'
#' @param iter_opt  Return value of iterative_param_optim.
#' @param param_names  Character vector of parameter names (default: auto).
#' @export
print_optim_results <- function(iter_opt,
                                 param_names = NULL) {
  n <- length(iter_opt$param_est)
  if (is.null(param_names))
    param_names <- paste0("param", seq_len(n))

  cat("Iterative Parameter Optimization Results\n")
  cat(strrep("=", 45), "\n")
  cat("Estimated Parameters:\n")
  for (i in seq_len(n))
    cat(sprintf("  %-10s : %g\n", param_names[i], iter_opt$param_est[i]))
  cat(sprintf("\nChi-Square (best) : %g\n", iter_opt$value))
  cat(sprintf("Iterations        : %d\n",   iter_opt$steps))
  cat(sprintf("Parameter L2 norm : %g\n",   iter_opt$param_norm))
  cat(strrep("=", 45), "\n")
}
