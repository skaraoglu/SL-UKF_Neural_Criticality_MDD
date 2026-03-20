# =============================================================================
# constants.R
# Central definition of all tuning constants used across the UKF pipeline.
# Edit values here; they propagate everywhere automatically.
# =============================================================================

UKF_CONSTANTS <- list(

  # --- Numerical stability (Cholesky / covariance conditioning) -------------
  JITTER_INIT  = 1e-8,   # First jitter added to Pxx/Pyy before Cholesky
  JITTER_MAX   = 1e-2,   # Maximum jitter before falling back to nearPD
  COND_NUM_MAX = 1e12,   # Maximum acceptable condition number for Pyy
  EIGVAL_MIN   = 1e-12,  # Minimum eigenvalue for positive-definiteness check

  # --- Parameter constraints ------------------------------------------------
  PARAM_MIN    = 1e-8,   # Smallest allowed value when forcePositive = TRUE

  # --- Physiological bounds for a, b (rad²/TR², TR = 2 s) ------------------
  # Oscillator natural frequency: f = sqrt(a) / (2π TR)
  # Lower: BOLD band lower edge  f = 0.01 Hz → a = (2π × 0.01 × 2)² ≈ 0.0158
  # Upper: 1.5 × BOLD upper edge f = 0.15 Hz → a = (2π × 0.15 × 2)² ≈ 3.553
  #        (50% slack above 0.10 Hz avoids over-constraining near-boundary pairs)
  # k is not bounded above — inter-regional coupling has no physiological ceiling.
  BOLD_AB_LOWER = (2 * pi * 0.01 * 2.0)^2,   # ≈ 0.0158 rad²/TR²
  BOLD_AB_UPPER = (2 * pi * 0.15 * 2.0)^2,   # ≈ 3.553 rad²/TR²

  # --- Iterative optimisation -----------------------------------------------
  PARAM_TOL_DEFAULT  = 1e-3,   # Default L2 convergence tolerance
  MAXSTEPS_DEFAULT   = 1000,   # Default maximum iterations

  # Chi-square plateau tolerance — static fallback value.
  #
  # History / bug fix:
  #   The original value of 1e-8 was ~8 orders of magnitude below the minimum
  #   achievable chi-square for this dataset (floor ≈ N_y × R_scale ≈ 0.17).
  #   It could never physically be reached, so the plateau criterion never fired
  #   and every pair ran until MAXSTEPS was exhausted.
  #
  # Static fix: 1e-5  (~0.006% of the chi-square floor — fires reliably)
  #
  # Dynamic version (preferred): computed in §3.3 of main_analysis.ipynb from
  #   the actual R_scale and N_y of the current dataset:
  #     CHISQ_PLATEAU_TOL_DYN = N_y * R_scale * 1e-3
  #   ≈ 2 × 0.085 × 0.001 = 1.7×10⁻⁴  (~5× the observed per-step improvement)
  #   This value is passed explicitly to iterative_param_optim() via the
  #   chisq_plateau_tol argument; the value here serves as the fallback when
  #   no dynamic value has been computed.
  CHISQ_PLATEAU_TOL  = 1e-5,   # static fallback — overridden by dynamic §3.3 value

  # --- RK4 integration ------------------------------------------------------
  DT_FRACTION        = 0.1,    # dt = DT_FRACTION * dT  (sub-step size)
  STIFFNESS_WARN     = 100     # Warn if any |parameter| exceeds this value
)
