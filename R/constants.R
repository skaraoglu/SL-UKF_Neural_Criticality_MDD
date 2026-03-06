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
  CHISQ_PLATEAU_TOL  = 1e-8,   # Chi-square change below which we declare plateau

  # --- RK4 integration ------------------------------------------------------
  DT_FRACTION        = 0.1,    # dt = DT_FRACTION * dT  (sub-step size)
  STIFFNESS_WARN     = 100     # Warn if any |parameter| exceeds this value
)
