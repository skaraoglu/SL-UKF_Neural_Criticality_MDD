# =============================================================================
# models.R
# Coupled oscillator ODE models used as the dynamical prior for UKF parameter
# estimation on fMRI resting-state time series.
#
# All models follow the interface:
#   ode_model(t, x, p)  ->  matrix of first derivatives (same dims as x)
# where:
#   t  : scalar dummy time (models have no explicit time dependence)
#   x  : (N_y x N_sigma) matrix of state variables (or (N_y x 1) vector)
#   p  : (N_p x N_sigma) matrix of parameters     (or (N_p x 1) vector)
#
# ── INTEGRATION-ORDER FIX (2026-03) ─────────────────────────────────────────
# The original models returned [theta1_ddot, theta2_ddot] into a state vector
# [theta1, theta2].  Because propagate_model() feeds model output directly into
# the RK4 accumulator as dy/dt, this caused the integrator to solve the gradient-
# flow system  dtheta/dt = -(gamma+k)*sin(theta)  rather than the second-order
# pendulum.  That system has time constant tau = 1/(gamma+k) ~= 0.22 TRs and
# collapses to theta = 0 within one TR, giving r2_ode ~= 0 for the remaining
# 259 TRs.  Confirmed by calibration experiment 2026-03-08 (r2_ode <= 0.029
# across all pairs and both models).
#
# The fix: extend the state to include angular velocities.
#   Old state (N_y = 2): [theta1,          theta2         ]
#   New state (N_y = 4): [theta1, theta1d, theta2, theta2d]
#
# The model now returns d/dt[theta1, theta1d, theta2, theta2d]
#                          = [theta1d, theta1dd, theta2d, theta2dd].
# The acceleration rows (theta1dd, theta2dd) are unchanged from original physics.
# ukf_engine.R, optim.R, and preprocessing.R are unaffected.
#
# Observation: only positions [theta1, theta2] are observed (rows 1 and 3).
# Set OBSERVED_ROWS <- c(1, 3) and N_Y_OBS <- 2 in main_analysis.ipynb.
# N_Y (state dim) = 4; N_Y_OBS (observation dim) = 2.
# =============================================================================


# -----------------------------------------------------------------------------
# coupled_osc_model_glk_2nd
#
# Second-order symmetric coupled oscillator with optional linear damping.
# State is [theta1, theta1d, theta2, theta2d]  (N_y = 4).
# Only positions [theta1, theta2] (rows 1 and 3) are observed via BOLD.
#
# EQUATIONS OF MOTION:
#   theta1_dd = -gamma*sin(theta1) - 2*zeta*sqrt(gamma)*theta1d - k*(sin(theta1)-sin(theta2))
#   theta2_dd = -gamma*sin(theta2) - 2*zeta*sqrt(gamma)*theta2d + k*(sin(theta1)-sin(theta2))
#
# Setting zeta = 0 recovers the conservative (undamped) pendulum.
# Setting zeta > 0 adds a velocity-proportional restoring force causing free
# oscillations to decay with envelope time constant tau = 1/(zeta*sqrt(gamma)) TRs.
#
# PARAMETERS (p rows):
#   p[1,] = gamma : shared squared natural frequency (rad^2/TR^2)
#                   INITIALISE per pair from spectral peak: (2*pi*f_dominant)^2
#                   Bounds: [BOLD_AB_LOWER, BOLD_AB_UPPER] = [0.0158, 3.553]
#   p[2,] = zeta  : damping ratio (dimensionless, >= 0)
#                   0       = conservative pendulum (oscillates forever)
#                   (0,1)   = underdamped decaying oscillation  <- BOLD regime
#                   1       = critically damped
#                   Physiological prior: zeta in [0.05, 0.30]
#                   Bounds: [0, 1]
#   p[3,] = k     : coupling strength (rad^2/TR^2)
#                   Bounds: [0, BOLD_K_UPPER] = [0, 10.659]
#
# IDENTIFIABILITY:
#   gamma -> oscillation frequency (spectral peak of trajectory)
#   zeta  -> envelope decay rate   (-zeta*sqrt(gamma) per TR)
#   k     -> mode splitting         (lambda2 - lambda1 = 2k)
#   Three distinct signatures -> all three parameters structurally identified.
#
# @param t  Scalar dummy time variable (unused).
# @param x  State matrix (4 x N_sigma): rows are [theta1, theta1d, theta2, theta2d].
# @param p  Parameter matrix (3 x N_sigma): rows are [gamma, zeta, k].
# @return   Matrix (4 x N_sigma): first derivatives [theta1d, theta1dd, theta2d, theta2dd].
coupled_osc_model_glk_2nd <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)

  gamma      <- p[1, ]
  zeta       <- p[2, ]
  k          <- p[3, ]

  theta1     <- x[1, ]     # position,  region 1  (observed)
  theta1_dot <- x[2, ]     # velocity,  region 1  (latent)
  theta2     <- x[3, ]     # position,  region 2  (observed)
  theta2_dot <- x[4, ]     # velocity,  region 2  (latent)

  # Linear damping coefficient c = 2*zeta*sqrt(gamma)
  # pmax guards against negative gamma from optimizer excursions
  c <- 2 * zeta * sqrt(pmax(gamma, 0))

  # Accelerations (same physics as original models, now correctly used as d(vel)/dt)
  theta1_ddot <- -gamma * sin(theta1) - c * theta1_dot - k * (sin(theta1) - sin(theta2))
  theta2_ddot <- -gamma * sin(theta2) - c * theta2_dot + k * (sin(theta1) - sin(theta2))

  # Return d/dt[theta1, theta1d, theta2, theta2d] = [theta1d, theta1dd, theta2d, theta2dd]
  rbind(theta1_dot,    # dtheta1/dt   = theta1d   (pass velocity through)
        theta1_ddot,   # dtheta1d/dt  = theta1dd  (computed from physics)
        theta2_dot,    # dtheta2/dt   = theta2d
        theta2_ddot)   # dtheta2d/dt  = theta2dd
}


# -----------------------------------------------------------------------------
# coupled_osc_model_abk_2nd
#
# Second-order ASYMMETRIC coupled oscillator with damping.
# Regions have separate intrinsic frequencies a, b.
# State is [theta1, theta1d, theta2, theta2d]  (N_y = 4).
# Use for pairs where the GLK constraint a = b is violated
# (Dchi2 = chi2_glk - chi2_abk above threshold, flagged in main pipeline S6).
#
# EQUATIONS OF MOTION:
#   theta1_dd = -a*sin(theta1) - 2*zeta*sqrt(a)*theta1d - k*(sin(theta1)-sin(theta2))
#   theta2_dd = -b*sin(theta2) - 2*zeta*sqrt(b)*theta2d + k*(sin(theta1)-sin(theta2))
#
# PARAMETERS (p rows):
#   p[1,] = a    : intrinsic squared frequency, region 1 (rad^2/TR^2)
#   p[2,] = b    : intrinsic squared frequency, region 2 (rad^2/TR^2)
#   p[3,] = zeta : shared damping ratio (dimensionless, >= 0)
#   p[4,] = k    : coupling strength (rad^2/TR^2)
#
# @param t  Scalar dummy time variable (unused).
# @param x  State matrix (4 x N_sigma).
# @param p  Parameter matrix (4 x N_sigma): rows are [a, b, zeta, k].
# @return   Matrix (4 x N_sigma) of first derivatives.
coupled_osc_model_abk_2nd <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)

  a          <- p[1, ]
  b          <- p[2, ]
  zeta       <- p[3, ]
  k          <- p[4, ]

  theta1     <- x[1, ]
  theta1_dot <- x[2, ]
  theta2     <- x[3, ]
  theta2_dot <- x[4, ]

  ca <- 2 * zeta * sqrt(pmax(a, 0))
  cb <- 2 * zeta * sqrt(pmax(b, 0))

  theta1_ddot <- -a * sin(theta1) - ca * theta1_dot - k * (sin(theta1) - sin(theta2))
  theta2_ddot <- -b * sin(theta2) - cb * theta2_dot + k * (sin(theta1) - sin(theta2))

  rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
}


# =============================================================================
# LEGACY MODELS  (first-order formulation - do NOT use for new analyses)
# Preserved to allow exact reproduction of pre-2026-03 results.
# =============================================================================

coupled_osc_model_abk <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)
  a <- p[1,]; b <- p[2,]; k <- p[3,]
  theta1 <- x[1,]; theta2 <- x[2,]
  rbind(-a*sin(theta1) - k*(sin(theta1)-sin(theta2)),
        -b*sin(theta2) + k*(sin(theta1)-sin(theta2)))
}

coupled_osc_model_glk <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)
  gamma <- p[1,]; k <- p[2,]
  theta1 <- x[1,]; theta2 <- x[2,]
  rbind(-gamma*sin(theta1) - k*(sin(theta1)-sin(theta2)),
        -gamma*sin(theta2) + k*(sin(theta1)-sin(theta2)))
}

coupled_osc_model_gLk <- function(t, x, p) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1)
  if (is.null(dim(p))) p <- matrix(p, ncol = 1)
  g <- p[1,]; L <- p[2,]; k <- p[3,]
  theta1 <- x[1,]; theta2 <- x[2,]
  rbind((-g*sin(theta1) - k*L*(sin(theta1)-sin(theta2)))/L,
        (-g*sin(theta2) + k*L*(sin(theta1)-sin(theta2)))/L)
}


# -----------------------------------------------------------------------------
# coupled_osc_model_glk2_fz  (fixed-zeta variant — PRIMARY MODEL from v4)
#
# Identical physics to coupled_osc_model_glk_2nd, but ζ is NOT a free parameter.
# It is baked in at construction time via make_glk2_fixed_zeta().
#
# RATIONALE (2026-03-09):
#   Three successive runs showed ζ piling up at whatever ZETA_UPPER was set to
#   (1.0 → 25% ceiling; 0.707 → 55% ceiling; p50 = exactly ceiling).
#   Bootstrap stability_k = −2,677,706 confirmed ζ has no identifiable gradient
#   in BOLD chi² at TR=2s / 260 TRs.  ζ must be fixed, not estimated.
#   Fixing ζ at 0.1 (lightly underdamped) restores the 2-parameter (γ, k)
#   identification proven in Phase 7, now with correct 2nd-order integration.
#
# PARAMETERS (p rows):
#   p[1,] = gamma : shared squared natural frequency (rad²/TR²)
#   p[2,] = k     : coupling strength (rad²/TR²)
#
# FIXED (closed over):
#   zeta_fixed    : damping ratio — use 0.1 for production; run sensitivity
#                   over {0.05, 0.10, 0.20} to confirm k is insensitive.
#
# @param zeta_fixed  Scalar damping ratio (default 0.1).
# @return  ODE function with signature (t, x, p) compatible with UKF engine.
make_glk2_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)   # capture in closure
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    gamma      <- p[1, ]
    k          <- p[2, ]
    zeta       <- zeta_fixed   # fixed — NOT taken from p

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    c_damp <- 2 * zeta * sqrt(pmax(gamma, 0))

    theta1_ddot <- -gamma * sin(theta1) - c_damp * theta1_dot - k * (sin(theta1) - sin(theta2))
    theta2_ddot <- -gamma * sin(theta2) - c_damp * theta2_dot + k * (sin(theta1) - sin(theta2))

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

# Production model: ζ fixed at 0.1 (lightly underdamped)
coupled_osc_model_glk2_fz <- make_glk2_fixed_zeta(0.1)


# -----------------------------------------------------------------------------
# make_abk2_fixed_zeta  (fixed-zeta ABK — for fair identifiability test, v5+)
#
# Identical physics to coupled_osc_model_abk_2nd, but ζ is closed over rather
# than taken from the parameter vector.  This gives the ABK model the same
# number of free parameters as glk2_fz (both have 2: either [γ,k] or [a,b,k]).
#
# With free ζ in abk but fixed ζ in gLk, the v4 Δchi² test was invalid:
# gLk beat abk purely through optimizer efficiency (2 vs 4 free params, same
# step budget).  With both models using fixed ζ, Δchi² isolates the frequency-
# asymmetry question (a ≠ b) cleanly.
#
# PARAMETERS (p rows) for abk2_fz:
#   p[1,] = a  : intrinsic squared frequency, region 1 (rad²/TR²)
#   p[2,] = b  : intrinsic squared frequency, region 2 (rad²/TR²)
#   p[3,] = k  : coupling strength (rad²/TR²)
#
# @param zeta_fixed  Scalar damping ratio to fix (default 0.1, matches gLk2_fz).
# @return  ODE function (t, x, p) compatible with UKF engine.
make_abk2_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    a          <- p[1, ]
    b          <- p[2, ]
    k          <- p[3, ]
    zeta       <- zeta_fixed   # fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    ca <- 2 * zeta * sqrt(pmax(a, 0))
    cb <- 2 * zeta * sqrt(pmax(b, 0))

    theta1_ddot <- -a * sin(theta1) - ca * theta1_dot - k * (sin(theta1) - sin(theta2))
    theta2_ddot <- -b * sin(theta2) - cb * theta2_dot + k * (sin(theta1) - sin(theta2))

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

# Production fixed-ζ ABK model (ζ=0.1, matching glk2_fz)
coupled_osc_model_abk2_fz <- make_abk2_fixed_zeta(0.1)


# -----------------------------------------------------------------------------
# make_lho_fixed_zeta  (Linearised Harmonic Oscillator — LHO, v5+)
#
# Replaces sin(θ) with θ in both the restoring force and coupling term of
# glk2_fz.  For BOLD fMRI, signal amplitudes are ~1-3% of baseline, so
# sin(θ) ≈ θ to 5 significant figures — the nonlinearity was physically
# vacuous and computationally harmful.
#
# The linearisation has two important consequences:
#
#   1. Coupling identifiability: k(θ1-θ2) in the state equations maps to
#      a specific off-diagonal cross-spectral density signature.  In the
#      frequency domain, k splits the single peak at ω=√γ into two normal
#      modes at ω₁=√γ and ω₂=√(γ+2k).  The chi² surface has a genuine
#      bowl at the true k (not a ramp toward the boundary).
#
#   2. State equations are linear conditional on γ and ζ.  Sigma-point
#      propagation (UKF) is exact for linear systems, eliminating
#      linearisation error from the Cholesky decomposition.
#
# EQUATIONS OF MOTION:
#   θ̈₁ = -γ·θ₁ - 2ζ√γ·θ̇₁ - k·(θ₁ - θ₂)
#   θ̈₂ = -γ·θ₂ - 2ζ√γ·θ̇₂ + k·(θ₁ - θ₂)
#
# PARAMETERS (p rows):
#   p[1,] = γ  : squared natural frequency, rad²/TR²  (shared — symmetric model)
#   p[2,] = k  : linear coupling strength, rad²/TR²
#
# FIXED: ζ closed over (unidentifiable from BOLD; fixed at 0.1).
#
# STATE (x rows): [θ₁, θ̇₁, θ₂, θ̇₂]  — identical layout to glk2_fz
# OBSERVED:       rows N_p+1, N_p+3  (θ₁, θ₂)
#
# @param zeta_fixed  Damping ratio (default 0.1).
# @return  ODE function (t, x, p) compatible with UKF engine.
make_lho_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    gamma      <- p[1, ]
    k          <- p[2, ]
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    c_damp <- 2 * zeta * sqrt(pmax(gamma, 0))

    # Linearised: sin(θ) → θ  and  sin(θ₁)-sin(θ₂) → (θ₁-θ₂)
    theta1_ddot <- -gamma * theta1 - c_damp * theta1_dot - k * (theta1 - theta2)
    theta2_ddot <- -gamma * theta2 - c_damp * theta2_dot + k * (theta1 - theta2)

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

# Production LHO model (ζ=0.1 fixed)
coupled_osc_model_lho_fz <- make_lho_fixed_zeta(0.1)

# Asymmetric-frequency LHO: p = [a, b, k]
# Allows region 1 and region 2 to have different natural frequencies.
# Used for identifiability comparison (matches abk2_fz role for gLk).
make_lho_abk_fixed_zeta <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    a          <- p[1, ]
    b          <- p[2, ]
    k          <- p[3, ]
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    ca <- 2 * zeta * sqrt(pmax(a, 0))
    cb <- 2 * zeta * sqrt(pmax(b, 0))

    theta1_ddot <- -a * theta1 - ca * theta1_dot - k * (theta1 - theta2)
    theta2_ddot <- -b * theta2 - cb * theta2_dot + k * (theta1 - theta2)

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}

coupled_osc_model_lho_abk_fz <- make_lho_abk_fixed_zeta(0.1)


# -----------------------------------------------------------------------------
# make_lho_single  (single-region LHO — Stage 1 of two-stage pipeline)
#
# Fits only γ to a single BOLD ROI time series.
# No coupling term. Simplest possible identifiable model.
#
# EQUATIONS:
#   θ̈ = -γ·θ - 2ζ√γ·θ̇
#
# PARAMETERS (p rows):
#   p[1,] = γ  : squared natural frequency, rad²/TR²
#
# FIXED: ζ closed over (default 0.1)
#
# STATE (x rows): [θ, θ̇]   N_y = 2
# OBSERVED:       row N_p+1 = θ  (position only)
#
# ts_data for this model has 3 columns: [time, θ, θ̇]
# (use prepare_ukf_input_single to build from smoothed matrix)
#
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
make_lho_single <- function(zeta_fixed = 0.1) {
  force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    gamma      <- p[1, ]
    zeta       <- zeta_fixed
    theta      <- x[1, ]
    theta_dot  <- x[2, ]

    c_damp     <- 2 * zeta * sqrt(pmax(gamma, 0))
    theta_ddot <- -gamma * theta - c_damp * theta_dot

    rbind(theta_dot, theta_ddot)
  }
}

lho_single_model <- make_lho_single(0.1)

# -----------------------------------------------------------------------------
# make_lho_fixed_ab  (two-stage coupled LHO — Stage 2)
#
# Both natural frequencies a (region 1) and b (region 2) are fixed from
# Stage 1 estimates.  Only coupling k is free.  This gives a strictly 1D
# chi² surface in k, which is necessary and sufficient for k to be
# identifiable via the normal-mode splitting mechanism.
#
# EQUATIONS:
#   θ̈₁ = -a·θ₁ - 2ζ√a·θ̇₁ - k·(θ₁ - θ₂)
#   θ̈₂ = -b·θ₂ - 2ζ√b·θ̇₂ + k·(θ₁ - θ₂)
#
# PARAMETERS (p rows):
#   p[1,] = k  : linear coupling strength, rad²/TR²
#
# FIXED: a, b, ζ closed over from Stage 1
#
# STATE (x rows): [θ₁, θ̇₁, θ₂, θ̇₂]  — same layout as all paired models
# OBSERVED:       rows N_p+1, N_p+3  (θ₁, θ₂)
#
# Normal mode frequencies:
#   ω₊ = √( (a+b)/2 + k + √( ((a-b)/2)² + k² ) )  / (2πTR)
#   ω₋ = √( (a+b)/2 + k - √( ((a-b)/2)² + k² ) )  / (2πTR)
# For a=b=γ: ω₋=√γ/(2πTR), ω₊=√(γ+2k)/(2πTR) — the bowl exists for
# k ∈ [0, ((0.1·2πTR)² - γ)/2], i.e. while both modes are in the BOLD band.
#
# @param a_fixed   Squared frequency for region 1 (rad²/TR², from Stage 1)
# @param b_fixed   Squared frequency for region 2 (rad²/TR², from Stage 1)
# @param zeta_fixed  Damping ratio (default 0.1)
# @return ODE function (t, x, p)
make_lho_fixed_ab <- function(a_fixed, b_fixed, zeta_fixed = 0.1) {
  force(a_fixed); force(b_fixed); force(zeta_fixed)
  function(t, x, p) {
    if (is.null(dim(x))) x <- matrix(x, ncol = 1)
    if (is.null(dim(p))) p <- matrix(p, ncol = 1)

    k          <- p[1, ]
    a          <- a_fixed
    b          <- b_fixed
    zeta       <- zeta_fixed

    theta1     <- x[1, ]
    theta1_dot <- x[2, ]
    theta2     <- x[3, ]
    theta2_dot <- x[4, ]

    ca <- 2 * zeta * sqrt(max(a, 0))
    cb <- 2 * zeta * sqrt(max(b, 0))

    theta1_ddot <- -a * theta1 - ca * theta1_dot - k * (theta1 - theta2)
    theta2_ddot <- -b * theta2 - cb * theta2_dot + k * (theta1 - theta2)

    rbind(theta1_dot, theta1_ddot, theta2_dot, theta2_ddot)
  }
}
