# =============================================================================
# preprocessing.R
# Data loading and physiologically-grounded kernel smoothing for MDD fMRI
# resting-state signals.
#
# Study context (from paper):
#   - Harvard-Oxford Atlas, 110 ROIs, 23 MDD subjects
#   - TR = 2000 ms (exactly), 260 usable volumes per run
#   - Signals already detrended and standardized by NiftiSpheresMasker
#   - BOLD resting-state band of interest: 0.01-0.1 Hz
#
# Smoothing strategy:
#   A Gaussian kernel with bandwidth h (time points) attenuates frequencies
#   above approximately f_cutoff = 0.44 / (h * TR) Hz.
#
#   Given TR = 2 s:
#     h = 0.44 / (f_cutoff_Hz * TR)
#
#   Preset bandwidths (TR = 2s):
#     "light"    -> h = 2.2  TRs  (cutoff ~ 0.10 Hz, full BOLD band preserved)
#     "moderate" -> h = 4.4  TRs  (cutoff ~ 0.05 Hz, mid-band)
#     "strong"   -> h = 8.8  TRs  (cutoff ~ 0.025 Hz, slow fluctuations only)
#
#   WHY LOO-CV WAS REMOVED:
#     LOO-CV on a kernel smoother minimizes in-sample MSE, which decreases
#     monotonically as bandwidth decreases (less smoothing = closer to raw).
#     It always returns bw_min, providing no useful selection. On time series
#     with temporal autocorrelation (which all BOLD signals have), true
#     leave-one-out cross-validation is also unreliable because adjacent points
#     are correlated — removing one point provides almost no independence.
#     With TR and acquisition parameters known exactly, a principled
#     frequency-derived bandwidth is both simpler and more meaningful.
#
# Functions:
#   fmri_bandwidth_bounds()  -- Compute bw_min / bw_max from TR and Hz limits
#   freq_to_bandwidth()      -- Convert a frequency cutoff (Hz) to bandwidth (TRs)
#   select_bandwidth()       -- Select bandwidth: preset string or numeric Hz cutoff
#   smooth_subject_data()    -- Smooth all ROI columns, return smoothed + report
#   diagnose_smoothing()     -- Multi-bandwidth overlay plot for one ROI
#   plot_bw_report()         -- Per-ROI bandwidth and variance retention plots
#   load_subject_data()      -- Load subject CSVs into a named list
#   build_region_pairs()     -- All unique 2-ROI combinations
#   prepare_ukf_input()      -- Extract one pair as UKF-ready matrix
# =============================================================================

library(KernSmooth)
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})


# =============================================================================
# ACQUISITION CONSTANTS — set once, used everywhere
# =============================================================================

FMRI_ACQ <- list(
  TR            = 2.0,    # repetition time in seconds (GE MR750, paper Sec 2.1)
  BOLD_FREQ_MIN = 0.01,   # Hz — lower bound of resting-state BOLD band
  BOLD_FREQ_MAX = 0.10,   # Hz — upper bound of resting-state BOLD band
  N_VOLS        = 260L,   # usable volumes per run (263 acquired - 3 discarded)
  VAR_THRESHOLD = 0.30    # warn (not error) if retained variance drops below this
)

# Named preset bandwidths — use these for smooth_subject_data(preset = ...)
BW_PRESETS <- list(
  light    = 0.44 / (0.10 * 2.0),   # 2.2  TRs — preserves up to 0.10 Hz
  moderate = 0.44 / (0.05 * 2.0),   # 4.4  TRs — preserves up to 0.05 Hz
  strong   = 0.44 / (0.025 * 2.0)   # 8.8  TRs — preserves up to 0.025 Hz
)


# -----------------------------------------------------------------------------
#' fmri_bandwidth_bounds
#'
#' Returns the physiological [bw_min, bw_max] window in TR units.
#' bw_min corresponds to f_max (least smoothing).
#' bw_max corresponds to f_min (most smoothing).
#'
#' @param TR           Repetition time in seconds.
#' @param freq_min_hz  Lower BOLD frequency bound (Hz).
#' @param freq_max_hz  Upper BOLD frequency bound (Hz).
#' @return Named numeric vector: bw_min, bw_max.
#' @export
fmri_bandwidth_bounds <- function(TR          = FMRI_ACQ$TR,
                                   freq_min_hz = FMRI_ACQ$BOLD_FREQ_MIN,
                                   freq_max_hz = FMRI_ACQ$BOLD_FREQ_MAX) {
  nyquist <- 1 / (2 * TR)
  if (freq_max_hz > nyquist)
    warning(sprintf("freq_max_hz (%.3f Hz) exceeds Nyquist (%.3f Hz) for TR=%.1f s.",
                    freq_max_hz, nyquist, TR))
  c(bw_min = 0.44 / (freq_max_hz * TR),
    bw_max = 0.44 / (freq_min_hz * TR))
}


# -----------------------------------------------------------------------------
#' freq_to_bandwidth
#'
#' Converts a frequency cutoff in Hz to a Gaussian kernel bandwidth in TRs.
#' Bandwidth h attenuates frequencies above f_cutoff ≈ 0.44 / (h * TR).
#'
#' @param freq_hz  Desired frequency cutoff in Hz (higher = less smoothing).
#' @param TR       Repetition time in seconds.
#' @return Scalar bandwidth in TR units.
#' @export
freq_to_bandwidth <- function(freq_hz, TR = FMRI_ACQ$TR) {
  bounds <- fmri_bandwidth_bounds(TR)
  h <- 0.44 / (freq_hz * TR)
  if (h < bounds["bw_min"] - 1e-9)
    warning(sprintf("freq_hz=%.4f Hz gives h=%.2f which is below bw_min=%.2f (above Nyquist range).",
                    freq_hz, h, bounds["bw_min"]))
  if (h > bounds["bw_max"] + 1e-9)
    warning(sprintf("freq_hz=%.4f Hz gives h=%.2f which is above bw_max=%.2f.",
                    freq_hz, h, bounds["bw_max"]))
  unname(h)
}


# -----------------------------------------------------------------------------
#' select_bandwidth
#'
#' Determines the Gaussian kernel bandwidth for one ROI signal using a
#' frequency-based rule rather than data-driven selection.
#'
#' The 'preset' argument accepts:
#'   "light"    — 2.2  TRs, preserves BOLD band up to 0.10 Hz (recommended start)
#'   "moderate" — 4.4  TRs, preserves up to 0.05 Hz
#'   "strong"   — 8.8  TRs, preserves slow fluctuations up to 0.025 Hz
#'   A numeric Hz value — converted to bandwidth via freq_to_bandwidth()
#'   A numeric TRs value > 1 (treated as direct bandwidth, clipped to window)
#'
#' After smoothing, variance retention is checked. A warning is issued if it
#' falls below VAR_THRESHOLD (but no automatic fallback — the user decides).
#'
#' @param signal        Numeric vector — one ROI time series.
#' @param preset        Bandwidth choice: "light", "moderate", "strong",
#'                      numeric Hz (≤ 0.5), or numeric TRs (> 0.5).
#' @param TR            Repetition time in seconds.
#' @param verbose       Print selection details if TRUE.
#' @param col_name      ROI name used in messages.
#' @return Named list: h (bandwidth TRs), method (string), var_ratio (numeric).
#' @export
select_bandwidth <- function(signal,
                               preset    = "moderate",
                               TR        = FMRI_ACQ$TR,
                               verbose   = FALSE,
                               col_name  = "") {

  bounds <- fmri_bandwidth_bounds(TR)

  # Resolve preset to a bandwidth value
  if (is.character(preset)) {
    if (!preset %in% names(BW_PRESETS))
      stop(sprintf("Unknown preset '%s'. Choose from: %s",
                   preset, paste(names(BW_PRESETS), collapse = ", ")))
    h      <- BW_PRESETS[[preset]]
    method <- sprintf("preset='%s' (h=%.2f TRs, cutoff~%.3f Hz)",
                      preset, h, 0.44 / (h * TR))

  } else if (is.numeric(preset) && preset <= 0.5) {
    # Treat as frequency cutoff in Hz
    h      <- freq_to_bandwidth(preset, TR)
    method <- sprintf("freq_cutoff=%.4f Hz -> h=%.2f TRs", preset, h)

  } else if (is.numeric(preset) && preset > 0.5) {
    # Treat as direct bandwidth in TRs — clip to physiological window
    h_raw  <- preset
    h      <- max(bounds["bw_min"], min(bounds["bw_max"], h_raw))
    method <- sprintf("direct h=%.2f TRs", h)
    if (abs(h - h_raw) > 1e-9)
      message(sprintf("  [%s] bandwidth %.2f clipped to physiological window [%.2f, %.2f]",
                       col_name, h_raw, bounds["bw_min"], bounds["bw_max"]))
  } else {
    stop("preset must be a character string or numeric value.")
  }

  # Compute variance retention as a diagnostic only
  n         <- length(signal)
  t_pts     <- seq_len(n)
  smoothed  <- ksmooth(t_pts, signal, kernel = "normal",
                        bandwidth = h, n.points = n)$y
  # Compute raw variance robustly (handle NA or constant signals)
  raw_var <- var(signal, na.rm = TRUE)
  if (!is.finite(raw_var) || raw_var <= 0) {
    var_ratio <- 1.0
  } else {
    var_ratio <- var(smoothed, na.rm = TRUE) / raw_var
  }

  if (var_ratio < FMRI_ACQ$VAR_THRESHOLD)
    warning(sprintf(
      "[%s] Low variance retention: %.1f%% < threshold %.0f%%. Consider using 'light' preset.",
      col_name, var_ratio * 100, FMRI_ACQ$VAR_THRESHOLD * 100))

  if (verbose)
    message(sprintf("  [%-12s] h=%.2f TRs  var_ratio=%.3f  (%s)",
                     col_name, h, var_ratio, method))

  list(h = h, method = method, var_ratio = var_ratio)
}


# -----------------------------------------------------------------------------
#' smooth_subject_data
#'
#' Applies Gaussian kernel smoothing to every ROI column of a subject data
#' frame using a single consistent bandwidth across all ROIs (chosen by preset
#' or explicit value).
#'
#' Returns both the smoothed data frame and a per-ROI bandwidth report.
#'
#' @param subject_df  Data frame: col 1 = time, remaining cols = ROI signals.
#' @param preset      Bandwidth choice passed to select_bandwidth().
#'                    Default "moderate" (4.4 TRs, cutoff ~0.05 Hz).
#' @param TR          Repetition time in seconds.
#' @param verbose     Print per-ROI details if TRUE.
#' @return List:
#'   $smoothed   — data.frame same structure as subject_df, values smoothed
#'   $bw_report  — data.frame: ROI, bandwidth, method, var_ratio
#' @export
smooth_subject_data <- function(subject_df,
                                 preset  = "moderate",
                                 TR      = FMRI_ACQ$TR,
                                 verbose = FALSE) {

  stopifnot(is.data.frame(subject_df) || is.matrix(subject_df),
            ncol(subject_df) >= 2)

  subject_df  <- as.data.frame(subject_df)
  time_vec    <- subject_df[[1]]
  n           <- nrow(subject_df)
  t_pts       <- seq_len(n)
  region_cols <- subject_df[, -1, drop = FALSE]
  roi_names   <- colnames(region_cols)

  smoothed_mat <- matrix(NA_real_, nrow = n, ncol = length(roi_names),
                          dimnames = list(NULL, roi_names))
  bw_records   <- vector("list", length(roi_names))

  for (i in seq_along(roi_names)) {
    col <- roi_names[i]
    sig <- region_cols[[col]]

    sel <- select_bandwidth(signal   = sig,
                             preset   = preset,
                             TR       = TR,
                             verbose  = verbose,
                             col_name = col)

    smoothed_mat[, i] <- ksmooth(t_pts, sig, kernel = "normal",
                                   bandwidth = sel$h, n.points = n)$y

    bw_records[[i]] <- data.frame(ROI       = col,
                                   bandwidth = sel$h,
                                   method    = sel$method,
                                   var_ratio = sel$var_ratio,
                                   stringsAsFactors = FALSE)
  }

  smoothed_df           <- data.frame(time_vec, smoothed_mat,
                                       stringsAsFactors = FALSE)
  colnames(smoothed_df) <- colnames(subject_df)

  list(smoothed  = smoothed_df,
       bw_report = do.call(rbind, bw_records))
}


# -----------------------------------------------------------------------------
#' diagnose_smoothing
#'
#' Overlays the raw signal with smoothed curves for all three presets
#' (light / moderate / strong) plus any user-specified bandwidth, so you can
#' visually pick the right level of smoothing for your data.
#'
#' @param signal       Numeric vector — one ROI time series.
#' @param bw_selected  Bandwidth (TRs) to highlight in red (e.g. the one in use).
#' @param TR           Repetition time in seconds.
#' @param signal_name  ROI label for the plot title.
#' @return ggplot object (printed as side-effect).
#' @export
diagnose_smoothing <- function(signal,
                                bw_selected = NULL,
                                TR          = FMRI_ACQ$TR,
                                signal_name = "ROI") {

  n     <- length(signal)
  t_pts <- seq_len(n)

  # Always show all three presets
  preset_bws <- unlist(BW_PRESETS)
  all_bws    <- sort(unique(c(preset_bws, bw_selected)))

  preset_labels <- setNames(
    sprintf("bw=%.1f TRs (~%.3f Hz)", preset_bws, 0.44 / (preset_bws * TR)),
    names(preset_bws)
  )

  make_label <- function(h) {
    nm <- names(preset_bws)[which(abs(preset_bws - h) < 1e-9)]
    if (length(nm) == 1)
      sprintf("%s: bw=%.1f (var=%.0f%%)", nm,
              h, var(ksmooth(t_pts, signal, "normal", h, n)$y) / var(signal) * 100)
    else
      sprintf("bw=%.1f (var=%.0f%%)",
              h, var(ksmooth(t_pts, signal, "normal", h, n)$y) / var(signal) * 100)
  }

  rows <- lapply(all_bws, function(h) {
    data.frame(Time  = t_pts,
               Value = ksmooth(t_pts, signal, "normal", h, n)$y,
               Label = make_label(h),
               stringsAsFactors = FALSE)
  })
  raw_df <- data.frame(Time = t_pts, Value = signal,
                        Label = "Raw signal", stringsAsFactors = FALSE)
  df     <- rbind(raw_df, do.call(rbind, rows))

  # Colours: grey for raw, blues for presets, red for selected
  preset_colours <- c(light = "#3498DB", moderate = "#1A5276", strong = "#AED6F1")
  curve_cols <- setNames(
    sapply(all_bws, function(h) {
      nm <- names(preset_bws)[which(abs(preset_bws - h) < 1e-9)]
      if (length(nm) == 1) preset_colours[[nm]] else "#7D3C98"
    }),
    sapply(all_bws, make_label)
  )
  all_cols <- c("Raw signal" = "grey65", curve_cols)

  if (!is.null(bw_selected)) {
    sel_lbl <- make_label(bw_selected)
    if (sel_lbl %in% names(all_cols)) all_cols[sel_lbl] <- "#E74C3C"
  }

  lw <- setNames(
    c(0.4, ifelse(abs(all_bws - ifelse(is.null(bw_selected), -1, bw_selected)) < 1e-9,
                   1.5, 0.9)),
    c("Raw signal", sapply(all_bws, make_label))
  )

  bounds <- fmri_bandwidth_bounds(TR)

  ggplot(df, aes(x = Time, y = Value, colour = Label, linewidth = Label)) +
    geom_line(alpha = 0.9) +
    scale_colour_manual(values = all_cols) +
    scale_linewidth_manual(values = lw) +
    labs(
      title    = paste("Smoothing Diagnosis —", signal_name),
      subtitle = sprintf(
        "Physiological window: [%.1f, %.1f] TRs  |  TR=%.1f s  |  Band: %.2f–%.2f Hz",
        bounds["bw_min"], bounds["bw_max"], TR,
        FMRI_ACQ$BOLD_FREQ_MIN, FMRI_ACQ$BOLD_FREQ_MAX),
      x = "Time (TR)", y = "Standardized BOLD signal",
      colour = NULL, linewidth = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "right",
          legend.text     = element_text(size = 9),
          plot.subtitle   = element_text(size = 9, colour = "grey40"))
}


# -----------------------------------------------------------------------------
#' plot_bw_report
#'
#' Plots per-ROI bandwidth and variance retention from smooth_subject_data().
#'
#' @param bw_report  $bw_report from smooth_subject_data().
#' @param subject_id String for plot title.
#' @return Invisibly returns list of two ggplot objects.
#' @export
plot_bw_report <- function(bw_report, subject_id = "") {
  bounds <- fmri_bandwidth_bounds()

  p1 <- ggplot(bw_report,
               aes(x = reorder(ROI, bandwidth), y = bandwidth)) +
    geom_col(fill = "#2980B9", width = 0.7) +
    geom_hline(yintercept = bounds["bw_min"], linetype = "dashed",
               colour = "steelblue", linewidth = 0.8) +
    geom_hline(yintercept = bounds["bw_max"], linetype = "dashed",
               colour = "firebrick", linewidth = 0.8) +
    # Preset reference lines
    geom_hline(yintercept = BW_PRESETS$light,    linetype = "dotted",
               colour = "#27AE60", linewidth = 0.7) +
    geom_hline(yintercept = BW_PRESETS$moderate, linetype = "dotted",
               colour = "#E67E22", linewidth = 0.7) +
    geom_hline(yintercept = BW_PRESETS$strong,   linetype = "dotted",
               colour = "#8E44AD", linewidth = 0.7) +
    annotate("text", x = nrow(bw_report) * 0.02,
             y = c(bounds["bw_min"] + 0.4, BW_PRESETS$moderate + 0.4,
                   bounds["bw_max"] - 0.8),
             label = c(sprintf("bw_min=%.1f (light)", bounds["bw_min"]),
                       sprintf("moderate=%.1f", BW_PRESETS$moderate),
                       sprintf("bw_max=%.1f", bounds["bw_max"])),
             colour = c("steelblue", "#E67E22", "firebrick"),
             size = 3, hjust = 0) +
    labs(title = paste("Selected bandwidth per ROI",
                        if (nchar(subject_id) > 0) paste0("— ", subject_id)),
         x = NULL, y = "Bandwidth (TRs)") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))

  p2 <- ggplot(bw_report,
               aes(x = reorder(ROI, var_ratio),
                   y = var_ratio * 100,
                   fill = var_ratio >= FMRI_ACQ$VAR_THRESHOLD)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = FMRI_ACQ$VAR_THRESHOLD * 100,
               linetype = "dashed", colour = "firebrick", linewidth = 0.8) +
    scale_fill_manual(
      values = c("TRUE" = "#2ECC71", "FALSE" = "#E74C3C"),
      labels = c("TRUE" = "Above threshold", "FALSE" = "Below threshold")) +
    labs(title    = "Variance retained after smoothing per ROI",
         subtitle = sprintf("Dashed = %.0f%% warning threshold",
                             FMRI_ACQ$VAR_THRESHOLD * 100),
         x = NULL, y = "Retained variance (%)", fill = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
          legend.position = "top")

  print(p1)
  print(p2)
  invisible(list(bandwidth_plot = p1, variance_plot = p2))
}


# =============================================================================
# DATA I/O
# =============================================================================

#' load_subject_data — load all subject CSVs into a named list
#' @export
load_subject_data <- function(folder,
                               pattern = "errts.*_rest(2)?_restroi_signals\\.csv$",
                               verbose = TRUE) {
  if (!dir.exists(folder)) stop("Folder not found: ", folder)
  files <- list.files(folder, pattern = pattern, full.names = TRUE)
  if (length(files) == 0) stop("No matching files found in: ", folder)
  var_names        <- gsub("errts\\.|_restroi_signals\\.csv", "", basename(files))
  data_list        <- vector("list", length(files))
  names(data_list) <- var_names
  for (i in seq_along(files)) {
    data_list[[i]] <- read.csv(files[[i]], header = TRUE)
    if (verbose)
      cat(sprintf("Loaded %-55s -> '%s'\n", basename(files[[i]]), var_names[[i]]))
  }
  data_list
}

#' build_region_pairs — all unique 2-ROI combinations
#' @export
build_region_pairs <- function(smoothed_df) {
  combn(colnames(smoothed_df)[-1], 2, simplify = FALSE)
}

#' prepare_ukf_input — extract one pair as UKF-ready matrix [time, ROI_A, ROI_B]
#' @export
prepare_ukf_input <- function(smoothed_df, region_pair) {
  stopifnot(length(region_pair) == 2,
            all(region_pair %in% colnames(smoothed_df)))
  as.matrix(smoothed_df[, c(colnames(smoothed_df)[1], region_pair)])
}
