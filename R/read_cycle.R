#' Read a single AquaResp cycle file
#'
#' AquaResp stores each experimental cycle (i.e., the measurement period when \eqn{MO_{2}} is calculated)
#' as a separate file. This function reads the data from a cycle file in the "All slopes" folder
#' of an AquaResp experiment.
#'
#' @param cycle_number Numeric. The cycle number to read.
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return A dataframe containing the data from the specified cycle.
#' @export
#'
#' @examples
#' # Locate the example AquaResp experiment directory shipped with the package
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' # Read cycle 1 from the example experiment
#' read_cycle(cycle_number = 1, path = exp_dir_path)

read_cycle <- function(cycle_number, path) {
  file_path <-
    file.path(path, "All slopes", paste0("Cycle_", cycle_number, ".txt"))

  if (!file.exists(file_path)) {
    warning(paste("Cycle file", file_path, "not found. Returning NULL."))
    return(NULL)
  }

  res <- utils::read.table(file_path, sep = ";", header = TRUE)

  if ("X" %in% names(res) && all(is.na(res$X))) {
    res <- res[, setdiff(names(res), "X"), drop = FALSE]
  }

  return(res)
}

#' Get summary PO2 values for each chamber in a cycle
#'
#' Extracts the average, median, minimum, maximum, and range of non-zero
#' \eqn{PO_{2}} values for each chamber from a cycle dataframe.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of summary \eqn{PO_{2}} values for each chamber.
#' Names are in the format `ch[chamber number].[statistic].po2`.
#' @importFrom stats median
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycle <- read_cycle(1, exp_dir_path)
#' get_cycle_min_max(cycle)
get_cycle_min_max <-
  function(cycle) {
    po2_cols <-
      grep("po2", names(cycle))

    lapply(cycle[po2_cols],
                  \(raw_x) {
      x <- raw_x[raw_x > 0 & !is.na(raw_x)]
      c(avg = mean(x), median = median(x), minimum = min(x),
        max = max(x), delta = max(x) - min(x))
      }) |> unlist() |>
      (\(x){ names(x) <- gsub(
          "^([^.]+)\\.([^.]+)\\.([^.]+)$",  # capture three dot-separated parts
          "\\1.\\3.\\2",                    # reorder as part1.part3.part2
          names(x)
        )
        return(x)})()
  }

#' Calculate percentage of missing/zero PO2 values in a cycle
#'
#' Computes the percentage of zero values in each \eqn{PO_{2}} column of a cycle
#' dataframe. Values greater than zero would have an incorrect \eqn{R^{2}} as
#' calculated by AquaResp. Missing values are defined as pO2 less than 0 or
#' NA/NaN etc. 0 is the common reading when a value is not recorded in one second
#' e.g. if the computer is busy, -300 is recorded when a firesting probe is not
#' in its housing.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @encoding UTF-8
#'
#' @return A named vector of proportion (0 to 1) indicating missing pO2 values in
#' the cycle.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycle <- read_cycle(1, exp_dir_path)
#' get_cycle_missingness(cycle)
get_cycle_missingness <-
  function(cycle) {
    po2_cols <-
      grep("po2", names(cycle))

    lapply(cycle[po2_cols], function(x)
      round(length(which(x <= 0 | !is.finite(x))) / length(x), 4) ) |>
      stats::setNames(gsub(".po2",
                           "",
                           paste0(names(cycle)[po2_cols], ".pct0"))
                      ) |> unlist()
  }

#' Calculate correlation between PO2 and time for each chamber
#'
#' Performs Pearson correlation (r) tests between \eqn{PO_{2}} values and time in seconds for each chamber in a cycle measurement period, returns r, r ^ 2, p value
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of correlation coefficients (both r and r^2 and p-values for each chamber.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycle <- read_cycle(1, exp_dir_path)
#' get_cycle_R2s(cycle)
get_cycle_R2s <-
  function(cycle) {
    po2_cols <-
      grep("po2", names(cycle))

    cors <-
      lapply(cycle[po2_cols], function(x, secs) {
        # was using 2nd column "Seconds..." renamed 'secs' above but at midnight
        # there's a bug
        secs.no0 <- secs[x != 0]
        x.no0 <- x[x != 0]
        cor_result <-
          stats::cor.test(x = x.no0, y = secs.no0, method = "pearson")
      }, secs = cycle$Unix.Time - min(cycle$Unix.Time))

    rs <- lapply(cors, function(x) {
      x$estimate |> unname()
    }) |> #returns r
      stats::setNames(gsub(".po2",
                           "",
                           paste0(names(cycle[po2_cols]), ".Pearson.R")
                           )
                      )

    r2s <- lapply(cors, function(x) {
      x$estimate^2 |> unname()
      }) |> #returns r ^ 2
      stats::setNames(gsub(".po2","",paste0(names(cycle[po2_cols]), ".R.2")))

    r2_pvals <- lapply(cors, function(x) {
    x$p.value
    }) |>
      stats::setNames(gsub(".po2","",paste0(names(cycle[po2_cols]), ".P")))

  unlist(c(rs, r2s, r2_pvals))
  }
#' Chamber selection utility
#' @keywords internal
.resolve_chambers <- function(all_chambers, chambers) {
  # If no selection, include all
  if (is.null(chambers)) return(all_chambers)

  # Normalize to integer
  chambers <- as.integer(chambers)
  if (any(is.na(chambers))) stop("`chambers` must be integer-like.")

  # Don't allow mixing positive & negative
  if (any(chambers > 0) && any(chambers < 0)) {
    stop("`chambers` cannot mix positive (include) and negative (exclude) values.")
  }

  if (all(chambers < 0)) {
    # Exclusion: drop abs(chambers)
    exclude <- abs(chambers)
    unknown <- setdiff(exclude, all_chambers)
    if (length(unknown)) stop("Unknown chambers to exclude: ", paste(unknown, collapse = ", "))
    return(setdiff(all_chambers, exclude))
  } else {
    # Inclusion: keep intersection
    include <- chambers
    unknown <- setdiff(include, all_chambers)
    if (length(unknown)) stop("Unknown chambers to include: ", paste(unknown, collapse = ", "))
    return(intersect(all_chambers, include))
  }
}


#' Calculate MO2 values for an AquaResp cycle
#'
#' Calculates mass-specific oxygen consumption (\eqn{MO_{2}}) for each chamber
#' in a given cycle using linear regression on \eqn{PO_{2}} data and experiment-level metadata. Also returns regression stats slope, standard error and intercept.
#'
#' Accepts either a **numeric cycle ID** (e.g., 1) or a **cycle dataframe** returned by `read_cycle()`.
#' `path` is always required to access metadata.

#' @param cycle Numeric or data.frame.
#'   - If \strong{numeric}: the cycle number to read (e.g., 1, 2, 3, ...).
#'   - If \strong{data.frame}: a cycle dataframe as returned by `read_cycle()`.
#' @param path Character. Path to the AquaResp experiment directory.
#' @param chambers Optional integer vector: subset of chambers to include,
#'   e.g. `c(1,2)`; or \strong{negative} to exclude, e.g. `-4` (means keep 1,2,3).
#'   Default `NULL` (all chambers).
#' @param long Logical. If `TRUE`, return a tibble in long format
#'   with columns `cycle`, `chamber`, `MO2`, `SLOPE`, `Std.Err`, `Intercept`.
#'   Default `FALSE` (wide named vector).
#' @param remove_missing Logical. If `TRUE` (default), missing \eqn{PO_{2}} measurements
#' are filtered out before calculating slopes used for MO2 estimation. `FALSE` would
#' give the same (sometimes spurious) results as AquaResp.
#'
#' @return If `long = FALSE`, a named numeric vector including
#'   `ch{N}.MO2`, `ch{N}.SLOPE`, `ch{N}.Std.Err`, `ch{N}.Intercept`.
#'   If `long = TRUE`, a tibble with columns `cycle`, `chamber`, `MO2`, `SLOPE`, `Std.Err`, `Intercept`.

calc_cycle_mo2s <- function(cycle, path, chambers = NULL, long = TRUE,
                            remove_missing = TRUE) {
  # --- Load/validate cycle data ------------------------------------------------
  if (is.numeric(cycle) && length(cycle) == 1L) {
    cycle_id   <- as.integer(cycle)
    cycle_data <- read_cycle(cycle_id, path)
  } else if (is.data.frame(cycle)) {
    cycle_data <- cycle
    # Try to detect cycle ID if present; otherwise NA
    cycle_id <- if ("Cycle" %in% names(cycle_data)) {
      unique(cycle_data$Cycle)
    } else if ("cycle" %in% names(cycle_data)) {
      unique(cycle_data$cycle)
    } else NA_integer_
    if (length(cycle_id) != 1L) cycle_id <- NA_integer_
  } else {
    stop("`cycle` must be a single numeric cycle number or a cycle data.frame as returned by `read_cycle()`.")
  }

  # --- Identify PO2 columns ----------------------------------------------------
  po2_cols <- grep("po2", names(cycle_data), value = TRUE)
  if (length(po2_cols) == 0) stop("`cycle` data contains no 'po2' columns.")

  # --- Experiment metadata & chamber selection --------------------------------
  all_chambers <- as.integer(get_chambers(path))
  sel_chambers <- .resolve_chambers(all_chambers, chambers)
  meta         <- get_exp_metadata(path)

  beta      <- meta$`Oxygen.solubilty..mg.O2...L`
  rRespFish <- meta$`Real.volume..vresp...vfish...neutrally.bouyant...L` / meta$`Mass.of.fish..kg`

  if (length(beta) != length(all_chambers) || length(rRespFish) != length(all_chambers)) {
    stop("Metadata lengths do not match number of chambers.")
  }

  # --- Time covariate: seconds since cycle start ------------------------------
  cycle_time_x <- as.numeric(cycle_data$Unix.Time) - as.numeric(min(cycle_data$Unix.Time))

  # Map each PO2 column to its chamber id and base name "chN"
  base_names <- sub("\\.po2$", "", po2_cols)                     # e.g., "ch1", "ch2"
  ch_ids     <- suppressWarnings(as.integer(sub("^ch(\\d+)$", "\\1", base_names)))

  # --- Fit PO2 ~ time, collect slope SE and intercept -------------------------
  fit_summaries <-
    lapply(seq_along(po2_cols), function(i) {
      po2s <- cycle_data[[po2_cols[i]]]
      if(remove_missing) {
        cycle_time_x <- cycle_time_x[po2s > 0 & !is.na(po2s)]
        po2s <- po2s[po2s > 0 & !is.na(po2s)]
      }
      fit_sum <- summary(stats::lm(po2s ~ cycle_time_x))
      fit_sum$coefficients  # rows: "(Intercept)", "cycle_time_x"
    })

  # --- Build outputs for selected chambers ------------------------------------
  out_list <- list()
  long_rows <- list()

  for (ch in sel_chambers) {
    idx <- which(ch_ids == ch)
    if (length(idx) != 1) next  # skip if unmatched

    coefs <- fit_summaries[[idx]]
    intercept_est <- unname(coefs["(Intercept)",   "Estimate"])
    slope_est     <- unname(coefs["cycle_time_x",  "Estimate"])
    slope_se      <- unname(coefs["cycle_time_x",  "Std. Error"])

    # Position of chamber in metadata arrays
    pos <- match(ch, all_chambers)
    if (is.na(pos)) stop("Chamber ", ch, " not found in experiment metadata order.")

    # --- MO2 (existing logic) --------------------------------------------------
    # slope_est is per-second change in PO2
    MO2_val <- -1 * (slope_est / 100) * beta[pos] * rRespFish[pos] * 3600.0

    base <- paste0("ch", ch)

    # Wide named vector pieces
    out_list[[length(out_list) + 1]] <- c(
      stats::setNames(MO2_val,       paste0(base, ".MO2")),
      stats::setNames(slope_est,     paste0(base, ".SLOPE")),
      stats::setNames(slope_se,      paste0(base, ".Std.Err")),
      stats::setNames(intercept_est, paste0(base, ".Intercept"))
    )

    # Long rows (keep extra stats available)
    long_rows[[length(long_rows) + 1]] <- data.frame(
      cycle     = cycle_id,
      chamber   = ch,
      MO2       = MO2_val,
      SLOPE     = slope_est,
      Std.Err   = slope_se,
      Intercept = intercept_est,
      check.names = FALSE
    )
  }

  if (!long) {
    # Return wide named vector
    return(unlist(out_list, use.names = TRUE))
  } else {
    # Return long tibble (with extra stats included)
    out_df <- do.call(rbind, long_rows)
    # If you want only MO2 per the docstring, select here:
    # out_df <- out_df[, c("cycle", "chamber", "MO2")]
    return(tibble::as_tibble(out_df))
  }
}

#' Pivot output from cycle_summary() to long format
#'
#' Helper to convert the wide summary into a tidy long format with columns
#' for `cycle`, `chamber`, and the metric columns: PO2 minimum/maximum,
#' cycle correlation R^2, R^2 p-value, and % missing PO2 data.
#'
#' @param df Dataframe. Output from cycle_summary().
#' @param cycle_ids Optional integer vector of cycle IDs (defaults to row order).
#'
#' @return A tibble in long format.
.pivot_cycle_summary_long <- function(df, cycle_ids = NULL) {
  if (is.null(cycle_ids)) cycle_ids <- seq_len(nrow(df))
  df$cycle <- cycle_ids

  df |>
    tidyr::pivot_longer(
      cols = -cycle,
      names_to = c("chamber", "metric"),
      names_pattern = "ch(\\d+)\\.(.+)"
    ) |>
    dplyr::mutate(chamber = as.integer(chamber)) |>
    tidyr::pivot_wider(names_from = metric, values_from = value)
}

#' Summarize all cycles in an AquaResp experiment
#'
#' Reads every cycle under `path/All slopes`, summarizes MO2, PO2 avg/median/min/max/delta,
#' R^2, R^2 p-value, and % missing values, with optional chamber subsetting and
#' long-format output.
#'
#' @param path Character. Experiment directory.
#' @param chambers Optional integer vector: subset of chambers to include
#'   (e.g. `c(1,2)`), or negative to exclude (e.g. `-4` means keep 1,2,3).
#'   Default `NULL` (all chambers).
#' @param long Logical. If `TRUE`, returns a long tibble; otherwise wide
#' data frame with `chX.` appended to the column names, where X is the chamber
#' number.
#'
#' @return A data.frame (wide) or tibble (long, default) summarizing all cycles.
#' @export
cycle_summary <- function(path, chambers = NULL, long = TRUE) {
  # Detect cycle IDs
  cycles <-
    gsub("Cycle_|\\.txt$", "", list.files(file.path(path, "All slopes"))) |>
    as.numeric() |>
    sort()

  all_chambers <- as.integer(get_chambers(path))
  sel_chambers <- .resolve_chambers(all_chambers, chambers)

  # Utility: filter a named vector to selected chambers
  filter_vec_by_ch <- function(v, sel) {
    rx <- sprintf("^ch(%s)\\.", paste(sel, collapse = "|"))
    keep <- grepl(rx, names(v))
    v[keep]
  }

  out <- sapply(
    cycles,
    function(cycle_id, path, sel_chambers) {
      cycle_data <- read_cycle(cycle_id, path)

      # Wide vectors with consistent names
      mo2   <- calc_cycle_mo2s(cycle_data, path,
                               chambers = sel_chambers,
                               long = FALSE)
      mm    <- filter_vec_by_ch(get_cycle_min_max(cycle_data), sel_chambers)
      r2    <- filter_vec_by_ch(get_cycle_R2s(cycle_data), sel_chambers)
      miss  <- filter_vec_by_ch(get_cycle_missingness(cycle_data), sel_chambers)

      c(mo2, mm, r2, miss)
    },
    path = path,
    sel_chambers = sel_chambers
  )

  wide <- t(out) |>
    as.data.frame()

  # wide[] <- lapply(wide, as.numeric)

  rownames(wide) <- NULL

  if (!long) {
    # Attach cycle IDs as a column for convenience
    wide$cycle <- cycles
    # Reorder to have cycle first
    wide <- wide[, c(ncol(wide), seq_len(ncol(wide) - 1))]
    return(wide)
  } else {
    # Long format using the helper
    .pivot_cycle_summary_long(df = wide, cycle_ids = cycles)
  }
}


#' Summarise experiment metadata
#'
#' Prints a summary of key experimental details including chambers, fish masses,
#' respirometer volumes, cycle timing, and environmental conditions.
#'
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return Invisibly returns NULL. Called for its side effect of printing to the console.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                              package = "flatheadresp")
#' summarise_experiment(path = exp_dir_path)
summarise_experiment <- function(path) {
  chambers <- get_chambers(path)
  meta <- get_exp_metadata(path)

  ## --- CLI Header ---
  cli::cli_h2("Experiment Summary")

  ## --- Experiment start ---
  cli::cli_text("Started: {.strong {as.character(meta$exp_start[1])}}")

  ## --- Chambers ---
  cli::cli_text("Number of chambers: {.strong {length(chambers)}}")

  ## --- Environmental conditions ---
  salinity <- meta$Salinity[1]
  temperature <- meta$Temperature[1]

  cli::cli_div(theme = list(
    span.sal = list(color = "cyan"),
    span.temp = list(color = "magenta")
  ))
  cli::cli_h3("Environmental Conditions")
  cli::cli_text("Salinity: {.sal {salinity}} ppt")
  cli::cli_text("Temperature: {.temp {temperature}} deg C")
  cli::cli_end()

  ## --- Fish masses ---
  cli::cli_h3("Fish Mass by Chamber")
  for (i in seq_along(chambers)) {
    chamber <- chambers[i]
    mass <- meta$`Mass.of.fish..kg`[i]
    cli::cli_text("Chamber {.strong {chamber}}: {mass} kg")
  }

  ## --- Respirometer volumes ---
  volumes <- meta$`Volume.respirometer..L`
  cli::cli_h3("Respirometer Volumes")
  if (length(unique(volumes)) == 1) {
    cli::cli_text("All chambers: {.strong {volumes[1]}} L")
  } else {
    for (i in seq_along(chambers)) {
      chamber <- chambers[i]
      vol <- volumes[i]
      cli::cli_text("Chamber {.strong {chamber}}: {vol} L")
    }
  }

  ## --- Number of cycles ---
  n_cycles <- get_n_cycles(path)

  ## --- Cycle timing ---
  flush_time <- meta$`Flush.time..s`[1]
  wait_time <- meta$`Wait.time..s`[1]
  measure_time <- meta$`Measurement.time..s`[1]

  cli::cli_h3("Cycle Timing")
  cli::cli_text("Flush: {flush_time} s")
  cli::cli_text("Wait: {wait_time} s")
  cli::cli_text("Measurement: {measure_time} s")
  cat("\n")
  cli::cli_text("{.strong {n_cycles}} total cycles.")

  invisible(NULL)
}


#' Plot PO2 values for a specific cycle
#'
#' Generates a scatter plot of \eqn{PO_{2}} values over time for each chamber in a given cycle.
#'
#' @param cycle_number Numeric. The cycle number to plot.
#' @param path Character. The path to the AquaResp experiment directory.
#' @param ylim Numeric vector. The y axis plot limits, if NULL (default) it will be the range of po2 values across the cycle. May need to change to c(0, 100) or e.g. c(80,100) to zoom enough to see if there were any 0s or missing data (-300) during the cycle.
#'
#' @return A base-R plot of time vs \eqn{PO_{2}}.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' plot_cycle_po2(cycle_number = 1, path = exp_dir_path)

plot_cycle_po2 <- function(cycle_number, path, ylim = NULL) {

  the_cycle <- read_cycle(cycle_number, path)
  color_palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728")
  po2cols <- grep("^ch[0-9]+\\.po2$", names(the_cycle))
  if(is.null(ylim))
    ylim <- range(the_cycle[po2cols]) +
    c(-0.05, 0.05) * diff(range(the_cycle[po2cols]))
  plot(the_cycle$Unix.Time, the_cycle[[po2cols[1]]],
       col = color_palette[1], ylim = ylim)
  for (i in 2:length(po2cols)) {
    graphics::points(the_cycle$Unix.Time, the_cycle[[po2cols[i]]],
                     col = color_palette[i])
  }
}

#' Get the number of AquaResp cycles in an experiment
#'
#' Counts the number of cycle files (named like \code{Cycle_<number>.txt}) in the
#' \code{"All slopes"} folder of an AquaResp experiment directory.
#'
#' @param path Character. The path to the AquaResp experiment directory.
#'
#' @return Integer. The number of cycles found.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' get_n_cycles(exp_dir_path)
get_n_cycles <- function(path) {
  slopes_dir <- file.path(path, "All slopes")
  if (!dir.exists(slopes_dir)) {
    warning(sprintf("Directory '%s' not found. Returning 0.", slopes_dir))
    return(0L)
  }

  all_files <- list.files(slopes_dir, full.names = FALSE)

  # Expected pattern: Cycle_<number>.txt
  pattern <- "^Cycle_[0-9]+\\.txt$"
  matching_files <- grep(pattern, all_files, value = TRUE)
  non_matching_files <- setdiff(all_files, matching_files)

  if (length(non_matching_files) > 0) {
    warning(sprintf(
      "Found %d file(s) in '%s' that do not match the expected pattern '%s': %s",
      length(non_matching_files),
      slopes_dir,
      pattern,
      paste(non_matching_files, collapse = ", ")
    ))
  }

  if (length(matching_files) == 0) return(0L)

  # Extract numeric cycle indices safely
  cycles <- as.integer(gsub("^Cycle_([0-9]+)\\.txt$", "\\1", matching_files))
  # Drop any NA (shouldn't occur with the strict pattern, but be defensive)
  length(stats::na.omit(cycles))
}


#' Get the start time of a cycle (local time of the first oxygen measurement recorded)
#'
#' @param cycle Data frame of cycle data as returned by read_cycle().
#' @keywords internal
#' @return POSIXct datetime
#'
.cycle_start_dt <- function(cycle) {
  # Validate input
  if (!"Unix.Time" %in% names(cycle)) {
    stop("The input must have a 'Unix.Time' column.")
  }

  # Coerce to numeric
  times <- as.numeric(cycle$Unix.Time)
  if (any(is.na(times))) {
    stop("'Unix.Time' contains non-numeric or NA values.")
  }

  # Check plausibility: Unix timestamps should be > 0 and < far future
  # 2008 - 2030
  if (times[1] < 1e9 || times[1] > 1.9e9) {
    warning("First Unix.Time value seems implausible as a timestamp.")
  }

  # Check if first element is the earliest
  if (times[1] != min(times)) {
    warning("First Unix.Time is not the earliest timestamp in the cycle.")
  }

  # Return as POSIXct
  as.POSIXct(times[1], origin = "1970-01-01")
}

#' Get start times for multiple cycles
#'
#' Returns the start datetime (local time of the first oxygen measurement recorded)
#' for each cycle in a list of cycles.
#'
#' @param cycles A list of cycle data frames as returned by `read_cycle()`.
#'
#' @return A POSIXct vector of start datetimes, named by cycle index.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment",
#'                             package = "flatheadresp")
#' cycles <- lapply(1:3, read_cycle, path = exp_dir_path)
#' flatheadresp:::.get_cycle_starts(cycles)
.get_cycle_starts <- function(cycles) {
  if (!is.list(cycles)) {
    stop("`cycles` must be a list of cycle data frames.")
  }

  starts <- vapply(seq_along(cycles), function(i) {
    cycle <- cycles[[i]]
    if (!"Unix.Time" %in% names(cycle)) {
      stop(sprintf("Cycle %d does not have a 'Unix.Time' column.", i))
    }

    times <- as.numeric(cycle$Unix.Time)
    if (any(is.na(times))) {
      stop(sprintf("Cycle %d has non-numeric or NA values in 'Unix.Time'.", i))
    }

    # Plausibility check (2008-2030)
    if (times[1] < 1e9 || times[1] > 1.9e9) {
      warning(sprintf("Cycle %d: First Unix.Time value seems implausible.", i), call. = FALSE)
    }

    # Warn if first element is not earliest
    if (times[1] != min(times)) {
      warning(sprintf("Cycle %d: First Unix.Time is not the earliest timestamp.", i), call. = FALSE)
    }

    as.POSIXct(times[1], origin = "1970-01-01", tz = "UTC")
  }, FUN.VALUE = as.POSIXct(NA))

  names(starts) <- paste0("cycle_", seq_along(cycles))
  starts
}

