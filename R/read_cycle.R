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

#' Calculate MO2 values for an AquaResp cycle
#'
#' Calculates mass-specific oxygen consumption (\eqn{MO_{2}}) for each chamber
#' in a given cycle using linear regression on \eqn{PO_{2}} data and experiment-level metadata.
#'
#' Accepts either a **numeric cycle ID** (e.g., 1) or a **cycle dataframe** returned by `read_cycle()`.
#' `path` is always required to access metadata.
#'
#' @param cycle Numeric or data.frame.
#'   - If **numeric**: the cycle number to read (e.g., 1, 2, 3, ...).
#'   - If **data.frame**: a cycle dataframe as returned by `read_cycle()`.
#' @param path Character. Path to the AquaResp experiment directory.
#' @param chambers Optional integer vector: subset of chambers to include,
#'   e.g. `c(1,2)`; or **negative** to exclude, e.g. `-4` (means keep 1,2,3).
#'   Default `NULL` (all chambers).
#' @param long Logical. If `TRUE`, return a tibble in long format
#'   with columns `cycle`, `chamber`, and `mo2`. Default `FALSE` (wide named vector).
#'
#' @return Either a named numeric vector (wide, default) with names like `"ch1.mo2"`,
#'   or a tibble (long) with columns `cycle`, `chamber`, `mo2`.
#' @export
#'
#' @examples
#' exp_dir_path <- system.file("extdata", "aquaresp_experiment", package = "flatheadresp")
#' # Using a cycle number (wide):
#' calc_cycle_mo2s(1, path = exp_dir_path)
#' # Using pre-read cycle data (long, chambers 1 & 2):
#' cy <- read_cycle(1, exp_dir_path)
#' calc_cycle_mo2s(cy, path = exp_dir_path, chambers = c(1,2), long = TRUE)
calc_cycle_mo2s <- function(cycle, path, chambers = NULL, long = FALSE) {
  # Validate path
  if (missing(path) || !is.character(path) || length(path) != 1L) {
    stop("`path` must be a single character string pointing to the experiment directory.")
  }

  # Resolve cycle data + cycle_id (if known)
  if (is.numeric(cycle) && length(cycle) == 1L && !is.na(cycle)) {
    cycle_id   <- as.integer(cycle)
    cycle_data <- read_cycle(cycle_id, path)
  } else if (is.data.frame(cycle)) {
    cycle_id   <- NA_integer_
    cycle_data <- cycle
  } else {
    stop("`cycle` must be either a single numeric cycle number or a cycle data.frame.")
  }

  # Basic checks
  if (!"Unix.Time" %in% names(cycle_data)) {
    stop("`cycle` data is missing required column `Unix.Time`.")
  }
  po2_cols <- grep("po2", names(cycle_data))
  if (length(po2_cols) == 0) stop("`cycle` data contains no 'po2' columns.")

  # Experiment level metadata
  all_chambers <- as.integer(get_chambers(path))
  sel_chambers <- .resolve_chambers(all_chambers, chambers)
  meta         <- get_exp_metadata(path)

  beta      <- meta$`Oxygen.solubilty..mg.O2...L`
  rRespFish <- meta$`Real.volume..vresp...vfish...neutrally.bouyant...L` /
    meta$`Mass.of.fish..kg`

  if (length(beta) != length(all_chambers) || length(rRespFish) != length(all_chambers)) {
    stop("Metadata lengths do not match number of chambers.")
  }

  # Time covariate: seconds since cycle start
  cycle_time_x <- cycle_data$Unix.Time - min(cycle_data$Unix.Time) + 1

  # Regress PO2 ~ time for each chamber; collect slopes
  cycle_slope_list <- lapply(po2_cols, function(i) {
    fit_sum <- stats::lm(cycle_data[[i]] ~ cycle_time_x) |> summary()
    fit_sum$coefficients[2, ]  # slope row (estimate, etc.)
  })

  # Map selected chambers to MO2
  mo2 <- sapply(sel_chambers, \(i) {
    slope <- cycle_slope_list[[i]][1]
    -1 * (slope / 100) * beta[i] * rRespFish[i] * 3600.0
  })

  if (!long) {
    # Wide: name as "chN.mo2"
    names(mo2) <- paste0("ch", sel_chambers, ".mo2")
    return(mo2)
  } else {
    # Long tibble
    tibble::tibble(
      cycle   = cycle_id,
      chamber = sel_chambers,
      mo2     = unname(mo2)
    )
  }
}


#' Get min and max PO2 values for each chamber in a cycle
#'
#' Extracts the minimum and maximum non-zero \eqn{PO_{2}} values for each chamber from a cycle dataframe.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of min and max \eqn{PO_{2}} values for each chamber.
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

    unlist(lapply(cycle[po2_cols], \(x) {
      c(min = min(x[x != 0]), max = max(x[x != 0]))
    }))
  }


#' Calculate percentage of zero PO2 values in a cycle
#'
#' Computes the percentage of zero values in each \eqn{PO_{2}} column of a cycle dataframe.
#' Values greater than zero would have an incorrect \eqn{R^{2}} as calculated by AquaResp.
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @encoding UTF-8
#'
#' @return A named vector of percentages (0 to 100%) indicating missingness per chamber.
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
      round(length(which(x == 0)) / length(x), 4) * 100) |>
      stats::setNames(gsub(".po2", "", paste0(names(cycle)[po2_cols], ".pct0"))) |>
      unlist()
  }

#' Calculate correlation between PO2 and time for each chamber
#'
#' Performs Pearson correlation (r) tests between \eqn{PO_{2}} values and Unix time for each chamber in a cycle, retyrns r ^ 2
#'
#' @param cycle Dataframe. A cycle dataframe as returned by `read_cycle()`.
#'
#' @return A named numeric vector of correlation coefficients and p-values for each chamber.
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
      }, secs = cycle$Unix.Time)

    r2s <- lapply(cors, function(x) {
      x$estimate^2 |> unname()
      }) |> #returns r ^ 2
      stats::setNames(gsub(".po2","",paste0(names(cycle[po2_cols]), ".r2")))

    r2_pvals <- lapply(cors, function(x) {
    x$p.value
    }) |>
      stats::setNames(gsub(".po2","",paste0(names(cycle[po2_cols]), ".p")))

  unlist(c(r2s, r2_pvals))
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

#' Summarize all cycles in an AquaResp experiment
#'
#' Reads every cycle under `path/All slopes`, summarizes MO2, PO2 min/max,
#' R^2, R^2 p-value, and % missing, with optional chamber subsetting and
#' long-format output.
#'
#' @param path Character. Experiment directory.
#' @param chambers Optional integer vector: subset of chambers to include
#'   (e.g. `c(1,2)`), or negative to exclude (e.g. `-4` means keep 1,2,3).
#'   Default `NULL` (all chambers).
#' @param long Logical. If `TRUE`, returns a long tibble; otherwise wide data frame.
#'
#' @return A data.frame (wide, default) or tibble (long) summarizing all cycles.
#' @export
cycle_summary <- function(path, chambers = NULL, long = FALSE) {
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
      mo2   <- calc_cycle_mo2s(cycle_data, path, chambers = sel_chambers, long = FALSE)
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

  rownames(wide) <- NULL

  if (!long) {
    # Attach cycle IDs as a column for convenience
    wide$cycle <- cycles
    # Reorder to have cycle first
    wide <- wide[, c(ncol(wide), seq_len(ncol(wide) - 1))]
    return(wide)
  } else {
    # Long format using the helper
    pivot_cycle_summary_long(wide, cycle_ids = cycles)
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
#' @export
pivot_cycle_summary_long <- function(df, cycle_ids = NULL) {
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

  cat("\n", "Experiment started",
      as.character(meta$exp_start[1]), "\n\n")
  cat(length(chambers), "chambers\n\n")

  # Salinity and temperature
  salinity <- meta$Salinity[1]
  temperature <- meta$Temperature[1]

  cat("Environmental conditions:\n")
  cat("  Salinity:", salinity, "ppt\n")
  cat("  Temperature:", temperature, "\u00B0C\n")

  # Fish masses
  cat("Mass of fish by chamber:\n")
  for (i in seq_along(chambers)) {
    chamber <- chambers[i]
    mass <- meta$`Mass.of.fish..kg`[i]
    cat("  Chamber", chamber, ":", mass, "kg\n")
  }

  # Respirometer volumes
  volumes <- meta$`Volume.respirometer..L`
  if (length(unique(volumes)) == 1) {
    cat("\nRespirometer volume is", volumes[1], "L for all chambers.\n")
  } else {
    cat("\nRespirometer volume by chamber:\n")
    for (i in seq_along(chambers)) {
      chamber <- chambers[i]
      vol <- volumes[i]
      cat("  Chamber", chamber, ":", vol, "L\n")
    }
  }

  # Cycle timing
  flush_time <- meta$`Flush.time..s`[1]
  wait_time <- meta$`Wait.time..s`[1]
  measure_time <- meta$`Measurement.time..s`[1]

  cat("\nCycle timing:\n")
  cat("  Flush time:", flush_time, "seconds\n")
  cat("  Wait time:", wait_time, "seconds\n")
  cat("  Measurement time:", measure_time, "seconds\n")

  # Number of cycles
  n_cycles <- get_n_cycles(path)
  cat("\n", n_cycles, "cycles")

  # # data QA/QC metrics
  # cycle_summary <- cycle_summary(path) |> pivot_cycle_summary_long()

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
#'
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
#' get_cycle_starts(cycles)
get_cycle_starts <- function(cycles) {
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

    # Plausibility check (2008–2030)
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

