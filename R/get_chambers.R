
#' Get the active chambers from an AquaResp experiment directory
#'
#' This function scans the experiment directory and returns the chamber numbers
#' for which summary data files are available. These correspond to the active chambers
#' used during the AquaResp experiment.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' This should be the folder created by AquaResp that contains folders like "All slopes",
#' "Experimental information", and "Oxygen data raw".
#'
#' @return A numeric vector of active chamber numbers.
#' @export
#'
#' @examples
#' \dontrun{
#' get_chambers(path = "~/data/Resp/Experiment_1")
#' }

get_chambers <- function(path) {
  gsub("Summary data resp |.txt",
       "",
       list.files(path, pattern = "Summary data resp .*\\.txt")) |>
    as.numeric() |>
    sort()

}

#
#' Get metadata for a specific chamber from an AquaResp experiment
#'
#' This function reads the first lines of the summary data file for a given chamber
#' and extracts experimental metadata such as fish mass, respirometer volume, and start time.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' @param chamber Numeric. The chamber number to retrieve metadata for. Use `get_chambers()` to list available chambers.
#'
#' @return A named list containing metadata for the specified chamber.
#' @export
#'
#' @examples
#' \dontrun{
#' get_chamber_metadata(path = "~/data/Resp/Experiment_1", chamber = 1)
#' }
get_chamber_metadata <- function(path, chamber) {

  sum_metadata <- list()

  sum_path <- file.path(path, paste0("Summary data resp ", chamber, ".txt"))

  con <- file(sum_path, "r")
  for (j in 1:15) {
    line <- readLines(con, n = 1)
    if (grepl(":", line)) {
      parts <- strsplit(line, ":")[[1]]
      key <- trimws(parts[1])
      value <- trimws(parts[2]) |> as.numeric() |> round(4)
      sum_metadata[[key]] <- value
    }
  }
  close(con)

  sum_metadata$exp_start <-
    sum_metadata$`Experiment start, UNIX time` |> as.POSIXct()

  sum_metadata <- c(list(chamber = chamber), sum_metadata)
  return(sum_metadata)
}

#' Get metadata for all chambers in an AquaResp experiment
#'
#' This function retrieves metadata for all active chambers in the experiment directory
#' and returns it as a dataframe, with one row per chamber.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#'
#' @return A dataframe where each row contains metadata for one chamber.
#' @export
#'
#' @examples
#' \dontrun{
#' get_exp_metadata(path = "~/data/Resp/Experiment_1")
#' }

get_exp_metadata <- function(path) {
  chambers <- get_chambers(path)

  metadata_list <-
    lapply(chambers, function(i) {
      get_chamber_metadata(path, chamber = i)

    })

  do.call(rbind, lapply(metadata_list, as.data.frame))
}


#' Get MO2 values from an AquaResp experiment
#'
#' This function reads the summary data files and returns mass-specific oxygen consumption (MO2)
#' values for one or more chambers. If no chamber is specified, data from all chambers is returned.
#'
#' @param path Character. The directory location of the AquaResp experiment on disk.
#' @param chamber Optional numeric. The chamber number(s) to retrieve MO2 values for. If omitted, all chambers are returned.
#'
#' @return A dataframe with columns for chamber number, cycle number, and MO2 values.
#' @export
#'
#' @examples
#' \dontrun{
#' get_exp_MO2s(path = "~/data/Resp/Experiment_1") # All chambers
#' get_exp_MO2s(path = "~/data/Resp/Experiment_1", chamber = 1) # Chamber 1 only
#' }

get_exp_MO2s <- function(path, chamber = NULL) {

  chambers <- get_chambers(path)

  if (is.null(chamber)) {
    summary_data_list <-
      lapply(chambers, function(i) {
        sum_path <-
          file.path(path, paste0("Summary data resp ", i, ".txt"))

        sum_data <-
          utils::read.table(sum_path,
                     skip = 15,
                     sep = ";",
                     header = TRUE) |>
          (\(res) res[, setdiff(names(res), "X"), drop = FALSE])()

        return(data.frame(chamber = i, cycle = 1:nrow(sum_data), sum_data))
      })

  } else {
    if (!all(chamber %in% chambers))
      stop(paste("chamber should be one of", paste(chambers, collapse = ", ")))

    summary_data_list <-
      lapply(chamber, function(i) {
        sum_path <-
          file.path(path, paste0("Summary data resp ", i, ".txt"))

        sum_data <-
          utils::read.table(sum_path,
                     skip = 15,
                     sep = ";",
                     header = TRUE) |>
          (\(res) res[, setdiff(names(res), "X"), drop = FALSE])()


        return(data.frame(chamber = i, cycle = 1:nrow(sum_data), sum_data))
      })
  }
  do.call('rbind', summary_data_list)
}
