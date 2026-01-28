#' Data for applying Q_10 correction to M_O2 values when temperature was not controlled.
#'
#' For use with `correct_for_temperature()`. This was recorded on a datalogger in
#' the TAF flathead respirometery water jacket during several trials, most of which had a chiller
#' set to 16 degrees C, but two of which were only temperature controlled with air AC
#' and the temperature increased across the course of the trials by ~3 degrees C
#' from about 15.5 to 18.5 if the author recalls correctly.
#'
#' The change is very similar for those two trials, so for the first few trials that
#' also did not have a chiller controlling temperature to 16 degrees these data can be matched
#' to estimate the temp across points in the trial and a Q_10 correction applied.
#'
#' The data is created by GAM- smoothing two trials for which temps were not chiller controlled,
#' and then extrapolating out to 25 hours in case some trials went longer than these two.
#'
#' @format A data frame with 300 rows and 2 variables:
#' \describe{
#'   \item{TIME.HOURS}{Decimal hours since experiment began}
#'   \item{est_temp}{Temperature in degrees Celsius}
#' }
#' @source Derived from internal data-raw pipeline.
#' @examples
#' data(temp_correct_data)
#' head(temp_correct_data)
"temp_correct_data"
