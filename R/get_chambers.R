get_chambers <- function(path) {

  gsub("Summary data resp |.txt", "",
       list.files(path, pattern = "Summary data resp .*\\.txt")) |>
    as.numeric() |>
    sort()

}
