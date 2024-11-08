#' Nicely formats plate reader data.
#'
#' @param loc Location of raw plate reader data.
#' @export
tidyplate <- function(loc) {

  dat <- read.csv(loc)

  # find plate labels
  labels <- dat[which(stringr::str_detect(dat[, 1], "Label: ")), 1]

  # find plate sizes
  plate.sizes <- list(
    rows = which(dat[, 1] == "End Time:") - which(dat[, 1] == "<>") - 5,
    cols = as.vector(apply(dat[which(dat[, 1] == "<>"), ], 1, function(x) {
      length(which(!is.na(suppressWarnings(as.numeric(x)))))
    }))
  )

  # run for loop to gather data for each plate and format
  plates <- data.frame()
  for (i in 1:length(labels)) {
    label <- labels[i]

    #read plate values
    plate <- data.frame(dat[(which(dat[, 1] == "<>")[i] + 1):(which(dat[, 1] == "<>")[i] + plate.sizes$rows[i]), 2:(2 + plate.sizes$cols[i] - 1)])
    colnames(plate) <- dat[which(dat[, 1] == "<>")[i], ][2:(2 + (plate.sizes$cols[i] - 1))]
    plate$rows <- dat[(which(dat[, 1] == "<>")[i] + 1):(which(dat[, 1] == "<>")[i] + plate.sizes$rows[i]), 1]

    #add in metadata
    plates <- tidyr::pivot_longer(plate, -rows, names_to = "cols", values_to = dat[which(dat[, 1] == label) + 1, 5]) |>
      dplyr::mutate(
        start.datetime = dat[which(dat[, 1] == label) + 7, 2],
        end.datetime = dat[which(dat[, 1] == label) + plate.sizes$rows[i] + 15, 2],
        wave.length = as.numeric(dat[which(dat[, 1] == label) + 2, 5]),
        bandwidth = as.numeric(dat[which(dat[, 1] == label) + 3, 5]),
        n.flash = as.numeric(dat[which(dat[, 1] == label) + 4, 5]),
        temp = as.numeric(stringr::str_extract(dat[which(dat[, 1] == label) + 9, 2], pattern = "[:digit:]..."))
      ) |>
      rbind(plates)
  }
  return(plates)
}
