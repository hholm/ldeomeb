#' @importFrom magrittr %>%
NULL

#' Import and format plate reader data from TECAN i-control.
#'
#' This function returns a data.frame containing formatted data from a from TECAN i-control csv output.
#'
#' @param loc A [character] string containing the file location of raw plate reader data to import.
#' @return A [data.frame] containing imported and formatted data.
#' @export
tidyplate <- function(loc) {
  dat <- read.csv(loc)

  # find plate labels
  labels <- dat[which(stringr::str_detect(dat[, 1], "Label: ")), 1]

  # find plate sizes
  plate.sizes <- list(
    rows = which(dat[, 1] == "End Time:") - which(dat[, 1] == "Start Time:") - 8,
    cols = as.vector(apply(dat[which(dat[, 1] == "Start Time:")+3, ], 1, function(x) {
      length(which(!is.na(suppressWarnings(as.numeric(x)))))
    }))
  )

  # run for loop to gather data for each plate and format
  plates <- data.frame()
  for (i in 1:length(labels)) {
    label <- labels[i]

    # read plate values
    plate <- data.frame(dat[(which(dat[, 1] == "Start Time:")[i] + 4):((which(dat[, 1] == "Start Time:")[i]+3) + plate.sizes$rows[i]), 2:(2 + plate.sizes$cols[i] - 1)]) %>%
      dplyr::mutate(across(everything(),as.numeric))
    colnames(plate) <- dat[(which(dat[, 1] == "Start Time:")[i]+3) , ][2:(2 + (plate.sizes$cols[i] - 1))]
    plate$rows <- dat[((which(dat[, 1] == "Start Time:")[i]+3)  + 1):((which(dat[, 1] == "Start Time:")[i]+3)  + plate.sizes$rows[i]), 1]

    # add in metadata
    meta.data <- data.frame(setting = dat[(which(dat[, 1] == label)+1):(which(dat[, 1] == "Start Time:")[i]-1), 1],
                            val = dat[(which(dat[, 1] == label)+1):(which(dat[, 1] == "Start Time:")[i]-1), 5])

    meta.data <- dat[(which(dat[, 1] == label)+1):(which(dat[, 1] == "Start Time:")[i]-1), 5]
    names(meta.data) <- dat[(which(dat[, 1] == label)+1):(which(dat[, 1] == "Start Time:")[i]-1), 1]
    
    if (nrow(plates) == 0) {
      plates <- plate %>%
        mutate(across(-rows,as.numeric)) %>%
        tidyr::pivot_longer(cols = -rows, names_to = "cols", values_to = "value") %>%
        mutate(!!!meta.data) %>%
        mutate(`Start Time` = dat[(which(dat[, 1] == "Start Time:")[i]),2]) %>%
        rbind()
    }else{
    plates <- plate %>%
      mutate(across(-rows,as.numeric)) %>%
      tidyr::pivot_longer(cols = -rows, names_to = "cols", values_to = "value") %>%
      mutate(!!!meta.data) %>%
      mutate(`Start Time` = dat[(which(dat[, 1] == "Start Time:")[i]),2]) %>%
      dplyr::full_join(plates)
    }
  }

  # add a name column if it exists
  if (any(dat[, 1] == "Names")) {
    name.frame <- dat[which(dat[, 1] == "Names"):(which(dat[, 1] == "Names") + plate.sizes$rows[1]), 1:(plate.sizes$cols[1] + 1)] %>%
      purrr::set_names(.[1, ]) %>%
      dplyr::slice(-1) %>%
      tidyr::pivot_longer(cols = !Names, names_to = "cols", values_to = "names")

    name.frame <- do.call(rbind, replicate(length(labels), name.frame, simplify = FALSE))

    # check alignment
    if (identical(name.frame$cols, plates$cols) & identical(name.frame$Names, plates$rows)) {
      plates$names <- name.frame$names
    }
  }
  return(plates)
}

#' Calculate pH values from m-cresol dye absorbance.
#'
#' This function returns a data.frame containing calculated pH values.
#'
#' @details
#' * A730_blank A [vector] or [numeric] of samples absorbance at 730nm with no dye.
#' * A578_blank A [vector] or [numeric] of samples absorbance at 578nm with no dye.
#' * A434_blank A [vector] or [numeric] of samples absorbance at 434nm with no dye.
#' * A730_dye A [vector] or [numeric] of samples absorbance at 730nm with m-cresol dye.
#' * A578_dye A [vector] or [numeric] of samples absorbance at 578nm with m-cresol dye.
#' * A434_dye A [vector] or [numeric] of samples absorbance at 434nm with m-cresol dye.
#' * vol.dye.L Volumn of dye added in liters (defaults to 0.01).
#' * salinity Sample salinity in PSU (defaults to 35).
#' * verbose Should only the pH values (False) be returned or metadata as well (True)?
#' @md
#'
#' @return A [data.frame] containing calculated pH values with metadata, or a [vector] containing just pH values if verbose == FALSE.
#' @export
calc_pH_spec <- function(A730_blank, A578_blank, A434_blank, A730_dye, A578_dye, A434_dye, vol.dye.L, salinity, verbose = TRUE) {
  # following hennon formula
  A1_A2 <- (A578_dye - A578_blank - (A730_dye - A730_blank)) / (A434_dye - A434_blank - (A730_dye - A730_blank))
  pK2 <- (1245.9 / 298) + 3.8275 + (0.00211 * (35 - salinity))
  #A1_A2_cor <- A1_A2 + (0.0218 - (0.0359 * A1_A2)) * vol.dye.L
  A1_A2_cor <- A1_A2 + (-0.00660453* vol.dye.L*1e6)
  # return all data if desired
  if (verbose) {
    return(data.frame(A1_A2, pK2, A1_A2_cor, pH = pK2 + log10((A1_A2_cor - 0.00691) / (2.222 - A1_A2_cor * 0.1331))))
  } else {
    return(pK2 + log10((A1_A2_cor - 0.00691) / (2.222 - A1_A2_cor * 0.1331)))
  }
}

#' Calculate pH values from m-cresol dye absorbance for an entire plate at once.
#'
#' This function returns a data.frame containing calculated pH values and metadata.
#'
#' @details
#' * tidyplate A [data.frame] containing formatted absorbance created using [tidyplate()].
#' * pairs A [data.frame] containing paired 'blank' and 'dye' rows. Defaults too: [data.frame(blanks = c("A", "C", "E", "G"), dyes = c("B", "D", "F", "H"))]
#' * vol.dye.L Volumn of dye added in liters (defaults to 0.01).
#' * salinity Sample salinity in PSU (defaults to 35).
#' * verbose Should only the pH values (False) be returned or metadata as well (True)?
#' @md
#'
#' @return A [data.frame] containing calculated pH values with metadata.
#' @export
calc_plate <- function(tidyplate, verbose = TRUE, pairs = data.frame(blanks = c("A", "C", "E", "G"), dyes = c("B", "D", "F", "H")),
                       salinty = 35, vol.dye.L = 0.01) {

  tidyplate <- tidyplate %>% subset(Mode == "Absorbance")

  # assign well type based on 'pairs" data.frame
  tidyplate[which(tidyplate$rows %in% pairs$blanks), "type"] <- "blank"
  tidyplate[which(tidyplate$rows %in% pairs$dyes), "type"] <- "dye"

  # add a position column
  tidyplate$position <- tidyplate$cols

  out <- list()
  for (i in 1:nrow(pairs)) {
    out[[i]] <- tidyplate %>%
      subset(tidyplate$rows %in% pairs[i, ]) %>%
      dplyr::mutate(position = paste0(position,paste0(pairs[i, ], collapse = "")))
  }

  tidyplate <- do.call(rbind, out)

  # if sample names exist, keep them, if not use only positions
  if (!suppressWarnings(is.null(tidyplate$names))) {
    # fix if names don't match across position
    tidyplate <- tidyplate %>%
      dplyr::group_by(position) %>%
      dplyr::mutate(names = list(unique(names))) %>%
      dplyr::ungroup()
    id_cols <- c("position", "names")
  } else {
    id_cols <- "position"
    #warning("No names column in data.frame. Using well positions as names.")
  }

  # time issues
  if (all(stringr::str_detect(tidyplate$`Start Time`, pattern = "AM|PM"))) {
    return <- tidyplate %>%
      mutate(`Start Time` = lubridate::parse_date_time(`Start Time`, orders = "%m/%d/%Y %h:%M:%S %p"))
  } else {
    return <- tidyplate %>%
      mutate(`Start Time` = lubridate::mdy_hm(`Start Time`))
  }

  # reformat
  return <- return %>%
    mutate(`Average Start Time` = mean(`Start Time`)) %>%
    select(-`Start Time`) %>%
    select(-c("rows","cols")) %>%
    tidyr::pivot_wider(names_from = c(type,`Measurement Wavelength`), values_from = value) %>%
    dplyr::mutate_at(dplyr::vars("blank_730", "dye_730", "blank_578", "dye_578", "blank_434", "dye_434"), as.numeric)

  # calculate pH
  pH <- calc_pH_spec(
    A730_blank = return$blank_730, A578_blank = return$blank_578, A434_blank = return$blank_434,
    A730_dye = return$dye_730, A578_dye = return$dye_578, A434_dye = return$dye_434,
    vol.dye.L = vol.dye.L, salinity = salinty,verbose = verbose
  )

  # return just pH column or all metadata
  if (verbose) {
    return <- cbind(return, pH)
  } else {
    return$pH <- pH
    return <- return[,c(id_cols,"pH")]
  }
  return(return)
}


