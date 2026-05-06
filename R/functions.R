#' @importFrom magrittr %>%
NULL

#' Import and format plate reader data from TECAN i-control.
#'
#' This function returns a data.frame containing formatted data from a from TECAN i-control csv output.
#'
#' @param loc A [character] string containing the file location of raw plate reader data to import.
#' @return A [data.frame] containing imported and formatted data.
#' @export
tidyplate <- function(dat) {
  # find plate labels
  labels <- dat[which(stringr::str_detect(dat[, 1], "Label: ")), 1]

  # find plate sizes
  plate.sizes <- list(
    rows = which(dat[, 1] == "End Time:") - which(dat[, 1] == "Start Time:") - 8,
    cols = as.vector(apply(dat[which(dat[, 1] == "Start Time:") + 3, ], 1, function(x) {
      length(which(!is.na(suppressWarnings(as.numeric(x)))))
    }))
  )

  # run for loop to gather data for each plate and format
  plates <- data.frame()
  for (i in 1:length(labels)) {
    label <- labels[i]

    # read plate values
    plate <- data.frame(dat[(which(dat[, 1] == "Start Time:")[i] + 4):((which(dat[, 1] == "Start Time:")[i] + 3) + plate.sizes$rows[i]), 2:(2 + plate.sizes$cols[i] - 1)]) %>%
      dplyr::mutate(across(everything(), as.numeric))
    colnames(plate) <- dat[(which(dat[, 1] == "Start Time:")[i] + 3), ][2:(2 + (plate.sizes$cols[i] - 1))]
    plate$rows <- dat[((which(dat[, 1] == "Start Time:")[i] + 3) + 1):((which(dat[, 1] == "Start Time:")[i] + 3) + plate.sizes$rows[i]), 1]

    # add in metadata
    meta.data <- data.frame(
      setting = dat[(which(dat[, 1] == label) + 1):(which(dat[, 1] == "Start Time:")[i] - 1), 1],
      val = dat[(which(dat[, 1] == label) + 1):(which(dat[, 1] == "Start Time:")[i] - 1), 5]
    )

    meta.data <- dat[(which(dat[, 1] == label) + 1):(which(dat[, 1] == "Start Time:")[i] - 1), 5]
    names(meta.data) <- dat[(which(dat[, 1] == label) + 1):(which(dat[, 1] == "Start Time:")[i] - 1), 1]

    if (nrow(plates) == 0) {
      plates <- plate %>%
        mutate(across(-rows, as.numeric)) %>%
        tidyr::pivot_longer(cols = -rows, names_to = "cols", values_to = "value") %>%
        mutate(!!!meta.data) %>%
        mutate(`Start Time` = dat[(which(dat[, 1] == "Start Time:")[i]), 2]) %>%
        rbind()
    } else {
      plates <- plate %>%
        mutate(across(-rows, as.numeric)) %>%
        tidyr::pivot_longer(cols = -rows, names_to = "cols", values_to = "value") %>%
        mutate(!!!meta.data) %>%
        mutate(`Start Time` = dat[(which(dat[, 1] == "Start Time:")[i]), 2]) %>%
        {suppressMessages(dplyr::full_join(.,plates))}
    }
  }

  # add a name column if it exists
  if (any(dat[, 1] == "Names",na.rm = T)) {
    name.frame <- dat[which(dat[, 1] == "Names"):(which(dat[, 1] == "Names") + plate.sizes$rows[1]), 1:(plate.sizes$cols[1] + 1)] %>%
      purrr::set_names(.[1, ]) %>%
      dplyr::slice(-1) %>%
      tidyr::pivot_longer(cols = !Names, names_to = "cols", values_to = "names")

    name.frame <- do.call(rbind, replicate(length(labels), name.frame, simplify = FALSE)) %>%
      rename(rows = Names) %>%
      distinct()

    plates <- suppressMessages(dplyr::left_join(plates, name.frame, by = c("rows", "cols")))
  }else{
    warning("No sample name data found in sheet.")
  }

  # add a name column if it exists
  if (any(dat[, 1] == "pH_pairs", na.rm = T)) {
    pairs <- dat[which(dat[, 1] == "pH_pairs"):nrow(dat), 1:3] %>%
      purrr::set_names(.[1, ]) %>%
      dplyr::slice(-1) %>%
      dplyr::select(-pH_pairs)

    # assign well type based on 'pairs" data.frame
    plates[which(plates$rows %in% pairs$blanks), "type"] <- "blank"
    plates[which(plates$rows %in% pairs$dyes), "type"] <- "dye"
  } else {
    pairs <- "No pair data found in sheet."
    warning("No pH pair data (i.e. dye/blank) found in sheet.")
  }
  return(list(plates = plates, pairs = pairs))
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
calc_pH_spec <- function(A730_blank, A578_blank, A434_blank, A730_dye, A578_dye,
                         A434_dye, vol.dye.uL, salinity, verbose = TRUE,
                         calibration = calibration) {
  # following hennon formula
  A1_A2 <- (A578_dye - A578_blank - (A730_dye - A730_blank)) / (A434_dye - A434_blank - (A730_dye - A730_blank))
  pK2 <- (1245.9 / 298) + 3.8275 + (0.00211 * (35 - salinity))
  #A1_A2_cor <- A1_A2 + (0.0218 - (0.0359 * A1_A2)) * vol.dye.L
  A1_A2_cor <- A1_A2 + (-1*calibration*vol.dye.uL)
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
calc_plate <- function(tidyplate, verbose = TRUE,
                       pairs = data.frame(blanks = c("A", "C", "E", "G"), dyes = c("B", "D", "F", "H")),
                       salinity = 35, vol.dye.uL = 3, calibration = -0.00660453) {

  # assign well type based on 'pairs" data.frame
  tidyplate[which(tidyplate$rows %in% pairs$blanks), "type"] <- "blank"
  tidyplate[which(tidyplate$rows %in% pairs$dyes), "type"] <- "dye"

  # drop any (assumed) accidental Ft readings on dye filled wells
  tidyplate <- tidyplate[!(tidyplate$type == "dye" & tidyplate$Mode == "Fluorescence Top Reading"),]

  # add a position column
  tidyplate$position <- tidyplate$cols

  if (!all(tidyplate$rows %in% unlist(pairs))) {
    warning("No cells measured not found in 'ph pairs' list.")
  }

  out <- list()
  for (i in 1:nrow(pairs)) {
    out[[i]] <- tidyplate %>%
      subset(tidyplate$rows %in% pairs[i, ]) %>%
      dplyr::mutate(position = paste0(position, paste0(pairs[i, ], collapse = "")))
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
    # warning("No names column in data.frame. Using well positions as names.")
  }

  # time issues
  if (all(stringr::str_detect(tidyplate$`Start Time`, pattern = "AM|PM"))) {
    return <- tidyplate %>%
      mutate(`Start Time` = lubridate::parse_date_time(`Start Time`, orders = "%m/%d/%Y %h:%M:%S %p"))
  } else {
    return <- tidyplate %>%
      mutate(`Start Time` = lubridate::mdy_hm(`Start Time`))
  }

  return <- return %>%
    mutate(`Average Start Time` = mean(`Start Time`)) %>%
    dplyr::select(-`Start Time`) %>%
    dplyr::select(-c("rows", "cols")) %>%
    tidyr::pivot_wider(names_from = Mode, values_from = value)

  if("Fluorescence Top Reading" %in% colnames(return)){
    return <- return %>%
    dplyr::select(c("Excitation Wavelength", "names", "type", "position", "Absorbance",
             "Fluorescence Top Reading", "Measurement Wavelength", "Average Start Time")) %>%
    tidyr::pivot_wider(names_from = `Excitation Wavelength`, values_from = `Fluorescence Top Reading`, names_prefix = "Ft_") %>%
    dplyr::select(!contains("Ft_NA"))
  }else{
    return <- return %>%
      dplyr::select(c( "names", "type", "position", "Absorbance",
               "Measurement Wavelength", "Average Start Time"))
  }

  # reformat
  return <- return %>%
    tidyr::pivot_wider(names_from = c(type, `Measurement Wavelength`), values_from = Absorbance)
  if (class(pull(return, "names")) == "list") {
    return$names <- pull(return, "names") %>% unlist()
  }

  return <- return %>%
    group_by(names, position, `Average Start Time`) %>%
    dplyr::select(!(contains("blank_NA") | contains("dye_NA"))) %>%
    summarise(across(everything(), function(x) {
      max(x, na.rm = T)
    }), .groups = "drop")


  # reformat ft values if they exist into a wide format, might be a simpler was to do this
  #  if(exists("Fluorescence Top Reading",where = return)){

  # get ft columns to pivot that are different accross 455 and 630
  #    return <- return %>%
  #    tidyr::pivot_wider(names_from = `Excitation Wavelength`, values_from = contains("Fluorescence Top Reading") | contains("Gain")  | contains("Z-Position")) %>%
  #    dplyr::select(!contains("_NA"))
  #
  #    # cut and rejoin the ft and abs parts
  #    ft_part <- return %>% subset(!is.na(`Emission Wavelength`))
  #    abs_part <- return %>% subset(is.na(`Emission Wavelength`))#

  #    ft_part <- ft_part %>% dplyr::select(where(~ !all(is.na(.))))
  #    abs_part <- abs_part %>% dplyr::select(where(~ !all(is.na(.))))

  # Identify columns with identical values
  #    common_cols <- intersect(names(ft_part), names(abs_part)) |>
  #      keep(~ identical(ft_part[[.x]], abs_part[[.x]]))

  #    return <- abs_part |>
  #        dplyr::select(-all_of(common_cols[common_cols != 'names'])) %>%
  #        full_join(ft_part,by = join_by(names), suffix = c(" Abs", " Ft"))
  #  }

  pH <- return %>%
    dplyr::select(
      starts_with("blank_7"), starts_with("dye_7"), starts_with("blank_5"),
      starts_with("dye_5"), starts_with("blank_4"), starts_with("dye_4")
    ) %>%
    dplyr::mutate_all(as.numeric)

  # calculate pH
  pH <- calc_pH_spec(
    A730_blank = pull(dplyr::select(pH, starts_with("blank_7"))),
    A578_blank = pull(dplyr::select(pH, starts_with("blank_5"))),
    A434_blank = pull(dplyr::select(pH, starts_with("blank_4"))),
    A730_dye = pull(dplyr::select(pH, starts_with("dye_7"))),
    A578_dye = pull(dplyr::select(pH, starts_with("dye_5"))),
    A434_dye = pull(dplyr::select(pH, starts_with("dye_4"))),
    vol.dye.uL = vol.dye.uL,
    salinity = salinity,
    verbose = verbose,
    calibration = calibration
  )

  # names fix idk why this is happening
  if (class(pull(return, "names")) == "list") {
    return$names <- pull(return, "names") %>% unlist()
  }

  # return just pH column or all metadata
  if (verbose) {
    return <- cbind(return, pH)
  } else {
    return$pH <- pH
    return <- return[, c(id_cols, "pH")]
  }
  return(return)
}

#' Convert fluorescence values between instruments or excitation colors
#'
#' Predicts equivalent fluorescence values across instruments or colour channels
#' using pre-fitted linear models stored in \code{\link{ft_convert_df}}, or
#' user-supplied model parameters. Prediction intervals are propagated via the
#' delta method and can be back-transformed when working on a log10 scale.
#'
#' @param ft_value Numeric vector of fluorescence values to convert.
#' @param strain Character. Phytoplankton strain values come from. Ignored when supplying manual model
#'   parameters. #' Should be one of:
#'   \itemize{
#'     \item \code{"Syn7002"}
#'     \item \code{"Syn8102"}
#'     \item \code{"Ehux"}
#'     \item \code{"Mcomm"}
#'     \item \code{"Tpsu"}
#'   }
#' @param from Character. Name of the input fluorescence (e.g.
#'   \code{"aqua_blue"}). Must be specified together with \code{to}. Should be one of:
#'   \itemize{
#'     \item \code{"aqua_red"}
#'     \item \code{"plate_blue"}
#'     \item \code{"plate_red"}
#'     \item \code{"aqua_blue"}
#'   }
#' @param to Character. Name of the output fluorescence (e.g.
#'   \code{"plate_blue"}). Must be specified together with \code{from}. Should be one of:
#'   \itemize{
#'     \item \code{"aqua_red"}
#'     \item \code{"plate_blue"}
#'     \item \code{"plate_red"}
#'     \item \code{"aqua_blue"}
#'   }
#' @param intercept Numeric. Model intercept. Required when not using
#'   \code{from}/\code{to} lookup.
#' @param slope Numeric. Model slope. Required when not using
#'   \code{from}/\code{to} lookup.
#' @param estimate_error Logical. If \code{TRUE} (default), prediction intervals
#'   are computed and returned as \code{lower} and \code{upper} columns. If
#'   \code{FALSE}, only \code{predicted} is returned. Automatically set to
#'   \code{FALSE} with a warning if any of \code{se_intercept}, \code{se_slope},
#'   or \code{rse} are missing.
#' @param se_intercept Numeric. Standard error of the intercept. Required when
#'   not using \code{from}/\code{to} lookup.
#' @param se_slope Numeric. Standard error of the slope. Required when not
#'   using \code{from}/\code{to} lookup.
#' @param rse Numeric. Residual standard error of the model. Required when not
#'   using \code{from}/\code{to} lookup.
#' @param level Numeric. Confidence level for the prediction interval.
#'   Defaults to \code{0.95}.
#' @param transform Character. Scale on which the model was fitted. One of
#'   \code{"raw"} (default) or \code{"log10"}. When using \code{from}/\code{to}
#'   lookup this selects the appropriate pre-fitted model; when using manual
#'   parameters this controls back-transformation of predictions.
#' @param direction Character. One of \code{"forward"} (default, predict
#'   \code{to} from \code{from}) or \code{"backward"} (invert the model to
#'   predict \code{from} from \code{to}). Set automatically when using
#'   \code{from}/\code{to} lookup.
#'
#' @return A \code{\link[tibble]{tibble}} with one row per element of
#'   \code{ft_value} and columns:
#'   \describe{
#'     \item{ft_value}{The original input values.}
#'     \item{predicted}{The predicted equivalent fluorescence values.}
#'     \item{lower}{Lower bound of the prediction interval. Only present when
#'       \code{estimate_error = TRUE}.}
#'     \item{upper}{Upper bound of the prediction interval. Only present when
#'       \code{estimate_error = TRUE}.}
#'   }
#'
#' @details
#' The function can be used in two ways:
#'
#' \strong{Lookup mode}: supply \code{from}, \code{to}, and \code{strain} to
#' automatically retrieve model parameters from \code{\link{ft_convert_df}}.
#' The conversion direction is inferred from the variable order in the stored
#' model.
#'
#' \strong{Manual mode}: supply all five model parameters (\code{intercept},
#' \code{slope}, \code{se_intercept}, \code{se_slope}, \code{rse}) directly,
#' along with \code{transform} and \code{direction}.
#'
#' Prediction intervals are computed using the delta method. On the log10
#' scale, uncertainty is propagated in log-space and back-transformed, so
#' intervals are asymmetric on the original scale.
#'
#' @examples
#' # Lookup mode: convert aqua_blue to plate_blue for Ehux
#' ft_convert(c(1000, 5000, 10000),
#'            strain = "Ehux",
#'            from   = "aqua_blue",
#'            to     = "plate_blue",
#'            transform = "log10")
#'
#' # Manual mode
#' ft_convert(c(1000, 5000, 10000),
#'            strain       = "Ehux",
#'            intercept    = 0.2665854,
#'            slope        = 1.035234,
#'            se_intercept = 0.04975249,
#'            se_slope     = 0.01626522,
#'            rse          = 5.459139e-02,
#'            transform    = "log10",
#'            direction = "reverse")
#'
#' @seealso \code{\link{ft_convert_df}} for the underlying model parameters.
#' @export
ft_convert <- function(ft_value, strain, from = NULL, to = NULL,
                       intercept = NULL, slope = NULL, estimate_error = TRUE,
                       se_intercept = NULL, se_slope = NULL, rse = NULL,
                       level = 0.95, transform = "raw", direction = "forward") {
  transform <- match.arg(transform, c("raw", "log10"))

  # --- Input validation ---
  using_lookup <- !is.null(from) || !is.null(to)
  using_manual <- !is.null(intercept) && !is.null(slope)

  if (!using_lookup && !using_manual) {
    stop("Provide either 'from'/'to' to use internal function values, or all model parameters manually.")
  }

  if (using_lookup) {
    if (is.null(from) || is.null(to)) {
      stop("Both 'from' and 'to' must be specified together.")
    }
    strain <- match.arg(strain, levels(ft_convert_df$species)) # throws informative error automatically
    from <- match.arg(from, unique(c(ft_convert_df$x_var, ft_convert_df$y_var)))
    to <- match.arg(to, unique(c(ft_convert_df$x_var, ft_convert_df$y_var)))

    # pull model params from lookup table
    params <- ft_convert_df %>%
      filter(species == strain) %>%
      filter(scale == transform) %>%
      filter(y_var %in% c(from, to) & x_var %in% c(from, to))

    if (nrow(params) != 1) {
      stop("Expected 1 model but found ", nrow(params), " — have Henry check ft_convert_df for duplicate entries.")
    }

    intercept <- params$intercept
    slope <- params$slope
    se_intercept <- params$intercept.std.err
    se_slope <- params$slope.std.err
    rse <- params$rse

    if (from == params$x_var) {
      direction <- "forward"
    } else if (from == params$y_var) {
      direction <- "backward"
    } else {
      stop("'from' does not match either variable in the selected model.")
    }
  }

  # --- Error estimation check ---
  has_error_params <- !is.null(se_intercept) && !is.null(se_slope) && !is.null(rse)
  if (estimate_error && !has_error_params) {
    warning("estimate_error = TRUE but se_intercept, se_slope, and/or rse are missing. Switching to estimate_error = FALSE.")
    estimate_error <- FALSE
  }

  t_val <- qt((1 + level) / 2, df = Inf)

  if (transform == "log10") {
    if (direction == "forward") {
      log_x <- log10(ft_value)
      log_pred <- intercept + slope * log_x
      if (estimate_error) {
        se_fit <- sqrt(se_intercept^2 + (log_x * se_slope)^2)
        se_pred <- sqrt(rse^2 + se_fit^2)
      }
    } else {
      # Invert: log_x = (log_y - intercept) / slope
      log_y <- log10(ft_value)
      log_pred <- (log_y - intercept) / slope
      if (estimate_error) {
        se_fit <- sqrt((se_intercept / slope)^2 + (log_pred * se_slope / slope)^2)
        se_pred <- sqrt((rse / slope)^2 + se_fit^2)
      }
    }
    if (estimate_error) {
      tibble(
        ft_value  = ft_value,
        predicted = 10^log_pred,
        lower     = 10^(log_pred - t_val * se_pred),
        upper     = 10^(log_pred + t_val * se_pred)
      )
    } else {
      tibble(
        ft_value  = ft_value,
        predicted = 10^log_pred
      )
    }
  } else {
    if (direction == "forward") {
      predicted <- intercept + slope * ft_value
      if (estimate_error) {
        se_fit <- sqrt(se_intercept^2 + (ft_value * se_slope)^2)
        se_pred <- sqrt(rse^2 + se_fit^2)
      }
    } else {
      predicted <- (ft_value - intercept) / slope
      if (estimate_error) {
        se_fit <- sqrt((se_intercept / slope)^2 + (predicted * se_slope / slope)^2)
        se_pred <- sqrt((rse / slope)^2 + se_fit^2)
      }
    }
    if (estimate_error) {
      tibble(ft_value, predicted,
        lower = predicted - t_val * se_pred,
        upper = predicted + t_val * se_pred
      )
    } else {
      tibble(
        ft_value,
        predicted
      )
    }
  }
}

