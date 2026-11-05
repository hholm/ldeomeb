#' @importFrom magrittr %>%
NULL

# ---------------------------------------------------------------------------
# Internal helpers (not exported)
# ---------------------------------------------------------------------------

# lower-case, alphanumerics only: "pH_pairs " -> "phpairs", "Start Time:" -> "starttime"
.tp_norm <- function(x) tolower(gsub("[^[:alnum:]]", "", as.character(x)))

# trim whitespace, empty strings -> NA
.tp_cell <- function(x) {
  x <- trimws(as.character(x))
  x[!nzchar(x)] <- NA_character_
  x
}

# "A1", "a1", " A 01 ", "A-1" -> "A1"; anything else -> NA
.tp_well <- function(x) {
  x <- .tp_cell(x)
  m <- regmatches(x, regexec("^([A-Za-z]{1,2})\\s*[-_:.]?\\s*0*([0-9]{1,2})$", x))
  vapply(m, function(z) if (length(z) == 3) paste0(toupper(z[2]), as.integer(z[3])) else NA_character_, "")
}

# Accept a path (xlsx/xls/csv), a data.frame/tibble/matrix; return an all-character data.frame
# with no header (row 1 = first line of the sheet), NA for empty cells.
.tp_as_df <- function(x, sheet = 1) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("File not found: ", x, call. = FALSE)
    ext <- tolower(tools::file_ext(x))
    if (ext %in% c("xlsx", "xls", "xlsm")) {
      x <- suppressMessages(readxl::read_excel(x, sheet = sheet, col_names = FALSE,
                                               col_types = "text", .name_repair = "minimal"))
    } else if (ext %in% c("csv", "txt")) {
      nc <- max(utils::count.fields(x, sep = ",", quote = "\"", blank.lines.skip = FALSE), na.rm = TRUE)
      x <- utils::read.csv(x, header = FALSE, colClasses = "character", fill = TRUE,
                           col.names = paste0("V", seq_len(nc)), na.strings = character(0),
                           blank.lines.skip = FALSE, check.names = FALSE)
    } else stop("Unsupported file type '", ext, "'. Use .xlsx, .xls or .csv.", call. = FALSE)
  }
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Input has fewer than 2 columns - is this a TECAN export?", call. = FALSE)
  x[] <- lapply(x, .tp_cell)
  names(x) <- paste0("V", seq_len(ncol(x)))
  x
}

# plate column header: whole numbers >= 1, strictly increasing (1,2,3 ... or 1,8,9,10 for a partial plate)
.tp_is_colheader <- function(v) length(v) >= 1 && !anyNA(v) && all(v == round(v)) && all(v >= 1) && all(diff(v) > 0)

# Parse TECAN time stamps, tolerant of 12 h / 24 h clocks and m/d vs d/m ordering.
.tp_parse_time <- function(x) {
  if (inherits(x, "POSIXt")) return(x)
  lubridate::parse_date_time(trimws(as.character(x)), quiet = TRUE,
                             orders = c("mdY IMS p", "mdY HMS", "mdY IM p", "mdY HM", "dmY HMS", "dmY HM"))
}

# Parse every "Label: ..." block of a TECAN i-control sheet into one long data.frame.
# Anchors on the labels in column 1 (not on fixed row offsets), so extra/missing
# lines (e.g. no temperature line) or any plate size are fine.
.tp_parse_blocks <- function(dat) {
  col1 <- dat[[1]]
  n <- nrow(dat)
  label_rows <- which(grepl("^label\\s*:", col1, ignore.case = TRUE))
  if (!length(label_rows))
    stop("No 'Label:' blocks found in column A - is this a TECAN i-control export?", call. = FALSE)
  block_end <- c(label_rows[-1] - 1L, n)
  reserved <- c("rows", "cols", "well", "value", "Label", "Start Time", "Temperature")
  out <- list(); bad_all <- character(0)
  
  for (i in seq_along(label_rows)) {
    r0 <- label_rows[i]; r1 <- block_end[i]
    lab <- trimws(sub("^label\\s*:\\s*", "", col1[r0], ignore.case = TRUE))
    r_start <- (r0:r1)[which(.tp_norm(col1[r0:r1]) == "starttime")[1]]
    if (is.na(r_start)) { warning("Block '", lab, "' has no 'Start Time:' line; skipped.", call. = FALSE); next }
    
    # settings between the label and 'Start Time:' (key in col A, first filled cell to its right = value)
    srows <- if (r_start - r0 > 1) (r0 + 1):(r_start - 1) else integer(0)
    srows <- srows[!is.na(col1[srows])]
    vals <- vapply(srows, function(r) { v <- unlist(dat[r, -1], use.names = FALSE); v <- v[!is.na(v)]
    if (length(v)) v[1] else NA_character_ }, "")
    keys <- col1[srows]; keys[keys %in% reserved] <- paste0("setting_", keys[keys %in% reserved])
    settings <- stats::setNames(as.list(vals), keys); settings <- settings[!duplicated(names(settings))]
    start_txt <- { v <- unlist(dat[r_start, -1], use.names = FALSE); v[!is.na(v)][1] }
    
    # plate header row: "<>" or a 1..k run of integers
    after <- if (r_start < r1) (r_start + 1):r1 else integer(0)
    is_hdr <- vapply(after, function(r) {
      if (identical(col1[r], "<>")) return(TRUE)
      v <- suppressWarnings(as.numeric(unlist(dat[r, -1], use.names = FALSE)))
      k <- sum(cumsum(is.na(v)) == 0)
      k >= 2 && is.na(col1[r]) && .tp_is_colheader(v[seq_len(k)])
    }, NA)
    r_h <- after[which(is_hdr)[1]]
    if (is.na(r_h)) { warning("Block '", lab, "': no plate grid found; skipped.", call. = FALSE); next }
    hv <- suppressWarnings(as.numeric(unlist(dat[r_h, -1], use.names = FALSE)))
    k <- sum(cumsum(is.na(hv)) == 0); colnums <- as.integer(hv[seq_len(k)])
    
    # temperature line (between Start Time and the grid)
    temp <- NA_real_
    if (r_h - r_start > 1) {
      cells <- unlist(dat[(r_start + 1):(r_h - 1), ], use.names = FALSE)
      tc <- grep("^temperature\\s*:", cells, ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(tc)) temp <- suppressWarnings(as.numeric(sub("^[^0-9-]*(-?[0-9.]+).*$", "\\1", tc)))
    }
    
    # grid rows: consecutive rows whose col A is 1-2 letters
    rr <- r_h; data_rows <- integer(0)
    while (rr + 1 <= r1 && !is.na(col1[rr + 1]) && grepl("^[A-Za-z]{1,2}$", col1[rr + 1])) {
      rr <- rr + 1; data_rows <- c(data_rows, rr) }
    if (!length(data_rows)) { warning("Block '", lab, "': empty plate grid; skipped.", call. = FALSE); next }
    m <- as.matrix(dat[data_rows, 1 + seq_len(k), drop = FALSE])
    num <- suppressWarnings(as.numeric(m))
    bad <- !is.na(m) & is.na(num)                      # e.g. "OVER"
    keep <- !is.na(m) | bad
    rl <- toupper(col1[data_rows])
    df <- data.frame(rows = rep(rl, times = k), cols = rep(colnums, each = length(data_rows)),
                     stringsAsFactors = FALSE)[keep, , drop = FALSE]
    df$well <- paste0(df$rows, df$cols); df$value <- num[keep]
    if (any(bad)) bad_all <- c(bad_all, paste0(lab, ": ", paste(df$well[is.na(df$value)], collapse = ",")))
    df$Label <- lab
    for (nm in names(settings)) df[[nm]] <- settings[[nm]]
    df[["Start Time"]] <- start_txt; df[["Temperature"]] <- temp
    out[[length(out) + 1]] <- df
  }
  if (!length(out)) stop("Found 'Label:' lines but no readable plate grids.", call. = FALSE)
  if (length(bad_all)) warning("Non-numeric readings (e.g. 'OVER') set to NA - ", paste(bad_all, collapse = " | "), call. = FALSE)
  nms <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(b) { for (z in setdiff(nms, names(b))) b[[z]] <- NA_character_; b[nms] })
  res <- do.call(rbind, out); rownames(res) <- NULL
  res[["Start Time"]] <- .tp_parse_time(res[["Start Time"]])
  if (any(is.na(res[["Start Time"]]))) warning("Some 'Start Time' values could not be parsed.", call. = FALSE)
  mode <- tolower(res[["Mode"]]); mode[is.na(mode)] <- ""
  res$wavelength <- suppressWarnings(as.numeric(ifelse(grepl("absorb", mode),
                                                       if ("Measurement Wavelength" %in% nms) res[["Measurement Wavelength"]] else NA,
                                                       if ("Excitation Wavelength" %in% nms) res[["Excitation Wavelength"]] else NA)))
  res
}

# Legacy "Names" plate grid (header "Names", 1..k across, row letters down) -> named vector well -> name
.tp_names_grid <- function(dat) {
  col1 <- dat[[1]]
  for (h in which(.tp_norm(col1) %in% c("names", "samplenames"))) {
    v <- suppressWarnings(as.numeric(unlist(dat[h, -1], use.names = FALSE)))
    k <- sum(cumsum(is.na(v)) == 0)
    if (k < 1 || !.tp_is_colheader(v[seq_len(k)])) next
    cn <- as.integer(v[seq_len(k)])
    rr <- h; rows <- integer(0)
    while (rr + 1 <= nrow(dat) && !is.na(col1[rr + 1]) && grepl("^[A-Za-z]{1,2}$", col1[rr + 1])) { rr <- rr + 1; rows <- c(rows, rr) }
    if (!length(rows)) next
    m <- as.matrix(dat[rows, 1 + seq_len(k), drop = FALSE])
    w <- paste0(rep(toupper(col1[rows]), times = k), rep(cn, each = length(rows)))
    nm <- trimws(gsub("\\s+", " ", as.vector(m)))
    ok <- !is.na(nm) & nzchar(nm)
    return(stats::setNames(nm[ok], w[ok]))
  }
  NULL
}

# Find a sample/pairs table anywhere in a sheet: a header row containing a 'blank(s)' and a 'dye(s)' cell.
# Returns data.frame(sample, blank, dye, salinity) of raw (uncleaned) strings, or NULL.
.tp_find_pairs <- function(dat) {
  p <- ncol(dat)
  blank_syn <- c("blank", "blanks", "blankwell", "blankwells")
  dye_syn   <- c("dye", "dyes", "dyewell", "dyewells")
  for (r in seq_len(nrow(dat))) {
    h <- .tp_norm(unlist(dat[r, seq_len(min(p, 10))], use.names = FALSE))
    bc <- which(h %in% blank_syn)[1]; dc <- which(h %in% dye_syn)[1]
    if (is.na(bc) || is.na(dc)) next
    sc <- which(h %in% c("sample", "samples", "samplename", "samplenames", "name", "names", "id", "sampleid"))[1]
    if (is.na(sc) && !(1 %in% c(bc, dc))) sc <- 1
    salc <- which(h %in% c("salinity", "sal", "salinitypsu", "psu"))[1]
    rows <- if (r < nrow(dat)) (r + 1):nrow(dat) else integer(0)
    stop_at <- which(.tp_norm(dat[[1]][rows]) %in% c("note", "notes", "comment", "comments", "remarks"))[1]
    if (!is.na(stop_at)) rows <- rows[seq_len(stop_at - 1)]
    sub <- dat[rows, , drop = FALSE]
    keep <- !is.na(sub[[bc]]) | !is.na(sub[[dc]])
    sub <- sub[keep, , drop = FALSE]
    if (!nrow(sub)) return(NULL)
    return(data.frame(sample = if (is.na(sc)) NA_character_ else sub[[sc]][seq_len(nrow(sub))],
                      blank = sub[[bc]], dye = sub[[dc]],
                      salinity = if (is.na(salc)) NA_character_ else sub[[salc]],
                      stringsAsFactors = FALSE)[seq_len(nrow(sub)), , drop = FALSE])
  }
  NULL
}

# Turn a user-supplied data.frame with column headers into the same raw form
.tp_pairs_from_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  raw <- rbind(as.character(names(df)), as.matrix(df))
  raw <- as.data.frame(raw, stringsAsFactors = FALSE); raw[] <- lapply(raw, .tp_cell); names(raw) <- paste0("V", seq_len(ncol(raw)))
  out <- .tp_find_pairs(raw)
  if (is.null(out)) stop("Sample table needs columns for the blank well and the dye well ",
                         "(e.g. 'sample', 'blank', 'dye').", call. = FALSE)
  out
}

# Raw pairs table -> canonical sample table (one row per blank/dye well pair).
# Accepts explicit wells ("A1","B1") and legacy whole-row pairs ("A","B", expanded across columns).
.tp_resolve_pairs <- function(pairs, wells, name_lookup = NULL, on_unmeasured = c("warn", "stop")) {
  on_unmeasured <- match.arg(on_unmeasured)
  res <- list(); mism <- character(0); skipped <- character(0)
  wr <- sub("[0-9]+$", "", wells); wc <- as.integer(sub("^[A-Z]+", "", wells))
  for (i in seq_len(nrow(pairs))) {
    b <- .tp_cell(pairs$blank[i]); d <- .tp_cell(pairs$dye[i])
    if (is.na(b) && is.na(d)) next
    if (is.na(b) || is.na(d)) {
      one <- if (is.na(b)) d else b
      if (!is.na(.tp_well(one)) || grepl("^[A-Za-z]{1,2}$", one))
        stop("Sample table row ", i, " ('", pairs$sample[i], "') has a ", if (is.na(b)) "dye" else "blank",
             " well but no ", if (is.na(b)) "blank" else "dye", " well.", call. = FALSE)
      next    # free text (a note), not a sample
    }
    bw <- .tp_well(b); dw <- .tp_well(d)
    sal <- suppressWarnings(as.numeric(pairs$salinity[i]))
    if (!is.na(bw) && !is.na(dw)) {
      nm <- .tp_cell(pairs$sample[i])
      if (is.na(nm) && !is.null(name_lookup)) nm <- unname(c(name_lookup[dw], name_lookup[bw]))[!is.na(c(name_lookup[dw], name_lookup[bw]))][1]
      if (is.na(nm)) nm <- paste0(bw, "_", dw)
      res[[length(res) + 1]] <- data.frame(sample = nm, blank = bw, dye = dw, salinity = sal, stringsAsFactors = FALSE)
    } else if (grepl("^[A-Za-z]{1,2}$", b) && grepl("^[A-Za-z]{1,2}$", d)) {
      rb <- toupper(b); rd <- toupper(d)
      cc <- sort(intersect(wc[wr == rb], wc[wr == rd]))
      if (!length(cc)) next
      wb <- paste0(rb, cc); wd <- paste0(rd, cc)
      nb <- if (is.null(name_lookup)) rep(NA_character_, length(cc)) else unname(name_lookup[wb])
      nd <- if (is.null(name_lookup)) rep(NA_character_, length(cc)) else unname(name_lookup[wd])
      if (!is.null(name_lookup) && any(!is.na(c(nb, nd)))) {
        unnamed <- is.na(nb) & is.na(nd)
        skipped <- c(skipped, paste0(wb, "/", wd)[unnamed]); keepc <- !unnamed
      } else keepc <- rep(TRUE, length(cc))
      nm <- ifelse(is.na(nd), nb, nd); nm[is.na(nm)] <- paste0(wb, "_", wd)[is.na(nm)]
      diff <- !is.na(nb) & !is.na(nd) & tolower(nb) != tolower(nd)
      if (any(diff & keepc)) mism <- c(mism, paste0(wb, "='", nb, "' vs ", wd, "='", nd, "'")[diff & keepc])
      res[[length(res) + 1]] <- data.frame(sample = nm, blank = wb, dye = wd, salinity = sal,
                                           stringsAsFactors = FALSE)[keepc, , drop = FALSE]
    } else stop("Sample table row ", i, ": could not read blank='", b, "' dye='", d,
                "'. Use well IDs like A1 / B1 (or whole-row letters like A / B).", call. = FALSE)
  }
  if (!length(res)) stop("The sample table has no usable rows.", call. = FALSE)
  s <- do.call(rbind, res); rownames(s) <- NULL
  auto <- s$sample == paste0(s$blank, "_", s$dye)
  if (any(auto)) warning(sum(auto), " sample(s) have no name and are labelled by well: ",
                         paste(utils::head(s$sample[auto], 4), collapse = ", "), if (sum(auto) > 4) " ..." else "",
                         ". Check the Names grid / sample column.", call. = FALSE)
  if (length(mism)) warning("Blank and dye wells have different names; using the dye-well name: ",
                            paste(utils::head(mism, 5), collapse = "; "), if (length(mism) > 5) " ...", call. = FALSE)
  if (length(skipped)) warning(length(skipped), " measured well pair(s) have no name and were skipped: ",
                               paste(utils::head(skipped, 6), collapse = ", "), if (length(skipped) > 6) " ...", call. = FALSE)
  allw <- c(s$blank, s$dye)
  if (anyDuplicated(allw)) stop("Well(s) used more than once in the sample table: ",
                                paste(unique(allw[duplicated(allw)]), collapse = ", "), call. = FALSE)
  miss <- setdiff(allw, wells)
  if (length(miss)) {
    msg <- paste0("Well(s) in the sample table were not read by the plate reader: ", paste(miss, collapse = ", "))
    if (on_unmeasured == "stop") stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  s$position <- paste(s$blank, s$dye, sep = "-")
  s
}

# ---------------------------------------------------------------------------
# Exported functions
# ---------------------------------------------------------------------------

#' Import and format plate reader data from TECAN i-control.
#'
#' Reads a TECAN i-control export and returns every plate read in long format, plus the
#' blank/dye sample table. The raw sheet can be left untouched: sample information may be
#' (in order of precedence) passed via `samples`, placed on a second sheet called
#' "Samples", or typed at the bottom of the data sheet.
#'
#' @details
#' The sample table has one row per sample and the columns `sample`, `blank`, `dye` and
#' (optionally) `salinity`; `blank` and `dye` are well IDs such as `A1` / `B1`. Any well can be
#' a blank or a dye. Column names are matched loosely (case, spaces and plural are ignored) and
#' well IDs may be written `a1`, `A01`, `A 1`. The older layout (a "Names" plate grid plus a
#' "pH_pairs" table of whole-row letters) is still understood.
#'
#' @param dat A path to an .xlsx/.xls/.csv TECAN export, or an already imported
#'   data.frame/tibble (read with `col_names = FALSE`).
#' @param sheet Sheet of the workbook holding the plate data (name or number). Default 1.
#' @param samples `NULL` (auto-detect, default), a data.frame or file path with the sample table,
#'   or `FALSE` to skip sample handling.
#' @return A list with `plates` (long data.frame: `rows`, `cols`, `well`, `value`, `Label`,
#'   reader settings, `Start Time`, `Temperature`, `wavelength`, `names`, `type`),
#'   `samples` (the sample table) and `pairs` (same as `samples`, kept for backwards compatibility).
#' @export
tidyplate <- function(dat, sheet = 1, samples = NULL) {
  path <- if (is.character(dat) && length(dat) == 1L) dat else NULL
  d <- .tp_as_df(dat, sheet)
  plates <- .tp_parse_blocks(d)
  wells <- unique(plates$well)
  
  # --- locate sample information ---
  tbl <- NULL; grid <- .tp_names_grid(d)
  if (isFALSE(samples)) {
    tbl <- NULL
  } else if (!is.null(samples)) {
    if (is.character(samples) && length(samples) == 1L) {
      s2 <- .tp_as_df(samples, sheet = if (grepl("\\.xlsx?$", samples, TRUE)) 1 else NULL)
      tbl <- .tp_find_pairs(s2); grid <- .tp_names_grid(s2) %||% grid
    } else tbl <- .tp_pairs_from_df(samples)
  } else {
    if (!is.null(path) && grepl("\\.xls[xm]?$", path, ignore.case = TRUE)) {
      sh <- readxl::excel_sheets(path)
      cand <- which(.tp_norm(sh) %in% c("samples", "samplemap", "platemap", "samplesheet", "samplelist", "phpairs", "pairs", "names"))
      cand <- setdiff(cand, if (is.numeric(sheet)) sheet else match(sheet, sh))
      if (length(cand)) {
        s2 <- .tp_as_df(path, sheet = cand[1]); tbl <- .tp_find_pairs(s2); grid <- .tp_names_grid(s2) %||% grid
      }
    }
    if (is.null(tbl)) tbl <- .tp_find_pairs(d)
  }
  smp <- NULL
  if (!is.null(tbl)) smp <- .tp_resolve_pairs(tbl, wells, name_lookup = grid, on_unmeasured = "warn")
  else if (!isFALSE(samples)) warning("No sample table (blank/dye wells) found. Add a 'Samples' sheet or table, ",
                                      "or pass `samples =` to calc_plate().", call. = FALSE)
  
  plates$names <- NA_character_; plates$type <- NA_character_
  if (!is.null(smp)) {
    ib <- match(plates$well, smp$blank); id <- match(plates$well, smp$dye)
    plates$type[!is.na(ib)] <- "blank"; plates$names[!is.na(ib)] <- smp$sample[ib[!is.na(ib)]]
    plates$type[!is.na(id)] <- "dye";   plates$names[!is.na(id)] <- smp$sample[id[!is.na(id)]]
  }
  attr(plates, "name_lookup") <- grid
  list(plates = plates, samples = smp, pairs = smp)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Calculate pH values from m-cresol dye absorbance.
#'
#' @param A730_blank,A578_blank,A434_blank Absorbance of the sample without dye at 730, 578, 434 nm.
#' @param A730_dye,A578_dye,A434_dye Absorbance of the sample with m-cresol dye at 730, 578, 434 nm.
#' @param vol.dye.uL Volume of dye added, in microlitres (default 3).
#' @param salinity Sample salinity in PSU; a single value or one per sample (default 35).
#' @param verbose Return a data.frame with intermediate values (`TRUE`) or only the pH vector (`FALSE`).
#' @param calibration Dye-addition correction, per microlitre of dye (default -0.00660453).
#' @return A data.frame (`verbose = TRUE`) or numeric vector of pH values.
#' @export
calc_pH_spec <- function(A730_blank, A578_blank, A434_blank, A730_dye, A578_dye,
                         A434_dye, vol.dye.uL = 3, salinity = 35, verbose = TRUE,
                         calibration = -0.00660453) {
  # following hennon formula
  A1_A2 <- (A578_dye - A578_blank - (A730_dye - A730_blank)) / (A434_dye - A434_blank - (A730_dye - A730_blank))
  pK2 <- (1245.9 / 298) + 3.8275 + (0.00211 * (35 - salinity))
  A1_A2_cor <- A1_A2 + (-1 * calibration * vol.dye.uL)
  pH <- pK2 + log10((A1_A2_cor - 0.00691) / (2.222 - A1_A2_cor * 0.1331))
  if (verbose) data.frame(A1_A2, pK2, A1_A2_cor, pH) else pH
}

#' Calculate pH values from m-cresol dye absorbance for an entire plate at once.
#'
#' Pairs each sample's blank and dye wells (from the sample table), picks the 434/578/730 nm
#' absorbance reads by wavelength, and returns one row per sample. Fluorescence reads
#' (`Ft_<excitation>`) are taken from the blank well.
#'
#' @param tidyplate The list returned by [tidyplate()] (preferred) or its `$plates` data.frame.
#' @param verbose `TRUE` for all intermediate columns, `FALSE` for name, position and pH only.
#' @param pairs,samples Optional sample table overriding the one found by [tidyplate()]
#'   (`pairs` is the old argument name). There is no default layout any more: the wells come
#'   from the sample table.
#' @param salinity Default salinity in PSU, used where the sample table has no `salinity` value.
#' @param vol.dye.uL Volume of dye added in microlitres (default 3).
#' @param calibration Dye-addition correction per microlitre (default -0.00660453).
#' @param wavelengths Wavelengths (nm) of the three absorbance reads, in the order 434, 578, 730.
#' @param wavelength_tol Tolerance in nm when matching reads to `wavelengths` (default 5).
#' @param duplicates What to do if a well was read twice at one wavelength: `"error"` (default) or `"mean"`.
#' @return A data.frame with one row per sample. Check the `flag` column for suspicious samples.
#' @export
calc_plate <- function(tidyplate, verbose = TRUE, pairs = NULL, salinity = 35, vol.dye.uL = 3,
                       calibration = -0.00660453, samples = NULL, wavelengths = c(434, 578, 730),
                       wavelength_tol = 5, duplicates = c("error", "mean")) {
  duplicates <- match.arg(duplicates)
  plates <- tidyplate; sm <- NULL
  if (is.list(tidyplate) && !is.data.frame(tidyplate)) { plates <- tidyplate$plates; sm <- tidyplate$samples }
  if (!is.data.frame(plates)) stop("`tidyplate` must be the output of tidyplate().", call. = FALSE)
  grid <- attr(plates, "name_lookup")
  plates <- as.data.frame(plates, stringsAsFactors = FALSE)
  if (!"well" %in% names(plates)) plates$well <- paste0(plates$rows, as.integer(plates$cols))
  if (is.null(sm)) sm <- attr(tidyplate, "samples")
  user <- if (!is.null(samples)) samples else pairs
  wells <- unique(plates$well)
  if (!is.null(user)) {
    lk <- grid
    if (is.null(lk) && "names" %in% names(plates)) { k <- !is.na(plates$names); lk <- stats::setNames(plates$names[k], plates$well[k]); lk <- lk[!duplicated(names(lk))] }
    sm <- .tp_resolve_pairs(.tp_pairs_from_df(user), wells, lk, "stop")
  } else if (is.null(sm)) {
    stop("No sample table (blank/dye wells). Add a 'Samples' sheet/table to the workbook or pass `samples =`.", call. = FALSE)
  } else sm <- .tp_resolve_pairs(.tp_pairs_from_df(sm), wells, NULL, "stop")
  n <- nrow(sm)
  
  mode <- tolower(plates$Mode); mode[is.na(mode)] <- ""
  if (!"wavelength" %in% names(plates))
    plates$wavelength <- suppressWarnings(as.numeric(ifelse(grepl("absorb", mode), plates[["Measurement Wavelength"]], plates[["Excitation Wavelength"]])))
  if (!"Label" %in% names(plates)) plates$Label <- NA_character_
  st <- .tp_parse_time(plates[["Start Time"]])
  
  # one value per well for a subset of reads
  well_values <- function(sel, what) {
    d <- plates[sel, ]
    if (anyDuplicated(d$well) && duplicates == "error")
      stop(what, ": well(s) ", paste(utils::head(unique(d$well[duplicated(d$well)]), 6), collapse = ", "),
           " were read more than once (labels: ", paste(unique(d$Label), collapse = ", "),
           "). Remove the extra read or use duplicates = \"mean\".", call. = FALSE)
    tapply(d$value, d$well, mean)
  }
  absn <- grepl("absorb", mode)
  get_abs <- function(wl) {
    sel <- absn & !is.na(plates$wavelength) & abs(plates$wavelength - wl) <= wavelength_tol
    if (!any(sel)) stop("No absorbance read within ", wavelength_tol, " nm of ", wl, " nm. Absorbance reads found: ",
                        paste(sort(unique(plates$wavelength[absn])), collapse = ", "), " nm.", call. = FALSE)
    v <- well_values(sel, paste0(wl, " nm"))
    miss <- setdiff(c(sm$blank, sm$dye), names(v))
    if (length(miss)) stop("Well(s) ", paste(miss, collapse = ", "), " have no ", wl, " nm reading.", call. = FALSE)
    v
  }
  v434 <- get_abs(wavelengths[1]); v578 <- get_abs(wavelengths[2]); v730 <- get_abs(wavelengths[3])
  g <- function(v, w) unname(v[w])
  
  sal <- sm$salinity; sal[is.na(sal)] <- salinity
  ph <- suppressWarnings(calc_pH_spec(
    A730_blank = g(v730, sm$blank), A578_blank = g(v578, sm$blank), A434_blank = g(v434, sm$blank),
    A730_dye = g(v730, sm$dye), A578_dye = g(v578, sm$dye), A434_dye = g(v434, sm$dye),
    vol.dye.uL = vol.dye.uL, salinity = sal, verbose = TRUE, calibration = calibration))
  
  # per-sample mean start time of the reads covering its wells
  blk <- unique(data.frame(well = plates$well, Label = plates$Label, t = st, stringsAsFactors = FALSE))
  avg_t <- do.call(c, lapply(seq_len(n), function(i) {
    b <- unique(blk[blk$well %in% c(sm$blank[i], sm$dye[i]), c("Label", "t")]); mean(b$t) }))
  
  out <- data.frame(names = sm$sample, position = sm$position, blank_well = sm$blank, dye_well = sm$dye,
                    `Average Start Time` = avg_t, salinity = sal, check.names = FALSE, stringsAsFactors = FALSE)
  
  # fluorescence: taken from the blank well
  fl <- grepl("fluoresc", mode) & !is.na(plates$wavelength) & plates$well %in% sm$blank
  for (wl in sort(unique(plates$wavelength[fl]))) {
    sel <- fl & plates$wavelength == wl
    pre <- if (any(grepl("bottom", mode[sel]))) "Fb_" else "Ft_"
    out[[paste0(pre, wl)]] <- g(well_values(sel, paste0("fluorescence ", wl, " nm")), sm$blank)
  }
  for (nm in c("730", "578", "434")) {
    v <- get(paste0("v", nm))
    out[[paste0("blank_", nm)]] <- g(v, sm$blank); out[[paste0("dye_", nm)]] <- g(v, sm$dye)
  }
  out <- cbind(out, ph)
  
  # sanity flags
  chk <- list(
    "missing/non-numeric absorbance" = !stats::complete.cases(out[, c("blank_730", "dye_730", "blank_578", "dye_578", "blank_434", "dye_434")]),
    "dye well not above blank at 434 nm (dye missing or wells swapped?)" = (out$dye_434 - out$blank_434) <= 0,
    "dye well not above blank at 578 nm (dye missing or wells swapped?)" = (out$dye_578 - out$blank_578) <= 0,
    "pH not computable" = is.na(out$pH))
  out$flag <- vapply(seq_len(n), function(i) paste(names(chk)[vapply(chk, function(z) isTRUE(z[i]), NA)], collapse = "; "), "")
  if (any(nzchar(out$flag)))
    warning(sum(nzchar(out$flag)), " of ", n, " samples flagged (see `flag` column): ",
            paste(utils::head(out$names[nzchar(out$flag)], 6), collapse = ", "), call. = FALSE)
  
  if (!verbose) return(out[, c("names", "position", "pH")])
  rownames(out) <- NULL
  out
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
      dplyr::filter(species == strain) %>%
      dplyr::filter(scale == transform) %>%
      dplyr::filter(y_var %in% c(from, to) & x_var %in% c(from, to))
    
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
      tibble::tibble(
        ft_value  = ft_value,
        predicted = 10^log_pred,
        lower     = 10^(log_pred - t_val * se_pred),
        upper     = 10^(log_pred + t_val * se_pred)
      )
    } else {
      tibble::tibble(
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
      tibble::tibble(ft_value, predicted,
                     lower = predicted - t_val * se_pred,
                     upper = predicted + t_val * se_pred
      )
    } else {
      tibble::tibble(
        ft_value,
        predicted
      )
    }
  }
}

#' Summarize Fluorometer CSV Files by Mean and Standard Deviation
#'
#' Reads multiple CSV files from a user-selected directory, extracts the first
#' 10 observations from each file, and computes summary statistics for
#' fluorometer datapoints. The function outputs a single CSV file containing
#' per-file means, standard deviations, and timestamps.
#'
#' @details
#' This is set up to expect data output by the Turner Fluorometer software. This function assumes each CSV file:
#' \itemize{
#'   \item Contains at least 10 rows of data after skipping the first 3 lines.
#'   \item Has columns named \code{Date}, \code{Time}, and \code{Datapoint.Value}.
#'   \item Stores \code{Date} in \code{"\%m/\%d/\%Y"} format.
#'   \item Stores \code{Time} in \code{"\%H:\%M:\%S"} format.
#' }
#'
#' Only the first 10 rows of each file are used to compute:
#' \itemize{
#'   \item Mean date (converted to POSIXct and averaged)
#'   \item Mean time (converted to POSIXct and averaged)
#'   \item Mean fluorometer value
#'   \item Standard deviation of fluorometer value
#'   \item Number of observations used (\code{n = 10}, unless truncated)
#' }
#'
#' Files generated by previous runs of this function (matching
#' \code{"file_means_"}) are excluded from processing.
#'
#' @return
#' Invisibly returns \code{NULL}. A CSV file is written to the selected directory
#' with the naming convention:
#' \preformatted{
#' file_means_YYYYMMDDHHMMSS.csv
#' }
#'
#' The output file contains the following columns:
#' \describe{
#'   \item{file}{Name of the source CSV file}
#'   \item{mean.date}{Mean date of first 10 entries}
#'   \item{mean.time}{Mean time of first 10 entries}
#'   \item{mean.ft}{Mean fluorometer value}
#'   \item{sd.ft}{Standard deviation of fluorometer value}
#'   \item{n.ft}{Number of observations used}
#' }
#'
#' @note
#' This function uses \code{choose.dir()}, which is only available on Windows.
#' It is not portable to macOS or Linux without modification.
#'
#' @seealso
#' \code{\link[utils]{read.csv}}, \code{\link[utils]{write.csv}}
#'
#' @examples
#' \dontrun{
#' average_turner()
#' }
#'
#' @export
average_turner <- function(wd = NA) {
  
  # Option 2: Let the user choose the directory if files are not in the current working directory
  if(is.na(wd)){
    wd <- choose.dir(default = getwd(), caption = "Select folder with CSV files.")
  }
  # List all CSV files in the directory (skip files that don't match pattern)
  csv_files <- list.files(wd, pattern = "\\.csv$", full.names = TRUE)
  
  # Drop any csv output by this script
  csv_files <- csv_files[!grepl(csv_files, pattern = "/file_means_")]
  
  cat("Summarizing ", length(csv_files), " files.\n")
  
  # Initialize an empty list to store the data from each CSV
  csv <- list()
  
  # Read each CSV file, skipping the first 3 rows
  for (file in csv_files) {
    csv[[file]] <- read.csv(file, skip = 3, colClasses = c("character", "character", "numeric", "NULL"))
  }
  
  # Initialize an empty dataframe to store the summary output
  out <- data.frame()
  
  # Loop through each CSV file and calculate summary statistics
  for (file in seq_along(csv)) {
    out <- rbind(
      out,
      data.frame(
        file = basename(csv_files[file]), # Use basename to avoid full path in the output
        mean.date = format(mean(as.POSIXct(csv[[file]][1:10, "Date"], format = "%m/%d/%Y")), "%m/%d/%Y"),
        mean.time = format(mean(as.POSIXct(csv[[file]][1:10, "Time"], format = "%H:%M:%S")), "%H:%M:%S"),
        mean.ft = mean(csv[[file]][1:10, "Datapoint.Value"], na.rm = TRUE),
        sd.ft = sd(csv[[file]][1:10, "Datapoint.Value"], na.rm = TRUE),
        n.ft = nrow(csv[[file]][1:10, ])
      )
    )
  }
  
  # Define the output file name with timestamp
  output_file <- paste0(wd, "/file_means_", format(Sys.time(), "%Y%m%d%H%M%S"), ".csv")
  
  # Write the output to a CSV file
  write.csv(out, file = output_file, row.names = FALSE)
  
  # Notify the user that the script has completed
  cat("Summary of fluormeter readings written to:", output_file, "\n")
  out
}