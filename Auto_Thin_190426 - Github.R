#### 00 USER SETTINGS ####
# ---- Paths ----
in_dir <- file.path("input", "Stand 91")
tls_in  <- file.path(in_dir, "TLS_Stand_91_pur.csv")                     # TLS 3DFin input CSV
mls_in  <- file.path(in_dir, "MLS_Friday_0045_sc_Katze_Stand91.csv") # MLS 3DFin input CSV
hand_in <- file.path(in_dir, "22012026_Stand91_Bäume_mit_XY_Hand.csv")                          # Hand/human input CSV or Excel (.xlsx/.xls)

hand_sheet <- "Sheet1"                                         # only used if hand_in is Excel
out_dir <- file.path("output", "Stand 91")                                             # output directory (no setwd)

# ---- Parsing ----
name_repair <- "minimal"                                       # readr name repair
guess_max   <- 200000                                          # readr guess_max
verbose     <- FALSE                                           # diagnostic messages
do_plots    <- TRUE                                            # write QC plots (local-rank maps)
do_write_gpkg_intermediate <- TRUE                             # write intermediate per-sensor GPKG

# ---- DBH imputation (DBH<=0 -> NA) ----
set.seed(42)
rf_predictors    <- c("TH", "X", "Y")                           # TH always; optionally X/Y
rf_num_trees     <- 500
rf_mtry          <- 2
rf_min_node_size <- 5
rf_importance    <- "impurity"
dbh_min_valid_m  <- 0.01                                        # values <= 1 cm treated as missing in meter-based DBH fields
mls_use_tls_transfer_model <- TRUE                              # use TLS TH->DBH model to impute MLS DBH gaps
mls_transfer_predictors <- c("TH")                              # transfer predictors (recommended: TH only)

# ---- TreeCompR thinning parameters ----
ci_methods     <- c("CI_Hegyi", "CI_Braathe", "CI_RK1", "CI_RK2", "CI_RK3", "CI_RK4")
compete_radius <- 7
target_radius  <- 7
target_source_treecomp <- "exclude_edge"                        # docs: use "all_trees" only for true full-census boundaries
local_radius   <- 4
min_neigh      <- 3
local_q        <- 0.70                                          # worst 30% => remove

# Which thinning decision is used downstream (optional):
# - set to a single thin_* column name (e.g. "thin_CI_Hegyi") to filter out "remove"
# - set to NULL or "" to carry all trees forward
THIN_RULE_COL <- NULL

# ---- Local->UTM32 conversion (add centerpoint) ----
E0 <- 612780.911
N0 <- 6618078.542
x_col <- "x_local"
y_col <- "y_local"
crs_out <- 32632

# ---- Hand conversion (only if needed) ----
# If hand is polar (retning/avstand) or Excel, these are used.
X_center <- 1226217
Y_center <- 8329487
hand_angle_col <- "retning"
hand_dist_col  <- "avstand"
hand_dist_divisor <- 100                                        # 100: cm->m; 1: already meters
crs_hand_raw <- 25833

# ---- Alignment (similarity ICP; geometry preserved in meters) ----
icp_max_iter   <- 60
icp_dmax_start <- 35
icp_dmax_end   <- 2.5
icp_min_pairs  <- 12
icp_trim_frac  <- 0.80
icp_tol        <- 1e-4

icp_fit_voxel_cell_m <- 1.0
icp_fit_max_points   <- 1500L
icp_scale_min        <- 0.85
icp_scale_max        <- 1.15

match_dmax_m   <- 5.0

# ---- Between-media diagnostics population ----
between_hand_frit_only <- FALSE
between_min_n <- 10L

# ---- Output naming ----
stand_id     <- "91"
thin_tag     <- "30p"
tls_tag      <- "TLS"
mls_tag      <- "MLS"
out_gpkg_base <- paste0("trees_aligned_density_Stand", stand_id, "_", thin_tag)

# ---- Visualization outputs (map + layered bar plot) ----
do_visualizations <- TRUE
vis_preview_plot <- FALSE
vis_interactive_preview <- FALSE
vis_crs_epsg <- crs_out

vis_overlap_csv <- file.path(out_dir, paste0("Stand_", stand_id, "_overlap_points.csv"))
vis_all_mls_csv <- file.path(out_dir, paste0("Stand_", stand_id, "_all_mls_registered.csv"))
vis_mls_traj_csv <- file.path(out_dir, paste0("Stand_", stand_id, "_mls_trajectory.csv"))
vis_map_png <- file.path(out_dir, paste0("map_thinned_overlap_Stand", stand_id, ".png"))
vis_bar_png <- file.path(out_dir, paste0("bar_thinned_overlap_Stand", stand_id, ".png"))

vis_dpi_out <- 300
vis_map_width_in <- 12
vis_map_height_in <- 8
vis_bar_width_in <- 13
vis_bar_height_in <- 6.5

vis_thin_rule_col <- NULL
vis_thin_remove_values <- c("remove", "removed", "cut", "thin", "x", "1", "true", "t", "yes", "y")

vis_map_title <- paste0("Comparision map of thinned trees in stand ", stand_id, "\nby recording media")
vis_author_text <- "Author: Janek Wuigk"
vis_date_text <- paste0("Date: ", format(Sys.Date(), "%d.%m.%Y"))
vis_crs_text <- "Coordinate Reference System:\nCRS: WGS 84 / UTM 32N\n(EPSG:32632)"

# ---- Proximity-merged thinning layer ----
do_write_proximity_thin_layer <- TRUE
proximity_max_match_m <- 1.0
proximity_primary_thin_col <- "thin_CI_Braathe"  # primary decision rule for TLS/MLS
proximity_k_remove <- 3L
proximity_thin_cols_5 <- c("thin_CI_Hegyi", "thin_CI_Braathe", "thin_CI_RK1", "thin_CI_RK2", "thin_CI_RK3")
proximity_hand_thin_candidates <- c("Frit", "friT", "frit", "thin_hand", "HandThin")
proximity_gpkg <- file.path(out_dir, paste0("Stand_", stand_id, "_merged_thin_by_proximity.gpkg"))
proximity_layer <- "THIN_merged"
proximity_csv <- file.path(out_dir, paste0("Stand_", stand_id, "_merged_thin_by_proximity.csv"))

#### 01 LIBRARIES ####

suppressPackageStartupMessages({
  library(TreeCompR)
  library(readr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(rlang)
  library(sf)
  library(ggplot2)
  library(RANN)
  library(ranger)
  library(readxl)
})

options(stringsAsFactors = FALSE)

#### 02 HELPERS ####

dir_create <- function(p) {
  if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
  invisible(p)
}

stamp <- function() format(Sys.time(), "%Y%m%d_%H%M%S")

fmt_txt_value <- function(x) {
  if (length(x) == 0) return("")
  if (is.null(x)) return("NULL")
  if (length(x) > 1) return(paste(vapply(as.list(x), fmt_txt_value, character(1)), collapse = ", "))
  if (is.logical(x)) return(ifelse(is.na(x), "NA", ifelse(x, "TRUE", "FALSE")))
  if (is.numeric(x)) return(ifelse(is.na(x), "NA", as.character(x)))
  if (inherits(x, "factor")) return(as.character(x))
  x_chr <- as.character(x)
  ifelse(is.na(x_chr), "NA", x_chr)
}

write_run_manifest_txt <- function(path, run_id, setting_names, output_named_list, extra_named_list = NULL) {
  settings_vals <- mget(setting_names, envir = environment(), inherits = TRUE, ifnotfound = list(NA))
  lines <- c(
    paste0("run_id=", run_id),
    paste0("created_at=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "[settings]"
  )
  for (nm in setting_names) {
    lines <- c(lines, paste0(nm, "=", fmt_txt_value(settings_vals[[nm]])))
  }
  lines <- c(lines, "", "[outputs]")
  for (nm in names(output_named_list)) {
    lines <- c(lines, paste0(nm, "=", fmt_txt_value(output_named_list[[nm]])))
  }
  if (!is.null(extra_named_list) && length(extra_named_list) > 0) {
    lines <- c(lines, "", "[extra]")
    for (nm in names(extra_named_list)) {
      lines <- c(lines, paste0(nm, "=", fmt_txt_value(extra_named_list[[nm]])))
    }
  }
  dir_create(dirname(path))
  writeLines(lines, con = path, useBytes = TRUE)
  invisible(path)
}

to_num <- function(v) {
  if (is.numeric(v)) return(as.numeric(v))
  x <- trimws(as.character(v))
  x[x %in% c("", "NA", "NaN", "nan", "null", "NULL")] <- NA_character_
  x <- gsub("\u00a0", "", x, fixed = TRUE)
  x <- gsub(" ", "", x, fixed = TRUE)

  parse_one <- function(s) {
    if (is.na(s) || !nzchar(s)) return(NA_real_)

    sign <- 1
    if (substr(s, 1, 1) == "-") {
      sign <- -1
      s <- substring(s, 2)
    } else if (substr(s, 1, 1) == "+") {
      s <- substring(s, 2)
    }
    if (!nzchar(s)) return(NA_real_)

    has_comma <- grepl(",", s, fixed = TRUE)
    has_dot   <- grepl(".", s, fixed = TRUE)

    if (has_comma && has_dot) {
      s <- gsub(".", "", s, fixed = TRUE)
      s <- sub(",", ".", s, fixed = TRUE)
      return(sign * suppressWarnings(as.numeric(s)))
    }

    if (has_comma && !has_dot) {
      s <- gsub(",", ".", s, fixed = TRUE)
      return(sign * suppressWarnings(as.numeric(s)))
    }

    # Dot-grouped export artifacts (e.g., 1.006.367.484, 15.028.731).
    # Require at least two grouped blocks so values like 0.001 stay decimal.
    if (!has_comma && grepl("^\\d{1,3}(\\.\\d{3}){2,}$", s)) {
      digits <- gsub(".", "", s, fixed = TRUE)
      n <- nchar(digits)
      if (n >= 10) return(sign * suppressWarnings(as.numeric(digits) / 1e9))
      if (n == 9)  return(sign * suppressWarnings(as.numeric(digits) / 1e8))
      if (n %in% c(7, 8)) return(sign * suppressWarnings(as.numeric(digits) / 1e6))
      return(sign * suppressWarnings(as.numeric(digits)))
    }

    # Plain-digit locale artifacts without separators.
    if (grepl("^\\d{7,}$", s)) {
      n <- nchar(s)
      if (n >= 9) return(sign * suppressWarnings(as.numeric(s) / 1e9))
      if (n == 8) return(sign * suppressWarnings(as.numeric(s) / 1e8))
      if (n == 7) return(sign * suppressWarnings(as.numeric(s) / 1e6))
    }

    sign * suppressWarnings(as.numeric(s))
  }

  as.numeric(vapply(x, parse_one, numeric(1)))
}

fix_names_case_insensitive_unique <- function(nm) {
  nm <- trimws(nm)
  nm <- gsub("^\ufeff", "", nm)
  nm <- gsub('^"|"$', "", nm)
  nm[nm == "" | is.na(nm)] <- "X"
  key <- tolower(nm)
  if (any(duplicated(key))) {
    nm_new <- character(length(nm))
    seen <- list()
    for (i in seq_along(nm)) {
      k <- key[i]
      n_seen <- if (!is.null(seen[[k]])) seen[[k]] + 1L else 0L
      seen[[k]] <- n_seen
      nm_new[i] <- if (n_seen == 0L) nm[i] else paste0(nm[i], "_", n_seen)
    }
    nm <- nm_new
  }
  make.unique(nm, sep = "_")
}

detect_delim <- function(path) {
  x <- readLines(path, n = 5, warn = FALSE, encoding = "UTF-8")
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return(";")
  s <- paste(x, collapse = "\n")
  counts <- c(`;` = str_count(s, fixed(";")),
              `,` = str_count(s, fixed(",")),
              `\t` = str_count(s, fixed("\t")),
              `|` = str_count(s, fixed("|")))
  names(counts)[which.max(counts)]
}

read_any_csv <- function(path, guess_max = 200000, name_repair = "minimal") {
  delim <- detect_delim(path)
  x <- suppressWarnings(readr::read_delim(
    file = path,
    delim = delim,
    col_types = cols(.default = col_character()),
    trim_ws = TRUE,
    progress = FALSE,
    guess_max = guess_max,
    name_repair = name_repair
  ))

  # Fallback for malformed quoted CSV headers that collapse everything into one column.
  if (ncol(x) == 1 && delim == ",") {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    lines <- lines[!is.na(lines) & nzchar(lines)]
    if (length(lines) >= 2) {
      header <- gsub('^"|"$', "", lines[1])
      header <- gsub('""', '"', header, fixed = TRUE)
      cols <- strsplit(header, ",", fixed = TRUE)[[1]]
      cols <- gsub('^"|"$', "", trimws(cols))
      rows <- strsplit(lines[-1], ",", fixed = TRUE)
      max_len <- max(length(cols), max(vapply(rows, length, integer(1))))
      if (length(cols) < max_len) cols <- c(cols, paste0("V", seq_len(max_len - length(cols))))
      m <- do.call(rbind, lapply(rows, function(r) {
        if (length(r) < max_len) r <- c(r, rep(NA_character_, max_len - length(r)))
        r[seq_len(max_len)]
      }))
      x <- as.data.frame(m, stringsAsFactors = FALSE)
      names(x) <- cols[seq_len(ncol(x))]
    }
  }

  nm <- names(x)
  if (length(nm) > 0 &&
      (!nzchar(trimws(nm[1])) || is.na(nm[1])) &&
      !any(tolower(nm) %in% c("treeid", "id", "tree_id", "tree.id"))) {
    nm[1] <- "TreeID"
  }
  nm <- fix_names_case_insensitive_unique(nm)
  names(x) <- nm
  as.data.frame(x, stringsAsFactors = FALSE)
}

find_col_ci <- function(df, candidates) {
  nms <- names(df)
  low <- tolower(nms)
  cand_low <- tolower(candidates)
  hit <- match(cand_low, low)
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0) return(NA_character_)
  nms[hit[1]]
}

ensure_cols <- function(df, cols, context = "") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop("Missing required columns", if (nzchar(context)) paste0(" (", context, ")") else "", ": ", paste(miss, collapse = ", "))
  df
}

standardize_3dfin <- function(df) {
  names(df) <- fix_names_case_insensitive_unique(names(df))

  id_col <- find_col_ci(df, c("TreeID", "ID", "id", "tree_id"))
  if (is.na(id_col)) {
    id_col <- names(df)[1]
    if (!nzchar(id_col)) id_col <- "TreeID"
  }

  th_col  <- find_col_ci(df, c("TH", "height", "hoyde", "H"))
  dbh_col <- find_col_ci(df, c("DBH", "dbh", "d"))
  x_col0  <- find_col_ci(df, c("X", "x", "x_m", "dx", "E", "E_utm32", "Easting"))
  y_col0  <- find_col_ci(df, c("Y", "y", "y_m", "dy", "N", "N_utm32", "Northing"))

  if (is.na(th_col) || is.na(dbh_col) || is.na(x_col0) || is.na(y_col0)) {
    stop("Could not standardize 3DFin columns. Found names: ", paste(names(df), collapse = ", "))
  }

  out <- df %>%
    mutate(
      TreeID = as.character(.data[[id_col]]),
      TH  = to_num(.data[[th_col]]),
      DBH = to_num(.data[[dbh_col]]),
      X   = to_num(.data[[x_col0]]),
      Y   = to_num(.data[[y_col0]])
    ) %>%
    mutate(
      tree_id = .data$TreeID,
      x_local = .data$X,
      y_local = .data$Y
    )

  out
}

write_csv_safe <- function(df, path) {
  dir_create(dirname(path))
  readr::write_csv(df, file = path, na = "")
  invisible(path)
}

sanitize_for_gpkg <- function(df) {
  names(df) <- fix_names_case_insensitive_unique(names(df))
  df
}

cleanup_gpkg <- function(gpkg_path) {
  sidecars <- c(
    gpkg_path,
    paste0(gpkg_path, "-journal"),
    paste0(gpkg_path, "-wal"),
    paste0(gpkg_path, "-shm")
  )
  suppressWarnings(unlink(sidecars, force = TRUE))
  invisible(gpkg_path)
}

ensure_sf_crs <- function(x, target_epsg = 32632, context = "layer", assume_if_missing = TRUE) {
  if (is.null(x) || nrow(x) == 0) return(x)
  target_crs <- sf::st_crs(as.integer(target_epsg))
  if (is.na(target_crs)) stop("Invalid target CRS EPSG in ", context, ": ", target_epsg)

  cur_crs <- sf::st_crs(x)
  if (is.na(cur_crs)) {
    if (!isTRUE(assume_if_missing)) stop("Missing CRS in ", context, ".")
    sf::st_crs(x) <- target_crs
    return(x)
  }

  if (!isTRUE(cur_crs == target_crs)) {
    x <- sf::st_transform(x, target_crs)
  }
  x
}

write_gpkg_points <- function(df, xcol, ycol, gpkg, layer, crs = 32632, overwrite = FALSE) {
  dir_create(dirname(gpkg))
  if (overwrite) cleanup_gpkg(gpkg)
  df <- sanitize_for_gpkg(df)
  pts <- st_as_sf(df, coords = c(xcol, ycol), crs = crs, remove = FALSE)
  pts <- ensure_sf_crs(pts, target_epsg = crs, context = paste0("write_gpkg_points:", layer), assume_if_missing = TRUE)
  gpkg_tmp <- tempfile(pattern = "tmp_gpkg_", fileext = ".gpkg")
  cleanup_gpkg(gpkg_tmp)
  st_write(pts, gpkg_tmp, layer = layer, driver = "GPKG", quiet = TRUE, delete_dsn = TRUE)
  ok_copy <- file.copy(gpkg_tmp, gpkg, overwrite = TRUE)
  cleanup_gpkg(gpkg_tmp)
  if (!isTRUE(ok_copy)) stop("Could not copy temporary GPKG to target path: ", gpkg)
  invisible(gpkg)
}

#### 03 DBH IMPUTATION / INTERPOLATION ####

impute_dbh_rf <- function(raw,
                          id_col = "TreeID",
                          th_col = "TH",
                          dbh_col = "DBH",
                          x_col = "X",
                          y_col = "Y",
                          predictors = c("TH","X","Y"),
                          num_trees = 500,
                          mtry = 2,
                          min_node_size = 5,
                          importance = "impurity",
                          seed = 42,
                          dbh_min_valid = 0,
                          external_model = NULL,
                          external_predictors = c("TH"),
                          external_label = "tls_transfer",
                          verbose = FALSE) {
  raw[[dbh_col]] <- to_num(raw[[dbh_col]])
  raw <- raw %>%
    mutate(
      DBH_meas = ifelse(.data[[dbh_col]] <= dbh_min_valid | is.na(.data[[dbh_col]]), NA_real_, to_num(.data[[dbh_col]]))
    )
  raw$DBH_amended <- NA_real_
  imp_source <- rep(NA_character_, nrow(raw))

  if (all(!is.na(raw$DBH_meas))) {
    raw <- raw %>%
      mutate(
        DBH_amended = NA_real_,
        DBH_source = "measured",
        DBH        = .data$DBH_meas
      )
    return(list(data = raw, model = NULL, metrics = NULL))
  }

  train_rf <- raw %>% filter(!is.na(.data$DBH_meas), !is.na(.data[[th_col]]))
  preds_use <- predictors[predictors %in% names(raw) & predictors != dbh_col]
  if (!("TH" %in% preds_use)) preds_use <- c("TH", preds_use)
  preds_use <- unique(preds_use[preds_use %in% names(raw)])

  idx_missing <- which(is.na(raw$DBH_meas) & !is.na(raw[[th_col]]))
  if (!is.null(external_model) && length(idx_missing) > 0) {
    preds_ext <- external_predictors[external_predictors %in% names(raw) & external_predictors != dbh_col]
    if (!("TH" %in% preds_ext)) preds_ext <- c("TH", preds_ext)
    preds_ext <- unique(preds_ext[preds_ext %in% names(raw)])
    if (length(preds_ext) > 0) {
      ext_new <- raw[idx_missing, , drop = FALSE]
      ext_ok <- stats::complete.cases(ext_new[, preds_ext, drop = FALSE])
      if (any(ext_ok)) {
        pred_ext <- rep(NA_real_, length(idx_missing))
        pred_ext[ext_ok] <- predict(external_model, data = ext_new[ext_ok, , drop = FALSE])$predictions
        raw$DBH_amended[idx_missing] <- pred_ext
        imp_source[idx_missing[!is.na(pred_ext)]] <- external_label
      }
    }
  }

  if (nrow(train_rf) < 5 || length(preds_use) == 0) {
    min_obs <- if (nrow(train_rf) > 0) min(train_rf$DBH_meas, na.rm = TRUE) else NA_real_
    raw <- raw %>%
      mutate(
        DBH_amended = ifelse(is.na(.data$DBH_meas) & is.na(.data$DBH_amended) & !is.na(.data[[th_col]]) & is.finite(min_obs), min_obs, .data$DBH_amended)
      )
    imp_source[is.na(imp_source) & !is.na(raw$DBH_amended) & is.na(raw$DBH_meas)] <- "fallback_min"
    raw <- raw %>%
      mutate(
        DBH_source = case_when(
          !is.na(.data$DBH_meas) ~ "measured",
          is.na(.data$DBH_meas) & !is.na(.data$DBH_amended) ~ imp_source,
          TRUE ~ "missing"
        ),
        DBH = coalesce(.data$DBH_meas, .data$DBH_amended)
      )
    return(list(data = raw, model = NULL, metrics = list(note = "Too few training rows; filled with min observed DBH")))
  }

  set.seed(seed)
  rf_formula <- as.formula(paste("DBH_meas ~", paste(preds_use, collapse = " + ")))

  rf_model <- ranger::ranger(
    formula        = rf_formula,
    data           = train_rf,
    num.trees      = num_trees,
    mtry           = min(max(1L, mtry), length(preds_use)),
    min.node.size  = min_node_size,
    importance     = importance,
    seed           = seed
  )

  oob_rmse <- sqrt(rf_model$prediction.error)
  mean_dbh <- mean(train_rf$DBH_meas, na.rm = TRUE)
  rel_rmse <- if (is.finite(mean_dbh) && mean_dbh > 0) 100 * oob_rmse / mean_dbh else NA_real_

  idx_missing <- which(is.na(raw$DBH_meas) & is.na(raw$DBH_amended) & !is.na(raw[[th_col]]))
  if (length(idx_missing) > 0) {
    rf_new <- raw[idx_missing, , drop = FALSE]
    ok_pred <- stats::complete.cases(rf_new[, preds_use, drop = FALSE])
    if (any(ok_pred)) {
      pred_rf <- rep(NA_real_, length(idx_missing))
      pred_rf[ok_pred] <- predict(rf_model, data = rf_new[ok_pred, , drop = FALSE])$predictions
      raw$DBH_amended[idx_missing] <- pred_rf
      imp_source[idx_missing[!is.na(pred_rf)]] <- "rf_imputed"
    }
  }

  if (any(is.na(raw$DBH_amended) & is.na(raw$DBH_meas))) {
    min_obs <- min(train_rf$DBH_meas, na.rm = TRUE)
    if (is.finite(min_obs)) {
      idx_fill <- which(is.na(raw$DBH_amended) & is.na(raw$DBH_meas))
      raw$DBH_amended[idx_fill] <- min_obs
      imp_source[idx_fill] <- "fallback_min"
    }
  }

  raw <- raw %>%
    mutate(
      DBH_source = case_when(
        !is.na(.data$DBH_meas) ~ "measured",
        is.na(.data$DBH_meas) & !is.na(.data$DBH_amended) ~ imp_source,
        TRUE ~ "missing"
      ),
      DBH = coalesce(.data$DBH_meas, .data$DBH_amended)
    )

  if (verbose) {
    message(sprintf("RF DBH imputation: OOB-RMSE=%.4f m (%.1f%% of mean DBH)", oob_rmse, rel_rmse))
    message("DBH_source table:\n", paste(capture.output(print(table(raw$DBH_source))), collapse = "\n"))
  }

  list(
    data = raw,
    model = rf_model,
    metrics = list(oob_rmse = oob_rmse, rel_rmse = rel_rmse, mean_dbh = mean_dbh, predictors = preds_use)
  )
}

fit_dbh_transfer_model <- function(raw,
                                   th_col = "TH",
                                   predictors = c("TH"),
                                   num_trees = 500,
                                   mtry = 1,
                                   min_node_size = 5,
                                   importance = "impurity",
                                   seed = 42,
                                   dbh_min_valid = 0,
                                   verbose = FALSE) {
  if (!"DBH_meas" %in% names(raw)) {
    raw$DBH <- to_num(raw$DBH)
    raw$DBH_meas <- ifelse(is.na(raw$DBH) | raw$DBH <= dbh_min_valid, NA_real_, raw$DBH)
  }

  preds_use <- predictors[predictors %in% names(raw)]
  if (!("TH" %in% preds_use)) preds_use <- c("TH", preds_use)
  preds_use <- unique(preds_use[preds_use %in% names(raw)])
  if (length(preds_use) == 0) return(NULL)

  train_rf <- raw %>% filter(!is.na(.data$DBH_meas), !is.na(.data[[th_col]]))
  if (nrow(train_rf) < 5) return(NULL)

  set.seed(seed)
  rf_formula <- as.formula(paste("DBH_meas ~", paste(preds_use, collapse = " + ")))
  model <- ranger::ranger(
    formula = rf_formula,
    data = train_rf,
    num.trees = num_trees,
    mtry = min(max(1L, mtry), length(preds_use)),
    min.node.size = min_node_size,
    importance = importance,
    seed = seed
  )

  oob_rmse <- sqrt(model$prediction.error)
  mean_dbh <- mean(train_rf$DBH_meas, na.rm = TRUE)
  rel_rmse <- if (is.finite(mean_dbh) && mean_dbh > 0) 100 * oob_rmse / mean_dbh else NA_real_
  if (verbose) {
    message(sprintf("TLS transfer model (TH->DBH): OOB-RMSE=%.4f m (%.1f%% of mean DBH)", oob_rmse, rel_rmse))
  }

  list(
    model = model,
    predictors = preds_use,
    metrics = list(oob_rmse = oob_rmse, rel_rmse = rel_rmse, mean_dbh = mean_dbh)
  )
}

#### 04 TREECOMP THINNING DECISIONS ####

add_id_safe <- function(df) {
  id_candidates <- c("ID", "TreeID", "id", "tree_id")
  has_id <- id_candidates[id_candidates %in% names(df)]
  if (length(has_id) == 0) df$.id <- seq_len(nrow(df)) else df$.id <- df[[has_id[1]]]
  df
}

pick_first_coords <- function(nms, base = c("x","y")) {
  out <- list()
  for (b in base) {
    cand <- c(b, paste0(b, c(".x",".y")), toupper(b), paste0(toupper(b), c(".x",".y")))
    out[[b]] <- intersect(cand, nms)
  }
  out
}

local_rank_one_index <- function(df_sf, idx_col, radius, min_n = 3) {
  val <- df_sf[[idx_col]]
  n   <- nrow(df_sf)
  rank_vec <- rep(NA_real_, n)
  neigh_list <- sf::st_is_within_distance(df_sf, df_sf, dist = radius)
  for (i in seq_len(n)) {
    neigh_idx <- as.integer(neigh_list[[i]])
    neigh_idx <- neigh_idx[!is.na(val[neigh_idx])]
    if (length(neigh_idx) < min_n || is.na(val[i])) {
      rank_vec[i] <- NA_real_
      next
    }
    neigh_vals <- val[neigh_idx]
    rank_vec[i] <- mean(neigh_vals <= val[i], na.rm = TRUE)
  }
  df_sf[[paste0(idx_col, "_rank_local")]] <- rank_vec
  df_sf
}

compute_treecomp_thinning <- function(raw,
                                      ci_methods,
                                      compete_radius,
                                      target_radius,
                                      target_source = "exclude_edge",
                                      local_radius,
                                      local_q,
                                      min_neigh,
                                      do_plots = FALSE,
                                      plot_dir = NULL,
                                      verbose = FALSE) {
  raw_tc <- raw
  id_like <- names(raw_tc)[tolower(names(raw_tc)) %in% c("id", "treeid", "tree_id", "tree.id")]
  if (length(id_like) > 1) {
    keep_id <- if ("TreeID" %in% id_like) "TreeID" else id_like[1]
    drop_id <- setdiff(id_like, keep_id)
    raw_tc <- raw_tc %>% select(-any_of(drop_id))
  }
  dbh_pref <- c("DBH_meas", "DBH_amended", "DBH")
  dbh_pref <- dbh_pref[dbh_pref %in% names(raw_tc)]
  if (length(dbh_pref) == 0) stop("No DBH columns available for TreeCompR (expected DBH_meas/DBH_amended/DBH).")
  raw_tc$DBH <- to_num(raw_tc[[dbh_pref[1]]])
  if (length(dbh_pref) > 1) {
    for (cc in dbh_pref[-1]) raw_tc$DBH <- dplyr::coalesce(raw_tc$DBH, to_num(raw_tc[[cc]]))
  }

  # TreeCompR expects x/y/dbh/height columns; keep standard names from chunk.
  inv <- TreeCompR::read_inv(
    inv_source  = raw_tc,
    x           = X,
    y           = Y,
    dbh         = DBH,
    height      = TH,
    dbh_unit    = "m",
    height_unit = "m",
    verbose     = verbose
  )

  targets <- TreeCompR::define_target(inv, target_source = target_source, radius = target_radius)

  ci <- TreeCompR::compete_inv(
    inv    = targets,
    radius = compete_radius,
    method = ci_methods
  )

  targets_df <- targets |>
    add_id_safe() |>
    as_tibble() |>
    select(any_of(c(".id","x","y")))

  ci_df <- ci |>
    add_id_safe() |>
    as_tibble()

  plot_df <- targets_df |>
    left_join(ci_df, by = ".id")

  coords <- pick_first_coords(names(plot_df), c("x","y"))
  if (length(coords$x) == 0 || length(coords$y) == 0) stop("Could not find x/y columns after TreeCompR join.")
  plot_df <- plot_df %>%
    mutate(
      x = dplyr::coalesce(!!!rlang::syms(coords$x)),
      y = dplyr::coalesce(!!!rlang::syms(coords$y))
    ) %>%
    select(-any_of(setdiff(unique(c(coords$x, coords$y)), c("x","y")))) %>%
    as_tibble()

  base_cols <- c(".id","x","y","ID","TreeID","id","tree_id")
  ci_cols <- setdiff(names(plot_df), base_cols)
  ci_cols <- ci_cols[grepl("^CI_", ci_cols)]
  ci_cols <- ci_cols[vapply(select(plot_df, all_of(ci_cols)), is.numeric, logical(1))]
  if (length(ci_cols) == 0) stop("No CI_* numeric columns found.")

  pts <- sf::st_as_sf(plot_df, coords = c("x","y"), remove = FALSE)

  pts_ranked <- purrr::reduce(
    ci_cols,
    .init = pts,
    .f = function(acc, idx) local_rank_one_index(acc, idx, radius = local_radius, min_n = min_neigh)
  )

  out <- sf::st_drop_geometry(pts_ranked)
  out <- out %>% select(-any_of(c("id", "dbh", "height", "ID", "TreeID", "tree_id", "DBH")))

  for (idx in ci_cols) {
    rank_col <- paste0(idx, "_rank_local")
    thin_col <- paste0("thin_", idx)
    if (!rank_col %in% names(out)) next
    out[[thin_col]] <- dplyr::case_when(
      is.na(out[[rank_col]])      ~ NA_character_,
      out[[rank_col]] >= local_q  ~ "remove",
      TRUE                        ~ "keep"
    )
  }

  # Join DBH info back.
  id_col <- find_col_ci(raw, c("TreeID","ID","tree_id","id"))
  if (is.na(id_col)) stop("Could not find ID column for joining TreeCompR output.")
  dbh_info <- raw %>%
    mutate(.id_join = as.character(.data[[id_col]])) %>%
    select(any_of(c(".id_join","TH","DBH_meas","DBH_amended","DBH_source")))

  out <- out %>%
    mutate(.id_join = as.character(.data$.id)) %>%
    left_join(dbh_info, by = ".id_join") %>%
    select(-".id_join")

  # Ensure required output columns exist, harmonize to requested names
  out <- out %>%
    mutate(
      tree_id = as.character(.data$.id),
      x_local = to_num(.data$x),
      y_local = to_num(.data$y)
    )

  # Optional plots (local rank maps)
  if (do_plots) {
    if (is.null(plot_dir)) plot_dir <- file.path(getwd(), "CI_plots_local_rank")
    dir_create(plot_dir)

    plot_ci_map_local_rank <- function(df, idx, q = 0.80) {
      rank_col <- paste0(idx, "_rank_local")
      df <- df %>%
        mutate(
          flag = dplyr::case_when(
            is.na(.data[[rank_col]]) ~ "NA",
            .data[[rank_col]] >= q   ~ paste0(">= local ", q),
            TRUE                     ~ paste0("< local ", q)
          )
        )
      lev_low  <- paste0("< local ", q)
      lev_high <- paste0(">= local ", q)
      df$flag  <- factor(df$flag, levels = c(lev_low, lev_high, "NA"))
      shape_vals <- c(16, 4, 1)
      names(shape_vals) <- c(lev_low, lev_high, "NA")

      ggplot(df, aes(x = .data$x, y = .data$y)) +
        geom_point(aes(shape = .data$flag, color = .data$flag), size = 2, stroke = 0.9, alpha = 0.9) +
        scale_shape_manual(values = shape_vals, na.translate = FALSE) +
        coord_equal() +
        labs(
          title    = paste0("Competition index (local rank): ", idx),
          subtitle = paste0("Neighbourhood r = ", local_radius, " m | Cross = worst ", round(100 * (1 - q)), "%"),
          x = "x (m)", y = "y (m)"
        )
    }

    for (idx in ci_cols) {
      rank_col <- paste0(idx, "_rank_local")
      if (!rank_col %in% names(out) || all(is.na(out[[rank_col]]))) next
      p <- plot_ci_map_local_rank(out, idx, q = local_q)
      ggsave(
        filename = file.path(plot_dir, paste0("plot_local_rank_", idx, ".png")),
        plot = p, width = 7, height = 6, dpi = 300
      )
    }
  }

  # Final column ordering (keep everything but lead with required)
  lead <- c("tree_id","x_local","y_local","TH","DBH_meas","DBH_amended","DBH_source")
  keep_lead <- lead[lead %in% names(out)]
  out <- out %>% select(all_of(keep_lead), everything())

  list(data = out, ci_cols = ci_cols)
}

#### 05 COORDINATE CONVERSION TO UTM32 ####

local_to_utm32 <- function(df, x_col = "x_local", y_col = "y_local", E0, N0) {
  ensure_cols(df, c(x_col, y_col), context = "local_to_utm32")
  df %>%
    mutate(
      E_utm32 = E0 + to_num(.data[[x_col]]),
      N_utm32 = N0 + to_num(.data[[y_col]])
    )
}

#### 06 HAND LOAD / CONVERSION ####

hand_polar_to_xy <- function(df,
                             angle_col = "retning",
                             dist_col = "avstand",
                             dist_divisor = 100,
                             X_center = 0,
                             Y_center = 0) {
  ensure_cols(df, c(angle_col, dist_col), context = "hand_polar_to_xy")
  ang <- to_num(df[[angle_col]])
  dist <- to_num(df[[dist_col]])
  dist_m <- dist / dist_divisor

  theta_rad <- ang * pi / 180
  dx <- dist_m * sin(theta_rad)
  dy <- dist_m * cos(theta_rad)

  df %>%
    mutate(
      avstand_m = dist_m,
      x_m = dx,
      y_m = dy,
      X = X_center + dx,
      Y = Y_center + dy
    )
}

load_hand <- function(hand_in,
                      hand_sheet = "Sheet1",
                      crs_out = 32632,
                      crs_hand_raw = 25833,
                      E0 = NA_real_, N0 = NA_real_,
                      X_center = 0, Y_center = 0,
                      angle_col = "retning",
                      dist_col = "avstand",
                      dist_divisor = 100) {
  ext <- tolower(tools::file_ext(hand_in))

  if (ext %in% c("xlsx","xls")) {
    raw <- readxl::read_excel(hand_in, sheet = hand_sheet)
    raw <- as.data.frame(raw, stringsAsFactors = FALSE)
    raw <- hand_polar_to_xy(raw, angle_col = angle_col, dist_col = dist_col, dist_divisor = dist_divisor, X_center = X_center, Y_center = Y_center)
    sf_raw <- st_as_sf(raw, coords = c("X","Y"), crs = crs_hand_raw, remove = FALSE)
    sf_out <- st_transform(sf_raw, crs_out)
    df_out <- sf::st_drop_geometry(sf_out)
    coords <- sf::st_coordinates(sf_out)
    df_out$E_utm32 <- coords[,1]
    df_out$N_utm32 <- coords[,2]
    df_out$tree_id <- if ("TreeID" %in% names(df_out)) as.character(df_out$TreeID) else if ("trenr" %in% names(df_out)) as.character(df_out$trenr) else as.character(seq_len(nrow(df_out)))
    return(df_out)
  }

  df <- read_any_csv(hand_in, guess_max = guess_max, name_repair = name_repair)

  # Try direct UTM columns
  e_col <- find_col_ci(df, c("E_utm32","E_base","E","Easting","east","e"))
  n_col <- find_col_ci(df, c("N_utm32","N_base","N","Northing","north","n"))

  if (!is.na(e_col) && !is.na(n_col)) {
    df <- df %>%
      mutate(
        E_utm32 = to_num(.data[[e_col]]),
        N_utm32 = to_num(.data[[n_col]])
      )
    df$tree_id <- if ("TreeID" %in% names(df)) as.character(df$TreeID) else if ("trenr" %in% names(df)) as.character(df$trenr) else as.character(seq_len(nrow(df)))
    return(df)
  }

  # Try local x/y in meters (x_m/y_m or X/Y)
  xm_col <- find_col_ci(df, c("x_m","x","X"))
  ym_col <- find_col_ci(df, c("y_m","y","Y"))

  if (!is.na(xm_col) && !is.na(ym_col) && is.finite(E0) && is.finite(N0)) {
    df <- df %>%
      mutate(
        x_local = to_num(.data[[xm_col]]),
        y_local = to_num(.data[[ym_col]]),
        E_utm32 = E0 + .data$x_local,
        N_utm32 = N0 + .data$y_local
      )
    df$tree_id <- if ("TreeID" %in% names(df)) as.character(df$TreeID) else if ("trenr" %in% names(df)) as.character(df$trenr) else as.character(seq_len(nrow(df)))
    return(df)
  }

  # Try polar columns in CSV
  if (angle_col %in% names(df) && dist_col %in% names(df)) {
    df <- hand_polar_to_xy(df, angle_col = angle_col, dist_col = dist_col, dist_divisor = dist_divisor, X_center = X_center, Y_center = Y_center)
    sf_raw <- st_as_sf(df, coords = c("X","Y"), crs = crs_hand_raw, remove = FALSE)
    sf_out <- st_transform(sf_raw, crs_out)
    df_out <- st_drop_geometry(sf_out)
    coords <- st_coordinates(sf_out)
    df_out$E_utm32 <- coords[,1]
    df_out$N_utm32 <- coords[,2]
    df_out$tree_id <- if ("TreeID" %in% names(df_out)) as.character(df_out$TreeID) else if ("trenr" %in% names(df_out)) as.character(df_out$trenr) else as.character(seq_len(nrow(df_out)))
    return(df_out)
  }

  stop("Could not load/convert hand input. Provide UTM columns, or local x/y with E0/N0 set, or polar retning/avstand.")
}

#### 07 ALIGNMENT / "AFFINE" TRANSFORMATION AGAINST HAND ####

rename_reserved <- function(df) {
  i <- which(tolower(names(df)) == "fid")
  if (length(i) > 0) names(df)[i] <- "fid_src"
  df
}

bbox_from_xy <- function(xy, crs = crs_out) {
  xy <- xy[is.finite(xy[,1]) & is.finite(xy[,2]), , drop = FALSE]
  st_bbox(c(
    xmin = min(xy[,1], na.rm = TRUE),
    ymin = min(xy[,2], na.rm = TRUE),
    xmax = max(xy[,1], na.rm = TRUE),
    ymax = max(xy[,2], na.rm = TRUE)
  ), crs = st_crs(crs))
}

clip_to_bbox_xy <- function(df, xy, bbox) {
  ok <- is.finite(xy[,1]) & is.finite(xy[,2]) &
    xy[,1] >= bbox["xmin"] & xy[,1] <= bbox["xmax"] &
    xy[,2] >= bbox["ymin"] & xy[,2] <= bbox["ymax"]
  list(df = df[ok, , drop = FALSE], xy = xy[ok, , drop = FALSE], keep = ok)
}

median_nn <- function(pts) {
  pts <- pts[is.finite(pts[,1]) & is.finite(pts[,2]), , drop = FALSE]
  if (nrow(pts) < 3) return(NA_real_)
  nn <- RANN::nn2(data = pts, query = pts, k = 2)
  stats::median(nn$nn.dists[,2], na.rm = TRUE)
}

make_T <- function(A, t) {
  rbind(
    c(A[1,1], A[1,2], t[1]),
    c(A[2,1], A[2,2], t[2]),
    c(0, 0, 1)
  )
}

apply_T <- function(pts, T) {
  X <- cbind(pts, 1)
  Y <- X %*% t(T)
  Y[, 1:2, drop = FALSE]
}

pca_angle <- function(pts) {
  pts <- pts[is.finite(pts[,1]) & is.finite(pts[,2]), , drop = FALSE]
  if (nrow(pts) < 3) return(0)
  S <- cov(pts)
  ev <- eigen(S, symmetric = TRUE)
  v1 <- ev$vectors[, 1]
  atan2(v1[2], v1[1])
}

voxel_thin_xy <- function(xy, cell) {
  xy <- as.matrix(xy)
  n <- nrow(xy)
  if (n < 1) return(logical(0))
  if (!is.finite(cell) || cell <= 0) return(rep(TRUE, n))
  ok <- is.finite(xy[,1]) & is.finite(xy[,2])
  keep <- rep(FALSE, n)
  if (!any(ok)) return(keep)
  kx <- floor(xy[ok,1] / cell)
  ky <- floor(xy[ok,2] / cell)
  key <- paste(kx, ky, sep = "_")
  keep[which(ok)] <- !duplicated(key)
  keep
}

thin_fit_idx <- function(xy, cell, max_points = Inf) {
  n <- nrow(xy)
  if (n < 1) return(integer(0))
  keep <- voxel_thin_xy(xy, cell = cell)
  idx <- which(keep & is.finite(xy[,1]) & is.finite(xy[,2]))
  if (length(idx) < 1) idx <- which(is.finite(xy[,1]) & is.finite(xy[,2]))
  if (!is.finite(max_points) || max_points < 1) max_points <- length(idx)
  if (length(idx) > max_points) {
    pick <- unique(round(seq(1, length(idx), length.out = max_points)))
    idx <- idx[pick]
  }
  idx
}

decompose_similarity <- function(T) {
  A <- T[1:2, 1:2, drop = FALSE]
  sc <- sqrt(abs(det(A)))
  if (!is.finite(sc) || sc <= 0) sc <- NA_real_
  R <- if (is.finite(sc) && sc > 0) A / sc else diag(2)
  ang <- atan2(R[2,1], R[1,1])
  list(
    scale = sc,
    rot_rad = ang,
    rot_deg = ang * 180 / pi,
    tx = T[1,3],
    ty = T[2,3]
  )
}

initial_similarity <- function(src, ref, allow_scale = TRUE, scale_min = -Inf, scale_max = Inf) {
  cs <- colMeans(src, na.rm = TRUE)
  cr <- colMeans(ref, na.rm = TRUE)

  ds <- sqrt(rowSums((src - matrix(cs, nrow(src), 2, byrow = TRUE))^2))
  dr <- sqrt(rowSums((ref - matrix(cr, nrow(ref), 2, byrow = TRUE))^2))
  ss <- stats::median(ds[is.finite(ds)], na.rm = TRUE); if (!is.finite(ss) || ss == 0) ss <- 1
  sr <- stats::median(dr[is.finite(dr)], na.rm = TRUE); if (!is.finite(sr) || sr == 0) sr <- 1
  sc <- if (isTRUE(allow_scale)) sr / ss else 1.0
  if (is.finite(scale_min)) sc <- max(sc, scale_min)
  if (is.finite(scale_max)) sc <- min(sc, scale_max)

  ang_s <- pca_angle(scale(src, center = cs, scale = FALSE))
  ang_r <- pca_angle(scale(ref, center = cr, scale = FALSE))
  ang <- ang_r - ang_s

  R <- matrix(c(cos(ang), -sin(ang),
                sin(ang),  cos(ang)), nrow = 2, byrow = TRUE)

  A <- sc * R
  t <- cr - as.numeric(A %*% cs)
  make_T(A, t)
}

fit_similarity <- function(src, dst, allow_scale = TRUE) {
  src <- as.matrix(src); dst <- as.matrix(dst)
  n <- nrow(src)
  if (n < 2) stop("Too few pairs for similarity fit.")

  mu_x <- colMeans(src)
  mu_y <- colMeans(dst)
  Xc <- src - matrix(mu_x, n, 2, byrow = TRUE)
  Yc <- dst - matrix(mu_y, n, 2, byrow = TRUE)

  C <- (t(Yc) %*% Xc) / n
  sv <- svd(C)
  U <- sv$u
  V <- sv$v
  D <- diag(c(1, sign(det(U %*% t(V)))))

  R <- U %*% D %*% t(V)
  var_x <- sum(Xc^2) / n
  s <- if (isTRUE(allow_scale) && is.finite(var_x) && var_x > 0) sum(diag(diag(sv$d) %*% D)) / var_x else 1.0
  if (!is.finite(s) || s <= 0) s <- 1.0

  A <- s * R
  tvec <- mu_y - as.numeric(A %*% mu_x)
  make_T(A, tvec)
}

fit_similarity_centered <- function(src_xy, dst_xy, allow_scale = TRUE) {
  src_xy <- as.matrix(src_xy)
  dst_xy <- as.matrix(dst_xy)
  c0 <- colMeans(dst_xy, na.rm = TRUE)
  src0 <- sweep(src_xy, 2, c0, "-")
  dst0 <- sweep(dst_xy, 2, c0, "-")
  T0 <- fit_similarity(src0, dst0, allow_scale = allow_scale)
  A <- T0[1:2,1:2,drop=FALSE]
  t0 <- T0[1:2,3,drop=TRUE]
  t <- t0 + c0 - as.vector(A %*% c0)
  T0[1:2,3] <- t
  T0
}

mutual_pairs <- function(src_t, ref, dmax) {
  nn_sr <- RANN::nn2(data = ref, query = src_t, k = 1)
  idx_r <- nn_sr$nn.idx[,1]
  dist  <- nn_sr$nn.dists[,1]

  keep <- which(is.finite(dist) & dist <= dmax)
  if (length(keep) == 0) return(NULL)

  src_keep <- keep
  ref_idx  <- idx_r[keep]

  nn_rs <- RANN::nn2(data = src_t, query = ref[ref_idx, , drop = FALSE], k = 1)
  back <- nn_rs$nn.idx[,1]

  ok <- which(back == src_keep)
  if (length(ok) == 0) return(NULL)

  data.frame(
    i_src = src_keep[ok],
    i_ref = ref_idx[ok],
    dist  = dist[src_keep[ok]]
  )
}

icp_similarity_core <- function(src, ref,
                                max_iter = 60,
                                dmax_start = 35, dmax_end = 2.5,
                                tol = 1e-4, min_pairs = 12,
                                trim_frac = 0.80,
                                allow_scale = TRUE,
                                init_T = NULL,
                                scale_min = -Inf,
                                scale_max = Inf) {
  src <- as.matrix(src)
  ref <- as.matrix(ref)
  if (nrow(src) < 2 || nrow(ref) < 2) {
    T_id <- diag(3)
    return(list(T = T_id, mean_err = NA_real_, n_pairs_final = 0L, n_iter = 0L))
  }

  T <- if (is.null(init_T)) {
    initial_similarity(src, ref, allow_scale = allow_scale, scale_min = scale_min, scale_max = scale_max)
  } else {
    init_T
  }

  prev_err <- Inf
  err_last <- NA_real_
  n_pairs_last <- 0L
  n_iter_used <- 0L

  for (it in seq_len(max_iter)) {
    dmax <- dmax_start + (dmax_end - dmax_start) * (it - 1) / max(1, (max_iter - 1))

    src_t <- apply_T(src, T)
    pairs <- mutual_pairs(src_t, ref, dmax = dmax)
    if (is.null(pairs) || nrow(pairs) < min_pairs) break

    if (is.finite(trim_frac) && trim_frac > 0 && trim_frac < 1 && nrow(pairs) > min_pairs) {
      k <- max(min_pairs, floor(trim_frac * nrow(pairs)))
      pairs <- pairs[order(pairs$dist), , drop = FALSE][seq_len(k), , drop = FALSE]
    }

    s_use <- src_t[pairs$i_src, , drop = FALSE]
    r_use <- ref[pairs$i_ref, , drop = FALSE]

    T_delta <- fit_similarity_centered(s_use, r_use, allow_scale = allow_scale)
    T <- T_delta %*% T

    err <- mean(pairs$dist, na.rm = TRUE)
    err_last <- err
    n_pairs_last <- nrow(pairs)
    n_iter_used <- it
    if (is.finite(prev_err) && is.finite(err) && abs(prev_err - err) < tol) break
    prev_err <- err
  }

  if (!is.finite(err_last)) {
    src_t <- apply_T(src, T)
    pairs <- mutual_pairs(src_t, ref, dmax = dmax_end)
    if (!is.null(pairs) && nrow(pairs) > 0) {
      err_last <- mean(pairs$dist, na.rm = TRUE)
      n_pairs_last <- nrow(pairs)
    }
  }

  list(
    T = T,
    mean_err = err_last,
    n_pairs_final = as.integer(n_pairs_last),
    n_iter = as.integer(n_iter_used)
  )
}

icp_similarity <- function(src, ref,
                           max_iter = 60,
                           dmax_start = 35, dmax_end = 2.5,
                           tol = 1e-4, min_pairs = 12,
                           trim_frac = 0.80,
                           scale_min = 0.85,
                           scale_max = 1.15) {
  src <- as.matrix(src)
  ref <- as.matrix(ref)
  ok_src <- is.finite(src[,1]) & is.finite(src[,2])
  ok_ref <- is.finite(ref[,1]) & is.finite(ref[,2])
  src <- src[ok_src, , drop = FALSE]
  ref <- ref[ok_ref, , drop = FALSE]

  if (nrow(src) < 2 || nrow(ref) < 2) {
    T_id <- diag(3)
    return(list(
      T = T_id,
      aligned = src,
      mean_err = NA_real_,
      n_pairs_final = 0L,
      n_iter = 0L,
      scale_total = 1.0,
      rotation_deg = 0.0,
      warn_scale_clamped = FALSE,
      fallback_rigid = FALSE,
      warn_flags = "warn_too_few_points"
    ))
  }

  c0 <- colMeans(ref, na.rm = TRUE)
  src0 <- sweep(src, 2, c0, "-")
  ref0 <- sweep(ref, 2, c0, "-")

  fit1 <- icp_similarity_core(
    src = src0, ref = ref0,
    max_iter = max_iter,
    dmax_start = dmax_start, dmax_end = dmax_end,
    tol = tol, min_pairs = min_pairs,
    trim_frac = trim_frac,
    allow_scale = TRUE,
    scale_min = scale_min,
    scale_max = scale_max
  )

  fit_use <- fit1
  info1 <- decompose_similarity(fit1$T)
  warn_flags <- character(0)
  fallback_rigid <- FALSE
  warn_scale <- !is.finite(info1$scale) || info1$scale < scale_min || info1$scale > scale_max

  if (warn_scale) {
    warn_flags <- c(warn_flags, "warn_scale_out_of_bounds")
    warning("ICP similarity scale outside guardrail [", scale_min, ", ", scale_max, "]; falling back to rigid (scale=1).")
    fit_use <- icp_similarity_core(
      src = src0, ref = ref0,
      max_iter = max_iter,
      dmax_start = dmax_start, dmax_end = dmax_end,
      tol = tol, min_pairs = min_pairs,
      trim_frac = trim_frac,
      allow_scale = FALSE
    )
    fallback_rigid <- TRUE
    warn_flags <- c(warn_flags, "fallback_rigid_scale_1")
  }

  T0 <- fit_use$T
  A <- T0[1:2,1:2,drop=FALSE]
  t0 <- T0[1:2,3,drop=TRUE]
  t_world <- t0 + c0 - as.vector(A %*% c0)
  T <- T0
  T[1:2,3] <- t_world

  info <- decompose_similarity(T)

  list(
    T = T,
    aligned = apply_T(src, T),
    mean_err = fit_use$mean_err,
    n_pairs_final = fit_use$n_pairs_final,
    n_iter = fit_use$n_iter,
    scale_total = info$scale,
    rotation_deg = info$rot_deg,
    warn_scale_clamped = warn_scale,
    fallback_rigid = fallback_rigid,
    warn_flags = paste(unique(warn_flags), collapse = "|")
  )
}

unique_match_pairs <- function(src_xy, hand_xy, dmax) {
  if (nrow(src_xy) < 1 || nrow(hand_xy) < 1) return(tibble(i_src = integer(), i_hand = integer(), dist = numeric()))
  nn <- RANN::nn2(data = src_xy, query = hand_xy, k = 1)
  pairs <- data.frame(
    i_hand = seq_len(nrow(hand_xy)),
    i_src  = nn$nn.idx[,1],
    dist   = nn$nn.dists[,1]
  )
  pairs <- pairs[is.finite(pairs$dist) & pairs$dist <= dmax, , drop = FALSE]
  if (nrow(pairs) == 0) return(tibble(i_src = integer(), i_hand = integer(), dist = numeric()))

  pairs <- pairs[order(pairs$dist), , drop = FALSE]
  used_hand <- rep(FALSE, nrow(hand_xy))
  used_src  <- rep(FALSE, nrow(src_xy))
  out <- vector("list", nrow(pairs))
  out_n <- 0L

  for (k in seq_len(nrow(pairs))) {
    h <- pairs$i_hand[k]
    s <- pairs$i_src[k]
    if (!used_hand[h] && !used_src[s]) {
      used_hand[h] <- TRUE
      used_src[s]  <- TRUE
      out_n <- out_n + 1L
      out[[out_n]] <- list(i_src = s, i_hand = h, dist = pairs$dist[k])
    }
  }
  if (out_n < 1L) return(tibble(i_src = integer(), i_hand = integer(), dist = numeric()))
  out <- bind_rows(out[seq_len(out_n)])
  out %>% arrange(.data$dist)
}

unique_match_keep <- function(src_xy, hand_xy, dmax) {
  pairs <- unique_match_pairs(src_xy, hand_xy, dmax)
  keep_src  <- rep(FALSE, nrow(src_xy))
  if (nrow(pairs) > 0) keep_src[pairs$i_src] <- TRUE
  keep_src
}

align_expand_cut <- function(src_df, src_xy, hand_xy, hand_bbox, hand_nn_med,
                            fit_voxel_cell_m = icp_fit_voxel_cell_m,
                            fit_max_points   = icp_fit_max_points,
                            max_iter         = icp_max_iter,
                            dmax_start       = icp_dmax_start,
                            dmax_end         = icp_dmax_end,
                            tol              = icp_tol,
                            min_pairs        = icp_min_pairs,
                            trim_frac        = icp_trim_frac,
                            scale_min        = icp_scale_min,
                            scale_max        = icp_scale_max,
                            match_dmax       = match_dmax_m) {
  warn_flags <- character(0)
  fit_cell <- fit_voxel_cell_m
  if (!is.finite(fit_cell) || fit_cell <= 0) {
    fit_cell <- if (is.finite(hand_nn_med) && hand_nn_med > 0) hand_nn_med else 1.0
  }

  src_fit_idx <- thin_fit_idx(src_xy, cell = fit_cell, max_points = fit_max_points)
  ref_fit_idx <- thin_fit_idx(hand_xy, cell = fit_cell, max_points = fit_max_points)
  src_fit <- src_xy[src_fit_idx, , drop = FALSE]
  ref_fit <- hand_xy[ref_fit_idx, , drop = FALSE]

  if (nrow(src_fit) < min_pairs || nrow(ref_fit) < min_pairs) {
    warn_flags <- c(warn_flags, "warn_fit_subset_too_small_fallback_all_points")
    src_fit <- src_xy
    ref_fit <- hand_xy
  }

  res <- icp_similarity(
    src = src_fit, ref = ref_fit,
    max_iter = max_iter,
    dmax_start = dmax_start, dmax_end = dmax_end,
    tol = tol, min_pairs = min_pairs,
    trim_frac = trim_frac,
    scale_min = scale_min, scale_max = scale_max
  )

  T <- res$T
  warn_flags <- c(warn_flags, res$warn_flags)

  xy_all <- apply_T(src_xy, T)
  cut2 <- clip_to_bbox_xy(src_df, xy_all, hand_bbox)
  df_cut <- cut2$df
  xy_cut <- cut2$xy

  if (nrow(xy_cut) < 1) {
    warn_flags <- c(warn_flags, "warn_bbox_empty_after_alignment_use_identity_clip")
    cut2 <- clip_to_bbox_xy(src_df, src_xy, hand_bbox)
    df_cut <- cut2$df
    xy_cut <- cut2$xy
    T <- diag(3)
    res$mean_err <- NA_real_
    res$n_pairs_final <- 0L
    res$scale_total <- 1.0
    res$rotation_deg <- 0.0
  }

  pairs <- unique_match_pairs(xy_cut, hand_xy, dmax = match_dmax)
  keep_src <- rep(FALSE, nrow(xy_cut))
  if (nrow(pairs) > 0) keep_src[pairs$i_src] <- TRUE
  if (!any(keep_src) && nrow(xy_cut) > 0) {
    warn_flags <- c(warn_flags, "warn_unique_match_zero_keep_bbox")
    keep_src <- rep(TRUE, nrow(xy_cut))
    pairs <- tibble(i_src = integer(), i_hand = integer(), dist = numeric())
  }

  df_out <- df_cut[keep_src, , drop = FALSE]
  xy_out <- xy_cut[keep_src, , drop = FALSE]
  d <- if (nrow(pairs) > 0) pairs$dist else numeric(0)
  warn_clean <- warn_flags[!is.na(warn_flags) & nzchar(warn_flags)]

  align_quality <- tibble(
    n_src_raw = nrow(src_xy),
    n_fit_src = nrow(src_fit),
    n_fit_ref = nrow(ref_fit),
    n_clip_bbox_after_align = nrow(xy_cut),
    n_unique_matches = sum(is.finite(d)),
    median_match_dist = if (sum(is.finite(d)) == 0) NA_real_ else stats::median(d[is.finite(d)], na.rm = TRUE),
    p90_match_dist = if (sum(is.finite(d)) == 0) NA_real_ else stats::quantile(d[is.finite(d)], probs = 0.90, na.rm = TRUE, names = FALSE),
    mean_err = res$mean_err,
    n_pairs_iter_final = res$n_pairs_final,
    scale_total = res$scale_total,
    rotation_deg = res$rotation_deg,
    warn_scale_clamped = isTRUE(res$warn_scale_clamped),
    fallback_rigid = isTRUE(res$fallback_rigid),
    warn_flags = paste(unique(warn_clean), collapse = "|")
  )

  list(df = df_out, xy = xy_out, T = T, align_quality = align_quality, match_pairs = pairs)
}

#### 08 TODO SECTIONS (PLACEHOLDERS) ####

relative_position_counter <- function(df, center_e = NULL, center_n = NULL) {
  message("NOTE: relative_position_counter() not implemented yet; returning input unchanged.")
  df
}

write_formatted_map_output <- function(hand_sf, tls_sf, mls_sf, out_path, example_map_path = NULL) {
  message("NOTE: write_formatted_map_output() not implemented yet; skipping formatted map output.")
  invisible(out_path)
}

to_removed <- function(x) {
  x <- trimws(tolower(as.character(x)))
  dplyr::case_when(
    is.na(x) ~ 0L,
    x %in% c("", "na", "nan") ~ 0L,
    x %in% c("true","t","1","yes","y","remove","removed","cut","thin","x") ~ 1L,
    x %in% c("false","f","0","no","n","keep","retained") ~ 0L,
    TRUE ~ 0L
  )
}

index_output_overlap_stub <- function(augmented_csv_path, ID_COL = "tree_id", THIN_PREFIX = "thin_") {
  dat <- read_any_csv(augmented_csv_path, guess_max = guess_max, name_repair = name_repair)
  if (!ID_COL %in% names(dat)) stop("ID_COL not found: ", ID_COL)

  thin_cols <- names(dat)[startsWith(names(dat), THIN_PREFIX)]
  if (length(thin_cols) == 0) stop("No thin_* columns found with prefix: ", THIN_PREFIX)

  thin_long <- dat %>%
    select(all_of(ID_COL), all_of(thin_cols)) %>%
    pivot_longer(cols = all_of(thin_cols), names_to = "index", values_to = "thin_val") %>%
    mutate(
      index = sub(paste0("^", THIN_PREFIX), "", .data$index),
      removed = to_removed(.data$thin_val)
    )

  removed_counts <- thin_long %>%
    group_by(.data$index) %>%
    summarise(
      n_removed = sum(.data$removed),
      n_total = n(),
      share_removed = n_removed / n_total,
      .groups = "drop"
    ) %>%
    arrange(desc(.data$n_removed))

  removed_sets <- thin_long %>%
    filter(.data$removed == 1) %>%
    group_by(.data$index) %>%
    summarise(removed_ids = list(unique(.data[[ID_COL]])), .groups = "drop")

  get_removed <- function(idx) {
    x <- removed_sets %>% filter(.data$index == idx) %>% pull(.data$removed_ids)
    if (length(x) == 0) character(0) else as.character(x[[1]])
  }

  compare_two <- function(idx_a, idx_b) {
    A <- get_removed(idx_a)
    B <- get_removed(idx_b)
    inter <- intersect(A, B)
    uni <- union(A, B)
    tibble(
      idx_a = idx_a,
      idx_b = idx_b,
      removed_A = length(A),
      removed_B = length(B),
      both_removed = length(inter),
      jaccard = ifelse(length(uni) == 0, NA_real_, length(inter) / length(uni))
    )
  }

  pairwise <- tidyr::expand_grid(a = unique(thin_long$index), b = unique(thin_long$index)) %>%
    filter(.data$a < .data$b) %>%
    mutate(res = purrr::map2(.data$a, .data$b, compare_two)) %>%
    tidyr::unnest(.data$res) %>%
    arrange(desc(.data$jaccard))

  per_tree_agreement <- thin_long %>%
    group_by(.data[[ID_COL]]) %>%
    summarise(
      n_indices_removed = sum(.data$removed),
      removed_by = paste(.data$index[.data$removed == 1], collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(desc(.data$n_indices_removed))

  list(
    removed_counts = removed_counts,
    pairwise = pairwise,
    per_tree_agreement = per_tree_agreement
  )
}

#### 08B VISUALIZATION HELPERS ####

vis_as_bool01 <- function(x) {
  if (is.logical(x)) return(as.integer(ifelse(is.na(x), 0L, x)))
  if (is.numeric(x)) return(as.integer(ifelse(is.na(x), 0L, x != 0)))
  z <- trimws(tolower(as.character(x)))
  z[z %in% c("", "na", "nan", "null")] <- NA_character_
  out <- rep(0L, length(z))
  out[z %in% c("1","true","t","yes","y","mls","tls","hand","present")] <- 1L
  out[is.na(z)] <- 0L
  out
}

vis_is_remove <- function(x) {
  z <- trimws(tolower(as.character(x)))
  z[z %in% c("", "na", "nan", "null")] <- NA_character_
  as.integer(!is.na(z) & z %in% vis_thin_remove_values)
}

vis_add_combo <- function(df) {
  df %>%
    mutate(
      Combination = dplyr::case_when(
        MLS == 1L & TLS == 1L & Hand == 1L ~ "MLS + TLS + Hand",
        MLS == 1L & TLS == 1L & Hand == 0L ~ "MLS + TLS",
        MLS == 1L & TLS == 0L & Hand == 1L ~ "MLS + Hand",
        MLS == 0L & TLS == 1L & Hand == 1L ~ "TLS + Hand",
        MLS == 1L & TLS == 0L & Hand == 0L ~ "only MLS",
        MLS == 0L & TLS == 1L & Hand == 0L ~ "only TLS",
        MLS == 0L & TLS == 0L & Hand == 1L ~ "only Hand",
        TRUE ~ "none"
      )
    ) %>%
    filter(.data$Combination != "none")
}

build_overlap_points <- function(hand_df, tls_df, mls_df, dmax = 8) {
  hx <- find_col_ci(hand_df, c("E_base", "E_utm32", "E_aligned", "E", "X", "x"))
  hy <- find_col_ci(hand_df, c("N_base", "N_utm32", "N_aligned", "N", "Y", "y"))
  tx <- find_col_ci(tls_df,  c("E_aligned", "E_utm32", "E", "X", "x"))
  ty <- find_col_ci(tls_df,  c("N_aligned", "N_utm32", "N", "Y", "y"))
  mx <- find_col_ci(mls_df,  c("E_aligned", "E_utm32", "E", "X", "x"))
  my <- find_col_ci(mls_df,  c("N_aligned", "N_utm32", "N", "Y", "y"))
  if (is.na(hx) || is.na(hy)) stop("Could not find Hand coordinate columns.")

  hand_xy <- cbind(to_num(hand_df[[hx]]), to_num(hand_df[[hy]]))
  keep_hand <- is.finite(hand_xy[,1]) & is.finite(hand_xy[,2])
  hand_xy <- hand_xy[keep_hand, , drop = FALSE]

  tls_xy <- if (!is.na(tx) && !is.na(ty) && nrow(tls_df) > 0) cbind(to_num(tls_df[[tx]]), to_num(tls_df[[ty]])) else matrix(numeric(0), ncol = 2)
  mls_xy <- if (!is.na(mx) && !is.na(my) && nrow(mls_df) > 0) cbind(to_num(mls_df[[mx]]), to_num(mls_df[[my]])) else matrix(numeric(0), ncol = 2)
  if (nrow(tls_xy) > 0) tls_xy <- tls_xy[is.finite(tls_xy[,1]) & is.finite(tls_xy[,2]), , drop = FALSE]
  if (nrow(mls_xy) > 0) mls_xy <- mls_xy[is.finite(mls_xy[,1]) & is.finite(mls_xy[,2]), , drop = FALSE]

  if (nrow(hand_xy) == 0) stop("No finite Hand points available for overlap calculation.")

  dist_tls <- rep(Inf, nrow(hand_xy))
  dist_mls <- rep(Inf, nrow(hand_xy))
  if (nrow(tls_xy) > 0) dist_tls <- RANN::nn2(data = tls_xy, query = hand_xy, k = 1)$nn.dists[,1]
  if (nrow(mls_xy) > 0) dist_mls <- RANN::nn2(data = mls_xy, query = hand_xy, k = 1)$nn.dists[,1]

  tibble::tibble(
    point_id = seq_len(nrow(hand_xy)),
    E_utm32 = hand_xy[,1],
    N_utm32 = hand_xy[,2],
    MLS = as.integer(is.finite(dist_mls) & dist_mls <= dmax),
    TLS = as.integer(is.finite(dist_tls) & dist_tls <= dmax),
    Hand = 1L,
    dist_mls = dist_mls,
    dist_tls = dist_tls
  )
}

write_mls_support_layers <- function(mls_df, all_mls_csv, traj_csv) {
  mx <- find_col_ci(mls_df, c("E_utm32", "E_aligned", "E", "X", "x"))
  my <- find_col_ci(mls_df, c("N_utm32", "N_aligned", "N", "Y", "y"))
  if (is.na(mx) || is.na(my)) return(invisible(NULL))

  out <- mls_df %>%
    mutate(E_utm32 = to_num(.data[[mx]]), N_utm32 = to_num(.data[[my]])) %>%
    filter(is.finite(.data$E_utm32), is.finite(.data$N_utm32))
  if (nrow(out) == 0) return(invisible(NULL))

  write_csv_safe(out %>% select(E_utm32, N_utm32), all_mls_csv)
  write_csv_safe(out %>% mutate(seq_id = row_number()) %>% select(seq_id, E_utm32, N_utm32), traj_csv)
  invisible(NULL)
}

write_overlap_map <- function(in_points_csv,
                              out_png,
                              in_all_mls_csv = NULL,
                              in_mls_traj_csv = NULL,
                              dpi_out = 300,
                              width_in = 12,
                              height_in = 8,
                              preview_plot = FALSE,
                              interactive_preview = FALSE,
                              crs_epsg = 32632,
                              thin_rule_col = NULL,
                              map_title  = vis_map_title,
                              author_text = vis_author_text,
                              date_text   = vis_date_text,
                              crs_text    = vis_crs_text) {
  if (!file.exists(in_points_csv)) {
    warning("Map skipped: overlap CSV not found: ", in_points_csv)
    return(invisible(NULL))
  }

  coord_candidates_x <- c("E_aligned", "E_utm32", "E", "X", "x", "easting", "east")
  coord_candidates_y <- c("N_aligned", "N_utm32", "N", "Y", "y", "northing", "north")
  mls_flag_candidates  <- c("MLS", "mls", "has_mls", "in_mls", "mls_present")
  tls_flag_candidates  <- c("TLS", "tls", "has_tls", "in_tls", "tls_present")
  hand_flag_candidates <- c("Hand", "hand", "has_hand", "in_hand", "hand_present")

  COL_ONLY_MLS           <- "#405DA6"
  COL_ONLY_TLS           <- "#FE7F00"
  COL_ONLY_HAND          <- "#339E2C"
  COL_MLS_TLS            <- "#6A3D9A"
  COL_TLS_HAND           <- "#FEDA31"
  COL_MLS_HAND           <- "#2087B9"
  COL_MLS_TLS_HAND       <- "#4A5C4F"
  COL_ALL_MLS_REGISTERED <- "#C5ABCB"
  COL_MLS_TRAJECTORY     <- "#D4D4D4"

  pts <- read_any_csv(in_points_csv, guess_max = guess_max, name_repair = name_repair)
  x_col <- find_col_ci(pts, coord_candidates_x)
  y_col <- find_col_ci(pts, coord_candidates_y)
  mls_col <- find_col_ci(pts, mls_flag_candidates)
  tls_col <- find_col_ci(pts, tls_flag_candidates)
  hand_col <- find_col_ci(pts, hand_flag_candidates)
  if (is.na(x_col) || is.na(y_col) || is.na(mls_col) || is.na(tls_col) || is.na(hand_col)) {
    warning("Map skipped: could not detect required coord/flag columns in overlap CSV.")
    return(invisible(NULL))
  }

  pts <- pts %>%
    mutate(
      x = to_num(.data[[x_col]]),
      y = to_num(.data[[y_col]]),
      MLS  = vis_as_bool01(.data[[mls_col]]),
      TLS  = vis_as_bool01(.data[[tls_col]]),
      Hand = vis_as_bool01(.data[[hand_col]])
    ) %>%
    filter(is.finite(.data$x), is.finite(.data$y))

  if (!is.null(thin_rule_col) && thin_rule_col %in% names(pts)) {
    pts <- pts %>%
      mutate(.remove = vis_is_remove(.data[[thin_rule_col]])) %>%
      filter(.data$.remove == 1L) %>%
      select(-.data$.remove)
  }

  pts <- vis_add_combo(pts)
  if (nrow(pts) == 0) {
    warning("Map skipped: no overlap points left after filtering.")
    return(invisible(NULL))
  }

  pts$cat <- factor(
    pts$Combination,
    levels = c("only MLS", "only TLS", "only Hand", "MLS + TLS", "TLS + Hand", "MLS + Hand", "MLS + TLS + Hand")
  )
  pts_sf <- st_as_sf(pts, coords = c("x", "y"), crs = crs_epsg, remove = FALSE)

  all_mls_sf <- NULL
  if (!is.null(in_all_mls_csv) && nzchar(in_all_mls_csv) && file.exists(in_all_mls_csv)) {
    all_mls <- read_any_csv(in_all_mls_csv, guess_max = guess_max, name_repair = name_repair)
    ax <- find_col_ci(all_mls, coord_candidates_x)
    ay <- find_col_ci(all_mls, coord_candidates_y)
    if (!is.na(ax) && !is.na(ay)) {
      all_mls <- all_mls %>% mutate(x = to_num(.data[[ax]]), y = to_num(.data[[ay]])) %>% filter(is.finite(.data$x), is.finite(.data$y))
      if (nrow(all_mls) > 0) all_mls_sf <- st_as_sf(all_mls, coords = c("x", "y"), crs = crs_epsg, remove = FALSE)
    }
  }

  traj_sf <- NULL
  if (!is.null(in_mls_traj_csv) && nzchar(in_mls_traj_csv) && file.exists(in_mls_traj_csv)) {
    tr <- read_any_csv(in_mls_traj_csv, guess_max = guess_max, name_repair = name_repair)
    tx <- find_col_ci(tr, coord_candidates_x)
    ty <- find_col_ci(tr, coord_candidates_y)
    if (!is.na(tx) && !is.na(ty)) {
      tr <- tr %>% mutate(x = to_num(.data[[tx]]), y = to_num(.data[[ty]])) %>% filter(is.finite(.data$x), is.finite(.data$y))
      if ("seq_id" %in% names(tr)) tr <- tr %>% arrange(to_num(.data$seq_id))
      if (nrow(tr) >= 2) {
        traj_sf <- st_sf(geometry = st_sfc(st_linestring(as.matrix(tr[, c("x", "y")])), crs = crs_epsg))
      }
    }
  }

  col_map <- c(
    "only MLS" = COL_ONLY_MLS,
    "only TLS" = COL_ONLY_TLS,
    "only Hand" = COL_ONLY_HAND,
    "MLS + TLS" = COL_MLS_TLS,
    "TLS + Hand" = COL_TLS_HAND,
    "MLS + Hand" = COL_MLS_HAND,
    "MLS + TLS + Hand" = COL_MLS_TLS_HAND
  )

  p <- ggplot() +
    {if (!is.null(all_mls_sf)) geom_sf(data = all_mls_sf, size = 2.3, alpha = 0.9, color = COL_ALL_MLS_REGISTERED)} +
    {if (!is.null(traj_sf)) geom_sf(data = traj_sf, linewidth = 1.2, alpha = 1, color = COL_MLS_TRAJECTORY)} +
    geom_point(data = pts, aes(x = .data$x, y = .data$y, color = .data$cat), size = 3.2) +
    scale_color_manual(values = col_map, name = NULL) +
    coord_sf(crs = sf::st_crs(crs_epsg), datum = NA) +
    labs(title = map_title, x = NULL, y = NULL) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 22, face = "bold", hjust = 0),
      legend.position = c(0.86, 0.74),
      legend.justification = c(0, 1),
      legend.box.background = element_rect(fill = "white", color = "black", linewidth = 0.6),
      legend.background = element_rect(fill = "white", color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      axis.text = element_text(color = "black")
    )

  if (requireNamespace("ggspatial", quietly = TRUE)) {
    p <- p +
      ggspatial::annotation_north_arrow(
        location = "tr",
        which_north = "true",
        pad_x = grid::unit(0.8, "cm"),
        pad_y = grid::unit(0.6, "cm"),
        style = ggspatial::north_arrow_fancy_orienteering
      ) +
      ggspatial::annotation_scale(
        location = "bl",
        width_hint = 0.22,
        pad_x = grid::unit(0.8, "cm"),
        pad_y = grid::unit(0.8, "cm"),
        text_cex = 0.9
      )
  }

  p <- p +
    annotate(
      "label",
      x = Inf, y = -Inf,
      label = paste(author_text, date_text, crs_text, sep = "\n"),
      hjust = 1, vjust = 0,
      size = 4.1,
      label.r = grid::unit(0, "pt"),
      fill = "white", color = "black",
      label.padding = grid::unit(0.6, "lines")
    )

  if (preview_plot) print(p)
  dir_create(dirname(out_png))
  ggsave(out_png, plot = p, width = width_in, height = height_in, dpi = dpi_out)

  if (interactive_preview) {
    if (!requireNamespace("mapview", quietly = TRUE)) {
      warning("interactive_preview=TRUE requires package 'mapview'.")
    } else {
      mv_layers <- list()
      if (!is.null(traj_sf)) mv_layers <- c(mv_layers, list(mapview::mapview(traj_sf, color = COL_MLS_TRAJECTORY, layer.name = "MLS trajectory")))
      if (!is.null(all_mls_sf)) mv_layers <- c(mv_layers, list(mapview::mapview(all_mls_sf, color = COL_ALL_MLS_REGISTERED, layer.name = "All MLS registered")))
      mv_layers <- c(mv_layers, list(mapview::mapview(pts_sf, zcol = "cat", col.regions = unname(col_map), layer.name = "Overlap points")))
      print(Reduce(`+`, mv_layers))
    }
  }
  invisible(out_png)
}

write_overlap_barplot <- function(in_points_csv, out_png, width_in = 13, height_in = 6.5, dpi_out = 300, thin_rule_col = NULL, bar_stand_id = stand_id) {
  if (!file.exists(in_points_csv)) {
    warning("Bar plot skipped: overlap CSV not found: ", in_points_csv)
    return(invisible(NULL))
  }

  dat <- read_any_csv(in_points_csv, guess_max = guess_max, name_repair = name_repair)
  mls_col <- find_col_ci(dat, c("MLS", "mls", "has_mls", "in_mls", "mls_present"))
  tls_col <- find_col_ci(dat, c("TLS", "tls", "has_tls", "in_tls", "tls_present"))
  hand_col <- find_col_ci(dat, c("Hand", "hand", "has_hand", "in_hand", "hand_present"))
  if (is.na(mls_col) || is.na(tls_col) || is.na(hand_col)) {
    warning("Bar plot skipped: could not detect MLS/TLS/Hand columns.")
    return(invisible(NULL))
  }

  dat <- dat %>%
    mutate(
      MLS  = vis_as_bool01(.data[[mls_col]]),
      TLS  = vis_as_bool01(.data[[tls_col]]),
      Hand = vis_as_bool01(.data[[hand_col]])
    )

  if (!is.null(thin_rule_col) && thin_rule_col %in% names(dat)) {
    dat <- dat %>%
      mutate(.remove = vis_is_remove(.data[[thin_rule_col]])) %>%
      filter(.data$.remove == 1L) %>%
      select(-.data$.remove)
  }

  dat <- vis_add_combo(dat)
  if (nrow(dat) == 0) {
    warning("Bar plot skipped: no overlap rows left after filtering.")
    return(invisible(NULL))
  }

  dat_long <- dat %>%
    select(MLS, TLS, Hand, Combination) %>%
    pivot_longer(cols = c("MLS", "TLS", "Hand"), names_to = "Method", values_to = "present") %>%
    filter(.data$present == 1L)
  if (nrow(dat_long) == 0) {
    warning("Bar plot skipped: no present-flag rows in overlap data.")
    return(invisible(NULL))
  }

  combo_levels_stack <- c("TLS + Hand", "only Hand", "only MLS", "only TLS", "MLS + TLS + Hand", "MLS + TLS", "MLS + Hand")
  method_levels <- c("Hand", "MLS", "TLS")
  color_map <- c(
    "only MLS" = "#405DA6",
    "only TLS" = "#FE7F00",
    "only Hand" = "#339E2C",
    "MLS + TLS" = "#6A3D9A",
    "TLS + Hand" = "#FEDA31",
    "MLS + Hand" = "#2087B9",
    "MLS + TLS + Hand" = "#4A5C4F"
  )

  summ <- dat_long %>%
    count(Method, Combination, name = "Count") %>%
    group_by(Method) %>%
    mutate(Percent = 100 * .data$Count / sum(.data$Count)) %>%
    ungroup() %>%
    mutate(
      Method = factor(.data$Method, levels = method_levels),
      Combination = factor(.data$Combination, levels = combo_levels_stack)
    )

  p <- ggplot(summ, aes(x = .data$Method, y = .data$Percent, fill = .data$Combination)) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = paste0(sprintf("%.1f", .data$Percent), "% (n=", .data$Count, ")")),
      position = position_stack(vjust = 0.5),
      size = 4.2,
      fontface = "bold"
    ) +
    scale_fill_manual(values = color_map, breaks = c("MLS + Hand","MLS + TLS","MLS + TLS + Hand","only Hand","only MLS","only TLS","TLS + Hand")) +
    scale_y_continuous(
      limits = c(0, 100),
      labels = scales::percent_format(scale = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    labs(
      title = paste0("Layered combination proportions for stand ", bar_stand_id),
      x = "Method",
      y = "Portion of the method (%)",
      fill = "Combination"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(face = "bold")
    )

  dir_create(dirname(out_png))
  ggsave(out_png, p, width = width_in, height = height_in, dpi = dpi_out)
  invisible(out_png)
}

#### 08C PROXIMITY-MERGED THINNING LAYER ####

merge_to_bool <- function(x) {
  x0 <- trimws(tolower(as.character(x)))
  out <- rep(FALSE, length(x0))
  out[is.na(x0) | x0 %in% c("", "na", "nan", "null")] <- FALSE
  out[x0 %in% c("1","true","t","yes","y","remove","removed","thin","cut","x")] <- TRUE
  out
}

merge_consensus_thin <- function(df, preferred_cols, k_remove = 3L, primary_col = "thin_CI_Braathe") {
  # Primary rule: if Braathe thinning exists, use it directly (replaces >=k consensus).
  if (!is.null(primary_col) && nzchar(primary_col)) {
    hit <- names(df)[tolower(names(df)) == tolower(primary_col)]
    if (length(hit) > 0) return(merge_to_bool(df[[hit[1]]]))
  }

  # Fallback: legacy >=k consensus across selected thin_* columns.
  cand <- preferred_cols[preferred_cols %in% names(df)]
  if (length(cand) == 0) cand <- names(df)[startsWith(names(df), "thin_")]
  if (length(cand) == 0) return(rep(FALSE, nrow(df)))
  mat <- sapply(cand, function(cc) merge_to_bool(df[[cc]]))
  if (is.vector(mat)) mat <- matrix(mat, ncol = 1)
  rowSums(mat, na.rm = TRUE) >= as.integer(k_remove)
}

merge_xy_from_sf <- function(x) {
  xy <- sf::st_coordinates(x)
  xy <- xy[, 1:2, drop = FALSE]
  colnames(xy) <- c("x", "y")
  xy
}

merge_greedy_unique_pairs <- function(src_xy, trg_xy, max_dist) {
  n_src <- nrow(src_xy); n_trg <- nrow(trg_xy)
  if (n_src < 1 || n_trg < 1) return(rep(NA_integer_, n_src))
  nn <- RANN::nn2(data = trg_xy, query = src_xy, k = 1)
  cand <- tibble(
    i_src = seq_len(n_src),
    i_trg = nn$nn.idx[, 1],
    dist  = nn$nn.dists[, 1]
  ) %>%
    filter(is.finite(.data$dist) & .data$dist <= max_dist) %>%
    arrange(.data$dist)

  used_src <- rep(FALSE, n_src)
  used_trg <- rep(FALSE, n_trg)
  map <- rep(NA_integer_, n_src)
  if (nrow(cand) == 0) return(map)

  for (k in seq_len(nrow(cand))) {
    s <- cand$i_src[k]; t <- cand$i_trg[k]
    if (!used_src[s] && !used_trg[t]) {
      used_src[s] <- TRUE
      used_trg[t] <- TRUE
      map[s] <- t
    }
  }
  map
}

build_proximity_merged_thin_layer <- function(hand_sf, tls_sf, mls_sf,
                                              max_match_m = 1.0,
                                              crs_target_epsg = 32632,
                                              primary_thin_col = "thin_CI_Braathe",
                                              k_remove = 3L,
                                              thin_cols_5 = c("thin_CI_Hegyi", "thin_CI_Braathe", "thin_CI_RK1", "thin_CI_RK2", "thin_CI_RK3"),
                                              hand_thin_candidates = c("Frit", "friT", "frit", "thin_hand", "HandThin")) {
  if (is.null(hand_sf) || is.null(tls_sf) || is.null(mls_sf)) return(NULL)
  if (nrow(hand_sf) == 0 && nrow(tls_sf) == 0 && nrow(mls_sf) == 0) return(NULL)

  gt_ok <- function(x) all(as.character(sf::st_geometry_type(x, by_geometry = TRUE)) %in% c("POINT", "MULTIPOINT"))
  if (!gt_ok(hand_sf) || !gt_ok(tls_sf) || !gt_ok(mls_sf)) stop("Proximity merge requires POINT/MULTIPOINT layers.")

  hand_sf <- ensure_sf_crs(hand_sf, target_epsg = crs_target_epsg, context = "proximity HAND", assume_if_missing = TRUE)
  tls_sf  <- ensure_sf_crs(tls_sf,  target_epsg = crs_target_epsg, context = "proximity TLS",  assume_if_missing = TRUE)
  mls_sf  <- ensure_sf_crs(mls_sf,  target_epsg = crs_target_epsg, context = "proximity MLS",  assume_if_missing = TRUE)

  crs_target <- sf::st_crs(hand_sf)
  if (is.na(crs_target)) stop("hand_sf has no CRS; set CRS before proximity matching.")
  if (sf::st_crs(tls_sf) != crs_target) tls_sf <- sf::st_transform(tls_sf, crs_target)
  if (sf::st_crs(mls_sf) != crs_target) mls_sf <- sf::st_transform(mls_sf, crs_target)

  hand_thin_col <- hand_thin_candidates[hand_thin_candidates %in% names(hand_sf)]
  hand_thin_col <- if (length(hand_thin_col) > 0) hand_thin_col[1] else NA_character_

  hand_sf <- hand_sf %>%
    mutate(
      hand_id = if ("tree_id" %in% names(.)) as.character(.data$tree_id) else as.character(seq_len(dplyr::n())),
      hand_frit01 = if (!is.na(hand_thin_col)) as.integer(merge_to_bool(.data[[hand_thin_col]])) else NA_integer_,
      thin_hand = if (!is.na(hand_thin_col)) merge_to_bool(.data[[hand_thin_col]]) else FALSE
    )
  tls_sf <- tls_sf %>%
    mutate(
      tls_id = if ("tree_id" %in% names(.)) as.character(.data$tree_id) else as.character(seq_len(dplyr::n())),
      thin_tls = merge_consensus_thin(., preferred_cols = thin_cols_5, k_remove = k_remove, primary_col = primary_thin_col)
    )
  mls_sf <- mls_sf %>%
    mutate(
      mls_id = if ("tree_id" %in% names(.)) as.character(.data$tree_id) else as.character(seq_len(dplyr::n())),
      thin_mls = merge_consensus_thin(., preferred_cols = thin_cols_5, k_remove = k_remove, primary_col = primary_thin_col)
    )

  hand_xy <- merge_xy_from_sf(hand_sf)
  tls_xy  <- merge_xy_from_sf(tls_sf)
  mls_xy  <- merge_xy_from_sf(mls_sf)

  hand_to_tls <- merge_greedy_unique_pairs(hand_xy, tls_xy, max_match_m)
  hand_to_mls <- merge_greedy_unique_pairs(hand_xy, mls_xy, max_match_m)

  used_tls <- rep(FALSE, nrow(tls_sf))
  used_mls <- rep(FALSE, nrow(mls_sf))
  used_tls[hand_to_tls[!is.na(hand_to_tls)]] <- TRUE
  used_mls[hand_to_mls[!is.na(hand_to_mls)]] <- TRUE

  dist_hand_tls <- rep(NA_real_, nrow(hand_sf))
  dist_hand_mls <- rep(NA_real_, nrow(hand_sf))
  if (any(!is.na(hand_to_tls))) {
    idx <- which(!is.na(hand_to_tls))
    dist_hand_tls[idx] <- sqrt(rowSums((hand_xy[idx, , drop = FALSE] - tls_xy[hand_to_tls[idx], , drop = FALSE])^2))
  }
  if (any(!is.na(hand_to_mls))) {
    idx <- which(!is.na(hand_to_mls))
    dist_hand_mls[idx] <- sqrt(rowSums((hand_xy[idx, , drop = FALSE] - mls_xy[hand_to_mls[idx], , drop = FALSE])^2))
  }

  tls_rem_idx <- which(!used_tls)
  mls_rem_idx <- which(!used_mls)
  tls_mls_pairs_src <- rep(NA_integer_, length(tls_rem_idx))
  dist_tls_mls_pair <- rep(NA_real_, length(tls_rem_idx))

  if (length(tls_rem_idx) > 0 && length(mls_rem_idx) > 0) {
    map_tls_to_mls_local <- merge_greedy_unique_pairs(
      src_xy = tls_xy[tls_rem_idx, , drop = FALSE],
      trg_xy = mls_xy[mls_rem_idx, , drop = FALSE],
      max_dist = max_match_m
    )
    has_pair <- !is.na(map_tls_to_mls_local)
    if (any(has_pair)) {
      tls_mls_pairs_src[has_pair] <- mls_rem_idx[map_tls_to_mls_local[has_pair]]
      dist_tls_mls_pair[has_pair] <- sqrt(rowSums(
        (tls_xy[tls_rem_idx[has_pair], , drop = FALSE] - mls_xy[tls_mls_pairs_src[has_pair], , drop = FALSE])^2
      ))
      used_mls[tls_mls_pairs_src[has_pair]] <- TRUE
      used_tls[tls_rem_idx[has_pair]] <- TRUE
    }
  }

  hand_rows <- tibble(
    group_id = paste0("H_", hand_sf$hand_id),
    x_out = hand_xy[, 1],
    y_out = hand_xy[, 2],
    hand_id = hand_sf$hand_id,
    hand_frit01 = hand_sf$hand_frit01,
    tls_id  = ifelse(is.na(hand_to_tls), NA_character_, tls_sf$tls_id[hand_to_tls]),
    mls_id  = ifelse(is.na(hand_to_mls), NA_character_, mls_sf$mls_id[hand_to_mls]),
    has_hand = TRUE,
    has_tls  = !is.na(hand_to_tls),
    has_mls  = !is.na(hand_to_mls),
    dist_hand_tls = dist_hand_tls,
    dist_hand_mls = dist_hand_mls,
    dist_tls_mls  = NA_real_,
    thin_hand = hand_sf$thin_hand,
    thin_tls  = ifelse(is.na(hand_to_tls), NA, tls_sf$thin_tls[hand_to_tls]),
    thin_mls  = ifelse(is.na(hand_to_mls), NA, mls_sf$thin_mls[hand_to_mls])
  ) %>%
    mutate(
      match_type = case_when(
        .data$has_hand & .data$has_tls & .data$has_mls ~ "H+T+M",
        .data$has_hand & .data$has_tls & !.data$has_mls ~ "H+T",
        .data$has_hand & !.data$has_tls & .data$has_mls ~ "H+M",
        TRUE ~ "H"
      )
    )

  tls_mls_pairs <- tibble()
  paired_tls_idx <- tls_rem_idx[!is.na(tls_mls_pairs_src)]
  paired_mls_idx <- tls_mls_pairs_src[!is.na(tls_mls_pairs_src)]
  paired_dist <- dist_tls_mls_pair[!is.na(tls_mls_pairs_src)]
  if (length(paired_tls_idx) > 0) {
    mid_xy <- (tls_xy[paired_tls_idx, , drop = FALSE] + mls_xy[paired_mls_idx, , drop = FALSE]) / 2
    tls_mls_pairs <- tibble(
      group_id = paste0("TM_", tls_sf$tls_id[paired_tls_idx], "__", mls_sf$mls_id[paired_mls_idx]),
      x_out = mid_xy[, 1],
      y_out = mid_xy[, 2],
      hand_id = NA_character_,
      hand_frit01 = NA_integer_,
      tls_id  = tls_sf$tls_id[paired_tls_idx],
      mls_id  = mls_sf$mls_id[paired_mls_idx],
      has_hand = FALSE,
      has_tls  = TRUE,
      has_mls  = TRUE,
      dist_hand_tls = NA_real_,
      dist_hand_mls = NA_real_,
      dist_tls_mls  = paired_dist,
      thin_hand = NA,
      thin_tls  = tls_sf$thin_tls[paired_tls_idx],
      thin_mls  = mls_sf$thin_mls[paired_mls_idx],
      match_type = "T+M"
    )
  }

  tls_only_idx <- which(!used_tls)
  mls_only_idx <- which(!used_mls)
  tls_only <- tibble()
  mls_only <- tibble()

  if (length(tls_only_idx) > 0) {
    tls_only <- tibble(
      group_id = paste0("T_", tls_sf$tls_id[tls_only_idx]),
      x_out = tls_xy[tls_only_idx, 1],
      y_out = tls_xy[tls_only_idx, 2],
      hand_id = NA_character_,
      hand_frit01 = NA_integer_,
      tls_id  = tls_sf$tls_id[tls_only_idx],
      mls_id  = NA_character_,
      has_hand = FALSE,
      has_tls  = TRUE,
      has_mls  = FALSE,
      dist_hand_tls = NA_real_,
      dist_hand_mls = NA_real_,
      dist_tls_mls  = NA_real_,
      thin_hand = NA,
      thin_tls  = tls_sf$thin_tls[tls_only_idx],
      thin_mls  = NA,
      match_type = "T"
    )
  }

  if (length(mls_only_idx) > 0) {
    mls_only <- tibble(
      group_id = paste0("M_", mls_sf$mls_id[mls_only_idx]),
      x_out = mls_xy[mls_only_idx, 1],
      y_out = mls_xy[mls_only_idx, 2],
      hand_id = NA_character_,
      hand_frit01 = NA_integer_,
      tls_id  = NA_character_,
      mls_id  = mls_sf$mls_id[mls_only_idx],
      has_hand = FALSE,
      has_tls  = FALSE,
      has_mls  = TRUE,
      dist_hand_tls = NA_real_,
      dist_hand_mls = NA_real_,
      dist_tls_mls  = NA_real_,
      thin_hand = NA,
      thin_tls  = NA,
      thin_mls  = mls_sf$thin_mls[mls_only_idx],
      match_type = "M"
    )
  }

  merged_df <- bind_rows(hand_rows, tls_mls_pairs, tls_only, mls_only) %>%
    mutate(
      population_definition = "AFTER_PROXIMITY_MERGE_THIN_ANY",
      thin_hand = as.logical(.data$thin_hand),
      thin_tls  = as.logical(.data$thin_tls),
      thin_mls  = as.logical(.data$thin_mls),
      thin_votes = rowSums(cbind(
        ifelse(is.na(.data$thin_hand), 0L, as.integer(.data$thin_hand)),
        ifelse(is.na(.data$thin_tls),  0L, as.integer(.data$thin_tls)),
        ifelse(is.na(.data$thin_mls),  0L, as.integer(.data$thin_mls))
      )),
      thin_any = .data$thin_votes >= 1L
    ) %>%
    filter(.data$thin_any)

  if (nrow(merged_df) == 0) return(NULL)
  merged_sf <- sf::st_as_sf(merged_df, coords = c("x_out", "y_out"), crs = crs_target, remove = FALSE)
  ensure_sf_crs(merged_sf, target_epsg = crs_target_epsg, context = "proximity merged layer", assume_if_missing = TRUE)
}

#### 09 RUN PIPELINE ####

dir_create(out_dir)

for (p in c(tls_in, mls_in, hand_in)) {
  if (!file.exists(p)) stop("Input file not found: ", p)
}
if (is.na(sf::st_crs(as.integer(crs_out)))) stop("Invalid crs_out EPSG: ", crs_out)

setting_names_manifest <- c(
  "in_dir", "tls_in", "mls_in", "hand_in", "hand_sheet", "out_dir",
  "name_repair", "guess_max", "verbose", "do_plots", "do_write_gpkg_intermediate",
  "rf_predictors", "rf_num_trees", "rf_mtry", "rf_min_node_size", "rf_importance",
  "dbh_min_valid_m", "mls_use_tls_transfer_model", "mls_transfer_predictors",
  "ci_methods", "compete_radius", "target_radius", "target_source_treecomp",
  "local_radius", "min_neigh", "local_q", "THIN_RULE_COL",
  "E0", "N0", "x_col", "y_col", "crs_out",
  "X_center", "Y_center", "hand_angle_col", "hand_dist_col", "hand_dist_divisor", "crs_hand_raw",
  "icp_max_iter", "icp_dmax_start", "icp_dmax_end", "icp_min_pairs", "icp_trim_frac", "icp_tol",
  "icp_fit_voxel_cell_m", "icp_fit_max_points", "icp_scale_min", "icp_scale_max", "match_dmax_m",
  "between_hand_frit_only", "between_min_n",
  "stand_id", "thin_tag", "tls_tag", "mls_tag", "out_gpkg_base",
  "do_visualizations", "vis_preview_plot", "vis_interactive_preview", "vis_crs_epsg",
  "vis_overlap_csv", "vis_all_mls_csv", "vis_mls_traj_csv", "vis_map_png", "vis_bar_png",
  "vis_dpi_out", "vis_map_width_in", "vis_map_height_in", "vis_bar_width_in", "vis_bar_height_in",
  "vis_thin_rule_col", "vis_thin_remove_values", "vis_map_title", "vis_author_text", "vis_date_text", "vis_crs_text",
  "do_write_proximity_thin_layer", "proximity_max_match_m", "proximity_primary_thin_col", "proximity_k_remove",
  "proximity_thin_cols_5", "proximity_hand_thin_candidates", "proximity_gpkg", "proximity_layer", "proximity_csv"
)

run_id <- stamp()

# ---- TLS ----
tls_raw <- read_any_csv(tls_in, guess_max = guess_max, name_repair = name_repair) %>% standardize_3dfin()
tls_imp <- impute_dbh_rf(
  raw = tls_raw,
  predictors = rf_predictors,
  num_trees = rf_num_trees,
  mtry = rf_mtry,
  min_node_size = rf_min_node_size,
  importance = rf_importance,
  seed = 42,
  dbh_min_valid = dbh_min_valid_m,
  verbose = verbose
)
tls_imputed <- tls_imp$data
tls_transfer_model <- NULL
if (isTRUE(mls_use_tls_transfer_model)) {
  tls_transfer_model <- fit_dbh_transfer_model(
    raw = tls_imputed,
    predictors = mls_transfer_predictors,
    num_trees = rf_num_trees,
    mtry = rf_mtry,
    min_node_size = rf_min_node_size,
    importance = rf_importance,
    seed = 42,
    dbh_min_valid = dbh_min_valid_m,
    verbose = verbose
  )
  if (is.null(tls_transfer_model) && isTRUE(verbose)) {
    message("TLS transfer model unavailable; MLS will use local RF/fallback only.")
  }
}

tls_tc <- compute_treecomp_thinning(
  raw = tls_imputed,
  ci_methods = ci_methods,
  compete_radius = compete_radius,
  target_radius = target_radius,
  target_source = target_source_treecomp,
  local_radius = local_radius,
  local_q = local_q,
  min_neigh = min_neigh,
  do_plots = do_plots,
  plot_dir = file.path(out_dir, paste0("CI_plots_local_rank_", tls_tag)),
  verbose = verbose
)
tls_aug <- tls_tc$data

tls_aug_path <- file.path(out_dir, paste0("Stand_", stand_id, "_", tls_tag, "_", thin_tag, "_augmented.csv"))
write_csv_safe(tls_aug, tls_aug_path)

tls_for_next <- tls_aug
if (!is.null(THIN_RULE_COL) && nzchar(THIN_RULE_COL) && THIN_RULE_COL %in% names(tls_for_next)) {
  tls_for_next <- tls_for_next %>% filter(is.na(.data[[THIN_RULE_COL]]) | .data[[THIN_RULE_COL]] != "remove")
}

tls_utm <- local_to_utm32(tls_for_next, x_col = x_col, y_col = y_col, E0 = E0, N0 = N0)
tls_utm_path <- file.path(out_dir, paste0("Stand_", stand_id, "_", tls_tag, "_", thin_tag, "_UTM32.csv"))
write_csv_safe(tls_utm, tls_utm_path)

if (do_write_gpkg_intermediate) {
  tls_gpkg <- file.path(out_dir, paste0("Stand_", stand_id, "_", tls_tag, "_", thin_tag, "_UTM32.gpkg"))
  tryCatch(
    suppressWarnings(write_gpkg_points(tls_utm, "E_utm32", "N_utm32", tls_gpkg, layer = paste0(tls_tag, "_UTM32"), crs = crs_out, overwrite = TRUE)),
    error = function(e) warning("Intermediate TLS GPKG write failed: ", conditionMessage(e))
  )
}

# ---- MLS ----
mls_raw <- read_any_csv(mls_in, guess_max = guess_max, name_repair = name_repair)
mls_raw <- standardize_3dfin(mls_raw)

mls_imp <- impute_dbh_rf(
  raw = mls_raw,
  predictors = rf_predictors,
  num_trees = rf_num_trees,
  mtry = rf_mtry,
  min_node_size = rf_min_node_size,
  importance = rf_importance,
  seed = 42,
  dbh_min_valid = dbh_min_valid_m,
  external_model = if (!is.null(tls_transfer_model)) tls_transfer_model$model else NULL,
  external_predictors = if (!is.null(tls_transfer_model)) tls_transfer_model$predictors else mls_transfer_predictors,
  external_label = "tls_transfer",
  verbose = verbose
)
mls_imputed <- mls_imp$data

mls_tc <- compute_treecomp_thinning(
  raw = mls_imputed,
  ci_methods = ci_methods,
  compete_radius = compete_radius,
  target_radius = target_radius,
  target_source = target_source_treecomp,
  local_radius = local_radius,
  local_q = local_q,
  min_neigh = min_neigh,
  do_plots = do_plots,
  plot_dir = file.path(out_dir, paste0("CI_plots_local_rank_", mls_tag)),
  verbose = verbose
)
mls_aug <- mls_tc$data

mls_aug_path <- file.path(out_dir, paste0("Stand_", stand_id, "_", mls_tag, "_", thin_tag, "_augmented.csv"))
write_csv_safe(mls_aug, mls_aug_path)

mls_for_next <- mls_aug
if (!is.null(THIN_RULE_COL) && nzchar(THIN_RULE_COL) && THIN_RULE_COL %in% names(mls_for_next)) {
  mls_for_next <- mls_for_next %>% filter(is.na(.data[[THIN_RULE_COL]]) | .data[[THIN_RULE_COL]] != "remove")
}

mls_utm <- local_to_utm32(mls_for_next, x_col = x_col, y_col = y_col, E0 = E0, N0 = N0)
mls_utm_path <- file.path(out_dir, paste0("Stand_", stand_id, "_", mls_tag, "_", thin_tag, "_UTM32.csv"))
write_csv_safe(mls_utm, mls_utm_path)

if (do_write_gpkg_intermediate) {
  mls_gpkg <- file.path(out_dir, paste0("Stand_", stand_id, "_", mls_tag, "_", thin_tag, "_UTM32.gpkg"))
  tryCatch(
    suppressWarnings(write_gpkg_points(mls_utm, "E_utm32", "N_utm32", mls_gpkg, layer = paste0(mls_tag, "_UTM32"), crs = crs_out, overwrite = TRUE)),
    error = function(e) warning("Intermediate MLS GPKG write failed: ", conditionMessage(e))
  )
}

# ---- HAND ----
hand_df <- load_hand(
  hand_in = hand_in,
  hand_sheet = hand_sheet,
  crs_out = crs_out,
  crs_hand_raw = crs_hand_raw,
  E0 = E0, N0 = N0,
  X_center = X_center, Y_center = Y_center,
  angle_col = hand_angle_col,
  dist_col = hand_dist_col,
  dist_divisor = hand_dist_divisor
) %>%
  rename_reserved()
hand_n_loaded <- nrow(hand_df)

hand_path <- file.path(out_dir, paste0("Stand_", stand_id, "_HAND_", thin_tag, "_UTM32.csv"))
write_csv_safe(hand_df, hand_path)

# ---- ALIGNMENT ----
hand_xy <- cbind(to_num(hand_df$E_utm32), to_num(hand_df$N_utm32))
tls_xy  <- cbind(to_num(tls_utm$E_utm32), to_num(tls_utm$N_utm32))
mls_xy  <- cbind(to_num(mls_utm$E_utm32), to_num(mls_utm$N_utm32))

ok_hand <- is.finite(hand_xy[,1]) & is.finite(hand_xy[,2])
ok_tls  <- is.finite(tls_xy[,1])  & is.finite(tls_xy[,2])
ok_mls  <- is.finite(mls_xy[,1])  & is.finite(mls_xy[,2])

hand_df <- hand_df[ok_hand, , drop = FALSE]; hand_xy <- hand_xy[ok_hand, , drop = FALSE]
tls_utm <- tls_utm[ok_tls,  , drop = FALSE]; tls_xy  <- tls_xy[ok_tls,  , drop = FALSE]
mls_utm <- mls_utm[ok_mls,  , drop = FALSE]; mls_xy  <- mls_xy[ok_mls,  , drop = FALSE]

hand_bbox   <- bbox_from_xy(hand_xy)
hand_nn_med <- median_nn(hand_xy)

tls_res <- align_expand_cut(tls_utm, tls_xy, hand_xy, hand_bbox, hand_nn_med)
mls_res <- align_expand_cut(mls_utm, mls_xy, hand_xy, hand_bbox, hand_nn_med)

hand_out <- hand_df
hand_out$E_base <- hand_xy[,1]
hand_out$N_base <- hand_xy[,2]
hand_out <- sanitize_for_gpkg(hand_out)

tls_out <- tls_res$df
tls_out$E_aligned <- tls_res$xy[,1]
tls_out$N_aligned <- tls_res$xy[,2]
tls_out <- sanitize_for_gpkg(tls_out)

mls_out <- mls_res$df
mls_out$E_aligned <- mls_res$xy[,1]
mls_out$N_aligned <- mls_res$xy[,2]
mls_out <- sanitize_for_gpkg(mls_out)

tls_aligned_path <- file.path(out_dir, paste0("Stand_", stand_id, "_", tls_tag, "_", thin_tag, "_aligned.csv"))
mls_aligned_path <- file.path(out_dir, paste0("Stand_", stand_id, "_", mls_tag, "_", thin_tag, "_aligned.csv"))
write_csv_safe(tls_out, tls_aligned_path)
write_csv_safe(mls_out, mls_aligned_path)

out_gpkg <- file.path(out_dir, paste0(out_gpkg_base, "_", run_id, ".gpkg"))
run_manifest_txt <- file.path(out_dir, paste0(out_gpkg_base, "_", run_id, "_settings.txt"))
cleanup_gpkg(out_gpkg)
out_gpkg_tmp <- tempfile(pattern = "aligned_layers_", fileext = ".gpkg")
cleanup_gpkg(out_gpkg_tmp)

hand_sf <- if (nrow(hand_out) > 0) st_as_sf(hand_out, coords = c("E_base", "N_base"), crs = crs_out, remove = FALSE) else NULL
tls_sf  <- if (nrow(tls_out)  > 0) st_as_sf(tls_out,  coords = c("E_aligned", "N_aligned"), crs = crs_out, remove = FALSE) else NULL
mls_sf  <- if (nrow(mls_out)  > 0) st_as_sf(mls_out,  coords = c("E_aligned", "N_aligned"), crs = crs_out, remove = FALSE) else NULL
hand_sf <- ensure_sf_crs(hand_sf, target_epsg = crs_out, context = "HAND_base", assume_if_missing = TRUE)
tls_sf  <- ensure_sf_crs(tls_sf,  target_epsg = crs_out, context = "TLS_aligned", assume_if_missing = TRUE)
mls_sf  <- ensure_sf_crs(mls_sf,  target_epsg = crs_out, context = "MLS_aligned", assume_if_missing = TRUE)

merged_sf <- NULL
if (isTRUE(do_write_proximity_thin_layer)) {
  merged_sf <- tryCatch(
    build_proximity_merged_thin_layer(
      hand_sf = hand_sf,
      tls_sf = tls_sf,
      mls_sf = mls_sf,
      max_match_m = proximity_max_match_m,
      crs_target_epsg = crs_out,
      primary_thin_col = proximity_primary_thin_col,
      k_remove = proximity_k_remove,
      thin_cols_5 = proximity_thin_cols_5,
      hand_thin_candidates = proximity_hand_thin_candidates
    ),
    error = function(e) {
      warning("Proximity merged-thin layer build failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(merged_sf) && nrow(merged_sf) > 0) {
    merged_sf <- ensure_sf_crs(merged_sf, target_epsg = crs_out, context = proximity_layer, assume_if_missing = TRUE)
    write_csv_safe(sf::st_drop_geometry(merged_sf), proximity_csv)
    tryCatch(
      {
        cleanup_gpkg(proximity_gpkg)
        merged_tmp <- tempfile(pattern = "merged_thin_", fileext = ".gpkg")
        cleanup_gpkg(merged_tmp)
        suppressWarnings(st_write(merged_sf, merged_tmp, layer = proximity_layer, driver = "GPKG", quiet = TRUE, delete_dsn = TRUE))
        ok_copy <- file.copy(merged_tmp, proximity_gpkg, overwrite = TRUE)
        cleanup_gpkg(merged_tmp)
        if (!isTRUE(ok_copy)) stop("Could not copy temporary proximity GPKG to output path.")
      },
      error = function(e) warning("Proximity merged-thin GPKG write failed: ", conditionMessage(e))
    )
  } else {
    message("NOTE: Proximity merged-thin layer has 0 rows (no medium voted to thin).")
  }
}

final_gpkg_ok <- TRUE
final_gpkg_err <- NULL
tryCatch(
  {
    if (!is.null(hand_sf) && nrow(hand_sf) > 0) {
      suppressWarnings(st_write(hand_sf, out_gpkg_tmp, layer = "HAND_base",   driver = "GPKG", quiet = TRUE))
    }
    if (!is.null(tls_sf) && nrow(tls_sf) > 0) {
      suppressWarnings(st_write(tls_sf,  out_gpkg_tmp, layer = "TLS_aligned", driver = "GPKG", quiet = TRUE))
    } else {
      message("NOTE: TLS_aligned layer is empty and was skipped in final GPKG.")
    }
    if (!is.null(mls_sf) && nrow(mls_sf) > 0) {
      suppressWarnings(st_write(mls_sf,  out_gpkg_tmp, layer = "MLS_aligned", driver = "GPKG", quiet = TRUE))
    } else {
      message("NOTE: MLS_aligned layer is empty and was skipped in final GPKG.")
    }
    if (!is.null(merged_sf) && nrow(merged_sf) > 0) {
      suppressWarnings(st_write(merged_sf, out_gpkg_tmp, layer = proximity_layer, driver = "GPKG", quiet = TRUE))
    }
    ok_copy <- file.copy(out_gpkg_tmp, out_gpkg, overwrite = TRUE)
    if (!isTRUE(ok_copy)) stop("Could not copy temporary aligned GPKG to output path.")
  },
  error = function(e) {
    final_gpkg_ok <<- FALSE
    final_gpkg_err <<- conditionMessage(e)
  }
)

if (!final_gpkg_ok) {
  cleanup_gpkg(out_gpkg)
  warning("Aligned GPKG write failed: ", final_gpkg_err, ". CSV outputs were written; GPKG was skipped.")
}
cleanup_gpkg(out_gpkg_tmp)

message("DONE")
message("Augmented TLS:  ", normalizePath(tls_aug_path, winslash = "/", mustWork = FALSE))
message("Augmented MLS:  ", normalizePath(mls_aug_path, winslash = "/", mustWork = FALSE))
message("UTM TLS:        ", normalizePath(tls_utm_path, winslash = "/", mustWork = FALSE))
message("UTM MLS:        ", normalizePath(mls_utm_path, winslash = "/", mustWork = FALSE))
message("HAND UTM:       ", normalizePath(hand_path, winslash = "/", mustWork = FALSE))
message("Aligned TLS:    ", normalizePath(tls_aligned_path, winslash = "/", mustWork = FALSE))
message("Aligned MLS:    ", normalizePath(mls_aligned_path, winslash = "/", mustWork = FALSE))
if (final_gpkg_ok) {
  message("Aligned GPKG:   ", normalizePath(out_gpkg, winslash = "/", mustWork = FALSE))
} else {
  message("Aligned GPKG:   skipped (write failed)")
}
if (isTRUE(do_write_proximity_thin_layer)) {
  if (!is.null(merged_sf) && nrow(merged_sf) > 0) {
    message("Merged Thin CSV:", normalizePath(proximity_csv, winslash = "/", mustWork = FALSE))
    message("Merged Thin GPKG:", normalizePath(proximity_gpkg, winslash = "/", mustWork = FALSE), " (layer=", proximity_layer, ")")
  } else {
    message("Merged Thin:    no output points (no medium voted to thin)")
  }
}

if (isTRUE(do_visualizations)) {
  vis_overlap <- build_overlap_points(hand_out, tls_out, mls_out, dmax = match_dmax_m)
  write_csv_safe(vis_overlap, vis_overlap_csv)
  write_mls_support_layers(mls_utm, vis_all_mls_csv, vis_mls_traj_csv)

  tryCatch(
    write_overlap_map(
      in_points_csv = vis_overlap_csv,
      out_png = vis_map_png,
      in_all_mls_csv = vis_all_mls_csv,
      in_mls_traj_csv = vis_mls_traj_csv,
      dpi_out = vis_dpi_out,
      width_in = vis_map_width_in,
      height_in = vis_map_height_in,
      preview_plot = vis_preview_plot,
      interactive_preview = vis_interactive_preview,
      crs_epsg = vis_crs_epsg,
      thin_rule_col = vis_thin_rule_col
    ),
    error = function(e) warning("Map visualization failed: ", conditionMessage(e))
  )

  tryCatch(
    write_overlap_barplot(
      in_points_csv = vis_overlap_csv,
      out_png = vis_bar_png,
      width_in = vis_bar_width_in,
      height_in = vis_bar_height_in,
      dpi_out = vis_dpi_out,
      thin_rule_col = vis_thin_rule_col
    ),
    error = function(e) warning("Barplot visualization failed: ", conditionMessage(e))
  )

  message("Overlap CSV:    ", normalizePath(vis_overlap_csv, winslash = "/", mustWork = FALSE))
  message("Map PNG:        ", normalizePath(vis_map_png, winslash = "/", mustWork = FALSE))
  message("Bar PNG:        ", normalizePath(vis_bar_png, winslash = "/", mustWork = FALSE))
}

manifest_outputs <- list(
  tls_augmented_csv = tls_aug_path,
  mls_augmented_csv = mls_aug_path,
  tls_utm_csv = tls_utm_path,
  mls_utm_csv = mls_utm_path,
  hand_utm_csv = hand_path,
  tls_aligned_csv = tls_aligned_path,
  mls_aligned_csv = mls_aligned_path,
  aligned_gpkg = if (final_gpkg_ok) out_gpkg else "skipped_write_failed",
  proximity_csv = if (isTRUE(do_write_proximity_thin_layer) && !is.null(merged_sf) && nrow(merged_sf) > 0) proximity_csv else "not_written_or_empty",
  proximity_gpkg = if (isTRUE(do_write_proximity_thin_layer) && !is.null(merged_sf) && nrow(merged_sf) > 0) proximity_gpkg else "not_written_or_empty",
  overlap_csv = if (isTRUE(do_visualizations)) vis_overlap_csv else "disabled",
  map_png = if (isTRUE(do_visualizations)) vis_map_png else "disabled",
  bar_png = if (isTRUE(do_visualizations)) vis_bar_png else "disabled"
)

write_run_manifest_txt(
  path = run_manifest_txt,
  run_id = run_id,
  setting_names = setting_names_manifest,
  output_named_list = manifest_outputs
)
message("Run Settings TXT:", normalizePath(run_manifest_txt, winslash = "/", mustWork = FALSE))

#### 10 DIAGNOSTICS / TEST CHUNKS ####
# Writes CSVs into out_dir/diagnostics_<run_id>/
# Uses only already-loaded packages (dplyr/tidyr/purrr/stringr/tibble/rlang/sf + base)

diag_dir <- file.path(out_dir, paste0("diagnostics_", run_id))
dir_create(diag_dir)

# ---------- helpers ----------
add_population_definition <- function(df, pop_def) {
  if (is.null(df)) return(df)
  if (!"population_definition" %in% names(df)) {
    df$population_definition <- pop_def
  } else {
    df$population_definition <- ifelse(is.na(df$population_definition) | !nzchar(df$population_definition), pop_def, df$population_definition)
  }
  df %>% relocate(population_definition, .before = 1)
}

collapse_warn_flags <- function(...) {
  x <- unlist(list(...), use.names = FALSE)
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0) return("")
  paste(unique(x), collapse = "|")
}

thin_bin01 <- function(x) {
  z <- trimws(tolower(as.character(x)))
  z[z %in% c("", "na", "nan", "null", "none")] <- NA_character_
  out <- rep(NA_integer_, length(z))
  out[z %in% c("remove","removed","thin","cut","x","1","true","t","yes","y")] <- 1L
  out[z %in% c("keep","retained","0","false","f","no","n")] <- 0L
  out
}

get_thin_cols_any <- function(df) {
  nms <- names(df)
  nms <- nms[grepl("^thin_", nms, ignore.case = TRUE)]
  # keep only thinning columns derived from indices (avoid "thin_hand" etc if present)
  nms
}

kappa_2cat <- function(a01, b01) {
  ok <- is.finite(a01) & is.finite(b01) & !is.na(a01) & !is.na(b01)
  if (!any(ok)) return(list(kappa = NA_real_, po = NA_real_, pe = NA_real_, n = 0L))
  a <- a01[ok]; b <- b01[ok]
  po <- mean(a == b)
  pa1 <- mean(a == 1); pa0 <- 1 - pa1
  pb1 <- mean(b == 1); pb0 <- 1 - pb1
  pe <- pa1 * pb1 + pa0 * pb0
  kap <- if (isTRUE(all.equal(1, pe))) NA_real_ else (po - pe) / (1 - pe)
  list(kappa = kap, po = po, pe = pe, n = length(a))
}

jaccard_remove <- function(a01, b01) {
  ok <- is.finite(a01) & is.finite(b01) & !is.na(a01) & !is.na(b01)
  if (!any(ok)) return(list(jaccard = NA_real_, inter = 0L, uni = 0L, n = 0L))
  a <- a01[ok]; b <- b01[ok]
  inter <- sum(a == 1 & b == 1)
  uni   <- sum(a == 1 | b == 1)
  jac <- if (uni == 0) NA_real_ else inter / uni
  list(jaccard = jac, inter = inter, uni = uni, n = length(a))
}

fleiss_kappa_binary_complete <- function(n_remove, k) {
  # complete rows only; each item has exactly k ratings, categories {keep, remove}
  if (length(n_remove) < 2) return(NA_real_)
  n_keep <- k - n_remove
  Pi <- (n_keep * (n_keep - 1) + n_remove * (n_remove - 1)) / (k * (k - 1))
  Pbar <- mean(Pi)
  p_remove <- mean(n_remove) / k
  p_keep   <- 1 - p_remove
  Pe <- p_remove^2 + p_keep^2
  if (isTRUE(all.equal(1, Pe))) return(NA_real_)
  (Pbar - Pe) / (1 - Pe)
}

boot_ci <- function(x, fun, B = 2000, seed = 42) {
  x <- x[is.finite(x)]
  if (length(x) < 5) return(c(lo = NA_real_, hi = NA_real_))
  set.seed(seed)
  n <- length(x)
  vals <- replicate(B, fun(sample(x, n, replace = TRUE)))
  stats::quantile(vals, probs = c(0.025, 0.975), na.rm = TRUE, names = FALSE) |> setNames(c("lo","hi"))
}

threshold_sweep <- function(n_remove, k) {
  tibble(
    k_threshold = 0:k,
    n = vapply(0:k, function(t) sum(n_remove >= t), integer(1)),
    pct = vapply(0:k, function(t) 100 * mean(n_remove >= t), numeric(1))
  )
}

expected_binom_baseline <- function(p, k) {
  # independence baseline with common p; for quick sanity comparisons only
  kk <- 0:k
  tibble(
    n_remove = kk,
    p_binom = stats::dbinom(kk, size = k, prob = p),
    pct_binom = 100 * stats::dbinom(kk, size = k, prob = p)
  )
}

# ---------- within-medium index diagnostics ----------
within_medium_index_diagnostics <- function(df, tag, out_dir_diag) {
  thin_cols <- get_thin_cols_any(df)
  thin_cols <- thin_cols[thin_cols %in% names(df)]
  if (length(thin_cols) < 2) {
    return(invisible(NULL))
  }

  mat <- sapply(thin_cols, function(cc) thin_bin01(df[[cc]]))
  if (is.vector(mat)) mat <- matrix(mat, ncol = 1)
  colnames(mat) <- thin_cols
  k <- ncol(mat)

  n_na <- rowSums(is.na(mat))
  complete <- (n_na == 0)
  any_present <- rowSums(!is.na(mat)) > 0
  n_complete <- sum(complete)
  n_any_present <- sum(any_present)

  n_remove_votes <- rowSums(mat == 1, na.rm = TRUE)
  any_remove <- (n_remove_votes >= 1)

  # summary with explicit denominators
  summ <- tibble(
    medium = tag,
    n_indices = k,
    n_rows_total = nrow(mat),
    n_rows_complete = n_complete,
    n_rows_any_decision = n_any_present,
    pct_complete = 100 * mean(complete),
    pct_any_remove_all = 100 * mean(any_remove, na.rm = TRUE),
    pct_any_remove_complete = ifelse(n_complete == 0, NA_real_, 100 * mean(any_remove[complete])),
    pct_all_keep_complete = ifelse(n_complete == 0, NA_real_, 100 * mean(n_remove_votes[complete] == 0)),
    pct_all_remove_complete = ifelse(n_complete == 0, NA_real_, 100 * mean(n_remove_votes[complete] == k)),
    mean_remove_votes_complete = ifelse(n_complete == 0, NA_real_, mean(n_remove_votes[complete])),
    fleiss_kappa_complete = ifelse(n_complete == 0, NA_real_, fleiss_kappa_binary_complete(n_remove_votes[complete], k))
  )
  summ <- add_population_definition(summ, "ALL_ROWS_WITH_WITHIN_THIN_COLUMNS")

  # distribution of remove-votes
  counts_complete <- as.integer(table(factor(n_remove_votes[complete], levels = 0:k)))
  counts_any_present <- as.integer(table(factor(n_remove_votes[any_present], levels = 0:k)))
  counts_complete_any_remove <- as.integer(table(factor(n_remove_votes[complete & n_remove_votes >= 1], levels = 0:k)))

  dist_tbl <- tibble(
    medium = tag,
    n_remove_votes = 0:k,
    count_complete = counts_complete,
    pct_of_complete = ifelse(n_complete == 0, NA_real_, 100 * count_complete / n_complete),
    count_any_decision = counts_any_present,
    pct_of_any_decision = ifelse(n_any_present == 0, NA_real_, 100 * count_any_decision / n_any_present),
    count_complete_any_remove = counts_complete_any_remove,
    pct_of_complete_any_remove = ifelse(sum(counts_complete_any_remove) == 0, NA_real_, 100 * count_complete_any_remove / sum(counts_complete_any_remove))
  ) %>%
    add_population_definition("COMPLETE_CASES_AND_ANY_DECISION")

  # threshold sweep on complete rows
  sweep <- if (n_complete == 0) {
    tibble()
  } else {
    threshold_sweep(n_remove_votes[complete], k) %>%
      mutate(medium = tag) %>%
      select(medium, everything()) %>%
      add_population_definition("COMPLETE_CASES_ONLY")
  }

  # per-index remove rates (complete rows)
  per_idx <- if (n_complete == 0) {
    tibble(
      medium = tag,
      index = thin_cols,
      remove_rate_complete = NA_real_,
      keep_rate_complete = NA_real_
    )
  } else {
    tibble(
      medium = tag,
      index = thin_cols,
      remove_rate_complete = 100 * colMeans(mat[complete, , drop = FALSE] == 1, na.rm = TRUE),
      keep_rate_complete   = 100 * colMeans(mat[complete, , drop = FALSE] == 0, na.rm = TRUE)
    )
  }
  per_idx <- add_population_definition(per_idx, "COMPLETE_CASES_ONLY")

  # pairwise kappa + jaccard on complete rows
  pairs <- combn(thin_cols, 2, simplify = FALSE)
  pw <- bind_rows(lapply(pairs, function(p) {
    a <- p[1]; b <- p[2]
    aa <- mat[, a]; bb <- mat[, b]
    if (n_complete > 0) {
      aa <- aa[complete]; bb <- bb[complete]
    } else {
      aa <- numeric(0); bb <- numeric(0)
    }
    kap <- kappa_2cat(aa, bb)
    jac <- jaccard_remove(aa, bb)
    tibble(
      medium = tag,
      i = a, j = b,
      n = kap$n,
      po = kap$po,
      kappa = kap$kappa,
      jaccard_remove = jac$jaccard,
      inter_remove = jac$inter,
      union_remove = jac$uni
    )
  }))
  pw <- add_population_definition(pw, "COMPLETE_CASES_ONLY")

  # average pairwise agreement (po) and kappa
  pw_summ <- if (nrow(pw) > 0) {
    pw %>%
      summarise(
        medium = first(.data$medium),
        mean_pairwise_po = mean(.data$po, na.rm = TRUE),
        mean_pairwise_kappa = mean(.data$kappa, na.rm = TRUE),
        mean_pairwise_jaccard_remove = mean(.data$jaccard_remove, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    tibble(
      medium = tag,
      mean_pairwise_po = NA_real_,
      mean_pairwise_kappa = NA_real_,
      mean_pairwise_jaccard_remove = NA_real_
    )
  }
  pw_summ <- add_population_definition(pw_summ, "COMPLETE_CASES_ONLY")

  # quick independence baseline using average p_remove (complete rows)
  base_tbl <- tibble()
  if (n_complete > 0) {
    p_bar <- mean(n_remove_votes[complete]) / k
    base_tbl <- expected_binom_baseline(p_bar, k) %>%
      mutate(medium = tag, p_bar = p_bar) %>%
      select(medium, p_bar, everything()) %>%
      add_population_definition("COMPLETE_CASES_ONLY")
  }

  # DBH_source sensitivity (if available): does imputation shift remove rates?
  dbh_sens <- tibble()
  if ("DBH_source" %in% names(df)) {
    src <- trimws(tolower(as.character(df$DBH_source)))
    src[is.na(src) | !nzchar(src)] <- "missing"
    src_class <- dplyr::case_when(
      src %in% c("measured", "field_measured", "meas") ~ "measured",
      grepl("imput|transfer|fallback|model", src) ~ "imputed_transfer_fallback",
      src %in% c("missing", "na", "nan", "") ~ "missing",
      TRUE ~ "other"
    )
    dbh_sens <- bind_rows(lapply(thin_cols, function(cc) {
      v <- thin_bin01(df[[cc]])
      ok <- !is.na(v)
      if (!any(ok)) return(NULL)
      tibble(
        medium = tag,
        index = cc,
        DBH_source_class = src_class[ok],
        removed01 = as.integer(v[ok] == 1L)
      ) %>%
        group_by(.data$medium, .data$index, .data$DBH_source_class) %>%
        summarise(
          n = dplyr::n(),
          remove_rate = 100 * mean(.data$removed01, na.rm = TRUE),
          .groups = "drop"
        )
    }))
    if (nrow(dbh_sens) > 0) {
      dbh_sens <- add_population_definition(dbh_sens, "ROWS_WITH_DECISION_BY_INDEX")
    }
  }

  # bootstrap CI for ">=3 remove" on complete rows (if k>=3)
  ge3_ci <- tibble()
  if (n_complete > 0 && k >= 3) {
    x <- n_remove_votes[complete]
    est <- 100 * mean(x >= 3)
    ci <- boot_ci(x, function(xx) 100 * mean(xx >= 3), B = 2000, seed = 42)
    ge3_ci <- tibble(medium = tag, k = k, pct_ge3_remove_complete = est, ci_lo = ci["lo"], ci_hi = ci["hi"]) %>%
      add_population_definition("COMPLETE_CASES_ONLY")
  }

  redundancy_clusters <- pw %>%
    filter(is.finite(.data$kappa), is.finite(.data$jaccard_remove)) %>%
    filter(.data$kappa > 0.80, .data$jaccard_remove > 0.80) %>%
    mutate(redundant_pair = TRUE)
  redundancy_clusters <- add_population_definition(redundancy_clusters, "COMPLETE_CASES_ONLY")

  # write outputs
  write_csv_safe(summ, file.path(out_dir_diag, paste0("within_", tag, "_summary.csv")))
  write_csv_safe(dist_tbl, file.path(out_dir_diag, paste0("within_", tag, "_removevote_distribution.csv")))
  write_csv_safe(sweep, file.path(out_dir_diag, paste0("within_", tag, "_threshold_sweep.csv")))
  write_csv_safe(per_idx, file.path(out_dir_diag, paste0("within_", tag, "_per_index_rates.csv")))
  write_csv_safe(pw, file.path(out_dir_diag, paste0("within_", tag, "_pairwise_kappa_jaccard.csv")))
  write_csv_safe(pw_summ, file.path(out_dir_diag, paste0("within_", tag, "_pairwise_summary.csv")))
  if (nrow(base_tbl) > 0) write_csv_safe(base_tbl, file.path(out_dir_diag, paste0("within_", tag, "_independence_baseline.csv")))
  if (nrow(dbh_sens) > 0) write_csv_safe(dbh_sens, file.path(out_dir_diag, paste0("within_", tag, "_dbhsource_sensitivity.csv")))
  if (nrow(ge3_ci) > 0) write_csv_safe(ge3_ci, file.path(out_dir_diag, paste0("within_", tag, "_ge3_boot_ci.csv")))
  write_csv_safe(redundancy_clusters, file.path(out_dir_diag, paste0("within_", tag, "_redundancy_clusters.csv")))

  invisible(list(summary = summ, dist = dist_tbl, per_idx = per_idx, pairwise = pw))
}

# ---------- alignment diagnostics ----------
decompose_T <- function(T) {
  z <- decompose_similarity(T)
  tibble(scale = z$scale, rot_rad = z$rot_rad, rot_deg = z$rot_deg, tx = z$tx, ty = z$ty)
}

alignment_diagnostics <- function(tag, src_xy_raw, align_res, hand_xy, hand_bbox, out_dir_diag) {
  if (nrow(src_xy_raw) < 1 || nrow(hand_xy) < 1) return(invisible(NULL))
  if (is.null(align_res) || is.null(align_res$T)) return(invisible(NULL))
  T <- align_res$T

  xy_aligned_all <- apply_T(src_xy_raw, T)
  keep_bbox <- is.finite(xy_aligned_all[,1]) & is.finite(xy_aligned_all[,2]) &
    xy_aligned_all[,1] >= hand_bbox["xmin"] & xy_aligned_all[,1] <= hand_bbox["xmax"] &
    xy_aligned_all[,2] >= hand_bbox["ymin"] & xy_aligned_all[,2] <= hand_bbox["ymax"]

  xy_bbox <- xy_aligned_all[keep_bbox, , drop = FALSE]
  n_raw <- nrow(src_xy_raw)
  n_bbox <- nrow(xy_bbox)

  pairs <- if (n_bbox > 0) unique_match_pairs(xy_bbox, hand_xy, dmax = match_dmax_m) else tibble(i_src = integer(), i_hand = integer(), dist = numeric())
  keep_match <- rep(FALSE, n_bbox)
  if (nrow(pairs) > 0) keep_match[pairs$i_src] <- TRUE
  n_match <- sum(keep_match)

  xy_kept <- if (n_bbox > 0 && any(keep_match)) xy_bbox[keep_match, , drop = FALSE] else matrix(numeric(0), ncol = 2)

  dist_to_hand <- if (nrow(pairs) > 0) pairs$dist else numeric(0)

  Tinfo <- decompose_T(T)
  aq <- align_res$align_quality
  warn_flags <- if (!is.null(aq) && "warn_flags" %in% names(aq)) as.character(aq$warn_flags[1]) else ""
  out <- tibble(
    population_definition = "AFTER_UNIQUE_MATCH_TO_HAND_REFERENCE",
    medium = tag,
    n_src_raw = n_raw,
    n_in_hand_bbox = n_bbox,
    n_unique_matched = n_match,
    pct_in_bbox = 100 * n_bbox / max(1, n_raw),
    pct_matched_of_bbox = ifelse(n_bbox == 0, NA_real_, 100 * n_match / n_bbox),
    pct_matched_of_raw = 100 * n_match / max(1, n_raw),
    dist_median = ifelse(length(dist_to_hand) == 0, NA_real_, stats::median(dist_to_hand)),
    dist_mean = ifelse(length(dist_to_hand) == 0, NA_real_, mean(dist_to_hand)),
    dist_p90 = ifelse(length(dist_to_hand) == 0, NA_real_, stats::quantile(dist_to_hand, 0.90, na.rm = TRUE, names = FALSE)),
    scale_total = Tinfo$scale,
    rotation_deg = Tinfo$rot_deg,
    tx = Tinfo$tx,
    ty = Tinfo$ty,
    mean_err = if (!is.null(aq) && "mean_err" %in% names(aq)) aq$mean_err[1] else NA_real_,
    n_pairs_iter_final = if (!is.null(aq) && "n_pairs_iter_final" %in% names(aq)) aq$n_pairs_iter_final[1] else NA_real_,
    n_clip_bbox_after_align = if (!is.null(aq) && "n_clip_bbox_after_align" %in% names(aq)) aq$n_clip_bbox_after_align[1] else n_bbox,
    warn_scale_clamped = if (!is.null(aq) && "warn_scale_clamped" %in% names(aq)) aq$warn_scale_clamped[1] else NA,
    fallback_rigid = if (!is.null(aq) && "fallback_rigid" %in% names(aq)) aq$fallback_rigid[1] else NA,
    warn_flags = warn_flags
  )

  write_csv_safe(out, file.path(out_dir_diag, paste0("alignment_", tag, "_diagnostics.csv")))

  # density check (median NN distance) for aligned-kept points
  dens <- tibble(
    population_definition = "HAND_REFERENCE_AND_ALIGNED_SUBSETS",
    medium = tag,
    hand_median_nn = median_nn(hand_xy),
    aligned_bbox_median_nn = ifelse(n_bbox < 3, NA_real_, median_nn(xy_bbox)),
    aligned_kept_median_nn = ifelse(nrow(xy_kept) < 3, NA_real_, median_nn(xy_kept))
  )
  write_csv_safe(dens, file.path(out_dir_diag, paste0("alignment_", tag, "_density.csv")))

  invisible(list(alignment = out, density = dens))
}

# ---------- between-media diagnostics (from proximity merged layer) ----------
confusion_stats <- function(ref, pred) {
  ok <- !is.na(ref) & !is.na(pred)
  if (!any(ok)) {
    return(tibble(
      n = 0L, tp = 0L, tn = 0L, fp = 0L, fn = 0L,
      acc = NA_real_, sens = NA_real_, spec = NA_real_, prec = NA_real_, kappa = NA_real_, mcnemar_p = NA_real_,
      warn_flags = "warn_no_valid_pairs"
    ))
  }
  r <- as.integer(ref[ok]); p <- as.integer(pred[ok])
  tp <- sum(r==1 & p==1); tn <- sum(r==0 & p==0); fp <- sum(r==0 & p==1); fn <- sum(r==1 & p==0)
  n <- tp+tn+fp+fn
  if (n < between_min_n) {
    return(tibble(
      n = n, tp = tp, tn = tn, fp = fp, fn = fn,
      acc = NA_real_, sens = NA_real_, spec = NA_real_, prec = NA_real_, kappa = NA_real_, mcnemar_p = NA_real_,
      warn_flags = paste0("warn_low_n_lt_", between_min_n)
    ))
  }
  acc <- (tp+tn)/n
  sens <- ifelse((tp+fn)==0, NA_real_, tp/(tp+fn))
  spec <- ifelse((tn+fp)==0, NA_real_, tn/(tn+fp))
  prec <- ifelse((tp+fp)==0, NA_real_, tp/(tp+fp))
  kap <- kappa_2cat(r, p)$kappa

  b10 <- fn # ref=1 pred=0
  b01 <- fp # ref=0 pred=1
  mcn_p <- if ((b10+b01)==0) NA_real_ else stats::binom.test(min(b10,b01), b10+b01, p=0.5)$p.value

  tibble(
    n = n, tp = tp, tn = tn, fp = fp, fn = fn,
    acc = acc, sens = sens, spec = spec, prec = prec, kappa = kap, mcnemar_p = mcn_p,
    warn_flags = ""
  )
}

between_media_diagnostics <- function(merged_sf_or_df, out_dir_diag) {
  if (is.null(merged_sf_or_df)) return(invisible(NULL))
  df <- if (inherits(merged_sf_or_df, "sf")) sf::st_drop_geometry(merged_sf_or_df) else merged_sf_or_df
  if (!all(c("has_hand","has_tls","has_mls","thin_hand","thin_tls","thin_mls") %in% names(df))) return(invisible(NULL))

  df2 <- df %>%
    mutate(
      match_type = if ("match_type" %in% names(df)) as.character(.data$match_type) else "unknown",
      hand_frit01 = if ("hand_frit01" %in% names(df)) as.integer(.data$hand_frit01) else NA_integer_,
      thin_hand = as.integer(.data$thin_hand),
      thin_tls  = as.integer(.data$thin_tls),
      thin_mls  = as.integer(.data$thin_mls),
      dist_hand_tls = if ("dist_hand_tls" %in% names(df)) to_num(.data$dist_hand_tls) else NA_real_,
      dist_hand_mls = if ("dist_hand_mls" %in% names(df)) to_num(.data$dist_hand_mls) else NA_real_,
      dist_tls_mls  = if ("dist_tls_mls"  %in% names(df)) to_num(.data$dist_tls_mls)  else NA_real_
    )

  # match-type counts + distance summaries on all merged rows
  match_levels <- c("H+T+M", "H+T", "H+M", "H", "T+M", "T", "M")
  mt <- df2 %>%
    count(.data$match_type, name="n") %>%
    right_join(tibble(match_type = match_levels), by = "match_type") %>%
    mutate(n = tidyr::replace_na(.data$n, 0L)) %>%
    mutate(
      pct = 100 * .data$n / sum(.data$n),
      population_definition = "ALL_ROWS_AFTER_PROXIMITY_MERGE_THIN_ANY"
    ) %>%
    relocate(population_definition, .before = 1)

  mt_dist <- df2 %>%
    group_by(.data$match_type) %>%
    summarise(
      n = dplyr::n(),
      dist_hand_tls_median = ifelse(sum(is.finite(.data$dist_hand_tls)) == 0, NA_real_, stats::median(.data$dist_hand_tls[is.finite(.data$dist_hand_tls)], na.rm = TRUE)),
      dist_hand_tls_p90 = ifelse(sum(is.finite(.data$dist_hand_tls)) == 0, NA_real_, stats::quantile(.data$dist_hand_tls[is.finite(.data$dist_hand_tls)], 0.90, na.rm = TRUE, names = FALSE)),
      dist_hand_mls_median = ifelse(sum(is.finite(.data$dist_hand_mls)) == 0, NA_real_, stats::median(.data$dist_hand_mls[is.finite(.data$dist_hand_mls)], na.rm = TRUE)),
      dist_hand_mls_p90 = ifelse(sum(is.finite(.data$dist_hand_mls)) == 0, NA_real_, stats::quantile(.data$dist_hand_mls[is.finite(.data$dist_hand_mls)], 0.90, na.rm = TRUE, names = FALSE)),
      dist_tls_mls_median = ifelse(sum(is.finite(.data$dist_tls_mls)) == 0, NA_real_, stats::median(.data$dist_tls_mls[is.finite(.data$dist_tls_mls)], na.rm = TRUE)),
      dist_tls_mls_p90 = ifelse(sum(is.finite(.data$dist_tls_mls)) == 0, NA_real_, stats::quantile(.data$dist_tls_mls[is.finite(.data$dist_tls_mls)], 0.90, na.rm = TRUE, names = FALSE)),
      .groups = "drop"
    ) %>%
    right_join(tibble(match_type = match_levels), by = "match_type") %>%
    mutate(n = tidyr::replace_na(.data$n, 0L)) %>%
    mutate(population_definition = "ALL_ROWS_AFTER_PROXIMITY_MERGE_THIN_ANY") %>%
    relocate(population_definition, .before = 1)

  write_csv_safe(mt, file.path(out_dir_diag, "between_media_matchtype_counts.csv"))
  write_csv_safe(mt_dist, file.path(out_dir_diag, "between_media_matchtype_distance_summary.csv"))

  # define between-media target population
  pop_def <- "HAND_REFERENCE_POINTS_ONLY"
  ref_df <- df2 %>% filter(.data$has_hand)
  if (isTRUE(between_hand_frit_only)) {
    ref_df <- ref_df %>% filter(.data$hand_frit01 == 1L)
    pop_def <- "HAND_REFERENCE_POINTS_FriT_EQ_1_ONLY"
  }

  pop_meta <- tibble(
    population_definition = pop_def,
    n_population_ref = nrow(ref_df),
    hand_frit_filter = between_hand_frit_only
  )
  write_csv_safe(pop_meta, file.path(out_dir_diag, "between_media_population_definition.csv"))

  stats_HT <- ref_df %>%
    filter(.data$has_tls) %>%
    summarise(out = list(confusion_stats(.data$thin_hand, .data$thin_tls)), .groups = "drop") %>%
    tidyr::unnest(.data$out) %>%
    mutate(
      comparison = "Hand_vs_TLS",
      dist_median = ifelse(sum(is.finite(ref_df$dist_hand_tls[ref_df$has_tls])) == 0, NA_real_, stats::median(ref_df$dist_hand_tls[ref_df$has_tls], na.rm = TRUE)),
      dist_p90 = ifelse(sum(is.finite(ref_df$dist_hand_tls[ref_df$has_tls])) == 0, NA_real_, stats::quantile(ref_df$dist_hand_tls[ref_df$has_tls], 0.90, na.rm = TRUE, names = FALSE))
    )
  stats_HM <- ref_df %>%
    filter(.data$has_mls) %>%
    summarise(out = list(confusion_stats(.data$thin_hand, .data$thin_mls)), .groups = "drop") %>%
    tidyr::unnest(.data$out) %>%
    mutate(
      comparison = "Hand_vs_MLS",
      dist_median = ifelse(sum(is.finite(ref_df$dist_hand_mls[ref_df$has_mls])) == 0, NA_real_, stats::median(ref_df$dist_hand_mls[ref_df$has_mls], na.rm = TRUE)),
      dist_p90 = ifelse(sum(is.finite(ref_df$dist_hand_mls[ref_df$has_mls])) == 0, NA_real_, stats::quantile(ref_df$dist_hand_mls[ref_df$has_mls], 0.90, na.rm = TRUE, names = FALSE))
    )
  stats_TM <- ref_df %>%
    filter(.data$has_tls & .data$has_mls) %>%
    summarise(out = list(confusion_stats(.data$thin_tls, .data$thin_mls)), .groups = "drop") %>%
    tidyr::unnest(.data$out) %>%
    mutate(
      comparison = "TLS_vs_MLS",
      dist_median = ifelse(sum(is.finite(ref_df$dist_tls_mls[ref_df$has_tls & ref_df$has_mls])) == 0, NA_real_, stats::median(ref_df$dist_tls_mls[ref_df$has_tls & ref_df$has_mls], na.rm = TRUE)),
      dist_p90 = ifelse(sum(is.finite(ref_df$dist_tls_mls[ref_df$has_tls & ref_df$has_mls])) == 0, NA_real_, stats::quantile(ref_df$dist_tls_mls[ref_df$has_tls & ref_df$has_mls], 0.90, na.rm = TRUE, names = FALSE))
    )

  out_stats <- bind_rows(stats_HT, stats_HM, stats_TM) %>%
    mutate(population_definition = pop_def) %>%
    relocate(population_definition, .before = 1) %>%
    select(population_definition, comparison, everything())

  write_csv_safe(out_stats, file.path(out_dir_diag, "between_media_consensus_confusion_stats.csv"))

  # distance-binned disagreement (where distances exist) on selected population
  make_bins <- function(x, nb=5) {
    x <- x[is.finite(x)]
    if (length(x) < between_min_n) return(NULL)
    qs <- stats::quantile(x, probs = seq(0,1,length.out=nb+1), na.rm=TRUE)
    unique(qs)
  }

  dist_bins_tbl <- tibble()
  # Hand-TLS
  dHT <- ref_df %>% filter(.data$has_tls & is.finite(.data$dist_hand_tls))
  if (nrow(dHT) >= between_min_n) {
    br <- make_bins(dHT$dist_hand_tls, nb=5)
    if (!is.null(br) && length(br) >= 3) {
      dHT <- dHT %>% mutate(bin = cut(.data$dist_hand_tls, breaks = br, include.lowest = TRUE, right = TRUE))
      dist_bins_tbl <- bind_rows(dist_bins_tbl,
        dHT %>% group_by(.data$bin) %>%
          summarise(
            comparison="Hand_vs_TLS",
            n=dplyr::n(),
            disagree_rate = 100 * mean(.data$thin_hand != .data$thin_tls, na.rm=TRUE),
            dist_med = stats::median(.data$dist_hand_tls, na.rm=TRUE),
            .groups="drop"
          )
      )
    } else {
      dist_bins_tbl <- bind_rows(dist_bins_tbl, tibble(comparison = "Hand_vs_TLS", bin = NA, n = nrow(dHT), disagree_rate = NA_real_, dist_med = NA_real_))
    }
  } else {
    dist_bins_tbl <- bind_rows(dist_bins_tbl, tibble(comparison = "Hand_vs_TLS", bin = NA, n = nrow(dHT), disagree_rate = NA_real_, dist_med = NA_real_))
  }
  # Hand-MLS
  dHM <- ref_df %>% filter(.data$has_mls & is.finite(.data$dist_hand_mls))
  if (nrow(dHM) >= between_min_n) {
    br <- make_bins(dHM$dist_hand_mls, nb=5)
    if (!is.null(br) && length(br) >= 3) {
      dHM <- dHM %>% mutate(bin = cut(.data$dist_hand_mls, breaks = br, include.lowest = TRUE, right = TRUE))
      dist_bins_tbl <- bind_rows(dist_bins_tbl,
        dHM %>% group_by(.data$bin) %>%
          summarise(
            comparison="Hand_vs_MLS",
            n=dplyr::n(),
            disagree_rate = 100 * mean(.data$thin_hand != .data$thin_mls, na.rm=TRUE),
            dist_med = stats::median(.data$dist_hand_mls, na.rm=TRUE),
            .groups="drop"
          )
      )
    } else {
      dist_bins_tbl <- bind_rows(dist_bins_tbl, tibble(comparison = "Hand_vs_MLS", bin = NA, n = nrow(dHM), disagree_rate = NA_real_, dist_med = NA_real_))
    }
  } else {
    dist_bins_tbl <- bind_rows(dist_bins_tbl, tibble(comparison = "Hand_vs_MLS", bin = NA, n = nrow(dHM), disagree_rate = NA_real_, dist_med = NA_real_))
  }
  # TLS-MLS
  dTM <- ref_df %>% filter(.data$has_tls & .data$has_mls & is.finite(.data$dist_tls_mls))
  if (nrow(dTM) >= between_min_n) {
    br <- make_bins(dTM$dist_tls_mls, nb=5)
    if (!is.null(br) && length(br) >= 3) {
      dTM <- dTM %>% mutate(bin = cut(.data$dist_tls_mls, breaks = br, include.lowest = TRUE, right = TRUE))
      dist_bins_tbl <- bind_rows(dist_bins_tbl,
        dTM %>% group_by(.data$bin) %>%
          summarise(
            comparison="TLS_vs_MLS",
            n=dplyr::n(),
            disagree_rate = 100 * mean(.data$thin_tls != .data$thin_mls, na.rm=TRUE),
            dist_med = stats::median(.data$dist_tls_mls, na.rm=TRUE),
            .groups="drop"
          )
      )
    } else {
      dist_bins_tbl <- bind_rows(dist_bins_tbl, tibble(comparison = "TLS_vs_MLS", bin = NA, n = nrow(dTM), disagree_rate = NA_real_, dist_med = NA_real_))
    }
  } else {
    dist_bins_tbl <- bind_rows(dist_bins_tbl, tibble(comparison = "TLS_vs_MLS", bin = NA, n = nrow(dTM), disagree_rate = NA_real_, dist_med = NA_real_))
  }

  dist_bins_tbl <- dist_bins_tbl %>%
    mutate(
      population_definition = pop_def,
      warn_flags = dplyr::case_when(
        .data$n < between_min_n ~ paste0("warn_low_n_lt_", between_min_n),
        is.na(.data$bin) ~ "warn_distance_bins_not_computable",
        TRUE ~ ""
      )
    ) %>%
    relocate(population_definition, .before = 1)
  write_csv_safe(dist_bins_tbl, file.path(out_dir_diag, "between_media_consensus_distance_bins.csv"))

  invisible(list(matchtype = mt, confusion = out_stats, distbins = dist_bins_tbl))
}

# ---------- run diagnostics ----------
diag_ok <- TRUE
diag_err <- NULL

tryCatch({
  # (A) within-medium index agreement (6 indices) for TLS/MLS augmented tables
  within_medium_index_diagnostics(tls_aug, "TLS", diag_dir)
  within_medium_index_diagnostics(mls_aug, "MLS", diag_dir)

  # (B) alignment diagnostics (losses + distances + transform + density)
  # use the finite-filtered UTM coords that the pipeline already created
  align_diag_tls <- alignment_diagnostics("TLS", src_xy_raw = tls_xy, align_res = tls_res, hand_xy = hand_xy, hand_bbox = hand_bbox, out_dir_diag = diag_dir)
  align_diag_mls <- alignment_diagnostics("MLS", src_xy_raw = mls_xy, align_res = mls_res, hand_xy = hand_xy, hand_bbox = hand_bbox, out_dir_diag = diag_dir)

  align_quality_raw <- bind_rows(
    if (!is.null(tls_res$align_quality)) tls_res$align_quality %>% mutate(medium = "TLS") else tibble(),
    if (!is.null(mls_res$align_quality)) mls_res$align_quality %>% mutate(medium = "MLS") else tibble()
  )
  names(align_quality_raw) <- fix_names_case_insensitive_unique(names(align_quality_raw))
  for (nm in c("medium", "scale_total", "rotation_deg", "mean_err", "n_pairs_iter_final",
               "n_clip_bbox_after_align", "n_unique_matches", "median_match_dist",
               "p90_match_dist", "warn_scale_clamped", "fallback_rigid", "warn_flags")) {
    if (!nm %in% names(align_quality_raw)) align_quality_raw[[nm]] <- NA
  }
  align_quality_flags <- align_quality_raw %>%
    transmute(
      population_definition = "AFTER_UNIQUE_MATCH_TO_HAND_REFERENCE",
      medium = .data$medium,
      scale_total = .data$scale_total,
      rotation_deg = .data$rotation_deg,
      mean_err = .data$mean_err,
      n_pairs_iter_final = .data$n_pairs_iter_final,
      n_clip_bbox_after_align = .data$n_clip_bbox_after_align,
      n_unique_matches = .data$n_unique_matches,
      median_match_dist = .data$median_match_dist,
      p90_match_dist = .data$p90_match_dist,
      warn_scale_clamped = .data$warn_scale_clamped,
      fallback_rigid = .data$fallback_rigid,
      WARN_FLAGS = .data$warn_flags
    )
  write_csv_safe(align_quality_flags, file.path(diag_dir, "alignment_quality_flags.csv"))

  # (C) between-media diagnostics on consensus-thin layer (proximity merge)
  if (!is.null(merged_sf) && nrow(merged_sf) > 0) {
    between_media_diagnostics(merged_sf, diag_dir)
  } else if (file.exists(proximity_csv)) {
    between_media_diagnostics(read_any_csv(proximity_csv, guess_max = guess_max, name_repair = name_repair), diag_dir)
  }

  # (D) sanity: raw input row counts & ID uniqueness (helps trace "trees disappearing")
  hand_frit_col <- find_col_ci(hand_df, c("Frit", "friT", "frit", "thin_hand", "HandThin"))
  n_hand_frit1 <- if (!is.na(hand_frit_col)) sum(merge_to_bool(hand_df[[hand_frit_col]]), na.rm = TRUE) else NA_integer_

  raw_counts <- bind_rows(
    tibble(medium = "TLS", population_definition = "ALL_RAW", n_rows = nrow(tls_raw)),
    tibble(medium = "MLS", population_definition = "ALL_RAW", n_rows = nrow(mls_raw)),
    tibble(medium = "HAND", population_definition = "ALL_RAW", n_rows = hand_n_loaded),
    tibble(medium = "TLS", population_definition = "AFTER_FINITE_XY", n_rows = nrow(tls_xy)),
    tibble(medium = "MLS", population_definition = "AFTER_FINITE_XY", n_rows = nrow(mls_xy)),
    tibble(medium = "HAND", population_definition = "AFTER_FINITE_XY", n_rows = nrow(hand_xy)),
    tibble(medium = "TLS", population_definition = "AFTER_BBOX_CLIP", n_rows = if (!is.null(tls_res$align_quality)) tls_res$align_quality$n_clip_bbox_after_align[1] else NA_integer_),
    tibble(medium = "MLS", population_definition = "AFTER_BBOX_CLIP", n_rows = if (!is.null(mls_res$align_quality)) mls_res$align_quality$n_clip_bbox_after_align[1] else NA_integer_),
    tibble(medium = "TLS", population_definition = "AFTER_UNIQUE_MATCH", n_rows = if (!is.null(tls_res$align_quality)) tls_res$align_quality$n_unique_matches[1] else NA_integer_),
    tibble(medium = "MLS", population_definition = "AFTER_UNIQUE_MATCH", n_rows = if (!is.null(mls_res$align_quality)) mls_res$align_quality$n_unique_matches[1] else NA_integer_),
    tibble(medium = "HAND", population_definition = "HAND_FriT_EQ_1_ONLY", n_rows = n_hand_frit1)
  )
  write_csv_safe(raw_counts, file.path(diag_dir, "sanity_input_rowcounts.csv"))

  # ID uniqueness (if present)
  id_sanity_one <- function(df, tag) {
    idc <- find_col_ci(df, c("tree_id","TreeID","ID","trenr","id"))
    if (is.na(idc)) return(tibble(medium=tag, id_col=NA_character_, n= nrow(df), n_unique=NA_integer_, n_dup=NA_integer_))
    ids <- as.character(df[[idc]])
    tibble(
      medium=tag,
      id_col=idc,
      n=nrow(df),
      n_unique=length(unique(ids)),
      n_dup=sum(duplicated(ids))
    )
  }
  id_sanity <- bind_rows(
    id_sanity_one(tls_aug, "TLS_aug"),
    id_sanity_one(mls_aug, "MLS_aug"),
    id_sanity_one(hand_df, "HAND_df"),
    id_sanity_one(tls_out, "TLS_aligned_out"),
    id_sanity_one(mls_out, "MLS_aligned_out")
  ) %>%
    mutate(population_definition = case_when(
      .data$medium %in% c("TLS_aug", "MLS_aug", "HAND_df") ~ "POST_AUGMENTED_TABLE",
      .data$medium %in% c("TLS_aligned_out", "MLS_aligned_out") ~ "AFTER_UNIQUE_MATCH",
      TRUE ~ "UNKNOWN"
    )) %>%
    relocate(population_definition, .before = 1)
  write_csv_safe(id_sanity, file.path(diag_dir, "sanity_id_uniqueness.csv"))

}, error = function(e) {
  diag_ok <<- FALSE
  diag_err <<- conditionMessage(e)
})

if (!diag_ok) {
  warning("Diagnostics failed: ", diag_err)
} else {
  message("Diagnostics written to: ", normalizePath(diag_dir, winslash="/", mustWork=FALSE))
}
#### END DIAGNOSTICS ####
