#' Build Xenium morphology images aligned to cell centroids
#'
#' Reads a single-channel morphology OME-TIFF pyramid level for each
#' Xenium sample (via reticulate + the Python \code{tifffile} package)
#' and converts it into a raster object aligned to cell centroid
#' coordinates in a monocle3 cell_data_set.
#'
#' Each sample's image is kept as a single \code{\link[grDevices]{as.raster}}
#' object (rendered via \code{\link[ggplot2]{annotation_raster}}) rather than
#' a pixel-per-row data.frame, since Xenium morphology images are typically
#' tens of thousands of pixels per side even at moderate pyramid levels --
#' far too large to explode into one row per pixel. Use \code{level > 0} to
#' load a downsampled pyramid level; this function assumes each pyramid
#' level downsamples by a factor of \code{2^level} relative to full
#' resolution (override with \code{downsample_factor} if a sample's pyramid
#' uses a different scheme).
#'
#' @param cds A monocle3 cell_data_set with \code{x_centroid}/\code{y_centroid}
#'   and \code{sample} columns in colData (see \code{build_xenium_cds()}).
#' @param sample_table A data.frame produced by discover_xenium_sample_table().
#' @param channel Integer morphology channel index, matching the
#'   \code{morphology_focus_XXXX.ome.tif} file name (0 = DAPI by default).
#' @param level Pyramid level to load (0 = full resolution). Even with
#'   raster-based rendering, very low levels can still be too large/slow;
#'   raise this if loading is slow or memory-heavy.
#' @param pixel_size_um Physical size, in microns, of one full-resolution
#'   (\code{level = 0}) pixel. Defaults to the standard Xenium value.
#' @param downsample_factor Optional override for the pyramid downsample
#'   factor at \code{level}. Defaults to \code{2^level}.
#' @param img_buffer Pixel buffer around cells used to compute a suggested
#'   cell-focused crop region (\code{cells_xlim}/\code{cells_ylim}) for
#'   plotting, without discarding any raster data.
#' @param conda_env Optional conda environment name containing
#'   \code{tifffile} (and \code{imagecodecs} for compressed OME-TIFFs).
#'
#' @return A list with elements:
#'   \describe{
#'     \item{images}{A named list, keyed by \code{sample_id}, each containing
#'       \code{raster} (a \code{grDevices::as.raster} object), \code{xlim}/
#'       \code{ylim} (full image pixel extent), \code{cells_xlim}/
#'       \code{cells_ylim} (suggested crop region around cells), and
#'       \code{sample_label}/\code{sample_group}.}
#'     \item{channel}{Channel index used}
#'     \item{level}{Pyramid level used}
#'     \item{pixel_size_um}{Full-resolution pixel size in microns}
#'     \item{downsample_factor}{Downsample factor applied at \code{level}}
#'     \item{effective_pixel_size_um}{Pixel size, in microns, at \code{level}}
#'   }
#'
#' @export
build_xenium_images <- function(
    cds,
    sample_table,
    channel           = 0,
    level             = 1,
    pixel_size_um     = 0.2125,
    downsample_factor = NULL,
    img_buffer        = 5,
    conda_env         = NULL
) {

  required_cols <- c("x_centroid", "y_centroid", "sample")

  if (!all(required_cols %in% colnames(SummarizedExperiment::colData(cds)))) {
    stop(
      "cds must contain cell centroid coordinates in colData: ",
      paste(required_cols, collapse = ", "),
      ". Did you run build_xenium_cds()?"
    )
  }

  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  }

  tifffile <- reticulate::import("tifffile")

  if (is.null(downsample_factor)) {
    downsample_factor <- 2^level
  }
  effective_pixel_size <- pixel_size_um * downsample_factor

  df_cells <- as.data.frame(SummarizedExperiment::colData(cds))
  df_cells$sample <- as.character(df_cells$sample)

  image_list <- list()

  for (i in seq_len(nrow(sample_table))) {

    sample_id <- as.character(sample_table$sample_id[i])
    df_samp   <- df_cells[df_cells$sample == sample_id, , drop = FALSE]

    if (nrow(df_samp) == 0) next

    morph_dir <- sample_table$morphology_dir[i]
    if (is.na(morph_dir) || !dir.exists(morph_dir)) next

    morph_file <- file.path(
      morph_dir,
      sprintf("morphology_focus_%04d.ome.tif", channel)
    )

    if (!file.exists(morph_file)) {
      warning("Morphology file not found, skipping sample ", sample_id,
              ": ", morph_file)
      next
    }

    message(
      "Reading Xenium morphology image: ", morph_file,
      " (channel ", channel, ", level ", level, ")"
    )

    img_py  <- tifffile$imread(morph_file, is_ome = FALSE, level = as.integer(level))
    img_mat <- as.matrix(img_py)

    n_pixels <- length(img_mat)
    if (n_pixels > 2e7) {
      warning(
        "Loaded image for sample ", sample_id, " has ", n_pixels,
        " pixels at level ", level, ". Consider a higher level if this is ",
        "slow or memory-heavy."
      )
    }

    value_range <- range(img_mat, na.rm = TRUE)
    img_norm <- if (diff(value_range) > 0) {
      (img_mat - value_range[1]) / diff(value_range)
    } else {
      matrix(0, nrow = nrow(img_mat), ncol = ncol(img_mat))
    }

    x_rng <- range(df_samp$x_centroid / effective_pixel_size, na.rm = TRUE)
    y_rng <- range(df_samp$y_centroid / effective_pixel_size, na.rm = TRUE)

    image_list[[sample_id]] <- list(
      raster       = grDevices::as.raster(img_norm),
      xlim         = c(0, ncol(img_mat)),
      ylim         = c(0, nrow(img_mat)),
      cells_xlim   = c(x_rng[1] - img_buffer, x_rng[2] + img_buffer),
      cells_ylim   = c(y_rng[1] - img_buffer, y_rng[2] + img_buffer),
      sample_label = sample_table$sample_label[i],
      sample_group = sample_table$sample_group[i]
    )
  }

  list(
    images                  = image_list,
    channel                 = channel,
    level                   = level,
    pixel_size_um           = pixel_size_um,
    downsample_factor       = downsample_factor,
    effective_pixel_size_um = effective_pixel_size
  )
}
