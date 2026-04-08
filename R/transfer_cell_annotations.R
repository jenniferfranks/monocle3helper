# --------------------------------------------------------------
# Transfer cell-level annotations between monocle3 CDS objects
# --------------------------------------------------------------
# Input: reference CDS object with celltype annotations
#        target CDS object
#
# Output: target CDS object with transferred celltype annotations
#
# Usage:
# cds2 <- transfer_cell_annotations(
#   cds_ref = cds1,
#   cds_target = cds2,
#   annotation_col = "cell_type",
#   new_annotation_col = "cell_type",
#   cell_id_col = "cell",   # set to NULL to use rownames
#   preserve_unmatched = TRUE
# )
#' @export

transfer_cell_annotations <- function(
    cds_ref,
    cds_target,
    annotation_col,
    new_annotation_col = NULL,
    cell_id_col = NULL,
    preserve_unmatched = TRUE,
    verbose = TRUE
) {
  # Extract pData
  pd_ref <- as.data.frame(pData(cds_ref), stringsAsFactors = FALSE)
  pd_target <- as.data.frame(pData(cds_target), stringsAsFactors = FALSE)

  # Determine output column name
  if (is.null(new_annotation_col)) {
    new_annotation_col <- annotation_col
  }

  # Check annotation column exists in reference
  if (!annotation_col %in% colnames(pd_ref)) {
    stop("annotation_col not found in cds_ref pData")
  }

  # Determine cell identifiers
  if (!is.null(cell_id_col)) {
    if (!cell_id_col %in% colnames(pd_ref)) {
      stop("cell_id_col not found in cds_ref pData")
    }
    if (!cell_id_col %in% colnames(pd_target)) {
      stop("cell_id_col not found in cds_target pData")
    }

    ref_ids <- as.character(pd_ref[[cell_id_col]])
    target_ids <- as.character(pd_target[[cell_id_col]])
  } else {
    ref_ids <- rownames(pd_ref)
    target_ids <- rownames(pd_target)
  }

  # Basic checks
  if (anyDuplicated(ref_ids) > 0) {
    dup_ids <- unique(ref_ids[duplicated(ref_ids)])
    stop(
      "Duplicate reference cell IDs found. Example duplicates: ",
      paste(head(dup_ids, 10), collapse = ", ")
    )
  }

  if (length(ref_ids) != nrow(pd_ref)) {
    stop("Length of reference IDs does not match number of reference cells")
  }

  if (length(target_ids) != nrow(pd_target)) {
    stop("Length of target IDs does not match number of target cells")
  }

  # Create mapping
  mapping <- setNames(as.character(pd_ref[[annotation_col]]), ref_ids)

  # Match target cells to reference
  match_idx <- match(target_ids, ref_ids)
  matched <- !is.na(match_idx)
  transferred <- rep(NA_character_, length(target_ids))
  transferred[matched] <- mapping[target_ids[matched]]

  # Build output vector
  if (new_annotation_col %in% colnames(pd_target)) {
    existing_values <- as.character(pd_target[[new_annotation_col]])
  } else {
    existing_values <- rep(NA_character_, length(target_ids))
  }

  output_values <- existing_values

  if (preserve_unmatched) {
    output_values[matched] <- transferred[matched]
  } else {
    output_values <- transferred
  }

  # Assign annotation back to CDS
  pData(cds_target)[[new_annotation_col]] <- output_values

  # Diagnostics
  if (verbose) {
    n_total <- length(target_ids)
    n_matched <- sum(matched)
    n_unmatched <- sum(!matched)
    n_ref_na <- sum(is.na(pd_ref[[annotation_col]]))

    message("Annotation transfer summary:")
    message("  Reference annotation: ", annotation_col)
    message("  Target column name:   ", new_annotation_col)
    message("  Total target cells:   ", n_total)
    message("  Matched cells:        ", n_matched)
    message("  Unmatched cells:      ", n_unmatched)
    message("  NA annotations in ref:", n_ref_na)
    message("  preserve_unmatched:   ", preserve_unmatched)

    if (preserve_unmatched) {
      message("  Unmatched target cells retained their previous values")
    } else {
      message("  Unmatched target cells assigned NA")
    }
  }

  return(cds_target)
}

# --------------------------------------------------------------
# Determine majority annotation per cluster
# --------------------------------------------------------------
# Input:  CDS object with transfered celltype annotations
#         must specify which columns contain cluster, celltype label
#
# Output: CDS object with column majority cluster label
#
#Usage:
# cluster_labels <- majority_label_by_cluster(
# cds,
# annotation_col = "cell_type",
# cluster_col = "cluster",
# min_prop = 0.6
# )
#
# # Or directly annotate the CDS
# cluster_labels <- majority_label_by_cluster(
#   cds,
#   annotation_col = "cell_type",
#   annotate_cds = TRUE
# )
#' @export

majority_label_by_cluster <- function(
    cds,
    annotation_col,
    cluster_col = "cluster",
    min_prop = 0,
    annotate_cds = FALSE
) {
  pd <- as.data.frame(pData(cds))

  if (!annotation_col %in% colnames(pd)) {
    stop("annotation_col not found in pData(cds)")
  }
  if (!cluster_col %in% colnames(pd)) {
    stop("cluster_col not found in pData(cds)")
  }

  # Remove NA annotations
  pd <- pd[!is.na(pd[[annotation_col]]), ]

  # Tabulate
  tab <- table(
    cluster = pd[[cluster_col]],
    annotation = pd[[annotation_col]]
  )

  # Compute majority label
  majority_df <- do.call(rbind, lapply(rownames(tab), function(cl) {
    counts <- tab[cl, ]
    prop <- counts / sum(counts)
    top_label <- names(which.max(counts))
    top_prop <- max(prop)

    if (top_prop < min_prop) {
      top_label <- NA
    }

    data.frame(
      cluster = cl,
      majority_label = top_label,
      majority_prop = top_prop,
      stringsAsFactors = FALSE
    )
  }))

  # Optionally annotate CDS
  if (annotate_cds) {
    cluster_to_label <- setNames(
      majority_df$majority_label,
      majority_df$cluster
    )
    pData(cds)$majority_label <- cluster_to_label[
      as.character(pData(cds)[[cluster_col]])
    ]
    # add a return to make sure that if annotate_cds is TRUE, what is returned is cds, not majority_df
    return(cds)
  }

  return(majority_df)
}
