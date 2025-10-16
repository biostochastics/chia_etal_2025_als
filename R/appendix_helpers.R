# appendix_helpers.R
# Helper functions for computational transparency appendix
#
# NOTE: tar_meta(), tar_manifest(), and tar_network() are called directly
# in Quarto chunks, not through a target, to avoid targets limitations.

#' Format Target Metadata Table
#'
#' @description
#' Formats target metadata into a clean table for display in appendix.
#'
#' @param meta_df Data frame from targets::tar_meta()
#' @param sort_by Character, one of "time", "size", "name" (default: "time")
#' @param top_n Integer, number of targets to show (default: all)
#'
#' @return Formatted tibble ready for knitr::kable()
#'
#' @export
format_target_table <- function(meta_df, sort_by = "time", top_n = NULL) {
  # Format into display table
  result <- meta_df |>
    dplyr::mutate(
      # Convert to numeric (may be character from tar_meta)
      seconds_num = as.numeric(seconds),
      bytes_num = as.numeric(bytes),
      time_min = seconds_num / 60,
      size_mb = bytes_num / 1e6,
      # warnings is a list column - count elements in each list
      warnings_count = purrr::map_int(warnings, function(w) {
        if (is.null(w) || length(w) == 0) 0L else length(w)
      }),
      has_error = !is.na(error)
    ) |>
    dplyr::select(name, time_min, size_mb, warnings_count, has_error, format)

  # Sort based on requested column
  result <- switch(sort_by,
    "time" = dplyr::arrange(result, dplyr::desc(time_min)),
    "size" = dplyr::arrange(result, dplyr::desc(size_mb)),
    "name" = dplyr::arrange(result, name),
    dplyr::arrange(result, dplyr::desc(time_min))
  )

  # Limit to top N if requested
  if (!is.null(top_n)) {
    result <- head(result, top_n)
  }

  result
}

#' Extract Target Command
#'
#' @description
#' Extracts command for a specific target from manifest.
#'
#' @param target_name Character, name of target
#' @param manifest_df Data frame from targets::tar_manifest()
#'
#' @return Character vector with target command
#'
#' @export
extract_target_command <- function(target_name, manifest_df) {
  # Get target command from manifest
  cmd <- manifest_df |>
    dplyr::filter(name == target_name) |>
    dplyr::pull(command)

  if (length(cmd) == 0) {
    return(paste("Target", target_name, "not found in manifest"))
  }

  # Return deparsed command
  as.character(cmd)
}

#' Create Dependency Subgraph
#'
#' @description
#' Creates filtered dependency graph for specific target prefixes.
#'
#' @param network_obj Network object from targets::tar_network()
#' @param prefix Character, target name prefix to filter (e.g., "lcv_")
#'
#' @return Filtered network object for tar_visnetwork()
#'
#' @export
create_subgraph <- function(network_obj, prefix) {
  # Filter nodes
  nodes <- network_obj$vertices |>
    dplyr::filter(stringr::str_starts(name, prefix))

  # Filter edges
  edges <- network_obj$edges |>
    dplyr::filter(
      from %in% nodes$name,
      to %in% nodes$name
    )

  # Return filtered network
  list(
    vertices = nodes,
    edges = edges
  )
}

#' Generate Execution Summary Statistics
#'
#' @description
#' Computes summary statistics for pipeline execution.
#'
#' @param meta_df Data frame from targets::tar_meta()
#'
#' @return Named list with summary statistics
#'
#' @export
summarize_execution <- function(meta_df) {
  # Convert seconds and bytes to numeric (may be character from tar_meta)
  seconds_numeric <- as.numeric(meta_df$seconds)
  bytes_numeric <- as.numeric(meta_df$bytes)

  # Count total warnings across all targets (warnings is a list column)
  total_warnings <- sum(purrr::map_int(meta_df$warnings, function(w) {
    if (is.null(w) || length(w) == 0) 0L else length(w)
  }))

  list(
    total_targets = nrow(meta_df),
    total_time_hours = sum(seconds_numeric, na.rm = TRUE) / 3600,
    total_size_gb = sum(bytes_numeric, na.rm = TRUE) / 1e9,
    avg_time_minutes = mean(seconds_numeric, na.rm = TRUE) / 60,
    median_time_seconds = median(seconds_numeric, na.rm = TRUE),
    max_time_minutes = max(seconds_numeric, na.rm = TRUE) / 60,
    total_warnings = total_warnings,
    total_errors = sum(!is.na(meta_df$error))
  )
}

#' Get Git Information
#'
#' @description
#' Retrieves git commit hash and repository status.
#'
#' @return List with git_commit, git_status, and git_clean flag
#'
#' @export
get_git_info <- function() {
  git_commit <- tryCatch(
    system("git rev-parse HEAD", intern = TRUE),
    error = function(e) "unknown"
  )

  git_status <- tryCatch(
    system("git status --short", intern = TRUE),
    error = function(e) character(0)
  )

  list(
    git_commit = git_commit,
    git_status = git_status,
    git_clean = length(git_status) == 0
  )
}
