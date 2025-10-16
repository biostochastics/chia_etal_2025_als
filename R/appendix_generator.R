#' Generate Computational Appendix Content
#'
#' @description
#' Creates markdown content for a computational appendix documenting the full pipeline:
#' - Complete code trace for all targets
#' - Function source code
#' - Output summaries
#' - Execution metadata
#'
#' This content is designed to be included in the main Quarto document.
#'
#' @param output_path Path where appendix content will be saved (child document)
#' @param include_outputs Logical; include sample outputs (default TRUE)
#' @param max_output_lines Maximum lines of output to include per target
#'
#' @return Path to generated appendix content file
#' @export
generate_computational_appendix <- function(
    output_path = "reports/_appendix_computational.qmd",
    include_outputs = TRUE,
    max_output_lines = 100
) {

  # Get pipeline manifest
  manifest <- targets::tar_manifest(fields = c("name", "command"))

  # Get execution metadata
  meta <- targets::tar_meta(fields = c("name", "time", "bytes", "seconds"))

  # Start building appendix content (no YAML - this is a child document)
  content <- c(
    "<!-- Auto-generated computational appendix -->",
    "<!-- DO NOT EDIT - regenerate with targets::tar_make() -->",
    "",
    "# Computational Appendix {.appendix}",
    "",
    "This appendix provides complete transparency into the computational pipeline that generated all results in this report.",
    "",
    "## Pipeline Overview",
    "",
    paste0("- **Total computational targets**: ", nrow(manifest)),
    paste0("- **Pipeline execution date**: ", Sys.time()),
    paste0("- **R version**: ", R.version.string),
    paste0("- **`targets` version**: ", packageVersion('targets')),
    "",
    "## Pipeline Dependency Graph {.unnumbered}",
    "",
    "```{r appendix-pipeline-graph}",
    "#| echo: false",
    "#| fig-width: 14",
    "#| fig-height: 12",
    "#| fig-cap: 'Interactive pipeline dependency graph showing all computational targets and their relationships'",
    "library(targets)",
    "tar_visnetwork(",
    "  targets_only = TRUE,",
    "  label = c('time', 'size'),",
    "  level_separation = 150",
    ")",
    "```",
    "",
    "## Execution Summary {.unnumbered}",
    "",
    "```{r appendix-execution-summary}",
    "#| echo: false",
    "library(dplyr)",
    "",
    "meta <- tar_meta(fields = c('name', 'time', 'seconds', 'bytes'))",
    "",
    "meta %>%",
    "  summarize(",
    "    `Total Targets` = n(),",
    "    `Total Time (hours)` = round(sum(seconds, na.rm = TRUE) / 3600, 2),",
    "    `Total Size (GB)` = round(sum(bytes, na.rm = TRUE) / 1e9, 2),",
    "    `Avg Time (minutes)` = round(mean(seconds, na.rm = TRUE) / 60, 2)",
    "  ) %>%",
    "  knitr::kable()",
    "```",
    "",
    "---",
    "",
    "## Complete Code Trace",
    "",
    "The following sections document every computational target in the pipeline, including the exact command executed, the source code of functions used, and output summaries.",
    ""
  )

  # Process each target
  for (i in seq_len(nrow(manifest))) {
    target_name <- manifest$name[i]
    target_cmd <- manifest$command[i]

    # Get metadata if available
    target_meta <- meta[meta$name == target_name, ]

    # Section header
    content <- c(
      content,
      paste0("### Target ", i, ": `", target_name, "` {.unnumbered}"),
      ""
    )

    # Execution metadata
    if (nrow(target_meta) > 0) {
      content <- c(
        content,
        "**Execution Metadata:**",
        "",
        paste0("- Time: ", target_meta$time),
        paste0("- Duration: ", round(target_meta$seconds, 2), " seconds"),
        paste0("- Size: ", format(target_meta$bytes, units = "auto")),
        ""
      )
    }

    # Command
    content <- c(
      content,
      "**Command:**",
      "",
      "```r",
      as.character(target_cmd),
      "```",
      ""
    )

    # Try to extract and show function source
    func_name <- gsub("\\(.*", "", as.character(target_cmd))
    func_name <- trimws(func_name)

    if (exists(func_name, mode = "function")) {
      func_obj <- get(func_name, mode = "function")
      func_source <- deparse(func_obj)

      content <- c(
        content,
        paste0("**Function Source (`", func_name, "`):**"),
        "",
        "```r",
        func_source,
        "```",
        ""
      )
    }

    # Include output summary if requested
    if (include_outputs) {
      tryCatch({
        target_data <- targets::tar_read_raw(target_name)
        output_summary <- capture.output(str(target_data, max.level = 2))

        content <- c(
          content,
          "**Output Structure:**",
          "",
          "```",
          head(output_summary, max_output_lines),
          if (length(output_summary) > max_output_lines) "... (truncated)" else NULL,
          "```",
          ""
        )
      }, error = function(e) {
        # Skip if target not built yet
        content <<- c(content, "*Target not yet built*", "")
      })
    }

    content <- c(content, "---", "")
  }

  # Write to file
  writeLines(content, output_path)

  message("✓ Computational appendix content generated: ", output_path)
  message("  Include in main report with: {{< include ", basename(output_path), " >}}")
  return(output_path)
}
