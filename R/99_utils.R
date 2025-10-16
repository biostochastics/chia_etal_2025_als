#' Utility Functions for ALS Biomarker Investigation
#'
#' Viridis-inspired dark theme with cohesive color palettes
#' Deep purples → cyans → greens → yellows

# Load viridis palette
library(viridisLite)

#' Viridis-inspired dark scientific theme for ggplot2
#'
#' @description
#' Sophisticated dark theme using viridis color palette for scientific figures.
#' Optimized for ultra-dark backgrounds with deep purples, cyans, greens, and lime accents.
#' Designed for both screen display and high-quality figure export.
#'
#' @param base_size Base font size in points
#' @param base_family Base font family
#' @return ggplot2 theme object
#' @export
#'
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_dark_scientific()
theme_dark_scientific <- function(base_size = 12, base_family = "") {
  # Design philosophy: Ultra-dark backgrounds reduce screen glare during extended
  # analysis sessions while maintaining WCAG AA contrast ratios for accessibility
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      # Backgrounds - dark purple-gray matching CSS --bg-surface: #151520
      plot.background = element_rect(fill = "#151520", color = NA),
      panel.background = element_rect(fill = "#151520", color = NA),

      # Grid lines - viridis teal (50%) provides subtle guidance without competing with data
      panel.grid.major = element_line(color = "#26828E", linewidth = 0.35, linetype = "dotted"),
      panel.grid.minor = element_line(color = "#1D1D2E", linewidth = 0.2, linetype = "dotted"),

      # Text elements - viridis lime (80%) for axis titles provides warmth without overwhelming
      text = element_text(color = "#F5F5FA", family = base_family),
      axis.text = element_text(color = "#C8C8D0", size = rel(0.9)),
      axis.title = element_text(color = "#6CCE59", face = "bold", size = rel(1.05)), # viridis lime
      plot.title = element_text(
        color = "#F5F5FA",
        face = "bold",
        size = rel(1.3),
        hjust = 0,
        margin = margin(b = 12),
        lineheight = 1.2
      ),
      plot.subtitle = element_text(
        color = "#9090A0",
        size = rel(1.0),
        hjust = 0,
        margin = margin(b = 10),
        lineheight = 1.3
      ),
      plot.caption = element_text(
        color = "#71717A",
        size = rel(0.85),
        hjust = 1,
        margin = margin(t = 12),
        lineheight = 1.2
      ),

      # Legend - refined with viridis tones and transparent backgrounds
      legend.background = element_rect(
        fill = "transparent",
        color = "#6CCE5933",
        linewidth = 0.75
      ),
      legend.key = element_rect(fill = "transparent", color = NA),
      legend.text = element_text(color = "#C8C8D0", size = rel(0.9)),
      legend.title = element_text(color = "#6CCE59", face = "bold", size = rel(1.0)),
      legend.position = "right",
      legend.justification = "center",

      # Facet strips - transparent backgrounds with viridis borders
      strip.background = element_rect(
        fill = "transparent",
        color = "#31688E",
        linewidth = 1
      ),
      strip.text = element_text(
        color = "#F5F5FA",
        face = "bold",
        size = rel(1.0),
        margin = margin(4, 4, 4, 4)
      ),

      # Panel borders - viridis cyan/teal borders
      panel.border = element_rect(
        color = "#31688E",
        fill = NA,
        linewidth = 1
      ),

      # Plot margins
      plot.margin = margin(15, 15, 15, 15),

      # Additional refinements
      axis.ticks = element_line(color = "#26828E", linewidth = 0.5),
      axis.ticks.length = unit(0.25, "cm"),
      axis.line = element_line(color = "#31688E", linewidth = 0.75)
    )
}

#' Viridis-inspired colorblind-friendly palette for dark backgrounds
#'
#' @description
#' Color palette directly derived from viridis, optimized for ultra-dark backgrounds.
#' Maintains perceptual uniformity and colorblind accessibility while providing
#' sophisticated aesthetic aligned with scientific visualization best practices.
#'
#' Colors span the full viridis gradient:
#' - Deep purples (dark end)
#' - Cyan-teal (mid-dark)
#' - Jade greens (mid-light)
#' - Bright lime-yellow (light end)
#'
#' @return Named vector of 10 hex color codes spanning viridis gradient
#' @export
viridis_dark_palette <- c(
  "#440154", # deep purple (viridis 0%)
  "#482677", # rich purple (viridis ~10%)
  "#3E4A89", # purple-blue (viridis ~20%)
  "#31688E", # deep cyan (viridis ~40%)
  "#26828E", # teal (viridis ~50%)
  "#1F9E89", # jade green (viridis ~60%)
  "#35B779", # green (viridis ~70%)
  "#6CCE59", # bright lime (viridis ~80%)
  "#B5DE2B", # yellow-green (viridis ~90%)
  "#FDE724" # vibrant yellow (viridis 100%)
)

#' Viridis-inspired continuous color scale (fill)
#'
#' @description
#' Continuous viridis fill scale optimized for dark backgrounds,
#' useful for heatmaps, density plots, and gradient visualizations.
#'
#' @param option Viridis palette option: "viridis" (default), "magma", "plasma", "inferno", "cividis"
#' @param ... Additional arguments passed to scale_fill_viridis_c
#' @return ggplot2 scale object
#' @export
scale_fill_viridis_dark <- function(option = "viridis", ...) {
  ggplot2::scale_fill_viridis_c(
    option = option,
    direction = 1,
    begin = 0.1, # Start slightly into palette for better contrast
    end = 0.95, # End slightly before white for visibility
    ...
  )
}

#' Viridis-inspired continuous color scale (color)
#'
#' @description
#' Continuous viridis color scale for lines, points, and edges
#'
#' @param option Viridis palette option
#' @param ... Additional arguments passed to scale_color_viridis_c
#' @return ggplot2 scale object
#' @export
scale_color_viridis_dark <- function(option = "viridis", ...) {
  ggplot2::scale_color_viridis_c(
    option = option,
    direction = 1,
    begin = 0.1,
    end = 0.95,
    ...
  )
}

#' Viridis-inspired discrete color scale (fill)
#'
#' @description
#' Discrete viridis fill scale for categorical data
#'
#' @param option Viridis palette option
#' @param ... Additional arguments passed to scale_fill_viridis_d
#' @return ggplot2 scale object
#' @export
scale_fill_viridis_discrete <- function(option = "viridis", ...) {
  ggplot2::scale_fill_viridis_d(
    option = option,
    direction = 1,
    begin = 0.15,
    end = 0.9,
    ...
  )
}

#' Viridis-inspired discrete color scale (color)
#'
#' @description
#' Discrete viridis color scale for categorical data
#'
#' @param option Viridis palette option
#' @param ... Additional arguments passed to scale_color_viridis_d
#' @return ggplot2 scale object
#' @export
scale_color_viridis_discrete <- function(option = "viridis", ...) {
  ggplot2::scale_color_viridis_d(
    option = option,
    direction = 1,
    begin = 0.15,
    end = 0.9,
    ...
  )
}

#' Custom diverging scale using viridis bookends
#'
#' @description
#' Diverging scale using viridis purple (low), gray (mid), and viridis yellow (high).
#' Useful for showing positive/negative values like log fold-changes.
#'
#' @param ... Additional arguments passed to scale_color_gradient2
#' @return ggplot2 scale object
#' @export
scale_color_diverging_viridis <- function(...) {
  ggplot2::scale_color_gradient2(
    low = "#31688E", # viridis cyan (cool end)
    mid = "#71717A", # neutral gray
    high = "#FDE724", # viridis yellow (warm end)
    midpoint = 0,
    ...
  )
}

#' Custom diverging fill scale using viridis bookends
#'
#' @description
#' Diverging fill scale with viridis palette endpoints
#'
#' @param ... Additional arguments passed to scale_fill_gradient2
#' @return ggplot2 scale object
#' @export
scale_fill_diverging_viridis <- function(...) {
  ggplot2::scale_fill_gradient2(
    low = "#31688E",
    mid = "#71717A",
    high = "#FDE724",
    midpoint = 0,
    ...
  )
}

#' Viridis-inspired color palette for diagnosis groups
#'
#' @description
#' Carefully selected viridis-inspired colors for ALS diagnosis categories,
#' ensuring visual distinction while maintaining cohesive aesthetic.
#'
#' Mapping philosophy:
#' - ALS: Yellow (bright, warm, attention-grabbing - disease focus)
#' - Healthy controls: Lime-green (positive, baseline reference)
#' - Neurological controls: Cyan (cool, distinct from both ALS and healthy)
#'
#' @return Named vector mapping diagnosis categories to viridis hex colors
#' @export
diagnosis_colors_viridis <- c(
  # ALS gets brightest color (yellow 100%) - primary research focus, maximum attention
  "ALS" = "#FDE724",
  # Healthy controls get warm lime (80%) - positive health state, visually distinct from disease
  "Healthy_control" = "#6CCE59",
  # Neurological controls get cool cyan (40%) - alternative disease state, distinct from both
  "Neurological_control" = "#31688E"
)

#' Viridis-inspired color palette for tube types
#'
#' @description
#' Viridis-derived colors emphasizing the technical artifact distinction.
#' Uses contrasting ends of viridis spectrum for maximum differentiation.
#'
#' Mapping philosophy:
#' - EDTA: Purple (dark viridis end - technical factor A)
#' - HEPARIN: Teal (mid viridis - technical factor B)
#'
#' @return Named vector mapping tube types to viridis hex colors
#' @export
tube_type_colors_viridis <- c(
  # EDTA gets dark purple (10%) - confounded with US samples, cooler tone
  "EDTA" = "#482677",
  # HEPARIN gets teal (50%) - confounded with Italy samples, maximum perceptual distance from EDTA
  "HEPARIN" = "#26828E"
)

#' Viridis-inspired color palette for country/geographic origin
#'
#' @description
#' Viridis-based colors for geographic stratification analysis
#'
#' @return Named vector mapping countries to viridis hex colors
#' @export
country_colors_viridis <- c(
  # Italy gets jade (60%) - mid-spectrum, associated with HEPARIN tubes in analysis
  "Italy" = "#1F9E89",
  # US gets yellow-green (90%) - warm end, associated with EDTA tubes, maximum separation from Italy
  "US" = "#B5DE2B"
)

#' Generate viridis palette with n colors
#'
#' @description
#' Utility function to generate n evenly-spaced colors from viridis palette
#'
#' @param n Number of colors to generate
#' @param option Viridis palette option
#' @param alpha Transparency (0-1)
#' @param begin Start point in palette (0-1)
#' @param end End point in palette (0-1)
#' @return Vector of hex colors
#' @export
viridis_palette <- function(n, option = "viridis", alpha = 1, begin = 0.1, end = 0.95) {
  viridisLite::viridis(
    n = n,
    option = option,
    alpha = alpha,
    begin = begin,
    end = end
  )
}

#' Okabe-Ito palette adapted with viridis tones for dark background
#'
#' @description
#' Colorblind-friendly Okabe-Ito palette adapted to harmonize with viridis
#' aesthetic while maintaining accessibility. Use when you need >10 distinct colors
#' or want to maintain Okabe-Ito's proven accessibility.
#'
#' @return Vector of 8 hex colors adapted from Okabe-Ito palette
#' @export
okabe_ito_viridis <- c(
  "#6CCE59", # viridis lime (replaces orange)
  "#31688E", # viridis cyan (replaces sky blue)
  "#1F9E89", # viridis jade (replaces bluish green)
  "#FDE724", # viridis yellow (replaces yellow)
  "#482677", # viridis purple (replaces blue)
  "#B5DE2B", # viridis yellow-green (replaces vermillion)
  "#26828E", # viridis teal (replaces reddish purple)
  "#C8C8D0" # neutral gray (replaces gray)
)

#' Apply viridis theme to existing ggplot
#'
#' @description
#' Convenience function to apply full viridis dark theme to a ggplot object
#'
#' @param p ggplot object
#' @param viridis_option Viridis palette option for color scales
#' @return Modified ggplot object
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, wt, color = cyl)) +
#'   geom_point()
#' apply_viridis_theme(p)
apply_viridis_theme <- function(p, viridis_option = "viridis") {
  p +
    theme_dark_scientific() +
    scale_color_viridis_discrete(option = viridis_option) +
    scale_fill_viridis_discrete(option = viridis_option)
}

#' Create viridis-inspired heatmap color function
#'
#' @description
#' Generate color mapping function for heatmaps using viridis palette
#'
#' @param option Viridis palette option
#' @param n Number of color breaks (default 100 for smooth gradients)
#' @return Function that maps values to colors
#' @export
viridis_heatmap_colors <- function(option = "viridis", n = 100) {
  colors <- viridisLite::viridis(n, option = option, begin = 0.1, end = 0.95)
  colorRampPalette(colors)
}

#' Viridis-inspired theme for correlation matrices
#'
#' @description
#' Specialized theme for correlation heatmaps with viridis colors
#' Optimized for corrplot package
#'
#' @return Named list of corrplot theme parameters (bg, col, tl.col, etc.)
#' @export
corrplot_viridis_theme <- function() {
  list(
    bg = "transparent",
    col = viridis_heatmap_colors(option = "viridis")(200),
    tl.col = "#6CCE59",
    tl.cex = 0.9,
    cl.cex = 0.8,
    addgrid.col = "#26828E",
    mar = c(0, 0, 2, 0)
  )
}

#' Set viridis-inspired default ggplot2 options
#'
#' @description
#' Configure global ggplot2 options to use viridis palette by default
#'
#' @export
set_viridis_defaults <- function() {
  # Set default discrete color scale
  options(
    ggplot2.discrete.colour = viridis_dark_palette,
    ggplot2.discrete.fill = viridis_dark_palette,
    ggplot2.continuous.colour = "viridis",
    ggplot2.continuous.fill = "viridis"
  )

  # Set default theme
  ggplot2::theme_set(theme_dark_scientific())

  invisible(TRUE)
}

#' Extract specific viridis colors by percentage
#'
#' @description
#' Get specific colors from viridis palette by percentage position (0-100)
#'
#' @param percentages Vector of percentages (0-100)
#' @param option Viridis palette option
#' @return Vector of hex colors
#' @export
#'
#' @examples
#' get_viridis_colors(c(0, 50, 100)) # Get start, middle, end colors
get_viridis_colors <- function(percentages, option = "viridis") {
  viridisLite::viridis(
    n = 100,
    option = option
  )[percentages + 1] # +1 because R is 1-indexed
}

#' Create viridis-inspired color gradient for specific value range
#'
#' @description
#' Generate smooth color gradient for specific data range
#'
#' @param values Numeric vector of values to map
#' @param option Viridis palette option
#' @param na_color Color for NA values
#' @return Vector of hex colors
#' @export
map_viridis_colors <- function(values, option = "viridis", na_color = "#71717A") {
  if (all(is.na(values))) {
    return(rep(na_color, length(values)))
  }

  # Normalize values to 0-1 range
  range_vals <- range(values, na.rm = TRUE)
  normalized <- (values - range_vals[1]) / diff(range_vals)

  # Map to viridis colors
  colors <- viridisLite::viridis(100, option = option, begin = 0.1, end = 0.95)
  color_indices <- pmax(1, pmin(100, round(normalized * 99) + 1))
  mapped_colors <- colors[color_indices]

  # Handle NAs
  mapped_colors[is.na(values)] <- na_color

  return(mapped_colors)
}

#' Save ggplot with proper dark background (no white artifacts)
#'
#' @description
#' Saves a ggplot object as PNG with proper handling of dark backgrounds to avoid
#' white padding/halos caused by premultiplied alpha in standard PNG devices.
#'
#' Uses ragg::agg_png device which handles alpha correctly and enforces zero
#' plot margins with solid dark background to eliminate white stripes.
#'
#' @param path Output file path
#' @param plot ggplot object to save
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution in dots per inch (default 300)
#' @param bg Background color (default matches theme_dark_scientific)
#' @return Invisibly returns the plot
#' @export
#'
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_dark_scientific()
#' save_dark_png("myplot.png", p, width = 10, height = 8)
save_dark_png <- function(path, plot, width, height, dpi = 300, bg = "#151520") {
  # Enforce dark backgrounds matching CSS --bg-surface for HTML integration
  plot_final <- plot +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = bg, colour = NA),
      panel.background = ggplot2::element_rect(fill = bg, colour = NA),
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )

  # Use ragg device which handles alpha correctly
  ggplot2::ggsave(
    filename = path,
    plot = plot_final,
    width = width,
    height = height,
    dpi = dpi,
    device = ragg::agg_png,
    bg = bg
  )

  invisible(plot)
}

#' Helper function for centered titles (summary panels and patchwork)
#'
#' @description
#' Returns theme modifications to center plot titles and subtitles.
#' Use this when combining plots with patchwork::plot_annotation() or
#' for standalone summary figures.
#'
#' @return ggplot2 theme object with centered title elements
#' @export
theme_centered_titles <- function() {
  ggplot2::theme(
    plot.title = ggplot2::element_text(hjust = 0.5),
    plot.subtitle = ggplot2::element_text(hjust = 0.5),
    plot.caption = ggplot2::element_text(hjust = 0.5)
  )
}

#' Helper function for zero margins (patchwork panels)
#'
#' @description
#' Returns theme modifications to remove plot margins.
#' Useful when combining plots with patchwork where outer annotation
#' provides spacing.
#'
#' @return ggplot2 theme object with zero margins
#' @export
theme_zero_margins <- function() {
  ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, 0, 0)
  )
}
