visualize_dim_selection <- function(seurat_obj,
                                    reduction = "pca",
                                    max_dims = 50,
                                    variance_thresholds = c(0.60, 0.70, 0.80, 0.90)) {

  # --- 1. Input Validation and Data Extraction ---
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(paste("Reduction '", reduction, "' not found in Seurat object. Please run the appropriate reduction method first.", sep=""))
  }

  if (!requireNamespace("ggplot2", quietly = TRUE) || !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'patchwork' are required but not installed.")
  }

  reduction_obj <- seurat_obj@reductions[[reduction]]

  # Check for importance scores (StDev for PCA, Singular Values for LSI)
  if (length(reduction_obj@stdev) > 0) {
    importance_scores <- reduction_obj@stdev
    y_label <- "Standard Deviation"
    plot_title <- "1. Elbow Plot (Standard Deviation)"
  } else if (length(reduction_obj@S) > 0) {
    importance_scores <- reduction_obj@S
    y_label <- "Singular Value"
    plot_title <- "1. Elbow Plot (Singular Values)"
  } else {
    stop(paste("Could not find importance scores (@stdev or @S) in reduction '", reduction, "'.", sep=""))
  }

  # Ensure max_dims is not greater than available dimensions
  max_dims <- min(max_dims, length(importance_scores))

  # Extract data for the specified number of dimensions
  scores_subset <- importance_scores[1:max_dims]
  variance <- scores_subset^2
  total_var <- sum(importance_scores^2) # Use all dimensions for total variance
  pct_var <- variance / total_var * 100
  cumulative_pct <- cumsum(variance) / total_var * 100

  # --- 2. Calculate Metrics ---
  marginal_gain <- c(pct_var[1], diff(cumulative_pct))
  diff1 <- diff(scores_subset)
  diff2 <- diff(diff1)
  curvature <- c(0, 0, diff2, rep(NA, max_dims - length(diff2) - 2))

  plot_data <- data.frame(
    Dimension = 1:max_dims,
    importance = scores_subset,
    individual_var = pct_var,
    cumulative_var = cumulative_pct,
    marginal_gain = marginal_gain,
    curvature = curvature
  )

  # --- 3. Annotations for Plots ---
  cumulative_pct_full <- cumsum(importance_scores^2) / total_var
  dim_at_thresholds <- sapply(variance_thresholds, function(thresh) {
    which(cumulative_pct_full >= thresh)[1]
  })

  annotation_data <- data.frame(
    dim = dim_at_thresholds,
    threshold = variance_thresholds * 100,
    label = paste0(dim_at_thresholds, " Dims\n(", variance_thresholds * 100, "%)")
  )

  colors <- c("red", "orange", "green", "blue", "purple")

  cat("\n=== Automatic Dimension Selection Based on Variance Thresholds ===\n")
  for (i in 1:nrow(annotation_data)) {
    cat(sprintf("  %.0f%% variance → Dimension %d\n",
                annotation_data$threshold[i],
                annotation_data$dim[i]))
  }
  cat("\n")

  # --- 4. Generate Plots ---
  # Plot 1: Elbow Plot
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dimension, y = importance)) +
    ggplot2::geom_line(color = "steelblue", size = 1) +
    ggplot2::geom_point(color = "steelblue", size = 2)

  for (i in 1:nrow(annotation_data)) {
    if (!is.na(annotation_data$dim[i]) && annotation_data$dim[i] <= max_dims) {
      p1 <- p1 +
        ggplot2::geom_vline(xintercept = annotation_data$dim[i], linetype = "dashed", color = colors[i], alpha = 0.6) +
        ggplot2::annotate("text", x = annotation_data$dim[i], y = max(scores_subset) * (1 - i * 0.1),
                          label = annotation_data$label[i], color = colors[i], size = 3, hjust = -0.1)
    }
  }
  p1 <- p1 + ggplot2::labs(title = plot_title, x = "Dimension", y = y_label) + ggplot2::theme_bw()

  # Plot 2: Cumulative Variance
  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dimension, y = cumulative_var)) +
    ggplot2::geom_line(color = "darkgreen", size = 1) +
    ggplot2::geom_point(color = "darkgreen", size = 1.5)
  for (i in 1:length(variance_thresholds)) {
    p2 <- p2 + ggplot2::geom_hline(yintercept = variance_thresholds[i] * 100, linetype = "dashed", color = colors[i], alpha = 0.6)
  }
  p2 <- p2 + ggplot2::scale_y_continuous(limits = c(0, 100)) +
    ggplot2::labs(title = "2. Cumulative Variance Explained", x = "Dimension", y = "Cumulative Variance (%)") + ggplot2::theme_bw()

  # Plot 3: Marginal Gain
  p3 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Dimension, y = marginal_gain)) +
    ggplot2::geom_bar(stat = "identity", fill = "coral", alpha = 0.7) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    ggplot2::labs(title = "3. Marginal Gain per Dimension", x = "Dimension", y = "Variance Gained (%)") + ggplot2::theme_bw()

  # Plot 4: Curvature
  p4 <- ggplot2::ggplot(plot_data[3:max_dims, ], ggplot2::aes(x = Dimension, y = abs(curvature))) +
    ggplot2::geom_line(color = "purple", size = 1) +
    ggplot2::geom_point(color = "purple", size = 1.5) +
    ggplot2::labs(title = "4. Curvature (2nd Derivative)", x = "Dimension", y = "Absolute Curvature") + ggplot2::theme_bw()

  # --- 5. Combine and Print ---
  combined <- (p1 | p2) / (p3 | p4)
  print(combined)

  # --- 6. Recommendations ---
  cat("\n=== Recommendations Based on Analysis ===\n")
  threshold_point <- which(plot_data$marginal_gain < 0.5)[1]
  if (!is.na(threshold_point)) {
    cat(sprintf("• Marginal gain <0.5%%: Dimension %d\n", threshold_point))
  }

  # Recommend elbow point within a reasonable range (e.g., first 30 dimensions)
  elbow_range <- 3:min(30, max_dims)
  elbow_point <- which.max(abs(plot_data$curvature[elbow_range])) + (elbow_range[1] - 1)
  cat(sprintf("• Maximum curvature (elbow point): Dimension %d\n", elbow_point))

  for (i in 1:nrow(annotation_data)) {
    cat(sprintf("• %.0f%% variance: Dimension %d\n",
                annotation_data$threshold[i],
                annotation_data$dim[i]))
  }
  cat("\n")

  return(plot_data)
}