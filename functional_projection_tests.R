# Run Hotelling's T2 test after functional observations have been projected.
hotelling_projection_test <- function(x_scores, y_scores, lam = 0) {
  x_scores <- as.matrix(x_scores)
  y_scores <- as.matrix(y_scores)

  n1 <- nrow(x_scores)
  n2 <- nrow(y_scores)
  d <- ncol(x_scores)

  pooled_cov <- ((n1 - 1) * stats::cov(x_scores) +
    (n2 - 1) * stats::cov(y_scores)) / (n1 + n2 - 2)
  pooled_cov <- matrix(pooled_cov, nrow = d, ncol = d)
  mean_diff <- colMeans(x_scores) - colMeans(y_scores)
  statistic <- as.numeric(n1 * n2 / (n1 + n2) *
    t(mean_diff) %*% MASS::ginv(lam * diag(d) + pooled_cov) %*% mean_diff)

  f_statistic <- statistic * (n1 + n2 - d - 1) / (d * (n1 + n2 - 2))
  p_value <- stats::pf(f_statistic, df1 = d, df2 = n1 + n2 - d - 1, lower.tail = FALSE)

  list(
    statistic = statistic,
    f_statistic = f_statistic,
    p.value = p_value,
    projected_dimension = d,
    mean_difference = mean_diff
  )
}

# Build conjugate-gradient directions using the same recursion as the simulation code.
build_cg_directions <- function(x_train, y_train, d = 5, n1_total = nrow(x_train), n2_total = nrow(y_train)) {
  x_train <- as.matrix(x_train)
  y_train <- as.matrix(y_train)

  mean_x <- colMeans(x_train)
  mean_y <- colMeans(y_train)
  mean_diff <- mean_x - mean_y
  x_centered <- sweep(x_train, 2, mean_x, "-")
  y_centered <- sweep(y_train, 2, mean_y, "-")
  covariance_hat <- (crossprod(x_centered) + crossprod(y_centered)) / (n1_total + n2_total - 2)

  phi <- rep(0, length(mean_diff))
  direction <- mean_diff
  residual <- mean_diff
  directions <- vector("list", d)

  for (i in seq_len(d)) {
    denominator <- as.numeric(t(direction) %*% covariance_hat %*% direction)
    step_size <- as.numeric(sum(direction * residual) / denominator)
    phi <- phi + step_size * direction
    residual <- mean_diff - as.vector(covariance_hat %*% phi)
    directions[[i]] <- phi

    cg_weight <- as.numeric(t(residual) %*% covariance_hat %*% direction / denominator)
    direction <- residual - cg_weight * direction
  }

  do.call(cbind, directions)
}

# Run the sample-splitting conjugate-gradient test.
cg_projection_test <- function(x, y, d = 5, n_splits = 100, split_fraction = 1 / 2, seed = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  n1 <- nrow(x)
  n2 <- nrow(y)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  split_statistics <- replicate(n_splits, {
    x_index <- sample(n1, floor(n1 * split_fraction), replace = FALSE)
    y_index <- sample(n2, floor(n2 * split_fraction), replace = FALSE)

    cg_directions <- build_cg_directions(
      x_train = x[-x_index, , drop = FALSE],
      y_train = y[-y_index, , drop = FALSE],
      d = d,
      n1_total = n1,
      n2_total = n2
    )
    x_scores <- x[x_index, , drop = FALSE] %*% cg_directions
    y_scores <- y[y_index, , drop = FALSE] %*% cg_directions

    vapply(seq_len(d), function(j) {
      x_score <- x_scores[, j]
      y_score <- y_scores[, j]
      mean_diff <- mean(x_score) - mean(y_score)
      variance_hat <- ((n1 - 1) * stats::var(x_score) + (n2 - 1) * stats::var(y_score)) /
        (n1 + n2 - 2)
      sqrt(n1 * n2 / (n1 + n2)) * abs(mean_diff) / sqrt(variance_hat)
    }, numeric(1))
  })

  statistics <- rowMeans(split_statistics)
  p_values <- stats::pnorm(statistics, lower.tail = FALSE)

  list(
    statistic = statistics[d],
    p.value = p_values[d],
    projected_dimension = d,
    statistics_by_dimension = statistics,
    p.values_by_dimension = p_values,
    method = "cg",
    selection = "cg_sample_splitting",
    n_splits = n_splits,
    split_fraction = split_fraction
  )
}

# Build FPCA projection functions using the overall covariance operator.
build_fpca_basis <- function(x, y,
                             d = NULL,
                             ranking = c("q", "lambda"),
                             variance_threshold = 0.85,
                             min_eigenvalue = 1e-8) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  ranking <- match.arg(ranking)

  xy <- rbind(x, y)
  overall_cov <- crossprod(scale(xy, center = colMeans(xy), scale = FALSE)) / (nrow(xy) - 1)
  eig <- eigen(overall_cov, symmetric = TRUE)
  keep <- which(eig$values > min_eigenvalue)

  values <- eig$values[keep]
  vectors <- eig$vectors[, keep, drop = FALSE]
  mean_diff <- colMeans(x) - colMeans(y)
  max_d <- min(length(values), nrow(x) + nrow(y) - 2)

  if (ranking == "q") {
    q_scores <- (drop(crossprod(vectors, mean_diff))^2) / values
    ordering <- order(q_scores, decreasing = TRUE)
    if (is.null(d)) {
      ordered_scores <- q_scores[ordering][seq_len(max_d)]
      ordered_scores <- ordered_scores[is.finite(ordered_scores) & ordered_scores > 0]
      d <- if (length(ordered_scores) <= 1) {
        max(1, length(ordered_scores))
      } else {
        which.max(ordered_scores[-length(ordered_scores)] / ordered_scores[-1])
      }
    }
  } else {
    ordering <- seq_along(values)
    if (is.null(d)) {
      variance_share <- cumsum(values[ordering]) / sum(values[ordering])
      d <- which(variance_share >= variance_threshold)[1]
    }
  }

  d <- min(d, max_d)
  selected <- ordering[seq_len(d)]
  list(
    basis = vectors[, selected, drop = FALSE],
    selected_dimension = d,
    criterion = paste("fpca", ranking, sep = "_"),
    eigenvalues = values[selected]
  )
}

# Build B-spline projection functions using the manuscript simulation setting.
build_bspline_basis <- function(grid, d = 6) {
  if (!requireNamespace("fda", quietly = TRUE)) {
    stop("Package 'fda' is required for the B-spline projection method.", call. = FALSE)
  }

  order <- min(d, 4)
  knots <- seq(0, 1, 1 / (d - order + 1))
  basis_object <- fda::create.bspline.basis(c(0, 1), nbasis = d, norder = order, breaks = knots)
  basis <- fda::eval.basis(grid, basis_object)

  list(
    basis = basis,
    selected_dimension = ncol(basis),
    criterion = "bspline"
  )
}

# Project functional observations onto a chosen set of projection functions.
project_curves <- function(curves, basis) {
  as.matrix(curves) %*% as.matrix(basis)
}

# Fit one projection-based functional two-sample test.
functional_two_sample_test <- function(x, y,
                                       method = c(
                                         "bspline",
                                         "cg",
                                         "fpca_q",
                                         "fpca_lambda"
                                       ),
                                       d = NULL,
                                       grid = NULL,
                                       variance_threshold = 0.85,
                                       lam = 0,
                                       seed = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  method <- match.arg(method)

  if (is.null(grid)) {
    grid <- seq(0, 1, length.out = ncol(x))
  }

  if (method == "cg") {
    return(cg_projection_test(x, y, d = if (is.null(d)) 5 else d, seed = seed))
  }

  projection <- switch(
    method,
    bspline = build_bspline_basis(grid, d = if (is.null(d)) 6 else d),
    fpca_q = build_fpca_basis(x, y, d = d, ranking = "q"),
    fpca_lambda = build_fpca_basis(
      x, y, d = d, ranking = "lambda",
      variance_threshold = variance_threshold
    )
  )

  x_scores <- project_curves(x, projection$basis)
  y_scores <- project_curves(y, projection$basis)
  test <- hotelling_projection_test(x_scores, y_scores, lam = lam)

  test$method <- method
  test$selection <- projection$criterion
  test$basis <- projection$basis
  test
}

# Plot generated two-sample functional curves and their sample means.
plot_functional_two_sample_data <- function(data,
                                            file = NULL,
                                            main = "Generated two-sample functional curves") {
  grid <- data$grid
  x <- as.matrix(data$x)
  y <- as.matrix(data$y)

  if (!is.null(file)) {
    grDevices::png(file, width = 1000, height = 650, res = 130)
    on.exit(grDevices::dev.off(), add = TRUE)
  }

  x_color <- grDevices::adjustcolor("#2C7FB8", alpha.f = 0.22)
  y_color <- grDevices::adjustcolor("#D95F02", alpha.f = 0.22)
  mean_x_color <- "#08519C"
  mean_y_color <- "#A63603"

  graphics::matplot(
    grid, t(x),
    type = "l",
    lty = 1,
    col = x_color,
    ylim = range(x, y, finite = TRUE),
    xlab = "Time",
    ylab = "Function value",
    main = main
  )
  graphics::matlines(grid, t(y), lty = 1, col = y_color)
  graphics::lines(grid, colMeans(x), col = mean_x_color, lwd = 3)
  graphics::lines(grid, colMeans(y), col = mean_y_color, lwd = 3)
  graphics::legend(
    "topright",
    legend = c("Group 1 mean", "Group 2 mean"),
    col = c(mean_x_color, mean_y_color),
    lwd = 3,
    bty = "n"
  )
  invisible(NULL)
}

# Generate two groups of smooth functional observations.
simulate_functional_two_sample_data <- function(n1 = 25,
                                                n2 = 25,
                                                grid = seq(0, 1, length.out = 50),
                                                signal_strength = 0.55,
                                                noise_sd = 0.35,
                                                seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  p <- length(grid)
  shared_basis <- cbind(
    sin(2 * pi * grid),
    cos(2 * pi * grid),
    sin(4 * pi * grid)
  )
  signal <- signal_strength * dnorm(grid, mean = 0.6, sd = 0.08)
  signal <- signal / max(signal)

  make_group <- function(n, mean_shift) {
    latent <- matrix(stats::rnorm(n * ncol(shared_basis)), nrow = n)
    latent <- sweep(latent, 2, c(1.2, 0.8, 0.5), `*`)
    smooth_noise <- matrix(stats::rnorm(n * p, sd = noise_sd), nrow = n)
    smooth_noise <- 0.6 * smooth_noise +
      0.4 * t(apply(smooth_noise, 1, stats::filter,
        filter = rep(1 / 5, 5), sides = 2, circular = TRUE
      ))
    latent %*% t(shared_basis) + matrix(mean_shift, nrow = n, ncol = p, byrow = TRUE) + smooth_noise
  }

  list(
    grid = grid,
    x = make_group(n1, mean_shift = rep(0, p)),
    y = make_group(n2, mean_shift = signal)
  )
}
