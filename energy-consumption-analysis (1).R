###############################################
# Appliance Energy Modelling in R
# Author: Priyanka Chandrahar Mane
# SIG718 assessment analysis and portfolio validation extension
###############################################

# Run from the folder containing this script, ENB.txt and AggWaFit718.R.
# In RStudio: Session > Set Working Directory > To Source File Location.
# The course helper and course-modified data must be obtained separately.
# Dependency: install.packages("lpSolve")
if (!requireNamespace("lpSolve", quietly = TRUE)) {
  stop('Install the required package first: install.packages("lpSolve")')
}
required_files <- c("ENB.txt", "AggWaFit718.R")
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop(paste("Missing required local files:", paste(missing_files, collapse = ", ")))
}
source("AggWaFit718.R")

#------------------------------------------------
# T1 – Understand the data
#------------------------------------------------

# Read the full ENB dataset as a numeric matrix
the.data <- as.matrix(read.table("ENB.txt"))

# Fixed sampling seed retained to reproduce the original analysis
set.seed(226162354)

# Required sample size for this assessment
num_row <- 650

# Randomly select 650 rows and keep only columns 2–7 (X1–X5 and Y)
my.data <- as.data.frame(
  the.data[sample(1:nrow(the.data), num_row), 2:7]
)

# Give meaningful names to the variables
colnames(my.data) <- c("X1_temp_LR",
                       "X2_humidity_LR",
                       "X3_temp_OFF",
                       "X4_humidity_OFF",
                       "X5_pressure",
                       "Y_energy")

# Quick structural check: types and first few values
str(my.data)

# Summary statistics for all variables (min, max, quartiles, mean)
summary(my.data)

# Scatter-plot matrix: shows pairwise relationships between X1–X5 and Y
pairs(my.data, main = "Scatter plot matrix – X1–X5 and Y")

# Histograms to inspect distribution of each variable
par(mfrow = c(2, 3))
for (i in 1:6) {
  hist(my.data[, i],
       main = paste("Histogram –", colnames(my.data)[i]),
       xlab  = colnames(my.data)[i])
}
par(mfrow = c(1, 1))   # reset plotting layout

#------------------------------------------------
# T2 – Transform the data
#------------------------------------------------

# Choose FOUR predictors for aggregation models (as required in T2)
selected_X <- c("X1_temp_LR",
                "X2_humidity_LR",
                "X3_temp_OFF",
                "X4_humidity_OFF")

# Keep the chosen Xs and the variable of interest Y
data_T2_raw <- my.data[, c(selected_X, "Y_energy")]

# Apply log transformation where appropriate
# Retain the transformations selected in the original assessment analysis.
data_T2_raw$X1_temp_LR <- log(data_T2_raw$X1_temp_LR)
data_T2_raw$Y_energy   <- log(data_T2_raw$Y_energy)

# Min–max normalisation function to scale each variable to [0, 1]
minmax <- function(x) { (x - min(x)) / (max(x) - min(x)) }

# Apply min–max scaling to all 5 variables (4 X + Y)
data_T2_scaled <- as.data.frame(lapply(data_T2_raw, minmax))

# Convert to matrix because AggWaFit functions expect a numeric matrix
data_T2_scaled <- as.matrix(data_T2_scaled)

# Store min/max of each transformed variable so we can invert scaling in T4
scale_params <- data.frame(
  variable = colnames(data_T2_raw),
  min = apply(data_T2_raw, 2, min),
  max = apply(data_T2_raw, 2, max)
)

# Save the transformed data to text file as required ("name-transformed.txt")
write.table(data_T2_scaled,
            file = "Priyanka-transformed.txt",
            row.names = FALSE,
            col.names = FALSE)

#------------------------------------------------
# T3 – Build models and compare
#------------------------------------------------

# Read back the transformed data from file for model fitting
data_T3 <- as.matrix(read.table("Priyanka-transformed.txt"))

# First 4 columns = inputs (chosen Xs), 5th column = output (Y)
X <- data_T3[, 1:4]
Y <- data_T3[, 5]

# Helper function to calculate RMSE between predicted and observed Y
rmse_fun <- function(a, b) sqrt(mean((a - b)^2))

# Number of input variables in the aggregation models
nvar <- 4

# --- 1) WAM (Weighted Arithmetic Mean) -----------------------
fit.QAM(data_T3,
        output.1 = "Priyanka_WAM_out.txt",
        stats.1  = "Priyanka_WAM_stats.txt",
        g = AM, g.inv = invAM)

# Read stats file and extract the learned weights w1–w4
wam_stats <- read.table("Priyanka_WAM_stats.txt",
                        header = FALSE, fill = TRUE)
wam_w <- as.numeric(wam_stats[6:(5 + nvar), 2])

# Get fitted Y values using the learned WAM weights
pred_WAM <- apply(X, 1, function(x_row)
  QAM(x_row, w = wam_w, g = AM, g.inv = invAM))

# Compute RMSE for WAM
rmse_WAM <- rmse_fun(pred_WAM, Y)

# --- 2) WPM(0.5) – Weighted Power Mean with p = 0.5 ----------
fit.QAM(data_T3,
        output.1 = "Priyanka_WPM05_out.txt",
        stats.1  = "Priyanka_WPM05_stats.txt",
        g = PM05, g.inv = invPM05)

# Extract WPM(0.5) weights
wpm05_stats <- read.table("Priyanka_WPM05_stats.txt",
                          header = FALSE, fill = TRUE)
wpm05_w <- as.numeric(wpm05_stats[6:(5 + nvar), 2])

# Predicted Y using WPM(0.5)
pred_WPM05 <- apply(X, 1, function(x_row)
  QAM(x_row, w = wpm05_w, g = PM05, g.inv = invPM05))

# RMSE for WPM(0.5)
rmse_WPM05 <- rmse_fun(pred_WPM05, Y)

# --- 3) WPM(2) – Weighted Power Mean with p = 2 --------------
fit.QAM(data_T3,
        output.1 = "Priyanka_WPM2_out.txt",
        stats.1  = "Priyanka_WPM2_stats.txt",
        g = QM, g.inv = invQM)

# Extract WPM(2) weights
wpm2_stats <- read.table("Priyanka_WPM2_stats.txt",
                         header = FALSE, fill = TRUE)
wpm2_w <- as.numeric(wpm2_stats[6:(5 + nvar), 2])

# Predicted Y using WPM(2)
pred_WPM2 <- apply(X, 1, function(x_row)
  QAM(x_row, w = wpm2_w, g = QM, g.inv = invQM))

# RMSE for WPM(2)
rmse_WPM2 <- rmse_fun(pred_WPM2, Y)

# --- 4) OWA – Ordered Weighted Average -----------------------
fit.OWA(data_T3,
        output.1 = "Priyanka_OWA_out.txt",
        stats.1  = "Priyanka_OWA_stats.txt")

# Extract OWA weights
owa_stats <- read.table("Priyanka_OWA_stats.txt",
                        header = FALSE, fill = TRUE)
owa_w <- as.numeric(owa_stats[7:(6 + nvar), 2])

# Predicted Y using OWA weights
pred_OWA <- apply(X, 1, function(x_row)
  OWA(x_row, w = owa_w))

# Fitting-sample RMSE for OWA
rmse_OWA <- rmse_fun(pred_OWA, Y)

# Collect the four models and their RMSE into a comparison table
T3_results <- data.frame(
  Model = c("WAM", "WPM(0.5)", "WPM(2)", "OWA"),
  RMSE  = c(rmse_WAM, rmse_WPM05, rmse_WPM2, rmse_OWA)
)

print(T3_results)   # Fitting-sample comparison, not holdout performance

# Identify best model based on minimum RMSE
best_model_index <- which.min(T3_results$RMSE)
best_model_name  <- T3_results$Model[best_model_index]

#------------------------------------------------
# T4 – Prediction
#------------------------------------------------

# New observation (original scale) provided in the question
new_original <- data.frame(
  X1_temp_LR      = 23.50,
  X2_humidity_LR  = 35.87125,
  X3_temp_OFF     = 24.89,
  X4_humidity_OFF = 32.93
  # Note: X5 is not used in the aggregation models, so we omit it here
)

# Apply the same log transformation used in T2
new_trans <- new_original
new_trans$X1_temp_LR <- log(new_trans$X1_temp_LR)

# Helper function: apply min–max scaling to a single value
scale_one <- function(x_raw, var_name) {
  minv <- scale_params$min[scale_params$variable == var_name]
  maxv <- scale_params$max[scale_params$variable == var_name]
  (x_raw - minv) / (maxv - minv)
}

# Build scaled input vector (length 4) in the same order as in T2/T3
new_scaled_vec <- c(
  scale_one(new_trans$X1_temp_LR,      "X1_temp_LR"),
  scale_one(new_trans$X2_humidity_LR,  "X2_humidity_LR"),
  scale_one(new_trans$X3_temp_OFF,     "X3_temp_OFF"),
  scale_one(new_trans$X4_humidity_OFF, "X4_humidity_OFF")
)

# Use the best model from T3 to predict scaled Y
if (best_model_name == "WAM") {
  pred_Y_scaled <- QAM(new_scaled_vec, w = wam_w, g = AM, g.inv = invAM)
} else if (best_model_name == "WPM(0.5)") {
  pred_Y_scaled <- QAM(new_scaled_vec, w = wpm05_w, g = PM05, g.inv = invPM05)
} else if (best_model_name == "WPM(2)") {
  pred_Y_scaled <- QAM(new_scaled_vec, w = wpm2_w, g = QM, g.inv = invQM)
} else if (best_model_name == "OWA") {
  pred_Y_scaled <- OWA(new_scaled_vec, w = owa_w)
}

# Undo min–max scaling for Y (back to log-scale)
Y_min_trans <- scale_params$min[scale_params$variable == "Y_energy"]
Y_max_trans <- scale_params$max[scale_params$variable == "Y_energy"]

pred_Y_trans    <- as.numeric(pred_Y_scaled) * (Y_max_trans - Y_min_trans) + Y_min_trans

# Undo log transformation to get prediction in original Wh units
pred_Y_original <- exp(pred_Y_trans)

# Given measured Y for this case
measured_Y <- 120

# Print a concise prediction summary (for T4 slide)
cat("\n----------- T4 prediction summary -----------\n")
cat("Best model:                ", best_model_name, "\n")
cat("Predicted Y (original):    ", pred_Y_original, "\n")
cat("Measured Y (given):        ", measured_Y, "\n")
cat("Absolute error |Y - Yhat|: ", abs(measured_Y - pred_Y_original), "\n")
cat("---------------------------------------------\n")
# Baseline comparison on the fitting sample
baseline_pred <- rep(mean(Y), length(Y))

baseline_comparison <- rbind(
  data.frame(
    Model = "Mean baseline",
    RMSE = rmse_fun(baseline_pred, Y)
  ),
  T3_results
)

print(baseline_comparison)
# Portfolio extension: training and validation split
set.seed(718)

train_rows <- sample(
  seq_len(nrow(my.data)),
  size = floor(0.8 * nrow(my.data))
)

train_raw <- my.data[train_rows, ]
validation_raw <- my.data[-train_rows, ]

cat("\nTraining rows:", nrow(train_raw))
cat("\nValidation rows:", nrow(validation_raw), "\n")
# Transform both groups consistently
model_columns <- c(selected_X, "Y_energy")

train_trans <- train_raw[, model_columns]
validation_trans <- validation_raw[, model_columns]

for (variable in c("X1_temp_LR", "Y_energy")) {
  stopifnot(
    all(train_trans[[variable]] > 0),
    all(validation_trans[[variable]] > 0)
  )
  train_trans[[variable]] <- log(train_trans[[variable]])
  validation_trans[[variable]] <- log(validation_trans[[variable]])
}

# Learn scaling parameters from training data only
train_min <- sapply(train_trans, min)
train_range <- sapply(train_trans, max) - train_min
stopifnot(all(train_range > 0))

scale_with_training <- function(data) {
  sweep(
    sweep(as.matrix(data), 2, train_min, "-"),
    2, train_range, "/"
  )
}

train_scaled <- scale_with_training(train_trans)
validation_scaled <- scale_with_training(validation_trans)

# Check validation values outside the training range
print(colSums(validation_scaled < 0 | validation_scaled > 1))
# Inspect validation inputs outside the training range
outside_counts <- data.frame(
  Variable = selected_X,
  Below_zero = colSums(
    validation_scaled[, selected_X, drop = FALSE] < 0
  ),
  Above_one = colSums(
    validation_scaled[, selected_X, drop = FALSE] > 1
  )
)

print(outside_counts)
print(range(validation_scaled[, "X3_temp_OFF"]))
# Fit WAM using training data only
fit.QAM(
  train_scaled,
  output.1 = "validation_WAM_fit.txt",
  stats.1 = "validation_WAM_stats.txt",
  g = AM, g.inv = invAM
)

validation_wam_stats <- read.table(
  "validation_WAM_stats.txt",
  header = FALSE, fill = TRUE
)

train_wam_weights <- as.numeric(validation_wam_stats[6:9, 2])
stopifnot(all(is.finite(train_wam_weights)))

# Predict validation observations
validation_wam_scaled <- apply(
  validation_scaled[, selected_X, drop = FALSE],
  1,
  function(x) QAM(x, w = train_wam_weights, g = AM, g.inv = invAM)
)

# Convert predictions back to Wh
validation_wam_wh <- exp(
  validation_wam_scaled * train_range["Y_energy"] +
    train_min["Y_energy"]
)

actual_wh <- validation_raw$Y_energy
baseline_wh <- rep(mean(train_raw$Y_energy), length(actual_wh))

validation_results <- data.frame(
  Model = c("Training-mean baseline", "WAM"),
  RMSE_Wh = c(
    rmse_fun(baseline_wh, actual_wh),
    rmse_fun(validation_wam_wh, actual_wh)
  ),
  MAE_Wh = c(
    mean(abs(baseline_wh - actual_wh)),
    mean(abs(validation_wam_wh - actual_wh))
  )
)

print(validation_results)
# Fit the remaining aggregation models on training data only
model_specs <- list(
  "WPM(0.5)" = list(g = PM05, inverse = invPM05),
  "WPM(2)" = list(g = QM, inverse = invQM),
  "OWA" = NULL
)

for (model_name in names(model_specs)) {
  model_spec <- model_specs[[model_name]]
  file_tag <- gsub("[^A-Za-z0-9]", "", model_name)
  stats_file <- paste0("validation_", file_tag, "_stats.txt")
  fit_file <- paste0("validation_", file_tag, "_fit.txt")
  
  if (model_name == "OWA") {
    fit.OWA(train_scaled, output.1 = fit_file, stats.1 = stats_file)
    weight_rows <- 7:10
  } else {
    fit.QAM(
      train_scaled, output.1 = fit_file, stats.1 = stats_file,
      g = model_spec$g, g.inv = model_spec$inverse
    )
    weight_rows <- 6:9
  }
  
  model_stats <- read.table(stats_file, header = FALSE, fill = TRUE)
  weights <- as.numeric(model_stats[weight_rows, 2])
  stopifnot(
    length(weights) == length(selected_X),
    all(is.finite(weights)),
    all(weights >= 0),
    abs(sum(weights) - 1) < 1e-6
  )
  
  predictions_scaled <- apply(
    validation_scaled[, selected_X, drop = FALSE], 1,
    function(x) {
      if (model_name == "OWA") {
        OWA(x, w = weights)
      } else {
        QAM(x, w = weights,
            g = model_spec$g, g.inv = model_spec$inverse)
      }
    }
  )
  
  predictions_wh <- exp(
    predictions_scaled * train_range["Y_energy"] +
      train_min["Y_energy"]
  )
  stopifnot(all(is.finite(predictions_wh)))
  
  validation_results <- rbind(
    validation_results,
    data.frame(
      Model = model_name,
      RMSE_Wh = rmse_fun(predictions_wh, actual_wh),
      MAE_Wh = mean(abs(predictions_wh - actual_wh))
    )
  )
}

print(validation_results)

# Export the exploratory holdout results in original energy units.
write.csv(validation_results, "validation-results.csv", row.names = FALSE)
