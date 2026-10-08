#' Add a weight variable using the sample frame (fast)
#'
#' data.table version of `analysistools::add_weights()`, same arguments,
#' checks and formula:
#' \deqn{w_h = \frac{N_h / \sum N}{n_h / \sum n}}{w_h = (N_h / sum(N)) / (n_h / sum(n))}
#' where \eqn{N_h} is the population of stratum \eqn{h} in the sample frame and
#' \eqn{n_h} the number of interviews in the dataset. Weights therefore sum to
#' the number of interviews (normalised weights). Row order is preserved.
#'
#' @param dataset The clean dataset.
#' @param sample_data Sample frame with the population per stratum.
#' @param strata_column_dataset Strata column in the dataset.
#' @param strata_column_sample Strata column in the sample frame.
#' @param population_column Population column in the sample frame.
#' @param weight_column Name of the weight column to add (default "weights").
#'
#' @return The dataset (data frame) with the weight column added.
#' @export
#' @examples
#' clean_data <- data.frame(uuid = 1:8, strata = c("s1", "s2", "s1", "s2", "s1", "s2", "s1", "s1"))
#' sample <- data.frame(strata = c("s1", "s2"), population = c(30000, 50000))
#' add_weights_fast(clean_data, sample, "strata", "strata", "population")
add_weights_fast <- function(dataset,
                             sample_data,
                             strata_column_dataset = NULL,
                             strata_column_sample = NULL,
                             population_column = NULL,
                             weight_column = "weights") {
  dataset <- as.data.frame(dataset)
  if (!strata_column_sample %in% names(sample_data)) {
    stop("Cannot find the defined strata column in the provided sample frame.")
  }
  if (!strata_column_dataset %in% names(dataset)) {
    stop("Cannot find the defined strata column in the provided dataset.")
  }
  if (!all(dataset[[strata_column_dataset]] %in% sample_data[[strata_column_sample]])) {
    stop("Not all strata from dataset are in sample frame")
  }
  if (!all(sample_data[[strata_column_sample]] %in% dataset[[strata_column_dataset]])) {
    stop("Not all strata from sample frame are in dataset")
  }
  if (!population_column %in% names(sample_data)) {
    stop("Cannot find the defined population_column column in the provided sample frame.")
  }
  if (weight_column %in% names(dataset)) {
    stop("Weight column already exists in the dataset. Please input another weights column")
  }
  cnt <- data.table(s = dataset[[strata_column_dataset]])[, list(count = .N), by = "s"]
  frame <- data.table(s = sample_data[[strata_column_sample]],
                      population = as.numeric(sample_data[[population_column]]))
  frame <- frame[, list(population = sum(population)), by = "s"]
  frame <- cnt[frame, on = "s"]
  frame[, wt := (population / sum(population)) / (as.numeric(count) / sum(as.numeric(count)))]
  dataset[[weight_column]] <- frame$wt[match(dataset[[strata_column_dataset]], frame$s)]
  if (anyNA(dataset[[weight_column]])) stop("There are NA values in the weights column")
  dataset
}
