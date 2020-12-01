
#' Application for Simulating Spikes
#'
#' Using runExample() produces a GUI for simulating spikes from our model. This allows for the easy simulation of spikes without using R directly, and the generated spike sequences can be locally downloaded.
#'
#' The application only allows for the implementation of 5 ISI types, namely: Exponential, Gamma, Inverse Gaussian, Log Normal and Weibull. To use alternatively distributions/parameterisations please use \code{\link{Spikes}} or \code{\link{MultiSpikes}}.
#'
#' @export
#'
#'
run_example <- function() {
  appDir <- base::system.file("shiny-examples", "SimulateSpikes", package = "SimSpikes")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `SimSpikes`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

