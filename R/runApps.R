
#' Application for Simulating Spikes
#'
#' Using runExample() produces a GUI for simulating spikes from our model. This allows for the easy simulation of spikes without using R directly, and the generated spike sequences can be locally downloaded.
#'
#' The application only allows for the implementation of 5 ISI types, namely: Exponential, Gamma, Inverse Gaussian, Log Normal and Weibull. To use alternatively distributions/parameterisations please use \code{\link{Spikes}} or \code{\link{MultiSpikes}}.
#'
#'@param app The type of application you want to open. The options are 'SimulatingSpikes' or 'ViewData'
#'
#' @export
#'
#'
run_app <- function(app) {
  if(app =='SimulateSpikes'){
    appDir <- base::system.file("shiny-examples", "SimulateSpikes", package = "SimSpikes")
  }
  else if(app == 'ViewData') {
    appDir <- base::system.file("shiny-examples", "ViewData", package = "SimSpikes")
  }
  else {
    stop("Could not find App. Valid options are 'SimulateSpikes' and 'ViewData'.")
  }

  shiny::runApp(appDir, display.mode = "normal")
}

