
#' Application for Simulating Spikes
#'
#' Using runExample() produces a GUI for simulating spikes from our model.
#' @return Shiny app of simulating spikes
#' @export
#'
#' @examples
runExample <- function() {
  appDir <- base::system.file("shiny-examples", "SimulateSpikes", package = "SimSpikes")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `SimSpikes`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}

