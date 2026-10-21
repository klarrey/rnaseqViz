#' Launch the rnaseqViz Shiny Application
#'
#' Launches the interactive breast cancer RNA-seq
#' visualization application.
#'
#' @return Invisibly returns the Shiny application result.
#' @export
run_app <- function() {

  if (!requireNamespace("shiny", quietly = TRUE)) {
    stop(
      "Please install the shiny package.",
      call. = FALSE
    )
  }

  app_dir <- system.file(
    "shiny",
    package = "rnaseqViz"
  )

  if (app_dir == "") {
    stop("Shiny application directory not found.")
  }

  shiny::runApp(app_dir)
}
