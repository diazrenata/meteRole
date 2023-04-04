#' package timestep
#'
#' @param row_from_modelSteps row from the modelSteps slot of a roleR object
#'
#' @return row tabled and tallied
#' @export
#'
#' @examples
#'
#' row_from_modelSteps <- c(1, 1, 1, 8, 10, 10)
#' package_timestep(row_from_modelSteps)
#'
package_timestep <- function(row_from_modelSteps) {

  table(row_from_modelSteps) |>
    as.data.frame()

}
