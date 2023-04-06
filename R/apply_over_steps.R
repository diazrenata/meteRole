#' Get z scores over steps
#'
#' @param a_role a role model
#'
#' @return data frame of step + z score
#' @export
#'
apply_over_steps <- function(a_role) {

  steps_dfs <- lapply(a_role@modelSteps, FUN = function(a_step) package_timestep(a_step@localComm@indSpecies))

  steps_zs <- lapply(steps_dfs, FUN = function(a_step_df)
    ifelse(nrow(a_step_df) > 10, meteZscore(a_step_df), NA))

  steps_results <- data.frame(step = 1:length(steps_dfs), z = unlist(steps_zs))

  steps_results

}
