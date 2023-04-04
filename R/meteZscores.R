#' meteZscore
#'
#' largely copied from meteR, calculates a Z score for obs SAD loglik compared to METE sads loglik
#'
#' @param packaged_timestep result of package_timestep
#' @param nrep n reps
#' @param return.sim keep sims? see meteR
#' @param ... addtl pars passed to meteR
#'
#' @return z score
#' @export
#'
#' @examples
#'row_from_modelSteps <- c(1, 1, 1, 8, 10, 10)
#' row_df <- package_timestep(row_from_modelSteps)
#' meteZscore(row_df)
#' @importFrom meteR meteESF sad
meteZscore <- function(packaged_timestep, nrep = 999, return.sim = FALSE, ...) {

  this_esf <- meteR::meteESF(packaged_timestep$row_from_modelSteps, packaged_timestep$Freq)
  x <- meteR::sad(this_esf)


    lik.obs <- meteR:::logLik.meteDist(x)
    state.var <- sum(x$data)
    lik.sim <- c()
    cat("simulating data that conform to state variables: \n")
    for (i in 1:10) {
      cat(sprintf("attempt %s \n", i))
      this.sim <- replicate(100 * nrep, {
        new.dat <- x$r(length(x$data))
        if (abs(sum(new.dat) - state.var) < 0.001 * state.var) {
          return(NA)
        }
        else {
          return(sum(x$d(new.dat, log = TRUE)))
        }
      })
      lik.sim <- c(lik.sim, this.sim[!is.na(this.sim)])
      if (length(lik.sim) >= nrep)
        break
    }
    if (length(lik.sim) >= nrep) {
      lik.sim <- c(lik.sim[1:nrep], lik.obs)
    }
    else {
      warning(sprintf("%s (not %s as desired) simulated replicates found that match the state variables",
                      length(lik.sim), nrep))
      lik.sim <- c(lik.sim, lik.obs)
    }
    #z <- ((lik.obs - mean(lik.sim))/sd(lik.sim))^2
    z <- ((lik.obs - mean(lik.sim))/sd(lik.sim))
    if (return.sim) {
      lik.sim <- ((lik.sim - mean(lik.sim))/sd(lik.sim))^2
    }
    else {
      lik.sim <- NULL
    }
    return(z)
}
