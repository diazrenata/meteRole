library(roleR)
library(dplyr)
library(ggplot2)
set.seed(3) #jbp
### Setting up a roleExperiment

#### Replicates of the same parameter settings

params1 <- untbParams(individuals_local = 100, individuals_meta = 1000,
                      species_meta = 50,
                      speciation = 0.2,
                      dispersal_prob = 0.1, init_type = 'oceanic_island',
                      niter = 10000, niterTimestep = 100)

paramsList1 <- list(a= params1, b= params1, c=params1, d= params1, e=params1)

trial1 <- roleExperiment(paramsList1)

results1 <- runRole(trial1)

row_from_modelSteps <- results1@modelRuns[[1]]@modelSteps[[10]]@localComm@indSpecies


####

results1@modelRuns[[1]]@modelSteps[[10]]@localComm@indSpecies # from here you can either run meteR on each item


sumStats1 <- getSumStats(results1, funs =  list(rich = richness,
                                                hill_abund= hillAbund,
                                                abund = totalN)) # or generate meteR for sumstats


library(meteR)


### meteR workflow

a_run_result <- results1@modelRuns[[1]]@modelSteps[[50]]@localComm@indSpecies

a_run_df <- table(a_run_result) %>% as.data.frame()

an_esf <- meteESF(spp = a_run_df$a_run_result, abund = a_run_df$Freq)
an_sad <- sad(an_esf)
ll_sad <- logLikZ(an_sad)
mseZ(an_sad)
meteLogLikZ(an_sad, nrep = 999)


### homebrew loglik
meteLogLikZ <- function (x, nrep = 999, return.sim = FALSE, ...) {
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
  z <- ((lik.obs - mean(lik.sim))/sd(lik.sim))^2
  if (return.sim) {
    lik.sim <- ((lik.sim - mean(lik.sim))/sd(lik.sim))^2
  }
  else {
    lik.sim <- NULL
  }
  return(list(z = z, sim = lik.sim))
}
