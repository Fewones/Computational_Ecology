# =============================================================================
# Computational Ecology I – Island Biogeography ("Stochasticity drives diversity")
# Baseline scenario: one island (N = 100), continent with 500 equally frequent
# species, immigration probability m, fully neutral dynamics.
# =============================================================================

#reset
rm(list=ls())
set.seed(123)

#parameter
N         <- 100    # island size
m         <- 0.01   # immigration probability
f         <- 50     # mean fecundity
S_cont    <- 500    # specices variety
n_years   <- 2000   # simulated years
burn_in   <- 1000   # years discarded before computing equilibrium values
n_reps    <- 20     # independent replicate runs
set.seed(42)

#summary statistics
# The island community is stored as an abundance vector of length S_cont:
# abund[i] = number of adults of species i on the island.
species_richness <- function(abund) sum(abund > 0)

gini_simpson <- function(abund) {
  p <- abund / sum(abund)
  1 - sum(p^2)
}

# One year sim
one_year <- function(abund, N, m, f, S_cont) {
  offspring <- rpois(S_cont, abund * f * (1 - m))

  # Immigration: on average N*f*m individuals from the continent, each
  # belonging to one of the S_cont equally frequent species.
  n_imm      <- rpois(1, N * f * m)
  immigrants <- tabulate(sample.int(S_cont, n_imm, replace = TRUE),
                         nbins = S_cont)

  pool <- offspring + immigrants

  #Regulation: exactly N of the offspring survive
  #sampling without replacement from the offspring pool.
  if (sum(pool) <= N)
    return(pool)
  survivors <- sample(rep.int(seq_len(S_cont), pool), N)
  tabulate(survivors, nbins = S_cont)
}

# full sim
simulate_island <- function(N, m, f, S_cont, n_years) {
  # Initial community N individuals sampled at random from the continent
  abund <- tabulate(sample.int(S_cont, N, replace = TRUE), nbins = S_cont)

  out <- data.frame(year = seq_len(n_years),
                    richness = NA_integer_, gini_simpson = NA_real_)
  for (t in seq_len(n_years)) {
    abund <- one_year(abund, N, m, f, S_cont)
    # 4. Store summary statistics at the end of each year
    out$richness[t]     <- species_richness(abund)
    out$gini_simpson[t] <- gini_simpson(abund)
  }
  out
}

#Run replicates
run_scenario <- function(N, m, f = 50, S_cont = 500,
                         n_years = 2000, burn_in = 1000, n_reps = 20) {
  reps <- lapply(seq_len(n_reps), function(r) {
    res <- simulate_island(N, m, f, S_cont, n_years)
    res$rep <- r
    res
  })
  ts  <- do.call(rbind, reps)
  eq  <- ts[ts$year > burn_in, ]
  list(timeseries = ts,
       equilibrium = c(richness     = mean(eq$richness),
                       richness_sd  = sd(eq$richness),
                       gini_simpson = mean(eq$gini_simpson),
                       gs_sd        = sd(eq$gini_simpson)))
}


# DEBUG RUNS TO CHECK FOR FAULTS


# No immigration (m = 0):
dbg0 <- run_scenario(N = 100, m = 0, n_years = 2000, burn_in = 1500, n_reps = 10)
cat("m = 0     :", round(dbg0$equilibrium, 3), "\n")

# drift alone -> one species fixes in the long run

# Baseline: N = 100, m = 0.01  (expected: richness ~33, Gini-Simpson ~0.94)
base <- run_scenario(N = N, m = m, f = f, S_cont = S_cont,
                     n_years = n_years, burn_in = burn_in, n_reps = n_reps)
cat("baseline  :", round(base$equilibrium, 3), "\n")


png("baseline_timeseries.png", width = 1600, height = 700, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
ts <- base$timeseries
plot(NA, xlim = c(1, n_years), ylim = c(0, max(ts$richness)),
     xlab = "Year", ylab = "Species richness",
     main = sprintf("Species richness (N = %d, m = %.2f)", N, m))
for (r in unique(ts$rep)) lines(ts$year[ts$rep == r], ts$richness[ts$rep == r],
                                col = adjustcolor("grey", 0.25))
lines(tapply(ts$richness, ts$year, mean), col = "blue", lwd = 2)
abline(v = burn_in, lty = 2)

plot(NA, xlim = c(1, n_years), ylim = c(0, 1),
     xlab = "Year", ylab = "Gini-Simpson index",
     main = "Gini-Simpson index")
for (r in unique(ts$rep)) lines(ts$year[ts$rep == r], ts$gini_simpson[ts$rep == r],
                                col = adjustcolor("grey", 0.25))
lines(tapply(ts$gini_simpson, ts$year, mean), col = "red", lwd = 2)
abline(v = burn_in, lty = 2)
dev.off()



plot_scenario <- function(res, N, m, file = NULL) {
  # file = NULL -> show in RStudio Plots pane; file = "x.png" -> save as PNG
  ts <- res$timeseries; eq <- res$equilibrium
  if (!is.null(file)) png(file, width = 1600, height = 700, res = 150)
  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
  plot(NA, xlim = c(1, max(ts$year)), ylim = c(0, max(ts$richness)),
       xlab = "Year", ylab = "Species richness",
       main = sprintf("Species richness (N = %d, m = %.2f)", N, m))
  for (r in unique(ts$rep)) lines(ts$year[ts$rep == r], ts$richness[ts$rep == r],
                                  col = adjustcolor("grey", 0.35))
  lines(tapply(ts$richness, ts$year, mean), col = "blue", lwd = 2)
  abline(v = burn_in, lty = 2)
  abline(h = eq["richness"], col = "blue", lty = 3)
  legend("topright", bty = "n", legend = sprintf("equilibrium mean = %.1f", eq["richness"]))

  plot(NA, xlim = c(1, max(ts$year)), ylim = c(0, 1),
       xlab = "Year", ylab = "Gini-Simpson index",
       main = sprintf("Gini-Simpson index (N = %d, m = %.2f)", N, m))
  for (r in unique(ts$rep)) lines(ts$year[ts$rep == r], ts$gini_simpson[ts$rep == r],
                                  col = adjustcolor("grey", 0.35))
  lines(tapply(ts$gini_simpson, ts$year, mean), col = "red", lwd = 2)
  abline(v = burn_in, lty = 2)
  abline(h = eq["gini_simpson"], col = "red", lty = 3)
  legend("bottomright", bty = "n", legend = sprintf("equilibrium mean = %.3f", eq["gini_simpson"]))
  if (!is.null(file)) dev.off()
}

base <- run_scenario(N = N, m = m, f = f, S_cont = S_cont,
                     n_years = n_years, burn_in = burn_in, n_reps = n_reps)
cat("m = 0.01:", round(base$equilibrium, 3), "\n")
plot_scenario(base, N, m)                                    # show on screen   # save as PNG

alt <- run_scenario(N = N, m = 0.08, f = f, S_cont = S_cont,
                    n_years = n_years, burn_in = burn_in, n_reps = n_reps)
cat("m = 0.08:", round(alt$equilibrium, 3), "\n")
plot_scenario(alt, N, 0.08)
plot_scenario(alt, N, 0.08, "baseline_timeseries_m0.08.png")