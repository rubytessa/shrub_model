# Description: test of "plant" package by D. Falster et al. 
# Author: Ruby An
# Date: 2022-03-14


library(plant)

# strategy
s <- FF16_Strategy()
str(s)

# individual
ind = FF16_Individual()

# grow plant 
env <- FF16_fixed_environment(1.0)
times <- seq(0,50, length.out = 101)
result <- grow_individual_to_time(ind, times, env)

## viz 
head(result$state)
plot(times, result$state[,"height"])

# grow patch 
params <- scm_base_parameters("FF16")
params$disturbance_mean_interval <- 30
patch <- expand_parameters(trait_matrix(0.0825, "lma"), params, mutant = FALSE)

result <- run_scm_collect(patch)

## viz
t <- result$time
h <- result$species[[1]]["height", ,]
matplot(t, h, lty =1, col = make_transparent("black", 0.25), type = "l",
        las = 1, xlab ="Time (years)", ylab="Height (m)")

## what environment do plants experience?
env <- result$env[[99]]
plot(env)

## what's the forest floor environment like? 
env_min <- lapply(result$env, function(e) e[1, "canopy_openness"])
plot(t, env_min)

## what's the abundance of plants at different heights over time? 
# Relativise the log densities onto (-4, max)
rel <- function(x, xmin = -4) {
  x[x < xmin] <- xmin
  xmax <- max(x, na.rm=TRUE)
  (x - xmin) / (xmax - xmin)
}

d <- result$species[[1]]["log_density", , ]

rd <- rel(d)

# R doesn't seem to offer a way to plot lines that vary in colour, so
# this is quite roundabout using `segments`, shaded by the density at
# the first part of the line segment:
n <- length(t)
x <- matrix(rep(t, ncol(h)), nrow(h))
col <- matrix(make_transparent("black", rd), nrow(d))

plot(NA, xlim=range(t), ylim=range(h, na.rm=TRUE),
     las=1, xlab="Time (years)", ylab="Cohort height (m)")
segments(x[-1, ], h[-1, ], x[-n, ], h[-n, ], col=col[-n, ], lend="butt")
