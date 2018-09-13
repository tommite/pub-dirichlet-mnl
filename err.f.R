library(evd)
FINAL.SCALE <- 0.3

## Centralized error function - EV-1
dce.err.f <- function(scale=FINAL.SCALE) {rgumbel(1, loc=0, scale=scale)}
