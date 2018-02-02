## Calculate dirichlet covariance
#' From https://en.wikipedia.org/wiki/Dirichlet_distribution
##
calc.covm <- function(dir.pars, names=c('pfs', 'mod', 'sev'), n=length(names)) {
    covm <- matrix(0, nrow=n, ncol=n)
    alpha0 <- sum(dir.pars$alpha)
    for (i in 1:n) {
        alphai <- dir.pars$alpha[i]
        for (j in 1:n) {
            alphaj <- dir.pars$alpha[j]
            if (i == j) {
                covm[i, j] <- (alphai * (alpha0 - alphai)) /
                    (alpha0^2 * (alpha0 + 1))
            }
            else {
                covm[i, j] <- (-alphai * alphaj) /
                    (alpha0^2 * (alpha0 + 1))
            }
        }
    }
    colnames(covm) <- names
    rownames(covm) <- names
    covm
}
