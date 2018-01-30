library(dplyr)
library(smaa)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(MCMCprecision)
source('load.dce.R')
set.seed(1911)


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

## Estimate a dirichlet distribution from the weight data ##
res <- llply(seq(from=10, to=nrow(df.w), by=20), function(nr.samples) {
    rlply(20, {
        w <- df.w[sample.int(nrow(df.w), nr.samples, replace=TRUE),c('pfs', 'mod', 'sev')]
        dir.pars <- dirichlet.mle(w)
        n <- length(dir.pars$alpha)
        covm <- calc.covm(dir.pars)
        list(alpha=dir.pars$alpha, cov=covm, nr.samples=nr.samples)
    })
})

dir.stats <- ldply(res, function(x) {
    ldply(x, function(y) {
        row <- c(y$alpha, diag(y$cov), y$nr.samples)
        names(row) <- c('a.PFS', 'a.mod', 'a.sev', 'var.PFS', 'var.mod', 'var.sev', 'n')
        row
    })
})

## normalize the weights
dir.stats.norm <- adply(dir.stats, 1, function(row) {
    w.sum <- sum(row[1:3])
    row[1:3] <- row[1:3] / w.sum
    row['w.sum'] <- w.sum
    colnames(row)[1:3] <- c('PFS.w', 'mod.w', 'sev.w')
    row
})

## Melt & plot
df.molten <- melt(dir.stats.norm,
                  measure.vars=c('PFS.w', 'mod.w', 'sev.w'))

plots <- dlply(df.molten, 'variable', function(df.plot) {
    w.name <- substring(unique(df.plot$variable), 0, 3)
    df.plot$n <- factor(df.plot$n,
                        labels=unique(df.plot$n))
    ggplot(df.plot, aes(x=n, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20) +
        ylab('weight') + theme_economist() + scale_colour_economist() +
        ggtitle(unique(df.plot$variable)) +
        geom_hline(yintercept=true.w[w.name], linetype='solid',
                   color='darkblue', size=1)

})
dev.new(width=10, height=6)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol=1)


## ## Visualize - TBD
## theme_set(theme_bw())
## f <- function(v) ddirichlet(v, dir.pars$alpha)
## mesh <- simplex_mesh(.0025) %>% as.data.frame %>% tbl_df
## mesh$f <- mesh %>% apply(1, function(v) f(bary2simp(v)))
## ggplot(mesh, aes(x, y)) +
##     geom_raster(aes(fill = f)) +
##     coord_equal(xlim = c(0,1), ylim = c(0, .85))
