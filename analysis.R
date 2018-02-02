library(plyr)
library(support.CEs)
library(smaa)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(MCMCprecision)
source('load.dce.R')
source('dirichlet.R')

set.seed(1911)

## Generate L^ma non-dominated design ##
attribute.names <- list('PFS'=sort(unique(df$level.PFS)),
                        'mod'=sort(unique(df$level.mod)),
                        'sev'=sort(unique(df$level.sev)))
design <- Lma.design(attribute.names=attribute.names, nalternatives=2, nblocks=1)
ok.qs <- laply(rownames(design$alternatives$alt.1), function(q.idx) {
    q1 <- as.numeric(as.matrix(design$alternatives$alt.1[q.idx,c('PFS', 'mod', 'sev')]))
    q2 <- as.numeric(as.matrix(design$alternatives$alt.2[q.idx,c('PFS', 'mod', 'sev')]))
    !((q1[1] >= q2[1] && q1[2] <= q2[2] && q1[3] <= q2[3]) ||
      (q1[1] <= q2[1] && q1[2] >= q2[2] && q1[3] >= q2[3]))
})
a1.qs <- design$alternatives$alt.1[ok.qs,]
a2.qs <- design$alternatives$alt.2[ok.qs,]
a1.qs$q.nr <- 1:nrow(a1.qs)
a2.qs$q.nr <- 1:nrow(a2.qs)
a1.qs$alt <- 'A'
a2.qs$alt <- 'B'
design.nondom <- rbind(a1.qs, a2.qs)[,c('q.nr', 'alt', 'PFS', 'mod', 'sev')]
design.nondom <- design.nondom[order(design.nondom$q.nr),]
design.nondom$PFS <- as.numeric(as.vector(design.nondom$PFS))
design.nondom$mod <- as.numeric(as.vector(design.nondom$mod))
design.nondom$sev <- as.numeric(as.vector(design.nondom$sev))

ranges <- laply(attribute.names, range)
rownames(ranges) <- names(attribute.names)

###
#' Simulates a DCE with the three models.
#'
#' @param n.questions Number of questions each respondent answers. These are
#' randomly selected from the set of all questions
#' @param n.respondents Number of respondents, sampled randomly from the
#' pool in the original study.
#' @return list of results from the 3 models
##
simulate.dce <- function(n.questions=6, n.respondents=50) {
    stopifnot(n.respondents <= length(unique(df$url))) # PRECOND
    stopifnot(n.respondents > 0)

    q.idx <- sample(unique(design.nondom$q.nr), n.questions, replace=FALSE)
    respondents <- sample(unique(df$url), n.respondents, replace=TRUE)

    qs <- design.nondom[design.nondom$q.nr %in% q.idx,]
    resp.w <- subset(df.w, url %in% respondents)[,c('url', 'pfs', 'mod', 'sev')]

    design.matrix <- adply(resp.w, 1, function(row) {
        rows <- qs
        rows$url <- row$url
        rows$question.no <- paste0(rownames(row), '.', rows$q.nr)

        ldply(unique(rows$q.nr), function(q) {
            r <- subset(rows, q.nr %in% q)

            data <- as.matrix(r[,c('PFS', 'mod', 'sev')])
            ## convert to partial values
            pvs <- cbind(smaa.pvf(data[,'PFS'], cutoffs=ranges['PFS',], values=c(0,1)),
                         smaa.pvf(data[,'mod'], cutoffs=ranges['mod',], values=c(1,0)),
                         smaa.pvf(data[,'sev'], cutoffs=ranges['sev',], values=c(1,0)))
            colnames(pvs) <- colnames(data)
            w <- as.matrix(subset(resp.w, url==unique(r$url))[,c('pfs', 'mod', 'sev')])
            vals <- w %*% t(pvs)
            r$choice <- if(vals[1] > vals[2]) c(1, 0) else c(0, 1)
            r
        })
    }, .expand=FALSE)
    design.matrix$idx <- paste0(design.matrix$url, '.', design.matrix$question.no)

    mdata <- mlogit.data(design.matrix, choice='choice',
                         ch.id='idx',
                         id.var='url',
                         shape='long',
                         alt.var='alt')

    res.rpl <- mlogit(choice ~ 0 + PFS + mod + sev,
                      rpar=c(PFS='n', mod='n', sev='n'),
                      data=mdata,
                      panel=TRUE,
                      halton=NA)
    res.mnl <- mlogit(choice ~ 0 + PFS + mod + sev,
                      data=mdata)
    res.dir <- dirichlet.mle(resp.w[,c('pfs', 'mod', 'sev')])

    list(mnl=res.mnl, rpl=res.rpl, dir=res.dir$alpha, respondents=respondents)
}

## Error handling routine to re-do the simulation in case of error
## due to singular design matrix (randomly bad set of questions)
error.catch.simulate.dce <- function(n.questions=6, n.respondents=50, n.simul=20) {
    n.errs <- 0
    n.ok <- 0
    resl <- list()
    while(n.ok < n.simul) {
        ## Somehow the same seed is being used in different calls (??)
        ## So need to manually set it
        seed <- n.errs * 10000 + n.questions*1000 + n.respondents*100 + n.simul*10 + n.ok
        set.seed(seed)
        tryCatch({
            res <- simulate.dce(n.questions, n.respondents)
            n.ok <- n.ok + 1
            resl[[n.ok]] <- res
        }, error=function(e) {
            n.errs <<- n.errs + 1
            cat('(seed ', seed, '): ', e$message, '\n', sep='')
        })
    }
    cat('[', n.questions, ' questions | ', n.respondents, ' respondents]: ',
        n.errs, ' errors\n', sep='')
    list(n.questions=n.questions, n.respondents=n.respondents,
         n.errors=n.errs, results=resl)
}

##
#' Constructs a normalized weight vector from DCE coefficients
##
coeff.to.w <- function(b) {
    rng.sizes <- aaply(ranges, 1, diff)
    w <- b * rng.sizes
    w / sum(rng.sizes)
}

##
#' function for computing the squared error
##
MSE <- function(x, y) {
    stopifnot(length(x) == length(y))
    sum((x - y)^2) / length(x)
}

## Fit models for the maximum possible data set
true.res <- simulate.dce(n.questions=16, n.respondents=560)
dir.fullsample <- dirichlet.mle(df.w[,c('pfs', 'mod', 'sev')])
dir.fullsample.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)

## vary number of respondents
res.vary.n <- llply(seq(from=20, to=500, by=20), error.catch.simulate.dce,
                    n.questions=6, n.simul=20)

## vary number of questions and respondents
f.pars <- expand.grid(n.questions=seq(from=5, to=nrow(design.nondom)/2, by=1),
                      n.respondents=c(50, 100, 200, 300))
res.vary.q <- mlply(f.pars, error.catch.simulate.dce, n.simul=20)

test.stats.p <- function(res, type) {
    ldply(res, function(y) {
        r <- laply(y$results, function(x) {
            x <- x[[type]]
            b <- x$coefficients[1:3]
            w <- coeff.to.w(b)
            p <- c(1, 1, 1)
            tryCatch({
                std.err <- sqrt(diag(solve(-x$hessian[1:3,1:3])))
                z <- b / std.err
                p <- 2 * (1 - pnorm(abs(z)))
            }, error=function(e) { })
            c(as.vector(p), as.vector(w), y$n.questions, y$n.respondents)
        })
        colnames(r) <- c(paste0(names(y$results[[1]]$mnl$coefficients), '.p'),
                         paste0(names(y$results[[1]]$mnl$coefficients), '.w'),
                         'n.quest', 'n.respondents')
        r
    })
}

test.stats.mse <- function(res) {
    ldply(res, function(y) {
        r <- laply(y$results, function(x) {
            dir.w <- x$dir / sum(x$dir)
            c(MSE(x$mnl$coefficients, true.res$mnl$coefficients),
              MSE(x$rpl$coefficients[1:3], true.res$rpl$coefficients[1:3]),
              MSE(dir.w, dir.fullsample.w),
              y$n.questions, y$n.respondents)
        })
        colnames(r) <- c('MSE.mnl', 'MSE.rpl', 'MSE.dir', 'n.quest', 'n.respondents')
        r
    })
}

test.stats.mnl <- test.stats.p(res.vary.n, 'mnl')
test.stats.rpl <- test.stats.p(res.vary.n, 'rpl')
test.stats.mse <- test.stats.mse(res.vary.n)

test.q.stats.mnl <- test.stats.p(res.vary.q, 'mnl')
test.q.stats.mse <- test.stats.mse(res.vary.q)

## Plot test stats vary n respondents ##
df.molten.mnl <- melt(as.data.frame(test.stats.mnl),
                      measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.rpl <- melt(as.data.frame(test.stats.rpl),
                      measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.mse <- melt(as.data.frame(test.stats.mse),
                      measure.vars=c('MSE.mnl', 'MSE.rpl', 'MSE.dir'))
df.molten.q <- melt(as.data.frame(test.q.stats.mse),
                    measure.vars=c('MSE.mnl', 'MSE.rpl', 'MSE.dir'))

do.plots <- function(df.molten, y.label) {
    plots <- dlply(df.molten, 'variable', function(df.plot) {
        df.plot$n.respondents <- factor(df.plot$n.respondents,
                                        labels=unique(df.plot$n.respondents))
        ggplot(df.plot, aes(x=n.respondents, y=value)) +
            geom_boxplot(outlier.colour='red', outlier.shape=20) +
            ylab(y.label) + theme_economist() + scale_colour_economist() +
            ggtitle(unique(df.plot$variable)) + ylim(0, 0.01)
    })
    dev.new(width=10, height=6)
    do.call(grid.arrange, c(plots, ncol=1))
}

do.plots.q <- function(df.molten, y.label) {
    plots <- dlply(df.molten, 'variable', function(df.plot) {
        df.plot$n.quest <- factor(df.plot$n.quest,
                                  labels=unique(df.plot$n.quest))
        ggplot(df.plot, aes(x=n.quest, y=value)) +
            geom_boxplot(outlier.colour='red', outlier.shape=20) +
            ylab(y.label) + theme_economist() + scale_colour_economist() +
            ggtitle(unique(df.plot$variable)) + ylim(0, 0.01)
    })
    dev.new(width=10, height=6)
    do.call(grid.arrange, c(plots, ncol=1))
}

do.plots(df.molten.mnl, 'p-value')
do.plots(df.molten.rpl, 'p-value')
do.plots(subset(df.molten.mse, variable %in% c('MSE.mnl', 'MSE.dir')), 'MSE')
do.plots.q(subset(df.molten.q, variable %in% c('MSE.mnl', 'MSE.rpl')), 'MSE')

