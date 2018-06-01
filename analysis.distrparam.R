library(plyr)
library(support.CEs)
library(smaa)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(MCMCprecision)
library(devEMF)
source('load.dce.R')
source('dirichlet.R')
source('simulate-cbm.R')
source('plotting.R')

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

n.simul <- 100

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

##
#' Constructs a normalized weight vector from DCE coefficients
##
coeff.to.w <- function(b) {
    rng.sizes <- aaply(ranges, 1, diff)
    w <- b * rng.sizes
    abs(w / sum(abs(w)))
}

eucl.dist <- function(x, y) {
    stopifnot(length(x) == length(y))
    sqrt(sum((x-y)^2))
}

## Binary MNL choice probability function
ch.prob <- function(u1, u2) {
    exp(u1) / (exp(u1) + exp(u2))
}

## Fit models for the maximum possible data set
rum.fullsample <- simulate.dce(n.questions=16, n.respondents=560)
mnl.fullsample.w <- coeff.to.w(rum.fullsample$mnl$coefficients)
dir.fullsample <- dirichlet.mle(df.w[,c('pfs', 'mod', 'sev')])
dir.fullsample.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)

## Calculate R2 for full-sample MNL
l0 <- (560 * 16 * log(0.5))
l1 <- rum.fullsample$mnl$logLik[1]
R2 <- 1 - ((l1 - 3) / l0)

n.dir.samples <- 1E3

## Error handling routine to re-do the simulation in case of error
## due to singular design matrix (randomly bad set of questions)
error.catch.simulate <- function(n.questions=6, n.respondents=50, n.simul=100) {
    n.errs <- 0
    n.ok <- 0
    resl <- list()
    while(n.ok < n.simul) {
        ## Somehow the same seed is being used in different calls (??)
        ## So need to manually set it
        seed <- n.errs * 10000 + n.questions*1000 + n.respondents*100 + n.simul*10 + n.ok
        set.seed(seed)
        tryCatch({
            res <- simulate.dce.distr(n.questions, n.respondents)
            n.ok <- n.ok + 1
            resl[[n.ok]] <- res
        }, error=function(e) {
            n.errs <<- n.errs + 1
            cat('(seed ', seed, '): ', e$message, '\n', sep='')
        })
    }
    ## Simulate from dirichlet
    res.dir <- raply(n.simul, {
      dir.w <- rdirichlet(n.respondents, dir.fullsample$alpha)
      dirichlet.mle(dir.w)$alpha
    })
    cat('[', n.questions, ' questions | ', n.respondents, ' respondents]: ',
        n.errs, ' errors\n', sep='')
    list(n.questions=n.questions, n.respondents=n.respondents,
         n.errors=n.errs, res.dce=resl, res.dir=res.dir)
}

## Add choice probabilities to the design
design.chprob <- ddply(design.nondom, ~q.nr, function(q) {
    utils <- aaply(q, 1, function(r) {
        a1 <- t(as.matrix(r[,c('PFS', 'mod', 'sev')]))
        a2 <- as.matrix(rum.fullsample$mnl$coefficients)
        t(a1) %*% a2
    }, .expand=FALSE)
    probs <- c(ch.prob(utils[1], utils[2]), ch.prob(utils[2], utils[1]))
    ## calculate dirichlet probabilities
    dir.w <- rdirichlet(n.dir.samples, dir.fullsample$alpha)

    pvs <- cbind(smaa.pvf(q[,'PFS'], cutoffs=ranges['PFS',], values=c(0,1)),
                 smaa.pvf(q[,'mod'], cutoffs=ranges['mod',], values=c(1,0)),
                 smaa.pvf(q[,'sev'], cutoffs=ranges['sev',], values=c(1,0)))
    dir.vals <- dir.w %*% t(pvs)
    a1.dir.prob <- sum(aaply(dir.vals, 1, which.max) == 1) / n.dir.samples

    cbind(q, u.mnl=utils, p.mnl=probs, p.dir=c(a1.dir.prob, 1-a1.dir.prob))
})

simulate.dce.distr <- function(n.questions=6, n.respondents=50) {
    stopifnot(n.respondents > 0 && n.questions > 0 && n.questions < (nrow(design.nondom)/2))

    q.idx.mat <- raply(n.respondents, sample(unique(design.nondom$q.nr), n.questions, replace=FALSE))
    q.idx <- as.vector(q.idx.mat)
    qs <- ldply(q.idx, function(x) {subset(design.chprob, q.nr ==  x)})

    design.matrix <- ldply(seq(1, nrow(qs), by=2), function(id) {
        rows <- qs[id:(id+1),]
        rows$idx <- paste0(id, '.', rows$q.nr)

        rnd <- runif(1)
        if (rnd <= rows[1,'p.mnl']) {
            rows$choice <- c(1, 0)
        } else {
            rows$choice <- c(0, 1)
        }
        rows
    })

    mdata <- mlogit.data(design.matrix, choice='choice',
                         ch.id='idx',
                         shape='long',
                         alt.var='alt')

    mlogit(choice ~ 0 + PFS + mod + sev, data=mdata)
}

## vary number of respondents
res.vary.n <- llply(seq(from=20, to=540, by=20), error.catch.simulate,
                    n.questions=6, n.simul=n.simul)

test.stats.p <- function(res) {
    ldply(res, function(y) {
        r <- laply(y$res.dce, function(x) {
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
        colnames(r) <- c(paste0(names(y$res.dce[[1]]$coefficients), '.p'),
                         paste0(names(y$res.dce[[1]]$coefficients), '.w'),
                         'n.quest', 'n.respondents')
        r
    })
}

test.stats.mse <- function(res, f=eucl.dist) {
    ldply(res, function(y) {
        r <- laply(y$res.dce, function(x) {
          f(coeff.to.w(x$coefficients),
            coeff.to.w(rum.fullsample$mnl$coefficients))
        })
        norm.dir <- aaply(y$res.dir, 1, function(row) {
          row / sum(row)})
        r.full <- cbind(r, aaply(norm.dir, 1, f, dir.fullsample.w),
                   y$n.questions, y$n.respondents)
        colnames(r.full) <- c('err.mnl', 'err.dir', 'n.quest', 'n.respondents')
        r.full
    })
}

test.p.stats.mnl <- test.stats.p(res.vary.n)
test.stats.mse <- test.stats.mse(res.vary.n)

## Plot test stats vary n respondents ##
df.molten.p <- melt(as.data.frame(test.p.stats.mnl),
                      measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.mse <- melt(as.data.frame(test.stats.mse),
                      measure.vars=c('err.mnl', 'err.dir'))

pdf('mod-ae-significance.pdf', width=15, height=8)
df.plot <- subset(df.molten.p, n.respondents <= 560 & n.respondents >=100 & variable == 'mod.p')
df.plot$n.respondents <- factor(df.plot$n.respondents,
                                labels=unique(df.plot$n.respondents))
p <- ggplot(df.plot, aes(x=n.respondents, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20) +
    ylab('p-value') + theme_economist() + scale_colour_economist() +
    ggtitle('Moderate AEs coefficient significance') + coord_cartesian(ylim=c(0, 0.2))
p + scale_y_continuous(breaks = sort(c(ggplot_build(p)$layout$panel_ranges[[1]]$y.major_source, 0.05)))
dev.off()

pdf('error-eucl.pdf', width=15, height=10)
## Revalue for having correct subplot titles ##
df.molten.mse$variable <- revalue(df.molten.mse$variable, c('err.mnl'='MNL', 'err.dir'='Dirichlet'))
plots <- dlply(subset(df.molten.mse, n.respondents <= 560), 'variable',
               function(df.plot) {
                   cut.off <- 0.15
                   df.plot[df.plot$value > cut.off, 'value'] <- cut.off
                   df.plot$n.respondents <- factor(df.plot$n.respondents,
                                                   labels=unique(df.plot$n.respondents))
                   ggplot(df.plot, aes(x=n.respondents, y=value)) +
                       geom_boxplot(outlier.colour='red', outlier.shape=20) +
                       ylab('Euclidean distance') + theme_economist() + scale_colour_economist() +
                       ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off))
})
do.call(grid.arrange, c(plots, ncol=1))
dev.off()

## Simulation #2 as per Douwe's email ##
simulate.dce.dir <- function(n.questions=6, n.respondents=50) {
    stopifnot(n.respondents > 0 && n.questions > 0 && n.questions < (nrow(design.nondom)/2))

    w.dir <- rdirichlet(n.respondents, dir.fullsample$alpha)
    q.idx.mat <- raply(n.respondents, sample(unique(design.nondom$q.nr), n.questions, replace=FALSE))
    qs.ws <- as.data.frame(cbind(w.dir, q.idx.mat))

    design.matrix <- adply(qs.ws, 1, function(row) {
        ws <- row[1:3]
        qs <- row[-(1:3)]
        q.rows <- subset(design.nondom, q.nr %in% qs)
        x <- ldply(unique(q.rows$q.nr), function(q) {
            r <- subset(q.rows, q.nr %in% q)

            data <- as.matrix(r[,c('PFS', 'mod', 'sev')])
            ## convert to partial values
            pvs <- cbind(smaa.pvf(data[,'PFS'], cutoffs=ranges['PFS',], values=c(0,1)),
                         smaa.pvf(data[,'mod'], cutoffs=ranges['mod',], values=c(1,0)),
                         smaa.pvf(data[,'sev'], cutoffs=ranges['sev',], values=c(1,0)))
            colnames(pvs) <- colnames(data)
            vals <- as.numeric(ws) %*% t(pvs)
            r$choice <- if(vals[1] > vals[2]) c(1, 0) else c(0, 1)
            r
        })
        x$id <- rownames(row)
        x
    }, .expand=FALSE)
    design.matrix$idx <- paste0(design.matrix$id, '.', design.matrix$q.nr)

    mdata <- mlogit.data(design.matrix, choice='choice',
                         ch.id='idx',
                         id.var='X1',
                         shape='long',
                         alt.var='alt')

    mlogit(choice ~ 0 + PFS + mod + sev, data=mdata)
}

error.catch.simulate.dce.dir <- function(n.questions=6, n.respondents=50, n.simul=100) {
    n.errs <- 0
    n.ok <- 0
    resl <- list()
    while(n.ok < n.simul) {
        ## Somehow the same seed is being used in different calls (??)
        ## So need to manually set it
        seed <- n.errs * 10000 + n.questions*1000 + n.respondents*100 + n.simul*10 + n.ok
        set.seed(seed)
        tryCatch({
            res <- simulate.dce.dir(n.questions, n.respondents)
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
         n.errors=n.errs, res.dce=resl)
}

res.dce.dir <- llply(seq(from=20, to=540, by=20), error.catch.simulate.dce.dir,
                    n.questions=6, n.simul=n.simul)

test.stats.dce.dir <- ldply(res.dce.dir, function(y) {
    r <- laply(y$res.dce, function(x) {
        eucl.dist(coeff.to.w(x$coefficients),
                  dir.fullsample.w)
    })
    r.full <- cbind(r, y$n.questions, y$n.respondents)
    colnames(r.full) <- c('err.mnl', 'n.quest', 'n.respondents')
    r.full
})

test.stats.dce.dir.p <- test.stats.p(res.dce.dir)

df.molten.dce.dir <- melt(as.data.frame(test.stats.dce.dir),
                          measure.vars=c('err.mnl'))
df.molten.dce.dir.p <- melt(as.data.frame(test.stats.dce.dir.p),
                            measure.vars=c('PFS.p', 'mod.p', 'sev.p'))

df.molten.dce.dir$variable <- revalue(df.molten.dce.dir$variable, c('err.mnl'='MNL (correct preference model Dirichlet)'))
plot.dce.dir <- dlply(subset(df.molten.dce.dir, n.respondents <= 550), 'variable',
               function(df.plot) {
                   cut.off <- 0.2
                   df.plot[df.plot$value > cut.off, 'value'] <- cut.off
                   df.plot$n.respondents <- factor(df.plot$n.respondents,
                                                   labels=unique(df.plot$n.respondents))
                   ggplot(df.plot, aes(x=n.respondents, y=value)) +
                       geom_boxplot(outlier.colour='red', outlier.shape=20) +
                       ylab('Euclidean distance') + theme_economist() + scale_colour_economist() +
                       ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off))
})

pdf('dce-dir.mod-ae-significance.pdf', width=15, height=8)
df.plot <- subset(df.molten.dce.dir.p, n.respondents <= 550 & variable == 'mod.p')
df.plot$n.respondents <- factor(df.plot$n.respondents,
                                labels=unique(df.plot$n.respondents))
p <- ggplot(df.plot, aes(x=n.respondents, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20) +
    ylab('p-value') + theme_economist() + scale_colour_economist() +
    ggtitle('Moderate AEs coefficient significance') + coord_cartesian(ylim=c(0, 0.7))
p + scale_y_continuous(breaks = sort(c(ggplot_build(p)$layout$panel_ranges[[1]]$y.major_source, 0.05)))
dev.off()

## Simulation #1 as per Douwe's email ##
res.dir.dce <- ldply(seq(from=20, to=540, by=20),function(n.respondents) {
    dir.params <- raply(n.simul, {
        survey.results <- gen.MNL.weights.survey(n.respondents, rum.fullsample$mnl$coefficients)
        dirichlet.mle(survey.results)$alpha
    })
    dir.norm <- aaply(dir.params, 1, function(row) {row / sum(row)})
    colnames(dir.norm) <- c('PFS', 'mod', 'sev')
    cbind(dir.norm, 'n.respondents'=n.respondents)
}, .progress='text')

res.dir.dce.err <- adply(res.dir.dce, 1, function(row) {
    cbind(row, err=eucl.dist(row[c('PFS', 'mod', 'sev')], mnl.fullsample.w))
})

df.molten.dir.dce <- melt(as.data.frame(res.dir.dce.err),
                          measure.vars=c('err'))

## Revalue for having correct subplot titles ##
df.molten.dir.dce$variable <- revalue(df.molten.dir.dce$variable, c('err'='Dirichlet (correct preference model MNL)'))
plot.dir.dce <- dlply(subset(df.molten.dir.dce, n.respondents <= 550), 'variable',
               function(df.plot) {
                   cut.off <- 0.2
                   df.plot[df.plot$value > cut.off, 'value'] <- cut.off
                   df.plot$n.respondents <- factor(df.plot$n.respondents,
                                                   labels=unique(df.plot$n.respondents))
                   ggplot(df.plot, aes(x=n.respondents, y=value)) +
                       geom_boxplot(outlier.colour='red', outlier.shape=20) +
                       ylab('Euclidean distance') + theme_economist() + scale_colour_economist() +
                       ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off))
               })
pdf('error-eucl-wrong-pmodel.pdf', width=15, height=10)
grid.arrange(plot.dce.dir[[1]], plot.dir.dce[[1]], ncol=1)
dev.off()

######### PLOT WEIGHT DISTRIBUTIONS ##########
w.mnl.distr <- gen.MNL.weights.survey(nrow(df.w), rum.fullsample$mnl$coefficients)

plot.mnl.surv <- pinkfloyd.fig(w.mnl.distr)
w.df <- df.w[,c('pfs', 'mod', 'sev')]
colnames(w.df)[1] <- 'PFS'
plot.w.orig <- pinkfloyd.fig(w.df)

## NEED TO write a function for sampling weights from the survey assuming DIRICHLET parameters
## GEN WEIGHTS
plot.dir <- pinkfloyd.fig(rdirichlet(nrow(w.mnl.distr), dir.fullsample$alpha))

## AND PLOT THAT
pdf('weightsamples.pdf', width=20, height=8)
grid.arrange(plot.w.orig, plot.dir, plot.mnl.surv, nrow=1)
dev.off()

## DO FIGURE: Add titles for subplots
