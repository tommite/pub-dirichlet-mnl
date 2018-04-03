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

##
#' Constructs a normalized weight vector from DCE coefficients
##
coeff.to.w <- function(b) {
    rng.sizes <- aaply(ranges, 1, diff)
    w <- b * rng.sizes
    abs(w / sum(abs(w)))
}

##
#' function for computing the error (euclidean distance)
##
err.f <- function(x, y) {
    stopifnot(length(x) == length(y))
    sqrt(sum((x-y)^2))
}

## Binary MNL choice probability function
ch.prob <- function(u1, u2) {
    exp(u1) / (exp(u1) + exp(u2))
}

## Fit models for the maximum possible data set
rum.fullsample <- simulate.dce(n.questions=16, n.respondents=560)
dir.fullsample <- dirichlet.mle(df.w[,c('pfs', 'mod', 'sev')])
dir.fullsample.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)

n.dir.samples <- 1E3

## Error handling routine to re-do the simulation in case of error
## due to singular design matrix (randomly bad set of questions)
error.catch.simulate <- function(n.questions=6, n.respondents=50, n.simul=20) {
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
    stopifnot(n.respondents > 0 && n.questions > 0 && n.questions < (nrow(design.matrix)/2))

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
res.vary.n <- llply(seq(from=20, to=550, by=10), error.catch.simulate,
                    n.questions=6, n.simul=20)

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
            c(err.f(x$mnl$coefficients, rum.fullsample$mnl$coefficients),
              err.f(x$rpl$coefficients[1:3], rum.fullsampl$rpl$coefficients[1:3]),
              err.f(dir.w, dir.fullsample.w),
              y$n.questions, y$n.respondents)
        })
        colnames(r) <- c('err.mnl', 'err.rpl', 'err.dir', 'n.quest', 'n.respondents')
        r
    })
}

test.stats.mnl <- test.stats.p(res.vary.n, 'mnl')
test.stats.rpl <- test.stats.p(res.vary.n, 'rpl')
test.n.stats.mse <- test.stats.mse(res.vary.n)

test.q.stats.mnl <- test.stats.p(res.vary.q, 'mnl')
test.q.stats.mse <- test.stats.mse(res.vary.q)

## Plot test stats vary n respondents ##
df.molten.mnl <- melt(as.data.frame(test.stats.mnl),
                      measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.rpl <- melt(as.data.frame(test.stats.rpl),
                      measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.mse <- melt(as.data.frame(test.n.stats.mse),
                      measure.vars=c('err.mnl', 'err.rpl', 'err.dir'))
df.molten.q <- melt(as.data.frame(test.q.stats.mse),
                    measure.vars=c('err.mnl', 'err.rpl', 'err.dir'))

do.plots <- function(df.molten, y.label, ymax=0.01) {
    plots <- dlply(df.molten, 'variable', function(df.plot) {
        df.plot$n.respondents <- factor(df.plot$n.respondents,
                                        labels=unique(df.plot$n.respondents))
        ggplot(df.plot, aes(x=n.respondents, y=value)) +
            geom_boxplot(outlier.colour='red', outlier.shape=20) +
            ylab(y.label) + theme_economist() + scale_colour_economist() +
            ggtitle(unique(df.plot$variable)) + ylim(0, ymax)
    })
    dev.new(width=10, height=8)
    do.call(grid.arrange, c(plots, ncol=1))
}

do.plots.q <- function(df.molten, y.label) {
    plots <- dlply(df.molten, 'n.respondents', function(df.plot) {
        df.plot$n.quest <- factor(df.plot$n.quest,
                                  labels=unique(df.plot$n.quest))
        ggplot(df.plot, aes(x=n.quest, y=value)) +
            geom_boxplot(outlier.colour='red', outlier.shape=20) +
            ylab(y.label) + theme_economist() + scale_colour_economist() +
            ggtitle(paste0('n.respondents=',unique(df.plot$n.respondents))) + ylim(0, 0.02)
    })
    dev.new(width=6, height=10)
    do.call(grid.arrange, c(plots, ncol=1))
}

do.plots(df.molten.mnl, 'p-value', 0.5)
do.plots(df.molten.rpl, 'p-value', 0.5)
do.plots(subset(df.molten.mse, variable %in% c('err.mnl', 'err.dir')), 'err', 0.05)
do.plots.q(subset(df.molten.q, variable=='err.mnl'), 'err')

## Do stats for the abstract ##
res.mse.stats <- ldply(c('mnl', 'rpl', 'dir'), function(model) {
    my.df <- subset(df.molten.q, variable==paste0('err.', model))
    cbind(ddply(my.df, c('n.questions', 'n.respondents'), summarise,
          mean=mean(value^2), sd=sd(value^2),
          upb=mean(value^2)+1.96*sd(value^2),
          lob=mean(value^2)-1.96*sd(value^2)),
          model=model)
})
res.dir.stats <- ddply(subset(df.molten.mse, variable=='err.dir'), 'n.respondents', summarise,
                       mean=mean(value^2), sd=sd(value^2),
                       upb=mean(value^2)+1.96*sd(value^2),
                       lob=mean(value^2)-1.96*sd(value^2))

## P-value stats
res.p.mnl <- ddply(subset(df.molten.mnl, variable=='mod.p'),
                   c('n.quest', 'n.respondents'), summarise,
                   meanmodp=mean(value),
                   modp02=sum(value>0.02))

res.p.rpl <- ddply(subset(df.molten.rpl, variable=='mod.p'),
                   c('n.quest', 'n.respondents'), summarise,
                   meanmodp=mean(value),
                   modp02=sum(value>0.02))

## for MNL/RPL
print(subset(res.mse.stats, model == 'mnl' & n.questions==6 & upb < 0.01)) # which are not ok
print(subset(res.mse.stats, model == 'mnl' & upb >= 0.01)) # which are not ok
print(subset(res.mse.stats, model == 'rpl' & upb < 0.01)) # which are ok

## for DIR
print(subset(res.dir.stats, upb < 0.01))
print(subset(res.dir.stats, upb < 0.005))
print(subset(res.dir.stats, upb < 0.001))

