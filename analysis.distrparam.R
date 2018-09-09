library(plyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(MCMCprecision)
library(devEMF)
source('dirichlet-cvm.R')
source('simulate-cbm.R')
source('pinkfloyd-plot.R')
source('load.dce.R')

memory.limit(size=16000)

n.simul <- 100

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

## Binary MNL choice probability function
ch.prob <- function(u1, u2) {
    exp(u1) / (exp(u1) + exp(u2))
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

plot.theme <- theme_economist(base_size=20)

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
    ylab('p-value') + xlab('Number of respondents') + plot.theme + scale_colour_economist() +
    coord_cartesian(ylim=c(0, 0.2))
##    ggtitle('Moderate AE coefficient significance')
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
                       ylab('Euclidean distance') + xlab('Number of respondents') +
                       plot.theme + scale_colour_economist() +
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

gen.RPL.weights <- function(coeff, n) {
    w <- cbind(rnorm(n, coeff['PFS'], abs(coeff['sd.PFS'])),
               rnorm(n, coeff['mod'], abs(coeff['sd.mod'])),
               rnorm(n, coeff['sev'], abs(coeff['sd.sev'])))
    colnames(w) <- c('pfs', 'mod', 'sev')
    w <- aaply(w, 1, coeff.to.w)
    w
}
w.rpl.distr <- gen.RPL.weights(rum.fullsample$rpl$coefficients, 560)

plot.mnl.surv <- pinkfloyd.fig(w.mnl.distr)
plot.rpl <- pinkfloyd.fig(w.rpl.distr)

w.df <- df.w[,c('pfs', 'mod', 'sev')]
colnames(w.df)[1] <- 'PFS'
plot.w.orig <- pinkfloyd.fig(w.df)

## NEED TO write a function for sampling weights from the survey assuming DIRICHLET parameters
## GEN WEIGHTS
plot.dir <- pinkfloyd.fig(rdirichlet(nrow(w.mnl.distr), dir.fullsample$alpha))

## AND PLOT THAT
pdf('weightsamples.pdf', width=20, height=8)
grid.arrange(plot.w.orig, plot.dir, plot.rpl, plot.mnl.surv, nrow=2)
dev.off()

## DO FIGURE: Add titles for subplots
