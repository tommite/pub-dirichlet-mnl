source('load.dce.R')
source('load.fullsample.res.R')

## Binary MNL choice probability function
ch.prob <- function(u1, u2) {
    exp(u1) / (exp(u1) + exp(u2))
}

simulate.dce.distr <- function(n.questions=6, n.respondents=50) {
    stopifnot(n.respondents > 0 && n.questions > 0 && n.questions < (nrow(design.nondom)/2))

    q.idx.mat <- raply(n.respondents, sample(unique(design.nondom$q.nr), n.questions, replace=FALSE))
    q.idx <- as.vector(q.idx.mat)

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

    qs <- ldply(q.idx, function(x) {subset(design.chprob, q.nr ==  x)})

    my.design.matrix <- ldply(seq(1, nrow(qs), by=2), function(id) {
        rows <- qs[id:(id+1),]
        rows$idx <- paste0(id, '.', rows$q.nr)

        if (rows[1,'u.mnl'] + dce.err.f() > rows[2, 'u.mnl'] + dce.err.f()) {
            rows$choice <- c(1, 0)
        } else {
            rows$choice <- c(0, 1)
        }
        rows
    })

    mdata <- mlogit.data(my.design.matrix, choice='choice',
                         ch.id='idx',
                         shape='long',
                         alt.var='alt')

    mlogit(choice ~ 0 + PFS + mod + sev, data=mdata)
}

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
##    cat('[', n.questions, ' questions | ', n.respondents, ' respondents]: ',
##        n.errs, ' errors\n', sep='')
    list(n.questions=n.questions, n.respondents=n.respondents,
         n.errors=n.errs, res.dce=resl, res.dir=res.dir)
}

cat('=== Convergence tests - varying number of respondents ===\n')
## vary number of respondents
res.vary.n <- llply(seq(from=20, to=540, by=20), error.catch.simulate,
                    n.questions=6, n.simul=n.simul, .progress='text')

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
    }, .progress='text')
}

test.stats.mar <- function(res) {
    rng.sizes <- aaply(ranges, 1, diff)

    ldply(res, function(y) {
        mar.dce <- laply(y$res.dce, function(x) {
            mars <- c(x$coefficients[1] / -x$coefficients[2],
                      x$coefficients[1] / -x$coefficients[3])
            names(mars) <- c('MAR.mod', 'MAR.sev')
            mars
        })
        mar.dir <- cbind((y$res.dir[,1] / y$res.dir[,2]) * (rng.sizes[2] / rng.sizes[1]),
        (y$res.dir[,1] / y$res.dir[,3]) * (rng.sizes[3] / rng.sizes[1]))
        colnames(mar.dir) <- c('MAR.mod', 'MAR.sev')
        r.full <- cbind(mar.dce, mar.dir,
                        y$n.questions, y$n.respondents)
        colnames(r.full) <- c('MAR.dce.mod', 'MAR.dce.sev',
                              'MAR.dir.mod', 'MAR.dir.sev',
                              'n.quest', 'n.respondents')
        r.full
    }, .progress='text')
}

## Save results for possibly re-doing the figures
saveRDS(res.vary.n, 'res.convergence.rds')

test.p.stats.mnl <- test.stats.p(res.vary.n)
test.stats.mse <- test.stats.mse(res.vary.n)
test.stats.mar <- test.stats.mar(res.vary.n)

## Plot test stats vary n respondents ##
df.molten.p <- melt(as.data.frame(test.p.stats.mnl),
                      measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.mse <- melt(as.data.frame(test.stats.mse),
                      measure.vars=c('err.mnl', 'err.dir'))
df.molten.mar <- melt(as.data.frame(test.stats.mar),
                      measure.vars=c('MAR.dce.mod', 'MAR.dce.sev', 'MAR.dir.mod', 'MAR.dir.sev'))

## Display percentages within 0.05
l_ply(seq(from=20, to=540, by=20), function(ss) {
    data <- subset(test.stats.mse, n.respondents == ss)
    cat('SS ', ss,
        ' MNL: ', sum(data$err.mnl < 0.05) / nrow(data),
        ' DIR: ', sum(data$err.dir < 0.05) / nrow(data), '\n')
})

## Revalue for having correct subplot titles ##
df.molten.mse$variable <- revalue(df.molten.mse$variable, c('err.mnl'='MNL', 'err.dir'='Dirichlet'))
pdf('error-eucl.pdf', width=15, height=10)
do.conv.plot <- function(df.plot, cut.off) {
    df.plot[df.plot$value > cut.off, 'value'] <- cut.off
    df.plot$n.respondents <- factor(df.plot$n.respondents,
                                    labels=unique(df.plot$n.respondents))
    ggplot(df.plot, aes(x=n.respondents, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20, outlier.size=2) +
        ylab('Euclidean distance') + xlab('Number of respondents') +
        plot.theme + scale_colour_economist() +
        ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off)) +
        geom_hline(aes(yintercept=cut.off), color='darkblue', linetype='dashed', size=1)
}
grid.arrange(do.conv.plot(subset(df.molten.mse, n.respondents <= 560 & variable == 'MNL'), 0.15),
             do.conv.plot(subset(df.molten.mse, n.respondents <= 560 & variable == 'Dirichlet'), 0.15), ncol=1)
dev.off()

df.plot <- subset(df.molten.p, n.respondents <= 100 & variable == 'mod.p')
pdf('mod-ae-significance.pdf', width=8, height=8)
df.plot$n.respondents <- factor(df.plot$n.respondents,
                                labels=unique(df.plot$n.respondents))
p <- ggplot(df.plot, aes(x=n.respondents, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20, outlier.size=2) +
    ylab('p-value') + xlab('Number of respondents') + plot.theme + scale_colour_economist() +
    coord_cartesian(ylim=c(0, 0.80))
p + scale_y_continuous(breaks = sort(c(ggplot_build(p)$layout$panel_ranges[[1]]$y.major_source, 0.05)))
dev.off()

## Revalue for having correct subplot titles ##
df.molten.mar$variable <- revalue(df.molten.mar$variable, c('MAR.dce.mod'='MNL - Moderate AEs',
                                                            'MAR.dce.sev'='MNL - Severe AEs',
                                                            'MAR.dir.mod'='Dirichlet - Moderate AEs',
                                                            'MAR.dir.sev'='Dirichlet - Severe AEs'
                                                            ))
do.mar.plot <- function(df.plot, cut.off) {
    df.plot[df.plot$value > cut.off, 'value'] <- cut.off
    df.plot$n.respondents <- factor(df.plot$n.respondents,
                                    labels=unique(df.plot$n.respondents))
    ggplot(df.plot, aes(x=n.respondents, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20, outlier.size=2) +
        ylab('MAR') + xlab('Number of respondents') +
        plot.theme + scale_colour_economist() +
        ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off)) +
        geom_hline(aes(yintercept=cut.off), color='darkblue', linetype='dashed', size=1)
}
pdf('error-mar-moderate.pdf', width=15, height=10)
grid.arrange(do.mar.plot(subset(df.molten.mar,
                                n.respondents <= 560 & variable == 'MNL - Moderate AEs'), 10.0),
             do.mar.plot(subset(df.molten.mar,
                                n.respondents <= 560 & variable == 'Dirichlet - Moderate AEs'), 10.0), ncol=1)
dev.off()
pdf('error-mar-severe.pdf', width=15, height=10)
grid.arrange(do.mar.plot(subset(df.molten.mar,
                                n.respondents <= 560 & variable == 'MNL - Severe AEs'), 5.0),
             do.mar.plot(subset(df.molten.mar,
                                n.respondents <= 560 & variable == 'Dirichlet - Severe AEs'), 5.0), ncol=1)
dev.off()
