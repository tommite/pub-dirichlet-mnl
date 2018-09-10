source('load.dce.R')
source('load.fullres.sample.R')

## Binary MNL choice probability function
ch.prob <- function(u1, u2) {
    exp(u1) / (exp(u1) + exp(u2))
}

## Add choice probabilities to the design
design.chprob <- ddply(design.nondom, ~q.nr, function(q) {
    utils <- aaply(q, 1, function(r) {
        a1 <- t(as.matrix(r[,c('PFS', 'mod', 'sev')]))
        a2 <- as.matrix(rum.fullsample$mnl$coefficients)
        (t(a1) %*% a2) + dce.err.f() # Add random error
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

    my.design.matrix <- ldply(seq(1, nrow(qs), by=2), function(id) {
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
    cat('[', n.questions, ' questions | ', n.respondents, ' respondents]: ',
        n.errs, ' errors\n', sep='')
    list(n.questions=n.questions, n.respondents=n.respondents,
         n.errors=n.errs, res.dce=resl, res.dir=res.dir)
}

cat('=== Convergence tests - varying number of respondents ===\n')
## vary number of respondents
res.vary.n <- llply(seq(from=20, to=540, by=20), error.catch.simulate,
                    n.questions=6, n.simul=n.simul, .progress='text')

## calculates p-values
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
