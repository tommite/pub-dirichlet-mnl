source('load.dce.R')
source('load.fullsample.res.R')

## Simulation #1: Dirichlet, MNL correct preference model ##
cat('=== Convergence tests -- MNL, correct preference model Dirichlet ===\n')
res.dir.dce <- ldply(seq(from=20, to=540, by=20),function(n.respondents) {
    dir.params <- raply(n.simul, {
        survey.results <- gen.MNL.weights.survey(n.respondents, rum.fullsample$mnl$coefficients)
        dirichlet.mle(survey.results)$alpha
    })
    dir.norm <- aaply(dir.params, 1, function(row) {row / sum(row)})
    colnames(dir.norm) <- c('PFS', 'mod', 'sev')
    cbind(dir.norm, 'n.respondents'=n.respondents)
}, .progress='text')

## Simulation #2: MNL, Dirichlet correct preference model
simulate.dce.dir <- function(n.questions=6, n.respondents=50) {
    stopifnot(n.respondents > 0 && n.questions > 0 && n.questions < (nrow(design.nondom)/2))

    w.dir <- rdirichlet(n.respondents, dir.fullsample$alpha)
    q.idx.mat <- raply(n.respondents, sample(unique(design.nondom$q.nr), n.questions, replace=FALSE))
    qs.ws <- as.data.frame(cbind(w.dir, q.idx.mat))

    my.design.matrix <- adply(qs.ws, 1, function(row) {
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
            r$choice <- if(vals[1] > vals[2] )
                            c(1, 0)
                        else c(0, 1)
            r
        })
        x$id <- rownames(row)
        x
    }, .expand=FALSE)
    my.design.matrix$idx <- paste0(my.design.matrix$id, '.', my.design.matrix$q.nr)

    mdata <- mlogit.data(my.design.matrix, choice='choice',
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
cat('=== Convergence tests -- Dirichlet, correct model MNL ===\n')
res.dce.dir <- llply(seq(from=20, to=540, by=20), error.catch.simulate.dce.dir,
                     n.questions=6, n.simul=n.simul, .progress='text')

## Save results for possibly re-doing the figures
saveRDS(res.dce.dir, file='res.incorrect-pref-model.dce.dir.rds')
saveRDS(res.dir.dce, file='res.incorrect-pref-model.dir.dce.rds')

## Compute test statistics ##
res.dir.dce.err <- adply(res.dir.dce, 1, function(row) {
    cbind(row, err=eucl.dist(row[c('PFS', 'mod', 'sev')], mnl.fullsample.w))
}, .progress='text')

test.stats.dce.dir <- ldply(res.dce.dir, function(y) {
    r <- laply(y$res.dce, function(x) {
        eucl.dist(coeff.to.w(x$coefficients),
                  dir.fullsample.w)
    })
    r.full <- cbind(r, y$n.questions, y$n.respondents)
    colnames(r.full) <- c('err.mnl', 'n.quest', 'n.respondents')
    r.full
}, .progress='text')
test.stats.dce.dir.p <- test.stats.p(res.dce.dir)

### Plotting ###
df.molten.dce.dir <- melt(as.data.frame(test.stats.dce.dir),
                          measure.vars=c('err.mnl'))
df.molten.dce.dir.p <- melt(as.data.frame(test.stats.dce.dir.p),
                            measure.vars=c('PFS.p', 'mod.p', 'sev.p'))
df.molten.dce.dir$variable <- revalue(df.molten.dce.dir$variable, c('err.mnl'='MNL (correct preference model Dirichlet)'))
df.molten.dir.dce <- melt(as.data.frame(res.dir.dce.err),
                          measure.vars=c('err'))
## Revalue for having correct subplot titles ##
df.molten.dir.dce$variable <- revalue(df.molten.dir.dce$variable, c('err'='Dirichlet (correct preference model MNL)'))

df.plot <- subset(df.molten.dce.dir.p, n.respondents <= 550 & variable == 'mod.p')
pdf('dce-dir.mod-ae-significance.pdf', width=15, height=8)
df.plot$n.respondents <- factor(df.plot$n.respondents,
                                labels=unique(df.plot$n.respondents))
p <- ggplot(df.plot, aes(x=n.respondents, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20, outlier.size-2) +
    ylab('p-value') + theme_economist() + scale_colour_economist() +
    ggtitle('Moderate AEs coefficient significance') + coord_cartesian(ylim=c(0, 0.7))
p + scale_y_continuous(breaks = sort(c(ggplot_build(p)$layout$panel_ranges[[1]]$y.major_source, 0.05)))
dev.off()

plot.dce.dir <- dlply(subset(df.molten.dce.dir, n.respondents <= 550), 'variable',
               function(df.plot) {
                   cut.off <- 0.2
                   df.plot[df.plot$value > cut.off, 'value'] <- cut.off
                   df.plot$n.respondents <- factor(df.plot$n.respondents,
                                                   labels=unique(df.plot$n.respondents))
                   ggplot(df.plot, aes(x=n.respondents, y=value)) +
                       geom_boxplot(outlier.colour='red', outlier.shape=20, outlier.size=2) +
                       ylab('Euclidean distance') + theme_economist() + scale_colour_economist() +
                       ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off))+
                       xlab('Number of respondents') +
                       geom_hline(aes(yintercept=cut.off), color='darkblue', linetype='dashed', size=1)
})
plot.dir.dce <- dlply(subset(df.molten.dir.dce, n.respondents <= 550), 'variable',
               function(df.plot) {
                   cut.off <- 0.2
                   df.plot[df.plot$value > cut.off, 'value'] <- cut.off
                   df.plot$n.respondents <- factor(df.plot$n.respondents,
                                                   labels=unique(df.plot$n.respondents))
                   ggplot(df.plot, aes(x=n.respondents, y=value)) +
                       geom_boxplot(outlier.colour='red', outlier.shape=20, outlier.size=2) +
                       ylab('Euclidean distance') + theme_economist() + scale_colour_economist() +
                       ggtitle(unique(df.plot$variable)) + scale_y_continuous(limits=c(0, cut.off)) +
                       xlab('Number of respondents') +
                       geom_hline(aes(yintercept=cut.off), color='darkblue', linetype='dashed', size=1)
               })
pdf('error-eucl-wrong-pmodel.pdf', width=15, height=10)
grid.arrange(plot.dce.dir[[1]], plot.dir.dce[[1]], ncol=1)
dev.off()
