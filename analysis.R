library(plyr)
library(support.CEs)
library(MCMCprecision)
library(smaa)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
source('load.dce.R')

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
ranges[2:3,] <- ranges[2:3,c(2,1)] # inverse mod,sev to get correct preference direction

## Estimate a dirichlet distribution from the weight data ##
dir.pars <- dirichlet.mle(df.w[,c('pfs', 'mod', 'sev')])

###
#' Simulates a DCE.
#'
#' @param n.questions Number of questions each respondent answers. These are
#' randomly selected from the set of all questions
#' @param n.dce.respondents Number of respondents, sampled randomly from the
#' pool in the original study.
#'  @param return DCE results as from mlogit
##
simulate.dce <- function(n.questions=6, n.dce.respondents=50) {
    stopifnot(n.dce.respondents <= length(unique(df$url))) # PRECOND
    stopifnot(n.dce.respondents > 0)

    q.idx <- sample(unique(design.nondom$q.nr), n.questions)
    respondents <- sample(unique(df$url), n.dce.respondents)

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
                         smaa.pvf(data[,'mod'], cutoffs=ranges['mod',], values=c(0,1)),
                         smaa.pvf(data[,'sev'], cutoffs=ranges['sev',], values=c(0,1)))
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
    mlogit(choice ~ 0 + PFS + mod + sev,
                  data=mdata)
}

## Error handling routine to re-do the simulation in case of error
## due to singular design matrix (randomly bad set of questions)
error.catch.simulate.dce <- function(n.questions=6, n.dce.respondents=50, n.dces=20) {
    n.errs <- 0
    n.ok <- 0
    resl <- list()
    while(n.ok < n.dces) {
        tryCatch({
            res <- simulate.dce(n.questions, n.dce.respondents)
            n.ok <- n.ok + 1
            resl[[n.ok]] <- res
        }, error=function(e) {
            n.errs <<- n.errs + 1
        })
    }
    cat(n.errs, ' errors\n')
    list(n.questions=n.questions, n.dce.respondents=n.dce.respondents,
         n.errors=n.errs, results=resl)
}

## vary number of respondents
res.vary.n <- llply(seq(from=10, to=300, by=10), error.catch.simulate.dce, n.questions=6, n.dces=20)
## vary number of questions
res.vary.q <- llply(seq(from=3, to=nrow(design.nondom)/2, by=1), error.catch.simulate.dce,
                    n.dce.respondents=200, n.dces=20)

test.stats.p <- function(res) {
    ldply(res, function(y) {
        r <- laply(y$results, function(x) {
            b <- x$coefficients
            std.err <- sqrt(diag(solve(-x$hessian)))
            z <- b / std.err
            p <- 2 * (1 - pnorm(abs(z)))
            c(as.vector(p), y$n.questions, y$n.dce.respondents)
        })
        colnames(r) <- c(names(y$results[[1]]$coefficients),
                         'n.quest', 'n.respondents')
        r
    })
}

test.stats.vary.n <- test.stats.p(res.vary.n)
test.stats.vary.q <- test.stats.p(res.vary.q)

## Plot test stats vary n respondents ##
df.molten <- melt(as.data.frame(test.stats.vary.n),
                  measure.vars=c(names(res.vary.n[[1]]$results[[1]]$coefficients)))

plots <- dlply(df.molten, 'variable', function(df.plot) {
    df.plot$n.respondents <- factor(df.plot$n.respondents,
                                    labels=unique(df.plot$n.respondents))
    ggplot(df.plot, aes(x=n.respondents, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20) +
        ylab('p-value') + theme_economist() + scale_colour_economist() +
        ggtitle(unique(df.plot$variable))
})
dev.new(width=10, height=6)
grid.arrange(plots[[1]], plots[[2]], plots[[3]], ncol=1)

## Plot test stats vary n questions
df.molten.q <- melt(as.data.frame(test.stats.vary.q),
                  measure.vars=c(names(res.vary.n[[1]]$results[[1]]$coefficients)))

plots.q <- dlply(df.molten.q, 'variable', function(df.plot) {
    df.plot$n.quest <- factor(df.plot$n.quest,
                              labels=unique(df.plot$n.quest))
    ggplot(df.plot, aes(x=n.quest, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20) +
        ylab('p-value') + theme_economist() + scale_colour_economist() +
        ggtitle(unique(df.plot$variable))
})
dev.new(width=8, height=6)
grid.arrange(plots.q[[1]], plots.q[[2]], plots.q[[3]], ncol=1)

## TODO: plot normalized weights per attribute (boxplots) + means from Douwe's study
## TODO: Plot dirichlet means + concentration parameter estimates
