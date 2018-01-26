library(plyr)
library(support.CEs)
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

###
#' Simulates a DCE.
#'
#' @param n.questions Number of questions each respondent answers. These are
#' randomly selected from the set of all questions
#' @param n.dce.respondents Number of respondents, sampled randomly from the
#' pool in the original study.
#' @param rpl if random parameter logit is to be used (TRUE) or not (FALSE)
#' @return DCE results as from mlogit
##
simulate.dce <- function(n.questions=6, n.dce.respondents=50, rpl=FALSE) {
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
    if (rpl) {
        mlogit(choice ~ 0 + PFS + mod + sev,
               rpar=c(PFS='n', mod='n', sev='n'),
               data=mdata)
    } else {
        mlogit(choice ~ 0 + PFS + mod + sev,
               data=mdata)
    }
}

## Error handling routine to re-do the simulation in case of error
## due to singular design matrix (randomly bad set of questions)
error.catch.simulate.dce <- function(n.questions=6, n.dce.respondents=50, n.dces=20, rpl) {
    n.errs <- 0
    n.ok <- 0
    resl <- list()
    while(n.ok < n.dces) {
        tryCatch({
            res <- simulate.dce(n.questions, n.dce.respondents, rpl)
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

##
#' Constructs a normalized weight vector from DCE coefficients
##
coeff.to.w <- function(b) {
    w <- b * aaply(ranges, 1, diff)
    w <- abs(w)
    w / sum(w)
}

### MNL ANALYSES ###
## vary number of respondents
res.vary.n <- llply(seq(from=10, to=300, by=10), error.catch.simulate.dce, n.questions=6, n.dces=20, rpl=FALSE)
## vary number of questions
res.vary.q <- llply(seq(from=3, to=nrow(design.nondom)/2, by=1), error.catch.simulate.dce,
                    n.dce.respondents=200, n.dces=20, rpl=FALSE)

### RPL ANALYSES ###
## vary number of respondents
rpl.res.vary.n <- llply(seq(from=10, to=300, by=10), error.catch.simulate.dce, n.questions=6, n.dces=20, rpl=TRUE)
## vary number of questions
rpl.res.vary.q <- llply(seq(from=3, to=nrow(design.nondom)/2, by=1), error.catch.simulate.dce,
                    n.dce.respondents=200, n.dces=20, rpl=TRUE)

test.stats.p <- function(res) {
    ldply(res, function(y) {
        r <- laply(y$results, function(x) {
            b <- x$coefficients
            std.err <- sqrt(diag(solve(-x$hessian)))
            z <- b / std.err
            p <- 2 * (1 - pnorm(abs(z)))
            w <- coeff.to.w(b)
            c(as.vector(p), as.vector(w), y$n.questions, y$n.dce.respondents)
        })
        colnames(r) <- c(paste0(names(y$results[[1]]$coefficients), '.p'),
                         paste0(names(y$results[[1]]$coefficients), '.w'),
                         'n.quest', 'n.respondents')
        r
    })
}

test.stats.vary.n <- test.stats.p(res.vary.n)
test.stats.vary.q <- test.stats.p(res.vary.q)

## Plot test stats vary n respondents ##
df.molten <- melt(as.data.frame(test.stats.vary.n),
                  measure.vars=c('PFS.p', 'mod.p', 'sev.p'))

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

## Plot weights ##
df.molten.w <- melt(as.data.frame(test.stats.vary.n),
                    measure.vars=c('PFS.w', 'mod.w', 'sev.w'))

plots.w <- dlply(df.molten.w, 'variable', function(df.plot) {
    w.name <- substring(unique(df.plot$variable), 0, 3)
    df.plot$n.respondents <- factor(df.plot$n.respondents,
                                    labels=unique(df.plot$n.respondents))
    ggplot(df.plot, aes(x=n.respondents, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20) +
        ylab('weight') + theme_economist() + scale_colour_economist() +
        ggtitle(w.name) +
        geom_hline(yintercept=true.w[w.name], linetype='solid',
                   color='darkblue', size=1)
})
dev.new(width=10, height=6)
grid.arrange(plots.w[[1]], plots.w[[2]], plots.w[[3]], ncol=1)

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
