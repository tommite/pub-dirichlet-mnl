library(plyr)
library(support.CEs)
library(MCMCprecision)
library(smaa)
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

n.dce.respondents <- length(unique(df$url))

## Simulate a DCE ##
n.questions <- 6
n.respondents <- 50

q.idx <- sample.int(nrow(design.nondom$alt.1), n.questions)
respondents <- sample(unique(df$url), n.dce.respondents)

qs <- design.nondom[design.nondom$q.nr %in% q.idx,]
resp.w <- subset(df.w, url %in% respondents)[,c('url', 'pfs', 'mod', 'sev')]

simul.dce <- adply(resp.w, 1, function(row) {
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
simul.dce$idx <- paste0(simul.dce$url, '.', simul.dce$question.no)

mdata <- mlogit.data(simul.dce, choice='choice',
                     ch.id='idx',
                     id.var='url',
                     shape='long',
                     alt.var='alt')
res <- mlogit(choice ~ 0 + PFS + mod + sev,
              data=mdata)

