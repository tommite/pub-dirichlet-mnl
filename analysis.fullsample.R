library(plyr)
library(smaa)
library(MCMCprecision)

source('load.dce.R')

memory.limit(8000) #8 gb should be enough

resp.w <- df.w[,c('url', 'pfs', 'mod', 'sev')]

design.matrix <- adply(resp.w, 1, function(row) {
    rows <- design.nondom
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
        vals <- as.matrix(row[-1]) %*% t(pvs)
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

res.rpl <- mlogit(choice ~ 1 + PFS + mod + sev,
                  rpar=c(PFS='n', mod='n', sev='n'),
                  data=mdata,
                  panel=TRUE,
                  R=5000,
                  halton=NA)

res.rpl.cov <- mlogit(choice ~ 1 + PFS + mod + sev,
                      rpar=c(PFS='n', mod='n', sev='n'),
                      data=mdata,
                      panel=TRUE,
                      R=5000,
                      halton=NA,
                      correlation=TRUE)

res.mnl <- mlogit(choice ~ 1 + PFS + mod + sev,
                  data=mdata)
res.dir <- dirichlet.mle(resp.w[,c('pfs', 'mod', 'sev')])

rum.fullsample <- list(mnl=res.mnl, rpl=res.rpl, rpl.cov=res.rpl.cov, dir=res.dir$alpha)
dir.fullsample <- dirichlet.mle(df.w[,c('pfs', 'mod', 'sev')])

mnl.fullsample.w <- coeff.to.w(rum.fullsample$mnl$coefficients)
rpl.fullsample.w <- coeff.to.w(rum.fullsample$rpl$coefficients[1:3])
rpl.cov.fullsample.w <- coeff.to.w(rum.fullsample$rpl.cov$coefficients[1:3])
dir.fullsample.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)

## Calculate R2's
#l0 <- (560 * 16 * log(0.5))
#l1 <- rum.fullsample$mnl$logLik[1]
#R2.mnl <- 1 - ((l1 - 3) / l0)

n.dir.samples <- 1E3

## Save results to user in other analyses
saveRDS('rum.fullsample.rds', rum.fullsample)
saveRDS('dir.fullsample.rds', dir.fullsample)
