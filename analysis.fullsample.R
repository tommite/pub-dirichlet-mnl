library(plyr)
library(smaa)
library(MCMCprecision)

source('load.dce.R')

memory.limit(8000) #8 gb should be enough

design.matrix <- make.design.matrix()

mdata <- mlogit.data(design.matrix, choice='choice',
                     ch.id='idx',
                     id.var='url',
                     shape='long',
                     alt.var='alt')

res.rpl <- mlogit(choice ~ 0 + PFS + mod + sev,
                  rpar=c(PFS='n', mod='n', sev='n'),
                  data=mdata,
                  panel=TRUE,
                  R=5000,
                  halton=NA)

res.mnl <- mlogit(choice ~ 0 + PFS + mod + sev,
                  data=mdata)

res.dir <- dirichlet.mle(resp.w[,c('pfs', 'mod', 'sev')])

rum.fullsample <- list(mnl=res.mnl, rpl=res.rpl, dir=res.dir$alpha)
dir.fullsample <- dirichlet.mle(df.w[,c('pfs', 'mod', 'sev')])

mnl.fullsample.w <- coeff.to.w(rum.fullsample$mnl$coefficients[1:3])
rpl.fullsample.w <- coeff.to.w(rum.fullsample$rpl$coefficients[1:3])
dir.fullsample.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)

n.dir.samples <- 1E3

## Calculate Adjusted R2's
LL0 <- (nrow(design.matrix) / 2) * log(0.5)
mnl.adj.r2 <- 1 - ((as.numeric(res.mnl$logLik) - length(res.mnl$coefficients)) / LL0)

cat("=== MNL model ===\n")
print(summary(res.mnl))
cat("MNL adjusted McFadden's R2: ", round(mnl.adj.r2, 2), '\n')
cat("MNL normalized weights: ", round(mnl.fullsample.w, 2), '\n')
cat("=============\n")

cat("=== DIR model ===\n")
print(dir.fullsample)
cat("DIR normalized weights: ", round(dir.fullsample.w, 2), '\n')
cat("=============\n")


## Save results to user in other analyses
saveRDS(file='rum.fullsample.rds', rum.fullsample)
saveRDS(file='dir.fullsample.rds', dir.fullsample)
