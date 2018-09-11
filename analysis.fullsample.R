library(plyr)
library(smaa)
library(DirichletReg)

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

df.w$Y <- DR_data(df.w[,c("pfs","mod","sev")]) # Add dirichlet response variable to the weight dataset
res.dir <- DirichReg(Y~1,df.w) # Fit null model (model without any covariates) to the weight data

rum.fullsample <- list(mnl=res.mnl, rpl=res.rpl, dir=exp(res.dir$coefficients))
dir.fullsample <- DirichReg(Y~1,df.w)

mnl.fullsample.w <- coeff.to.w(rum.fullsample$mnl$coefficients[1:3])
rpl.fullsample.w <- coeff.to.w(rum.fullsample$rpl$coefficients[1:3])
dir.fullsample.w <- exp(res.dir$coefficients) / sum(exp(res.dir$coefficients))

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
