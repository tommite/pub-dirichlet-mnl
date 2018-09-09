rum.fullsample <- readRDS('rum.fullsample.rds')
dir.fullsample <- readRDS('dir.fullsample.rds')

mnl.fullsample.w <- coeff.to.w(rum.fullsample$mnl$coefficients)
rpl.fullsample.w <- coeff.to.w(rum.fullsample$rpl$coefficients[1:3])
dir.fullsample.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)
