library(plyr)
library(smaa)
library(DirichletReg)

source('load.dce.R')
source('dirichlet-cvm.R')

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


## Perform bootstrapping to obtain confidence intervals for the sample mean of the weights
bootstrapWeights <- function(weight.data,n.samples) { 
  weights <- c()
  for (i in 1:n.samples) {
    weights <- rbind(weights,colMeans(weight.data[sample(1:nrow(weight.data),nrow(weight.data),replace=T),]))
  }
  weights
}
  
weight.data <- df.w[,c("pfs","mod","sev")]
bootstrap.samples <- bootstrapWeights(weight.data,1e4)
round(apply(bootstrap.samples,MARGIN=2,quantile,probs=c(0.025,0.975)),2) # Bootstap 95% confidence intervals

## Fit models to full-sample data
df.w$Y <- DR_data(df.w[,c("pfs","mod","sev")]) # Add dirichlet response variable to the weight dataset
res.dir <- DirichReg(Y~1,df.w) # Fit null model (model without any covariates) to the weight data

rum.fullsample <- list(mnl=res.mnl, rpl=res.rpl, dir=exp(res.dir$coefficients))
dir.fullsample <- list(alpha=as.numeric(exp(res.dir$coefficients)),
                       sum=sum(exp(as.numeric(res.dir$coefficients))))

mnl.fullsample.w <- coeff.to.w(rum.fullsample$mnl$coefficients[1:3])
rpl.fullsample.w <- coeff.to.w(rum.fullsample$rpl$coefficients[1:3])
dir.fullsample.w <- exp(res.dir$coefficients) / sum(exp(res.dir$coefficients))

## Calculate Adjusted R2's
adj.r2 <- function(res) {
    LL0 <- (nrow(design.matrix) / 2) * log(0.5)
    1 - ((as.numeric(res$logLik) - length(res$coefficients)) / LL0)
}
mnl.adj.r2 <- adj.r2(rum.fullsample$mnl)
mxl.adj.r2 <- adj.r2(rum.fullsample$rpl)

## Calculate SE's for dirichlet
dir.cvm <- calc.covm(dir.fullsample)
dir.se <- sqrt(diag(dir.cvm))

cat("=== MNL model ===\n")
print(summary(res.mnl))
cat("MNL adjusted McFadden's R2: ", round(mnl.adj.r2, 2), '\n')
cat("MNL normalized weights: ", round(mnl.fullsample.w, 2), '\n')
cat("=============\n")

cat("=== MXL model ===\n")
print(summary(res.rpl))
#cat("MXL adjusted McFadden's R2: ", round(mxl.adj.r2, 2), '\n')
cat("MXL normalized weights: ", round(rpl.fullsample.w, 2), '\n')
cat("=============\n")

cat("=== DIR model ===\n")
print(dir.fullsample)
cat("DIR normalized weights: ", round(dir.fullsample.w, 2), '\n')
cat("DIR parameters (alpha): ", round(dir.fullsample.w, 2), '\n')
cat("DIR standard errors: ", round(res.dir$se), 3), '\n')
cat("DIR confidence intervals:\n")
cat(round(exp(res.dir$coefficients - (res.dir$se * 1.96)), 3), '\n')
cat(round(exp(res.dir$coefficients + (res.dir$se * 1.96)), 3), '\n')
cat("=============\n")

## Save results to use in other analyses
saveRDS(file='rum.fullsample.rds', rum.fullsample)
saveRDS(file='dir.fullsample.rds', dir.fullsample)


