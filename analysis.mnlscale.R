library(plyr)
library(smaa)
library(MCMCprecision)
source('load.dce.R')

scale.seq <- seq(from=0.1, to=1.0, by=0.1)

## Run teh simulations, save results to files
llply(scale.seq, function(scale) {
    my.res <- rlply(1000, {
        design.matrix <- make.design.matrix(scale)

        mdata <- mlogit.data(design.matrix, choice='choice',
                             ch.id='idx',
                             id.var='url',
                             shape='long',
                             alt.var='alt')

        mod.res <- mlogit(choice ~ 0 + PFS + mod + sev,
                          data=mdata)
        mod.res$coefficients
    })
    saveRDS(my.res, file=paste0('results-scale-', scale, '.rds'))
}, .progress='text')

## Load the results from files
res <- llply(scale.seq, function(scale) {
    readRDS(paste0('results-scale-', scale, '.rds'))
})

err.res <- as.data.frame(laply(res, function(scale.res) {
    norm.ws <- laply(scale.res, function(xres) {coeff.to.w(xres[1:3])})
    aaply(norm.ws, 1, function(x) eucl.dist(x, as.vector(colMeans(resp.w[,2:4]))))
}, .progress='text'))
err.res$scale <- scale.seq

df.molten.mse <- melt(err.res, measure.vars=1:1000)
df.molten.mse$scale <- as.factor(err.res$scale)

pdf('mnl-errors-per-scale.pdf', width=15, height=10)
ggplot(df.molten.mse, aes(x=scale, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20) +
    ylab('Euclidean distance') + xlab(expression(Scale~(beta))) +
    plot.theme + scale_colour_economist() +
    ggtitle("Full-sample MNL distance of sample mean per Gumbel scale")
dev.off()




