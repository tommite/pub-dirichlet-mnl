library(plyr)
library(smaa)
library(MCMCprecision)

source('load.dce.R')

memory.limit(8000) #8 gb should be enough

scale.seq <- seq(from=0.1, to=1.0, by=0.1)

res <- llply(scale.seq, function(scale) {
    rlply(100, {

        design.matrix <- make.design.matrix(scale)

        mdata <- mlogit.data(design.matrix, choice='choice',
                             ch.id='idx',
                             id.var='url',
                             shape='long',
                             alt.var='alt')

        mlogit(choice ~ 0 + PFS + mod + sev,
               data=mdata)
        })
}, .progress='text')

## Save raw results for possibly changing plot later on
saveRDS(res, file='results-mnl-scale.rds')

err.res <- as.data.frame(laply(res, function(scale.res) {
    norm.ws <- laply(scale.res, function(res) {coeff.to.w(res$coefficients[1:3])})
    aaply(norm.ws, 1, function(x) eucl.dist(x, as.vector(colMeans(resp.w[,2:4]))))
}))
err.res$scale <- scale.seq

df.molten.mse <- melt(err.res, measure.vars=1:20)
df.molten.mse$scale <- as.factor(err.res$scale)

pdf('mnl-errors-per-scale.pdf', width=15, height=10)
ggplot(df.molten.mse, aes(x=scale, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20) +
    ylab('Euclidean distance') + xlab(expression(Scale~(beta))) +
    plot.theme + scale_colour_economist()
#    ggtitle("Full-sample MNL distance of sample mean per Gumbel scale")
dev.off()




