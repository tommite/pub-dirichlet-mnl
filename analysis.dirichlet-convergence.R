source('load.dce.R')
source('load.fullres.sample.R')

n.simul <- 100

test.stats.mse <- function(res, f=eucl.dist) {
    laply(res, function(y) {
        norm.dir <- aaply(y, 1, function(row) {
            row / sum(row)
        })
        r.full <- aaply(norm.dir, 1, f, dir.fullsample.w)
        r.full
    })
}

n.seq <- seq(from=200, to=3000, by=200)

res.vary.n <- llply(n.seq, function(n.respondents) {
    raply(n.simul, {
        dir.w <- rdirichlet(n.respondents, dir.fullsample$alpha)
        dirichlet.mle(dir.w)$alpha
    })
})

t.stats <- as.data.frame(test.stats.mse(res.vary.n))
t.stats$n.respondents <- n.seq

df.molten.mse <- melt(t.stats, measure.vars=1:100)
df.molten.mse$n.respondents <- as.factor(df.molten.mse$n.respondents)

pdf('error-dirichlet-convergence.pdf', width=15, height=10)
ggplot(df.molten.mse, aes(x=n.respondents, y=value)) +
    geom_boxplot(outlier.colour='red', outlier.shape=20) +
    ylab('Euclidean distance') + xlab('Number of respondents') +
    plot.theme + scale_colour_economist() +
    ggtitle("Dirichlet convergence")
dev.off()

