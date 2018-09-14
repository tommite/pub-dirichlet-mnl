source('load.dce.R')
source('load.fullsample.res.R')

######### PLOT WEIGHT DISTRIBUTIONS ##########
gen.RPL.weights <- function(coeff, n) {
    w <- cbind(rnorm(n, coeff['PFS'], abs(coeff['sd.PFS'])),
               rnorm(n, coeff['mod'], abs(coeff['sd.mod'])),
               rnorm(n, coeff['sev'], abs(coeff['sd.sev'])))
    colnames(w) <- c('pfs', 'mod', 'sev')
    w <- aaply(w, 1, coeff.to.w)
    w
}
w.rpl.distr <- gen.RPL.weights(rum.fullsample$rpl$coefficients, 560)

plot.rpl <- pinkfloyd.fig(w.rpl.distr, coeff.to.w(rum.fullsample$rpl$coefficients[1:3]))

w.df <- df.w[,c('pfs', 'mod', 'sev')]
colnames(w.df)[1] <- 'PFS'
plot.w.orig <- pinkfloyd.fig(w.df, colMeans(w.df))

dir.w <- dir.fullsample$alpha / sum(dir.fullsample$alpha)
plot.dir <- pinkfloyd.fig(rdirichlet(nrow(w.df), dir.fullsample$alpha), dir.w)

pdf('weightsamples.pdf', width=20, height=8)
grid.arrange(plot.w.orig + ggtitle("Source data"),
             plot.rpl + ggtitle("MXL"),
             plot.dir + ggtitle("Dirichlet"), nrow=1)
dev.off()

## DO FIGURE: Add titles for subplots
