source('load.dce.R')
source('load.fullres.sample.R')

######### PLOT WEIGHT DISTRIBUTIONS ##########
w.mnl.distr <- gen.MNL.weights.survey(nrow(df.w), rum.fullsample$mnl$coefficients)

gen.RPL.weights <- function(coeff, n) {
    w <- cbind(rnorm(n, coeff['PFS'], abs(coeff['sd.PFS'])),
               rnorm(n, coeff['mod'], abs(coeff['sd.mod'])),
               rnorm(n, coeff['sev'], abs(coeff['sd.sev'])))
    colnames(w) <- c('pfs', 'mod', 'sev')
    w <- aaply(w, 1, coeff.to.w)
    w
}
w.rpl.distr <- gen.RPL.weights(rum.fullsample$rpl$coefficients, 560)

plot.mnl.surv <- pinkfloyd.fig(w.mnl.distr)
plot.rpl <- pinkfloyd.fig(w.rpl.distr)

w.df <- df.w[,c('pfs', 'mod', 'sev')]
colnames(w.df)[1] <- 'PFS'
plot.w.orig <- pinkfloyd.fig(w.df)

plot.dir <- pinkfloyd.fig(rdirichlet(nrow(w.mnl.distr), dir.fullsample$alpha))
pdf('weightsamples.pdf', width=20, height=8)
grid.arrange(plot.w.orig, plot.dir, plot.rpl, plot.mnl.surv, nrow=2)
dev.off()

## DO FIGURE: Add titles for subplots
