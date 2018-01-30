library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(MCMCprecision)

source('load.dce.R')
set.seed(1911)

true.w <- colMeans(df.w[,c('pfs', 'mod', 'sev')])

sample.mean <- function(nr.samples) {
  estimated.w <- colMeans(df.w[sample.int(nrow(df.w), nr.samples, replace=TRUE),c('pfs', 'mod', 'sev')])
  c(estimated.w,SE(estimated.w))
}

# function for computing the squared error
SE <- function(estimated.w) {
  sum((estimated.w - true.w)^2)
}

### Conduct experiment ###
samples <- seq(from=10, to=300, by=20)
n.experiments <- 100

true.w <- colMeans(df.w[,c('pfs', 'mod', 'sev')])
results.SM <- c()
for (nr.samples in samples) {
  results.SM <- rbind(results.SM,cbind(t(sapply(rep(nr.samples,n.experiments),sample.mean)),rep(nr.samples,n.experiments)))
}  

colnames(results.SM) <- c("pfs","mod","sev","MSE","n")
results.SM <- as.data.frame(results.SM)
results.SM$n <- as.factor(results.SM$n)

## Visualize results ###
p.MSE <- ggplot(results.SM) +
  geom_boxplot(mapping=aes(x=n,y=MSE)) +
  ylab("squared error") + ggtitle("Squared error") + theme_economist() + scale_colour_economist() 

p.pfs <- ggplot(results.SM) +
  geom_boxplot(mapping=aes(x=n,y=pfs)) +
  geom_hline(yintercept=true.w[1],linetype='solid',color='darkblue',size=1) +
  ylab("Weight") + ggtitle("PFS.w") + theme_economist() + scale_colour_economist() 

p.mod <- ggplot(results.SM) +
  geom_boxplot(mapping=aes(x=n,y=mod)) +
  geom_hline(yintercept=true.w[2],linetype='solid',color='darkblue',size=1) +
  ylab("Weight") + ggtitle("mod.w") + theme_economist() + scale_colour_economist()

p.sev <- ggplot(results.SM) +
  geom_boxplot(mapping=aes(x=n,y=sev)) +
  geom_hline(yintercept=true.w[3],linetype='solid',color='darkblue',size=1) +
  ylab("Weight") + ggtitle("sev.w") + theme_economist() + scale_colour_economist()

dev.new(width=20, height=14)
grid.arrange(p.pfs,p.mod,p.sev,p.MSE,ncol=1)

