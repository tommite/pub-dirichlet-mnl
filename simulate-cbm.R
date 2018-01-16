### Simulate the outcome from a choice-based matching experiment ###

library(plyr)
library(smaa)
library(ggplot2)
library(reshape2)
library(MCMCprecision)
library(hitandrun)
library(gridExtra)
library(ggthemes)

# Simulates a choice-based matching procedure
# Assumption: additive value function with linear partial value functions
#
#' @param w respondent's true weight vector
#' @param n.steps number of bisection steps in the choice-based matching procedure (per pairsie comparison)
#' @param har.samples number of har samples based on which the respondent's weight vector is estimated 

simulate.cbm <- function(w,n.steps=2,har.samples=1e5) {
  
  n.crit <- length(w)
  ranking <- order(w,decreasing=T) # Ordinal ranking of the criteria weights
  
  har.constr <- simplexConstraints(n.crit) # Initial constraint set 
  for (index in 1:(n.crit-1)) { 
    
    exact.ratio <- w[ranking[index]]/w[ranking[index+1]] # Exact ratio contraint for the current pairwise comparison
    bounds <- seq(0,1,length.out=(2^n.steps + 1)) # Possible lower- and upper bounds for the ratio contraints for the given number of bisecion steps
    
    # Ratio bound constraint values resulting from the choice-based matching procedure
    upper.bound <- 1/bounds[sum(bounds < 1/exact.ratio)]
    lower.bound <- 1/bounds[sum(bounds < 1/exact.ratio) + 1]
    
    # Generate linear constraints and add them to the constraint set
    har.constr <- mergeConstraints(har.constr,lowerRatioConstraint(n.crit,ranking[index],ranking[index+1],lower.bound))
    if (upper.bound!=Inf) {
      har.constr <- mergeConstraints(har.constr,upperRatioConstraint(n.crit,ranking[index],ranking[index+1],upper.bound))
    }
    
  }
    
  # Generate representative weights by applying HAR sampling
  colMeans(hitandrun(har.constr,n.samples=har.samples))
  
}

cbm.survey <- function(n.respondents,dirichlet.alpha,...) {
  
  resp.w <- rdirichlet(n.respondents,dirichlet.alpha)
  weights.cbm <- adply(resp.w,1,simulate.cbm,...,.id=NULL)
  colMeans(weights.cbm)
  
}

cbm.survey(200,c(2.96,0.97,1.79))

### Conduct experiment ###

respondents <- c(10,50,100,250,500) 
n.experiments <- 20
dirichlet.alpha <- c(2.96,0.97,1.79)
true.w <- dirichlet.alpha/sum(dirichlet.alpha)

results.cbm <- c()
for (n.respondents in respondents) {
  print(n.respondents)
  for (e in 1:n.experiments) {
    print(e)
    results.cbm <- rbind(results.cbm,c(cbm.survey(n.respondents,dirichlet.alpha),n.respondents))
  }
}  
colnames(results.cbm) <- c("pfs","mod","sev","n")
results.cbm <- as.data.frame(results.cbm)
results.cbm$n <- as.factor(results.cbm$n)

## Visualize results 
p.pfs <- ggplot(results.cbm) +
geom_boxplot(mapping=aes(x=n,y=pfs)) +
geom_hline(yintercept=true.w[1],linetype='solid',color='darkblue',size=1) +
ylab("Weight") + ggtitle("PFS.w") + theme_economist() + scale_colour_economist() 

p.mod <- ggplot(results.cbm) +
geom_boxplot(mapping=aes(x=n,y=mod)) +
geom_hline(yintercept=true.w[2],linetype='solid',color='darkblue',size=1) +
ylab("Weight") + ggtitle("mod.w") + theme_economist() + scale_colour_economist()

p.sev <- ggplot(results.cbm) +
geom_boxplot(mapping=aes(x=n,y=sev)) +
geom_hline(yintercept=true.w[3],linetype='solid',color='darkblue',size=1) +
ylab("Weight") + ggtitle("sev.w") + theme_economist() + scale_colour_economist()

dev.new(width=8, height=6)
grid.arrange(p.pfs,p.mod,p.sev,ncol=1)

