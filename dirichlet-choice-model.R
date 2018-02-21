# Estimate the parameters of the dirichlet distribution from discrete choice data

library("smaa")
library("plyr")
library("MCMCprecision")

partialValue <- function(best, worst, cutoffs=numeric(), values=numeric()) {
  if (best > worst) {
    # Increasing
    v <- c(0, values, 1)
    y <- c(worst, cutoffs, best)
  } else {
    # Decreasing
    v <- c(1, values, 0)
    y <- c(best, cutoffs, worst)
  }
  function(x) {
    smaa.pvf(x, y, v, outOfBounds="interpolate")
  }
}

loglik.question <- function(alpha,meas,selected,N=1e2) {
  
  weights <- rdirichlet(N,alpha)
  
  total.values <- weights %*% t(meas)
  prob.selected <- sum(aaply(total.values,1,which.max)==which.max(selected))/N
  log(max(prob.selected,1e-12)) # Log-likelihood contribution
  
}

# Fixed effects model
loglik.survey <- function(alpha,meas,selected,question.id,N=1e2) {
  
  set.seed(prod(alpha))
  
  loglik.sum <- 0
  for (question in unique(question.id)) {
    loglik.sum <- loglik.sum + loglik.question(alpha,meas[question.id==question,],selected[question.id==question],N)
  }
  
  loglik.sum
  
}

# Random effects model
loglik.sample.re <- function(theta,meas,selected,question.id,subject.id,n.outer=10,n.inner=10) {
  
  set.seed(prod(theta))
  
  # theta <- c(alpha,C)
  alpha <- theta[1:(length(theta)-1)]
  C <- theta[length(theta)]
  
  loglik.sum <- 0
  for (subject in unique(subject.id)) {
    meas.subject <- meas[subject.id==subject,]
    selected.subject <- selected[subject.id==subject]
    question.id.subject <- question.id[subject.id==subject]
    prob.selected.subject <- c()
    for (i.outer in 1:n.outer) {
      alpha.inner <- C*as.vector(rdirichlet(1,alpha)) 
      prob.selected.iter <- 1
      for (question in unique(question.id.subject)) {
        weights <- rdirichlet(n.inner,alpha.inner)
        total.values <- weights %*% t(meas.subject[question.id.subject==question,])
        prob.selected <-sum(aaply(total.values,1,which.max)==which.max(selected.subject[question.id.subject==question]))/n.inner
        prob.selected.iter <- prob.selected.iter*prob.selected
      }
      prob.selected.subject <- c(prob.selected.subject,prob.selected.iter)
    }
    loglik.sum <- loglik.sum + log(max(mean(prob.selected.subject),1e-12)) # Log-likelihood contribution
  } 
   
  loglik.sum
  
}

###############

PVFs <- list(pfs=partialValue(90,50),mod=partialValue(45,85),sev=partialValue(20,80))

source("load.dce.R")
df <- df[df$url %in% unique(df$url)[1:5],]
#df <- df[df$url=="00ncn5es",]

meas <- data.frame(pfs=sapply(df$level.PFS,PVFs[["pfs"]]),mod=sapply(df$level.mod,PVFs[["mod"]]),sev=sapply(df$level.sev,PVFs[["sev"]])) 
selected <- df$selected.by.subject
question.id <- df$idx
subject.id <- df$url

alpha <- c(2,2,2)
theta <- c(1,1,1,5)

loglik.sample.re(theta,meas,selected,question.id,subject.id,n.outer=10,n.inner=10)
constrOptim(theta=theta,f=loglik.sample.re,grad=NULL,ui=diag(4),ci=rep(0,4),meas=meas,selected=selected,question.id=question.id,subject.id=subject.id,
            control=list(fnscale=-1,trace=1))
  
loglik.survey(alpha,meas.full,selected.full,question.id)
constrOptim(theta=alpha,f=loglik.survey,grad=NULL,ui=diag(3),ci=rep(0,3),meas=meas,selected=selected,question.id=question.id,
            control=list(fnscale=-1,trace=1))


