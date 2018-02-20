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
  log(max(prob.selected,(1/N)^2)) # Log-likelihood contribution
  
}

loglik.survey <- function(alpha,meas,selected,question.id,N=1e2) {
  
  set.seed(sum(alpha))
  
  loglik.sum <- 0
  for (question in unique(question.id)) {
    loglik.sum <- loglik.sum + loglik.question(alpha,meas[question.id==question,],selected[question.id==question],N)
  }
  
  loglik.sum
  
}

loglik.sample <- function(alpha,meas,selected,question.id,subject.id,N=1e2) {
  
  set.seed(sum(alpha))
  n.subjects <- length(unique(subject.id))
  
  loglik.sum <- 0
  for (subject in unique(subject.id)) {
    weights <-rdirichlet(N,alpha)
    total.values <- weights %*% t(meas)
    prob.selected <- rep(1,N)
    for (question in unique(question.id[subject.id==subject])) {
      prob.selected <- prob.selected*as.numeric(aaply(total.values[,question.id==question & subject.id==subject],1,which.max)==which.max(selected[question.id==question & subject.id==subject]))
    }
    loglik.sum <- loglik.sum + log(max(mean(prob.selected),(1/N)^2)) # Log-likelihood contribution
  }
  
  loglik.sum
  
}


###############

PVFs <- list(pfs=partialValue(90,50),mod=partialValue(45,85),sev=partialValue(20,80))

source("load.dce.R")
df <- df[df$url %in% unique(df$url)[1:10],]
#df <- df[df$url=="00ncn5es",]

meas.full <- data.frame(pfs=sapply(df$level.PFS,PVFs[["pfs"]]),mod=sapply(df$level.mod,PVFs[["mod"]]),sev=sapply(df$level.sev,PVFs[["sev"]])) 
selected.full <- df$selected.by.subject
question.id <- df$idx
subject.id <- df$url

alpha <- c(1,1,1)

loglik.survey(alpha,meas.full,selected.full,question.id)

constrOptim(theta=alpha,f=loglik.survey,grad=NULL,ui=diag(3),ci=rep(0,3),meas=meas.full,selected=selected.full,question.id=question.id,
            control=list(fnscale=-1,trace=1))

loglik.sample(alpha,meas.full,selected.full,question.id,df$url)

constrOptim(theta=alpha,f=loglik.sample,grad=NULL,ui=diag(3),ci=rep(0,3),meas=meas.full,selected=selected.full,question.id=question.id,subject.id=subject.id,
            control=list(fnscale=-1,trace=1))


