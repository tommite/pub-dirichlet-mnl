# Estimate the parameters of the dirichlet distribution from discrete choice data

library("smaa")
library("plyr")

#library("MCMCprecision")
library("DirichletReg")

########### Numerical integration test ##################################

marginal.w1 <- function(w1,alpha=c(1,1,1)) {
  
  loc.w2 <- Vectorize(function(x) {
    
    w <- c(w1,x,1-w1-x)
    
    product <- prod(mapply(function(w,a){w^(a-1)},w,alpha))
    norm.constant <- prod(sapply(alpha,gamma))/gamma(sum(alpha))
    product/norm.constant
    
  })
  
  loc.w3 <- Vectorize(function(x) {
    
    w <- c(w1,1-w1-1,x)
    
    product <- prod(mapply(function(w,a){w^(a-1)},w,alpha))
    norm.constant <- prod(sapply(alpha,gamma))/gamma(sum(alpha))
    product/norm.constant
    
  })
  
  integrate(loc.w2,lower=0,upper=1-w1,subdivisions=5e3,rel.tol=1e-2)$value + integrate(loc.w3,lower=0,upper=1-w1,subdivisions=5e3,rel.tol=1e-2)$value
  
}

integrate(Vectorize(marginal.w1),lower=0,upper=1,alpha=rep(1,3))

###########################

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


loglik.question <- function(alpha,meas,selected,N) {
  
  weights <- rdirichlet(N,alpha)
  
  total.values <- weights %*% t(meas)
  prob.selected <- sum(aaply(total.values,1,which.max)==which.max(selected))/N
  log(max(prob.selected,1e-30)) # Log-likelihood contribution
  
}

# Fixed effects model
loglik.survey <- function(alpha,meas,selected,question.id,R=round(length(unique(question.id))^(5/9))) {
  
  N <- length(unique(question.id))
  
  set.seed(prod(alpha))
  weights <- rdirichlet(R*N,alpha)
  
  loglik.sum <- 0
  weights.index <- 0
  
  for (question in unique(question.id)) {
    total.values <- weights[1:R+weights.index,] %*% t(meas[question.id==question,])
    prob.selected <- sum(aaply(total.values,1,which.max)==which.max(selected[question.id==question]))/R
    loglik.sum <- loglik.sum + log(max(prob.selected,1e-30)) 
    weights.index <- weights.index + R
  }
  
  loglik.sum
  
}

loglik.survey.2 <- function(alpha,meas,selected,question.id,R=round(length(unique(question.id))^(5/9))) {
  
  set.seed(prod(alpha))
  
  loglik.sum <- 0
  for (question in unique(question.id)) {
    loglik.sum <- loglik.sum + loglik.question(alpha,meas[question.id==question,],selected[question.id==question],R)
  }
  
  loglik.sum
  
}

loglik.survey.3 <- function(alpha,meas,selected,question.id,R=round(length(unique(question.id))^(5/9))) {
  
  set.seed(prod(alpha))
  sum(sapply(unique(question.id),function(x) {loglik.question(alpha,meas[question.id==x,],selected[question.id==x],R)}))
  
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
df <- df[df$url %in% unique(df$url)[1],]
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
  
loglik.survey(alpha,meas,selected,question.id)
loglik.survey.2(alpha,meas,selected,question.id)
loglik.survey.3(alpha,meas,selected,question.id)

constrOptim(theta=alpha,f=loglik.survey.3,grad=NULL,ui=diag(3),ci=rep(0,3),meas=meas,selected=selected,question.id=question.id,
            control=list(fnscale=-1,trace=1))


