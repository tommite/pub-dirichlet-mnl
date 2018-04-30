##### Simulate the outcome from a choice-based matching experiment #####

library(smaa)
library(hitandrun)

#### Choice-based matching with MNL population model ####

# Generate MNL choice-probabilities to the questions for the choice-based matching procedure
add.choice.prob <- function(coefficients,questions) {
  questions$xbeta <- as.matrix(questions[,c("PFS","mod","sev")]) %*% coefficients 
  choice.prob <- c()
  for (i in 1:dim(questions)[1]) {
    choice.prob <- c(choice.prob,exp(questions$xbeta[i])/sum(exp(questions$xbeta[questions$q.nr==questions$q.nr[i]])))
  }
  questions$choice.prob <- choice.prob
  questions
}

# Ordinal swing weighting
ordinal.swing.MNL <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,3),alt=c("A","B","C"),PFS=c(90,50,50),mod=c(85,45,85),sev=c(80,80,20))
  question.1 <- add.choice.prob(coefficients,question.1)
  choice.1 <- runif(1)
  if (choice.1<=question.1$choice.prob[1]) { # PFS most important
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(90,90),mod=c(45,85),sev=c(80,20))
    question.2 <- add.choice.prob(coefficients,question.2) 
    if (runif(1)<=question.2$choice.prob[1]) { # mod second most important
      results.ordinal.swing <- list(first="PFS",second="mod",third="sev",prob=question.1$choice.prob[1]*question.2$choice.prob[1])
    } else { # sev second most important
      results.ordinal.swing <- list(first="PFS",second="sev",third="mod",prob=question.1$choice.prob[1]*question.2$choice.prob[2])
    }
  } else {
    if (choice.1<=question.1$choice.prob[1]+question.1$choice.prob[2]) { # mod most important
      question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(90,50),mod=c(45,45),sev=c(80,20))
      question.2 <- add.choice.prob(coefficients,question.2) 
      if (runif(1)<=question.2$choice.prob[1]) { # PFS second most important
        results.ordinal.swing <- list(first="mod",second="PFS",third="sev",prob=question.1$choice.prob[2]*question.2$choice.prob[1])
      } else { # sev second most important
        results.ordinal.swing <- list(first="mod",second="sev",third="PFS",prob=question.1$choice.prob[2]*question.2$choice.prob[2])
      }
    } else { # sev most important
      question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(90,50),mod=c(85,45),sev=c(20,20))
      question.2 <- add.choice.prob(coefficients,question.2) 
      if (runif(1)<=question.2$choice.prob[1]) { # PFS second most important
        results.ordinal.swing <- list(first="sev",second="PFS",third="mod",prob=question.1$choice.prob[3]*question.2$choice.prob[1])
      } else { # modv second most important
        results.ordinal.swing <- list(first="sev",second="mod",third="PFS",prob=question.1$choice.prob[3]*question.2$choice.prob[2])
      }
    }
  }
  
  results.ordinal.swing
  
}

# Choice-based matching to trade-off PFS against sev
cbm.PFS.sev <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,2),alt=c("A","B"),PFS=c(70,50),mod=c(65,65),sev=c(80,20)) # First bisection step
  question.1 <- add.choice.prob(coefficients,question.1)
  
  if(runif(1)<=question.1$choice.prob[1]) { 
    
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(60,50),mod=c(65,65),sev=c(80,20)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) {
      constr <- lowerRatioConstraint(3,1,3,4)
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,1,3,2),upperRatioConstraint(3,1,3,4))
    }
    
  } else { 
    
    question.2 <- data.frame(q.nr=rep(3,2),alt=c("A","B"),PFS=c(80,50),mod=c(65,65),sev=c(80,20)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- mergeConstraints(lowerRatioConstraint(3,1,3,1+1/3),upperRatioConstraint(3,1,3,2))
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,1,3,1),upperRatioConstraint(3,1,3,1+1/3))
    }
    
  }
  
  constr
  
}

# Choice-based matching to trade-off PFS against mod
cbm.PFS.mod <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,2),alt=c("A","B"),PFS=c(70,50),mod=c(85,45),sev=c(50,50)) # First bisection step
  question.1 <- add.choice.prob(coefficients,question.1)
  
  if(runif(1)<=question.1$choice.prob[1]) { 
    
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(60,50),mod=c(85,45),sev=c(50,50)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- lowerRatioConstraint(3,1,2,4)
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,1,2,2),upperRatioConstraint(3,1,2,4))
    }
    
  } else { 
    
    question.2 <- data.frame(q.nr=rep(3,2),alt=c("A","B"),PFS=c(80,50),mod=c(85,45),sev=c(50,50)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- mergeConstraints(lowerRatioConstraint(3,1,2,1+1/3),upperRatioConstraint(3,1,2,2))
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,1,2,1),upperRatioConstraint(3,1,2,1+1/3))
    }
    
  }
  
  constr
  
}

# Choice-based matching to trade-off mod against sev
cbm.mod.sev <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,2),alt=c("A","B"),PFS=c(70,70),mod=c(65,85),sev=c(80,20)) # First bisection step
  question.1 <- add.choice.prob(coefficients,question.1)
  
  if(runif(1)<=question.1$choice.prob[1]) { 
    
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(70,70),mod=c(75,85),sev=c(80,20)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) {
      constr <- lowerRatioConstraint(3,2,3,4)
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,2,3,2),upperRatioConstraint(3,2,3,4))
    }
    
  } else { 
    
    question.2 <- data.frame(q.nr=rep(3,2),alt=c("A","B"),PFS=c(70,70),mod=c(55,65),sev=c(80,20)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- mergeConstraints(lowerRatioConstraint(3,2,3,1+1/3),upperRatioConstraint(3,2,3,2))
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,2,3,1),upperRatioConstraint(3,2,3,1+1/3))
    }
    
  }
  
  constr
  
}

# Choice-based matching to trade-off mod against PFS
cbm.mod.PFS <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,2),alt=c("A","B"),PFS=c(50,90),mod=c(65,85),sev=c(50,50)) # First bisection step
  question.1 <- add.choice.prob(coefficients,question.1)
  
  if(runif(1)<=question.1$choice.prob[1]) { 
    
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(50,90),mod=c(75,85),sev=c(50,50)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) {
      constr <- lowerRatioConstraint(3,2,1,4)
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,2,1,2),upperRatioConstraint(3,2,1,4))
    }
    
  } else { 
    
    question.2 <- data.frame(q.nr=rep(3,2),alt=c("A","B"),PFS=c(50,90),mod=c(55,85),sev=c(50,50)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- mergeConstraints(lowerRatioConstraint(3,2,1,1+1/3),upperRatioConstraint(3,2,1,2))
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,2,1,1),upperRatioConstraint(3,2,1,1+1/3))
    }
    
  }
  
  constr
  
}

# Choice-based matching to trade-off sev against PFS
cbm.sev.PFS <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,2),alt=c("A","B"),PFS=c(50,90),mod=c(65,65),sev=c(50,80)) # First bisection step
  question.1 <- add.choice.prob(coefficients,question.1)
  
  if(runif(1)<=question.1$choice.prob[1]) { 
    
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(50,90),mod=c(65,65),sev=c(65,80)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) {
      constr <- lowerRatioConstraint(3,3,1,4)
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,3,1,2),upperRatioConstraint(3,3,1,4))
    }
    
  } else { 
    
    question.2 <- data.frame(q.nr=rep(3,2),alt=c("A","B"),PFS=c(50,90),mod=c(65,65),sev=c(35,80)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- mergeConstraints(lowerRatioConstraint(3,3,1,1+1/3),upperRatioConstraint(3,3,1,2))
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,3,1,1),upperRatioConstraint(3,3,1,1+1/3))
    }
    
  }
  
  constr
  
}

# Choice-based matching to trade-off sev against mod
cbm.sev.mod <- function(coefficients) {
  
  question.1 <- data.frame(q.nr=rep(1,2),alt=c("A","B"),PFS=c(70,70),mod=c(85,45),sev=c(50,80)) # First bisection step
  question.1 <- add.choice.prob(coefficients,question.1)
  
  if(runif(1)<=question.1$choice.prob[1]) { 
    
    question.2 <- data.frame(q.nr=rep(2,2),alt=c("A","B"),PFS=c(70,70),mod=c(45,85),sev=c(65,80)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) {
      constr <- lowerRatioConstraint(3,3,2,4)
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,3,2,2),upperRatioConstraint(3,3,2,4))
    }
    
  } else { 
    
    question.2 <- data.frame(q.nr=rep(3,2),alt=c("A","B"),PFS=c(70,70),mod=c(85,45),sev=c(35,80)) # Second bisection step
    question.2 <- add.choice.prob(coefficients,question.2)
    
    if(runif(1)<=question.2$choice.prob[1]) { 
      constr <- mergeConstraints(lowerRatioConstraint(3,3,2,1+1/3),upperRatioConstraint(3,3,2,2))
    } else {
      constr <- mergeConstraints(lowerRatioConstraint(3,3,2,1),upperRatioConstraint(3,3,2,1+1/3))
    }
    
  }
  
  constr
  
}

qen.MNL.weights.patient <- function(coefficients) {
  
  ord.swing <- ordinal.swing.MNL(coefficients)
  constr <- mergeConstraints(simplexConstraints(3),do.call(paste0("cbm.",ord.swing$first,".",ord.swing$second),list(coefficients=coefficients)))
  constr <- mergeConstraints(constr,do.call(paste0("cbm.",ord.swing$second,".",ord.swing$third),list(coefficients=coefficients)))
  
  # Generate representative weights by applying HAR sampling (slack-LP with random slack)
  weights <- colMeans(hitandrun(constr,n.samples=1,x0.randomize=T))
  names(weights) <- c("PFS","mod","sev")
  
  weights
  
}

qen.MNL.weights.survey <- function(n.subjects,coefficients) {
  
  weights <- c()
  for (i in 1:n.subjects) {
    weights <- rbind(weights,qen.MNL.weights.patient(coefficients))
  }
  
  weights
  
}


################### Example #########################################

#set.seed(1911)
#rum.fullsample <- simulate.dce(n.questions=16, n.respondents=560)
#coefficients <- rum.fullsample$mnl$coefficients
#survey.results <- qen.MNL.weights.survey(200,coefficients) 

####################################################################


#### Choice-based matching with Dirichlet population model ####

# Simulates a choice-based matching procedure
# Assumption: additive value function with linear partial value functions
#
#' @param w respondent's true weight vector
#' @param n.steps number of bisection steps in the choice-based matching procedure (per pairwise comparison)
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
