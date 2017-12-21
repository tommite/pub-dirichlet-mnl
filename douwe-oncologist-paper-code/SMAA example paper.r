### SMAA based on an input effects table with exact attribute measurements and the distribution of preferences in the study population 
### Author: Douwe Postmus (d.postmus@umcg.nl)

library("smaa")

### Additive value function 
# alt: the vector of attribute measurements associated with an alternative (PFS,G2,G34)
# differences: the vector of d_ij parameters used to parametrize the value function

value.function <- function(alt,differences) {
  
  partial.value.pfs <- smaa.pvf(alt[1],seq(50,90,10),c(0,differences[1],sum(differences[1:2]),sum(differences[1:3]),sum(differences[1:4])),outOfBounds="interpolate")
  partial.value.mod <- smaa.pvf(alt[2],seq(45,85,10),c(sum(differences[5:8]),sum(differences[5:7]),sum(differences[5:6]),differences[5],0),outOfBounds="interpolate")
  partial.value.sev <- smaa.pvf(alt[3],seq(20,80,15),c(sum(differences[9:12]),sum(differences[9:11]),sum(differences[9:10]),differences[9],0),outOfBounds="interpolate")
  
  overall.value <- partial.value.pfs + partial.value.mod + partial.value.sev
  
}

### Conduct SMAA ###

# Import preference data
data <- read.csv("preference_data.csv")

# Effect size estimates for the example included in the manuscript
treatment <- c(66,71,60) # Ixazomib
control <- c(59,69,53) # placebo

# Compute first rank acceptability index for Ixazomib
n <- dim(data)[1]
treatment.best <- rep(0,n)
for (i in 1:n) {
  if (value.function(treatment,data[i,])>value.function(control,data[i,])) {
    treatment.best[i] <- 1
  }
}
mean(treatment.best) # First rank acceptability index

