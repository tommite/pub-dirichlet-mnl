library(dplyr)
library(smaa)
library(ggplot2)
library(reshape2)
library(MCMCprecision)
library(hitandrun)

source('load.dce.R')

## Rotate points 90% clockwise
rotate <- function(x) {
  c(x[2],-1*x[1])
}

## Visualize density in fancyplot style

visualize.density <- function(data,params) {
  
  mean <- params$alpha/params$sum
  
  # Transform weights to n-1 dimensions
  eq.constr <- list(constr = t(rep(1, 3)), dir = "=", rhs = 1)
  basis <- solution.basis(eq.constr)
  data.transformed <- c()
  for (i in 1:dim(data)[1]) {
    data.transformed <- rbind(data.transformed,rotate(as.vector(t(basis$basis) %*% (as.numeric(data[i,]) - basis$translate))))
  }
  colnames(data.transformed) <- c("x1","x2")
  data.transformed <- data.frame(data.transformed)
  
  # Obtain the coordinates of the simplex's vertices
  extreme.1 <- rotate(as.vector(t(basis$basis) %*% (c(1,0,0) - basis$translate)))
  extreme.2 <- rotate(as.vector(t(basis$basis) %*% (c(0,1,0) - basis$translate)))
  extreme.3 <- rotate(as.vector(t(basis$basis) %*% (c(0,0,1) - basis$translate)))
  
  p <- ggplot(data.transformed,aes(x=x1,y=x2))
  
  p <- p + scale_x_continuous(limits=c(-0.75,0.75),name=NULL,breaks=NULL)
  p <- p + scale_y_continuous(limits=c(-0.6,0.9),name=NULL,breaks=NULL)
  
  # Add coloured polygons
  polygon.1 <- data.frame(x=c(extreme.1[1],(extreme.1[1]+extreme.2[1])/2,0),y=c(extreme.1[2],(extreme.1[2]+extreme.2[2])/2,0))
  polygon.2 <- data.frame(x=c(extreme.1[1],(extreme.1[1]+extreme.3[1])/2,0),y=c(extreme.1[2],(extreme.1[2]+extreme.3[2])/2,0))
  polygon.3 <- data.frame(x=c(extreme.3[1],(extreme.1[1]+extreme.3[1])/2,0),y=c(extreme.3[2],(extreme.1[2]+extreme.3[2])/2,0))
  polygon.4 <- data.frame(x=c(extreme.3[1],(extreme.2[1]+extreme.3[1])/2,0),y=c(extreme.3[2],(extreme.2[2]+extreme.3[2])/2,0))
  polygon.5 <- data.frame(x=c(extreme.2[1],(extreme.2[1]+extreme.3[1])/2,0),y=c(extreme.2[2],(extreme.2[2]+extreme.3[2])/2,0))
  polygon.6 <- data.frame(x=c(extreme.2[1],(extreme.1[1]+extreme.2[1])/2,0),y=c(extreme.2[2],(extreme.1[2]+extreme.2[2])/2,0))
  
  data.polygon <- rbind(polygon.1,polygon.2,polygon.3,polygon.4,polygon.5,polygon.6)
  data.polygon$ranking <- c(rep("PFS>mod>sev (n=215)",3),rep("PFS>sev>mod (n=197)",3),rep("sev>PFS>mod (n=109)",3),rep("sev>mod>PFS (n=9)",3),rep("mod>sev>PFS (n=11)",3),rep("mod>PFS>sev (n=19)",3))
  
  p <- p + geom_polygon(aes(x=x,y=y,group=ranking,fill=ranking),data=data.polygon,alpha=0.8)
  p <- p + theme(panel.background=element_rect(fill="white"),legend.position="bottom",legend.title=element_blank(),legend.direction="horizontal",
                 legend.margin=margin(t=-0.75, unit='cm'),aspect.ratio=1,legend.text=element_text(size=7), plot.margin=unit(c(-0.5,0,0.5,0),unit="cm"))
  
  # Add coordinate lines
  ticks <- seq(0,1,0.1)
  
  data.grid.1 <- data.frame(y=rep(extreme.2[2],length(ticks)) + ticks*(extreme.1[2] - extreme.2[2]),yend=rep(extreme.2[2],length(ticks)) + ticks*(extreme.1[2] - extreme.2[2]),
                            x=extreme.2[1]*(1-ticks),xend=(extreme.3[1]*(1-ticks)))
  p <- p + geom_segment(data=data.grid.1,mapping=aes(x=x,xend=xend,y=y,yend=yend),colour="white")
  p <- p + geom_text(data=data.frame(data.grid.1,label=as.character(ticks)),mapping=aes(x=xend,y=yend,label=label),nudge_x=0.02,nudge_y=0.02,hjust="left",size=2)
  
  data.grid.2 <- data.frame(y=rep(extreme.1[2],length(ticks)) + ticks*(extreme.2[2] - extreme.1[2]),yend=rep(extreme.2[2],length(ticks)),
                            x=rep(extreme.1[1],length(ticks)) + ticks*(extreme.2[1] - extreme.1[1]),xend=rep(extreme.3[1],length(ticks)) + ticks*(extreme.2[1] - extreme.3[1]))
  p <- p + geom_segment(data=data.grid.2,mapping=aes(x=x,xend=xend,y=y,yend=yend),colour="white")
  p <- p + geom_text(data=data.frame(data.grid.2,label=as.character(ticks)),mapping=aes(x=x,y=y,label=label),nudge_x=-0.03,nudge_y=-0.02,hjust="right",angle=300,size=2)
  
  data.grid.3 <- data.frame(y=rep(extreme.1[2],length(ticks)) + ticks*(extreme.3[2] - extreme.1[2]),yend=rep(extreme.3[2],length(ticks)),
                            x=rep(extreme.1[1],length(ticks)) + ticks*(extreme.3[1] - extreme.1[1]),xend=rep(extreme.2[1],length(ticks)) + ticks*(extreme.3[1] - extreme.2[1]))
  p <- p + geom_segment(data=data.grid.3,mapping=aes(x=x,xend=xend,y=y,yend=yend),colour="white")
  p <- p + geom_text(data=data.frame(data.grid.3,label=as.character(ticks)),mapping=aes(x=xend,y=yend,label=label),nudge_x=0.02,nudge_y=-0.02,hjust="right",angle=60,size=2)
  
  # Add axis labels
  p <- p +geom_text(data=NULL,x=(extreme.1[1]+extreme.3[1])/2 + 0.1,y=(extreme.1[2]+extreme.2[2])/2 + 0.075,label="1-year progression-free survival (50% -> 90%)",angle=300,size=3)
  p <- p +geom_text(data=NULL,x=(extreme.1[1]+extreme.2[1])/2 - 0.1,y=(extreme.1[2]+extreme.2[2])/2 + 0.075,label="Moderate chronic toxicity (85% -> 45%)",angle=60,size=3)
  p <- p +geom_text(data=NULL,x=0,y=extreme.2[2] - 0.12,label="Severe toxicity (80% -> 20%)",size=3)
  
  # Add simplex boundary
  p <- p + geom_segment(data=NULL,x=extreme.1[1],xend=extreme.2[1],y=extreme.1[2],yend=extreme.2[2],colour="black")
  p <- p + geom_segment(data=NULL,x=extreme.1[1],xend=extreme.3[1],y=extreme.1[2],yend=extreme.3[2],colour="black")
  p <- p + geom_segment(data=NULL,x=extreme.2[1],xend=extreme.3[1],y=extreme.2[2],yend=extreme.3[2],colour="black")
  
  # Plot individual samples
  p <- p + geom_point(size=1,alpha=0.8)
  
  # Plot mean of the fitted Dirichlet distribution
  mean.transformed <- rotate(as.vector(t(basis$basis) %*% (as.numeric(mean) - basis$translate)))
  p <- p + geom_point(x=mean.transformed[1],y=mean.transformed[2],size=2,colour="red",shape=23,fill="red")
  
  p
  
}

## Generate random sample from the Dirichlet distribution fitted to our dataset
params <- dirichlet.mle(df.w[,c("pfs","mod","sev")])
samples <- rdirichlet(500,params$alpha)
colnames(samples) <- c("pfs","mod","sev")

mean <- params$alpha/params$sum
names(mean) <- c("pfs","mod","sev")

## Create plot 
p <- visualize.density(samples,params)
ggsave("samples_500.pdf",plot=p,width=14.8,height=14.8,units="cm")

