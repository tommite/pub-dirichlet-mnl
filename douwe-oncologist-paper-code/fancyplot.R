library(ggplot2)
library(hitandrun)

df.w <- read.csv2(data.file.w, header=TRUE, sep=',', stringsAsFactors=FALSE)

crit.names <- c('pfs', 'mod', 'sev')

weights <- aaply(df.w, 1, function(x) {
    laply(paste0('d', 1:3), function(c) {
        sum(as.numeric(x[grep(c, names(df.w))]))
    })
}, .expand=FALSE)

colnames(weights) <- crit.names

data.PFS <- weights

# Function for rotating points 90% clockwise
rotate <- function(x) {
  c(x[2],-1*x[1])
}

# Transform weights to n-1 dimensions
data <- data.PFS[,c("pfs","mod","sev")]
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

# Add additional line segments (inner triangles)
#p <- p + geom_segment(data=NULL,x=(extreme.1[1]+extreme.2[1])/2,xend=(extreme.1[1]+extreme.3[1])/2,y=(extreme.1[2]+extreme.2[2])/2,yend=(extreme.1[2]+extreme.2[2])/2,linetype="solid",size=1,colour="darkgray")
#p <- p + geom_segment(data=NULL,x=(extreme.1[1]+extreme.2[1])/2,xend=(extreme.2[1]+extreme.3[1])/2,y=(extreme.1[2]+extreme.2[2])/2,yend=extreme.2[2],linetype="solid",size=1,colour="darkgray")
#p <- p + geom_segment(data=NULL,x=(extreme.1[1]+extreme.3[1])/2,xend=(extreme.2[1]+extreme.3[1])/2,y=(extreme.1[2]+extreme.3[2])/2,yend=extreme.3[2],linetype="solid",size=1,colour="darkgray")

#p <- p + annotate("text",x=extreme.1[1],y=extreme.1[2] - 0.15,label="PFS",colour="darkgray",size=5,fontface="bold")
#p <- p + annotate("text",x=(extreme.2[1]+0.15),y=(extreme.2[2]+0.075),label="mod",colour="darkgray",size=5,fontface="bold",angle=300)
#p <- p + annotate("text",x=(extreme.3[1]-0.15),y=(extreme.3[2]+0.075),label="sev",colour="darkgray",size=5,fontface="bold",angle=60)

# Add simplex boundary
p <- p + geom_segment(data=NULL,x=extreme.1[1],xend=extreme.2[1],y=extreme.1[2],yend=extreme.2[2],colour="black")
p <- p + geom_segment(data=NULL,x=extreme.1[1],xend=extreme.3[1],y=extreme.1[2],yend=extreme.3[2],colour="black")
p <- p + geom_segment(data=NULL,x=extreme.2[1],xend=extreme.3[1],y=extreme.2[2],yend=extreme.3[2],colour="black")

# Add individual weights
p <- p + geom_point(size=1,alpha=0.8)
p <- p + geom_point(x=0.12801,y=0.24789,size=2,colour="red",shape=23,fill="red")

#ggsave("simplex_revised4.pdf",plot=p,width=14.8,height=14.8,units="cm")
ggsave("Figure 3.tiff",plot=p,dpi=300,width=11.4,height=11.4,units="cm")
