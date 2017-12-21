### CLuster analysis of the survey data ###

library("cluster")
library("colorspace")
library("ggplot2")
library("hitandrun")
library("car")
library("compositions")
library("grid")
library("smaa")

# Compare distribution of baseline characteristics across subgroups
compare.baseline <- function(data,characteristic,group) {
  cur.table <- table(data[,characteristic],data[,group])
  print(100*round(prop.table(cur.table,margin=2),digits=2))
  print(chisq.test(cur.table,simulate.p.value=TRUE))
}

compare.baseline.cont <- function(data,characteristic,outcome) {
  for (j in levels(data[,characteristic])) {
    print(j)
    cur.data <- data[,outcome][data[,characteristic]==j]
    print(round(quantile(cur.data,probs=c(0.25,0.50,0.75),na.rm=T),digits=2))
  }
  test <- kruskal.test(data[,outcome],data[,characteristic])
  print(test)
}

covariates <- c("sex","age","caucasian","working","alone","dependent_family",
                "diagnosis","current_status","treatment_hist","side_effects",
                "patient_organisation","chronic_health")

data.PFS <- read.csv("data_PFS_har.csv")
data.PFS$diagnosis <- ordered(recode(data.PFS$diagnosis,"c('Between 4 and 6 years ago','More than 6 years ago')='More than 4 years ago'"),levels=c("In the last 2 years","Between 2 and 4 years ago","More than 4 years ago"))
data.PFS$treatment_hist <- recode(data.PFS$treatment_hist,"c('None','1 line of treatment')='Less than 2 lines of treatment'")

feedback.coded <- read.xlsx("feedback_coded_final.xlsx",1)
feedback.coded$Code <- ifelse(is.na(feedback.coded$Code),"U",feedback.coded$Code)
feedback.coded$Code <- factor(feedback.coded$Code,labels=c("C","D","U"))
summary(feedback.coded$Code)

data.PFS$feedback.code <- feedback.coded$Code

### Assessment of the partial value functions ###

tiff("Figure 2.tiff",width=14.8,height=4,units="cm",res=300)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,3)))

data.boxplot.PFS <- data.frame(y=c(rep(0,560),data.PFS$X1,data.PFS$X1+data.PFS$X2,data.PFS$X1+data.PFS$X2+data.PFS$X3,data.PFS$pfs),group=c(rep(50,560),rep(60,560),rep(70,560),rep(80,560),rep(90,560)))
ggp <- ggplot(data.boxplot.PFS,aes(group,y))
ggp <- ggp + geom_boxplot(aes(group=group),size=0.25,outlier.size=0.25)
ggp <- ggp + labs(x="1-year progression-free survival (%)",y="Part-worth")
ggp <- ggp + ylim(0,1)

means.ext <- c(0,mean(data.PFS$X1),mean(data.PFS$X1+data.PFS$X2),mean(data.PFS$X1+data.PFS$X2+data.PFS$X3),mean(data.PFS$pfs))
ggp <- ggp + geom_line(aes(x=x,y=y),data=data.frame(x=(seq(50,90,10)),y=means.ext),size=0.25,col="red")
ggp <- ggp + geom_point(aes(x=x,y=y),data=data.frame(x=(seq(50,90,10)),y=means.ext),size=0.5,col="red")

ggp <- ggp + theme(title=element_text(size=5),axis.text=element_text(size=3))

print(ggp,vp=viewport(layout.pos.row = 1,layout.pos.col = 1))

#tiff("pvf_pfs.tiff",res=280,width=8,height=8,units="cm")
#print(ggp)
#dev.off()

data.boxplot.mod <- data.frame(y=c(data.PFS$mod,data.PFS$X5+data.PFS$X6+data.PFS$X7,data.PFS$X5+data.PFS$X6,data.PFS$X5,rep(0,560)),group=c(rep(45,560),rep(55,560),rep(65,560),rep(75,560),rep(85,560)))
ggp <- ggplot(data.boxplot.mod,aes(group,y))
ggp <- ggp + geom_boxplot(aes(group=group),size=0.25,outlier.size=0.25)
ggp <- ggp + labs(x="Moderate chronic toxicity (%)",y="Part-worth")
ggp <- ggp + ylim(0,1)
ggp <- ggp + scale_x_continuous(breaks=seq(45,85,10))

means.ext <- sort(c(0,mean(data.PFS$X5),mean(data.PFS$X5+data.PFS$X6),mean(data.PFS$X5+data.PFS$X6+data.PFS$X7),mean(data.PFS$mod)),decreasing = T)
ggp <- ggp + geom_line(aes(x=x,y=y),data=data.frame(x=(seq(45,85,10)),y=means.ext),size=0.25,col="red")
ggp <- ggp + geom_point(aes(x=x,y=y),data=data.frame(x=(seq(45,85,10)),y=means.ext),size=0.5,col="red")

ggp <- ggp + theme(title=element_text(size=5),axis.text=element_text(size=3))

print(ggp,vp=viewport(layout.pos.row = 1,layout.pos.col = 2))

#tiff("pvf_mod.tiff",res=280,width=8,height=8,units="cm")
#print(ggp)
#dev.off()

data.boxplot.sev <- data.frame(y=c(data.PFS$sev,data.PFS$X9+data.PFS$X10+data.PFS$X11,data.PFS$X9+data.PFS$X10,data.PFS$X9,rep(0,560)),group=c(rep(20,560),rep(35,560),rep(50,560),rep(65,560),rep(80,560)))
ggp <- ggplot(data.boxplot.sev,aes(group,y))
ggp <- ggp + geom_boxplot(aes(group=group),size=0.25,outlier.size=0.25)
ggp <- ggp + labs(x="Severe toxicity (%)",y="Part-worth")
ggp <- ggp + ylim(0,1)
ggp <- ggp + scale_x_continuous(breaks=seq(20,80,15))

means.ext <- sort(c(0,mean(data.PFS$X9),mean(data.PFS$X9+data.PFS$X10),mean(data.PFS$X9+data.PFS$X10+data.PFS$X11),mean(data.PFS$sev)),decreasing = T)
ggp <- ggp + geom_line(aes(x=x,y=y),data=data.frame(x=(seq(20,80,15)),y=means.ext),size=0.25,col="red")
ggp <- ggp + geom_point(aes(x=x,y=y),data=data.frame(x=(seq(20,80,15)),y=means.ext),size=0.5,col="red")

ggp <- ggp + theme(title=element_text(size=5),axis.text=element_text(size=3))

print(ggp,vp=viewport(layout.pos.row = 1,layout.pos.col = 3))

#tiff("pvf_sev.tiff",res=280,width=8,height=8,units="cm")
#print(ggp)
#dev.off()

dev.off()

### Assessment of the attribute weights ###

data.PFS$ord_pfs <- ifelse(data.PFS$pfs>data.PFS$sev & data.PFS$pfs>data.PFS$mod,1,0)
data.PFS$ord_sev <- ifelse(data.PFS$sev>data.PFS$pfs & data.PFS$sev>data.PFS$mod,1,0)
data.PFS$ord_mod <- ifelse(data.PFS$mod>data.PFS$pfs & data.PFS$mod>data.PFS$sev,1,0)
data.PFS$ordinal <- as.factor(ifelse(data.PFS$ord_pfs==1,"pfs",ifelse(data.PFS$sev>data.PFS$mod,"sev","mod")))

#data.PFS$group <- ifelse(data.PFS$ordinal=="pfs",ifelse(data.PFS$sev>data.PFS$mod,"pfs_sev","pfs_mod"),data.PFS$ordinal)
data.PFS$group <- ifelse(data.PFS$sev>data.PFS$mod,"sev","mod")

#data.PFS$group <- ifelse(data.PFS$sev_os>data.PFS$mod_os,"sev","mod")

for (i in covariates) {
  print(i)
  compare.baseline(data.PFS[data.PFS$consistent==1,],i,"group")
}

### Compositional data analysis ###

weights <- data.PFS[,c("pfs","sev","mod")]
comps <- acomp(weights)

mean(comps)
var <- variation(comps)
total.var <- sum(var)/6

transform.1 <- alr(weights)

summary(manova(transform.1~sex,data=data.PFS))
summary(manova(transform.1~age,data=data.PFS))
summary(manova(transform.1~working,data=data.PFS))
summary(manova(transform.1~alone,data=data.PFS))
summary(manova(transform.1~dependent_family,data=data.PFS))
summary(manova(transform.1~diagnosis,data=data.PFS))
summary(manova(transform.1~current_status,data=data.PFS))
summary(manova(transform.1~treatment_hist,data=data.PFS))
summary(manova(transform.1~side_effects,data=data.PFS))
summary(manova(transform.1~patient_organisation,data=data.PFS))
summary(manova(transform.1~chronic_health,data=data.PFS))
summary(manova(transform.1~caucasian,data=data.PFS))

### Visualization ###

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

### Consistency checkning and sensitivity analysis ###

# Consistency checking
data.OS <- read.csv("data_OS_har.csv")
data.PFS$ord_os <- ifelse(data.OS$os>data.OS$sev & data.OS$os>data.OS$mod,1,0)
data.PFS$os <- data.OS$os
data.PFS$sev_os <- data.OS$sev
data.PFS$mod_os <- data.OS$mod

data.PFS$sev.mod_os <- log(data.PFS$sev_os/data.PFS$mod_os)
data.PFS$sev.mod_pfs <- log(data.PFS$sev/data.PFS$mod)

data.PFS$sev.gr.mod_pfs <- ifelse(data.PFS$sev>data.PFS$mod,1,0)
data.PFS$sev.gr.mod_os <- ifelse(data.PFS$sev_os>data.PFS$mod_os,1,0)
table(data.PFS$sev.gr.mod_pfs,data.PFS$sev.gr.mod_os)

data.PFS$consistent <- as.numeric(data.PFS$sev.gr.mod_pfs==data.PFS$sev.gr.mod_os)
prop.table(table(data.PFS$consistent,data.PFS$feedback.code),margin=2)

#pdf("consistency.pdf")
plot(log(data.PFS$sev.mod_pfs),log(data.PFS$sev.mod_os),xlim=c(-3,5),ylim=c(-3,5),xlab="model for progression-free survival",ylab="model for overall survival")
abline(0,1)
abline(v=0)
abline(h=0)
#dev.off()

a <- ggplot(data=data.PFS,mapping=aes(x=sev.mod_pfs,y=sev.mod_os))
a <- a + geom_point()
a <- a + xlim(-3,5)
a <- a + ylim(-3,5)
a <- a + geom_hline(yintercept=0)
a <- a + geom_vline(xintercept=0)
a <- a + xlab("model for progression-free survival")
a <- a + ylab("model for overall survival")

pdf("consistency.pdf")
print(a)
dev.off()

cor(data.PFS$sev.mod_pfs,data.PFS$sev.mod_os,method="s")
cor(data.PFS$sev.mod_pfs[data.PFS$consistent==1],data.PFS$sev.mod_os[data.PFS$consistent==1],method="s")

cor(data.PFS$sev.mod_pfs[data.PFS$sev.gr.mod_pfs==0],data.PFS$sev.mod_os[data.PFS$sev.gr.mod_pfs==0],method="s")

for (i in covariates) {
  print(i)
  compare.baseline(data.PFS,i,"consistent")
}

# Sensitivity analysis
data.PFS$pfs.tox <- log(data.PFS$pfs/(1-data.PFS$pfs))
data.PFS$os.tox <- log(data.PFS$os/(1-data.PFS$os))

plot(data.PFS$pfs.tox,data.PFS$os.tox)
abline(0,1)
abline(v=0)
abline(h=0)

cor(data.PFS$pfs.tox,data.PFS$os.tox,method="s")

plot(data.PFS$sev+data.PFS$mod,data.PFS$sev_os+data.PFS$mod_os)
abline(0,1)
abline(v=0.5)
abline(h=0.5)

cor(data.PFS$sev+data.PFS$mod,data.PFS$sev_os+data.PFS$mod_os,method="s")

summary(1-data.PFS$pfs)
summary(1-data.PFS$os)

boxplot(1-data.PFS$pfs,1-data.PFS$os,ylab="Cumulative weight of the toxicities",names=c("PFS","OS"))
wilcox.test(1-data.PFS$pfs,1-data.PFS$os,paired=T)
