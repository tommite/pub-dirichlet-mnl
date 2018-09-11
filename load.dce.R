library(mlogit)
library(support.CEs)
library(plyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(ggthemes)
library(MCMCprecision)
library(devEMF)
library(evd)
source('dirichlet-cvm.R')
source('simulate-cbm.R')
source('pinkfloyd-plot.R')

set.seed(1911)
memory.limit(size=16000)
plot.theme <- theme_economist(base_size=20)
n.simul <- 100
n.dir.samples <- 1E3

data.file <- 'data/data_DCE.csv'
data.file.w <- 'data/data_PFS_har.csv'

df.raw <- read.csv(data.file, header=TRUE, sep=',', stringsAsFactors=FALSE)
df <- df.raw
df[,'idx'] <- paste0(df$url, '.', df$question.no)

## Load weight data
df.w <- read.csv(data.file.w, header=TRUE, sep=',', stringsAsFactors=FALSE)
df.w.highpfs <- subset(df.w, pfs > sev)
df.w.highpfs <- subset(df.w.highpfs,  sev > mod)


## Transform first question to 2 preference statements
df.q1 <- subset(df, question.no == 1)

df.q1.2alt <- ldply(unique(df.q1$url), function(responder) {
    slice <- df.q1[df.q1$url == responder,]

    idx.selected <- which(slice$selected.by.subject == 1)
    idx.others <- setdiff(1:3, idx.selected)
    q1 <- slice[c(idx.selected, idx.others[1]),]
    q2 <- slice[c(idx.selected, idx.others[2]),]
    q1$idx <- paste0(q1$idx, '.', 1)
    q2$idx <- paste0(q2$idx, '.', 2)
    q1$question.no <- paste0(q1$question.no, '.', 1)
    q2$question.no <- paste0(q2$question.no, '.', 2)
    rbind(q1, q2)
})

## Combine
df <- rbind(df.q1.2alt, subset(df, question.no != 1))

## Add alt.var
alt.vars <- c('A', 'B')
df[,'alt'] <- rep(alt.vars, times=nrow(df) / length(alt.vars))

df <- df[,c(2:7, 22, 23)]

## Calculate population weight average
true.w <- colMeans(df.w[,c('pfs', 'mod', 'sev')])
names(true.w) <- c('PFS', 'mod', 'sev') # compatible names

## Generate L^ma non-dominated design ##
attribute.names <- list('PFS'=sort(unique(df$level.PFS)),
                        'mod'=sort(unique(df$level.mod)),
                        'sev'=sort(unique(df$level.sev)))
design <- Lma.design(attribute.names=attribute.names, nalternatives=2, nblocks=1)
ok.qs <- laply(rownames(design$alternatives$alt.1), function(q.idx) {
    q1 <- as.numeric(as.matrix(design$alternatives$alt.1[q.idx,c('PFS', 'mod', 'sev')]))
    q2 <- as.numeric(as.matrix(design$alternatives$alt.2[q.idx,c('PFS', 'mod', 'sev')]))
    !((q1[1] >= q2[1] && q1[2] <= q2[2] && q1[3] <= q2[3]) ||
      (q1[1] <= q2[1] && q1[2] >= q2[2] && q1[3] >= q2[3]))
})
a1.qs <- design$alternatives$alt.1[ok.qs,]
a2.qs <- design$alternatives$alt.2[ok.qs,]
a1.qs$q.nr <- 1:nrow(a1.qs)
a2.qs$q.nr <- 1:nrow(a2.qs)
a1.qs$alt <- 'A'
a2.qs$alt <- 'B'
design.nondom <- rbind(a1.qs, a2.qs)[,c('q.nr', 'alt', 'PFS', 'mod', 'sev')]
design.nondom <- design.nondom[order(design.nondom$q.nr),]
design.nondom$PFS <- as.numeric(as.vector(design.nondom$PFS))
design.nondom$mod <- as.numeric(as.vector(design.nondom$mod))
design.nondom$sev <- as.numeric(as.vector(design.nondom$sev))

ranges <- laply(attribute.names, range)
rownames(ranges) <- names(attribute.names)

##
#' Constructs a normalized weight vector from DCE coefficients
##
coeff.to.w <- function(b) {
    rng.sizes <- aaply(ranges, 1, diff)
    w <- b * rng.sizes
    abs(w / sum(abs(w)))
}

eucl.dist <- function(x, y) {
    stopifnot(length(x) == length(y))
    sqrt(sum((x-y)^2))
}

FINAL.SCALE <- 0.3

## Centralized error function - EV-1
dce.err.f <- function(scale=FINAL.SCALE) {rgumbel(1, loc=0, scale=scale)}

resp.w <- df.w[,c('url', 'pfs', 'mod', 'sev')]

make.design.matrix <- function(scale=FINAL.SCALE) {
    design.matrix <- adply(resp.w, 1, function(row) {
        rows <- design.nondom
        rows$url <- row$url
        rows$question.no <- paste0(rownames(row), '.', rows$q.nr)

        ldply(unique(rows$q.nr), function(q) {
            r <- subset(rows, q.nr %in% q)

            data <- as.matrix(r[,c('PFS', 'mod', 'sev')])
            ## convert to partial values
            pvs <- cbind(smaa.pvf(data[,'PFS'], cutoffs=ranges['PFS',], values=c(0,1)),
                         smaa.pvf(data[,'mod'], cutoffs=ranges['mod',], values=c(1,0)),
                         smaa.pvf(data[,'sev'], cutoffs=ranges['sev',], values=c(1,0)))
            colnames(pvs) <- colnames(data)
            vals <- as.matrix(row[-1]) %*% t(pvs)
            ## ADD error
            vals[1] <- vals[1] + dce.err.f(scale)
            vals[2] <- vals[2] + dce.err.f(scale)
            r$choice <- if(vals[1] > vals[2]) c(1, 0) else c(0, 1)
            r
        })
    }, .expand=FALSE)
    design.matrix$idx <- paste0(design.matrix$url, '.', design.matrix$question.no)
    design.matrix
}

