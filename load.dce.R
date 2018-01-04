library(mlogit)
library(plyr)

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
