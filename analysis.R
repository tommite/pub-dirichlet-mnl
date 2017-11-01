library(mlogit)
library(plyr)

data.file <- 'data/data_DCE.csv'

df.raw <- read.csv2(data.file, header=TRUE, sep=',', stringsAsFactors=FALSE)
df <- df.raw
df[,'idx'] <- paste0(df$url, '.', df$question.no)

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
alt.vars <- c('a', 'b')
df[,'alt'] <- rep(alt.vars, times=nrow(df) / length(alt.vars))

## Fit MNL
mdata <- mlogit.data(df, choice='selected.by.subject',
                     ch.id='idx',
                     id.var='url',
                     shape='long',
                     alt.var='alt')
res <- mlogit(selected.by.subject ~ 0 + level.PFS + level.mod + level.sev, data=mdata)

## Normalize weights to scale
scales <- c(diff(range(df$level.PFS)), diff(range(df$level.mod)), diff(range(df$level.sev)))
norm.to.scale <- abs(as.matrix(res$coefficients) * scales)
norm.weights <- norm.to.scale / sum(norm.to.scale)

summary(res)
print(norm.weights)
