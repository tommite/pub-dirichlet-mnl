library(mlogit)
library(plyr)
library(colorspace)
library(support.CEs)

## Dirichlet to Douwe's data

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

attribute.names <- list('PFS'=sort(unique(df$level.PFS)),
                        'mod'=sort(unique(df$level.mod)),
                        'sev'=sort(unique(df$level.sev)))

design <- Lma.design(attribute.names=attribute.names, nalternatives=2, nblocks=1)
## filter out dominated alternative questions
ok.qs <- laply(rownames(design$alternatives$alt.1), function(q.idx) {
    q1 <- as.numeric(as.matrix(design$alternatives$alt.1[q.idx,c('PFS', 'mod', 'sev')]))
    q2 <- as.numeric(as.matrix(design$alternatives$alt.2[q.idx,c('PFS', 'mod', 'sev')]))
    !((q1[1] >= q2[1] && q1[2] <= q2[2] && q1[3] <= q2[3]) ||
      (q1[1] <= q2[1] && q1[2] >= q2[2] && q1[3] >= q2[3]))
})

design.nondom <- list(alt.1=design$alternatives$alt.1[ok.qs,],
                      alt.2=design$alternatives$alt.2[ok.qs,])

##
##df <- subset(df, url %in% df.w.highpfs$url)

## Fit MNL
mdata <- mlogit.data(df, choice='selected.by.subject',
                     ch.id='idx',
                     id.var='url',
                     shape='long',
                     alt.var='alt')
res <- mlogit(selected.by.subject ~ 0 + level.PFS + level.mod + level.sev,
#              rpar=c(level.PFS='n', level.mod='n', level.sev='n'),
              data=mdata)

## Normalize weights to scale
scales <- c(diff(range(df$level.PFS)), diff(range(df$level.mod)), diff(range(df$level.sev)))
norm.to.scale <- abs(as.matrix(res$coefficients) * scales)
norm.weights <- norm.to.scale / sum(norm.to.scale)

print(summary(res))
print(norm.weights)

## L^ma Design
design <- Lma.design(attribute.names=list(
                         'pfs'=c('50', '60', '70', '80', '90'),
                         'mod'=c('45', '55', '65', '75', '85'),
                         'sev'=c('20', '35', '50', '65', '80')),
                     continuous.attributes=c('pfs', 'mod', 'sev'),
                     nalternatives=2,
                     nblocks=1,
                     seed=1911)

