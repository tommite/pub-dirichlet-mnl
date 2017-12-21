library(mlogit)

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
