## vary number of questions
res.vary.q <- llply(seq(from=3, to=nrow(design.nondom)/2, by=1), error.catch.simulate.dce,
                    n.respondents=200, n.simul=20)

## Plot test stats vary n questions
df.molten.q <- melt(as.data.frame(test.stats.vary.q),
                  measure.vars=c(names(res.vary.n[[1]]$results[[1]]$coefficients)))

plots.q <- dlply(df.molten.q, 'variable', function(df.plot) {
    df.plot$n.quest <- factor(df.plot$n.quest,
                              labels=unique(df.plot$n.quest))
    ggplot(df.plot, aes(x=n.quest, y=value)) +
        geom_boxplot(outlier.colour='red', outlier.shape=20) +
        ylab('p-value') + theme_economist() + scale_colour_economist() +
        ggtitle(unique(df.plot$variable))
})
dev.new(width=8, height=6)
grid.arrange(plots.q[[1]], plots.q[[2]], plots.q[[3]], ncol=1)


test.stats.vary.q <- test.stats.p(res.vary.q)
