distribution_plot = function(
    scores, to_compare, ncol=2, n=5000, top_lim=1.5, scales='fixed', wrap='~ func', grid=F,
    p_annotations=NULL, color='group', mark=NULL, mark_shape=4, mark_name=NA
) {
    levels(scores$func) = sapply(levels(scores$func), latex2exp::TeX, 'expression')

    if (any(scores$group == 'unassigned') && n > 1) {
        # due to memory limitations, only n randomly selected unassigned profiles are plotted.
        unassigned_rows <- which(scores$group=='unassigned')
        unassigned <- scores[sample(unassigned_rows, n),]
        managable_subset <- rbind(scores[scores$group != 'unassigned',], unassigned)
    } else {
        managable_subset <- scores
    }

    if(grid) {
        grid = facet_grid
        grid_args = c(space='free')
    } else {
        grid = facet_wrap
        grid_args = c(ncol=ncol)
    }

    g <- (
        ggplot(scores, aes(x=group, y=score, color=group))
        + do.call(grid, c(wrap, labeller=label_parsed, scales=scales, grid_args))
    )

    if(!is.null(mark)) {
        g = (
            g
            + ggbeeswarm::geom_quasirandom(
                data=managable_subset,
                aes(
                    color=ifelse(get(mark), mark_name, group),
                    shape=ifelse(get(mark) == T, mark_shape, 19)
                )
            ) + scale_shape_identity()
        )
    } else {
        g = (
            g
            + ggbeeswarm::geom_quasirandom(data=managable_subset)
        )
    }

    g = (
        g
        + geom_boxplot(notch=T, alpha=0.5, outlier.shape = NA)
        + colors
        + theme(
            legend.position='bottom',
            legend.box='horizontal',
            text=element_text(size=15),
            #axis.text.x=element_text(angle=45, hjust=1)
            axis.ticks.x=element_blank(),
            axis.text.x=element_blank(),
            legend.margin=margin(t=-1, unit='cm')
        )
        + guides(color=guide_legend(title=''))
        + xlab('')
        + ylab('Normalized anticorrelation score')
    )

    if(!is.null(to_compare)) {
        g = (
            g
            + ggpubr::stat_compare_means(
                comparisons=to_compare,
                method="wilcox.test",
                method.args=list(alternative="greater")
            )
        )
    }

    if(!is.null(p_annotations)) {
        levels(p_annotations$func) = sapply(levels(p_annotations$func), latex2exp::TeX, 'expression')

        g = (
            g
            + ggsignif::geom_signif(
                data=p_annotations,
                aes(
                    xmin=start, xmax=end,
                    annotations=label,
                    y_position=y,
                    inherit.aes=F
                ),
                textsize=3.5, vjust=-0.2,
                manual=TRUE,
                color='black',
                size=0.25
            )
        )
    }
    if(is.numeric(top_lim)) {
        g = g + ylim(-1, top_lim)
    }
    g
}
