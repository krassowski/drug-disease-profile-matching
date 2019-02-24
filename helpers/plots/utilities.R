trans <- function(values, reference, margin=0.5) {
    mi = min(reference)
    (
        (values - mi)
        /
        (max(reference) - mi)
    ) + margin
}

trans_seq <- function(reference, side, zoom, margin=0.5, start_adjustment=0, end_shift=0) {
    start = round(max(reference) + end_shift, floor(log10(zoom)))
    range = max(reference) - min(reference)
    seq(
        side * trans(start, reference, margin),
        to=margin * side * (1 - start_adjustment),
        by=1/(range * zoom) * side * -1
    )
}

labels_seq <- function(reference, side, zoom, end_shift=0) {
    level = floor(log10(zoom))
    end = round(max(reference) + end_shift, level)
    start = round(min(reference), level)
    if(side == 1) {
        e = end
        end = start
        start = e
    }
    s = seq(
        end,
        to=start,
        by=1/(zoom) * side
    )
    if (side == 1) {
        s = rev(s)
    }
    s
}

subplot <- function(letter) {
    annotate('text', x=-Inf, y=Inf, label=letter, hjust=-1, vjust=1)
}
