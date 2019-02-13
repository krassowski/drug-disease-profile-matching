from matplotlib import pyplot as plt


def y_line(p, score_mean, color='#bf4045', width=1.5, label='', n=1):
    x, y = p.get_lines()[n].get_data()
    y_by_dist_from_mean = {}
    for xi, yi in zip(x, y):
        dist = abs(xi - score_mean)
        y_by_dist_from_mean[dist] = yi
    closest = min(y_by_dist_from_mean.keys())

    plt.vlines(score_mean, 0, y_by_dist_from_mean[closest], linewidth=width, colors=color, label=label)


