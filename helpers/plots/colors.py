from itertools import combinations
from math import sqrt

from IPython.core.display import display
from matplotlib.colors import to_rgba, to_hex
from plotnine import scales

from helpers.gui import NeatNamespace, HorizontalNamespace


def generate_colors(labels, scale=scales.scale_fill_hue()):
    return dict(zip(labels, map(Color.from_hex, scale.palette(len(labels)))))


def mix_colors(colors, to_mix, joined_label=lambda labels: ' and '.join(sorted(labels))):
    # make the mixed "clusters" have averaged colors
    for a, b in combinations(to_mix, 2):
        label = joined_label([a, b])
        a_color = colors[a]
        b_color = colors[b]
        colors[label] = ((pow(a_color, 2) + pow(b_color, 2)) / 2).sqrt()
    return colors


def propagate_colors(colors, order, data, transitive=False):
    # TODO: generalise?
    for i, current in enumerate(order[1:], 1):
        previous = order[i - 1 if transitive else 0]

        previous = data[data.group == previous]
        current = data[data.group == current]
        for cluster in current.cluster.unique():
            cluster_participants = set(current[current.cluster == cluster].participant)
            base = previous[previous.participant.isin(cluster_participants)]
            counts = base.groupby('cluster').participant.count().to_frame(name='participants_count')
            counts = counts[counts.participants_count > 0]
            counts['color'] = counts.index.map(colors)

            color = (
                    sum(pow(counts.color, 2) * counts.participants_count, Color()
            ) / sum(counts.participants_count)).sqrt()
            colors[cluster] = color
    return colors


def display_scales(scales, n=5):
    display(NeatNamespace({
        score: HorizontalNamespace(dict(enumerate(map(Color.from_hex, scale.palette(n)))))
        for scale, score in scales.items()
    }))


def display_scale(scale):
    display(HorizontalNamespace(dict(enumerate(map(Color.from_hex, scale.palette(5))))))


class Color:
    def __init__(self, components=(0, 0, 0, 0.0)):
        self.components = components

    @classmethod
    def from_hex(cls, components):
        return cls(to_rgba(components))

    def __add__(self, other):
        return Color([a + b for a, b in zip(self.components, other.components)])

    def __mul__(self, n: float):
        return Color([c * n for c in self.components])

    __rmul__ = __mul__

    def __truediv__(self, n: float):
        return Color([c / n for c in self.components])

    def __repr__(self):
        return f'<Color {self.to_hex()}>'

    def __pow__(self, n=2):
        return Color([c * c for c in self.components])

    def sqrt(self):
        return Color([sqrt(c) for c in self.components])

    def to_hex(self, alpha=False):
        return to_hex(self.components, keep_alpha=alpha)

    def _repr_html_(self):
        return f'<div style="background:{self.to_hex()}">{self.to_hex()}</div>'


class CustomPalette:

    def __init__(self, colors):
        self.colors = colors

    def palette(self, n):
        return self.colors

    def __repr__(self):
        return f'{self.colors}'
