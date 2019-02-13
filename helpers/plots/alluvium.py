from itertools import combinations, permutations

from pandas import Series
from scipy.stats import spearmanr
from sklearn.metrics import homogeneity_score
from tqdm import tqdm

from colorzero import Color as ColorZero


def rank_by_similarity(data, scoring_function=homogeneity_score):
    groups = data.group.unique()
    similarity_ranking = {}
    for a_name, b_name in combinations(groups, 2):
        a = data[data.group == a_name]
        b = data[data.group == b_name]
        common_participants = list(set(a.participant) & set(b.participant))
        available_participants = max([len(a.participant.unique()), len(b.participant.unique())])
        shared_ratio = len(common_participants) / available_participants
        similarity_ranking[frozenset({a_name, b_name})] = (
            shared_ratio * scoring_function(
                a.set_index('participant').loc[common_participants].cluster,
                b.set_index('participant').loc[common_participants].cluster
            )
        )
    return similarity_ranking


def suggest_groups_ordering(data):
    similarity_ranking = rank_by_similarity(data)
    groups = data.group.unique()

    ordering_ranking = {}
    for permutation in permutations(groups):
        similarity = 0
        for i in range(len(groups) - 1):
            a = permutation[i - 1]
            b = permutation[i]
            similarity += similarity_ranking[frozenset({a, b})]
        ordering_ranking[permutation] = similarity

    return Series(ordering_ranking).sort_values(ascending=False)


def determine_order_for_clusters_in_groups(data, ordered_groups, reference_group_order):
    """First group is used as a reference"""
    ordered_participants_ranks = []
    participants_order = {}
    rank = 0

    reference_group = data[data.group == ordered_groups[0]]
    for cluster in reference_group_order:
        for participant in reference_group[reference_group.cluster == cluster].participant:
            ordered_participants_ranks.append(rank)
            participants_order[participant] = rank
        rank += 1

    all_group_orders = []
    for group in tqdm(ordered_groups[1:]):
        ranked_permutations = {}
        group = data[data.group == group]
        all_clusters = group.cluster.unique()
        tested_clusters = [
            c
            for c in all_clusters
            if len(group[group.cluster == c].participant) > 10
        ]
        for permutation in permutations(tested_clusters):
            reference = []
            permutation_participants_ranks = []
            for rank, cluster in enumerate(permutation):
                for participant in group[group.cluster == cluster].participant:
                    reference.append(participants_order.get(participant, 0))
                    permutation_participants_ranks.append(rank)

            score = spearmanr(reference, permutation_participants_ranks).correlation

            ranked_permutations[permutation] = score

        chosen_group_order = [
            *Series(ranked_permutations).sort_values(ascending=False).index[0],
            *list(set(all_clusters) - set(tested_clusters))
        ]
        all_group_orders.extend(chosen_group_order)
    return all_group_orders


def suggest_contrastive_colors(scales, colors_generation_piepeline, groups):
    scales_ranking = {}
    for scale in tqdm(scales):
        colors = colors_generation_piepeline(scale=scale)
        total_difference = 0
        differences = []
        for group_name, group in groups.items():
            colors_in_group = [
                (cluster, colors[cluster])
                for cluster in group.cluster.unique()
            ]
            n_comb = sum(1 for _ in combinations(colors_in_group, 2))
            difference = sum(
                ColorZero.from_string(a.to_hex()).difference(ColorZero.from_string(b.to_hex()), method='ciede2000') * (
                    len(group[group.cluster == c_b]) * len(group[group.cluster == c_a])
                )
                for (c_a, a), (c_b, b) in combinations(colors_in_group, 2)
            ) / n_comb
            differences.append(difference)
            total_difference += pow(difference, 2)
        scales_ranking[scale] = max(differences)
    return scales_ranking
