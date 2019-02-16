from collections import defaultdict
from functools import partial
from math import sqrt
from statistics import mean
from traceback import print_exc
from typing import Dict

import numpy as np
from pandas import DataFrame, Series
from sklearn.metrics import mean_squared_error

from helpers.source import source_for_table
from helpers.r import r_ks_test
from helpers import on_division_by_zero

from .scores_models import ProcessedScores
from . import calculation_utilities as calc


class MetricsManager:
    registry: Dict[str, Dict[str, 'EvaluationMetric']] = defaultdict(dict)

    def metrics_by_objective(self, convert_to=dict, exclude=None):
        grouped_by_objective = defaultdict(lambda: defaultdict(set))

        for category, metrics in self.registry.items():
            for metric in metrics.values():
                group = grouped_by_objective[metric.objective]
                group[category].add(metric.name)

        return convert_to({
            objective: convert_to(metrics_by_category)
            for objective, metrics_by_category in grouped_by_objective.items()
        })

    def defined_metrics_table(self):
        defined_metrics = []

        for category, metrics in self.registry.items():
            for metric in metrics.values():
                metric_data = {
                    'category': category,
                    'objective': metric.objective,
                    'code': source_for_table(metric.function, trim_first_decorator=True, type_hints=False),
                    'name': metric.name,
                    'combine': metric.combine.__name__
                }
                defined_metrics.append(metric_data)

        defined_metrics = DataFrame(defined_metrics).set_index(['category', 'name'])
        return defined_metrics

    def best_scores(self, category, name, scores) -> Series:
        objective_functions = {
            None: lambda column: np.nan,
            'maximize': max,
            'minimize': min
        }
        try:
            metric = self.registry[category][name]
            choose_best = objective_functions[metric.objective]
        except KeyError:
            choose_best = objective_functions[None]

        return scores == choose_best(scores)


metrics_manager = MetricsManager()


class EvaluationMetric:

    all_metrics_categories = {'overall', 'indications', 'contraindications', 'controls'}

    def __init__(self, function, category='overall', name=None, objective='maximize', combine=mean):
        assert objective in {'maximize', 'minimize', None}
        assert category in self.all_metrics_categories

        if not name:
            name = function.__name__.replace('_', ' ').title()

        # set up the metadata
        self.category = category
        self.objective = objective
        self.name = name
        self.function = function
        self.combine = combine

        # add to registry
        metrics_manager.registry[category][self.name] = self

    def __call__(self, score: ProcessedScores):
        try:
            return self.function(score)
        except Exception as e:
            print(f'Evaluation metric {self.category}:{self.name} failed:')
            print_exc()
            print(f'Returning NaN')
            return np.nan


def evaluation_metric(category='overall', name=None, objective='maximize', combine=mean):

    def decorator(metric_function):
        return EvaluationMetric(metric_function, category=category, objective=objective, name=name, combine=combine)

    return decorator


indications_metric = partial(evaluation_metric, category='indications')
contraindications_metric = partial(evaluation_metric, category='contraindications')
controls_metric = partial(evaluation_metric, category='controls')


# Precision, recall

@indications_metric(objective=None)
def precision(scores: ProcessedScores):
    return calc.precision(true_positives=scores.top.indications, all_selected=scores.top.all)


@indications_metric(objective=None)
def recall(scores: ProcessedScores):
    return calc.recall(true_positives=scores.top.indications, all_positives=scores.indications)


# F1 Score

@indications_metric()
def f1_score(scores: ProcessedScores):
    selected = scores.top
    return calc.f1(
        calc.precision(true_positives=selected.indications, all_selected=selected.all),
        calc.recall(true_positives=selected.indications, all_positives=scores.indications)
    )


@contraindications_metric(objective='minimize')
def f1_score(scores: ProcessedScores):
    selected = scores.top
    return calc.f1(
        calc.precision(true_positives=selected.contraindications, all_selected=selected.all),
        calc.recall(true_positives=selected.contraindications, all_positives=scores.contraindications)
    )


# Mean Square Error

@evaluation_metric(objective='minimize', name='RMSE')
def rmse(scores: ProcessedScores):
    results = scores.vector_overall
    return (
        sqrt(mean_squared_error(results.expected, results.observed))
        if not results.observed.dropna().empty else
        np.nan
    )


# Means comparison

@controls_metric()
def is_mean_better(scores: ProcessedScores):
    return scores.indications.mean() > scores.controls.mean()


@contraindications_metric()
def is_mean_better(scores: ProcessedScores):
    return scores.indications.mean() > scores.contraindications.mean()


# Means

@indications_metric(objective=None)
def mean(scores: ProcessedScores):
    return scores.indications.mean()


@contraindications_metric(objective=None)
def mean(scores: ProcessedScores):
    return scores.contraindications.mean()


@controls_metric(objective=None)
def mean(scores: ProcessedScores):
    return scores.controls.mean()


# Kolmogorov–Smirnov

@controls_metric(objective='minimize', name='KS p-value', combine=calc.fisher_method)
def ks_p(scores: ProcessedScores):
    # From R docs:
    # "Thus in the two-sample case alternative = "greater" includes distributions for which x is
    # stochastically *smaller* than y (the CDF of x lies above and hence to the left of that for y),
    # in contrast to t.test or wilcox.test."
    ks = r_ks_test(scores.indications, scores.controls, alternative='less')
    return ks['p.value'][0]


@contraindications_metric(objective='minimize', name='KS p-value', combine=calc.fisher_method)
def ks_p(scores: ProcessedScores):
    # See ks_p controls metric for explanation of alternative='less'
    ks = r_ks_test(scores.indications, scores.contraindications, alternative='less')
    return ks['p.value'][0]


# ROC AUC

@controls_metric(name='AUC ROC')
def roc(scores: ProcessedScores):
    return calc.generalized_roc_auc_score(scores.vector_controls)


@contraindications_metric(name='AUC ROC')
def roc(scores: ProcessedScores):
    return calc.generalized_roc_auc_score(scores.vector_contraindications)


@controls_metric(name='AUC ROC classification')
def roc_binary(scores: ProcessedScores):
    return calc.generalized_roc_auc_score(scores.vector_controls_binary)


@contraindications_metric(name='AUC ROC classification')
def roc_binary(scores: ProcessedScores):
    return calc.generalized_roc_auc_score(scores.vector_contraindications_binary)


@evaluation_metric(name='AUC ROC classification')
def roc_binary(scores: ProcessedScores):
    return calc.generalized_roc_auc_score(scores.vector_overall_binary)


# Ratio prioritized

@contraindications_metric()
@on_division_by_zero(fill_with=np.nan)
def indications_prioritized(scores: ProcessedScores):
    return (
        len(scores.indications[scores.indications > scores.contraindications.max()])
        /
        len(scores.contraindications)
    )


@controls_metric()
@on_division_by_zero(fill_with=np.nan)
def indications_prioritized(scores: ProcessedScores):
    return (
        len(scores.indications[scores.indications > scores.controls.max()])
        /
        len(scores.controls)
    )


# Normalized means difference

@contraindications_metric()
def normalized_means_difference(scores: ProcessedScores):
    return calc.normalized_means_difference(scores.indications, scores.contraindications)


@controls_metric()
def normalized_means_difference(scores: ProcessedScores):
    return calc.normalized_means_difference(scores.indications, scores.controls)


# metrics from Cheng 2014

@indications_metric(name='AUC0.1')
def partial_retrieval_auc_01(scores: ProcessedScores):
    """partial retrieval area under the ROC curve (AUC0.1) at false positive rate 0.1"""
    return calc.generalized_roc_auc_score(scores.vector_indications_over_non_indications, max_fpr=0.1)


@indications_metric(name='AUC0.01')
def partial_retrieval_auc_001(scores: ProcessedScores):
    """partial retrieval area under the ROC curve (AUC0.01) at false positive rate 0.01"""
    return calc.generalized_roc_auc_score(scores.vector_indications_over_non_indications, max_fpr=0.01)
