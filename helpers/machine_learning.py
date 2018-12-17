from pandas import Series
from sklearn.metrics import roc_curve, auc
from sklearn.model_selection import StratifiedKFold, ShuffleSplit
import numpy as np
from scipy import interp
import matplotlib.pyplot as plt

from sklearn import metrics
from sklearn.preprocessing import LabelEncoder


def transform_labels(y) -> Series:
    if type(next(iter(y))) is str:
        le = LabelEncoder()
        le.fit(y)
        y = le.transform(y)
    return Series(y)


def calc_auc(clf, test_x, test_y):
    y_pred = clf.predict(test_x)
    return metrics.roc_auc_score(
        transform_labels(test_y),
        transform_labels(y_pred.tolist())
    )


def roc_plot(classifier, X, y, n_splits=3, title='', labeller=None):
    cv = StratifiedKFold(n_splits=n_splits)
    #if labeller:
    #    y = [labeller(i) for i in y]
    y = transform_labels(y)
    #cv = ShuffleSplit(n_splits=n_splits)
    tprs = []
    aucs = []
    mean_fpr = np.linspace(0, 1, 100)

    i = 0
    for train, test in cv.split(X, y):
        probas_ = classifier.fit(X.iloc[train], y.iloc[train]).predict_proba(X.iloc[test])
        # Compute ROC curve and area the curve
        fpr, tpr, thresholds = roc_curve(y.iloc[test], probas_[:, 1])
        tprs.append(interp(mean_fpr, fpr, tpr))
        tprs[-1][0] = 0.0
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
        plt.plot(fpr, tpr, lw=1, alpha=0.3,
                 label='ROC fold %d (AUC = %0.2f)' % (i, roc_auc))

        i += 1
    plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='r',
             label='Chance', alpha=.8)

    mean_tpr = np.mean(tprs, axis=0)
    mean_tpr[-1] = 1.0
    mean_auc = auc(mean_fpr, mean_tpr)
    std_auc = np.std(aucs)
    plt.plot(mean_fpr, mean_tpr, color='b',
             label=r'Mean ROC (AUC = %0.2f $\pm$ %0.2f)' % (mean_auc, std_auc),
             lw=2, alpha=.8)

    std_tpr = np.std(tprs, axis=0)
    tprs_upper = np.minimum(mean_tpr + std_tpr, 1)
    tprs_lower = np.maximum(mean_tpr - std_tpr, 0)
    plt.fill_between(mean_fpr, tprs_lower, tprs_upper, color='grey', alpha=.2,
                     label=r'$\pm$ 1 std. dev.')

    plt.xlim([-0.05, 1.05])
    plt.ylim([-0.05, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver operating characteristic ' + title)
    plt.legend(loc="lower right")
    plt.show()
    return plt