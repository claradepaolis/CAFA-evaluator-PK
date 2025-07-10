import numpy as np


def test_norm_metric(metric, name=None):
    """
    Test if metric is between 0 and 1
    :param metric: array of metric to test
    :param name: string of metric name for error message if necessary
    """
    name = 'Metric' if name is None else name

    if not (metric.max() <= 1.0 and metric.min() >= 0):
        raise Exception(f'{name} error: max: {metric.max():.4f} min: {metric.min():.4f}')


def test_intersection(intersection, pred_counts, true_counts):
    """
    Test if intersection of counts is always smaller than counts from one method
    :param intersection: array of positive counts in intersection
    :param pred_counts: array of positive counts from prediction
    :param true_counts: array of positive counts from ground truth
    """
    if not np.all(pred_counts >= intersection):
        raise Exception(f'Count of positive prediction smaller than intersection')

    if not np.all(true_counts >= intersection):
        raise Exception(f'Count of positive annotations smaller than intersection')
