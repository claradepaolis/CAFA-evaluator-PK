import os
import numpy as np
import pandas as pd
import multiprocessing as mp
import random
from parser import obo_parser, gt_parser, pred_parser, gt_exclude_parser, update_toi
import logging
logging.getLogger(__name__).addHandler(logging.NullHandler())


# Return a mask for all the predictions (matrix) >= tau
def solidify_prediction(pred, tau):
    return pred >= tau


# computes the f metric for each precision and recall in the input arrays
def compute_f(pr, rc):
    n = 2 * pr * rc
    d = pr + rc
    return np.divide(n, d, out=np.zeros_like(n, dtype=float), where=d != 0)


def compute_s(ru, mi):
    return np.sqrt(ru**2 + mi**2)
    # return np.where(np.isnan(ru), mi, np.sqrt(ru + np.nan_to_num(mi)))


def compute_confusion_matrix(tau_arr, g, pred, toi, n_gt, ic_arr=None, B_ind = None):
    """
    Perform the evaluation at the matrix level for all tau thresholds
    The calculation is
    """
    # n, tp, fp, fn, pr, rc (fp = misinformation, fn = remaining uncertainty)
    metrics = np.zeros((len(tau_arr), 6), dtype='float')
    metrics_B_tau = {}
    metrics_B = []

    for i, tau in enumerate(tau_arr):

        # Filter predictions based on tau threshold
        p = solidify_prediction(pred.matrix[:, toi], tau)

        # Terms subsets
        intersection = np.logical_and(p, g)  # TP
        mis = np.logical_and(p, np.logical_not(g))  # FP, predicted but not in the ground truth
        remaining = np.logical_and(np.logical_not(p), g)  # FN, not predicted but in the ground truth

        # Weighted evaluation
        if ic_arr is not None:
            p = p * ic_arr[toi]
            intersection = intersection * ic_arr[toi]  # TP
            mis = mis * ic_arr[toi]  # FP, predicted but not in the ground truth
            remaining = remaining * ic_arr[toi]  # FN, not predicted but in the ground truth

        n_pred = p.sum(axis=1)  # TP + FP
        n_intersection = intersection.sum(axis=1)  # TP

        # Number of proteins with at least one term predicted with score >= tau
        metrics[i, 0] = (p.sum(axis=1) > 0).sum()

        # Sum of confusion matrices
        metrics[i, 1] = n_intersection.sum()  # TP
        metrics[i, 2] = mis.sum(axis=1).sum()  # FP
        metrics[i, 3] = remaining.sum(axis=1).sum()  # FN

        # Macro-averaging
        metrics[i, 4] = np.divide(n_intersection, n_pred, out=np.zeros_like(n_intersection, dtype='float'), where=n_pred > 0).sum()  # Precision
        metrics[i, 5] = np.divide(n_intersection, n_gt, out=np.zeros_like(n_gt, dtype='float'), where=n_gt > 0).sum()  # Recall

        if B_ind is not None:
            metrics_B_tau[i] = bootstrap(p, intersection, mis, remaining, n_gt, B_ind)
    #if B_ind is not None:
    #    metrics_B = get_metrics_B(metrics_B_tau)
    return metrics, metrics_B_tau

def compute_confusion_matrix_exclude(tau_arr, g_perprotein, pred, toi_perprotein, n_gt, ic_arr=None, B_ind = None):
    """
    Perform the evaluation at the matrix level for all tau thresholds
    The calculation is

    Here, g is the full ground truth matrix without filtering terms of interest (toi).
    Instead,
    """
    # n, tp, fp, fn, pr, rc (fp = misinformation, fn = remaining uncertainty)
    metrics = np.zeros((len(tau_arr), 6), dtype='float')
    metrics_B_tau = {}
    metrics_B = []

    for i, tau in enumerate(tau_arr):

        # Filter predictions based on tau threshold
        p_perprotein = [solidify_prediction(pred.matrix[p_idx, tois], tau) for p_idx, tois in enumerate(toi_perprotein)]

        # Terms subsets
        intersection = [np.logical_and(p_i, g_i) for p_i, g_i in zip(p_perprotein, g_perprotein)]  # TP
        mis = [np.logical_and(p_i, np.logical_not(g_i)) for p_i, g_i in zip(p_perprotein, g_perprotein)]  # FP, predicted but not in the ground truth
        remaining = [np.logical_and(np.logical_not(p_i), g_i) for p_i, g_i in zip(p_perprotein, g_perprotein)]  # FN, not predicted but in the ground truth

        # Weighted evaluation
        if ic_arr is not None:
            p_perprotein = [p_i * ic_arr[tois] for p_i, tois in zip(p_perprotein, toi_perprotein)]
            intersection = [inter * ic_arr[tois] for inter, tois in zip(intersection, toi_perprotein)]  # TP
            mis = [misinf * ic_arr[tois] for misinf, tois in zip(mis, toi_perprotein)]  # FP, predicted but not in the ground truth
            remaining = [rem * ic_arr[tois] for rem, tois in zip(remaining, toi_perprotein)]  # FN, not predicted but in the ground truth

        n_pred = np.array([p_i.sum() for p_i in p_perprotein])  # TP + FP
        n_intersection = np.array([inter.sum() for inter in intersection])  # TP

        # Number of proteins with at least one term predicted with score >= tau
        metrics[i, 0] = (n_pred > 0).sum()

        # Sum of confusion matrices
        metrics[i, 1] = n_intersection.sum()  # TP
        metrics[i, 2] = np.sum([m.sum() for m in mis])  # FP
        metrics[i, 3] = np.sum([r.sum() for r in remaining])  # FN

        # Macro-averaging
        metrics[i, 4] = np.divide(n_intersection, n_pred, out=np.zeros_like(n_intersection, dtype='float'), where=n_pred > 0).sum()  # Precision
        metrics[i, 5] = np.divide(n_intersection, n_gt, out=np.zeros_like(n_gt, dtype='float'), where=n_gt > 0).sum()  # Recall
        if B_ind is not None:
            metrics_B_tau[tau] = bootstrap_exclude(p_perprotein, intersection, mis, remaining, n_gt, B_ind)
    #if B_ind is not None:
    #    metrics_B = get_metrics_B(metrics_B_tau)

    return metrics, metrics_B_tau

# Input-> metrics_B_tau : a dict where thresholds are the keys, and a metrics array per threshold containing B rows, is in the values
# output-> metrics_B: a dict where a b index corresponding to each bootstrap round is the key and a metrics array containing one row per tau (threshold) is the output
def get_metrics_B(metrics_B_tau):
    taus = sorted(list(metrics_B_tau.keys()))
    B = len(metrics_B_tau[taus[0]]) #B = number of rows in the dict at the first key (threshold)
    metrics_B = {}
    columns = ["n", "tp", "fp", "fn", "pr", "rc"]
    for b in range(B):
        rows = []
        metrics_b = np.zeros((len(metrics_B_tau.keys()), 6), dtype='float')
        for i, tau in enumerate(taus):
            metrics_b[i] = metrics_B_tau[tau][b]
        metrics_B[b] = pd.DataFrame(metrics_b, columns=columns)
    return metrics_B
def bootstrap(p, intersection, mis, remaining, n_gt, B_ind):
    metrics_B_tau = np.zeros((len(B_ind), 6), dtype='float')
    for b, ind in enumerate(B_ind):
        n_gt_b = n_gt[ind]
        p_b = p[ind]
        intersection_b = intersection[ind]
        mis_b = mis[ind]
        remaining_b = remaining[ind]

        n_pred_b = p_b.sum(axis=1)  # TP + FP
        n_intersection_b = intersection_b.sum(axis=1)  # TP

        # Number of proteins with at least one term predicted with score >= tau
        metrics_B_tau[b, 0] = (p_b.sum(axis=1) > 0).sum()

        # Sum of confusion matrices
        metrics_B_tau[b, 1] = n_intersection_b.sum()  # TP
        metrics_B_tau[b, 2] = mis_b.sum(axis=1).sum()  # FP
        metrics_B_tau[b, 3] = remaining_b.sum(axis=1).sum()  # FN

        # Macro-averaging
        metrics_B_tau[b, 4] = np.divide(n_intersection_b, n_pred_b, out=np.zeros_like(n_intersection_b, dtype='float'),
                                  where=n_pred_b > 0).sum()  # Precision
        metrics_B_tau[b, 5] = np.divide(n_intersection_b, n_gt_b, out=np.zeros_like(n_gt_b, dtype='float'),
                                  where=n_gt_b > 0).sum()  # Recall

    return metrics_B_tau

def bootstrap_exclude(p_perprotein, intersection, mis, remaining, n_gt, B_ind):
    metrics_B_tau = np.zeros((len(B_ind), 6), dtype='float')
    for b, ind in enumerate(B_ind):
        n_gt_b = n_gt[ind]

        p_perprotein_b = [p_perprotein[p] for p in ind]
        intersection_b = [intersection[p] for p in ind]
        mis_b = [mis[p] for p in ind]
        remaining_b = [remaining[p] for p in ind]

        n_pred_b = np.array([p_i.sum() for p_i in p_perprotein_b])  # TP + FP
        n_intersection_b = np.array([inter.sum() for inter in intersection_b])  # TP

        # Sum of confusion matrices
        metrics_B_tau[b, 1] = n_intersection_b.sum()  # TP
        metrics_B_tau[b, 2] = np.sum([m.sum() for m in mis_b])  # FP
        metrics_B_tau[b, 3] = np.sum([r.sum() for r in remaining_b])  # FN

        # Macro-averaging
        metrics_B_tau[b, 4] = np.divide(n_intersection_b, n_pred_b, out=np.zeros_like(n_intersection_b, dtype='float'),
                                  where=n_pred_b > 0).sum()  # Precision
        metrics_B_tau[b, 5] = np.divide(n_intersection_b, n_gt_b, out=np.zeros_like(n_gt_b, dtype='float'),
                                  where=n_gt_b > 0).sum()  # Recall

    return metrics_B_tau
def compute_metrics(pred, gt, tau_arr, toi, gt_exclude=None, ic_arr=None, n_cpu=0, B_ind = None):
    """
    Takes the prediction and the ground truth and for each threshold in tau_arr
    calculates the confusion matrix and returns the coverage,
    precision, recall, remaining uncertainty and misinformation.
    Toi is the list of terms (indexes) to be considered
    """
    # Parallelization
    if n_cpu == 0:
        n_cpu = mp.cpu_count()

    columns = ["n", "tp", "fp", "fn", "pr", "rc"]
    g = gt.matrix[:, toi]

    if gt_exclude is not None:
        g_exclude = gt_exclude.matrix[:, toi]
        toi_perprotein = [np.setdiff1d(toi, gt_exclude.matrix[p, :].nonzero()[0], assume_unique=True) for p in
                          range(g.shape[0])]
        gt_perprotein = [gt.matrix[p_idx, tois] for p_idx, tois in enumerate(toi_perprotein)]
        # The number of GT annotations per proteins will change to exclude the set from g_exclude
        count_g = np.logical_and(np.logical_not(g_exclude), g)  # count terms in g only if they are not in exclude list
    else:
        count_g = g


    # Simple metrics
    if ic_arr is None:
        n_gt = count_g.sum(axis=1)
    # Weighted metrics
    else:
        n_gt = (count_g * ic_arr[toi]).sum(axis=1)

    if gt_exclude is None:
        arg_lists = [[tau_arr, g, pred, toi, n_gt, ic_arr, B_ind] for tau_arr in np.array_split(tau_arr, n_cpu)]
        with mp.Pool(processes=n_cpu) as pool:
            #metrics = np.concatenate(pool.starmap(compute_confusion_matrix, arg_lists), axis=0)
            results = pool.starmap(compute_confusion_matrix, arg_lists)
            metrics = [results[i][0] for i in range(len(results))]
            metrics = pd.DataFrame(np.concatenate(metrics), columns=columns)
            metrics_B_tau = {}
            if results[0][1]:  # If metrics from the bootstrapping were calculated
                for thread in range(len(results)):
                    for tau, metrics_b in results[thread][1].items():
                        metrics_B_tau[tau] = metrics_b
                metrics_B = get_metrics_B(metrics_B_tau)
    else:
        arg_lists = [[tau_arr, gt_perprotein, pred, toi_perprotein, n_gt, ic_arr, B_ind] for tau_arr in np.array_split(tau_arr, n_cpu)]
        with mp.Pool(processes=n_cpu) as pool:
            #metrics = np.concatenate(pool.starmap(compute_confusion_matrix_exclude, arg_lists), axis=0)
            results = pool.starmap(compute_confusion_matrix_exclude, arg_lists)
            metrics = [results[i][0] for i in range(len(results))]
            metrics = pd.DataFrame(np.concatenate(metrics), columns=columns)
            metrics_B_tau = {}
            if results[0][1]:  # If metrics from the bootstrapping were calculated
                for thread in range(len(results)):
                    for tau, metrics_b in results[thread][1].items():
                        metrics_B_tau[tau] = metrics_b
                metrics_B = get_metrics_B(metrics_B_tau)

    return metrics, metrics_B


def normalize(metrics, ns, tau_arr, ne, normalization):

    # Normalize columns
    for column in metrics.columns:
        if column != "n":
            # By default normalize by gt
            denominator = ne
            # Otherwise normalize by pred
            if normalization == 'pred' or (normalization == 'cafa' and column == "pr"):
                denominator = metrics["n"]
            metrics[column] = np.divide(metrics[column], denominator,
                                        out=np.zeros_like(metrics[column], dtype='float'),
                                        where=denominator > 0)

    metrics['ns'] = [ns] * len(tau_arr)
    metrics['tau'] = tau_arr
    metrics['cov'] = metrics['n'] / ne
    metrics['mi'] = metrics['fp']
    metrics['ru'] = metrics['fn']

    metrics['f'] = compute_f(metrics['pr'], metrics['rc'])
    metrics['s'] = compute_s(metrics['ru'], metrics['mi'])

    # Micro-average, calculation is based on the average of the confusion matrices
    metrics['pr_micro'] = np.divide(metrics['tp'], metrics['tp'] + metrics['fp'],
                                    out=np.zeros_like(metrics['tp'], dtype='float'),
                                    where=(metrics['tp'] + metrics['fp']) > 0)
    metrics['rc_micro'] = np.divide(metrics['tp'], metrics['tp'] + metrics['fn'],
                                    out=np.zeros_like(metrics['tp'], dtype='float'),
                                    where=(metrics['tp'] + metrics['fn']) > 0)
    metrics['f_micro'] = compute_f(metrics['pr_micro'], metrics['rc_micro'])

    return metrics


def evaluate_prediction(prediction, gt, ontologies, tau_arr, gt_exclude=None, normalization='cafa', n_cpu=0, B = 0, B_pct = 0):

    dfs = []
    dfs_w = []
    metrics_B_dfs = []
    metrics_B_w_dfs = []

    # Unweighted metrics
    for ns in prediction:
        if gt_exclude is None:
            exclude = None
            # number of proteins with positive predicitons
            num_pred_prots = (gt[ns].matrix[:, ontologies[ns].toi].sum(1) > 0).sum()
        else:
            exclude = gt_exclude[ns]
            # number of proteins with positive predicitons
            num_proteins = gt[ns].matrix.shape[0]
            toi_perprotein = [
                np.setdiff1d(ontologies[ns].toi, gt_exclude[ns].matrix[p, :].nonzero()[0], assume_unique=True) for p in
                range(num_proteins)]
            num_pred_prots = sum([gt[ns].matrix[p, toi_perprotein[p]].sum()>0 for p in range(num_proteins)])

        ne = np.full(len(tau_arr), num_pred_prots)

        # Generate B sets of indices
        B_ind = []
        N = len(gt[ns].ids)  # Number of proteins
        nB = 0
        if B and B_pct > 0:
            nB = round((B_pct / 100) * N)  # Number of proteins to be included in each bootstrap round
            for b in range(B):
                B_ind.append(random.choices(range(0, N), k=nB))

        metrics, metrics_B = compute_metrics(prediction[ns], gt[ns], tau_arr, ontologies[ns].toi, exclude, None, n_cpu, B_ind = B_ind)
        dfs.append(normalize(metrics, ns, tau_arr, ne, normalization))
        metrics_B_df = []
        if metrics_B:
            for b in metrics_B.keys():
                metrics_b = normalize(metrics_B[b], ns, tau_arr, ne, normalization)
                metrics_b["b"] = b
                metrics_B_df.append(metrics_b)
            metrics_B_dfs.append(pd.concat(metrics_B_df)) # Concats dfs for each b iteration together

        # Weighted metrics
        if ontologies[ns].ia is not None:
            exclude = None if gt_exclude is None else gt_exclude[ns]
            ne = np.full(len(tau_arr), gt[ns].matrix[:, ontologies[ns].toi_ia].shape[0])
            metrics_w, metrics_B_w = compute_metrics(prediction[ns], gt[ns], tau_arr, ontologies[ns].toi_ia, exclude, ontologies[ns].ia, n_cpu, B_ind = B_ind)
            dfs_w.append(normalize(metrics_w, ns, tau_arr, ne, normalization))
            metrics_B_w_df = []
            if metrics_B_w:
                for b in metrics_B_w.keys():
                    metrics_b_w = normalize(metrics_B_w[b], ns, tau_arr, ne, normalization)
                    metrics_b_w["b"] = b
                    metrics_B_w_df.append(metrics_b_w)
                metrics_B_w_dfs.append(pd.concat(metrics_B_w_df))  # Concats dfs for each b iteration together

    dfs = pd.concat(dfs) # Concats df from each ns

    # Merge weighted and unweighted dataframes
    if dfs_w:
        dfs_w = pd.concat(dfs_w)
        dfs = pd.merge(dfs, dfs_w, on=['ns', 'tau'], suffixes=('', '_w'))

    if metrics_B_dfs:
        metrics_B_dfs = pd.concat(metrics_B_dfs)
        if metrics_B_w_dfs:
            metrics_B_w_dfs = pd.concat(metrics_B_w_dfs)
            metrics_B_dfs = pd.merge(metrics_B_dfs, metrics_B_w_dfs, on=['ns', 'tau', 'b'], suffixes=('', '_w'))

    return dfs, metrics_B_dfs


def cafa_eval(obo_file, pred_dir, gt_file, ia=None, no_orphans=False, norm='cafa', prop='max',
              exclude=None, toi_file=None, max_terms=None, th_step=0.01, n_cpu=1, B = 0, B_pct = 50):

    # Tau array, used to compute metrics at different score thresholds
    tau_arr = np.arange(th_step, 1, th_step)

    # Parse the OBO file and creates a different graphs for each namespace
    ontologies = obo_parser(obo_file, ("is_a", "part_of"), ia, not no_orphans)
    if toi_file is not None:
        ontologies = update_toi(ontologies, toi_file)

    # Parse ground truth file
    gt = gt_parser(gt_file, ontologies)
    if exclude is not None:
        gt_exclude = gt_exclude_parser(exclude, gt, ontologies)
    else:
        gt_exclude = None

    # Set prediction files looking recursively in the prediction folder
    pred_folder = os.path.normpath(pred_dir) + "/"  # add the tailing "/"
    pred_files = []
    for root, dirs, files in os.walk(pred_folder):
        for file in files:
            pred_files.append(os.path.join(root, file))
    logging.debug("Prediction paths {}".format(pred_files))

    # Parse prediction files and perform evaluation
    dfs = []
    for file_name in pred_files:
        prediction = pred_parser(file_name, ontologies, gt, prop, max_terms)
        if not prediction:
            logging.warning("Prediction: {}, not evaluated".format(file_name))
        else:
            df_pred = evaluate_prediction(prediction, gt, ontologies, tau_arr, gt_exclude,
                                          normalization=norm, n_cpu=n_cpu, B = B, B_pct= B_pct)
            df_pred['filename'] = file_name.replace(pred_folder, '').replace('/', '_')
            dfs.append(df_pred)
            logging.info("Prediction: {}, evaluated".format(file_name))

    # Concatenate all dataframes and save them
    df = None
    dfs_best = {}
    if dfs:
        df = pd.concat(dfs)

        # Remove rows with no coverage
        df = df[df['cov'] > 0].reset_index(drop=True)
        df.set_index(['filename', 'ns', 'tau'], inplace=True)

        # Calculate the best index for each namespace and each evaluation metric
        for metric, cols in [('f', ['rc', 'pr']), ('f_w', ['rc_w', 'pr_w']), ('s', ['ru', 'mi']), ('f_micro', ['rc_micro', 'pr_micro']), ('f_micro_w', ['rc_micro_w', 'pr_micro_w'])]:
            if metric in df.columns:
                index_best = df.groupby(level=['filename', 'ns'])[metric].idxmax() if metric in ['f', 'f_w', 'f_micro', 'f_micro_w'] else df.groupby(['filename', 'ns'])[metric].idxmin()
                df_best = df.loc[index_best]
                if metric[-2:] != '_w':
                    df_best['cov_max'] = df.reset_index('tau').loc[[ele[:-1] for ele in index_best]].groupby(level=['filename', 'ns'])['cov'].max()
                else:
                    df_best['cov_max'] = df.reset_index('tau').loc[[ele[:-1] for ele in index_best]].groupby(level=['filename', 'ns'])['cov_w'].max()
                dfs_best[metric] = df_best
    else:
        logging.info("No predictions evaluated")

    return df, dfs_best


def write_results(df, dfs_best, out_dir='results', th_step=0.01):

    # Create output folder here in order to store the log file
    out_folder = os.path.normpath(out_dir) + "/"
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)

    # Set the number of decimals to write in the output files based on the threshold step size
    decimals = int(np.ceil(-np.log10(th_step))) + 1

    df.to_csv('{}/evaluation_all.tsv'.format(out_folder), float_format="%.{}f".format(decimals), sep="\t")

    for metric in dfs_best:
        dfs_best[metric].to_csv('{}/evaluation_best_{}.tsv'.format(out_folder, metric), float_format="%.{}f".format(decimals), sep="\t")
