#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import scipy.stats
from multiprocessing import Pool, cpu_count
import networkx as nx
import pandas as pd
import tqdm


def _cc_pp_(pairs, method, ncpus):

    coeffs = []

    pool = Pool(ncpus)

    if len(pairs) > ncpus > 1:
        pairs = [pairs[i: i + int(len(pairs) / (ncpus - 1))] for i in range(0, len(pairs), int(len(pairs) / (ncpus - 1)))]
    else:
        pairs = [pairs]

    if method == "pearson":
        results = pool.map(_pearsonr, pairs)
    elif method == "spearman":
        results = pool.map(_spearmanr, pairs)
    else:
        raise ValueError("Method {} does not exist".format(method))

    pool.close()
    pool.join()

    for result in results:
        coeffs.extend(result)

    return coeffs


def _pearsonr(pairs):
    temp = []
    for pair in pairs:
        out = scipy.stats.pearsonr(pair[0], pair[1])
        temp.append([out[0], out[1]])
    return temp


def _spearmanr(pairs):
    temp = []
    for pair in pairs:
        out = scipy.stats.spearmanr(pair[0], pair[1])
        temp.append([out[0], out[1]])
    return temp


def correlation_coefficients(df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=0.05, method="pearson", block=5000, ncpus=None):
    if ["name", "mz", "rt"] == list(df.columns.values[0:3]):
        ncols = 4
        df = df.sort_values(['rt', 'mz']).reset_index(drop=True)
    elif "mz" == df.columns.values[0] and "rt" not in df.columns.values:
        ncols = 1
    else:
        raise ValueError("Incorrect column names: [name, mz, rt] or [mz, intensity]")

    if ncpus is None:
        ncpus = cpu_count()
        if ncpus > 1:
            ncpus -= 1

    column_names = ["name_a", "name_b", "r_value", "p_value"]
    df_coeffs = pd.DataFrame(columns=column_names)

    pairs, peaks = [], []
    n = len(df.iloc[:, 0])

    if n >= 100:
        disable_tqdm = False
    else:
        disable_tqdm = True

    for i in tqdm.trange(n, disable=disable_tqdm):

        intens_i = df.iloc[i, ncols:].values

        if pd.notnull(intens_i).sum() < 4:
            continue

        for j in range(i + 1, n):

            if max_rt_diff is not None:
                rt_diff = abs(float(df.loc[j, "rt"] - df.loc[i, "rt"]))
            else:
                rt_diff = 0.0  # Direct Infusion - no retention time available

            if rt_diff <= max_rt_diff and max_rt_diff is not None:

                intens_j = df.iloc[j, ncols:].values
                nas = np.logical_or(pd.isnull(intens_i), pd.isnull(intens_j))
                intens_filt_i, intens_filt_j = intens_i[~nas], intens_j[~nas]

                if len(intens_filt_i) > 3 and len(intens_filt_j) > 3:
                    peaks.append([df.iloc[i, 0], df.iloc[j, 0], rt_diff])
                    pairs.append([intens_filt_i, intens_filt_j])
                    if len(pairs) == block * ncpus:
                        #print("Calculating correlations for {} pairs (subset)".format(len(pairs)))
                        coeffs = _cc_pp_(pairs, method, ncpus)
                        for k in range(len(coeffs)):
                            if abs(coeffs[k][0]) > coeff_thres and (abs(coeffs[k][1]) < pvalue_thres or pvalue_thres is None):
                                s = pd.Series([peaks[k][0], peaks[k][1], round(coeffs[k][0], 2), coeffs[k][1]], index=column_names)
                                df_coeffs = df_coeffs.append(s, ignore_index=True) # pandas append
                        pairs, peaks = [], []
            else:
                break

    if len(pairs) > 0:
        #print("Calculating correlations for {} pairs (subset)".format(len(pairs)))
        coeffs = _cc_pp_(pairs, method, ncpus)
        for k in range(len(coeffs)):
            if abs(coeffs[k][0]) > coeff_thres and (abs(coeffs[k][1]) < pvalue_thres or pvalue_thres is None):
                s = pd.Series([peaks[k][0], peaks[k][1], round(coeffs[k][0], 2), coeffs[k][1]], index=column_names)
                df_coeffs = df_coeffs.append(s, ignore_index=True)
    return df_coeffs


def correlation_graphs(df_coeffs, df):
    df_coeffs = df_coeffs.merge(df[["name", "mz", "intensity", "rt"]], how='left', left_on=['name_a'], right_on=['name'])
    df_coeffs = df_coeffs.merge(df[["name", "mz", "intensity", "rt"]], how='left', left_on=['name_b'], right_on=['name'])
    from decimal import Decimal
    graphs = nx.OrderedDiGraph()
    for index, row in df_coeffs.iterrows():
        graphs.add_node(str(row["name_a"]), mz=row["mz_x"], intensity=row["intensity_x"], rt=row["rt_x"])
        graphs.add_node(str(row["name_b"]), mz=row["mz_y"], intensity=row["intensity_y"], rt=row["rt_y"])

        mz_diff = row["mz_x"] - row["mz_y"]

        if mz_diff < 0:
            graphs.add_edge(str(row["name_a"]), str(row["name_b"]), rvalue=row["r_value"], pvalue=row["p_value"],
                            mzdiff=abs(row["mz_x"]-row["mz_y"]), rtdiff=abs(row["rt_x"]-row["rt_y"]))
        else:
            graphs.add_edge(str(row["name_b"]), str(row["name_a"]), rvalue=row["r_value"], pvalue=row["p_value"],
                            mzdiff=abs(row["mz_x"] - row["mz_y"]), rtdiff=abs(row["rt_x"] - row["rt_y"]))
    return graphs
