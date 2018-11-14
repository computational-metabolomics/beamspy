#!/usr/bin/python
# -*- coding: utf-8 -*-

from sys import platform
import matplotlib

if platform != "win32":
    gui_env = ['Qt5Agg', 'TkAgg']
    for gui in gui_env:
        try:
            matplotlib.use(gui, warn=False, force=True)
            break
        except:
            continue

import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


def report(df, column_ppm_error, column_adducts, fn_pdf):

    fig = plt.figure()

    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 5])

    ax_box = plt.subplot(gs[0])
    ax_hist = plt.subplot(gs[2], sharex=ax_box)
    ax_count = plt.subplot(gs[3])
    #ax = plt.subplot(gs[1])

    ppm_errors = df[column_ppm_error].dropna()

    sns.boxplot(ppm_errors, ax=ax_box)
    sns.distplot(ppm_errors, ax=ax_hist)

    std = ppm_errors.std()
    mean = ppm_errors.mean()
    median = ppm_errors.median()
    Q1 = ppm_errors.quantile(0.25)
    Q3 = ppm_errors.quantile(0.75)

    # Remove x axis name for the boxplot
    ax_box.set(xlabel="")
    #ax_box.set_xticks([])
    ax_box.set_title("Q1={}; median={}; Q3={}".format(round(Q1, 2), round(median, 2), round(Q3, 2)))

    ax_hist.set_title("mean={}; std={}".format(round(mean, 2), round(std, 2)))
    ax_hist.set(xlabel="ppm error")

    sns.set(style="whitegrid")
    sns.countplot(df[column_adducts].dropna(), ax=ax_count)

    plt.setp(ax_box.get_xticklabels(), visible=False)

    fig.suptitle('Summary - BEAMS', fontsize=20)
    fig.set_size_inches(11.69, 8.27)
    fig.savefig(fn_pdf, format="pdf")
