#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import matplotlib

if "linux" in sys.platform:
    gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']
elif sys.platform == "darwin":
    try:
        import PySide2
        gui_env = ['Qt5Agg']
    except ImportError:
        gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']
else:
    pass

if sys.platform != "win32":
    for gui in gui_env:
        try:
            matplotlib.use(gui, warn=False, force=True)
            break
        except:
            continue


import sqlite3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


def plot_correlations(column_corr, column_pvalue, df):

    fig = plt.figure(figsize=(8, 8))
    fig.set_size_inches(8.27, 11.69)

    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_x_dist = plt.subplot(gs[0, :2], sharex=ax_main)
    ax_y_dist = plt.subplot(gs[1:3, 2], sharey=ax_main)

    ax_main.grid(linestyle='dashed')
    ax_x_dist.grid(linestyle='dashed')
    ax_y_dist.grid(linestyle='dashed')

    ax_main.set_axisbelow(True)
    ax_x_dist.set_axisbelow(True)
    ax_y_dist.set_axisbelow(True)

    max_pvalue = df[column_pvalue].max()
    bin_size_pvalue = max_pvalue / 10.0

    hb = ax_main.hexbin(x=column_corr, y=column_pvalue, data=df, gridsize=(40, 40), mincnt=1, extent=[-1, 1.0, 0, max_pvalue])
    ax_main.set(xlabel="Correlation coefficient (R)", ylabel="P-value",
                xticks=np.arange(-1, 1.1, 0.1), yticks=np.arange(0.0,  max_pvalue * 1.1, bin_size_pvalue))
    ax_main.ticklabel_format(style='sci', axis='y', scilimits=(0, 0), useMathText=True)

    bins = np.arange(-1, 1.05, 0.05)
    ax_x_dist.hist(x=column_corr, data=df, bins=bins, align='mid', color="lightblue")
    ax_x_dist.set(ylabel='Frequency', xlim=(-1.1, 1.1))
    ax_x_dist.axvline(0, color='k', linestyle='dashed', linewidth=1)
    ax_x_dist.tick_params(axis="x", labelsize=7.5)

    ax_xcum_dist = ax_x_dist.twinx()
    ax_xcum_dist.hist(x=column_corr, data=df, bins=bins, cumulative=True, histtype='step',
                      density=True, color='darkblue', align='mid')
    ax_xcum_dist.set(xlim=(-1.1, 1.1))
    ax_xcum_dist.tick_params(axis="y", colors='darkblue')
    ax_xcum_dist.set_ylabel('cumulative', color='darkblue')
    ax_xcum_dist.set(yticks=np.arange(0.0, 1.2, 0.2))

    bins = np.arange(0, max_pvalue + bin_size_pvalue, bin_size_pvalue)
    ax_y_dist.hist(x=column_pvalue, data=df, bins=bins, orientation='horizontal',
                   align='mid', color="lightblue")
    ax_y_dist.set(xlabel='Frequency')
    ax_ycum_dist = ax_y_dist.twiny()
    ax_ycum_dist.hist(x=column_pvalue, data=df, bins=bins, cumulative=True, histtype='step',
                      density=True, color='darkblue', align='mid', orientation='horizontal')
    ax_ycum_dist.tick_params(axis="x", colors='darkblue')
    ax_ycum_dist.set_xlabel('cumulative', color='darkblue')
    ax_ycum_dist.set(xticks=np.arange(0.0, 1.2, 0.2), ylim=(-bin_size_pvalue, max_pvalue * 1.1))

    #plt.setp(ax_x_dist.get_xticklabels(), visible=False)
    plt.setp(ax_y_dist.get_yticklabels(), visible=False)
    plt.setp(ax_x_dist.get_xticklabels(), rotation=90)
    plt.setp(ax_main.get_xticklabels(), rotation=90)
    plt.setp(ax_y_dist.get_xticklabels(), rotation=90)
    plt.setp(ax_ycum_dist.get_xticklabels(), rotation=90)

    fig.subplots_adjust(top=0.85, right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.4])

    cb = plt.colorbar(hb, cax=cbar_ax)
    cb.set_label('Frequency')

    return plt


def plot_annotations(column_ppm_error, column_adducts, df):

    fig = plt.figure()
    fig.set_size_inches(8.27, 11.69)

    gs = gridspec.GridSpec(5, 2, height_ratios=[1, 1, 5, 1, 1])

    ax_box = plt.subplot(gs[2])
    ax_hist = plt.subplot(gs[4], sharex=ax_box)
    ax_count = plt.subplot(gs[5])
    # ax = plt.subplot(gs[1])

    ppm_errors = df[column_ppm_error].dropna()

    sns.boxplot(x=ppm_errors, ax=ax_box)

    bin_size = 0.1
    bins = np.arange(np.floor(ppm_errors.min()) - bin_size, np.ceil(ppm_errors.max()) + bin_size, bin_size).round(3)
    ax_hist.hist(x=column_ppm_error, data=df, bins=bins, align='mid', color="lightblue")

    ax_hist.grid(False)
    ax_box.grid(False)

    std = ppm_errors.std()
    mean = ppm_errors.mean()
    median = ppm_errors.median()
    Q1 = ppm_errors.quantile(0.25)
    Q3 = ppm_errors.quantile(0.75)

    # Remove x axis name for the boxplot
    ax_box.set(xlabel="")
    # ax_box.set_xticks([])
    ax_box.set_title("Q1={}; median={}; Q3={}".format(round(Q1, 2), round(median, 2), round(Q3, 2)))

    ax_hist.set_title("mean={}; std={}".format(round(mean, 2), round(std, 2)))
    ax_hist.set(xlabel="Ppm error", ylabel="Frequency")

    sns.countplot(x=df[column_adducts].dropna(), ax=ax_count)
    ax_count.set(xlabel="Adduct", ylabel="Frequency")

    plt.setp(ax_box.get_xticklabels(), visible=False)
    plt.setp(ax_count.get_xticklabels(), rotation=90)

    return plt


def report(db, pdf_out, column_corr, column_pvalue, column_ppm_error, column_adducts):

    with PdfPages(pdf_out) as pdf:

        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        cursor.execute("""SELECT name FROM sqlite_master WHERE type='table';""")
        title = "Summary - BEAMSpy\n\n\n"
        for i, table in enumerate(cursor.fetchall()):
            if str(table[0]) == "groups":

                df = pd.read_sql_query("SELECT {}, {} FROM groups".format(column_corr, column_pvalue), conn)

                plt = plot_correlations(column_corr, column_pvalue, df)
                plt.suptitle('{}Grouping features'.format(title), fontsize=20)
                title = "\n\n\n"
                pdf.savefig(dpi=300)
                plt.close()

            elif table[0][0:10] == "compounds_":

                df = pd.read_sql_query("SELECT {}, {} FROM {}".format(column_ppm_error, column_adducts, table[0]), conn)

                plt = plot_annotations("ppm_error", "adduct", df)
                plt.suptitle('{}Compound Annotation\nDatabase: {}'.format(title, table[0].replace("compounds_", "")), fontsize=20)
                title = "\n\n\n"
                pdf.savefig(dpi=300)
                plt.close()
        conn.close()


# if __name__ == '__main__':
#
#     report("../tests/test_results/results_annotation.sqlite", "test_report_01.pdf",
#            "r_value", "p_value", "ppm_error", "adduct")
#     statinfo = os.stat("test_report_01.pdf")
#
#     report("../tests/test_results/results_pearson_all.sqlite", "test_report_02.pdf",
#            "r_value", "p_value", "ppm_error", "adduct")
