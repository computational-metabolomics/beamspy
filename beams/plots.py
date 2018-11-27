#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import matplotlib

if "linux" in sys.platform:
    gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']
elif sys.platform == "darwin":
    try:
        import PyQt5
        gui_env = ['Qt5Agg']
    except ImportError:
        gui_env = ['TkAgg', 'GTKAgg', 'Qt5Agg', 'WXAgg']
else:
    pass

if sys.platform != "win32":
    for gui in gui_env:
        try:
            matplotlib.use(gui, warn=False, force=True)
            print(gui)
            break
        except:
            continue

import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sns


import sys
import sqlite3
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages


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

    hb = ax_main.hexbin(x=column_corr, y=column_pvalue, data=df, gridsize=(20, 20), mincnt=1, extent=[-1, 1.0, 0, 0.1])
    ax_main.set(xlabel="Correlation coefficient (R)", ylabel="P-value", xticks=np.arange(-1, 1.1, 0.1), yticks=np.arange(0.0,  0.105, 0.005))

    ax_x_dist.hist(x=column_corr, data=df, bins=40, align='mid', color="lightblue")
    ax_x_dist.set(ylabel='count', xlim=(-1, 1))
    ax_x_dist.axvline(0, color='k', linestyle='dashed', linewidth=1)

    ax_xcum_dist = ax_x_dist.twinx()
    ax_xcum_dist.hist(x=column_corr, data=df, bins=40, cumulative=True, histtype='step',
                      density=True, color='darkblue', align='mid')
    ax_xcum_dist.set(xlim=(-1, 1))
    ax_xcum_dist.tick_params(column_corr, colors='darkblue')
    ax_xcum_dist.set_ylabel('cumulative', color='darkblue')
    ax_xcum_dist.set(yticks=np.arange(0.0, 1.2, 0.2))

    ax_y_dist.hist(x=column_pvalue, data=df, bins=200, orientation='horizontal',
                   align='mid', color="lightblue")
    ax_y_dist.set(xlabel='count')
    ax_ycum_dist = ax_y_dist.twiny()
    ax_ycum_dist.hist(x=column_pvalue, data=df, bins=200, cumulative=True, histtype='step',
                      density=True, color='darkblue', align='mid', orientation='horizontal')
    ax_ycum_dist.tick_params(column_pvalue, colors='darkblue')
    ax_ycum_dist.set_xlabel('cumulative', color='darkblue')
    ax_ycum_dist.set(xticks=np.arange(0.0, 1.2, 0.2), ylim=(0, 0.1))

    plt.setp(ax_x_dist.get_xticklabels(), visible=False)
    plt.setp(ax_y_dist.get_yticklabels(), visible=False)

    plt.setp(ax_main.get_xticklabels(), rotation=90)
    plt.setp(ax_y_dist.get_xticklabels(), rotation=90)

    fig.subplots_adjust(top=0.85, right=0.85)
    cbar_ax = fig.add_axes([0.87, 0.15, 0.03, 0.4])

    cb = plt.colorbar(hb, cax=cbar_ax)
    cb.set_label('counts')

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

    sns.boxplot(ppm_errors, ax=ax_box)
    sns.distplot(ppm_errors, ax=ax_hist)

    ax_hist.grid(False)
    ax_box.grid(False)
    ax_hist.grid(False)

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
    ax_hist.set(xlabel="ppm error")

    sns.set(style="whitegrid")
    sns.countplot(df[column_adducts].dropna(), ax=ax_count)

    plt.setp(ax_box.get_xticklabels(), visible=False)

    return plt


def report(db, pdf_out, column_corr, column_pvalue, column_ppm_error, column_adducts):

    with PdfPages(pdf_out) as pdf:

        conn = sqlite3.connect(db)
        cursor = conn.cursor()
        cursor.execute("""SELECT name FROM sqlite_master WHERE type='table';""")

        for table in cursor.fetchall():
            if str(table[0]) == "groups":

                df = pd.read_sql_query("SELECT {}, {} FROM groups".format(column_corr, column_pvalue), conn)

                plt = plot_correlations(column_corr, column_pvalue, df)
                plt.suptitle('Summary - BEAMS\n\n\nGrouping features', fontsize=20)
                pdf.savefig(dpi=300)
                plt.close()

            elif table[0][0:10] == "compounds_":

                df = pd.read_sql_query("SELECT {}, {} FROM {}".format(column_ppm_error, column_adducts, table[0]), conn)

                plt = plot_annotations("ppm_error", "adduct", df)
                plt.suptitle('Annotation\n\n{}'.format(table[0].replace("compounds_", "")), fontsize=20)
                pdf.savefig(dpi=300)
                plt.close()
        conn.close()
