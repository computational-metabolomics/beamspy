#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
from collections import OrderedDict
import numpy as np
import pandas as pd
from utils import to_test_data, create_matrix
from beams.statistics import correlation_coefficients, correlation_graphs
from beams.in_out import read_peaklist


class StatisticsTestCase(unittest.TestCase):

    def setUp(self):
        self.df_peaklist = read_peaklist(to_test_data("peaklist_dims_pos_theoretical.txt"))
        self.df_matrix = create_matrix(self.df_peaklist, 50)
        self.df_matrix = self.df_matrix.set_index("name")
        self.df_matrix.ix["336_9720324"] = self.df_matrix.ix["168_9896544"] * 0.5
        self.df_matrix.ix["504_9544104"] = self.df_matrix.ix["168_9896544"] * 0.4
        self.df_matrix.ix["492_0604104"] = self.df_matrix.ix["493_0637654"] * 1.1
        self.df_matrix.ix["230_9901644"] = self.df_matrix.ix["215_0162264"] * 2.0
        self.df_matrix = self.df_matrix.reset_index(drop=False)
        self.df_peaklist = self.df_peaklist.assign(intensity=pd.Series(self.df_matrix.median(axis=1, skipna=True).values))
        self.df = pd.merge(self.df_peaklist, self.df_matrix, how='left', left_on="name", right_on="name")

    def test_correlation_coefficients(self):
        df_coeffs_comp = pd.DataFrame({"name_a": ["168_9896544", "168_9896544", "215_0162264", "336_9720324", "492_0604104"],
                            "name_b": ["336_9720324", "504_9544104", "230_9901644", "504_9544104", "493_0637654"],
                            "r_value": [np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0)],
                            "p_value": [np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)]}, columns=["name_a", "name_b", "r_value", "p_value"])
        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=None, method="pearson", block=5000, ncpus=None)
        self.assertTrue(df_coeffs.equals(df_coeffs_comp))

        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=None, method="spearman", block=5000, ncpus=None)
        self.assertTrue(df_coeffs.equals(df_coeffs_comp))

        df_coeffs = correlation_coefficients(self.df, max_rt_diff=50000.0, coeff_thres=0.0, pvalue_thres=None, method="pearson", block=5000, ncpus=None)
        self.assertEqual(df_coeffs.shape, (136, 4))

    def test_correlation_graphs(self):
        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.0, pvalue_thres=None, method="pearson", block=5000, ncpus=None)
        graph = correlation_graphs(df_coeffs, self.df)

        # print list(graph.nodes(data=True))[0]
        # print list(graph.nodes(data=True))[-1]
        # print list(graph.edges(data=True))[0]
        # print list(graph.edges(data=True))[-1]

        self.assertEqual(list(graph.nodes(data=True))[0],
                         ('126_9792044', {'rt': 0.0, 'intensity': 481.78681200439746, 'mz': 126.9792044}))
        self.assertEqual(list(graph.nodes(data=True))[-1],
                         ('550_0658904', {'rt': 0.0, 'intensity': 512.2308133856002, 'mz': 550.0658904000001}))
        self.assertEqual(list(graph.edges(data=True))[0],
                         ('126_9792044', '135_0288014', OrderedDict([('rtdiff', 0.0), ('mzdiff', 8.049596999999991), ('rvalue', -0.02), ('pvalue', 0.8764648560841006)])))
        self.assertEqual(list(graph.edges(data=True))[-1],
                         ('504_9544104', '550_0658904', OrderedDict([('rtdiff', 0.0), ('mzdiff', 45.11148000000014), ('rvalue', 0.06), ('pvalue', 0.6882762403055659)])))

        #nx.write_gml(graph, to_test_data("graph_simple.gml"))
        #graph = nx.read_gml(to_test_data("graph_simple.gml"))


if __name__ == '__main__':
    unittest.main()