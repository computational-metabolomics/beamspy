#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
from collections import OrderedDict

import numpy as np
import pandas as pd

from beams.in_out import combine_peaklist_matrix
from beams.statistics import correlation_coefficients, correlation_graphs
from tests.utils import to_test_data


class StatisticsTestCase(unittest.TestCase):

    def setUp(self):
        self.df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_theoretical.txt"))

    def test_correlation_coefficients(self):
        df_coeffs_comp = pd.DataFrame({"name_a": ["M169T120", "M169T120", "M337T121", "M215T170", "M492T190"],
                            "name_b": ["M337T121", "M505T122", "M505T122", "M231T174", "M493T192"],
                            "r_value": [np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0)],
                            "p_value": [np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(5.854150141280045e-157)]}, columns=["name_a", "name_b", "r_value", "p_value"])
        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)
        pd.testing.assert_frame_equal(df_coeffs, df_coeffs_comp, check_exact=True)

        df_coeffs_comp = pd.DataFrame({"name_a": ["M169T120", "M169T120", "M337T121", "M215T170", "M492T190"],
                            "name_b": ["M337T121", "M505T122", "M505T122", "M231T174", "M493T192"],
                            "r_value": [np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0)],
                            "p_value": [np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)]}, columns=["name_a", "name_b", "r_value", "p_value"])

        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="spearman", block=5000, ncpus=None)
        pd.testing.assert_frame_equal(df_coeffs, df_coeffs_comp, check_exact=True)

        df_coeffs = correlation_coefficients(self.df, max_rt_diff=50000.0, coeff_thres=0.0, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)
        self.assertEqual(df_coeffs.shape, (136, 4))

    def test_correlation_graphs(self):
        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)
        graph = correlation_graphs(df_coeffs, self.df)

        # print list(graph.nodes(data=True))[0]
        # print list(graph.nodes(data=True))[-1]
        # print list(graph.edges(data=True))[0]
        # print list(graph.edges(data=True))[-1]

        self.assertEqual(list(graph.nodes(data=True))[0],
                         ('M169T120', {'rt': 120.0, 'intensity': 520.0, 'mz': 168.9896544}))
        self.assertEqual(list(graph.nodes(data=True))[-1],
                         ('M493T192', {'rt': 192.5, 'intensity': 163.33, 'mz': 493.0637654}))
        self.assertEqual(list(graph.edges(data=True))[0],
                         ('M169T120', 'M337T121', OrderedDict([('rtdiff', 1.0), ('mzdiff', 167.98237799999993), ('rvalue', 1.0), ('pvalue', 0.0)])))
        self.assertEqual(list(graph.edges(data=True))[-1],
                         ('M492T190', 'M493T192', OrderedDict([('rtdiff', 2.5), ('mzdiff',  1.0033550000001128), ('rvalue', 1.0), ('pvalue', 5.854150141280045e-157)])))

        #nx.write_gml(graph, to_test_data("graph_simple.gml"))
        #graph = nx.read_gml(to_test_data("graph_simple.gml"))


if __name__ == '__main__':
    unittest.main()