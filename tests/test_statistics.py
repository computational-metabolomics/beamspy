#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest

import numpy as np
import pandas as pd

from beamspy.in_out import combine_peaklist_matrix
from beamspy.statistics import correlation_coefficients, correlation_graphs
from tests.utils import to_test_data


class StatisticsTestCase(unittest.TestCase):

    def setUp(self):
        self.df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_lcms_theoretical.txt"))

    def test_correlation_coefficients(self):
        df_coeffs_comp = pd.DataFrame({"name_a": ["M169T120", "M169T120", "M337T121", "M215T170", "M492T190"],
                            "name_b": ["M337T121", "M505T122", "M505T122", "M231T174", "M493T192"],
                            "r_value": [np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0)],
                            "p_value": [np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(5.85415087865495e-157)]}, columns=["name_a", "name_b", "r_value", "p_value"])
        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)
        pd.testing.assert_frame_equal(df_coeffs, df_coeffs_comp)

        df_coeffs_comp = pd.DataFrame({"name_a": ["M169T120", "M169T120", "M337T121", "M215T170", "M492T190"],
                            "name_b": ["M337T121", "M505T122", "M505T122", "M231T174", "M493T192"],
                            "r_value": [np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0), np.float64(1.0)],
                            "p_value": [np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0), np.float64(0.0)]}, columns=["name_a", "name_b", "r_value", "p_value"])

        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="spearman", block=5000, ncpus=None)
        pd.testing.assert_frame_equal(df_coeffs, df_coeffs_comp)

        df_coeffs = correlation_coefficients(self.df, max_rt_diff=50000.0, coeff_thres=0.0, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)
        self.assertEqual(df_coeffs.shape, (136, 4))

    def test_correlation_graphs(self):
        df_coeffs = correlation_coefficients(self.df, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)
        graph = correlation_graphs(df_coeffs, self.df)

        n0 = list(graph.nodes(data=True))[0]
        n1 = list(graph.nodes(data=True))[-1]

        e0 = list(graph.edges(data=True))[0]
        e1 = list(graph.edges(data=True))[-1]

        # order is different between python 2 and 3
        np.testing.assert_almost_equal([n0[1]["mz"], n0[1]["intensity"], n0[1]["rt"]], [168.989654, 520.0, 120.0])
        np.testing.assert_almost_equal([n1[1]["mz"], n1[1]["intensity"], n1[1]["rt"]], [493.063765, 163.33, 192.5])
        np.testing.assert_almost_equal([e0[2]["rvalue"], e0[2]["pvalue"], e0[2]["mzdiff"], e0[2]["rtdiff"]], [1.0, 0.0, 167.982378, 1.0])
        np.testing.assert_almost_equal([e1[2]["rvalue"], e1[2]["pvalue"], e1[2]["mzdiff"], e1[2]["rtdiff"]], [1.0, 0.0, 1.003355, 2.5])

        #nx.write_gml(graph, to_test_data("graph_simple.gml"))
        #graph = nx.read_gml(to_test_data("graph_simple.gml"))


if __name__ == '__main__':
    unittest.main()