#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import numpy as np
from tests.utils import to_test_data, to_test_results, sqlite_records
from beamspy.in_out import combine_peaklist_matrix
from beamspy.grouping import group_features


class GroupFeaturesTestCase(unittest.TestCase):

    def setUp(self):
        self.df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_lcms_theoretical.txt"))

    def test_group_features(self):
        fn_sql = "results_pearson.sqlite"
        db_out = to_test_results(fn_sql)
        group_features(self.df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)

        records = sqlite_records(to_test_results(fn_sql), "groups")
        records_comp = sqlite_records(to_test_data(fn_sql), "groups")
        for i in range(len(records)):
            self.assertEqual(records[i][0:6], records_comp[i][0:6])
            np.testing.assert_almost_equal(records[i][6:], records_comp[i][6:])

        fn_sql = "results_pearson_all.sqlite"
        db_out = to_test_results(fn_sql)
        group_features(self.df, db_out, max_rt_diff=200.0, coeff_thres=0.0, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)

        records = sqlite_records(to_test_results(fn_sql), "groups")
        records_comp = sqlite_records(to_test_data(fn_sql), "groups")
        for i in range(len(records)):
            self.assertEqual(records[i][0:6], records_comp[i][0:6])
            np.testing.assert_almost_equal(records[i][6:], records_comp[i][6:])

        fn_sql = "results_pearson_all.sqlite"
        db_out = to_test_results(fn_sql)
        group_features(self.df, db_out, max_rt_diff=200.0, coeff_thres=0.0, pvalue_thres=1.0, method="pearson", block=20, ncpus=1)

        records = sqlite_records(to_test_results(fn_sql), "groups")
        records_comp = sqlite_records(to_test_data(fn_sql), "groups")
        for i in range(len(records)):
            self.assertEqual(records[i][0:6], records_comp[i][0:6])
            np.testing.assert_almost_equal(records[i][6:], records_comp[i][6:])

        fn_sql = "results_spearman.sqlite"
        db_out = to_test_results(fn_sql)
        group_features(self.df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="spearman", block=5000, ncpus=None)
        records = sqlite_records(to_test_results(fn_sql), "groups")
        records_comp = sqlite_records(to_test_data(fn_sql), "groups")

        for i in range(len(records)):
            self.assertEqual(records[i][0:6], records_comp[i][0:6])
            np.testing.assert_almost_equal(records[i][6:], records_comp[i][6:])


if __name__ == '__main__':
    unittest.main()