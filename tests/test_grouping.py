#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import os
import unittest
import pandas as pd
from utils import to_test_data, to_test_results, sqlite_records, create_matrix
from beams.in_out import read_peaklist
from beams.grouping import group_features


class GroupFeaturesTest(unittest.TestCase):

    def setUp(self):
        self.df_peaklist = read_peaklist(to_test_data("peaklist_dims_pos_theoretical.txt"))
        self.df_matrix = create_matrix(self.df_peaklist, 50)
        self.df_matrix = self.df_matrix.set_index("name")
        self.df_matrix.ix["492_0604104"] = self.df_matrix.ix["493_0637654"] * 1.1
        self.df_matrix.ix["230_9901644"] = self.df_matrix.ix["215_0162264"] * 2.0
        self.df_matrix = self.df_matrix.reset_index(drop=False)
        self.df_peaklist = self.df_peaklist.assign(intensity=pd.Series(self.df_matrix.median(axis=1, skipna=True).values))
        self.df = pd.merge(self.df_peaklist, self.df_matrix, how='left', left_on="name", right_on="name")

    def test_group_features(self):
        db_out = to_test_results("db_out_pearson.sqlite")
        group_features(self.df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=None, method="pearson", block=5000, ncpus=None)
        self.assertEqual(sqlite_records(to_test_data("db_out_pearson.sqlite"), "groups"), sqlite_records(db_out, "groups"))

        db_out = to_test_results("db_out_spearman.sqlite")
        group_features(self.df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=None, method="spearman", block=5000, ncpus=None)
        self.assertEqual(sqlite_records(to_test_data("db_out_spearman.sqlite"), "groups"), sqlite_records(db_out, "groups"))
