#!/usr/bin/env python
#  -*- coding: utf-8 -*-
import unittest
from utils import to_test_data
from beams.in_out import read_peaklist, combine_peaklist_matrix


class InOutTestCase(unittest.TestCase):

    def test_read_peaklist(self):

        self.df_peaklist = read_peaklist(to_test_data("peaklist_lcms_pos_theoretical.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "M127T60")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "M550T200")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.9792044)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.0658904000001)

        self.assertEqual(self.df_peaklist["rt"].iloc[0], 60)
        self.assertEqual(self.df_peaklist["rt"].iloc[-1], 200)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

        self.df_peaklist = read_peaklist(to_test_data("peaklist_dims_pos_theoretical.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "126_9792044")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "550_0658904")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.9792044)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.0658904000001)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

    def test_combine_peaklist_matrix(self):
        df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_theoretical.txt"))

        self.assertEqual(df["name"].iloc[0], "M127T60")
        self.assertEqual(df["name"].iloc[-1], "M550T200")

        self.assertEqual(df["mz"].iloc[0], 126.9792044)
        self.assertEqual(df["mz"].iloc[-1], 550.0658904000001)

        self.assertEqual(df["rt"].iloc[0], 60)
        self.assertEqual(df["rt"].iloc[-1], 200)

        self.assertEqual(df["intensity"].iloc[0], 1421.775)
        self.assertEqual(df["intensity"].iloc[-1], 4549.65)

if __name__ == '__main__':
    unittest.main()