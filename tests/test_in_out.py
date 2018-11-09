#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
from collections import OrderedDict
from beams.in_out import *
from tests.utils import to_test_data
import numpy as np

class InOutTestCase(unittest.TestCase):

    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_read_peaklist(self):

        self.df_peaklist = read_peaklist(to_test_data("peaklist_lcms_pos_theoretical.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "M127T60")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "M550T200")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.9792044)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.0658904)

        self.assertEqual(self.df_peaklist["rt"].iloc[0], 60)
        self.assertEqual(self.df_peaklist["rt"].iloc[-1], 200)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

        self.df_peaklist = read_peaklist(to_test_data("peaklist_dims_pos_theoretical.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "126_9792044")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "550_0658904")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.9792044)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.0658904)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

    def test_combine_peaklist_matrix(self):
        df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_theoretical.txt"))

        self.assertEqual(df["name"].iloc[0], "M127T60")
        self.assertEqual(df["name"].iloc[-1], "M550T200")

        self.assertEqual(df["mz"].iloc[0], 126.9792044)
        self.assertEqual(df["mz"].iloc[-1], np.float64(550.0658904))

        self.assertEqual(df["rt"].iloc[0], 60)
        self.assertEqual(df["rt"].iloc[-1], 200)

        self.assertEqual(df["intensity"].iloc[0], 1421.775)
        self.assertEqual(df["intensity"].iloc[-1], 4549.65)

    def test_read_molecular_formulae(self):

        db_molecular_formula = os.path.join(self.path, "beams", "data", "db_mf.txt")
        records = read_molecular_formulae(db_molecular_formula, separator="\t")
        self.assertEqual(len(records), 13061)

        record_01 = [("composition", OrderedDict([('C', 1), ('H', 2), ('O', 1)])), ('CHNOPS', True),
                     ('exact_mass', 30.010565),
                     ('HC', 1), ('NOPSC', 1), ('lewis', 1), ('senior', 1), ('double_bond_equivalents', 1.0)]
        record_02 = [("composition", OrderedDict([('C', 17), ('H', 19), ('Cl', 1), ('N', 2), ('O', 1), ('S', 1)])), ('CHNOPS', False),
                     ('exact_mass', 334.090662),
                     ('HC', 1), ('NOPSC', 1), ('lewis', 0), ('senior', 1), ('double_bond_equivalents', 9.0)]
        record_03 = [("composition", OrderedDict([('C', 48), ('H', 86), ('O', 18), ('P', 2)])), ('CHNOPS', True),
                     ('exact_mass', 1012.528942),
                     ('HC', 1), ('NOPSC', 1), ('lewis', 1), ('senior', 1), ('double_bond_equivalents', 6.0)]

        self.assertEqual(records[0], OrderedDict(record_01))
        self.assertEqual(records[5000], OrderedDict(record_02))
        self.assertEqual(records[-1], OrderedDict(record_03))

    def test_read_compounds(self):

        db_compounds = os.path.join(self.path, "beams", "data", "db_compounds.txt")
        records = read_compounds(db_compounds, separator="\t")
        self.assertEqual(len(records), 31644)
        record_01 = [("composition", OrderedDict([('C', 10), ('Cl', 10), ('O', 1)])), ('CHNOPS', False),
                     ('exact_mass', 485.683442),
                     ('compound_id', 1638),
                     ('compound_name', 'Chlordecone'),
                     ('molecular_formula', 'C10Cl10O')]
        record_02 = [("composition", OrderedDict([('C', 24), ('H', 42), ('O', 21)])), ('CHNOPS', True),
                     ('exact_mass', 666.221865),
                     ('compound_id', 17543),
                     ('compound_name', '6G,6-kestotetraose'),
                     ('molecular_formula', 'C24H42O21')]
        record_03 = [("composition", OrderedDict([('H', 1), ('N', 1), ('O', 3)])), ('CHNOPS', True),
                     ('exact_mass', 62.995644),
                     ('compound_id', 40762),
                     ('compound_name', 'Peroxynitrite'),
                     ('molecular_formula', 'HNO3')]

        self.assertEqual(records[0], OrderedDict(record_01))
        self.assertEqual(records[14000], OrderedDict(record_02))
        self.assertEqual(records[-1], OrderedDict(record_03))


    def test_read_adducts(self):
        adducts_lib = os.path.join(self.path, "beams", "data", "adducts.txt")
        records_pos = read_adducts(adducts_lib, "pos")
        records_pos_comp = [('[M+H]+', 1.0072764), ('[M+Na]+', 22.9892214),
                            ('[M+K]+', 38.9631594)]
        self.assertEqual(records_pos.lib, OrderedDict(records_pos_comp))
        records_neg = read_adducts(adducts_lib, "neg")
        records_neg_comp = [('[M-H]-', -1.0072764), ('[M+Na-2H]-', 20.9746686),
                            ('[M+Cl]-', 34.9694016), ('[M+K-2H]-', 36.9486066),
                            ('[M+Hac-H]-', 59.0138536)]
        self.assertEqual(records_neg.lib, OrderedDict(records_neg_comp))

    def test_read_isotopes(self):
        isotopes_lib = os.path.join(self.path, "beams", "data", "isotopes.txt")
        records_pos = read_isotopes(isotopes_lib, "pos")
        records_pos_comp = [OrderedDict([('C', {'abundance': 100.0}), ('(13C)', {'abundance': 1.1}), ('mass_difference', 1.003355)]),
                            OrderedDict([('S', {'abundance': 100.0}), ('(34S)', {'abundance': 4.21}), ('mass_difference', 1.995796)]),
                            OrderedDict([('K', {'abundance': 100.0}), ('(41K)', {'abundance': 6.73}), ('mass_difference', 1.998117)])]
        self.assertEqual(records_pos.lib, records_pos_comp)
        records_neg = read_isotopes(isotopes_lib, "neg")
        records_neg_comp = [OrderedDict([('C', {'abundance': 100.0}), ('(13C)', {'abundance': 1.1}), ('mass_difference', 1.003355)]),
                            OrderedDict([('S', {'abundance': 100.0}), ('(34S)', {'abundance': 4.21}), ('mass_difference', 1.995796)]),
                            OrderedDict([('Cl', {'abundance': 100.0}), ('(37Cl)', {'abundance': 24.23}), ('mass_difference', 1.99705)])]
        self.assertEqual(records_neg.lib, records_neg_comp)


if __name__ == '__main__':
    unittest.main()
