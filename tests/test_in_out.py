#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
from collections import OrderedDict
from beamspy.in_out import *
from tests.utils import to_test_data
import numpy as np


class InOutTestCase(unittest.TestCase):

    def setUp(self):
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

    def test_read_peaklist(self):

        self.df_peaklist = read_peaklist(to_test_data("peaklist_lcms_pos_theoretical.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "M127T60")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "M550T200")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.979204)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.065890)

        self.assertEqual(self.df_peaklist["rt"].iloc[0], 60)
        self.assertEqual(self.df_peaklist["rt"].iloc[-1], 200)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

        self.df_peaklist = read_peaklist(to_test_data("peaklist_lcms_pos_theoretical_no_name.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "M127T60")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "M550T200")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.979204)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.065890)

        self.assertEqual(self.df_peaklist["rt"].iloc[0], 60)
        self.assertEqual(self.df_peaklist["rt"].iloc[-1], 200)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

        self.df_peaklist = read_peaklist(to_test_data("peaklist_dims_pos_theoretical.txt"))

        self.assertEqual(self.df_peaklist["name"].iloc[0], "126_979204")
        self.assertEqual(self.df_peaklist["name"].iloc[-1], "550_065890")

        self.assertEqual(self.df_peaklist["mz"].iloc[0], 126.979204)
        self.assertEqual(self.df_peaklist["mz"].iloc[-1], 550.065890)

        self.assertEqual(self.df_peaklist["intensity"].iloc[0], 1421.78)
        self.assertEqual(self.df_peaklist["intensity"].iloc[-1], 4549.65)

    def test_combine_peaklist_matrix(self):
        df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_lcms_theoretical.txt"))

        self.assertEqual(df["name"].iloc[0], "M127T60")
        self.assertEqual(df["name"].iloc[-1], "M550T200")

        self.assertEqual(df["mz"].iloc[0], 126.979204)
        self.assertEqual(df["mz"].iloc[-1], 550.065890)

        self.assertEqual(df["rt"].iloc[0], 60)
        self.assertEqual(df["rt"].iloc[-1], 200)

        self.assertEqual(df["intensity"].iloc[0], 1421.775)
        self.assertEqual(df["intensity"].iloc[-1], 4549.65)

        df = combine_peaklist_matrix(to_test_data("peaklist_dims_pos_theoretical.txt"), to_test_data("dataMatrix_dims_theoretical.txt"))

        self.assertEqual(df["name"].iloc[0], "126_979204")
        self.assertEqual(df["name"].iloc[-1], "550_065890")
        
        self.assertEqual(df["mz"].iloc[0], 126.979204)
        self.assertEqual(df["mz"].iloc[-1], 550.065890)

        self.assertEqual(df["intensity"].iloc[0], 1421.78)
        self.assertEqual(df["intensity"].iloc[-1], 4549.65)
    
    def test_read_xset_matrix(self):
        df = read_xset_matrix(to_test_data("xset_matrix.txt"), "sample01")

        self.assertEqual(df["name"].iloc[0], "M127T60")
        self.assertEqual(df["name"].iloc[-1], "M550T200")

        np.testing.assert_almost_equal(df["mz"].iloc[0], 126.979204, 8)
        np.testing.assert_almost_equal(df["mz"].iloc[-1], 550.065890, 8)

        self.assertEqual(df["rt"].iloc[0], 60)
        self.assertEqual(df["rt"].iloc[-1], 200)

        np.testing.assert_almost_equal(df["intensity"].iloc[0], 1421.775, 8)
        np.testing.assert_almost_equal(df["intensity"].iloc[-1], 4549.65, 8)

    def test_read_molecular_formulae(self):

        db_molecular_formula = os.path.join(self.path, "beamspy", "data", "db_mf.txt")
        records = read_molecular_formulae(db_molecular_formula, separator="\t")
        self.assertEqual(len(records), 13061)

        record_01 = [("composition", OrderedDict([('C', 1), ('H', 2), ('O', 1)])), ('CHNOPS', True),
                     ('exact_mass', 30.010565),
                     ('HC', 1), ('NOPSC', 1), ('lewis', 1), ('senior', 1), ('double_bond_equivalents', 1.0)]
        record_02 = [("composition", OrderedDict([('C', 17), ('H', 19), ('Cl', 1), ('N', 2), ('O', 1), ('S', 1)])), ('CHNOPS', False),
                     ('exact_mass', 334.090662),
                     ('HC', 1), ('NOPSC', 1), ('lewis', 0), ('senior', 1), ('double_bond_equivalents', 9.0)]
        record_03 = [("composition", OrderedDict([('C', 48), ('H', 86), ('O', 18), ('P', 2)])), ('CHNOPS', True),
                     ('exact_mass', 1012.528940),
                     ('HC', 1), ('NOPSC', 1), ('lewis', 1), ('senior', 1), ('double_bond_equivalents', 6.0)]

        self.assertEqual(records[0], OrderedDict(record_01))
        self.assertEqual(records[5000], OrderedDict(record_02))
        self.assertEqual(records[-1], OrderedDict(record_03))

    def test_read_compounds(self):

        db_compounds = os.path.join(self.path, "beamspy", "data", "db_compounds.txt")
        records = read_compounds(db_compounds, separator="\t")
        self.assertEqual(len(records), 31644)
        record_01 = [("composition", OrderedDict([('C', 10), ('Cl', 10), ('O', 1)])), ('CHNOPS', False),
                     ('exact_mass', 485.683441),
                     ('compound_id', 1638),
                     ('compound_name', 'Chlordecone'),
                     ('molecular_formula', 'C10Cl10O')]
        record_02 = [("composition", OrderedDict([('C', 24), ('H', 42), ('O', 21)])), ('CHNOPS', True),
                     ('exact_mass', 666.221858),
                     ('compound_id', 17543),
                     ('compound_name', '6G,6-kestotetraose'),
                     ('molecular_formula', 'C24H42O21')]
        record_03 = [("composition", OrderedDict([('H', 1), ('N', 1), ('O', 3)])), ('CHNOPS', True),
                     ('exact_mass', 62.995643),
                     ('compound_id', 40762),
                     ('compound_name', 'Peroxynitrite'),
                     ('molecular_formula', 'HNO3')]

        self.assertEqual(records[0], OrderedDict(record_01))
        self.assertEqual(records[14000], OrderedDict(record_02))
        self.assertEqual(records[-1], OrderedDict(record_03))


    def test_read_adducts(self):
        adducts_lib = os.path.join(self.path, "beamspy", "data", "adducts.txt")
        records_pos = read_adducts(adducts_lib, "pos")
        records_pos_comp = OrderedDict([('[M+H]+', OrderedDict([('mass', 1.007276), ('charge', 1)])),
                                        ('[M+Na]+', OrderedDict([('mass', 22.989221), ('charge', 1)])),
                                        ('[M+K]+', OrderedDict([('mass', 38.963158), ('charge', 1)]))])
        self.assertEqual(records_pos.lib, OrderedDict(records_pos_comp))
        records_neg = read_adducts(adducts_lib, "neg")
        records_neg_comp = OrderedDict([('[M-H]-', OrderedDict([('mass', -1.007276), ('charge', 1)])),
                                        ('[M+Na-2H]-', OrderedDict([('mass', 20.974668), ('charge', 1)])),
                                        ('[M+Cl]-', OrderedDict([('mass', 34.969401), ('charge', 1)])),
                                        ('[M+K-2H]-', OrderedDict([('mass', 36.948605), ('charge', 1)])),
                                        ('[M+Hac-H]-', OrderedDict([('mass', 59.013853), ('charge', 1)]))])
        self.assertEqual(records_neg.lib, OrderedDict(records_neg_comp))

    def test_read_isotopes(self):
        isotopes_lib = os.path.join(self.path, "beamspy", "data", "isotopes.txt")
        records_pos = read_isotopes(isotopes_lib, "pos")
        records_pos_comp = [OrderedDict([('C', {'abundance': 98.93}), ('(13C)', {'abundance': 1.07}), ('mass_difference', 0.5016775), ('charge', 2)]),
                            OrderedDict([('C', {'abundance': 98.93}), ('(13C)', {'abundance': 1.07}), ('mass_difference', 1.003355), ('charge', 1)]),
                            OrderedDict([('S', {'abundance': 94.99}), ('(34S)', {'abundance': 4.25}), ('mass_difference', 1.995796), ('charge', 1)]),
                            OrderedDict([('K', {'abundance': 93.25}), ('(41K)', {'abundance': 6.73}), ('mass_difference', 1.998119), ('charge', 1)])]
        self.assertEqual(records_pos.lib, records_pos_comp)
        records_neg = read_isotopes(isotopes_lib, "neg")
        records_neg_comp = [OrderedDict([('C', {'abundance': 98.93}), ('(13C)', {'abundance': 1.07}), ('mass_difference', 0.5016775), ('charge', 2)]),
                            OrderedDict([('C', {'abundance': 98.93}), ('(13C)', {'abundance': 1.07}), ('mass_difference', 1.003355), ('charge', 1)]),
                            OrderedDict([('S', {'abundance': 94.99}), ('(34S)', {'abundance': 4.25}), ('mass_difference', 1.995796), ('charge', 1)]),
                            OrderedDict([('Cl', {'abundance': 75.76}), ('(37Cl)', {'abundance': 24.24}), ('mass_difference', 1.997050), ('charge', 1)])]
        self.assertEqual(records_neg.lib, records_neg_comp)

    # def test_read_mass_differences(self):
    #     differences_lib = os.path.join(self.path, "beamspy", "data", "adducts_differences.txt")
    #     records = read_mass_differences(differences_lib, ion_mode="pos")
    #     self.assertEqual(records.lib, [OrderedDict([('[M+H]+', {'charge': 1.0}),
    #                                                 ('[M+Na]+', {'charge': 1.0}),
    #                                                 ('mass_difference', 21.981945)])])
    #     records = read_mass_differences(differences_lib, ion_mode="neg")
    #     self.assertEqual(records.lib, [])
    #     records = read_mass_differences(differences_lib, ion_mode="both")
    #     self.assertEqual(records.lib, [])


if __name__ == '__main__':
    unittest.main()
