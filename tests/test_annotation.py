#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import unittest
import gzip
import shutil
import pandas as pd
import filecmp

from beamspy.annotation import *
from beamspy.grouping import group_features
from beamspy.in_out import *
from tests.utils import *


class AnnotationTestCase(unittest.TestCase):

    def setUp(self):

        self.df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical.txt"), to_test_data("dataMatrix_lcms_theoretical.txt"))
        self.path, f = os.path.split(os.path.dirname(os.path.abspath(__file__)))

        self.lib_isotopes = read_isotopes(os.path.join(self.path, "beamspy", "data", "isotopes.txt"), "pos")
        self.lib_adducts = read_adducts(os.path.join(self.path, "beamspy", "data", "adducts.txt"), "pos")
        self.lib_multiple_charged_ions = read_multiple_charged_ions(os.path.join(self.path, "beamspy", "data", "multiple_charged_ions.txt"), "pos")
        # lib_mass_differences = read_mass_differences(os.path.join(self.path, "beamspy", "data", "multiple_charged_differences.txt"), "pos")

        self.db_results = "results_annotation.sqlite"
        self.db_results_graph = "results_annotation_graph.sqlite"
        self.graph = group_features(self.df, to_test_results(self.db_results_graph), max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)

        self.ppm = 2.0

    #def tearDown(self):
    #    os.remove(to_test_results("hmdb_full_v4_0_v1.sqlite"))

    def test_annotate_adducts(self):
        annotate_adducts(self.df, to_test_results(self.db_results), self.ppm, self.lib_adducts)
        self.assertEqual(sqlite_records(to_test_results(self.db_results), "adduct_pairs"), sqlite_records(to_test_data(self.db_results), "adduct_pairs"))

        annotate_adducts(self.graph, to_test_results(self.db_results_graph), self.ppm, self.lib_adducts)
        self.assertEqual(sqlite_records(to_test_results(self.db_results_graph), "adduct_pairs"), sqlite_records(to_test_data(self.db_results_graph), "adduct_pairs"))

    def test_annotate_isotopes(self):
        annotate_isotopes(self.df, to_test_results(self.db_results), self.ppm, self.lib_isotopes)
        self.assertEqual(sqlite_records(to_test_results(self.db_results), "isotopes"), sqlite_records(to_test_data(self.db_results), "isotopes"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "isotopes"), 1)

        annotate_isotopes(self.graph, to_test_results(self.db_results_graph), self.ppm, self.lib_isotopes)
        self.assertEqual(sqlite_records(to_test_results(self.db_results_graph), "isotopes"), sqlite_records(to_test_data(self.db_results_graph), "isotopes"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results_graph), "isotopes"), 1)

    def test_annotate_oligomers(self):
        annotate_oligomers(self.df, to_test_results(self.db_results), self.ppm, self.lib_adducts, maximum=5)
        self.assertEqual(sqlite_records(to_test_results(self.db_results), "oligomers"),
                         sqlite_records(to_test_data(self.db_results), "oligomers"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "oligomers"), 2)

        annotate_oligomers(self.graph, to_test_results(self.db_results_graph), self.ppm, self.lib_adducts)
        self.assertEqual(sqlite_records(to_test_results(self.db_results_graph), "oligomers"),
                         sqlite_records(to_test_data(self.db_results_graph), "oligomers"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results_graph), "oligomers"), 2)

    # def test_annotate_drug_products(self):
    #     df = pd.DataFrame({"name": pd.Series(["M152T100", "M188T100", "M310T200", "M348T200"]),
    #                        "mz": pd.Series([152.0706054, 188.0682004, 310.1413254, 348.0972084], dtype=np.float64),
    #                        "rt": pd.Series([100.0, 100, 0, 200.0, 200.0], dtype=np.float64),
    #                        "intensity": pd.Series([1234.45, 2345.67, 3456.78, 4567.89], dtype=np.float64)},
    #                        columns=["name", "mz", "rt", "intensity"],
    #                        index=range(0, 4))
    #     smiles = ["CC(=O)NC1=CC=C(C=C1)O", "CNCCC(OC1=CC=C(C=C1)C(F)(F)F)C1=CC=CC=C1"]
    #     annotate_drug_products(df, to_test_results(self.db_results), smiles, self.lib_adducts, self.ppm,
    #                            phase1_cycles=1, phase2_cycles=1)
    #     self.assertEqual(sqlite_records(to_test_results(self.db_results), "drug_products"),
    #                      sqlite_records(to_test_data(self.db_results), "drug_products"))
    #     self.assertEqual(sqlite_count(to_test_results(self.db_results), "drug_products"), 4)

    def test_annotate_compounds(self):

        db_name = "hmdb_full_v4_0_v1"

        path_hmdb_sql_gz = os.path.join(os.getcwd(), "beamspy", "data", "databases", db_name + ".sql.gz")
        path_hmdb_sqlite = to_test_results("{}.sqlite".format(db_name))

        if os.path.isfile(path_hmdb_sqlite):
            os.remove(path_hmdb_sqlite)

        with gzip.GzipFile(path_hmdb_sql_gz, mode='rb') as db_dump:
            conn = sqlite3.connect(path_hmdb_sqlite)
            cursor = conn.cursor()
            cursor.executescript(db_dump.read().decode('utf-8'))
            conn.commit()
            conn.close()

        # sqlite profile provided
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), db_name, path_hmdb_sqlite)
        self.assertEqual(sqlite_records(to_test_results(self.db_results), "compounds_{}".format(db_name)),
                         sqlite_records(to_test_data(self.db_results), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 57)

        # internal sqlite databases
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), db_name)
        self.assertEqual(sqlite_records(to_test_results(self.db_results), "compounds_{}".format(db_name)),
                         sqlite_records(to_test_data(self.db_results), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 57)

        path_db_txt = os.path.join(os.getcwd(), "beamspy", "data", "db_compounds.txt")
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), "test", path_db_txt)
        #self.assertEqual(sqlite_records(to_test_results(self.db_results), "compounds_{}".format(db_name)), sqlite_records(to_test_data(self.db_results), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 57)

    def test_annotate_molecular_formulae(self):
        fn_mf = os.path.join(self.path, "beamspy", "data", "db_mf.txt")
        annotate_molecular_formulae(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), fn_mf)
        self.assertEqual(sqlite_records(to_test_results(self.db_results), "molecular_formulae"), sqlite_records(to_test_data(self.db_results), "molecular_formulae"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "molecular_formulae"), 16)

        db_mfdb_results = "results_mfdb.sqlite"
        annotate_molecular_formulae(self.df, self.lib_adducts, self.ppm, to_test_results(db_mfdb_results))
        self.assertEqual(sqlite_records(to_test_results(db_mfdb_results), "molecular_formulae"), sqlite_records(to_test_data(db_mfdb_results), "molecular_formulae"))
        self.assertEqual(sqlite_count(to_test_results(db_mfdb_results), "molecular_formulae"), 586)

    def test_summary(self):

        def _assert(summary_test_data, summary_result):
            with open(summary_result) as result:
                with open(summary_test_data) as test_data:
                    lines_results = result.read().splitlines()
                    lines_test_data = test_data.read().splitlines()
                    for i in range(len(lines_results)):
                        self.assertEqual(lines_results[i], lines_test_data[i])
                        self.assertEqual(sqlite_records(to_test_results(self.db_results), "summary"), sqlite_records(to_test_data(self.db_results), "summary"))

        fn_summary = "summary_mr_mc.txt"
        df_summary = summary(self.df, to_test_results(self.db_results), single_row=False, single_column=False, convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_sr_mc.txt"
        df_summary = summary(self.df, to_test_results(self.db_results), single_row=True, single_column=False, convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_sr_sc.txt"
        df_summary = summary(self.df, to_test_results(self.db_results), single_row=True, single_column=True, convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_mr_mc_graphs.txt"
        df_summary = summary(self.df, to_test_results(self.db_results_graph), single_row=False, single_column=False, convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

if __name__ == '__main__':
    unittest.main()