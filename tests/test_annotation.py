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
        # self.lib_multiple_charged_ions = read_multiple_charged_ions(os.path.join(self.path, "beamspy", "data", "multiple_charged_ions.txt"), "pos")
        # lib_mass_differences = read_mass_differences(os.path.join(self.path, "beamspy", "data", "multiple_charged_differences.txt"), "pos")

        self.db_results = "results_annotation.sqlite"
        self.db_results_graph = "results_annotation_graph.sqlite"
        self.graph = group_features(self.df, to_test_results(self.db_results_graph), max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)

        self.ppm = 2.0

    def test_neutral_losses(self):

        df_nls = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical_nls.txt"),
                                     to_test_data("dataMatrix_lcms_theoretical_nls.txt"))

        db_nls = "results_annotation_nls.sqlite"

        graph = group_features(df_nls, to_test_results(db_nls), max_rt_diff=5.0, coeff_thres=0.7,
                               pvalue_thres=1.0, method="pearson", block=5000, ncpus=None)

        lib_nls = read_neutral_losses(os.path.join(self.path, "beamspy", "data", "neutral_losses.txt"))

        annotate_neutral_losses(graph, to_test_results(db_nls), self.ppm, lib_nls)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_nls), "neutral_losses"),
                                 sqlite_records(to_test_data(db_nls), "neutral_losses"))

        annotate_neutral_losses(df_nls, to_test_results(db_nls), self.ppm, lib_nls)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_nls), "neutral_losses"),
                                 sqlite_records(to_test_data(db_nls), "neutral_losses"))

        lib_adducts = read_adducts(os.path.join(self.path, "beamspy", "data", "multiple_charged_ions.txt"), "pos")

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

        annotate_compounds(df_nls, lib_adducts, self.ppm, to_test_results(db_nls), db_name,
                           filter=True, db_in=path_hmdb_sqlite)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_nls), "compounds_{}".format(db_name)),
                                 sqlite_records(to_test_data(db_nls), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(db_nls), "compounds_{}".format(db_name)), 12)

        annotate_molecular_formulae(df_nls, self.lib_adducts, self.ppm, to_test_results(db_nls),
                                    filter=True, rules=True)

        self.assertSequenceEqual(sqlite_records(to_test_results(db_nls), "molecular_formulae"),
                                 sqlite_records(to_test_data(db_nls), "molecular_formulae"))
        self.assertEqual(sqlite_count(to_test_results(db_nls), "molecular_formulae"), 6)

    #def tearDown(self):
    #    os.remove(to_test_results("hmdb_full_v4_0_v1.sqlite"))

    def test_annotate_multiple_charged_adducts(self):
        df = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical_mc_o.txt"),
                                          to_test_data("dataMatrix_lcms_theoretical_mc_o.txt"))

        lib_adducts = read_adducts(os.path.join(self.path, "beamspy", "data", "multiple_charged_ions.txt"), "pos")
        db_mc = "results_annotation_mc_o.sqlite"

        annotate_adducts(df, to_test_results(db_mc), self.ppm, lib_adducts)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_mc), "adduct_pairs"),
                                 sqlite_records(to_test_data(db_mc), "adduct_pairs"))

        annotate_isotopes(df, to_test_results(db_mc), self.ppm, self.lib_isotopes)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_mc), "isotopes"),
                                 sqlite_records(to_test_data(db_mc), "isotopes"))
        self.assertEqual(sqlite_count(to_test_results(db_mc), "isotopes"), 2)

        annotate_oligomers(df, to_test_results(db_mc), self.ppm, lib_adducts)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_mc), "oligomers"),
                                 sqlite_records(to_test_data(db_mc), "oligomers"))
        self.assertEqual(sqlite_count(to_test_results(db_mc), "oligomers"), 4)

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

        annotate_compounds(df, lib_adducts, self.ppm, to_test_results(db_mc), db_name,
                           filter=True, db_in=path_hmdb_sqlite)

        self.assertSequenceEqual(sqlite_records(to_test_results(db_mc), "compounds_{}".format(db_name)),
                                 sqlite_records(to_test_data(db_mc), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(db_mc), "compounds_{}".format(db_name)), 41)

        annotate_molecular_formulae(df, lib_adducts, self.ppm, to_test_results(db_mc),
                                    filter=True, rules=True)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_mc), "molecular_formulae"),
                                 sqlite_records(to_test_data(db_mc), "molecular_formulae"))
        self.assertEqual(sqlite_count(to_test_results(db_mc), "molecular_formulae"), 4257)

    def test_annotate_adducts(self):
        annotate_adducts(self.df, to_test_results(self.db_results), self.ppm, self.lib_adducts)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "adduct_pairs"),
                                 sqlite_records(to_test_data(self.db_results), "adduct_pairs"))

        annotate_adducts(self.graph, to_test_results(self.db_results_graph), self.ppm, self.lib_adducts)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results_graph), "adduct_pairs"),
                                 sqlite_records(to_test_data(self.db_results_graph), "adduct_pairs"))

    def test_annotate_isotopes(self):
        annotate_isotopes(self.df, to_test_results(self.db_results), self.ppm, self.lib_isotopes)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "isotopes"),
                                 sqlite_records(to_test_data(self.db_results), "isotopes"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "isotopes"), 1)

        annotate_isotopes(self.graph, to_test_results(self.db_results_graph), self.ppm, self.lib_isotopes)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results_graph), "isotopes"),
                                 sqlite_records(to_test_data(self.db_results_graph), "isotopes"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results_graph), "isotopes"), 1)

    def test_annotate_oligomers(self):
        annotate_oligomers(self.df, to_test_results(self.db_results), self.ppm, self.lib_adducts)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "oligomers"),
                                 sqlite_records(to_test_data(self.db_results), "oligomers"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "oligomers"), 2)

        annotate_oligomers(self.graph, to_test_results(self.db_results_graph), self.ppm, self.lib_adducts)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results_graph), "oligomers"),
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

        annotate_adducts(self.df, to_test_results(self.db_results), self.ppm, self.lib_adducts)
        annotate_isotopes(self.df, to_test_results(self.db_results), self.ppm, self.lib_isotopes)
        annotate_oligomers(self.df, to_test_results(self.db_results), self.ppm, self.lib_adducts)

        annotate_adducts(self.df, to_test_results(self.db_results_graph), self.ppm, self.lib_adducts)
        annotate_isotopes(self.df, to_test_results(self.db_results_graph), self.ppm, self.lib_isotopes)
        annotate_oligomers(self.df, to_test_results(self.db_results_graph), self.ppm, self.lib_adducts)

        if os.path.isfile(path_hmdb_sqlite):
            os.remove(path_hmdb_sqlite)

        with gzip.GzipFile(path_hmdb_sql_gz, mode='rb') as db_dump:
            conn = sqlite3.connect(path_hmdb_sqlite)
            cursor = conn.cursor()
            cursor.executescript(db_dump.read().decode('utf-8'))
            conn.commit()
            conn.close()

        # sqlite file provided
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), db_name,
                           filter=True, db_in=path_hmdb_sqlite)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "compounds_{}".format(db_name)),
                                 sqlite_records(to_test_data(self.db_results), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 51)

        # internal sqlite databases
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), db_name,
                           filter=True, db_in="")
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "compounds_{}".format(db_name)),
                                 sqlite_records(to_test_data(self.db_results), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 51)

        # internal sqlite databases, including grouping
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results_graph), db_name,
                           filter=True, db_in="")
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results_graph), "compounds_{}".format(db_name)),
                                 sqlite_records(to_test_data(self.db_results_graph), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results_graph), "compounds_{}".format(db_name)), 51)

        # internal sqlite databases (excl. patterns)
        db_results_excl_patterns = self.db_results.replace(".sqlite", "_excl_pattern.sqlite")
        annotate_adducts(self.df, to_test_results(db_results_excl_patterns), self.ppm, self.lib_adducts)
        annotate_isotopes(self.df, to_test_results(db_results_excl_patterns), self.ppm, self.lib_isotopes)

        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(db_results_excl_patterns), db_name,
                           filter=False, db_in="")
        self.assertEqual(sqlite_count(to_test_results(db_results_excl_patterns), "compounds_{}".format(db_name)), 58)

        # text file provided
        path_db_txt = os.path.join(os.getcwd(), "beamspy", "data", "db_compounds.txt")
        db_name = "test"
        annotate_compounds(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results), db_name,
                           filter=True, db_in=path_db_txt)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "compounds_{}".format(db_name)),
                                sqlite_records(to_test_data(self.db_results), "compounds_{}".format(db_name)))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 66)

        path_db_txt = to_test_results("db_compounds_rt.txt")
        db_name = "test_rt"
        with open(path_db_txt, "w") as out:
            out.write("compound_id\tmolecular_formula\tcompound_name\tadduct\tretention_time\n")
            out.write("HMDB0000263\tC3H5O6P\tPhosphoenolpyruvic acid\t[M+H]+\t118.0\n")

        annotate_compounds(self.df, self.lib_adducts, 100.0, to_test_results(self.db_results), db_name, filter=True,
                           db_in=path_db_txt, rt_tol=5.0)
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "compounds_{}".format(db_name)), 1)

    def test_annotate_molecular_formulae(self):

        fn_mf = os.path.join(self.path, "beamspy", "data", "db_mf.txt")
        annotate_molecular_formulae(self.df, self.lib_adducts, self.ppm, to_test_results(self.db_results),
                                    fn_mf, filter=True)

        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "molecular_formulae"),
                                 sqlite_records(to_test_data(self.db_results), "molecular_formulae"))
        self.assertEqual(sqlite_count(to_test_results(self.db_results), "molecular_formulae"), 15)

        db_mfdb_results = "results_mfdb.sqlite"
        annotate_molecular_formulae(self.df, self.lib_adducts, self.ppm, to_test_results(db_mfdb_results),
                                    filter=False, rules=True)
        self.assertSequenceEqual(sqlite_records(to_test_results(db_mfdb_results), "molecular_formulae"),
                                 sqlite_records(to_test_data(db_mfdb_results), "molecular_formulae"))
        self.assertEqual(sqlite_count(to_test_results(db_mfdb_results), "molecular_formulae"), 586)

        db_mfdb_results = "results_mfdb_excl_hrules.sqlite"
        annotate_adducts(self.df, to_test_results(db_mfdb_results), self.ppm, self.lib_adducts)
        annotate_isotopes(self.df, to_test_results(db_mfdb_results), self.ppm, self.lib_isotopes)

        annotate_molecular_formulae(self.df, self.lib_adducts, self.ppm, to_test_results(db_mfdb_results),
                                    filter=True, rules=False)
        records = sqlite_records(to_test_results(db_mfdb_results), "molecular_formulae")
        self.assertEqual(len(records), 3869)
        self.assertSequenceEqual(records[897],
                                 ('M493T192', 493.063765, 493.06376, 0.010140676303965665, '[M+Na]+', '(13C)', '',
                                  14, 23, 4, 8, 2, 1, 1, 'C14H23N4O8P2S', 1, 1, 0, 1, 5.0))

    def test_summary(self):

        def _assert(summary_test_data, summary_result):
            with open(summary_result) as result:
                with open(summary_test_data) as test_data:
                    lines_results = result.read().splitlines()
                    lines_test_data = test_data.read().splitlines()
                    for i in range(len(lines_results)):
                        self.assertEqual(lines_results[i], lines_test_data[i])

        fn_summary = "summary_mr_mc.txt"
        df_summary = summary(self.df, to_test_results(self.db_results), single_row=False, single_column=False,
                             convert_rt=None, ndigits_mz=None)

        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results), "summary"),
                                 sqlite_records(to_test_data(self.db_results), "summary"))

        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        self.assertSequenceEqual(df_summary.shape, (118, 28))
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_sr_mc.txt"
        df_summary = summary(self.df, to_test_results(self.db_results), single_row=True, single_column=False,
                             convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        self.assertSequenceEqual(df_summary.shape, (17, 17))
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_sr_sc.txt"
        df_summary = summary(self.df, to_test_results(self.db_results), single_row=True, single_column=True,
                             convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        self.assertSequenceEqual(df_summary.shape, (17, 12))
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_mr_mc_graphs.txt"
        df_summary = summary(self.df, to_test_results(self.db_results_graph), single_row=False, single_column=False,
                             convert_rt=None, ndigits_mz=None)
        self.assertSequenceEqual(sqlite_records(to_test_results(self.db_results_graph), "summary"),
                                 sqlite_records(to_test_data(self.db_results_graph), "summary"))
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        self.assertSequenceEqual(df_summary.shape, ((51, 32)))
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))

        fn_summary = "summary_mr_mc_graphs_nls.txt"
        db_nls = "results_annotation_nls.sqlite"
        df_nls = combine_peaklist_matrix(to_test_data("peaklist_lcms_pos_theoretical_nls.txt"),
                                         to_test_data("dataMatrix_lcms_theoretical_nls.txt"))
        df_summary = summary(df_nls, to_test_results(db_nls), single_row=False, single_column=False,
                             convert_rt=None, ndigits_mz=None)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)
        df_summary.to_csv(to_test_results(fn_summary), sep="\t", index=False)

        self.assertSequenceEqual(df_summary.shape, ((18, 27)))
        _assert(to_test_data(fn_summary), to_test_results(fn_summary))


if __name__ == '__main__':
    unittest.main()