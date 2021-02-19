#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from beamspy import in_out
import networkx as nx
from beamspy.grouping import group_features
from beamspy.annotation import annotate_adducts
from beamspy.annotation import annotate_isotopes
from beamspy.annotation import annotate_oligomers
from beamspy.annotation import annotate_compounds
from beamspy.annotation import annotate_molecular_formulae
from beamspy.annotation import summary
from beamspy import plots


def main():

    path = "../tests/test_data/"
    fn_peaklist = os.path.join(path, "peaklist_lcms_pos_theoretical.txt")
    fn_matrix = os.path.join(path, "dataMatrix_lcms_theoretical.txt")

    df = in_out.combine_peaklist_matrix(fn_peaklist, fn_matrix)

    ion_mode = "pos"
    db_out = "results.sqlite".format(ion_mode)

    graphs = group_features(df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=0.01, method="pearson")

    nx.write_gml(graphs, "graphs.gml")
    # graphs = nx.read_gml("graphs.gml")

    path = "../beamspy/data"
    lib_isotopes = in_out.read_isotopes(os.path.join(path, "isotopes.txt"), ion_mode)
    lib_adducts = in_out.read_adducts(os.path.join(path, "adducts.txt"), ion_mode)

    print(lib_isotopes)
    print(lib_adducts)

    ppm = 5.0

    annotate_adducts(graphs, db_out, ppm, lib_adducts)
    annotate_isotopes(graphs, db_out, ppm, lib_isotopes)

    # annotate_molecular_formulae(df, lib_adducts, ppm, db_out)
    annotate_compounds(df, lib_adducts, ppm, db_out, "hmdb_full_v4_0_20200909_v1")

    df_out = summary(df, db_out)
    fn_out = "summary.txt"
    df_out.to_csv(fn_out, sep="\t", index=False, encoding="utf-8")

    pdf_out = "report.pdf"
    plots.report(db=db_out, pdf_out=pdf_out, column_corr="r_value", column_pvalue="p_value",
                 column_ppm_error="ppm_error", column_adducts="adduct")


if __name__ == '__main__':
    main()
