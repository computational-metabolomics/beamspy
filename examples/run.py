#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from beamspy import in_out
import networkx as nx
from beamspy.grouping import group_features
from beamspy.annotation import annotate_adducts
from beamspy.annotation import annotate_isotopes
from beamspy.annotation import annotate_oligomers
from beamspy.annotation import annotate_multiple_charged_ions
from beamspy.annotation import annotate_compounds
from beamspy.annotation import annotate_molecular_formulae
from beamspy.annotation import summary
from beamspy import plots


def main():

    path = "../tests/test_data/"
    fn_peaklist = os.path.join(path, "variableMetadata.txt")
    fn_matrix = os.path.join(path, "dataMatrix.txt")

    df = in_out.combine_peaklist_matrix(fn_peaklist, fn_matrix)

    ion_mode = "pos"

    db_out = "results_{}.sqlite".format(ion_mode)

    graphs = group_features(df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=0.01, method="pearson")

    nx.write_gml(graphs, "graphs.gml")
    # graphs = nx.read_gml("graphs.gml")

    path = "../beamspy/data"
    lib_isotopes = in_out.read_isotopes(os.path.join(path, "isotopes.txt"), ion_mode)
    lib_adducts = in_out.read_adducts(os.path.join(path, "adducts.txt"), ion_mode)
    lib_multiple_charged_ions = in_out.read_multiple_charged_ions(os.path.join(path, "multiple_charged_ions.txt"), ion_mode)
    lib_mass_differences = in_out.read_mass_differences(os.path.join(path, "multiple_charged_differences.txt"), ion_mode)

    print(lib_isotopes)
    print(lib_adducts)

    ppm = 5.0

    annotate_adducts(graphs, db_out, ppm, lib_adducts)
    annotate_isotopes(graphs, db_out, ppm, lib_isotopes)
    annotate_oligomers(graphs, db_out, ppm, lib_adducts)
    annotate_multiple_charged_ions(graphs, db_out, ppm, lib_multiple_charged_ions)

    # annotate_molecular_formulae(df, lib_adducts, ppm, db_out)
    annotate_compounds(df, lib_adducts, ppm, db_out, "lipidmaps_full_20181217_v1")

    df_out = summary(df, db_out)
    fn_out = "summary_{}.txt".format(ion_mode)
    df_out.to_csv(fn_out, sep="\t", index=False, encoding="utf-8")

    pdf_out = "report_{}.pdf".format(ion_mode)
    plots.report(db=db_out, pdf_out=pdf_out, column_corr="r_value", column_pvalue="p_value",
                 column_ppm_error="ppm_error", column_adducts="adduct")


if __name__ == '__main__':
    main()
