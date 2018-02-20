#!/usr/bin/python
# -*- coding: utf-8 -*-

from beams import in_out
from beams.annotation import annotate_adducts
from beams.annotation import annotate_isotopes
from beams.annotation import annotate_compounds
from beams.annotation import annotate_molecular_formulae
from beams.annotation import summary

def main():

    fn_peaklist = "../tests/test_data/variableMetadata.txt"
    fn_matrix = "../tests/test_data/dataMatrix.txt"

    #df = in_out.combine_peaklist_matrix(fn_peaklist, fn_matrix)

    fn_matrix = "../tests/test_data/Study_pospeaks_both_parts_no_MPA.txt"
    #df = in_out.read_xset_matrix(fn_matrix, first_sample="s1")

    fn_peaklist = "../tests/test_data/peaklist_dims_pos_theoretical.txt"
    df = in_out.read_peaklist(fn_peaklist)

    db_out = "../tests/test_data/results.sqlite"

    #df = df[0:700]
    #print df
    #graph = group_features(df, db_out, max_rt_diff=5.0, coeff_thres=0.7, pvalue_thres=None, method="pearson")

    #nx.write_gml(graph, "../tests/test_data/graph.gml")
    #graph = nx.read_gml("../tests/test_data/graph.gml")

    lib_isotopes = in_out.read_isotopes("../beams/data/isotopes.txt", "pos")
    lib_adducts = in_out.read_adducts("../beams/data/adducts.txt", "pos")
    lib_multiple_charged_ions = in_out.read_multiple_charged_ions("../beams/data/multiple_charged_ions.txt", "pos")
    #lib_mass_differences = in_out.read_mass_differences("data/multiple_charged_differences.txt", "pos")
    fn_mf = "../beams/data/trimMMD_sortAmass.txt"
    fn_cpds = "../beams/data/trimMMD_sortMF.txt"
    fn_sql_db = "../beams/data/MIDB_260916.sqlite"

    print lib_isotopes
    print lib_adducts
    print lib_multiple_charged_ions

    ppm = 5.0

    #annotate_adducts(graph, db_out, ppm, lib_adducts)
    annotate_adducts(df, db_out, ppm, lib_adducts)
    #annotate_isotopes(graph, db_out, ppm, lib_isotopes)
    #annotate_isotopes(df, db_out, ppm, lib_isotopes)
    #annotate_oligomers(graphs, db_out, ppm, lib_adducts)
    #annotate_artifacts(testdata_no_MPA, db_out, 0.005)
    #annotate_multiple_charged_ions(graph, db_out, ppm, lib_multiple_charged_ions)

    annotate_molecular_formulae(df, lib_adducts, ppm, db_out, fn_mf)
    #annotate_molecular_formulae(df, lib_adducts, ppm, db_out)
    annotate_compounds(df, lib_adducts, ppm, db_out, fn_sql_db, "HMDB__25_09_2016")

    df_out = summary(df, db_out)
    fn_out = "../tests/test_data/summary.txt"
    df_out.to_csv(fn_out, sep="\t", index=False)

if __name__ == '__main__':
    main()