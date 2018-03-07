#!/usr/bin/python
# -*- coding: utf-8 -*-

import argparse
import networkx as nx
from . import __version__
import in_out
import grouping
import annotation


def map_delimiter(delimiter):
    seps = {"comma": ",", "tab": "\t"}
    if delimiter in seps:
        return seps[delimiter]
    else:
        return delimiter


def main():
    print("Executing BEAMS version %s." % __version__)

    parser = argparse.ArgumentParser(description='Annotation package of LC-MS and DIMS data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')


    parser_fi = subparsers.add_parser('format-inputs', help='Convert and format input files.')

    parser_gf = subparsers.add_parser('group-features', help='Group features.')

    parser_app = subparsers.add_parser('annotate-peak-patterns', help='Annotate peak patterns, molecular formulae and metabolites')

    parser_amf = subparsers.add_parser('annotate-mf', help='Annotate molecular formulae')

    parser_am = subparsers.add_parser('annotate-compounds', help='Annotate metabolites')


    parser_sr = subparsers.add_parser('summary-results', help='Summarise results')


    #################################
    # FORMAT INPUTS
    #################################

    parser_fi.add_argument('-l', '--peaklist',
                           type=str, required=True, help="Tab-delimited peaklist")

    parser_fi.add_argument('-m', '--intensity-matrix',
                           type=str, required=True, help="Tab-delimited intensity matrix")

    parser_fi.add_argument('-o', '--output', type=str, required=True,
                           help="Filename / Path to save the merged peaklist and intensity matrix.")

    parser_fi.add_argument('-s', '--sep', default="tab", choices=["tab", "comma"], required=True,
                           help="Values on each line of the file are separated by this character.")


    #################################
    # GROUP FEATURES
    #################################

    parser_gf.add_argument('-l', '--peaklist',
                           type=str, required=True, help="Tab-delimited peaklist")

    parser_gf.add_argument('-i', '--intensity-matrix',
                           type=str, required=False, help="Tab-delimited intensity matrix")

    #parser_gf.add_argument('-x', '--xset-matrix',
    #                       type=str, required=False, help="Tab-delimited intensity matrix")

    parser_gf.add_argument('-d', '--db', type=str, required=True,
                           help="Sqlite database to write results")

    parser_gf.add_argument('-r', '--max-rt-diff', default=5.0, type=float, required=True,
                           help="Maximum difference in retention time between two peaks")

    parser_gf.add_argument('-m', '--method', default="pearson", choices=["pearson", "spearman"], required=True,
                           help="Method to apply for grouping features")

    parser_gf.add_argument('-c', '--coeff-threshold', default=0.7, type=float, required=True,
                           help="Threshold for correlation coefficient")

    parser_gf.add_argument('-p', '--pvalue-threshold', default=0.01, type=float, required=True,
                           help="Threshold for p-value")

    parser_gf.add_argument('-g', '--gml-file', type=str, required=True,
                           help="Tab-delimited intensity matrix")


    #################################
    # ANNOTATE PEAK PATTERS
    #################################

    parser_app.add_argument('-l', '--peaklist', type=str, required=True,
                             help="Tab-delimited peaklist")

    parser_app.add_argument('-i', '--intensity-matrix', type=str, required=False,
                             help="Tab-delimited intensity matrix")

    parser_app.add_argument('-g', '--gml-file', type=str, required=False,
                             help="Tab-delimited intensity matrix")

    parser_app.add_argument('-d', '--db', type=str, required=True,
                             help="Sqlite database to write results")

    parser_app.add_argument('-a', '--adducts', action='store_true', required=False,
                             help="Annotate adducts")

    parser_app.add_argument('-b', '--adducts-library', action='append', required=False,
                             help="Annotate adducts")

    parser_app.add_argument('-e', '--isotopes', action='store_true', required=False,
                             help="Annotate isotopes")

    parser_app.add_argument('-f', '--isotopes-library', action='append', required=False,
                             help="List of isotopes")

    parser_app.add_argument('-r', '--multiple-charged-ions', action='store_true', required=False,
                             help="Annotate multiple-charged ions")

    parser_app.add_argument('-s', '--multiple-charged-ions-library', action='append', required=False,
                             help="List of multiple charged ions")

    parser_app.add_argument('-o', '--oligomers', action='store_true', required=False,
                             help="Annotate oligomers")

    parser_app.add_argument('-m', '--ion-mode', choices=["pos", "neg"], required=True,
                             help="Define the ion mode of the libraries")

    parser_app.add_argument('-p', '--ppm', default=3.0, type=float, required=True,
                             help="Mass tolerance in parts per million.")


    #################################
    # ANNOTATE MOLECULAR FORMULAE
    #################################

    parser_amf.add_argument('-l', '--peaklist', type=str, required=True,
                            help="Tab-delimited peaklist")

    parser_amf.add_argument('-i', '--intensity-matrix', type=str, required=False,
                            help="Tab-delimited intensity matrix")

    parser_amf.add_argument('-d', '--db', type=str, required=True,
                            help="Sqlite database to write results")

    parser_amf.add_argument('-c', '--db-mf', type=str, required=True, default="http://multiomics-int.cs.bham.ac.uk",
                            help="Molecular formulae database (reference)")

    parser_amf.add_argument('-a', '--adducts-library', required=True,
                            help="List of adducts to search for")

    parser_amf.add_argument('-m', '--ion-mode', choices=["pos", "neg"], required=True,
                             help="Define the ion mode of the libraries")

    parser_amf.add_argument('-p', '--ppm', default=3.0, type=float, required=True,
                            help="Mass tolerance in parts per million.")

    parser_amf.add_argument('-z', '--max-mz', type=float, required=False, default=None,
                            help="Maximum m/z value to assign molecular formula(e).")


    #################################
    # ANNOTATE METABOLITES
    #################################

    parser_am.add_argument('-l', '--peaklist', type=str, required=True,
                           help="Tab-delimited peaklist")

    parser_am.add_argument('-i', '--intensity-matrix', type=str, required=False,
                           help="Tab-delimited intensity matrix")

    parser_am.add_argument('-d', '--db', type=str, required=True,
                           help="Sqlite database to write results")

    parser_am.add_argument('-c', '--db-compounds', type=str, required=True, help="Metabolite database (reference)")

    parser_am.add_argument('-n', '--db-name', type=str, default="", required=False,
                           help="Name compound / metabolite database (within --db-compounds)")

    parser_am.add_argument('-a', '--adducts-library', required=True,
                           help="List of adducts to search for")

    parser_am.add_argument('-m', '--ion-mode', choices=["pos", "neg"], required=True,
                             help="Define the ion mode of the libraries")

    parser_am.add_argument('-p', '--ppm', default=3.0, type=float, required=True,
                           help="Mass tolerance in parts per million.")

    #################################
    # SUMMARY RESULTS
    #################################

    parser_sr.add_argument('-l', '--peaklist', type=str, required=True,
                           help="Tab-delimited peaklist")

    parser_sr.add_argument('-i', '--intensity-matrix', type=str, required=False,
                           help="Tab-delimited intensity matrix")

    parser_sr.add_argument('-o', '--output', type=str, required=False,
                           help="Summary file")

    parser_sr.add_argument('-d', '--db', type=str, required=True,
                           help="Sqlite database to write results")

    parser_sr.add_argument('-s', '--sep', default="tab", choices=["tab", "comma"], required=True,
                           help="Values on each line of the file are separated by this character.")

    args = parser.parse_args()

    print args

    separators = {"tab": "\t", "comma": ","}

    if args.step == "format-inputs":

        df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        df.to_csv(args.output, sep=separators[args.sep], index=False)

    if args.step == "group-features":
        df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        graph = grouping.group_features(df, db_out=args.db, max_rt_diff=args.max_rt_diff, coeff_thres=args.coeff_threshold,
                                        pvalue_thres=args.pvalue_threshold, method=args.method)
        print graph
        nx.write_gml(graph, str(args.gml_file))

    if args.step == "annotate-peak-patterns":

        if args.gml_file:
            inp = nx.read_gml(args.gml_file)
        else:
            inp = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)

        if args.adducts:
            for i, a in enumerate(args.adducts_library):
                try:
                    lib = in_out.read_adducts(a, args.ion_mode)
                except:
                    lib = in_out.read_mass_differences(a, args.ion_mode)
                if i > 0:
                    add = True
                else:
                    add = False
                annotation.annotate_adducts(inp, db_out=args.db, ppm=args.ppm, lib=lib, add=add)


        if args.isotopes:
            for i in args.isotopes_library:
                lib = in_out.read_isotopes(i, args.ion_mode)
                annotation.annotate_isotopes(inp, db_out=args.db, ppm=args.ppm, lib=lib)

        if args.multiple_charged_ions:
            for i, m in enumerate(args.multiple_charged_ions_library):
                try:
                    lib = in_out.read_multiple_charged_ions(m, args.ion_mode)
                except:
                    lib = in_out.read_mass_differences(m, args.ion_mode)

                if i > 0:
                    add = True
                else:
                    add = False

                annotation.annotate_multiple_charged_ions(inp, db_out=args.db, ppm=args.ppm, lib=lib, add=add)

        if args.oligomers:
            annotation.annotate_oligomers(inp, db_out=args.db, ppm=args.ppm, lib=lib)

    if args.step == "annotate-mf":
        df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        lib = in_out.read_adducts(args.adducts_library, args.ion_mode)
        annotation.annotate_molecular_formulae(df, ppm=args.ppm, lib_adducts=lib, db_out=args.db, db_in=args.db_mf, max_mz=args.max_mz)

    if args.step == "annotate-compounds":
        df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        lib = in_out.read_adducts(args.adducts_library, args.ion_mode)
        annotation.annotate_compounds(df, lib_adducts=lib, ppm=args.ppm, db_out=args.db, db_in=args.db_compounds, db_name=args.db_name)

    if args.step == "summary-results":
        df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        df_out = annotation.summary(df, db=args.db)
        df_out.to_csv(args.output, sep=separators[args.sep], index=False)


"""
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


    #annotate_oligomers(graphs, db_out, ppm, lib_adducts)
    #annotate_artifacts(testdata_no_MPA, db_out, 0.005)
    #annotate_multiple_charged_ions(graph, db_out, ppm, lib_multiple_charged_ions)

    annotate_molecular_formulae(df, lib_adducts, ppm, db_out, fn_mf)
    #annotate_molecular_formulae(df, lib_adducts, ppm, db_out)
    annotate_compounds(df, lib_adducts, ppm, db_out, fn_sql_db, "HMDB__25_09_2016")


"""


if __name__ == "__main__":
    main()