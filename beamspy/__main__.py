#!/usr/bin/env python
# -*- coding: utf-8 -*-

from beamspy import __version__
import argparse
import sys
import os
import networkx as nx
from beamspy import in_out
from beamspy import grouping
from beamspy import annotation
from beamspy import plots


def map_delimiter(delimiter):
    seps = {"comma": ",", "tab": "\t"}
    if delimiter in seps:
        return seps[delimiter]
    else:
        return delimiter


def main():
    print("Executing BEAMSpy version {}.".format(__version__))

    parser = argparse.ArgumentParser(description='Annotation package of LC-MS and DIMS data',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
                                     # formatter_class=RawTextHelpFormatter)

    subparsers = parser.add_subparsers(dest='step')


    parser_gf = subparsers.add_parser('group-features', help='Group features.')

    parser_app = subparsers.add_parser('annotate-peak-patterns', help='Annotate peak patterns, molecular formulae and metabolites.')

    parser_amf = subparsers.add_parser('annotate-mf', help='Annotate molecular formulae.')

    parser_am = subparsers.add_parser('annotate-compounds', help='Annotate metabolites.')

    parser_sr = subparsers.add_parser('summary-results', help='Summarise results.')

    parser_gui = subparsers.add_parser('start-gui', help='Start GUI.')


    #################################
    # GROUP FEATURES
    #################################

    parser_gf.add_argument('-l', '--peaklist',
                           type=str, required=True, help="Tab-delimited peaklist.")

    parser_gf.add_argument('-i', '--intensity-matrix',
                           type=str, required=True, help="Tab-delimited intensity matrix.")

    #parser_gf.add_argument('-x', '--xset-matrix',
    #                       type=str, required=False, help="Tab-delimited intensity matrix")

    parser_gf.add_argument('-d', '--db', type=str, required=True,
                           help="Sqlite database to write results.")

    parser_gf.add_argument('-r', '--max-rt-diff', default=5.0, type=float, required=True,
                           help="Maximum difference in retention time between two peaks.")

    parser_gf.add_argument('-m', '--method', default="pearson", choices=["pearson", "spearman"], required=True,
                           help="Method to apply for grouping features.")

    parser_gf.add_argument('-c', '--coeff-threshold', default=0.7, type=float, required=True,
                           help="Threshold for correlation coefficient.")

    parser_gf.add_argument('-p', '--pvalue-threshold', default=0.01, type=float, required=True,
                           help="Threshold for p-value.")

    parser_gf.add_argument('-g', '--gml-file', type=str, required=True,
                           help="Write graph to GraphML format.")

    parser_gf.add_argument('-n', '--ncpus', type=int, required=False,
                           help="Number of central processing units (CPUs).")

    #################################
    # ANNOTATE PEAK PATTERS
    #################################

    parser_app.add_argument('-l', '--peaklist', type=str, required=True,
                             help="Tab-delimited peaklist.")

    parser_app.add_argument('-i', '--intensity-matrix', type=str, required=False,
                             help="Tab-delimited intensity matrix.")

    parser_app.add_argument('-g', '--gml-file', type=str, required=False,
                             help="Correlation graph in GraphML format.")

    parser_app.add_argument('-d', '--db', type=str, required=True,
                             help="Sqlite database to write results.")

    parser_app.add_argument('-a', '--adducts', action='store_true', required=False,
                             help="Annotate adducts.")

    parser_app.add_argument('-b', '--adducts-library', type=str, default=None, required=False,
                             help="List of adducts.")

    parser_app.add_argument('-e', '--isotopes', action='store_true', required=False,
                             help="Annotate isotopes.")

    parser_app.add_argument('-f', '--isotopes-library', required=False,
                             help="List of isotopes.")

    parser_app.add_argument('-o', '--oligomers', action='store_true', required=False,
                             help="Annotate oligomers.")

    parser_app.add_argument('-n', '--neutral-losses', action='store_true', required=False,
                             help="Annotate neutral losses.")

    parser_app.add_argument('-s', '--neutral-losses-library', required=False,
                             help="List of neutral losses.")

    parser_app.add_argument('-m', '--ion-mode', choices=["pos", "neg"], required=True,
                             help="Ion mode of the libraries.")

    parser_app.add_argument('-p', '--ppm', default=3.0, type=float, required=True,
                             help="Mass tolerance in parts per million.")

    parser_app.add_argument('-u', '--max-monomer-units', default=2, type=int, required=False,
                             help="Maximum number of monomer units.")


    #################################
    # ANNOTATE MOLECULAR FORMULAE
    #################################

    parser_amf.add_argument('-l', '--peaklist', type=str, required=True,
                            help="Tab-delimited peaklist.")

    parser_amf.add_argument('-i', '--intensity-matrix', type=str, required=False,
                            help="Tab-delimited intensity matrix.")

    parser_amf.add_argument('-d', '--db', type=str, required=True,
                            help="Sqlite database to write results.")

    parser_amf.add_argument('-c', '--db-mf', type=str, default="http://mfdb.bham.ac.uk",
                            help="Molecular formulae database (reference).")

    parser_amf.add_argument('-a', '--adducts-library', type=str, default=None, required=False,
                            help="List of adducts to search for.")

    parser_amf.add_argument('-m', '--ion-mode', choices=["pos", "neg"], required=True,
                             help="Ion mode of the libraries.")

    parser_amf.add_argument('-p', '--ppm', default=3.0, type=float, required=True,
                            help="Mass tolerance in parts per million.")

    parser_amf.add_argument('-e', '--skip-patterns', action="store_false",
                            help="Skip applying/using peak patterns (e.g. adduct and isotope patterns) to filter annotations.")

    parser_amf.add_argument('-r', '--skip-rules', action="store_false",
                            help="Skip heuritic rules to filter annotations.")

    parser_amf.add_argument('-z', '--max-mz', type=float, required=False, default=500.0,
                            help="Maximum m/z value to assign molecular formula(e).")


    #################################
    # ANNOTATE METABOLITES
    #################################

    parser_am.add_argument('-l', '--peaklist', type=str, required=True,
                           help="Tab-delimited peaklist.")

    parser_am.add_argument('-i', '--intensity-matrix', type=str, required=False,
                           help="Tab-delimited intensity matrix.")

    parser_am.add_argument('-d', '--db', type=str, required=True,
                           help="Sqlite database to write results.")

    parser_am.add_argument('-c', '--db-compounds', type=str, default="", required=False,
                           help="Metabolite database (reference).")

    parser_am.add_argument('-n', '--db-name', type=str, default="", required=True,
                           help="Name compound / metabolite database (within --db-compounds).")

    parser_am.add_argument('-a', '--adducts-library', type=str, default=None, required=False,
                           help="List of adducts to search for.")

    parser_am.add_argument('-m', '--ion-mode', choices=["pos", "neg"], required=True,
                             help="Ion mode of the libraries.")

    parser_am.add_argument('-p', '--ppm', default=3.0, type=float, required=True,
                           help="Mass tolerance in parts per million.")

    parser_am.add_argument('-e', '--skip-patterns', action="store_false",
                            help="Skip applying/using peak patterns (e.g. adduct and isotope patterns) to filter annotations.")

    parser_am.add_argument('-r', '--rt', default=None, type=float,
                           help="Retention time tolerance in seconds.")

    #################################
    # SUMMARY RESULTS
    #################################

    parser_sr.add_argument('-l', '--peaklist', type=str, required=True,
                           help="Tab-delimited peaklist")

    parser_sr.add_argument('-i', '--intensity-matrix', type=str, required=False,
                           help="Tab-delimited intensity matrix.")

    parser_sr.add_argument('-o', '--output', type=str, required=True,
                           help="Output file for the summary")

    parser_sr.add_argument('-p', '--pdf', type=str, required=False,
                           help="Output pdf file for the summary plots")

    parser_sr.add_argument('-d', '--db', type=str, required=True,
                           help="Sqlite database that contains the results from the previous steps.")

    parser_sr.add_argument('-s', '--sep', default="tab", choices=["tab", "comma"], required=True,
                           help="Values on each line of the output are separated by this character.")

    parser_sr.add_argument('-r', '--single-row', action="store_true",
                           help="Concatenate the annotations for each spectral feature and represent in a single row.")

    parser_sr.add_argument('-c', '--single-column', action="store_true",
                           help="Concatenate the annotations for each spectral feature and keep seperate columns for molecular formula, adduct, name, etc.")

    parser_sr.add_argument('-n', '--ndigits-mz', default=None, type=int, required=False,
                           help="Digits after the decimal point for m/z values.")

    parser_sr.add_argument('-t', '--convert-rt', default=None, choices=["sec", "min", None],
                           required=False, help="Covert the retention time to seconds or minutes. An additional column will be added.")

    args = parser.parse_args()

    print(args)

    separators = {"tab": "\t", "comma": ","}

    if args.step == "group-features":
        df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        graph = grouping.group_features(df, db_out=args.db, max_rt_diff=args.max_rt_diff,
                                        coeff_thres=args.coeff_threshold, pvalue_thres=args.pvalue_threshold,
                                        method=args.method, ncpus=args.ncpus)
        nx.write_gml(graph, str(args.gml_file))

    if args.step == "annotate-peak-patterns":

        if args.gml_file:
            inp = nx.read_gml(args.gml_file)
        elif args.intensity_matrix:
            inp = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        else:
            inp = in_out.read_peaklist(args.peaklist)

        if args.adducts:
            if args.adducts_library:
                lib = in_out.read_adducts(args.adducts_library, args.ion_mode)
            else:
                path = 'data/adducts.txt'
                p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                lib = in_out.read_adducts(p, args.ion_mode)
            annotation.annotate_adducts(inp, db_out=args.db, ppm=args.ppm, lib=lib, add=False)

        if args.isotopes:
            if args.isotopes_library:
                lib = in_out.read_isotopes(args.isotopes_library, args.ion_mode)
            else:
                path = 'data/isotopes.txt'
                p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                lib = in_out.read_isotopes(p, args.ion_mode)
            annotation.annotate_isotopes(inp, db_out=args.db, ppm=args.ppm, lib=lib)

        if args.neutral_losses:
            if args.neutral_losses_library:
                lib = in_out.read_neutral_losses(args.neutral_losses_library)
            else:
                path = 'data/neutral_losses.txt'
                p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                lib = in_out.read_neutral_losses(p)
            annotation.neutral_losses(inp, db_out=args.db, ppm=args.ppm, lib=lib)

        if args.oligomers:
            if args.adducts_library:
                lib = in_out.read_adducts(args.adducts_library, args.ion_mode)
            else:
                path = 'data/adducts.txt'
                p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                lib = in_out.read_adducts(p, args.ion_mode)

            annotation.annotate_oligomers(inp, db_out=args.db, ppm=args.ppm, lib=lib, maximum=args.max_monomer_units)

    if args.step == "annotate-mf":

        if args.intensity_matrix:
            df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        else:
            df = in_out.read_peaklist(args.peaklist)

        if args.adducts_library:
            lib = in_out.read_adducts(args.adducts_library, args.ion_mode)
        else:
            path = 'data/adducts.txt'
            p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
            lib = in_out.read_adducts(p, args.ion_mode)
        annotation.annotate_molecular_formulae(df, ppm=args.ppm, lib_adducts=lib, db_out=args.db, db_in=args.db_mf,
                                               patterns=args.skip_patterns, rules=args.skip_rules, max_mz=args.max_mz)

    if args.step == "annotate-compounds":
        
        if args.intensity_matrix:
            df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        else:
            df = in_out.read_peaklist(args.peaklist)

        if args.adducts_library:
            lib = in_out.read_adducts(args.adducts_library, args.ion_mode)
        else:
            path = 'data/adducts.txt'
            p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
            lib = in_out.read_adducts(p, args.ion_mode)
        annotation.annotate_compounds(df, lib_adducts=lib, ppm=args.ppm, db_out=args.db, db_name=args.db_name, patterns=args.skip_patterns, db_in=args.db_compounds, rt_tol=args.rt)

    if args.step == "summary-results":
        
        if args.intensity_matrix:
            df = in_out.combine_peaklist_matrix(args.peaklist, args.intensity_matrix)
        else:
            df = in_out.read_peaklist(args.peaklist)

        df_out = annotation.summary(df, db=args.db, single_row=args.single_row, single_column=args.single_column, convert_rt=args.convert_rt, ndigits_mz=args.ndigits_mz)
        df_out.to_csv(args.output, sep=separators[args.sep], index=False, encoding="utf-8")
        if args.pdf:
            plots.report(db=args.db, pdf_out=args.pdf,
                         column_corr="r_value", column_pvalue="p_value",
                         column_ppm_error="ppm_error", column_adducts="adduct")

    if args.step == "start-gui":
        from PySide2 import QtWidgets
        from beamspy.gui import BeamsApp
        app = QtWidgets.QApplication(sys.argv)
        # app.setStyle("Fusion")
        form = BeamsApp()
        form.show()
        sys.exit(app.exec_())


if __name__ == "__main__":
    main()
