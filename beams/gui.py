#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from functools import partial
import networkx as nx
from beams import in_out
from beams import grouping
from beams import annotation
from PyQt5 import QtCore, QtGui, QtWidgets
from beams.qt import form
import sqlite3


class BeamsApp(QtWidgets.QMainWindow, form.Ui_MainWindow):
    def __init__(self, parent=None):
        super(BeamsApp, self).__init__(parent)
        self.setupUi(self)

        self.pushButton_cancel.clicked.connect(QtCore.QCoreApplication.instance().quit)

        self.pushButton_peaklist.clicked.connect(partial(self.open_file, self.lineEdit_peaklist))
        self.pushButton_peak_matrix.clicked.connect(partial(self.open_file, self.lineEdit_intensity_matrix))

        self.pushButton_sql_database.clicked.connect(partial(self.save_file, self.lineEdit_sql_database, "database.sqlite"))
        self.pushButton_graph.clicked.connect(partial(self.save_file, self.lineEdit_graph, "graph.gml"))

        self.pushButton_default_adduct_library.clicked.connect(partial(self.open_file, self.lineEdit_default_adduct_library))
        self.pushButton_adduct_library.clicked.connect(partial(self.open_file, self.lineEdit_adduct_library))
        self.pushButton_isotopes.clicked.connect(partial(self.open_file, self.lineEdit_isotopes))
        self.pushButton_multiple_charged.clicked.connect(partial(self.open_file, self.lineEdit_multiple_charged))

        self.checkBox_filename_reference.clicked.connect(self.source_compounds)
        self.pushButton_filename_reference.clicked.connect(partial(self.open_file, self.lineEdit_filename_reference))
        self.pushButton_filename_mf.clicked.connect(partial(self.open_file, self.lineEdit_filename_mf))

        self.pushButton_summary_filename.clicked.connect(partial(self.save_file, self.lineEdit_summary_filename, "summary.txt"))

        self.checkBox_group_features.clicked.connect(self.group_features)
        self.checkBox_annotate_peak_patterns.clicked.connect(self.annotate_peak_patterns)
        self.checkBox_annotate_molecular_formulae.clicked.connect(self.annotate_molecular_formulae)
        self.checkBox_annotate_compounds.clicked.connect(self.annotate_compounds)
        self.checkBox_create_summary.clicked.connect(self.create_summary)

        self.comboBox_source_mf.activated.connect(self.source_mf)
        self.checkBox_adduct_library.clicked.connect(self.source_peak_patterns)
        self.checkBox_isotopes.clicked.connect(self.source_peak_patterns)
        self.checkBox_multiple_charged.clicked.connect(self.source_peak_patterns)
        self.checkBox_oligomers.clicked.connect(self.source_peak_patterns)

        self.checkBox_mz_digits.clicked.connect(self.create_summary)
        self.checkBox_convert_rt.clicked.connect(self.create_summary)

        self.add_databases()

        self.pushButton_start.clicked.connect(self.run)  # When the button is pressed


    def open_file(self, field):

        d = QtWidgets.QFileDialog.getOpenFileName(self, 'Select File', "")
        if d:
            if str(d[0]) == "":
                QtWidgets.QMessageBox.critical(None, "Select File", "No file selected", QtWidgets.QMessageBox.Ok)
            else:
                field.setText(d[0])
        return

    def save_file(self, field, filename):
        d = QtWidgets.QFileDialog.getSaveFileName(self, 'Save File', filename)
        if d:
            if str(d[0]) == "":
                QtWidgets.QMessageBox.critical(None, "Save File", "Provide a valid filename", QtWidgets.QMessageBox.Ok)
            else:
                field.setText(d[0])
        return

    def source_graph_file(self):
        if not self.checkBox_group_features.isChecked() and not self.checkBox_annotate_peak_patterns.isChecked():
            self.lineEdit_graph.setEnabled(False)
            self.pushButton_graph.setEnabled(False)
            self.label_graph.setEnabled(False)
        else:
            self.lineEdit_graph.setEnabled(True)
            self.pushButton_graph.setEnabled(True)
            self.label_graph.setEnabled(True)

    def source_default_adducts(self):
        self.lineEdit_default_adduct_library.setEnabled(True)
        self.pushButton_default_adduct_library.setEnabled(True)

    def source_mf(self):
        if self.comboBox_source_mf.currentText() == "Tab-delimited text file":
            self.label_filename_mf.setEnabled(True)
            self.lineEdit_filename_mf.setEnabled(True)
            self.pushButton_filename_mf.setEnabled(True)
            self.label_max_mz.setEnabled(False)
            self.spinBox_max_mz.setEnabled(False)
            self.checkBox_heuristic_rules.setEnabled(False)
        else:
            self.label_filename_mf.setEnabled(False)
            self.lineEdit_filename_mf.setEnabled(False)
            self.pushButton_filename_mf.setEnabled(False)
            self.label_max_mz.setEnabled(True)
            self.spinBox_max_mz.setEnabled(True)
            self.checkBox_heuristic_rules.setEnabled(True)

    def source_peak_patterns(self):
        if not self.checkBox_adduct_library.isChecked():
            self.lineEdit_adduct_library.setEnabled(False)
            self.pushButton_adduct_library.setEnabled(False)
        else:
            self.lineEdit_adduct_library.setEnabled(True)
            self.pushButton_adduct_library.setEnabled(True)
        if not self.checkBox_isotopes.isChecked():
            self.lineEdit_isotopes.setEnabled(False)
            self.pushButton_isotopes.setEnabled(False)
        else:
            self.lineEdit_isotopes.setEnabled(True)
            self.pushButton_isotopes.setEnabled(True)
        if not self.checkBox_multiple_charged.isChecked():
            self.pushButton_multiple_charged.setEnabled(False)
            self.lineEdit_multiple_charged.setEnabled(False)
        else:
            self.pushButton_multiple_charged.setEnabled(True)
            self.lineEdit_multiple_charged.setEnabled(True)
        if not self.checkBox_oligomers.isChecked():
            self.spinBox_max_monomer_units.setEnabled(False)
            self.label_max_monomer_units.setEnabled(False)
        else:
            self.spinBox_max_monomer_units.setEnabled(True)
            self.label_max_monomer_units.setEnabled(True)

    def source_compounds(self):
        if self.checkBox_filename_reference.isChecked():
            self.listWidget_databases.setEnabled(False)
            self.listWidget_categories.setEnabled(False)
            self.label_databases.setEnabled(False)
            self.pushButton_filename_reference.setEnabled(True)
            self.lineEdit_filename_reference.setEnabled(True)
        else:
            self.label_databases.setEnabled(True)
            self.listWidget_databases.setEnabled(True)
            self.label_categories.setEnabled(False)
            self.listWidget_categories.setEnabled(False)
            self.pushButton_filename_reference.setEnabled(False)
            self.lineEdit_filename_reference.setEnabled(False)

    def group_features(self):
        if not self.checkBox_group_features.isChecked():
            self.label_max_rt.setEnabled(False)
            self.doubleSpinBox_max_rt.setEnabled(False)
            self.comboBox_grouping_method.setEnabled(False)
            self.doubleSpinBox_p_value.setEnabled(False)
            self.label_tool_p_value.setEnabled(False)
            self.label_grouping_method.setEnabled(False)
            self.doubleSpinBox_coefficent.setEnabled(False)
            self.label_tool_coefficient.setEnabled(False)
        else:
            self.label_max_rt.setEnabled(True)
            self.doubleSpinBox_max_rt.setEnabled(True)
            self.comboBox_grouping_method.setEnabled(True)
            self.doubleSpinBox_p_value.setEnabled(True)
            self.label_tool_p_value.setEnabled(True)
            self.label_grouping_method.setEnabled(True)
            self.doubleSpinBox_coefficent.setEnabled(True)
            self.label_tool_coefficient.setEnabled(True)
        self.source_graph_file()

    def annotate_peak_patterns(self):
        if not self.checkBox_annotate_peak_patterns.isChecked():
            self.lineEdit_adduct_library.setEnabled(False)
            self.pushButton_adduct_library.setEnabled(False)
            self.checkBox_adduct_library.setEnabled(False)
            self.checkBox_isotopes.setEnabled(False)
            self.pushButton_multiple_charged.setEnabled(False)
            self.lineEdit_multiple_charged.setEnabled(False)
            self.lineEdit_isotopes.setEnabled(False)
            self.pushButton_isotopes.setEnabled(False)
            self.checkBox_multiple_charged.setEnabled(False)
            self.checkBox_oligomers.setEnabled(False)
            self.label_max_monomer_units.setEnabled(False)
            self.spinBox_max_monomer_units.setEnabled(False)
        else:
            self.checkBox_adduct_library.setEnabled(True)
            self.checkBox_isotopes.setEnabled(True)
            self.checkBox_multiple_charged.setEnabled(True)
            self.checkBox_oligomers.setEnabled(True)
            self.source_peak_patterns()
        self.source_graph_file()

    def annotate_molecular_formulae(self):
        if not self.checkBox_annotate_molecular_formulae.isChecked():
            self.label_filename_mf.setEnabled(False)
            self.lineEdit_filename_mf.setEnabled(False)
            self.comboBox_source_mf.setEnabled(False)
            self.label_source_mf.setEnabled(False)
            self.pushButton_filename_mf.setEnabled(False)
            self.label_max_mz.setEnabled(False)
            self.spinBox_max_mz.setEnabled(False)
            self.checkBox_heuristic_rules.setEnabled(False)
        else:
            self.comboBox_source_mf.setEnabled(True)
            self.label_source_mf.setEnabled(True)
            self.source_mf()
        return

    def annotate_compounds(self):
        if not self.checkBox_annotate_compounds.isChecked():
            self.listWidget_databases.setEnabled(False)
            self.listWidget_categories.setEnabled(False)
            self.label_databases.setEnabled(False)
            self.checkBox_filename_reference.setEnabled(False)
            self.pushButton_filename_reference.setEnabled(False)
            self.lineEdit_filename_reference.setEnabled(False)
        else:
            self.checkBox_filename_reference.setEnabled(True)
            self.source_compounds()

    def create_summary(self):
        if not self.checkBox_create_summary.isChecked():
            self.label_summary_filename.setEnabled(False)
            self.label_separator.setEnabled(False)
            self.checkBox_mz_digits.setEnabled(False)
            self.checkBox_convert_rt.setEnabled(False)
            self.label_annotations_format.setEnabled(False)
            self.comboBox_annotations_format.setEnabled(False)
            self.pushButton_summary_filename.setEnabled(False)
            self.comboBox_separator.setEnabled(False)
            self.comboBox_convert_rt.setEnabled(False)
            self.spinBox_mz_digits.setEnabled(False)
        else:
            self.label_summary_filename.setEnabled(True)
            self.label_separator.setEnabled(True)
            self.checkBox_mz_digits.setEnabled(True)
            self.checkBox_convert_rt.setEnabled(True)
            self.label_annotations_format.setEnabled(True)
            self.comboBox_annotations_format.setEnabled(True)
            self.pushButton_summary_filename.setEnabled(True)
            self.comboBox_separator.setEnabled(True)

            if self.checkBox_convert_rt.isChecked():
                self.comboBox_convert_rt.setEnabled(True)
            else:
                self.comboBox_convert_rt.setEnabled(False)
            if self.checkBox_mz_digits.isChecked():
                self.spinBox_mz_digits.setEnabled(True)
            else:
                self.spinBox_mz_digits.setEnabled(False)

    def add_databases(self):
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
        dbs = [fn for fn in os.listdir(path) if fn.endswith('.sqlite')]
        conn = sqlite3.connect(os.path.join(path, dbs[0]))
        cursor = conn.cursor()
        cursor.execute("""SELECT name FROM sqlite_master where type='table'""")
        db_names = [str(name[0]) for name in cursor.fetchall()]
        for i, db in enumerate(dbs[1:]):
            cursor.execute("ATTACH DATABASE ? AS db?", (os.path.join(path, dbs[i+1]), i, ))
            cursor.execute("""SELECT name FROM db?.sqlite_master where type='table'""", (i, ))
            db_names.extend([str(name[0]) for name in cursor.fetchall()])
        self.listWidget_databases.setSelectionMode(QtWidgets.QListWidget.MultiSelection)
        for db_name in db_names:
            item = QtWidgets.QListWidgetItem(db_name)
            self.listWidget_databases.addItem(item)
        return db_names

    def run(self):

        if not os.path.isfile(self.lineEdit_peaklist.text()) or not os.path.isfile(self.lineEdit_intensity_matrix.text()):
            QtWidgets.QMessageBox.critical(None, "Select file", "Select file(s) for Peaklist and/or Intensity Matrix", QtWidgets.QMessageBox.Ok)
            return

        elif self.lineEdit_sql_database.text() == "":
            QtWidgets.QMessageBox.critical(None, "Select File", "Select file for SQLite database to save output", QtWidgets.QMessageBox.Ok)
            return

        if self.checkBox_annotate_compounds.isChecked():
            if len(self.listWidget_databases.selectedItems()) == 0 and not self.checkBox_filename_reference.isChecked():
                QtWidgets.QMessageBox.critical(None, "Select File", "Select database or file for 'Annotate Compounds / Metabolites'", QtWidgets.QMessageBox.Ok)
                return

        if self.checkBox_create_summary.isChecked() and self.lineEdit_summary_filename.text() == "":
            QtWidgets.QMessageBox.critical(None, "Save File As", "Select file to save summary", QtWidgets.QMessageBox.Ok)
            return

        self.hide()

        lib_ion_mode = {"Positive": "pos", "Negative": "neg"}

        if self.checkBox_group_features.isChecked():
            print("Grouping features....")
            if self.comboBox_grouping_method.currentText() == "Pearson correlation":
                method = "pearson"
            else:
                method = "spearman"

            df = in_out.combine_peaklist_matrix(self.lineEdit_peaklist.text(), self.lineEdit_intensity_matrix.text())
            graph = grouping.group_features(df,
                                            db_out=self.lineEdit_sql_database.text(),
                                            max_rt_diff=self.doubleSpinBox_max_rt.value(),
                                            coeff_thres=self.doubleSpinBox_coefficent.value(),
                                            pvalue_thres=self.doubleSpinBox_p_value.value(),
                                            method=method,
                                            ncpus=None)
            nx.write_gml(graph, str(self.lineEdit_graph.text()))
            print("Done")
            print("")
        if self.checkBox_annotate_peak_patterns.isChecked():
            print("Annotating peak patterns....")
            if str(self.lineEdit_graph.text()) != "":
                inp = nx.read_gml(str(self.lineEdit_graph.text()))
            else:
                inp = in_out.combine_peaklist_matrix(self.lineEdit_peaklist.text(), self.lineEdit_intensity_matrix.text())

            if self.checkBox_adduct_library.isChecked():

                print("Adducts...."),

                if self.lineEdit_adduct_library.text() == "Use default":
                    path = 'data/adducts.txt'
                    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                    lib = in_out.read_adducts(p, lib_ion_mode[self.comboBox_ion_mode.currentText()])

                elif os.path.isfile(self.lineEdit_adduct_library.text()):
                    try:
                        lib = in_out.read_adducts(self.lineEdit_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                    except:
                        lib = in_out.read_mass_differences(self.lineEdit_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                    else:
                        QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)
                else:
                    QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts or 'Use default'", QtWidgets.QMessageBox.Ok)
                annotation.annotate_adducts(inp, db_out=self.lineEdit_sql_database.text(), ppm=self.doubleSpinBox_ppm_error.value(), lib=lib)
                print("Done")

            if self.checkBox_isotopes.isChecked():
                print("Isotopes...."),
                if self.lineEdit_isotopes.text() == "Use default":
                    path = 'data/isotopes.txt'
                    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                    lib = in_out.read_isotopes(p, lib_ion_mode[self.comboBox_ion_mode.currentText()])

                elif os.path.isfile(self.lineEdit_isotopes.text()):
                    lib = in_out.read_isotopes(self.lineEdit_isotopes.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                else:
                    QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for isotopes or 'Use default'", QtWidgets.QMessageBox.Ok)
                annotation.annotate_isotopes(inp, db_out=self.lineEdit_sql_database.text(), ppm=self.doubleSpinBox_ppm_error.value(), lib=lib)
                print("Done")

            if self.checkBox_multiple_charged.isChecked():
                print("Multiple charged ions...."),
                if self.lineEdit_multiple_charged.text() == "Use default":
                    path = 'data/multiple_charged_ions.txt'
                    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                    lib = in_out.read_multiple_charged_ions(p, lib_ion_mode[self.comboBox_ion_mode.currentText()])
                elif os.path.isfile(self.lineEdit_multiple_charged.text()):
                    lib = in_out.read_multiple_charged_ions(self.lineEdit_multiple_charged.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                else:
                    QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for multiple charged ions or 'Use default'", QtWidgets.QMessageBox.Ok)
                annotation.annotate_multiple_charged_ions(inp, db_out=self.lineEdit_sql_database.text(), ppm=self.doubleSpinBox_ppm_error.value(), lib=lib)
                print("Done")

            if self.checkBox_oligomers.isChecked():
                print("Oligomers...."),
                if self.lineEdit_default_adduct_library.text() == "Use default":
                    path = 'data/adducts.txt'
                    p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                    lib = in_out.read_adducts(p, lib_ion_mode[self.comboBox_ion_mode.currentText()])
                elif os.path.isfile(self.lineEdit_default_adduct_library.text()):
                    try:
                        lib = in_out.read_adducts(self.lineEdit_default_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                    except:
                        lib = in_out.read_mass_differences(self.lineEdit_default_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                    else:
                        QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)
                else:
                    QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)
                inp = in_out.combine_peaklist_matrix(self.lineEdit_peaklist.text(), self.lineEdit_intensity_matrix.text())
                annotation.annotate_oligomers(inp, db_out=self.lineEdit_sql_database.text(), ppm=self.doubleSpinBox_ppm_error.value(), lib=lib, maximum=self.spinBox_max_monomer_units.value())
                print("Done")
            print

        if self.checkBox_annotate_molecular_formulae.isChecked():
            print("Annotating molecular formulae...."),
            df = in_out.combine_peaklist_matrix(self.lineEdit_peaklist.text(), self.lineEdit_intensity_matrix.text())
            if self.lineEdit_default_adduct_library.text() == "Use default":
                path = 'data/adducts.txt'
                p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                lib = in_out.read_adducts(p, lib_ion_mode[self.comboBox_ion_mode.currentText()])
            elif os.path.isfile(self.lineEdit_default_adduct_library.text()):
                try:
                    lib = in_out.read_adducts(self.lineEdit_default_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                except:
                    lib = in_out.read_mass_differences(self.lineEdit_default_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                else:
                    QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)
            else:
                QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)

            if self.comboBox_source_mf.currentText() == "Tab-delimited text file":
                db_in = self.lineEdit_filename_mf.text()
                rules = None
                max_mz = None
            else:
                db_in = "http://multiomics-int.cs.bham.ac.uk"
                rules = self.checkBox_heuristic_rules.isChecked()
                max_mz = self.spinBox_max_mz.value()

            annotation.annotate_molecular_formulae(df,
                                                   lib_adducts=lib,
                                                   ppm=self.doubleSpinBox_ppm_error.value(),
                                                   db_out=self.lineEdit_sql_database.text(),
                                                   db_in=db_in,
                                                   rules=rules,
                                                   max_mz=max_mz)
            print("Done")
            print("")
        if self.checkBox_annotate_compounds.isChecked():
            print("Annotating compounds...."),
            df = in_out.combine_peaklist_matrix(self.lineEdit_peaklist.text(), self.lineEdit_intensity_matrix.text())

            if self.lineEdit_default_adduct_library.text() == "Use default":
                path = 'data/adducts.txt'
                p = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                lib = in_out.read_adducts(p, lib_ion_mode[self.comboBox_ion_mode.currentText()])
            elif os.path.isfile(self.lineEdit_default_adduct_library.text()):
                try:
                    lib = in_out.read_adducts(self.lineEdit_default_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                except:
                    lib = in_out.read_mass_differences(self.lineEdit_default_adduct_library.text(), lib_ion_mode[self.comboBox_ion_mode.currentText()])
                else:
                    QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)
            else:
                QtWidgets.QMessageBox.critical(None, "Select file", "Provide a valid filename for adducts", QtWidgets.QMessageBox.Ok)

            if self.checkBox_filename_reference.isChecked():
                annotation.annotate_compounds(df, lib_adducts=lib, ppm=self.doubleSpinBox_ppm_error.value(), db_out=self.lineEdit_sql_database.text(), db_in=self.lineEdit_filename_reference.text(), db_name=None)
            else:
                path = 'data/BEAMS_DB.sqlite'
                path_db = os.path.join(os.path.dirname(os.path.abspath(__file__)), path)
                print(path_db)
                for db_name in self.listWidget_databases.selectedItems():
                    annotation.annotate_compounds(df, lib_adducts=lib, ppm=self.doubleSpinBox_ppm_error.value(), db_out=self.lineEdit_sql_database.text(), db_in=path_db, db_name=db_name.text())
            print("Done")
            print

        if self.checkBox_create_summary.isChecked():
            print("Creating summary...."),
            if self.checkBox_convert_rt.isChecked():
                lib = {"Seconds": "sec", "Minutes": "min"}
                convert_rt = lib[self.comboBox_convert_rt.currentText()]
            else:
                convert_rt = None

            if self.checkBox_mz_digits.isChecked():
                ndigits_mz = self.spinBox_mz_digits.value()
            else:
                ndigits_mz = None

            df = in_out.combine_peaklist_matrix(self.lineEdit_peaklist.text(),
                                                self.lineEdit_intensity_matrix.text())

            if self.comboBox_annotations_format.currentText() == "Single row for each feature and separate columns":
                single_row = True
                single_column = False
            elif self.comboBox_annotations_format.currentText() == "Single row for each feature and merged columns":
                single_row = True
                single_column = True
            else:
                single_row = False
                single_column = False

            df_out = annotation.summary(df,
                                        db=self.lineEdit_sql_database.text(),
                                        single_row=single_row,
                                        single_column=single_column,
                                        convert_rt=convert_rt,
                                        ndigits_mz=ndigits_mz)

            separators = {"tab": "\t", "comma": ","}
            df_out.to_csv(self.lineEdit_summary_filename.text(), sep=separators[self.comboBox_separator.currentText()], index=False)
            print("Done")
            print("")

        self.close()