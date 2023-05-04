# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'form.ui'
##
## Created by: Qt User Interface Compiler version 5.15.6
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide2.QtCore import *  # type: ignore
from PySide2.QtGui import *  # type: ignore
from PySide2.QtWidgets import *  # type: ignore


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        if not MainWindow.objectName():
            MainWindow.setObjectName(u"MainWindow")
        MainWindow.setEnabled(True)
        MainWindow.resize(974, 775)
        MainWindow.setAnimated(False)
        self.actionExampleData = QAction(MainWindow)
        self.actionExampleData.setObjectName(u"actionExampleData")
        self.actionAbout = QAction(MainWindow)
        self.actionAbout.setObjectName(u"actionAbout")
        self.centralwidget = QWidget(MainWindow)
        self.centralwidget.setObjectName(u"centralwidget")
        self.scrollArea = QScrollArea(self.centralwidget)
        self.scrollArea.setObjectName(u"scrollArea")
        self.scrollArea.setGeometry(QRect(6, 5, 961, 731))
        self.scrollArea.setWidgetResizable(True)
        self.scrollAreaWidgetContents = QWidget()
        self.scrollAreaWidgetContents.setObjectName(u"scrollAreaWidgetContents")
        self.scrollAreaWidgetContents.setGeometry(QRect(0, 0, 959, 729))
        self.verticalLayout_2 = QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout_2.setObjectName(u"verticalLayout_2")
        self.groupBox_general = QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_general.setObjectName(u"groupBox_general")
        self.gridLayout_3 = QGridLayout(self.groupBox_general)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.gridLayout_3.setVerticalSpacing(5)
        self.gridLayout_3.setContentsMargins(10, 5, 10, 5)
        self.pushButton_peaklist = QPushButton(self.groupBox_general)
        self.pushButton_peaklist.setObjectName(u"pushButton_peaklist")
        self.pushButton_peaklist.setEnabled(True)
        self.pushButton_peaklist.setMinimumSize(QSize(95, 23))
        self.pushButton_peaklist.setMaximumSize(QSize(95, 16777215))

        self.gridLayout_3.addWidget(self.pushButton_peaklist, 5, 2, 1, 1)

        self.pushButton_graph = QPushButton(self.groupBox_general)
        self.pushButton_graph.setObjectName(u"pushButton_graph")
        self.pushButton_graph.setEnabled(True)
        self.pushButton_graph.setMinimumSize(QSize(95, 23))
        self.pushButton_graph.setMaximumSize(QSize(95, 16777215))

        self.gridLayout_3.addWidget(self.pushButton_graph, 5, 6, 1, 1)

        self.pushButton_default_adduct_library = QPushButton(self.groupBox_general)
        self.pushButton_default_adduct_library.setObjectName(u"pushButton_default_adduct_library")
        self.pushButton_default_adduct_library.setEnabled(True)
        self.pushButton_default_adduct_library.setMinimumSize(QSize(95, 23))
        self.pushButton_default_adduct_library.setMaximumSize(QSize(95, 16777215))

        self.gridLayout_3.addWidget(self.pushButton_default_adduct_library, 7, 6, 1, 1)

        self.lineEdit_intensity_matrix = QLineEdit(self.groupBox_general)
        self.lineEdit_intensity_matrix.setObjectName(u"lineEdit_intensity_matrix")
        self.lineEdit_intensity_matrix.setMinimumSize(QSize(150, 0))
        self.lineEdit_intensity_matrix.setMaximumSize(QSize(150, 16777215))
        self.lineEdit_intensity_matrix.setReadOnly(True)

        self.gridLayout_3.addWidget(self.lineEdit_intensity_matrix, 7, 1, 1, 1)

        self.lineEdit_graph = QLineEdit(self.groupBox_general)
        self.lineEdit_graph.setObjectName(u"lineEdit_graph")
        self.lineEdit_graph.setEnabled(True)
        self.lineEdit_graph.setMinimumSize(QSize(150, 0))
        self.lineEdit_graph.setMaximumSize(QSize(150, 16777215))
        self.lineEdit_graph.setReadOnly(True)

        self.gridLayout_3.addWidget(self.lineEdit_graph, 5, 5, 1, 1)

        self.label_graph = QLabel(self.groupBox_general)
        self.label_graph.setObjectName(u"label_graph")
        self.label_graph.setEnabled(True)

        self.gridLayout_3.addWidget(self.label_graph, 5, 3, 1, 1)

        self.lineEdit_wd = QLineEdit(self.groupBox_general)
        self.lineEdit_wd.setObjectName(u"lineEdit_wd")
        self.lineEdit_wd.setEnabled(True)
        self.lineEdit_wd.setMinimumSize(QSize(150, 0))
        self.lineEdit_wd.setMaximumSize(QSize(150, 16777215))
        self.lineEdit_wd.setReadOnly(True)

        self.gridLayout_3.addWidget(self.lineEdit_wd, 4, 1, 1, 1)

        self.pushButton_sql_database = QPushButton(self.groupBox_general)
        self.pushButton_sql_database.setObjectName(u"pushButton_sql_database")
        self.pushButton_sql_database.setEnabled(True)
        self.pushButton_sql_database.setMinimumSize(QSize(95, 23))
        self.pushButton_sql_database.setMaximumSize(QSize(95, 16777215))

        self.gridLayout_3.addWidget(self.pushButton_sql_database, 4, 6, 1, 1)

        self.lineEdit_sql_database = QLineEdit(self.groupBox_general)
        self.lineEdit_sql_database.setObjectName(u"lineEdit_sql_database")
        self.lineEdit_sql_database.setEnabled(True)
        self.lineEdit_sql_database.setMinimumSize(QSize(150, 0))
        self.lineEdit_sql_database.setMaximumSize(QSize(150, 16777215))
        self.lineEdit_sql_database.setReadOnly(True)

        self.gridLayout_3.addWidget(self.lineEdit_sql_database, 4, 5, 1, 1)

        self.lineEdit_peaklist = QLineEdit(self.groupBox_general)
        self.lineEdit_peaklist.setObjectName(u"lineEdit_peaklist")
        self.lineEdit_peaklist.setEnabled(True)
        self.lineEdit_peaklist.setMinimumSize(QSize(150, 0))
        self.lineEdit_peaklist.setMaximumSize(QSize(150, 16777215))
        self.lineEdit_peaklist.setReadOnly(True)

        self.gridLayout_3.addWidget(self.lineEdit_peaklist, 5, 1, 1, 1)

        self.pushButton_wd = QPushButton(self.groupBox_general)
        self.pushButton_wd.setObjectName(u"pushButton_wd")
        self.pushButton_wd.setEnabled(True)
        self.pushButton_wd.setMinimumSize(QSize(95, 23))
        self.pushButton_wd.setMaximumSize(QSize(95, 16777215))

        self.gridLayout_3.addWidget(self.pushButton_wd, 4, 2, 1, 1)

        self.label_default_adduct_library = QLabel(self.groupBox_general)
        self.label_default_adduct_library.setObjectName(u"label_default_adduct_library")
        self.label_default_adduct_library.setEnabled(True)

        self.gridLayout_3.addWidget(self.label_default_adduct_library, 7, 3, 1, 1)

        self.label_intensity_matrix = QLabel(self.groupBox_general)
        self.label_intensity_matrix.setObjectName(u"label_intensity_matrix")

        self.gridLayout_3.addWidget(self.label_intensity_matrix, 7, 0, 1, 1)

        self.lineEdit_default_adduct_library = QLineEdit(self.groupBox_general)
        self.lineEdit_default_adduct_library.setObjectName(u"lineEdit_default_adduct_library")
        self.lineEdit_default_adduct_library.setEnabled(True)
        self.lineEdit_default_adduct_library.setMinimumSize(QSize(150, 0))
        self.lineEdit_default_adduct_library.setMaximumSize(QSize(150, 16777215))
        self.lineEdit_default_adduct_library.setReadOnly(True)

        self.gridLayout_3.addWidget(self.lineEdit_default_adduct_library, 7, 5, 1, 1)

        self.pushButton_peak_matrix = QPushButton(self.groupBox_general)
        self.pushButton_peak_matrix.setObjectName(u"pushButton_peak_matrix")
        self.pushButton_peak_matrix.setEnabled(True)
        self.pushButton_peak_matrix.setMinimumSize(QSize(95, 23))
        self.pushButton_peak_matrix.setMaximumSize(QSize(95, 16777215))

        self.gridLayout_3.addWidget(self.pushButton_peak_matrix, 7, 2, 1, 1)

        self.label_ion_mode = QLabel(self.groupBox_general)
        self.label_ion_mode.setObjectName(u"label_ion_mode")
        self.label_ion_mode.setMaximumSize(QSize(75, 16777215))

        self.gridLayout_3.addWidget(self.label_ion_mode, 7, 7, 1, 1)

        self.comboBox_ion_mode = QComboBox(self.groupBox_general)
        self.comboBox_ion_mode.addItem("")
        self.comboBox_ion_mode.addItem("")
        self.comboBox_ion_mode.setObjectName(u"comboBox_ion_mode")
        self.comboBox_ion_mode.setMinimumSize(QSize(100, 0))
        self.comboBox_ion_mode.setMaximumSize(QSize(100, 16777215))

        self.gridLayout_3.addWidget(self.comboBox_ion_mode, 7, 8, 1, 1)

        self.label_sql_database = QLabel(self.groupBox_general)
        self.label_sql_database.setObjectName(u"label_sql_database")
        self.label_sql_database.setEnabled(True)

        self.gridLayout_3.addWidget(self.label_sql_database, 4, 3, 1, 1)

        self.label_peaklist = QLabel(self.groupBox_general)
        self.label_peaklist.setObjectName(u"label_peaklist")
        self.label_peaklist.setEnabled(True)

        self.gridLayout_3.addWidget(self.label_peaklist, 5, 0, 1, 1)

        self.label_wd = QLabel(self.groupBox_general)
        self.label_wd.setObjectName(u"label_wd")
        self.label_wd.setEnabled(True)

        self.gridLayout_3.addWidget(self.label_wd, 4, 0, 1, 1)

        self.label_data_files = QLabel(self.groupBox_general)
        self.label_data_files.setObjectName(u"label_data_files")
        self.label_data_files.setEnabled(True)
        font = QFont()
        font.setFamily(u".SF NS Text")
        font.setPointSize(13)
        font.setBold(True)
        font.setWeight(75)
        self.label_data_files.setFont(font)

        self.gridLayout_3.addWidget(self.label_data_files, 3, 0, 1, 2)


        self.verticalLayout_2.addWidget(self.groupBox_general)

        self.groupBox_group_features = QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_group_features.setObjectName(u"groupBox_group_features")
        self.gridLayout_4 = QGridLayout(self.groupBox_group_features)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.gridLayout_4.setHorizontalSpacing(6)
        self.gridLayout_4.setVerticalSpacing(5)
        self.gridLayout_4.setContentsMargins(10, 5, 10, 5)
        self.doubleSpinBox_max_rt = QDoubleSpinBox(self.groupBox_group_features)
        self.doubleSpinBox_max_rt.setObjectName(u"doubleSpinBox_max_rt")
        self.doubleSpinBox_max_rt.setMaximum(9999999.000000000000000)
        self.doubleSpinBox_max_rt.setValue(5.000000000000000)

        self.gridLayout_4.addWidget(self.doubleSpinBox_max_rt, 3, 1, 1, 1)

        self.label_grouping_ncpus = QLabel(self.groupBox_group_features)
        self.label_grouping_ncpus.setObjectName(u"label_grouping_ncpus")
        self.label_grouping_ncpus.setEnabled(True)

        self.gridLayout_4.addWidget(self.label_grouping_ncpus, 3, 4, 1, 1)

        self.label_grouping_block = QLabel(self.groupBox_group_features)
        self.label_grouping_block.setObjectName(u"label_grouping_block")
        self.label_grouping_block.setEnabled(True)

        self.gridLayout_4.addWidget(self.label_grouping_block, 10, 4, 1, 1)

        self.label_max_rt = QLabel(self.groupBox_group_features)
        self.label_max_rt.setObjectName(u"label_max_rt")
        self.label_max_rt.setEnabled(True)

        self.gridLayout_4.addWidget(self.label_max_rt, 3, 0, 1, 1)

        self.label_tool_coefficient = QLabel(self.groupBox_group_features)
        self.label_tool_coefficient.setObjectName(u"label_tool_coefficient")
        self.label_tool_coefficient.setEnabled(True)

        self.gridLayout_4.addWidget(self.label_tool_coefficient, 10, 0, 1, 1)

        self.doubleSpinBox_ncpus = QDoubleSpinBox(self.groupBox_group_features)
        self.doubleSpinBox_ncpus.setObjectName(u"doubleSpinBox_ncpus")
        self.doubleSpinBox_ncpus.setDecimals(0)
        self.doubleSpinBox_ncpus.setMaximum(10000.000000000000000)
        self.doubleSpinBox_ncpus.setSingleStep(1.000000000000000)
        self.doubleSpinBox_ncpus.setValue(1.000000000000000)

        self.gridLayout_4.addWidget(self.doubleSpinBox_ncpus, 3, 5, 1, 1)

        self.doubleSpinBox_coefficent = QDoubleSpinBox(self.groupBox_group_features)
        self.doubleSpinBox_coefficent.setObjectName(u"doubleSpinBox_coefficent")
        self.doubleSpinBox_coefficent.setMaximum(1.000000000000000)
        self.doubleSpinBox_coefficent.setSingleStep(0.100000000000000)
        self.doubleSpinBox_coefficent.setValue(0.700000000000000)

        self.gridLayout_4.addWidget(self.doubleSpinBox_coefficent, 10, 1, 1, 1)

        self.checkBox_group_features = QCheckBox(self.groupBox_group_features)
        self.checkBox_group_features.setObjectName(u"checkBox_group_features")
        self.checkBox_group_features.setEnabled(True)
        self.checkBox_group_features.setFont(font)
        self.checkBox_group_features.setChecked(True)

        self.gridLayout_4.addWidget(self.checkBox_group_features, 0, 0, 1, 2)

        self.doubleSpinBox_p_value = QDoubleSpinBox(self.groupBox_group_features)
        self.doubleSpinBox_p_value.setObjectName(u"doubleSpinBox_p_value")
        self.doubleSpinBox_p_value.setDecimals(10)
        self.doubleSpinBox_p_value.setMaximum(1.000000000000000)
        self.doubleSpinBox_p_value.setSingleStep(0.010000000000000)
        self.doubleSpinBox_p_value.setValue(0.010000000000000)

        self.gridLayout_4.addWidget(self.doubleSpinBox_p_value, 10, 3, 1, 1)

        self.doubleSpinBox_block = QDoubleSpinBox(self.groupBox_group_features)
        self.doubleSpinBox_block.setObjectName(u"doubleSpinBox_block")
        self.doubleSpinBox_block.setDecimals(0)
        self.doubleSpinBox_block.setMaximum(100000000.000000000000000)
        self.doubleSpinBox_block.setSingleStep(1000.000000000000000)
        self.doubleSpinBox_block.setValue(5000.000000000000000)

        self.gridLayout_4.addWidget(self.doubleSpinBox_block, 10, 5, 1, 1)

        self.label_grouping_method = QLabel(self.groupBox_group_features)
        self.label_grouping_method.setObjectName(u"label_grouping_method")
        self.label_grouping_method.setEnabled(True)
        self.label_grouping_method.setIndent(-1)

        self.gridLayout_4.addWidget(self.label_grouping_method, 3, 2, 1, 1)

        self.comboBox_grouping_method = QComboBox(self.groupBox_group_features)
        self.comboBox_grouping_method.addItem("")
        self.comboBox_grouping_method.addItem("")
        self.comboBox_grouping_method.setObjectName(u"comboBox_grouping_method")

        self.gridLayout_4.addWidget(self.comboBox_grouping_method, 3, 3, 1, 1)

        self.label_tool_p_value = QLabel(self.groupBox_group_features)
        self.label_tool_p_value.setObjectName(u"label_tool_p_value")
        self.label_tool_p_value.setEnabled(True)

        self.gridLayout_4.addWidget(self.label_tool_p_value, 10, 2, 1, 1)

        self.checkBox_positive_correlation = QCheckBox(self.groupBox_group_features)
        self.checkBox_positive_correlation.setObjectName(u"checkBox_positive_correlation")
        self.checkBox_positive_correlation.setEnabled(True)
        self.checkBox_positive_correlation.setChecked(True)

        self.gridLayout_4.addWidget(self.checkBox_positive_correlation, 11, 0, 1, 1)


        self.verticalLayout_2.addWidget(self.groupBox_group_features)

        self.groupBox_annotate_peak_patterns = QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_annotate_peak_patterns.setObjectName(u"groupBox_annotate_peak_patterns")
        self.gridLayout = QGridLayout(self.groupBox_annotate_peak_patterns)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout.setVerticalSpacing(5)
        self.gridLayout.setContentsMargins(10, 5, 10, 5)
        self.pushButton_isotopes = QPushButton(self.groupBox_annotate_peak_patterns)
        self.pushButton_isotopes.setObjectName(u"pushButton_isotopes")
        self.pushButton_isotopes.setEnabled(True)
        self.pushButton_isotopes.setMinimumSize(QSize(95, 23))
        self.pushButton_isotopes.setMaximumSize(QSize(95, 16777215))

        self.gridLayout.addWidget(self.pushButton_isotopes, 6, 4, 1, 1)

        self.label_max_monomer_units = QLabel(self.groupBox_annotate_peak_patterns)
        self.label_max_monomer_units.setObjectName(u"label_max_monomer_units")
        self.label_max_monomer_units.setEnabled(True)
        self.label_max_monomer_units.setMinimumSize(QSize(110, 0))

        self.gridLayout.addWidget(self.label_max_monomer_units, 6, 7, 1, 1)

        self.checkBox_neutral_losses = QCheckBox(self.groupBox_annotate_peak_patterns)
        self.checkBox_neutral_losses.setObjectName(u"checkBox_neutral_losses")
        self.checkBox_neutral_losses.setEnabled(True)
        self.checkBox_neutral_losses.setChecked(True)

        self.gridLayout.addWidget(self.checkBox_neutral_losses, 4, 5, 1, 1)

        self.pushButton_neutral_losses = QPushButton(self.groupBox_annotate_peak_patterns)
        self.pushButton_neutral_losses.setObjectName(u"pushButton_neutral_losses")
        self.pushButton_neutral_losses.setEnabled(True)
        self.pushButton_neutral_losses.setMinimumSize(QSize(63, 23))
        self.pushButton_neutral_losses.setMaximumSize(QSize(95, 16777215))

        self.gridLayout.addWidget(self.pushButton_neutral_losses, 6, 6, 1, 1)

        self.checkBox_isotopes = QCheckBox(self.groupBox_annotate_peak_patterns)
        self.checkBox_isotopes.setObjectName(u"checkBox_isotopes")
        self.checkBox_isotopes.setEnabled(True)
        self.checkBox_isotopes.setChecked(True)

        self.gridLayout.addWidget(self.checkBox_isotopes, 4, 3, 1, 1)

        self.lineEdit_isotopes = QLineEdit(self.groupBox_annotate_peak_patterns)
        self.lineEdit_isotopes.setObjectName(u"lineEdit_isotopes")
        self.lineEdit_isotopes.setEnabled(True)
        self.lineEdit_isotopes.setMinimumSize(QSize(120, 0))
        self.lineEdit_isotopes.setMaximumSize(QSize(120, 16777215))
        self.lineEdit_isotopes.setReadOnly(True)

        self.gridLayout.addWidget(self.lineEdit_isotopes, 6, 3, 1, 1)

        self.checkBox_annotate_peak_patterns = QCheckBox(self.groupBox_annotate_peak_patterns)
        self.checkBox_annotate_peak_patterns.setObjectName(u"checkBox_annotate_peak_patterns")
        self.checkBox_annotate_peak_patterns.setEnabled(True)
        self.checkBox_annotate_peak_patterns.setFont(font)
        self.checkBox_annotate_peak_patterns.setChecked(True)

        self.gridLayout.addWidget(self.checkBox_annotate_peak_patterns, 0, 0, 1, 2)

        self.checkBox_adduct_library = QCheckBox(self.groupBox_annotate_peak_patterns)
        self.checkBox_adduct_library.setObjectName(u"checkBox_adduct_library")
        self.checkBox_adduct_library.setEnabled(True)
        self.checkBox_adduct_library.setChecked(True)

        self.gridLayout.addWidget(self.checkBox_adduct_library, 4, 0, 1, 1)

        self.checkBox_oligomers = QCheckBox(self.groupBox_annotate_peak_patterns)
        self.checkBox_oligomers.setObjectName(u"checkBox_oligomers")
        self.checkBox_oligomers.setEnabled(True)
        self.checkBox_oligomers.setChecked(False)

        self.gridLayout.addWidget(self.checkBox_oligomers, 4, 7, 1, 1)

        self.label_pp_ppm_tolerance = QLabel(self.groupBox_annotate_peak_patterns)
        self.label_pp_ppm_tolerance.setObjectName(u"label_pp_ppm_tolerance")
        self.label_pp_ppm_tolerance.setEnabled(True)

        self.gridLayout.addWidget(self.label_pp_ppm_tolerance, 10, 0, 1, 2)

        self.lineEdit_neutral_losses = QLineEdit(self.groupBox_annotate_peak_patterns)
        self.lineEdit_neutral_losses.setObjectName(u"lineEdit_neutral_losses")
        self.lineEdit_neutral_losses.setEnabled(True)
        self.lineEdit_neutral_losses.setMinimumSize(QSize(120, 0))
        self.lineEdit_neutral_losses.setMaximumSize(QSize(120, 16777215))
        self.lineEdit_neutral_losses.setReadOnly(True)

        self.gridLayout.addWidget(self.lineEdit_neutral_losses, 6, 5, 1, 1)

        self.spinBox_max_monomer_units = QSpinBox(self.groupBox_annotate_peak_patterns)
        self.spinBox_max_monomer_units.setObjectName(u"spinBox_max_monomer_units")
        self.spinBox_max_monomer_units.setMinimum(2)
        self.spinBox_max_monomer_units.setMaximum(1000000)
        self.spinBox_max_monomer_units.setValue(2)
        self.spinBox_max_monomer_units.setDisplayIntegerBase(10)

        self.gridLayout.addWidget(self.spinBox_max_monomer_units, 6, 8, 1, 1)

        self.lineEdit_adduct_library = QLineEdit(self.groupBox_annotate_peak_patterns)
        self.lineEdit_adduct_library.setObjectName(u"lineEdit_adduct_library")
        self.lineEdit_adduct_library.setEnabled(True)
        self.lineEdit_adduct_library.setMinimumSize(QSize(120, 0))
        self.lineEdit_adduct_library.setMaximumSize(QSize(120, 16777215))
        self.lineEdit_adduct_library.setReadOnly(True)

        self.gridLayout.addWidget(self.lineEdit_adduct_library, 6, 0, 1, 1)

        self.pushButton_adduct_library = QPushButton(self.groupBox_annotate_peak_patterns)
        self.pushButton_adduct_library.setObjectName(u"pushButton_adduct_library")
        self.pushButton_adduct_library.setEnabled(True)
        self.pushButton_adduct_library.setMinimumSize(QSize(95, 23))
        self.pushButton_adduct_library.setMaximumSize(QSize(95, 16777215))

        self.gridLayout.addWidget(self.pushButton_adduct_library, 6, 1, 1, 2)

        self.doubleSpinBox_pp_ppm_error = QDoubleSpinBox(self.groupBox_annotate_peak_patterns)
        self.doubleSpinBox_pp_ppm_error.setObjectName(u"doubleSpinBox_pp_ppm_error")
        self.doubleSpinBox_pp_ppm_error.setMaximum(100000.000000000000000)
        self.doubleSpinBox_pp_ppm_error.setValue(5.000000000000000)

        self.gridLayout.addWidget(self.doubleSpinBox_pp_ppm_error, 10, 3, 1, 1)


        self.verticalLayout_2.addWidget(self.groupBox_annotate_peak_patterns)

        self.groupBox_annotate_molecular_formulae = QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_annotate_molecular_formulae.setObjectName(u"groupBox_annotate_molecular_formulae")
        self.gridLayout_7 = QGridLayout(self.groupBox_annotate_molecular_formulae)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.gridLayout_7.setVerticalSpacing(5)
        self.gridLayout_7.setContentsMargins(10, 5, 10, 5)
        self.lineEdit_filename_mf = QLineEdit(self.groupBox_annotate_molecular_formulae)
        self.lineEdit_filename_mf.setObjectName(u"lineEdit_filename_mf")
        self.lineEdit_filename_mf.setEnabled(False)
        self.lineEdit_filename_mf.setReadOnly(True)

        self.gridLayout_7.addWidget(self.lineEdit_filename_mf, 2, 1, 1, 1)

        self.doubleSpinBox_mf_ppm_error = QDoubleSpinBox(self.groupBox_annotate_molecular_formulae)
        self.doubleSpinBox_mf_ppm_error.setObjectName(u"doubleSpinBox_mf_ppm_error")
        self.doubleSpinBox_mf_ppm_error.setEnabled(False)
        self.doubleSpinBox_mf_ppm_error.setMaximum(100000.000000000000000)
        self.doubleSpinBox_mf_ppm_error.setValue(5.000000000000000)

        self.gridLayout_7.addWidget(self.doubleSpinBox_mf_ppm_error, 2, 4, 1, 1)

        self.pushButton_filename_mf = QPushButton(self.groupBox_annotate_molecular_formulae)
        self.pushButton_filename_mf.setObjectName(u"pushButton_filename_mf")
        self.pushButton_filename_mf.setEnabled(False)
        self.pushButton_filename_mf.setMinimumSize(QSize(63, 23))
        self.pushButton_filename_mf.setFlat(False)

        self.gridLayout_7.addWidget(self.pushButton_filename_mf, 2, 2, 1, 1)

        self.label_max_mz = QLabel(self.groupBox_annotate_molecular_formulae)
        self.label_max_mz.setObjectName(u"label_max_mz")
        self.label_max_mz.setEnabled(False)

        self.gridLayout_7.addWidget(self.label_max_mz, 1, 3, 1, 1)

        self.checkBox_mf_pp_rules = QCheckBox(self.groupBox_annotate_molecular_formulae)
        self.checkBox_mf_pp_rules.setObjectName(u"checkBox_mf_pp_rules")
        self.checkBox_mf_pp_rules.setEnabled(False)
        self.checkBox_mf_pp_rules.setChecked(True)

        self.gridLayout_7.addWidget(self.checkBox_mf_pp_rules, 1, 5, 1, 1)

        self.spinBox_max_mz = QSpinBox(self.groupBox_annotate_molecular_formulae)
        self.spinBox_max_mz.setObjectName(u"spinBox_max_mz")
        self.spinBox_max_mz.setEnabled(False)
        self.spinBox_max_mz.setMaximum(1000000)
        self.spinBox_max_mz.setValue(500)
        self.spinBox_max_mz.setDisplayIntegerBase(10)

        self.gridLayout_7.addWidget(self.spinBox_max_mz, 1, 4, 1, 1)

        self.label_mf_ppm_tolerance = QLabel(self.groupBox_annotate_molecular_formulae)
        self.label_mf_ppm_tolerance.setObjectName(u"label_mf_ppm_tolerance")
        self.label_mf_ppm_tolerance.setEnabled(False)

        self.gridLayout_7.addWidget(self.label_mf_ppm_tolerance, 2, 3, 1, 1)

        self.comboBox_source_mf = QComboBox(self.groupBox_annotate_molecular_formulae)
        self.comboBox_source_mf.addItem("")
        self.comboBox_source_mf.addItem("")
        self.comboBox_source_mf.setObjectName(u"comboBox_source_mf")
        self.comboBox_source_mf.setEnabled(False)

        self.gridLayout_7.addWidget(self.comboBox_source_mf, 1, 1, 1, 1)

        self.label_filename_mf = QLabel(self.groupBox_annotate_molecular_formulae)
        self.label_filename_mf.setObjectName(u"label_filename_mf")
        self.label_filename_mf.setEnabled(False)

        self.gridLayout_7.addWidget(self.label_filename_mf, 2, 0, 1, 1)

        self.checkBox_heuristic_rules = QCheckBox(self.groupBox_annotate_molecular_formulae)
        self.checkBox_heuristic_rules.setObjectName(u"checkBox_heuristic_rules")
        self.checkBox_heuristic_rules.setEnabled(False)
        self.checkBox_heuristic_rules.setChecked(True)

        self.gridLayout_7.addWidget(self.checkBox_heuristic_rules, 2, 5, 1, 1)

        self.label_source_mf = QLabel(self.groupBox_annotate_molecular_formulae)
        self.label_source_mf.setObjectName(u"label_source_mf")
        self.label_source_mf.setEnabled(False)

        self.gridLayout_7.addWidget(self.label_source_mf, 1, 0, 1, 1)

        self.checkBox_annotate_molecular_formulae = QCheckBox(self.groupBox_annotate_molecular_formulae)
        self.checkBox_annotate_molecular_formulae.setObjectName(u"checkBox_annotate_molecular_formulae")
        self.checkBox_annotate_molecular_formulae.setEnabled(True)
        self.checkBox_annotate_molecular_formulae.setFont(font)
        self.checkBox_annotate_molecular_formulae.setCheckable(True)
        self.checkBox_annotate_molecular_formulae.setChecked(False)

        self.gridLayout_7.addWidget(self.checkBox_annotate_molecular_formulae, 0, 0, 1, 2)


        self.verticalLayout_2.addWidget(self.groupBox_annotate_molecular_formulae)

        self.groupBox_annotate_compounds = QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_annotate_compounds.setObjectName(u"groupBox_annotate_compounds")
        self.gridLayout_6 = QGridLayout(self.groupBox_annotate_compounds)
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.gridLayout_6.setVerticalSpacing(5)
        self.gridLayout_6.setContentsMargins(10, 5, 10, 5)
        self.pushButton_filename_reference = QPushButton(self.groupBox_annotate_compounds)
        self.pushButton_filename_reference.setObjectName(u"pushButton_filename_reference")
        self.pushButton_filename_reference.setEnabled(False)
        self.pushButton_filename_reference.setMinimumSize(QSize(63, 23))

        self.gridLayout_6.addWidget(self.pushButton_filename_reference, 2, 2, 1, 1)

        self.lineEdit_filename_reference = QLineEdit(self.groupBox_annotate_compounds)
        self.lineEdit_filename_reference.setObjectName(u"lineEdit_filename_reference")
        self.lineEdit_filename_reference.setEnabled(False)
        self.lineEdit_filename_reference.setReadOnly(True)

        self.gridLayout_6.addWidget(self.lineEdit_filename_reference, 2, 1, 1, 1)

        self.checkBox_cpds_pp_rules = QCheckBox(self.groupBox_annotate_compounds)
        self.checkBox_cpds_pp_rules.setObjectName(u"checkBox_cpds_pp_rules")
        self.checkBox_cpds_pp_rules.setEnabled(True)
        self.checkBox_cpds_pp_rules.setChecked(True)

        self.gridLayout_6.addWidget(self.checkBox_cpds_pp_rules, 3, 3, 1, 1)

        self.checkBox_filename_reference = QCheckBox(self.groupBox_annotate_compounds)
        self.checkBox_filename_reference.setObjectName(u"checkBox_filename_reference")
        self.checkBox_filename_reference.setEnabled(True)
        self.checkBox_filename_reference.setChecked(False)

        self.gridLayout_6.addWidget(self.checkBox_filename_reference, 1, 1, 1, 1)

        self.listWidget_databases = QListWidget(self.groupBox_annotate_compounds)
        self.listWidget_databases.setObjectName(u"listWidget_databases")
        self.listWidget_databases.setMinimumSize(QSize(510, 0))
        self.listWidget_databases.setMaximumSize(QSize(16777215, 95))

        self.gridLayout_6.addWidget(self.listWidget_databases, 1, 0, 4, 1)

        self.doubleSpinBox_cpds_ppm_error = QDoubleSpinBox(self.groupBox_annotate_compounds)
        self.doubleSpinBox_cpds_ppm_error.setObjectName(u"doubleSpinBox_cpds_ppm_error")
        self.doubleSpinBox_cpds_ppm_error.setMaximum(100000.000000000000000)
        self.doubleSpinBox_cpds_ppm_error.setValue(5.000000000000000)

        self.gridLayout_6.addWidget(self.doubleSpinBox_cpds_ppm_error, 3, 2, 1, 1)

        self.label_cpds_ppm_tolerance = QLabel(self.groupBox_annotate_compounds)
        self.label_cpds_ppm_tolerance.setObjectName(u"label_cpds_ppm_tolerance")
        self.label_cpds_ppm_tolerance.setEnabled(True)

        self.gridLayout_6.addWidget(self.label_cpds_ppm_tolerance, 3, 1, 1, 1)

        self.checkBox_annotate_compounds = QCheckBox(self.groupBox_annotate_compounds)
        self.checkBox_annotate_compounds.setObjectName(u"checkBox_annotate_compounds")
        self.checkBox_annotate_compounds.setEnabled(True)
        self.checkBox_annotate_compounds.setFont(font)
        self.checkBox_annotate_compounds.setChecked(True)

        self.gridLayout_6.addWidget(self.checkBox_annotate_compounds, 0, 0, 1, 1)


        self.verticalLayout_2.addWidget(self.groupBox_annotate_compounds)

        self.groupBox_create_summary = QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_create_summary.setObjectName(u"groupBox_create_summary")
        self.gridLayout_2 = QGridLayout(self.groupBox_create_summary)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.gridLayout_2.setVerticalSpacing(5)
        self.gridLayout_2.setContentsMargins(10, 5, 10, 5)
        self.pushButton_summary_filename = QPushButton(self.groupBox_create_summary)
        self.pushButton_summary_filename.setObjectName(u"pushButton_summary_filename")
        self.pushButton_summary_filename.setMinimumSize(QSize(63, 23))

        self.gridLayout_2.addWidget(self.pushButton_summary_filename, 1, 2, 1, 1)

        self.spinBox_mz_digits = QSpinBox(self.groupBox_create_summary)
        self.spinBox_mz_digits.setObjectName(u"spinBox_mz_digits")
        self.spinBox_mz_digits.setEnabled(False)
        self.spinBox_mz_digits.setMaximum(1000000)
        self.spinBox_mz_digits.setValue(5)
        self.spinBox_mz_digits.setDisplayIntegerBase(10)

        self.gridLayout_2.addWidget(self.spinBox_mz_digits, 3, 5, 1, 1)

        self.checkBox_convert_rt = QCheckBox(self.groupBox_create_summary)
        self.checkBox_convert_rt.setObjectName(u"checkBox_convert_rt")
        self.checkBox_convert_rt.setEnabled(True)
        self.checkBox_convert_rt.setChecked(False)

        self.gridLayout_2.addWidget(self.checkBox_convert_rt, 3, 6, 1, 1)

        self.comboBox_annotations_format = QComboBox(self.groupBox_create_summary)
        self.comboBox_annotations_format.addItem("")
        self.comboBox_annotations_format.addItem("")
        self.comboBox_annotations_format.addItem("")
        self.comboBox_annotations_format.setObjectName(u"comboBox_annotations_format")

        self.gridLayout_2.addWidget(self.comboBox_annotations_format, 3, 1, 1, 1)

        self.label_separator = QLabel(self.groupBox_create_summary)
        self.label_separator.setObjectName(u"label_separator")
        self.label_separator.setEnabled(True)

        self.gridLayout_2.addWidget(self.label_separator, 1, 3, 1, 1)

        self.lineEdit_summary_filename = QLineEdit(self.groupBox_create_summary)
        self.lineEdit_summary_filename.setObjectName(u"lineEdit_summary_filename")
        self.lineEdit_summary_filename.setReadOnly(True)

        self.gridLayout_2.addWidget(self.lineEdit_summary_filename, 1, 1, 1, 1)

        self.comboBox_convert_rt = QComboBox(self.groupBox_create_summary)
        self.comboBox_convert_rt.addItem("")
        self.comboBox_convert_rt.addItem("")
        self.comboBox_convert_rt.setObjectName(u"comboBox_convert_rt")
        self.comboBox_convert_rt.setEnabled(False)

        self.gridLayout_2.addWidget(self.comboBox_convert_rt, 3, 7, 1, 1)

        self.checkBox_create_summary = QCheckBox(self.groupBox_create_summary)
        self.checkBox_create_summary.setObjectName(u"checkBox_create_summary")
        self.checkBox_create_summary.setEnabled(True)
        self.checkBox_create_summary.setFont(font)
        self.checkBox_create_summary.setChecked(True)

        self.gridLayout_2.addWidget(self.checkBox_create_summary, 0, 0, 1, 2)

        self.checkBox_mz_digits = QCheckBox(self.groupBox_create_summary)
        self.checkBox_mz_digits.setObjectName(u"checkBox_mz_digits")
        self.checkBox_mz_digits.setEnabled(True)
        self.checkBox_mz_digits.setChecked(False)

        self.gridLayout_2.addWidget(self.checkBox_mz_digits, 3, 2, 1, 2)

        self.comboBox_separator = QComboBox(self.groupBox_create_summary)
        self.comboBox_separator.addItem("")
        self.comboBox_separator.addItem("")
        self.comboBox_separator.setObjectName(u"comboBox_separator")

        self.gridLayout_2.addWidget(self.comboBox_separator, 1, 5, 1, 1)

        self.label_summary_filename = QLabel(self.groupBox_create_summary)
        self.label_summary_filename.setObjectName(u"label_summary_filename")
        self.label_summary_filename.setEnabled(True)

        self.gridLayout_2.addWidget(self.label_summary_filename, 1, 0, 1, 1)

        self.label_annotations_format = QLabel(self.groupBox_create_summary)
        self.label_annotations_format.setObjectName(u"label_annotations_format")
        self.label_annotations_format.setEnabled(True)

        self.gridLayout_2.addWidget(self.label_annotations_format, 3, 0, 1, 1)


        self.verticalLayout_2.addWidget(self.groupBox_create_summary)

        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.pushButton_start = QPushButton(self.centralwidget)
        self.pushButton_start.setObjectName(u"pushButton_start")
        self.pushButton_start.setGeometry(QRect(860, 740, 100, 32))
        self.pushButton_cancel = QPushButton(self.centralwidget)
        self.pushButton_cancel.setObjectName(u"pushButton_cancel")
        self.pushButton_cancel.setGeometry(QRect(760, 740, 100, 32))
        self.pushButton_cancel.setMaximumSize(QSize(16777204, 16777215))
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)

        QMetaObject.connectSlotsByName(MainWindow)
    # setupUi

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QCoreApplication.translate("MainWindow", u"BEAMSpy - Birmingham mEtabolite Annotation for Mass Spectrometry (Python package)", None))
        self.actionExampleData.setText(QCoreApplication.translate("MainWindow", u"Add example data", None))
        self.actionAbout.setText(QCoreApplication.translate("MainWindow", u"About", None))
        self.groupBox_general.setTitle("")
        self.pushButton_peaklist.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.pushButton_graph.setText(QCoreApplication.translate("MainWindow", u"Save as...", None))
        self.pushButton_default_adduct_library.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.lineEdit_intensity_matrix.setText("")
        self.lineEdit_graph.setText("")
        self.label_graph.setText(QCoreApplication.translate("MainWindow", u"Graph:", None))
        self.lineEdit_wd.setText("")
        self.pushButton_sql_database.setText(QCoreApplication.translate("MainWindow", u"Save as...", None))
        self.lineEdit_sql_database.setText("")
        self.lineEdit_peaklist.setText("")
        self.pushButton_wd.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.label_default_adduct_library.setText(QCoreApplication.translate("MainWindow", u"Adduct library:", None))
        self.label_intensity_matrix.setText(QCoreApplication.translate("MainWindow", u"Intensity matrix:", None))
        self.lineEdit_default_adduct_library.setText(QCoreApplication.translate("MainWindow", u"Use default", None))
        self.pushButton_peak_matrix.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.label_ion_mode.setText(QCoreApplication.translate("MainWindow", u"Ion mode:", None))
        self.comboBox_ion_mode.setItemText(0, QCoreApplication.translate("MainWindow", u"Positive", None))
        self.comboBox_ion_mode.setItemText(1, QCoreApplication.translate("MainWindow", u"Negative", None))

        self.label_sql_database.setText(QCoreApplication.translate("MainWindow", u"Database:", None))
        self.label_peaklist.setText(QCoreApplication.translate("MainWindow", u"Peaklist:", None))
        self.label_wd.setText(QCoreApplication.translate("MainWindow", u"Working directory:", None))
        self.label_data_files.setText(QCoreApplication.translate("MainWindow", u"Data Files & General Settings", None))
        self.groupBox_group_features.setTitle("")
        self.label_grouping_ncpus.setText(QCoreApplication.translate("MainWindow", u"cpus:", None))
        self.label_grouping_block.setText(QCoreApplication.translate("MainWindow", u"Block size:", None))
        self.label_max_rt.setText(QCoreApplication.translate("MainWindow", u"Maximum RT difference (sec):", None))
        self.label_tool_coefficient.setText(QCoreApplication.translate("MainWindow", u"Coefficent threshold:", None))
        self.checkBox_group_features.setText(QCoreApplication.translate("MainWindow", u"Group Features", None))
        self.label_grouping_method.setText(QCoreApplication.translate("MainWindow", u"Grouping method:", None))
        self.comboBox_grouping_method.setItemText(0, QCoreApplication.translate("MainWindow", u"Pearson correlation", None))
        self.comboBox_grouping_method.setItemText(1, QCoreApplication.translate("MainWindow", u"Spearman-rank correlation", None))

        self.label_tool_p_value.setText(QCoreApplication.translate("MainWindow", u"P-value threshold:", None))
        self.checkBox_positive_correlation.setText(QCoreApplication.translate("MainWindow", u"Positive Correlation", None))
        self.groupBox_annotate_peak_patterns.setTitle("")
        self.pushButton_isotopes.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.label_max_monomer_units.setText(QCoreApplication.translate("MainWindow", u"Monomer Units:", None))
        self.checkBox_neutral_losses.setText(QCoreApplication.translate("MainWindow", u"Neutral losses", None))
        self.pushButton_neutral_losses.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.checkBox_isotopes.setText(QCoreApplication.translate("MainWindow", u"Isotopes", None))
        self.lineEdit_isotopes.setText(QCoreApplication.translate("MainWindow", u"Use default", None))
        self.checkBox_annotate_peak_patterns.setText(QCoreApplication.translate("MainWindow", u"Annotate Peak Patterns", None))
        self.checkBox_adduct_library.setText(QCoreApplication.translate("MainWindow", u"Adducts", None))
        self.checkBox_oligomers.setText(QCoreApplication.translate("MainWindow", u"Oligomers", None))
        self.label_pp_ppm_tolerance.setText(QCoreApplication.translate("MainWindow", u"Mass tolerance (ppm):", None))
        self.lineEdit_neutral_losses.setText(QCoreApplication.translate("MainWindow", u"Use default", None))
        self.lineEdit_adduct_library.setText(QCoreApplication.translate("MainWindow", u"Use default", None))
        self.pushButton_adduct_library.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.groupBox_annotate_molecular_formulae.setTitle("")
        self.lineEdit_filename_mf.setText("")
        self.pushButton_filename_mf.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.label_max_mz.setText(QCoreApplication.translate("MainWindow", u"Maximum m/z:", None))
        self.checkBox_mf_pp_rules.setText(QCoreApplication.translate("MainWindow", u"use peak patterns", None))
        self.label_mf_ppm_tolerance.setText(QCoreApplication.translate("MainWindow", u"Mass tolerance (ppm):", None))
        self.comboBox_source_mf.setItemText(0, QCoreApplication.translate("MainWindow", u"https://mfdb.bham.ac.uk", None))
        self.comboBox_source_mf.setItemText(1, QCoreApplication.translate("MainWindow", u"Tab-delimited text file", None))

        self.label_filename_mf.setText(QCoreApplication.translate("MainWindow", u"Reference file:", None))
        self.checkBox_heuristic_rules.setText(QCoreApplication.translate("MainWindow", u"Heuristic rules", None))
        self.label_source_mf.setText(QCoreApplication.translate("MainWindow", u"Source:", None))
        self.checkBox_annotate_molecular_formulae.setText(QCoreApplication.translate("MainWindow", u"Annotate Molecular Formulae", None))
        self.groupBox_annotate_compounds.setTitle("")
        self.pushButton_filename_reference.setText(QCoreApplication.translate("MainWindow", u"Browse...", None))
        self.lineEdit_filename_reference.setText("")
        self.checkBox_cpds_pp_rules.setText(QCoreApplication.translate("MainWindow", u"use peak patterns", None))
        self.checkBox_filename_reference.setText(QCoreApplication.translate("MainWindow", u"Reference file", None))
        self.label_cpds_ppm_tolerance.setText(QCoreApplication.translate("MainWindow", u"Mass tolerance (ppm):", None))
        self.checkBox_annotate_compounds.setText(QCoreApplication.translate("MainWindow", u"Annotate Compounds / Metabolites", None))
        self.groupBox_create_summary.setTitle("")
        self.pushButton_summary_filename.setText(QCoreApplication.translate("MainWindow", u"Save as...", None))
        self.checkBox_convert_rt.setText(QCoreApplication.translate("MainWindow", u"Convert RT:", None))
        self.comboBox_annotations_format.setItemText(0, QCoreApplication.translate("MainWindow", u"Multiple rows for each feature", None))
        self.comboBox_annotations_format.setItemText(1, QCoreApplication.translate("MainWindow", u"Single row for each feature and separate columns", None))
        self.comboBox_annotations_format.setItemText(2, QCoreApplication.translate("MainWindow", u"Single row for each feature and merged columns", None))

        self.label_separator.setText(QCoreApplication.translate("MainWindow", u"Separator:", None))
        self.lineEdit_summary_filename.setText("")
        self.comboBox_convert_rt.setItemText(0, QCoreApplication.translate("MainWindow", u"Minutes", None))
        self.comboBox_convert_rt.setItemText(1, QCoreApplication.translate("MainWindow", u"Seconds", None))

        self.checkBox_create_summary.setText(QCoreApplication.translate("MainWindow", u"Create summary", None))
        self.checkBox_mz_digits.setText(QCoreApplication.translate("MainWindow", u"Number of digits m/z:", None))
        self.comboBox_separator.setItemText(0, QCoreApplication.translate("MainWindow", u"tab", None))
        self.comboBox_separator.setItemText(1, QCoreApplication.translate("MainWindow", u"comma", None))

        self.label_summary_filename.setText(QCoreApplication.translate("MainWindow", u"Summary:", None))
        self.label_annotations_format.setText(QCoreApplication.translate("MainWindow", u"Annotations:", None))
        self.pushButton_start.setText(QCoreApplication.translate("MainWindow", u"Start", None))
        self.pushButton_cancel.setText(QCoreApplication.translate("MainWindow", u"Cancel", None))
    # retranslateUi

