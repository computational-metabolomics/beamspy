# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'databases.ui'
#
# Created by: PyQt5 UI code generator 5.6
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        Dialog.setObjectName("Dialog")
        Dialog.resize(531, 323)
        self.listWidget_datbase_names = QtWidgets.QListWidget(Dialog)
        self.listWidget_datbase_names.setGeometry(QtCore.QRect(10, 10, 250, 301))
        self.listWidget_datbase_names.setObjectName("listWidget_datbase_names")
        item = QtWidgets.QListWidgetItem()
        self.listWidget_datbase_names.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listWidget_datbase_names.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listWidget_datbase_names.addItem(item)
        self.listWidget_categories = QtWidgets.QListWidget(Dialog)
        self.listWidget_categories.setGeometry(QtCore.QRect(270, 10, 250, 301))
        self.listWidget_categories.setObjectName("listWidget_categories")
        item = QtWidgets.QListWidgetItem()
        self.listWidget_categories.addItem(item)
        item = QtWidgets.QListWidgetItem()
        self.listWidget_categories.addItem(item)

        self.retranslateUi(Dialog)
        QtCore.QMetaObject.connectSlotsByName(Dialog)

    def retranslateUi(self, Dialog):
        _translate = QtCore.QCoreApplication.translate
        Dialog.setWindowTitle(_translate("Dialog", "Dialog"))
        __sortingEnabled = self.listWidget_datbase_names.isSortingEnabled()
        self.listWidget_datbase_names.setSortingEnabled(False)
        item = self.listWidget_datbase_names.item(0)
        item.setText(_translate("Dialog", "HMDB"))
        item = self.listWidget_datbase_names.item(1)
        item.setText(_translate("Dialog", "LIPIDMAPS"))
        item = self.listWidget_datbase_names.item(2)
        item.setText(_translate("Dialog", "BIOCYC"))
        self.listWidget_datbase_names.setSortingEnabled(__sortingEnabled)
        __sortingEnabled = self.listWidget_categories.isSortingEnabled()
        self.listWidget_categories.setSortingEnabled(False)
        item = self.listWidget_categories.item(0)
        item.setText(_translate("Dialog", "CAT1"))
        item = self.listWidget_categories.item(1)
        item.setText(_translate("Dialog", "CAT2"))
        self.listWidget_categories.setSortingEnabled(__sortingEnabled)

