# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgImportDataLC.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgImportLC(object):
    def setupUi(self, DlgImportLC):
        DlgImportLC.setObjectName("DlgImportLC")
        DlgImportLC.setWindowModality(QtCore.Qt.ApplicationModal)
        DlgImportLC.resize(434, 196)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgImportLC.sizePolicy().hasHeightForWidth())
        DlgImportLC.setSizePolicy(sizePolicy)
        DlgImportLC.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgImportLC)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox_lc_agg = QtWidgets.QGroupBox(DlgImportLC)
        self.groupBox_lc_agg.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_lc_agg.sizePolicy().hasHeightForWidth())
        self.groupBox_lc_agg.setSizePolicy(sizePolicy)
        self.groupBox_lc_agg.setObjectName("groupBox_lc_agg")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_lc_agg)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.btn_agg_edit_def = QtWidgets.QPushButton(self.groupBox_lc_agg)
        self.btn_agg_edit_def.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_agg_edit_def.sizePolicy().hasHeightForWidth())
        self.btn_agg_edit_def.setSizePolicy(sizePolicy)
        self.btn_agg_edit_def.setMinimumSize(QtCore.QSize(0, 30))
        self.btn_agg_edit_def.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.btn_agg_edit_def.setObjectName("btn_agg_edit_def")
        self.verticalLayout_2.addWidget(self.btn_agg_edit_def)
        self.checkBox_use_sample = QtWidgets.QCheckBox(self.groupBox_lc_agg)
        self.checkBox_use_sample.setChecked(True)
        self.checkBox_use_sample.setObjectName("checkBox_use_sample")
        self.verticalLayout_2.addWidget(self.checkBox_use_sample)
        self.label = QtWidgets.QLabel(self.groupBox_lc_agg)
        self.label.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(0, 0))
        self.label.setBaseSize(QtCore.QSize(0, 0))
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.verticalLayout_2.addWidget(self.label)
        self.verticalLayout.addWidget(self.groupBox_lc_agg)
        self.btnBox = QtWidgets.QDialogButtonBox(DlgImportLC)
        self.btnBox.setOrientation(QtCore.Qt.Horizontal)
        self.btnBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.btnBox.setObjectName("btnBox")
        self.verticalLayout.addWidget(self.btnBox)

        self.retranslateUi(DlgImportLC)
        self.btnBox.accepted.connect(DlgImportLC.accept)
        self.btnBox.rejected.connect(DlgImportLC.reject)
        QtCore.QMetaObject.connectSlotsByName(DlgImportLC)

    def retranslateUi(self, DlgImportLC):
        _translate = QtCore.QCoreApplication.translate
        DlgImportLC.setWindowTitle(_translate("DlgImportLC", "Load a Custom Land Cover Dataset"))
        self.groupBox_lc_agg.setTitle(_translate("DlgImportLC", "Choose a land cover aggregation method"))
        self.btn_agg_edit_def.setText(_translate("DlgImportLC", "Edit definition"))
        self.checkBox_use_sample.setText(_translate("DlgImportLC", "Use sample when reading cover classes from input file"))
        self.label.setText(_translate("DlgImportLC", "Note: If reading a large file it is recommended that the above option be checked, as it will singificantly speed the process of reading the input classes from the dataset. However, if you find that Trends.Earth is not identifying all of the classes in the input file, it may be necessary to turn off this option."))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgImportLC = QtWidgets.QDialog()
    ui = Ui_DlgImportLC()
    ui.setupUi(DlgImportLC)
    DlgImportLC.show()
    sys.exit(app.exec_())