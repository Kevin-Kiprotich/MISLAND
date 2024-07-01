# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgImportDataSOC.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgImportSOC(object):
    def setupUi(self, DlgImportSOC):
        DlgImportSOC.setObjectName("DlgImportSOC")
        DlgImportSOC.setWindowModality(QtCore.Qt.ApplicationModal)
        DlgImportSOC.resize(700, 82)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgImportSOC.sizePolicy().hasHeightForWidth())
        DlgImportSOC.setSizePolicy(sizePolicy)
        DlgImportSOC.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgImportSOC)
        self.verticalLayout.setObjectName("verticalLayout")
        self.btnBox = QtWidgets.QDialogButtonBox(DlgImportSOC)
        self.btnBox.setOrientation(QtCore.Qt.Horizontal)
        self.btnBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.btnBox.setObjectName("btnBox")
        self.verticalLayout.addWidget(self.btnBox)

        self.retranslateUi(DlgImportSOC)
        self.btnBox.accepted.connect(DlgImportSOC.accept)
        self.btnBox.rejected.connect(DlgImportSOC.reject)
        QtCore.QMetaObject.connectSlotsByName(DlgImportSOC)

    def retranslateUi(self, DlgImportSOC):
        _translate = QtCore.QCoreApplication.translate
        DlgImportSOC.setWindowTitle(_translate("DlgImportSOC", "Load a Custom Soil Organic Carbon (SOC) dataset"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgImportSOC = QtWidgets.QDialog()
    ui = Ui_DlgImportSOC()
    ui.setupUi(DlgImportSOC)
    DlgImportSOC.show()
    sys.exit(app.exec_())