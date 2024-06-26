# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgImportDataProd.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgImportProd(object):
    def setupUi(self, DlgImportProd):
        DlgImportProd.setObjectName("DlgImportProd")
        DlgImportProd.setWindowModality(QtCore.Qt.ApplicationModal)
        DlgImportProd.resize(238, 230)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgImportProd.sizePolicy().hasHeightForWidth())
        DlgImportProd.setSizePolicy(sizePolicy)
        DlgImportProd.setModal(True)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgImportProd)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox_lc_agg = QtWidgets.QGroupBox(DlgImportProd)
        self.groupBox_lc_agg.setEnabled(True)
        self.groupBox_lc_agg.setObjectName("groupBox_lc_agg")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_lc_agg)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.textBrowser = QtWidgets.QTextBrowser(self.groupBox_lc_agg)
        self.textBrowser.setMinimumSize(QtCore.QSize(200, 150))
        self.textBrowser.setAutoFillBackground(False)
        self.textBrowser.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.textBrowser.setObjectName("textBrowser")
        self.gridLayout_2.addWidget(self.textBrowser, 0, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_lc_agg)
        self.btnBox = QtWidgets.QDialogButtonBox(DlgImportProd)
        self.btnBox.setOrientation(QtCore.Qt.Horizontal)
        self.btnBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.btnBox.setObjectName("btnBox")
        self.verticalLayout.addWidget(self.btnBox)

        self.retranslateUi(DlgImportProd)
        self.btnBox.accepted.connect(DlgImportProd.accept)
        self.btnBox.rejected.connect(DlgImportProd.reject)
        QtCore.QMetaObject.connectSlotsByName(DlgImportProd)

    def retranslateUi(self, DlgImportProd):
        _translate = QtCore.QCoreApplication.translate
        DlgImportProd.setWindowTitle(_translate("DlgImportProd", "Load a Custom Land Productivity Dataset"))
        self.groupBox_lc_agg.setTitle(_translate("DlgImportProd", "Productivity class definition"))
        self.textBrowser.setHtml(_translate("DlgImportProd", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:7.875pt; font-weight:400; font-style:normal;\">\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt; font-weight:600;\">Productivity classes in the input data must be coded as follows:</span></p>\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:8pt;\"><br /></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">1: Declining</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">2: Early signs of decline</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">3: Stable but stressed</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">4: Stable</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">5: Increasing</span></p>\n"
"<p style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-size:8pt;\">0 or -32768: No data</span></p></body></html>"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgImportProd = QtWidgets.QDialog()
    ui = Ui_DlgImportProd()
    ui.setupUi(DlgImportProd)
    DlgImportProd.show()
    sys.exit(app.exec_())
