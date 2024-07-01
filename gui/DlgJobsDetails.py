# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgJobsDetails.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgJobsDetails(object):
    def setupUi(self, DlgJobsDetails):
        DlgJobsDetails.setObjectName("DlgJobsDetails")
        DlgJobsDetails.setWindowModality(QtCore.Qt.WindowModal)
        DlgJobsDetails.resize(418, 359)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgJobsDetails.sizePolicy().hasHeightForWidth())
        DlgJobsDetails.setSizePolicy(sizePolicy)
        self.gridLayout = QtWidgets.QGridLayout(DlgJobsDetails)
        self.gridLayout.setObjectName("gridLayout")
        self.comments = QtWidgets.QTextBrowser(DlgJobsDetails)
        self.comments.setMinimumSize(QtCore.QSize(400, 0))
        self.comments.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.comments.setOpenExternalLinks(True)
        self.comments.setObjectName("comments")
        self.gridLayout.addWidget(self.comments, 3, 0, 1, 2)
        self.commentsLabel = QtWidgets.QLabel(DlgJobsDetails)
        self.commentsLabel.setObjectName("commentsLabel")
        self.gridLayout.addWidget(self.commentsLabel, 2, 0, 1, 1)
        self.taskNameLabel = QtWidgets.QLabel(DlgJobsDetails)
        self.taskNameLabel.setObjectName("taskNameLabel")
        self.gridLayout.addWidget(self.taskNameLabel, 0, 0, 1, 1)
        self.outputLabel = QtWidgets.QLabel(DlgJobsDetails)
        self.outputLabel.setObjectName("outputLabel")
        self.gridLayout.addWidget(self.outputLabel, 6, 0, 1, 1)
        self.inputLabel = QtWidgets.QLabel(DlgJobsDetails)
        self.inputLabel.setObjectName("inputLabel")
        self.gridLayout.addWidget(self.inputLabel, 4, 0, 1, 1)
        self.statusLabel = QtWidgets.QLabel(DlgJobsDetails)
        self.statusLabel.setObjectName("statusLabel")
        self.gridLayout.addWidget(self.statusLabel, 0, 1, 1, 1)
        self.output = QtWidgets.QTextBrowser(DlgJobsDetails)
        self.output.setMinimumSize(QtCore.QSize(400, 0))
        self.output.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.output.setObjectName("output")
        self.gridLayout.addWidget(self.output, 7, 0, 1, 2)
        self.input = QtWidgets.QTextBrowser(DlgJobsDetails)
        self.input.setMinimumSize(QtCore.QSize(400, 0))
        self.input.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.input.setObjectName("input")
        self.gridLayout.addWidget(self.input, 5, 0, 1, 2)
        self.task_name = QtWidgets.QTextBrowser(DlgJobsDetails)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.task_name.sizePolicy().hasHeightForWidth())
        self.task_name.setSizePolicy(sizePolicy)
        self.task_name.setMinimumSize(QtCore.QSize(0, 30))
        self.task_name.setMaximumSize(QtCore.QSize(16777215, 40))
        self.task_name.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.task_name.setObjectName("task_name")
        self.gridLayout.addWidget(self.task_name, 1, 0, 1, 1)
        self.task_status = QtWidgets.QTextBrowser(DlgJobsDetails)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.task_status.sizePolicy().hasHeightForWidth())
        self.task_status.setSizePolicy(sizePolicy)
        self.task_status.setMinimumSize(QtCore.QSize(0, 30))
        self.task_status.setMaximumSize(QtCore.QSize(16777215, 40))
        self.task_status.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.task_status.setObjectName("task_status")
        self.gridLayout.addWidget(self.task_status, 1, 1, 1, 1)

        self.retranslateUi(DlgJobsDetails)
        QtCore.QMetaObject.connectSlotsByName(DlgJobsDetails)

    def retranslateUi(self, DlgJobsDetails):
        _translate = QtCore.QCoreApplication.translate
        DlgJobsDetails.setWindowTitle(_translate("DlgJobsDetails", "Google Earth Engine Task Detail"))
        self.comments.setHtml(_translate("DlgJobsDetails", "<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"</style></head><body style=\" font-family:\'MS Shell Dlg 2\'; font-size:7.875pt; font-weight:400; font-style:normal;\">\n"
"<p style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-size:8pt;\"><br /></p></body></html>"))
        self.commentsLabel.setText(_translate("DlgJobsDetails", "Comments:"))
        self.taskNameLabel.setText(_translate("DlgJobsDetails", "Task name:"))
        self.outputLabel.setText(_translate("DlgJobsDetails", "Output:"))
        self.inputLabel.setText(_translate("DlgJobsDetails", "Input:"))
        self.statusLabel.setText(_translate("DlgJobsDetails", "Status:"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgJobsDetails = QtWidgets.QDialog()
    ui = Ui_DlgJobsDetails()
    ui.setupUi(DlgJobsDetails)
    DlgJobsDetails.show()
    sys.exit(app.exec_())