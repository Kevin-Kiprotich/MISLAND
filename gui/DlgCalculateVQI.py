# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgCalculateVQI.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculateVQI(object):
    def setupUi(self, DlgCalculateVQI):
        DlgCalculateVQI.setObjectName("DlgCalculateVQI")
        DlgCalculateVQI.resize(350, 160)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgCalculateVQI.sizePolicy().hasHeightForWidth())
        DlgCalculateVQI.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgCalculateVQI)
        self.verticalLayout.setObjectName("verticalLayout")
        self.TabBox = QtWidgets.QTabWidget(DlgCalculateVQI)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.TabBox.sizePolicy().hasHeightForWidth())
        self.TabBox.setSizePolicy(sizePolicy)
        self.TabBox.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.TabBox.setObjectName("TabBox")
        self.verticalLayout.addWidget(self.TabBox)
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.button_next = QtWidgets.QPushButton(DlgCalculateVQI)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_next.sizePolicy().hasHeightForWidth())
        self.button_next.setSizePolicy(sizePolicy)
        self.button_next.setMinimumSize(QtCore.QSize(0, 40))
        self.button_next.setObjectName("button_next")
        self.gridLayout_3.addWidget(self.button_next, 0, 1, 1, 1)
        self.button_calculate = QtWidgets.QPushButton(DlgCalculateVQI)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_calculate.sizePolicy().hasHeightForWidth())
        self.button_calculate.setSizePolicy(sizePolicy)
        self.button_calculate.setMinimumSize(QtCore.QSize(0, 40))
        self.button_calculate.setObjectName("button_calculate")
        self.gridLayout_3.addWidget(self.button_calculate, 1, 0, 1, 2)
        self.button_prev = QtWidgets.QPushButton(DlgCalculateVQI)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_prev.sizePolicy().hasHeightForWidth())
        self.button_prev.setSizePolicy(sizePolicy)
        self.button_prev.setMinimumSize(QtCore.QSize(0, 40))
        self.button_prev.setObjectName("button_prev")
        self.gridLayout_3.addWidget(self.button_prev, 0, 0, 1, 1)
        self.verticalLayout.addLayout(self.gridLayout_3)

        self.retranslateUi(DlgCalculateVQI)
        self.TabBox.setCurrentIndex(-1)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculateVQI)
        DlgCalculateVQI.setTabOrder(self.TabBox, self.button_prev)
        DlgCalculateVQI.setTabOrder(self.button_prev, self.button_next)
        DlgCalculateVQI.setTabOrder(self.button_next, self.button_calculate)

    def retranslateUi(self, DlgCalculateVQI):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculateVQI.setWindowTitle(_translate("DlgCalculateVQI", "Calculate Vegetation Quality Index"))
        self.button_next.setText(_translate("DlgCalculateVQI", "Next"))
        self.button_calculate.setText(_translate("DlgCalculateVQI", "Calculate"))
        self.button_prev.setText(_translate("DlgCalculateVQI", "Previous"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculateVQI = QtWidgets.QDialog()
    ui = Ui_DlgCalculateVQI()
    ui.setupUi(DlgCalculateVQI)
    DlgCalculateVQI.show()
    sys.exit(app.exec_())