# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgCalculateLCSetAggregation.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculateLCSetAggregation(object):
    def setupUi(self, DlgCalculateLCSetAggregation):
        DlgCalculateLCSetAggregation.setObjectName("DlgCalculateLCSetAggregation")
        DlgCalculateLCSetAggregation.setWindowModality(QtCore.Qt.WindowModal)
        DlgCalculateLCSetAggregation.resize(1065, 501)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgCalculateLCSetAggregation.sizePolicy().hasHeightForWidth())
        DlgCalculateLCSetAggregation.setSizePolicy(sizePolicy)
        DlgCalculateLCSetAggregation.setMinimumSize(QtCore.QSize(0, 0))
        self.gridLayout = QtWidgets.QGridLayout(DlgCalculateLCSetAggregation)
        self.gridLayout.setObjectName("gridLayout")
        self.btn_save = QtWidgets.QPushButton(DlgCalculateLCSetAggregation)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_save.sizePolicy().hasHeightForWidth())
        self.btn_save.setSizePolicy(sizePolicy)
        self.btn_save.setMinimumSize(QtCore.QSize(0, 40))
        self.btn_save.setObjectName("btn_save")
        self.gridLayout.addWidget(self.btn_save, 3, 1, 1, 1)
        self.btn_load = QtWidgets.QPushButton(DlgCalculateLCSetAggregation)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_load.sizePolicy().hasHeightForWidth())
        self.btn_load.setSizePolicy(sizePolicy)
        self.btn_load.setMinimumSize(QtCore.QSize(0, 40))
        self.btn_load.setObjectName("btn_load")
        self.gridLayout.addWidget(self.btn_load, 3, 0, 1, 1)
        self.remap_view = QtWidgets.QTableView(DlgCalculateLCSetAggregation)
        self.remap_view.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectItems)
        self.remap_view.setSortingEnabled(True)
        self.remap_view.setObjectName("remap_view")
        self.remap_view.horizontalHeader().setVisible(True)
        self.remap_view.horizontalHeader().setCascadingSectionResizes(True)
        self.remap_view.horizontalHeader().setStretchLastSection(True)
        self.remap_view.verticalHeader().setVisible(False)
        self.remap_view.verticalHeader().setSortIndicatorShown(True)
        self.gridLayout.addWidget(self.remap_view, 0, 0, 1, 2)
        self.btn_reset = QtWidgets.QPushButton(DlgCalculateLCSetAggregation)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_reset.sizePolicy().hasHeightForWidth())
        self.btn_reset.setSizePolicy(sizePolicy)
        self.btn_reset.setMinimumSize(QtCore.QSize(0, 40))
        self.btn_reset.setObjectName("btn_reset")
        self.gridLayout.addWidget(self.btn_reset, 1, 0, 1, 2)
        self.btn_close = QtWidgets.QPushButton(DlgCalculateLCSetAggregation)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_close.sizePolicy().hasHeightForWidth())
        self.btn_close.setSizePolicy(sizePolicy)
        self.btn_close.setMinimumSize(QtCore.QSize(0, 40))
        self.btn_close.setObjectName("btn_close")
        self.gridLayout.addWidget(self.btn_close, 4, 0, 1, 2)

        self.retranslateUi(DlgCalculateLCSetAggregation)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculateLCSetAggregation)

    def retranslateUi(self, DlgCalculateLCSetAggregation):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculateLCSetAggregation.setWindowTitle(_translate("DlgCalculateLCSetAggregation", "Setup aggregation of data"))
        self.btn_save.setText(_translate("DlgCalculateLCSetAggregation", "Save definition to file"))
        self.btn_load.setText(_translate("DlgCalculateLCSetAggregation", "Load definition from file"))
        self.btn_reset.setText(_translate("DlgCalculateLCSetAggregation", "Reset to default"))
        self.btn_close.setText(_translate("DlgCalculateLCSetAggregation", "Save"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculateLCSetAggregation = QtWidgets.QDialog()
    ui = Ui_DlgCalculateLCSetAggregation()
    ui.setupUi(DlgCalculateLCSetAggregation)
    DlgCalculateLCSetAggregation.show()
    sys.exit(app.exec_())
