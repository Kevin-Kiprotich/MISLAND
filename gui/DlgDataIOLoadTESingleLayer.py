# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgDataIOLoadTESingleLayer.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgDataIOLoadTESingleLayer(object):
    def setupUi(self, DlgDataIOLoadTESingleLayer):
        DlgDataIOLoadTESingleLayer.setObjectName("DlgDataIOLoadTESingleLayer")
        DlgDataIOLoadTESingleLayer.setWindowModality(QtCore.Qt.ApplicationModal)
        DlgDataIOLoadTESingleLayer.resize(294, 272)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgDataIOLoadTESingleLayer.sizePolicy().hasHeightForWidth())
        DlgDataIOLoadTESingleLayer.setSizePolicy(sizePolicy)
        self.gridLayout = QtWidgets.QGridLayout(DlgDataIOLoadTESingleLayer)
        self.gridLayout.setObjectName("gridLayout")
        self.buttonBox = QtWidgets.QDialogButtonBox(DlgDataIOLoadTESingleLayer)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 3, 0, 1, 1)
        self.groupBox = QtWidgets.QGroupBox(DlgDataIOLoadTESingleLayer)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.layers_view = QtWidgets.QListView(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.layers_view.sizePolicy().hasHeightForWidth())
        self.layers_view.setSizePolicy(sizePolicy)
        self.layers_view.setMinimumSize(QtCore.QSize(200, 150))
        self.layers_view.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
        self.layers_view.setObjectName("layers_view")
        self.verticalLayout_2.addWidget(self.layers_view)
        self.gridLayout.addWidget(self.groupBox, 0, 0, 1, 1)

        self.retranslateUi(DlgDataIOLoadTESingleLayer)
        QtCore.QMetaObject.connectSlotsByName(DlgDataIOLoadTESingleLayer)

    def retranslateUi(self, DlgDataIOLoadTESingleLayer):
        _translate = QtCore.QCoreApplication.translate
        DlgDataIOLoadTESingleLayer.setWindowTitle(_translate("DlgDataIOLoadTESingleLayer", "Open a MISLAND file"))
        self.groupBox.setTitle(_translate("DlgDataIOLoadTESingleLayer", "Select a layer"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgDataIOLoadTESingleLayer = QtWidgets.QDialog()
    ui = Ui_DlgDataIOLoadTESingleLayer()
    ui.setupUi(DlgDataIOLoadTESingleLayer)
    DlgDataIOLoadTESingleLayer.show()
    sys.exit(app.exec_())
