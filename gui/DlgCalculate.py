# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DlgCalculate.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculate(object):
    def setupUi(self, DlgCalculate):
        DlgCalculate.setObjectName("DlgCalculate")
        DlgCalculate.setWindowModality(QtCore.Qt.ApplicationModal)
        DlgCalculate.resize(570, 385)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgCalculate)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox = QtWidgets.QGroupBox(DlgCalculate)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.pushButton_ld = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_ld.sizePolicy().hasHeightForWidth())
        self.pushButton_ld.setSizePolicy(sizePolicy)
        self.pushButton_ld.setMinimumSize(QtCore.QSize(0, 48))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_ld.setFont(font)
        self.pushButton_ld.setObjectName("pushButton_ld")
        self.verticalLayout_2.addWidget(self.pushButton_ld)
        self.pushButton_timeseries = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_timeseries.sizePolicy().hasHeightForWidth())
        self.pushButton_timeseries.setSizePolicy(sizePolicy)
        self.pushButton_timeseries.setMinimumSize(QtCore.QSize(0, 48))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_timeseries.setFont(font)
        self.pushButton_timeseries.setObjectName("pushButton_timeseries")
        self.verticalLayout_2.addWidget(self.pushButton_timeseries)
        self.pushButton_forest = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_forest.sizePolicy().hasHeightForWidth())
        self.pushButton_forest.setSizePolicy(sizePolicy)
        self.pushButton_forest.setMinimumSize(QtCore.QSize(0, 48))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_forest.setFont(font)
        self.pushButton_forest.setObjectName("pushButton_forest")
        self.verticalLayout_2.addWidget(self.pushButton_forest)
        self.pushButton_medalus = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_medalus.sizePolicy().hasHeightForWidth())
        self.pushButton_medalus.setSizePolicy(sizePolicy)
        self.pushButton_medalus.setMinimumSize(QtCore.QSize(0, 48))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_medalus.setFont(font)
        self.pushButton_medalus.setObjectName("pushButton_medalus")
        self.verticalLayout_2.addWidget(self.pushButton_medalus)
        self.pushbutton_CE = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushbutton_CE.sizePolicy().hasHeightForWidth())
        self.pushbutton_CE.setSizePolicy(sizePolicy)
        self.pushbutton_CE.setMinimumSize(QtCore.QSize(0, 48))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushbutton_CE.setFont(font)
        self.pushbutton_CE.setObjectName("pushbutton_CE")
        self.verticalLayout_2.addWidget(self.pushbutton_CE)
        self.pushButton_SE = QtWidgets.QPushButton(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_SE.sizePolicy().hasHeightForWidth())
        self.pushButton_SE.setSizePolicy(sizePolicy)
        self.pushButton_SE.setMinimumSize(QtCore.QSize(0, 48))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_SE.setFont(font)
        self.pushButton_SE.setObjectName("pushButton_SE")
        self.verticalLayout_2.addWidget(self.pushButton_SE)
        self.verticalLayout.addWidget(self.groupBox)

        self.retranslateUi(DlgCalculate)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculate)

    def retranslateUi(self, DlgCalculate):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculate.setWindowTitle(_translate("DlgCalculate", "Calculate Indicators"))
        self.groupBox.setTitle(_translate("DlgCalculate", "MISLAND tools"))
        self.pushButton_ld.setText(_translate("DlgCalculate", "Land degradation indicator\n"
"(SDG indicator 15.3.1)"))
        self.pushButton_timeseries.setText(_translate("DlgCalculate", "Vegetation Indices Time-Series\n"
"(NDVI, MSAVI2 and SAVI)"))
        self.pushButton_forest.setText(_translate("DlgCalculate", "Forest Degradation Hotspots\n"
"(Forest Loss, Gain, Cover and Fires)"))
        self.pushButton_medalus.setText(_translate("DlgCalculate", "MEDALUS Indicators for vulnerability to desertification\n"
"(CQI, SQI, LMSQI and VQI) (*BETA Version)"))
        self.pushbutton_CE.setText(_translate("DlgCalculate", "Coastal Erosion (Coastal Vulnerability Index)"))
        self.pushButton_SE.setText(_translate("DlgCalculate", "Soil Erosion (ILSWE, RUSSLE)"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculate = QtWidgets.QDialog()
    ui = Ui_DlgCalculate()
    ui.setupUi(DlgCalculate)
    DlgCalculate.show()
    sys.exit(app.exec_())
