# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'DlgCalculateCE.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculateCE(object):
    def setupUi(self, DlgCalculateCE):
        DlgCalculateCE.setObjectName("DlgCalculateCE")
        DlgCalculateCE.resize(487, 428)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgCalculateCE.sizePolicy().hasHeightForWidth())
        DlgCalculateCE.setSizePolicy(sizePolicy)
        DlgCalculateCE.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(DlgCalculateCE)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox = QtWidgets.QGroupBox(DlgCalculateCE)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.groupBox)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.verticalLayout = QtWidgets.QVBoxLayout()
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName("verticalLayout")
        self.label = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.verticalLayout.addWidget(self.label)
        self.cvi_country_combobox = QtWidgets.QComboBox(self.groupBox)
        self.cvi_country_combobox.setMinimumSize(QtCore.QSize(0, 32))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cvi_country_combobox.setFont(font)
        self.cvi_country_combobox.setObjectName("cvi_country_combobox")
        self.verticalLayout.addWidget(self.cvi_country_combobox)
        self.horizontalLayout.addLayout(self.verticalLayout)
        self.verticalLayout_3 = QtWidgets.QVBoxLayout()
        self.verticalLayout_3.setSpacing(0)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.label_2 = QtWidgets.QLabel(self.groupBox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Maximum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_2.setFont(font)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_3.addWidget(self.label_2)
        self.cvi_region_combobox = QtWidgets.QComboBox(self.groupBox)
        self.cvi_region_combobox.setMinimumSize(QtCore.QSize(0, 32))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cvi_region_combobox.setFont(font)
        self.cvi_region_combobox.setObjectName("cvi_region_combobox")
        self.verticalLayout_3.addWidget(self.cvi_region_combobox)
        self.horizontalLayout.addLayout(self.verticalLayout_3)
        self.verticalLayout_2.addWidget(self.groupBox)
        self.groupBox_4 = QtWidgets.QGroupBox(DlgCalculateCE)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_4.setFont(font)
        self.groupBox_4.setObjectName("groupBox_4")
        self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.groupBox_4)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.cvi_dateEdit = QtWidgets.QDateEdit(self.groupBox_4)
        self.cvi_dateEdit.setMinimumSize(QtCore.QSize(0, 32))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cvi_dateEdit.setFont(font)
        self.cvi_dateEdit.setObjectName("cvi_dateEdit")
        self.horizontalLayout_4.addWidget(self.cvi_dateEdit)
        self.verticalLayout_2.addWidget(self.groupBox_4)
        self.groupBox_5 = QtWidgets.QGroupBox(DlgCalculateCE)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_5.setFont(font)
        self.groupBox_5.setObjectName("groupBox_5")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox_5)
        self.horizontalLayout_2.setSpacing(16)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.cvi_lineEdit = QtWidgets.QLineEdit(self.groupBox_5)
        self.cvi_lineEdit.setMinimumSize(QtCore.QSize(0, 32))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.cvi_lineEdit.setFont(font)
        self.cvi_lineEdit.setObjectName("cvi_lineEdit")
        self.horizontalLayout_2.addWidget(self.cvi_lineEdit)
        self.pushButton_cvi_getFolder = QtWidgets.QPushButton(self.groupBox_5)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_cvi_getFolder.sizePolicy().hasHeightForWidth())
        self.pushButton_cvi_getFolder.setSizePolicy(sizePolicy)
        self.pushButton_cvi_getFolder.setMinimumSize(QtCore.QSize(32, 32))
        self.pushButton_cvi_getFolder.setMaximumSize(QtCore.QSize(32, 16777215))
        self.pushButton_cvi_getFolder.setObjectName("pushButton_cvi_getFolder")
        self.horizontalLayout_2.addWidget(self.pushButton_cvi_getFolder)
        self.verticalLayout_2.addWidget(self.groupBox_5)
        self.groupBox_3 = QtWidgets.QGroupBox(DlgCalculateCE)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_3.setFont(font)
        self.groupBox_3.setObjectName("groupBox_3")
        self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.groupBox_3)
        self.horizontalLayout_5.setObjectName("horizontalLayout_5")
        self.pushButton_cvi_sub = QtWidgets.QPushButton(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_cvi_sub.sizePolicy().hasHeightForWidth())
        self.pushButton_cvi_sub.setSizePolicy(sizePolicy)
        self.pushButton_cvi_sub.setMinimumSize(QtCore.QSize(124, 32))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_cvi_sub.setFont(font)
        self.pushButton_cvi_sub.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.pushButton_cvi_sub.setObjectName("pushButton_cvi_sub")
        self.horizontalLayout_5.addWidget(self.pushButton_cvi_sub)
        self.verticalLayout_2.addWidget(self.groupBox_3)
        self.groupBox_2 = QtWidgets.QGroupBox(DlgCalculateCE)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.groupBox_2)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.pushButton_cvi = QtWidgets.QPushButton(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.pushButton_cvi.sizePolicy().hasHeightForWidth())
        self.pushButton_cvi.setSizePolicy(sizePolicy)
        self.pushButton_cvi.setMinimumSize(QtCore.QSize(124, 32))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.pushButton_cvi.setFont(font)
        self.pushButton_cvi.setLayoutDirection(QtCore.Qt.RightToLeft)
        self.pushButton_cvi.setObjectName("pushButton_cvi")
        self.horizontalLayout_3.addWidget(self.pushButton_cvi)
        self.verticalLayout_2.addWidget(self.groupBox_2)

        self.retranslateUi(DlgCalculateCE)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculateCE)

    def retranslateUi(self, DlgCalculateCE):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculateCE.setWindowTitle(_translate("DlgCalculateCE", "Calculate Indicators"))
        self.groupBox.setTitle(_translate("DlgCalculateCE", "Select AOI"))
        self.label.setText(_translate("DlgCalculateCE", "Select Country"))
        self.label_2.setText(_translate("DlgCalculateCE", "Select Region"))
        self.groupBox_4.setTitle(_translate("DlgCalculateCE", "Select Year"))
        self.cvi_dateEdit.setDisplayFormat(_translate("DlgCalculateCE", "yyyy"))
        self.groupBox_5.setTitle(_translate("DlgCalculateCE", "Select Folder"))
        self.pushButton_cvi_getFolder.setText(_translate("DlgCalculateCE", "..."))
        self.groupBox_3.setTitle(_translate("DlgCalculateCE", "Option 1: Compute CVI Sub-indicators"))
        self.pushButton_cvi_sub.setText(_translate("DlgCalculateCE", "Compute CVI Sub-indicators"))
        self.groupBox_2.setTitle(_translate("DlgCalculateCE", "Option 2: Compute Coastal Vulnerability Index (CVI)"))
        self.pushButton_cvi.setText(_translate("DlgCalculateCE", "Compute CVI"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculateCE = QtWidgets.QDialog()
    ui = Ui_DlgCalculateCE()
    ui.setupUi(DlgCalculateCE)
    DlgCalculateCE.show()
    sys.exit(app.exec_())
