# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgDataIO.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgDataIO(object):
    def setupUi(self, DlgDataIO):
        DlgDataIO.setObjectName("DlgDataIO")
        DlgDataIO.setWindowModality(QtCore.Qt.ApplicationModal)
        DlgDataIO.resize(509, 668)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgDataIO.sizePolicy().hasHeightForWidth())
        DlgDataIO.setSizePolicy(sizePolicy)
        self.gridLayout = QtWidgets.QGridLayout(DlgDataIO)
        self.gridLayout.setObjectName("gridLayout")
        self.groupBox = QtWidgets.QGroupBox(DlgDataIO)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox.setFont(font)
        self.groupBox.setCursor(QtGui.QCursor(QtCore.Qt.ArrowCursor))
        self.groupBox.setObjectName("groupBox")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox_3 = QtWidgets.QGroupBox(self.groupBox)
        self.groupBox_3.setObjectName("groupBox_3")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox_3)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.btn_prod = QtWidgets.QPushButton(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_prod.sizePolicy().hasHeightForWidth())
        self.btn_prod.setSizePolicy(sizePolicy)
        self.btn_prod.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_prod.setFont(font)
        self.btn_prod.setStyleSheet("background-color: rgb(225, 213, 255);")
        self.btn_prod.setObjectName("btn_prod")
        self.horizontalLayout_2.addWidget(self.btn_prod)
        self.btn_soc = QtWidgets.QPushButton(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_soc.sizePolicy().hasHeightForWidth())
        self.btn_soc.setSizePolicy(sizePolicy)
        self.btn_soc.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_soc.setFont(font)
        self.btn_soc.setStyleSheet("background-color: rgb(225, 213, 255);")
        self.btn_soc.setObjectName("btn_soc")
        self.horizontalLayout_2.addWidget(self.btn_soc)
        self.btn_lc = QtWidgets.QPushButton(self.groupBox_3)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_lc.sizePolicy().hasHeightForWidth())
        self.btn_lc.setSizePolicy(sizePolicy)
        self.btn_lc.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_lc.setFont(font)
        self.btn_lc.setStyleSheet("background-color: rgb(225, 213, 255);")
        self.btn_lc.setObjectName("btn_lc")
        self.horizontalLayout_2.addWidget(self.btn_lc)
        self.verticalLayout_2.addWidget(self.groupBox_3)
        self.groupBox_5 = QtWidgets.QGroupBox(self.groupBox)
        self.groupBox_5.setObjectName("groupBox_5")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_5)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.groupBox_4 = QtWidgets.QGroupBox(self.groupBox_5)
        self.groupBox_4.setObjectName("groupBox_4")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.btn_ppt = QtWidgets.QPushButton(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_ppt.sizePolicy().hasHeightForWidth())
        self.btn_ppt.setSizePolicy(sizePolicy)
        self.btn_ppt.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_ppt.setFont(font)
        self.btn_ppt.setStyleSheet("background-color: rgb(213, 255, 252);")
        self.btn_ppt.setObjectName("btn_ppt")
        self.gridLayout_3.addWidget(self.btn_ppt, 0, 0, 1, 1)
        self.btn_pet = QtWidgets.QPushButton(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_pet.sizePolicy().hasHeightForWidth())
        self.btn_pet.setSizePolicy(sizePolicy)
        self.btn_pet.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_pet.setFont(font)
        self.btn_pet.setStyleSheet("background-color: rgb(213, 255, 252);")
        self.btn_pet.setObjectName("btn_pet")
        self.gridLayout_3.addWidget(self.btn_pet, 0, 1, 1, 1)
        self.btn_cqi = QtWidgets.QPushButton(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_cqi.sizePolicy().hasHeightForWidth())
        self.btn_cqi.setSizePolicy(sizePolicy)
        self.btn_cqi.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_cqi.setFont(font)
        self.btn_cqi.setStyleSheet("background-color: rgb(213, 255, 252);")
        self.btn_cqi.setObjectName("btn_cqi")
        self.gridLayout_3.addWidget(self.btn_cqi, 1, 0, 1, 2)
        self.verticalLayout_3.addWidget(self.groupBox_4)
        self.groupBox_6 = QtWidgets.QGroupBox(self.groupBox_5)
        self.groupBox_6.setObjectName("groupBox_6")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_6)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.btn_pm = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_pm.sizePolicy().hasHeightForWidth())
        self.btn_pm.setSizePolicy(sizePolicy)
        self.btn_pm.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_pm.setFont(font)
        self.btn_pm.setStyleSheet("background-color: rgb(247, 234, 211);")
        self.btn_pm.setObjectName("btn_pm")
        self.gridLayout_2.addWidget(self.btn_pm, 0, 0, 1, 1)
        self.btn_drain = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_drain.sizePolicy().hasHeightForWidth())
        self.btn_drain.setSizePolicy(sizePolicy)
        self.btn_drain.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_drain.setFont(font)
        self.btn_drain.setStyleSheet("background-color: rgb(247, 234, 211);")
        self.btn_drain.setObjectName("btn_drain")
        self.gridLayout_2.addWidget(self.btn_drain, 0, 1, 1, 1)
        self.btn_rock = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_rock.sizePolicy().hasHeightForWidth())
        self.btn_rock.setSizePolicy(sizePolicy)
        self.btn_rock.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_rock.setFont(font)
        self.btn_rock.setStyleSheet("background-color: rgb(247, 234, 211);")
        self.btn_rock.setObjectName("btn_rock")
        self.gridLayout_2.addWidget(self.btn_rock, 1, 0, 1, 1)
        self.btn_texture = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_texture.sizePolicy().hasHeightForWidth())
        self.btn_texture.setSizePolicy(sizePolicy)
        self.btn_texture.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_texture.setFont(font)
        self.btn_texture.setStyleSheet("background-color: rgb(247, 234, 211);")
        self.btn_texture.setObjectName("btn_texture")
        self.gridLayout_2.addWidget(self.btn_texture, 1, 1, 1, 1)
        self.btn_sqi = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_sqi.sizePolicy().hasHeightForWidth())
        self.btn_sqi.setSizePolicy(sizePolicy)
        self.btn_sqi.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_sqi.setFont(font)
        self.btn_sqi.setStyleSheet("background-color: rgb(247, 234, 211);")
        self.btn_sqi.setObjectName("btn_sqi")
        self.gridLayout_2.addWidget(self.btn_sqi, 2, 0, 1, 2)
        self.verticalLayout_3.addWidget(self.groupBox_6)
        self.groupBox_7 = QtWidgets.QGroupBox(self.groupBox_5)
        self.groupBox_7.setObjectName("groupBox_7")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_7)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.btn_vqi = QtWidgets.QPushButton(self.groupBox_7)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_vqi.sizePolicy().hasHeightForWidth())
        self.btn_vqi.setSizePolicy(sizePolicy)
        self.btn_vqi.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_vqi.setFont(font)
        self.btn_vqi.setStyleSheet("background-color: rgb(229, 255, 233);")
        self.btn_vqi.setObjectName("btn_vqi")
        self.gridLayout_4.addWidget(self.btn_vqi, 0, 0, 1, 1)
        self.verticalLayout_3.addWidget(self.groupBox_7)
        self.groupBox_8 = QtWidgets.QGroupBox(self.groupBox_5)
        self.groupBox_8.setObjectName("groupBox_8")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_8)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.btn_mqi = QtWidgets.QPushButton(self.groupBox_8)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_mqi.sizePolicy().hasHeightForWidth())
        self.btn_mqi.setSizePolicy(sizePolicy)
        self.btn_mqi.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_mqi.setFont(font)
        self.btn_mqi.setStyleSheet("background-color: rgb(254, 255, 234);")
        self.btn_mqi.setObjectName("btn_mqi")
        self.gridLayout_5.addWidget(self.btn_mqi, 0, 0, 1, 1)
        self.verticalLayout_3.addWidget(self.groupBox_8)
        self.verticalLayout_2.addWidget(self.groupBox_5)
        self.gridLayout.addWidget(self.groupBox, 1, 0, 1, 1)
        self.groupBox_2 = QtWidgets.QGroupBox(DlgDataIO)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        font = QtGui.QFont()
        font.setBold(True)
        font.setWeight(75)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout.setObjectName("verticalLayout")
        self.btn_te = QtWidgets.QPushButton(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_te.sizePolicy().hasHeightForWidth())
        self.btn_te.setSizePolicy(sizePolicy)
        self.btn_te.setMinimumSize(QtCore.QSize(0, 30))
        self.btn_te.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.btn_te.setFont(font)
        self.btn_te.setStyleSheet("")
        self.btn_te.setObjectName("btn_te")
        self.verticalLayout.addWidget(self.btn_te)
        self.gridLayout.addWidget(self.groupBox_2, 0, 0, 1, 1)

        self.retranslateUi(DlgDataIO)
        QtCore.QMetaObject.connectSlotsByName(DlgDataIO)

    def retranslateUi(self, DlgDataIO):
        _translate = QtCore.QCoreApplication.translate
        DlgDataIO.setWindowTitle(_translate("DlgDataIO", "Load data"))
        self.groupBox.setTitle(_translate("DlgDataIO", "Import a custom input dataset"))
        self.groupBox_3.setTitle(_translate("DlgDataIO", "SDG 15.3.1"))
        self.btn_prod.setText(_translate("DlgDataIO", "Productivity"))
        self.btn_soc.setText(_translate("DlgDataIO", "Soil organic carbon"))
        self.btn_lc.setText(_translate("DlgDataIO", "Land cover"))
        self.groupBox_5.setTitle(_translate("DlgDataIO", "MEDALUS"))
        self.groupBox_4.setTitle(_translate("DlgDataIO", "Climate Quality Index"))
        self.btn_ppt.setText(_translate("DlgDataIO", "Precipitation"))
        self.btn_pet.setText(_translate("DlgDataIO", "Evapotranspiration"))
        self.btn_cqi.setText(_translate("DlgDataIO", "Climate Quality Index Dataset"))
        self.groupBox_6.setTitle(_translate("DlgDataIO", "Soil Quality Index"))
        self.btn_pm.setText(_translate("DlgDataIO", "Parent Material"))
        self.btn_drain.setText(_translate("DlgDataIO", "Drainage"))
        self.btn_rock.setText(_translate("DlgDataIO", "Rock Fragment"))
        self.btn_texture.setText(_translate("DlgDataIO", "Soil Texture"))
        self.btn_sqi.setText(_translate("DlgDataIO", "Soil Quality Index Dataset"))
        self.groupBox_7.setTitle(_translate("DlgDataIO", "Vegetation Quality Index"))
        self.btn_vqi.setText(_translate("DlgDataIO", "Vegetation Quality Index Dataset"))
        self.groupBox_8.setTitle(_translate("DlgDataIO", "Management Quality Index"))
        self.btn_mqi.setText(_translate("DlgDataIO", "Management Quality Index Dataset"))
        self.groupBox_2.setTitle(_translate("DlgDataIO", "Load a dataset produced by MISLAND"))
        self.btn_te.setText(_translate("DlgDataIO", "Load an existing MISLAND output file"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgDataIO = QtWidgets.QDialog()
    ui = Ui_DlgDataIO()
    ui.setupUi(DlgDataIO)
    DlgDataIO.show()
    sys.exit(app.exec_())
