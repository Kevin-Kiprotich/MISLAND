# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgCalculateSDISummaryTableAdmin.ui'
#
# Created by: PyQt5 UI code generator 5.15.2
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculateSDISummaryTableAdmin(object):
    def setupUi(self, DlgCalculateSDISummaryTableAdmin):
        DlgCalculateSDISummaryTableAdmin.setObjectName("DlgCalculateSDISummaryTableAdmin")
        DlgCalculateSDISummaryTableAdmin.resize(719, 517)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgCalculateSDISummaryTableAdmin.sizePolicy().hasHeightForWidth())
        DlgCalculateSDISummaryTableAdmin.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgCalculateSDISummaryTableAdmin)
        self.verticalLayout.setObjectName("verticalLayout")
        self.TabBox = QtWidgets.QTabWidget(DlgCalculateSDISummaryTableAdmin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.TabBox.sizePolicy().hasHeightForWidth())
        self.TabBox.setSizePolicy(sizePolicy)
        self.TabBox.setTabShape(QtWidgets.QTabWidget.Rounded)
        self.TabBox.setObjectName("TabBox")
        self.InputTab = QtWidgets.QWidget()
        self.InputTab.setObjectName("InputTab")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.InputTab)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.groupBox = QtWidgets.QGroupBox(self.InputTab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.groupBox.setFont(font)
        self.groupBox.setObjectName("groupBox")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox)
        self.gridLayout.setObjectName("gridLayout")
        self.combo_layer_sqi = WidgetDataIOSelectTELayerExisting(self.groupBox)
        self.combo_layer_sqi.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.combo_layer_sqi.sizePolicy().hasHeightForWidth())
        self.combo_layer_sqi.setSizePolicy(sizePolicy)
        self.combo_layer_sqi.setMinimumSize(QtCore.QSize(0, 30))
        self.combo_layer_sqi.setProperty("layer_type", "Soil Quality Index (cm deep)")
        self.combo_layer_sqi.setObjectName("combo_layer_sqi")
        self.gridLayout.addWidget(self.combo_layer_sqi, 2, 0, 1, 1)
        self.verticalLayout_2.addWidget(self.groupBox)
        self.groupBox_4 = QtWidgets.QGroupBox(self.InputTab)
        self.groupBox_4.setObjectName("groupBox_4")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_4)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.combo_layer_mqi = WidgetDataIOSelectTELayerExisting(self.groupBox_4)
        self.combo_layer_mqi.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.combo_layer_mqi.sizePolicy().hasHeightForWidth())
        self.combo_layer_mqi.setSizePolicy(sizePolicy)
        self.combo_layer_mqi.setMinimumSize(QtCore.QSize(0, 30))
        self.combo_layer_mqi.setProperty("layer_type", "Vegetation Quality Index")
        self.combo_layer_mqi.setObjectName("combo_layer_mqi")
        self.gridLayout_4.addWidget(self.combo_layer_mqi, 0, 0, 1, 1)
        self.verticalLayout_2.addWidget(self.groupBox_4)
        self.groupBox_3 = QtWidgets.QGroupBox(self.InputTab)
        self.groupBox_3.setObjectName("groupBox_3")
        self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_3)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.combo_layer_cqi = WidgetDataIOSelectTELayerExisting(self.groupBox_3)
        self.combo_layer_cqi.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.combo_layer_cqi.sizePolicy().hasHeightForWidth())
        self.combo_layer_cqi.setSizePolicy(sizePolicy)
        self.combo_layer_cqi.setMinimumSize(QtCore.QSize(0, 30))
        self.combo_layer_cqi.setProperty("layer_type", "Climate Quality Index (year)")
        self.combo_layer_cqi.setObjectName("combo_layer_cqi")
        self.gridLayout_2.addWidget(self.combo_layer_cqi, 0, 0, 1, 1)
        self.verticalLayout_2.addWidget(self.groupBox_3)
        self.groupBox_2 = QtWidgets.QGroupBox(self.InputTab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_2.sizePolicy().hasHeightForWidth())
        self.groupBox_2.setSizePolicy(sizePolicy)
        self.groupBox_2.setMinimumSize(QtCore.QSize(0, 0))
        self.groupBox_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.groupBox_2.setFont(font)
        self.groupBox_2.setObjectName("groupBox_2")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.groupBox_2)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.combo_layer_vqi = WidgetDataIOSelectTELayerExisting(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.combo_layer_vqi.sizePolicy().hasHeightForWidth())
        self.combo_layer_vqi.setSizePolicy(sizePolicy)
        self.combo_layer_vqi.setMinimumSize(QtCore.QSize(0, 30))
        self.combo_layer_vqi.setProperty("layer_type", "Management Quality Index")
        self.combo_layer_vqi.setObjectName("combo_layer_vqi")
        self.verticalLayout_6.addWidget(self.combo_layer_vqi)
        self.verticalLayout_2.addWidget(self.groupBox_2)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.TabBox.addTab(self.InputTab, "")
        self.verticalLayout.addWidget(self.TabBox)
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.button_prev = QtWidgets.QPushButton(DlgCalculateSDISummaryTableAdmin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_prev.sizePolicy().hasHeightForWidth())
        self.button_prev.setSizePolicy(sizePolicy)
        self.button_prev.setMinimumSize(QtCore.QSize(0, 30))
        self.button_prev.setObjectName("button_prev")
        self.gridLayout_3.addWidget(self.button_prev, 0, 0, 1, 1)
        self.button_next = QtWidgets.QPushButton(DlgCalculateSDISummaryTableAdmin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_next.sizePolicy().hasHeightForWidth())
        self.button_next.setSizePolicy(sizePolicy)
        self.button_next.setMinimumSize(QtCore.QSize(0, 30))
        self.button_next.setObjectName("button_next")
        self.gridLayout_3.addWidget(self.button_next, 0, 1, 1, 1)
        self.button_calculate = QtWidgets.QPushButton(DlgCalculateSDISummaryTableAdmin)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_calculate.sizePolicy().hasHeightForWidth())
        self.button_calculate.setSizePolicy(sizePolicy)
        self.button_calculate.setMinimumSize(QtCore.QSize(0, 30))
        self.button_calculate.setObjectName("button_calculate")
        self.gridLayout_3.addWidget(self.button_calculate, 1, 0, 1, 2)
        self.verticalLayout.addLayout(self.gridLayout_3)

        self.retranslateUi(DlgCalculateSDISummaryTableAdmin)
        self.TabBox.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculateSDISummaryTableAdmin)
        DlgCalculateSDISummaryTableAdmin.setTabOrder(self.TabBox, self.combo_layer_vqi)
        DlgCalculateSDISummaryTableAdmin.setTabOrder(self.combo_layer_vqi, self.button_prev)
        DlgCalculateSDISummaryTableAdmin.setTabOrder(self.button_prev, self.button_next)
        DlgCalculateSDISummaryTableAdmin.setTabOrder(self.button_next, self.button_calculate)
        DlgCalculateSDISummaryTableAdmin.setTabOrder(self.button_calculate, self.combo_layer_sqi)

    def retranslateUi(self, DlgCalculateSDISummaryTableAdmin):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculateSDISummaryTableAdmin.setWindowTitle(_translate("DlgCalculateSDISummaryTableAdmin", "Calculate Sensitivity Desertification Index Indicator"))
        self.groupBox.setTitle(_translate("DlgCalculateSDISummaryTableAdmin", "Soil Quality Index (SQI):"))
        self.groupBox_4.setTitle(_translate("DlgCalculateSDISummaryTableAdmin", "Vegetation Quality Index (VQI):"))
        self.groupBox_3.setTitle(_translate("DlgCalculateSDISummaryTableAdmin", "Climate Quality Index (CQI):"))
        self.groupBox_2.setTitle(_translate("DlgCalculateSDISummaryTableAdmin", "Management Quality Index (MQI):"))
        self.TabBox.setTabText(self.TabBox.indexOf(self.InputTab), _translate("DlgCalculateSDISummaryTableAdmin", "Input"))
        self.button_prev.setText(_translate("DlgCalculateSDISummaryTableAdmin", "Previous"))
        self.button_next.setText(_translate("DlgCalculateSDISummaryTableAdmin", "Next"))
        self.button_calculate.setText(_translate("DlgCalculateSDISummaryTableAdmin", "Calculate"))
from MISLAND.data_io import WidgetDataIOSelectTELayerExisting


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculateSDISummaryTableAdmin = QtWidgets.QDialog()
    ui = Ui_DlgCalculateSDISummaryTableAdmin()
    ui.setupUi(DlgCalculateSDISummaryTableAdmin)
    DlgCalculateSDISummaryTableAdmin.show()
    sys.exit(app.exec_())