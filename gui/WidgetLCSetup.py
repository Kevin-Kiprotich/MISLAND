# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\WidgetLCSetup.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_WidgetLCSetup(object):
    def setupUi(self, WidgetLCSetup):
        WidgetLCSetup.setObjectName("WidgetLCSetup")
        WidgetLCSetup.resize(512, 461)
        self.verticalLayout = QtWidgets.QVBoxLayout(WidgetLCSetup)
        self.verticalLayout.setObjectName("verticalLayout")
        self.use_esa = QtWidgets.QRadioButton(WidgetLCSetup)
        self.use_esa.setChecked(True)
        self.use_esa.setObjectName("use_esa")
        self.verticalLayout.addWidget(self.use_esa)
        self.groupBox_esa_period = QtWidgets.QGroupBox(WidgetLCSetup)
        self.groupBox_esa_period.setObjectName("groupBox_esa_period")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.groupBox_esa_period)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.use_esa_tg_year = QtWidgets.QDateEdit(self.groupBox_esa_period)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_esa_tg_year.sizePolicy().hasHeightForWidth())
        self.use_esa_tg_year.setSizePolicy(sizePolicy)
        self.use_esa_tg_year.setMinimumSize(QtCore.QSize(0, 30))
        self.use_esa_tg_year.setMaximumSize(QtCore.QSize(120, 35))
        self.use_esa_tg_year.setMinimumDateTime(QtCore.QDateTime(QtCore.QDate(1992, 1, 1), QtCore.QTime(0, 0, 0)))
        self.use_esa_tg_year.setMaximumDate(QtCore.QDate(2030, 12, 31))
        self.use_esa_tg_year.setMinimumDate(QtCore.QDate(1992, 1, 1))
        self.use_esa_tg_year.setDisplayFormat("yyyy")
        self.use_esa_tg_year.setDate(QtCore.QDate(2015, 1, 1))
        self.use_esa_tg_year.setObjectName("use_esa_tg_year")
        self.gridLayout_8.addWidget(self.use_esa_tg_year, 1, 1, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBox_esa_period)
        self.label_10.setMinimumSize(QtCore.QSize(0, 30))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.gridLayout_8.addWidget(self.label_10, 0, 0, 1, 1)
        self.label_3 = QtWidgets.QLabel(self.groupBox_esa_period)
        self.label_3.setMinimumSize(QtCore.QSize(0, 30))
        self.label_3.setObjectName("label_3")
        self.gridLayout_8.addWidget(self.label_3, 0, 1, 1, 1)
        self.use_esa_bl_year = QtWidgets.QDateEdit(self.groupBox_esa_period)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_esa_bl_year.sizePolicy().hasHeightForWidth())
        self.use_esa_bl_year.setSizePolicy(sizePolicy)
        self.use_esa_bl_year.setMinimumSize(QtCore.QSize(0, 30))
        self.use_esa_bl_year.setMaximumSize(QtCore.QSize(120, 35))
        self.use_esa_bl_year.setMaximumDate(QtCore.QDate(2015, 12, 31))
        self.use_esa_bl_year.setMinimumDate(QtCore.QDate(1992, 1, 1))
        self.use_esa_bl_year.setDisplayFormat("yyyy")
        self.use_esa_bl_year.setDate(QtCore.QDate(2001, 1, 1))
        self.use_esa_bl_year.setObjectName("use_esa_bl_year")
        self.gridLayout_8.addWidget(self.use_esa_bl_year, 1, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_esa_period)
        self.groupBox_esa_agg = QtWidgets.QGroupBox(WidgetLCSetup)
        self.groupBox_esa_agg.setCheckable(False)
        self.groupBox_esa_agg.setChecked(False)
        self.groupBox_esa_agg.setObjectName("groupBox_esa_agg")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox_esa_agg)
        self.gridLayout.setObjectName("gridLayout")
        self.use_esa_agg_edit = QtWidgets.QPushButton(self.groupBox_esa_agg)
        self.use_esa_agg_edit.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_esa_agg_edit.sizePolicy().hasHeightForWidth())
        self.use_esa_agg_edit.setSizePolicy(sizePolicy)
        self.use_esa_agg_edit.setMinimumSize(QtCore.QSize(100, 30))
        self.use_esa_agg_edit.setMaximumSize(QtCore.QSize(150, 16777215))
        self.use_esa_agg_edit.setCheckable(False)
        self.use_esa_agg_edit.setObjectName("use_esa_agg_edit")
        self.gridLayout.addWidget(self.use_esa_agg_edit, 1, 0, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_esa_agg)
        self.use_custom = QtWidgets.QRadioButton(WidgetLCSetup)
        self.use_custom.setObjectName("use_custom")
        self.verticalLayout.addWidget(self.use_custom)
        self.groupBox_custom_bl = QtWidgets.QGroupBox(WidgetLCSetup)
        self.groupBox_custom_bl.setEnabled(True)
        self.groupBox_custom_bl.setObjectName("groupBox_custom_bl")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_custom_bl)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.use_custom_initial = WidgetDataIOSelectTELayerImport(self.groupBox_custom_bl)
        self.use_custom_initial.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_custom_initial.sizePolicy().hasHeightForWidth())
        self.use_custom_initial.setSizePolicy(sizePolicy)
        self.use_custom_initial.setMinimumSize(QtCore.QSize(0, 30))
        self.use_custom_initial.setProperty("layer_type", "Land cover (7 class)")
        self.use_custom_initial.setObjectName("use_custom_initial")
        self.gridLayout_4.addWidget(self.use_custom_initial, 0, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_custom_bl)
        self.groupBox_custom_tg = QtWidgets.QGroupBox(WidgetLCSetup)
        self.groupBox_custom_tg.setEnabled(True)
        self.groupBox_custom_tg.setObjectName("groupBox_custom_tg")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.groupBox_custom_tg)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.use_custom_final = WidgetDataIOSelectTELayerImport(self.groupBox_custom_tg)
        self.use_custom_final.setEnabled(True)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.use_custom_final.sizePolicy().hasHeightForWidth())
        self.use_custom_final.setSizePolicy(sizePolicy)
        self.use_custom_final.setMinimumSize(QtCore.QSize(0, 30))
        self.use_custom_final.setProperty("layer_type", "Land cover (7 class)")
        self.use_custom_final.setObjectName("use_custom_final")
        self.gridLayout_6.addWidget(self.use_custom_final, 0, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_custom_tg)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)

        self.retranslateUi(WidgetLCSetup)
        QtCore.QMetaObject.connectSlotsByName(WidgetLCSetup)

    def retranslateUi(self, WidgetLCSetup):
        _translate = QtCore.QCoreApplication.translate
        WidgetLCSetup.setWindowTitle(_translate("WidgetLCSetup", "Form"))
        self.use_esa.setText(_translate("WidgetLCSetup", "European Space Agency CCI-LC (default land cover dataset for UNCCD reporting)"))
        self.groupBox_esa_period.setTitle(_translate("WidgetLCSetup", "Period"))
        self.label_10.setText(_translate("WidgetLCSetup", "Initial year:"))
        self.label_3.setText(_translate("WidgetLCSetup", "Target year:"))
        self.groupBox_esa_agg.setTitle(_translate("WidgetLCSetup", "Customize land cover aggregation method"))
        self.use_esa_agg_edit.setText(_translate("WidgetLCSetup", "Edit definition"))
        self.use_custom.setText(_translate("WidgetLCSetup", "Custom land cover dataset"))
        self.groupBox_custom_bl.setTitle(_translate("WidgetLCSetup", "Initial layer (initial year)"))
        self.groupBox_custom_tg.setTitle(_translate("WidgetLCSetup", "Final layer (target year)"))
from MISLAND.data_io import WidgetDataIOSelectTELayerImport


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    WidgetLCSetup = QtWidgets.QWidget()
    ui = Ui_WidgetLCSetup()
    ui.setupUi(WidgetLCSetup)
    WidgetLCSetup.show()
    sys.exit(app.exec_())
