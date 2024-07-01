# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgCalculateRestBiomassSummaryTable.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculateRestBiomassSummaryTable(object):
    def setupUi(self, DlgCalculateRestBiomassSummaryTable):
        DlgCalculateRestBiomassSummaryTable.setObjectName("DlgCalculateRestBiomassSummaryTable")
        DlgCalculateRestBiomassSummaryTable.resize(293, 285)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgCalculateRestBiomassSummaryTable.sizePolicy().hasHeightForWidth())
        DlgCalculateRestBiomassSummaryTable.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgCalculateRestBiomassSummaryTable)
        self.verticalLayout.setObjectName("verticalLayout")
        self.TabBox = QtWidgets.QTabWidget(DlgCalculateRestBiomassSummaryTable)
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
        self.combo_layer_biomass_diff = WidgetDataIOSelectTELayerExisting(self.groupBox_2)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.combo_layer_biomass_diff.sizePolicy().hasHeightForWidth())
        self.combo_layer_biomass_diff.setSizePolicy(sizePolicy)
        self.combo_layer_biomass_diff.setMinimumSize(QtCore.QSize(0, 30))
        self.combo_layer_biomass_diff.setProperty("layer_type", "Biomass (tonnes CO2e per ha)")
        self.combo_layer_biomass_diff.setObjectName("combo_layer_biomass_diff")
        self.verticalLayout_6.addWidget(self.combo_layer_biomass_diff)
        self.verticalLayout_2.addWidget(self.groupBox_2)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_2.addItem(spacerItem)
        self.TabBox.addTab(self.InputTab, "")
        self.OutputTab = QtWidgets.QWidget()
        self.OutputTab.setObjectName("OutputTab")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.OutputTab)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.groupBox_6 = QtWidgets.QGroupBox(self.OutputTab)
        self.groupBox_6.setObjectName("groupBox_6")
        self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.groupBox_6)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.output_file_layer = QtWidgets.QLineEdit(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.output_file_layer.sizePolicy().hasHeightForWidth())
        self.output_file_layer.setSizePolicy(sizePolicy)
        self.output_file_layer.setMinimumSize(QtCore.QSize(150, 35))
        self.output_file_layer.setReadOnly(True)
        self.output_file_layer.setObjectName("output_file_layer")
        self.horizontalLayout_3.addWidget(self.output_file_layer)
        self.browse_output_file_layer = QtWidgets.QPushButton(self.groupBox_6)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.browse_output_file_layer.sizePolicy().hasHeightForWidth())
        self.browse_output_file_layer.setSizePolicy(sizePolicy)
        self.browse_output_file_layer.setMinimumSize(QtCore.QSize(0, 35))
        self.browse_output_file_layer.setObjectName("browse_output_file_layer")
        self.horizontalLayout_3.addWidget(self.browse_output_file_layer)
        self.verticalLayout_3.addWidget(self.groupBox_6)
        self.groupBox_7 = QtWidgets.QGroupBox(self.OutputTab)
        self.groupBox_7.setObjectName("groupBox_7")
        self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.groupBox_7)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.output_file_table = QtWidgets.QLineEdit(self.groupBox_7)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.output_file_table.sizePolicy().hasHeightForWidth())
        self.output_file_table.setSizePolicy(sizePolicy)
        self.output_file_table.setMinimumSize(QtCore.QSize(75, 35))
        self.output_file_table.setReadOnly(True)
        self.output_file_table.setObjectName("output_file_table")
        self.horizontalLayout_2.addWidget(self.output_file_table)
        self.browse_output_file_table = QtWidgets.QPushButton(self.groupBox_7)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.browse_output_file_table.sizePolicy().hasHeightForWidth())
        self.browse_output_file_table.setSizePolicy(sizePolicy)
        self.browse_output_file_table.setMinimumSize(QtCore.QSize(0, 35))
        self.browse_output_file_table.setObjectName("browse_output_file_table")
        self.horizontalLayout_2.addWidget(self.browse_output_file_table)
        self.verticalLayout_3.addWidget(self.groupBox_7)
        spacerItem1 = QtWidgets.QSpacerItem(20, 0, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_3.addItem(spacerItem1)
        self.TabBox.addTab(self.OutputTab, "")
        self.verticalLayout.addWidget(self.TabBox)
        self.gridLayout_3 = QtWidgets.QGridLayout()
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.button_prev = QtWidgets.QPushButton(DlgCalculateRestBiomassSummaryTable)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_prev.sizePolicy().hasHeightForWidth())
        self.button_prev.setSizePolicy(sizePolicy)
        self.button_prev.setMinimumSize(QtCore.QSize(0, 30))
        self.button_prev.setObjectName("button_prev")
        self.gridLayout_3.addWidget(self.button_prev, 0, 0, 1, 1)
        self.button_next = QtWidgets.QPushButton(DlgCalculateRestBiomassSummaryTable)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_next.sizePolicy().hasHeightForWidth())
        self.button_next.setSizePolicy(sizePolicy)
        self.button_next.setMinimumSize(QtCore.QSize(0, 30))
        self.button_next.setObjectName("button_next")
        self.gridLayout_3.addWidget(self.button_next, 0, 1, 1, 1)
        self.button_calculate = QtWidgets.QPushButton(DlgCalculateRestBiomassSummaryTable)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_calculate.sizePolicy().hasHeightForWidth())
        self.button_calculate.setSizePolicy(sizePolicy)
        self.button_calculate.setMinimumSize(QtCore.QSize(0, 30))
        self.button_calculate.setObjectName("button_calculate")
        self.gridLayout_3.addWidget(self.button_calculate, 1, 0, 1, 2)
        self.verticalLayout.addLayout(self.gridLayout_3)

        self.retranslateUi(DlgCalculateRestBiomassSummaryTable)
        self.TabBox.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculateRestBiomassSummaryTable)
        DlgCalculateRestBiomassSummaryTable.setTabOrder(self.TabBox, self.combo_layer_biomass_diff)
        DlgCalculateRestBiomassSummaryTable.setTabOrder(self.combo_layer_biomass_diff, self.button_prev)
        DlgCalculateRestBiomassSummaryTable.setTabOrder(self.button_prev, self.button_next)
        DlgCalculateRestBiomassSummaryTable.setTabOrder(self.button_next, self.button_calculate)
        DlgCalculateRestBiomassSummaryTable.setTabOrder(self.button_calculate, self.browse_output_file_table)
        DlgCalculateRestBiomassSummaryTable.setTabOrder(self.browse_output_file_table, self.output_file_table)

    def retranslateUi(self, DlgCalculateRestBiomassSummaryTable):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculateRestBiomassSummaryTable.setWindowTitle(_translate("DlgCalculateRestBiomassSummaryTable", "Calculate Carbon Change Summary Table"))
        self.groupBox_2.setTitle(_translate("DlgCalculateRestBiomassSummaryTable", "Biomass change"))
        self.TabBox.setTabText(self.TabBox.indexOf(self.InputTab), _translate("DlgCalculateRestBiomassSummaryTable", "Input"))
        self.groupBox_6.setTitle(_translate("DlgCalculateRestBiomassSummaryTable", "Output file for biomass difference layers"))
        self.output_file_layer.setPlaceholderText(_translate("DlgCalculateRestBiomassSummaryTable", "Click \"Browse\" to choose a file..."))
        self.browse_output_file_layer.setText(_translate("DlgCalculateRestBiomassSummaryTable", "Browse"))
        self.groupBox_7.setTitle(_translate("DlgCalculateRestBiomassSummaryTable", "Output file for summary table"))
        self.output_file_table.setPlaceholderText(_translate("DlgCalculateRestBiomassSummaryTable", "Click \"Browse\" to choose a file..."))
        self.browse_output_file_table.setText(_translate("DlgCalculateRestBiomassSummaryTable", "Browse"))
        self.TabBox.setTabText(self.TabBox.indexOf(self.OutputTab), _translate("DlgCalculateRestBiomassSummaryTable", "Output"))
        self.button_prev.setText(_translate("DlgCalculateRestBiomassSummaryTable", "Previous"))
        self.button_next.setText(_translate("DlgCalculateRestBiomassSummaryTable", "Next"))
        self.button_calculate.setText(_translate("DlgCalculateRestBiomassSummaryTable", "Calculate"))
from MISLAND.data_io import WidgetDataIOSelectTELayerExisting


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculateRestBiomassSummaryTable = QtWidgets.QDialog()
    ui = Ui_DlgCalculateRestBiomassSummaryTable()
    ui.setupUi(DlgCalculateRestBiomassSummaryTable)
    DlgCalculateRestBiomassSummaryTable.show()
    sys.exit(app.exec_())
