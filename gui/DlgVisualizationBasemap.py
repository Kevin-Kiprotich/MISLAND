# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgVisualizationBasemap.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgVisualizationBasemap(object):
    def setupUi(self, DlgVisualizationBasemap):
        DlgVisualizationBasemap.setObjectName("DlgVisualizationBasemap")
        DlgVisualizationBasemap.resize(400, 275)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgVisualizationBasemap.sizePolicy().hasHeightForWidth())
        DlgVisualizationBasemap.setSizePolicy(sizePolicy)
        DlgVisualizationBasemap.setMinimumSize(QtCore.QSize(400, 275))
        self.verticalLayout = QtWidgets.QVBoxLayout(DlgVisualizationBasemap)
        self.verticalLayout.setObjectName("verticalLayout")
        self.area_groupbox = QtWidgets.QGroupBox(DlgVisualizationBasemap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.area_groupbox.sizePolicy().hasHeightForWidth())
        self.area_groupbox.setSizePolicy(sizePolicy)
        self.area_groupbox.setObjectName("area_groupbox")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.area_groupbox)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.checkBox_mask = QtWidgets.QCheckBox(self.area_groupbox)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.checkBox_mask.sizePolicy().hasHeightForWidth())
        self.checkBox_mask.setSizePolicy(sizePolicy)
        self.checkBox_mask.setObjectName("checkBox_mask")
        self.verticalLayout_3.addWidget(self.checkBox_mask)
        self.label_maskselect = QtWidgets.QLabel(self.area_groupbox)
        self.label_maskselect.setEnabled(False)
        self.label_maskselect.setObjectName("label_maskselect")
        self.verticalLayout_3.addWidget(self.label_maskselect)
        self.area_gridlayout = QtWidgets.QGridLayout()
        self.area_gridlayout.setObjectName("area_gridlayout")
        self.label_admin1 = QtWidgets.QLabel(self.area_groupbox)
        self.label_admin1.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_admin1.sizePolicy().hasHeightForWidth())
        self.label_admin1.setSizePolicy(sizePolicy)
        self.label_admin1.setMinimumSize(QtCore.QSize(150, 0))
        self.label_admin1.setMaximumSize(QtCore.QSize(110, 16777215))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_admin1.setFont(font)
        self.label_admin1.setObjectName("label_admin1")
        self.area_gridlayout.addWidget(self.label_admin1, 1, 0, 1, 1)
        self.area_admin_1 = QtWidgets.QComboBox(self.area_groupbox)
        self.area_admin_1.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.area_admin_1.sizePolicy().hasHeightForWidth())
        self.area_admin_1.setSizePolicy(sizePolicy)
        self.area_admin_1.setMinimumSize(QtCore.QSize(0, 30))
        self.area_admin_1.setObjectName("area_admin_1")
        self.area_gridlayout.addWidget(self.area_admin_1, 1, 1, 1, 2)
        self.label_admin0 = QtWidgets.QLabel(self.area_groupbox)
        self.label_admin0.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_admin0.sizePolicy().hasHeightForWidth())
        self.label_admin0.setSizePolicy(sizePolicy)
        self.label_admin0.setMinimumSize(QtCore.QSize(150, 0))
        self.label_admin0.setMaximumSize(QtCore.QSize(110, 16777215))
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_admin0.setFont(font)
        self.label_admin0.setObjectName("label_admin0")
        self.area_gridlayout.addWidget(self.label_admin0, 0, 0, 1, 1)
        self.area_admin_0 = QtWidgets.QComboBox(self.area_groupbox)
        self.area_admin_0.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.area_admin_0.sizePolicy().hasHeightForWidth())
        self.area_admin_0.setSizePolicy(sizePolicy)
        self.area_admin_0.setMinimumSize(QtCore.QSize(0, 30))
        self.area_admin_0.setObjectName("area_admin_0")
        self.area_gridlayout.addWidget(self.area_admin_0, 0, 1, 1, 2)
        self.verticalLayout_3.addLayout(self.area_gridlayout)
        self.verticalLayout.addWidget(self.area_groupbox)
        self.label_disclaimer = QtWidgets.QLabel(DlgVisualizationBasemap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_disclaimer.sizePolicy().hasHeightForWidth())
        self.label_disclaimer.setSizePolicy(sizePolicy)
        self.label_disclaimer.setWordWrap(True)
        self.label_disclaimer.setOpenExternalLinks(True)
        self.label_disclaimer.setObjectName("label_disclaimer")
        self.verticalLayout.addWidget(self.label_disclaimer)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)
        self.buttonBox = QtWidgets.QDialogButtonBox(DlgVisualizationBasemap)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout.addWidget(self.buttonBox)

        self.retranslateUi(DlgVisualizationBasemap)
        QtCore.QMetaObject.connectSlotsByName(DlgVisualizationBasemap)

    def retranslateUi(self, DlgVisualizationBasemap):
        _translate = QtCore.QCoreApplication.translate
        DlgVisualizationBasemap.setWindowTitle(_translate("DlgVisualizationBasemap", "Add basemap"))
        self.area_groupbox.setTitle(_translate("DlgVisualizationBasemap", "Mask "))
        self.checkBox_mask.setText(_translate("DlgVisualizationBasemap", "Use a mask"))
        self.label_maskselect.setText(_translate("DlgVisualizationBasemap", "Mask all areas outside of:"))
        self.label_admin1.setText(_translate("DlgVisualizationBasemap", "Second level:"))
        self.label_admin0.setText(_translate("DlgVisualizationBasemap", "First level:"))
        self.label_disclaimer.setText(_translate("DlgVisualizationBasemap", "<html><head/><body><p>Disclaimer: The provided boundaries are from <a href=\"http://www.naturalearthdata.com/\"><span style=\" text-decoration: underline; color:#0000ff;\">Natural Earth</span></a>, and are in the <a href=\"http://creativecommons.org/publicdomain\"><span style=\" text-decoration: underline; color:#0000ff;\">public domain</span></a>. The boundaries and names used, and the designations used, in trends.earth do not imply official endorsement or acceptance by Conservation International Foundation, or by its partner organizations and contributors.</p></body></html>"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgVisualizationBasemap = QtWidgets.QDialog()
    ui = Ui_DlgVisualizationBasemap()
    ui.setupUi(DlgVisualizationBasemap)
    DlgVisualizationBasemap.show()
    sys.exit(app.exec_())
