# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgVisualizationCreateMap.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgVisualizationCreateMap(object):
    def setupUi(self, DlgVisualizationCreateMap):
        DlgVisualizationCreateMap.setObjectName("DlgVisualizationCreateMap")
        DlgVisualizationCreateMap.resize(293, 314)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgVisualizationCreateMap.sizePolicy().hasHeightForWidth())
        DlgVisualizationCreateMap.setSizePolicy(sizePolicy)
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(DlgVisualizationCreateMap)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.layer_combo_label = QtWidgets.QLabel(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.layer_combo_label.sizePolicy().hasHeightForWidth())
        self.layer_combo_label.setSizePolicy(sizePolicy)
        self.layer_combo_label.setMinimumSize(QtCore.QSize(0, 30))
        self.layer_combo_label.setMaximumSize(QtCore.QSize(16777215, 33))
        self.layer_combo_label.setObjectName("layer_combo_label")
        self.verticalLayout_2.addWidget(self.layer_combo_label)
        self.combo_layers = QtWidgets.QComboBox(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.combo_layers.sizePolicy().hasHeightForWidth())
        self.combo_layers.setSizePolicy(sizePolicy)
        self.combo_layers.setMinimumSize(QtCore.QSize(275, 30))
        self.combo_layers.setObjectName("combo_layers")
        self.verticalLayout_2.addWidget(self.combo_layers)
        self.label = QtWidgets.QLabel(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setMinimumSize(QtCore.QSize(0, 30))
        self.label.setMaximumSize(QtCore.QSize(16777215, 33))
        self.label.setObjectName("label")
        self.verticalLayout_2.addWidget(self.label)
        self.title = QtWidgets.QLineEdit(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.title.sizePolicy().hasHeightForWidth())
        self.title.setSizePolicy(sizePolicy)
        self.title.setMinimumSize(QtCore.QSize(200, 30))
        self.title.setMaxLength(40)
        self.title.setObjectName("title")
        self.verticalLayout_2.addWidget(self.title)
        self.label_2 = QtWidgets.QLabel(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy)
        self.label_2.setMinimumSize(QtCore.QSize(0, 30))
        self.label_2.setMaximumSize(QtCore.QSize(16777215, 33))
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.authors = QtWidgets.QLineEdit(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.authors.sizePolicy().hasHeightForWidth())
        self.authors.setSizePolicy(sizePolicy)
        self.authors.setMinimumSize(QtCore.QSize(200, 30))
        self.authors.setMaxLength(175)
        self.authors.setObjectName("authors")
        self.verticalLayout_2.addWidget(self.authors)
        self.groupBox = QtWidgets.QGroupBox(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox.sizePolicy().hasHeightForWidth())
        self.groupBox.setSizePolicy(sizePolicy)
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayout = QtWidgets.QHBoxLayout(self.groupBox)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.portrait_layout = QtWidgets.QRadioButton(self.groupBox)
        self.portrait_layout.setChecked(True)
        self.portrait_layout.setObjectName("portrait_layout")
        self.horizontalLayout.addWidget(self.portrait_layout)
        self.landscape_layout = QtWidgets.QRadioButton(self.groupBox)
        self.landscape_layout.setObjectName("landscape_layout")
        self.horizontalLayout.addWidget(self.landscape_layout)
        self.verticalLayout_2.addWidget(self.groupBox)
        self.buttonBox = QtWidgets.QDialogButtonBox(DlgVisualizationCreateMap)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.buttonBox.sizePolicy().hasHeightForWidth())
        self.buttonBox.setSizePolicy(sizePolicy)
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.verticalLayout_2.addWidget(self.buttonBox)

        self.retranslateUi(DlgVisualizationCreateMap)
        QtCore.QMetaObject.connectSlotsByName(DlgVisualizationCreateMap)

    def retranslateUi(self, DlgVisualizationCreateMap):
        _translate = QtCore.QCoreApplication.translate
        DlgVisualizationCreateMap.setWindowTitle(_translate("DlgVisualizationCreateMap", "Create Map"))
        self.layer_combo_label.setText(_translate("DlgVisualizationCreateMap", "Layer to display:"))
        self.label.setText(_translate("DlgVisualizationCreateMap", "Map name:"))
        self.label_2.setText(_translate("DlgVisualizationCreateMap", "Author(s):"))
        self.groupBox.setTitle(_translate("DlgVisualizationCreateMap", "Layout"))
        self.portrait_layout.setText(_translate("DlgVisualizationCreateMap", "Portrait"))
        self.landscape_layout.setText(_translate("DlgVisualizationCreateMap", "Landscape"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgVisualizationCreateMap = QtWidgets.QDialog()
    ui = Ui_DlgVisualizationCreateMap()
    ui.setupUi(DlgVisualizationCreateMap)
    DlgVisualizationCreateMap.show()
    sys.exit(app.exec_())
