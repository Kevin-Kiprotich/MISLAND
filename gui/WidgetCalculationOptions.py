# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\WidgetCalculationOptions.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_WidgetCalculationOptions(object):
    def setupUi(self, WidgetCalculationOptions):
        WidgetCalculationOptions.setObjectName("WidgetCalculationOptions")
        WidgetCalculationOptions.resize(535, 416)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(WidgetCalculationOptions.sizePolicy().hasHeightForWidth())
        WidgetCalculationOptions.setSizePolicy(sizePolicy)
        self.verticalLayout = QtWidgets.QVBoxLayout(WidgetCalculationOptions)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox_where_to_run = QtWidgets.QGroupBox(WidgetCalculationOptions)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_where_to_run.sizePolicy().hasHeightForWidth())
        self.groupBox_where_to_run.setSizePolicy(sizePolicy)
        self.groupBox_where_to_run.setObjectName("groupBox_where_to_run")
        self.gridLayout = QtWidgets.QGridLayout(self.groupBox_where_to_run)
        self.gridLayout.setObjectName("gridLayout")
        self.radioButton_run_in_cloud = QtWidgets.QRadioButton(self.groupBox_where_to_run)
        self.radioButton_run_in_cloud.setChecked(True)
        self.radioButton_run_in_cloud.setObjectName("radioButton_run_in_cloud")
        self.gridLayout.addWidget(self.radioButton_run_in_cloud, 0, 0, 1, 1)
        self.radioButton_run_locally = QtWidgets.QRadioButton(self.groupBox_where_to_run)
        self.radioButton_run_locally.setObjectName("radioButton_run_locally")
        self.gridLayout.addWidget(self.radioButton_run_locally, 1, 0, 1, 1)
        self.lineEdit_local_data_folder = QtWidgets.QLineEdit(self.groupBox_where_to_run)
        self.lineEdit_local_data_folder.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lineEdit_local_data_folder.sizePolicy().hasHeightForWidth())
        self.lineEdit_local_data_folder.setSizePolicy(sizePolicy)
        self.lineEdit_local_data_folder.setMinimumSize(QtCore.QSize(60, 30))
        self.lineEdit_local_data_folder.setReadOnly(True)
        self.lineEdit_local_data_folder.setObjectName("lineEdit_local_data_folder")
        self.gridLayout.addWidget(self.lineEdit_local_data_folder, 2, 0, 1, 1)
        self.btn_local_data_folder_browse = QtWidgets.QPushButton(self.groupBox_where_to_run)
        self.btn_local_data_folder_browse.setEnabled(False)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_local_data_folder_browse.sizePolicy().hasHeightForWidth())
        self.btn_local_data_folder_browse.setSizePolicy(sizePolicy)
        self.btn_local_data_folder_browse.setMinimumSize(QtCore.QSize(0, 30))
        self.btn_local_data_folder_browse.setMaximumSize(QtCore.QSize(80, 16777215))
        self.btn_local_data_folder_browse.setObjectName("btn_local_data_folder_browse")
        self.gridLayout.addWidget(self.btn_local_data_folder_browse, 2, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_where_to_run)
        self.groupBox_metadata = QtWidgets.QGroupBox(WidgetCalculationOptions)
        self.groupBox_metadata.setObjectName("groupBox_metadata")
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.groupBox_metadata)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.label = QtWidgets.QLabel(self.groupBox_metadata)
        self.label.setObjectName("label")
        self.verticalLayout_2.addWidget(self.label)
        self.task_name = QtWidgets.QLineEdit(self.groupBox_metadata)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.task_name.sizePolicy().hasHeightForWidth())
        self.task_name.setSizePolicy(sizePolicy)
        self.task_name.setMinimumSize(QtCore.QSize(0, 30))
        self.task_name.setMaxLength(60)
        self.task_name.setObjectName("task_name")
        self.verticalLayout_2.addWidget(self.task_name)
        self.label_2 = QtWidgets.QLabel(self.groupBox_metadata)
        self.label_2.setObjectName("label_2")
        self.verticalLayout_2.addWidget(self.label_2)
        self.task_notes = QtWidgets.QTextEdit(self.groupBox_metadata)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.task_notes.sizePolicy().hasHeightForWidth())
        self.task_notes.setSizePolicy(sizePolicy)
        self.task_notes.setMinimumSize(QtCore.QSize(0, 0))
        self.task_notes.setMaximumSize(QtCore.QSize(16777215, 100))
        self.task_notes.setObjectName("task_notes")
        self.verticalLayout_2.addWidget(self.task_notes)
        self.task_notes.raise_()
        self.task_name.raise_()
        self.label_2.raise_()
        self.label.raise_()
        self.verticalLayout.addWidget(self.groupBox_metadata)
        spacerItem = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout.addItem(spacerItem)

        self.retranslateUi(WidgetCalculationOptions)
        QtCore.QMetaObject.connectSlotsByName(WidgetCalculationOptions)

    def retranslateUi(self, WidgetCalculationOptions):
        _translate = QtCore.QCoreApplication.translate
        WidgetCalculationOptions.setWindowTitle(_translate("WidgetCalculationOptions", "Form"))
        self.groupBox_where_to_run.setTitle(_translate("WidgetCalculationOptions", "Where to run calculations"))
        self.radioButton_run_in_cloud.setText(_translate("WidgetCalculationOptions", "Run in cloud"))
        self.radioButton_run_locally.setText(_translate("WidgetCalculationOptions", "Run locally"))
        self.lineEdit_local_data_folder.setPlaceholderText(_translate("WidgetCalculationOptions", "Click \"Browse\" to select where local data is saved..."))
        self.btn_local_data_folder_browse.setText(_translate("WidgetCalculationOptions", "Browse"))
        self.groupBox_metadata.setTitle(_translate("WidgetCalculationOptions", "Metadata"))
        self.label.setText(_translate("WidgetCalculationOptions", "Task name:"))
        self.label_2.setText(_translate("WidgetCalculationOptions", "Notes:"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    WidgetCalculationOptions = QtWidgets.QWidget()
    ui = Ui_WidgetCalculationOptions()
    ui.setupUi(WidgetCalculationOptions)
    WidgetCalculationOptions.show()
    sys.exit(app.exec_())
