# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MISLAND/gui\DlgCalculateProd.ui'
#
# Created by: PyQt5 UI code generator 5.15.0
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_DlgCalculateProd(object):
    def setupUi(self, DlgCalculateProd):
        DlgCalculateProd.setObjectName("DlgCalculateProd")
        DlgCalculateProd.resize(608, 950)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(DlgCalculateProd.sizePolicy().hasHeightForWidth())
        DlgCalculateProd.setSizePolicy(sizePolicy)
        DlgCalculateProd.setMinimumSize(QtCore.QSize(0, 450))
        self.verticalLayout_2 = QtWidgets.QVBoxLayout(DlgCalculateProd)
        self.verticalLayout_2.setObjectName("verticalLayout_2")
        self.TabBox = QtWidgets.QTabWidget(DlgCalculateProd)
        self.TabBox.setMinimumSize(QtCore.QSize(0, 0))
        self.TabBox.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.TabBox.setObjectName("TabBox")
        self.IndicatorsTab = QtWidgets.QWidget()
        self.IndicatorsTab.setObjectName("IndicatorsTab")
        self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.IndicatorsTab)
        self.verticalLayout_4.setObjectName("verticalLayout_4")
        self.groupBox_lpd_mode = QtWidgets.QGroupBox(self.IndicatorsTab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_lpd_mode.sizePolicy().hasHeightForWidth())
        self.groupBox_lpd_mode.setSizePolicy(sizePolicy)
        self.groupBox_lpd_mode.setObjectName("groupBox_lpd_mode")
        self.gridLayout_6 = QtWidgets.QGridLayout(self.groupBox_lpd_mode)
        self.gridLayout_6.setObjectName("gridLayout_6")
        self.mode_lpd_jrc = QtWidgets.QRadioButton(self.groupBox_lpd_mode)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mode_lpd_jrc.sizePolicy().hasHeightForWidth())
        self.mode_lpd_jrc.setSizePolicy(sizePolicy)
        self.mode_lpd_jrc.setText("")
        self.mode_lpd_jrc.setChecked(False)
        self.mode_lpd_jrc.setObjectName("mode_lpd_jrc")
        self.gridLayout_6.addWidget(self.mode_lpd_jrc, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(self.groupBox_lpd_mode)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy)
        self.label.setWordWrap(True)
        self.label.setObjectName("label")
        self.gridLayout_6.addWidget(self.label, 1, 1, 1, 1)
        self.mode_te_prod = QtWidgets.QRadioButton(self.groupBox_lpd_mode)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.mode_te_prod.sizePolicy().hasHeightForWidth())
        self.mode_te_prod.setSizePolicy(sizePolicy)
        self.mode_te_prod.setChecked(True)
        self.mode_te_prod.setObjectName("mode_te_prod")
        self.gridLayout_6.addWidget(self.mode_te_prod, 0, 0, 1, 2)
        self.verticalLayout_4.addWidget(self.groupBox_lpd_mode)
        self.groupBox_ndvi_dataset = QtWidgets.QGroupBox(self.IndicatorsTab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_ndvi_dataset.sizePolicy().hasHeightForWidth())
        self.groupBox_ndvi_dataset.setSizePolicy(sizePolicy)
        self.groupBox_ndvi_dataset.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.groupBox_ndvi_dataset.setFlat(True)
        self.groupBox_ndvi_dataset.setObjectName("groupBox_ndvi_dataset")
        self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.groupBox_ndvi_dataset)
        self.verticalLayout_3.setObjectName("verticalLayout_3")
        self.dataset_ndvi = QtWidgets.QComboBox(self.groupBox_ndvi_dataset)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.dataset_ndvi.sizePolicy().hasHeightForWidth())
        self.dataset_ndvi.setSizePolicy(sizePolicy)
        self.dataset_ndvi.setMinimumSize(QtCore.QSize(0, 40))
        self.dataset_ndvi.setObjectName("dataset_ndvi")
        self.verticalLayout_3.addWidget(self.dataset_ndvi)
        self.verticalLayout_4.addWidget(self.groupBox_ndvi_dataset)
        spacerItem = QtWidgets.QSpacerItem(20, 0, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_4.addItem(spacerItem)
        self.TabBox.addTab(self.IndicatorsTab, "")
        self.AdvancedTab = QtWidgets.QWidget()
        self.AdvancedTab.setObjectName("AdvancedTab")
        self.verticalLayout_7 = QtWidgets.QVBoxLayout(self.AdvancedTab)
        self.verticalLayout_7.setObjectName("verticalLayout_7")
        self.scrollArea = QtWidgets.QScrollArea(self.AdvancedTab)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.scrollArea.sizePolicy().hasHeightForWidth())
        self.scrollArea.setSizePolicy(sizePolicy)
        self.scrollArea.setFrameShape(QtWidgets.QFrame.NoFrame)
        self.scrollArea.setFrameShadow(QtWidgets.QFrame.Plain)
        self.scrollArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 568, 793))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.verticalLayout_9 = QtWidgets.QVBoxLayout(self.scrollAreaWidgetContents)
        self.verticalLayout_9.setObjectName("verticalLayout_9")
        self.groupBox_traj = QtWidgets.QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_traj.setCheckable(True)
        self.groupBox_traj.setChecked(True)
        self.groupBox_traj.setObjectName("groupBox_traj")
        self.verticalLayout = QtWidgets.QVBoxLayout(self.groupBox_traj)
        self.verticalLayout.setObjectName("verticalLayout")
        self.groupBox_4 = QtWidgets.QGroupBox(self.groupBox_traj)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.groupBox_4.sizePolicy().hasHeightForWidth())
        self.groupBox_4.setSizePolicy(sizePolicy)
        self.groupBox_4.setFlat(True)
        self.groupBox_4.setObjectName("groupBox_4")
        self.verticalLayout_6 = QtWidgets.QVBoxLayout(self.groupBox_4)
        self.verticalLayout_6.setObjectName("verticalLayout_6")
        self.traj_indic = QtWidgets.QComboBox(self.groupBox_4)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.traj_indic.sizePolicy().hasHeightForWidth())
        self.traj_indic.setSizePolicy(sizePolicy)
        self.traj_indic.setMinimumSize(QtCore.QSize(0, 40))
        self.traj_indic.setObjectName("traj_indic")
        self.verticalLayout_6.addWidget(self.traj_indic)
        self.verticalLayout.addWidget(self.groupBox_4)
        self.groupBox_traj_period = QtWidgets.QGroupBox(self.groupBox_traj)
        self.groupBox_traj_period.setFlat(True)
        self.groupBox_traj_period.setObjectName("groupBox_traj_period")
        self.gridLayout_5 = QtWidgets.QGridLayout(self.groupBox_traj_period)
        self.gridLayout_5.setObjectName("gridLayout_5")
        self.traj_year_start = QtWidgets.QDateEdit(self.groupBox_traj_period)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.traj_year_start.sizePolicy().hasHeightForWidth())
        self.traj_year_start.setSizePolicy(sizePolicy)
        self.traj_year_start.setMinimumSize(QtCore.QSize(0, 40))
        self.traj_year_start.setMaximumSize(QtCore.QSize(120, 35))
        self.traj_year_start.setMaximumDate(QtCore.QDate(2030, 12, 31))
        self.traj_year_start.setMinimumDate(QtCore.QDate(1980, 1, 1))
        self.traj_year_start.setDisplayFormat("yyyy")
        self.traj_year_start.setDate(QtCore.QDate(2001, 1, 1))
        self.traj_year_start.setObjectName("traj_year_start")
        self.gridLayout_5.addWidget(self.traj_year_start, 1, 0, 1, 1)
        self.label_7 = QtWidgets.QLabel(self.groupBox_traj_period)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.gridLayout_5.addWidget(self.label_7, 0, 0, 1, 1)
        self.label_6 = QtWidgets.QLabel(self.groupBox_traj_period)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.gridLayout_5.addWidget(self.label_6, 0, 1, 1, 1)
        self.traj_year_end = QtWidgets.QDateEdit(self.groupBox_traj_period)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.traj_year_end.sizePolicy().hasHeightForWidth())
        self.traj_year_end.setSizePolicy(sizePolicy)
        self.traj_year_end.setMinimumSize(QtCore.QSize(0, 40))
        self.traj_year_end.setMaximumSize(QtCore.QSize(120, 35))
        self.traj_year_end.setMaximumDateTime(QtCore.QDateTime(QtCore.QDate(2030, 12, 31), QtCore.QTime(23, 59, 59)))
        self.traj_year_end.setMinimumDate(QtCore.QDate(1980, 1, 1))
        self.traj_year_end.setDisplayFormat("yyyy")
        self.traj_year_end.setDate(QtCore.QDate(2015, 12, 31))
        self.traj_year_end.setObjectName("traj_year_end")
        self.gridLayout_5.addWidget(self.traj_year_end, 1, 1, 1, 1)
        self.verticalLayout.addWidget(self.groupBox_traj_period)
        self.groupBox_traj_climate = QtWidgets.QGroupBox(self.groupBox_traj)
        self.groupBox_traj_climate.setFlat(True)
        self.groupBox_traj_climate.setObjectName("groupBox_traj_climate")
        self.verticalLayout_8 = QtWidgets.QVBoxLayout(self.groupBox_traj_climate)
        self.verticalLayout_8.setObjectName("verticalLayout_8")
        self.traj_climate = QtWidgets.QComboBox(self.groupBox_traj_climate)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Preferred)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.traj_climate.sizePolicy().hasHeightForWidth())
        self.traj_climate.setSizePolicy(sizePolicy)
        self.traj_climate.setMinimumSize(QtCore.QSize(0, 40))
        self.traj_climate.setObjectName("traj_climate")
        self.verticalLayout_8.addWidget(self.traj_climate)
        self.verticalLayout.addWidget(self.groupBox_traj_climate)
        self.verticalLayout_9.addWidget(self.groupBox_traj)
        self.groupBox_perf = QtWidgets.QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_perf.setCheckable(True)
        self.groupBox_perf.setObjectName("groupBox_perf")
        self.gridLayout_8 = QtWidgets.QGridLayout(self.groupBox_perf)
        self.gridLayout_8.setObjectName("gridLayout_8")
        self.perf_year_start = QtWidgets.QDateEdit(self.groupBox_perf)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.perf_year_start.sizePolicy().hasHeightForWidth())
        self.perf_year_start.setSizePolicy(sizePolicy)
        self.perf_year_start.setMinimumSize(QtCore.QSize(0, 40))
        self.perf_year_start.setMaximumSize(QtCore.QSize(120, 35))
        self.perf_year_start.setMaximumDate(QtCore.QDate(2030, 12, 31))
        self.perf_year_start.setMinimumDate(QtCore.QDate(1980, 1, 1))
        self.perf_year_start.setDisplayFormat("yyyy")
        self.perf_year_start.setDate(QtCore.QDate(2001, 1, 1))
        self.perf_year_start.setObjectName("perf_year_start")
        self.gridLayout_8.addWidget(self.perf_year_start, 1, 0, 1, 1)
        self.label_10 = QtWidgets.QLabel(self.groupBox_perf)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_10.setFont(font)
        self.label_10.setObjectName("label_10")
        self.gridLayout_8.addWidget(self.label_10, 0, 0, 1, 1)
        self.label_14 = QtWidgets.QLabel(self.groupBox_perf)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_14.setFont(font)
        self.label_14.setObjectName("label_14")
        self.gridLayout_8.addWidget(self.label_14, 0, 1, 1, 1)
        self.perf_year_end = QtWidgets.QDateEdit(self.groupBox_perf)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.perf_year_end.sizePolicy().hasHeightForWidth())
        self.perf_year_end.setSizePolicy(sizePolicy)
        self.perf_year_end.setMinimumSize(QtCore.QSize(0, 40))
        self.perf_year_end.setMaximumSize(QtCore.QSize(120, 35))
        self.perf_year_end.setMaximumDateTime(QtCore.QDateTime(QtCore.QDate(2030, 12, 31), QtCore.QTime(23, 59, 59)))
        self.perf_year_end.setMinimumDate(QtCore.QDate(1980, 1, 1))
        self.perf_year_end.setDisplayFormat("yyyy")
        self.perf_year_end.setDate(QtCore.QDate(2015, 12, 31))
        self.perf_year_end.setObjectName("perf_year_end")
        self.gridLayout_8.addWidget(self.perf_year_end, 1, 1, 1, 1)
        self.verticalLayout_9.addWidget(self.groupBox_perf)
        self.groupBox_state = QtWidgets.QGroupBox(self.scrollAreaWidgetContents)
        self.groupBox_state.setCheckable(True)
        self.groupBox_state.setObjectName("groupBox_state")
        self.verticalLayout_5 = QtWidgets.QVBoxLayout(self.groupBox_state)
        self.verticalLayout_5.setObjectName("verticalLayout_5")
        self.groupBox_state_baseline = QtWidgets.QGroupBox(self.groupBox_state)
        self.groupBox_state_baseline.setFlat(True)
        self.groupBox_state_baseline.setObjectName("groupBox_state_baseline")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_state_baseline)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.state_year_bl_start = QtWidgets.QDateEdit(self.groupBox_state_baseline)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.state_year_bl_start.sizePolicy().hasHeightForWidth())
        self.state_year_bl_start.setSizePolicy(sizePolicy)
        self.state_year_bl_start.setMinimumSize(QtCore.QSize(0, 40))
        self.state_year_bl_start.setMaximumSize(QtCore.QSize(120, 35))
        self.state_year_bl_start.setMaximumDate(QtCore.QDate(2030, 12, 31))
        self.state_year_bl_start.setMinimumDate(QtCore.QDate(1980, 1, 1))
        self.state_year_bl_start.setDisplayFormat("yyyy")
        self.state_year_bl_start.setDate(QtCore.QDate(2001, 1, 1))
        self.state_year_bl_start.setObjectName("state_year_bl_start")
        self.gridLayout_3.addWidget(self.state_year_bl_start, 1, 0, 1, 1)
        self.label_8 = QtWidgets.QLabel(self.groupBox_state_baseline)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.gridLayout_3.addWidget(self.label_8, 0, 0, 1, 1)
        self.state_year_bl_end = QtWidgets.QDateEdit(self.groupBox_state_baseline)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.state_year_bl_end.sizePolicy().hasHeightForWidth())
        self.state_year_bl_end.setSizePolicy(sizePolicy)
        self.state_year_bl_end.setMinimumSize(QtCore.QSize(0, 40))
        self.state_year_bl_end.setMaximumSize(QtCore.QSize(120, 35))
        self.state_year_bl_end.setMaximumDateTime(QtCore.QDateTime(QtCore.QDate(2030, 12, 31), QtCore.QTime(23, 59, 59)))
        self.state_year_bl_end.setMinimumDateTime(QtCore.QDateTime(QtCore.QDate(1980, 1, 1), QtCore.QTime(0, 0, 0)))
        self.state_year_bl_end.setDisplayFormat("yyyy")
        self.state_year_bl_end.setDate(QtCore.QDate(2012, 12, 31))
        self.state_year_bl_end.setObjectName("state_year_bl_end")
        self.gridLayout_3.addWidget(self.state_year_bl_end, 1, 1, 1, 1)
        self.label_12 = QtWidgets.QLabel(self.groupBox_state_baseline)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_12.setFont(font)
        self.label_12.setObjectName("label_12")
        self.gridLayout_3.addWidget(self.label_12, 0, 1, 1, 1)
        self.verticalLayout_5.addWidget(self.groupBox_state_baseline)
        self.groupBox_state_comparison = QtWidgets.QGroupBox(self.groupBox_state)
        self.groupBox_state_comparison.setFlat(True)
        self.groupBox_state_comparison.setObjectName("groupBox_state_comparison")
        self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_state_comparison)
        self.gridLayout_4.setObjectName("gridLayout_4")
        self.label_16 = QtWidgets.QLabel(self.groupBox_state_comparison)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_16.setFont(font)
        self.label_16.setObjectName("label_16")
        self.gridLayout_4.addWidget(self.label_16, 3, 1, 1, 1)
        self.state_year_tg_start = QtWidgets.QDateEdit(self.groupBox_state_comparison)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.state_year_tg_start.sizePolicy().hasHeightForWidth())
        self.state_year_tg_start.setSizePolicy(sizePolicy)
        self.state_year_tg_start.setMinimumSize(QtCore.QSize(0, 40))
        self.state_year_tg_start.setMaximumSize(QtCore.QSize(120, 35))
        self.state_year_tg_start.setMaximumDateTime(QtCore.QDateTime(QtCore.QDate(2015, 12, 31), QtCore.QTime(23, 59, 59)))
        self.state_year_tg_start.setMinimumDateTime(QtCore.QDateTime(QtCore.QDate(2011, 1, 1), QtCore.QTime(0, 0, 0)))
        self.state_year_tg_start.setMaximumDate(QtCore.QDate(2015, 12, 31))
        self.state_year_tg_start.setMinimumDate(QtCore.QDate(2011, 1, 1))
        self.state_year_tg_start.setDisplayFormat("yyyy")
        self.state_year_tg_start.setDate(QtCore.QDate(2013, 1, 1))
        self.state_year_tg_start.setObjectName("state_year_tg_start")
        self.gridLayout_4.addWidget(self.state_year_tg_start, 4, 0, 1, 1)
        self.state_year_tg_end = QtWidgets.QDateEdit(self.groupBox_state_comparison)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Fixed, QtWidgets.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.state_year_tg_end.sizePolicy().hasHeightForWidth())
        self.state_year_tg_end.setSizePolicy(sizePolicy)
        self.state_year_tg_end.setMinimumSize(QtCore.QSize(0, 40))
        self.state_year_tg_end.setMaximumSize(QtCore.QSize(120, 35))
        self.state_year_tg_end.setMaximumDateTime(QtCore.QDateTime(QtCore.QDate(2030, 12, 31), QtCore.QTime(23, 59, 59)))
        self.state_year_tg_end.setMinimumDate(QtCore.QDate(1980, 1, 1))
        self.state_year_tg_end.setDisplayFormat("yyyy")
        self.state_year_tg_end.setDate(QtCore.QDate(2015, 12, 31))
        self.state_year_tg_end.setObjectName("state_year_tg_end")
        self.gridLayout_4.addWidget(self.state_year_tg_end, 4, 1, 1, 1)
        self.label_15 = QtWidgets.QLabel(self.groupBox_state_comparison)
        font = QtGui.QFont()
        font.setBold(False)
        font.setWeight(50)
        self.label_15.setFont(font)
        self.label_15.setObjectName("label_15")
        self.gridLayout_4.addWidget(self.label_15, 3, 0, 1, 1)
        self.verticalLayout_5.addWidget(self.groupBox_state_comparison)
        spacerItem1 = QtWidgets.QSpacerItem(20, 40, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.verticalLayout_5.addItem(spacerItem1)
        self.verticalLayout_9.addWidget(self.groupBox_state)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.verticalLayout_7.addWidget(self.scrollArea)
        self.TabBox.addTab(self.AdvancedTab, "")
        self.verticalLayout_2.addWidget(self.TabBox)
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.button_calculate = QtWidgets.QPushButton(DlgCalculateProd)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_calculate.sizePolicy().hasHeightForWidth())
        self.button_calculate.setSizePolicy(sizePolicy)
        self.button_calculate.setMinimumSize(QtCore.QSize(0, 40))
        self.button_calculate.setObjectName("button_calculate")
        self.gridLayout.addWidget(self.button_calculate, 1, 0, 1, 2)
        self.button_prev = QtWidgets.QPushButton(DlgCalculateProd)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_prev.sizePolicy().hasHeightForWidth())
        self.button_prev.setSizePolicy(sizePolicy)
        self.button_prev.setMinimumSize(QtCore.QSize(0, 40))
        self.button_prev.setObjectName("button_prev")
        self.gridLayout.addWidget(self.button_prev, 0, 0, 1, 1)
        self.button_next = QtWidgets.QPushButton(DlgCalculateProd)
        sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.button_next.sizePolicy().hasHeightForWidth())
        self.button_next.setSizePolicy(sizePolicy)
        self.button_next.setMinimumSize(QtCore.QSize(0, 40))
        self.button_next.setObjectName("button_next")
        self.gridLayout.addWidget(self.button_next, 0, 1, 1, 1)
        self.verticalLayout_2.addLayout(self.gridLayout)

        self.retranslateUi(DlgCalculateProd)
        self.TabBox.setCurrentIndex(1)
        QtCore.QMetaObject.connectSlotsByName(DlgCalculateProd)
        DlgCalculateProd.setTabOrder(self.TabBox, self.button_prev)
        DlgCalculateProd.setTabOrder(self.button_prev, self.button_next)
        DlgCalculateProd.setTabOrder(self.button_next, self.button_calculate)
        DlgCalculateProd.setTabOrder(self.button_calculate, self.traj_indic)
        DlgCalculateProd.setTabOrder(self.traj_indic, self.traj_climate)

    def retranslateUi(self, DlgCalculateProd):
        _translate = QtCore.QCoreApplication.translate
        DlgCalculateProd.setWindowTitle(_translate("DlgCalculateProd", "Calculate Productivity"))
        self.groupBox_lpd_mode.setTitle(_translate("DlgCalculateProd", "Indicators to calculate"))
        self.label.setText(_translate("DlgCalculateProd", "UNCCD default data (Land Productivity Dynamics (LPD) from Joint Research Commission)"))
        self.mode_te_prod.setText(_translate("DlgCalculateProd", "MISLAND land productivity"))
        self.groupBox_ndvi_dataset.setTitle(_translate("DlgCalculateProd", "NDVI dataset"))
        self.TabBox.setTabText(self.TabBox.indexOf(self.IndicatorsTab), _translate("DlgCalculateProd", "Indicators"))
        self.groupBox_traj.setTitle(_translate("DlgCalculateProd", "Trajectory (related to rate of change of  productivity over time)"))
        self.groupBox_4.setTitle(_translate("DlgCalculateProd", "Trajectory indicator"))
        self.groupBox_traj_period.setTitle(_translate("DlgCalculateProd", "Period"))
        self.label_7.setText(_translate("DlgCalculateProd", "Starting year:"))
        self.label_6.setText(_translate("DlgCalculateProd", "Ending year:"))
        self.groupBox_traj_climate.setTitle(_translate("DlgCalculateProd", "Climate dataset"))
        self.groupBox_perf.setTitle(_translate("DlgCalculateProd", "Performance (a measure of how productivity in an area compares to that of similar areas)"))
        self.label_10.setText(_translate("DlgCalculateProd", "Starting year:"))
        self.label_14.setText(_translate("DlgCalculateProd", "Ending year:"))
        self.groupBox_state.setTitle(_translate("DlgCalculateProd", "State (compares current productivity in an area to past productivity in the same area)"))
        self.groupBox_state_baseline.setTitle(_translate("DlgCalculateProd", "Initial period"))
        self.label_8.setText(_translate("DlgCalculateProd", "Starting year:"))
        self.label_12.setText(_translate("DlgCalculateProd", "Ending year:"))
        self.groupBox_state_comparison.setTitle(_translate("DlgCalculateProd", "Comparison period"))
        self.label_16.setText(_translate("DlgCalculateProd", "Ending year:"))
        self.label_15.setText(_translate("DlgCalculateProd", "Starting year:"))
        self.TabBox.setTabText(self.TabBox.indexOf(self.AdvancedTab), _translate("DlgCalculateProd", "Advanced"))
        self.button_calculate.setText(_translate("DlgCalculateProd", "Calculate"))
        self.button_prev.setText(_translate("DlgCalculateProd", "Previous"))
        self.button_next.setText(_translate("DlgCalculateProd", "Next"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    DlgCalculateProd = QtWidgets.QDialog()
    ui = Ui_DlgCalculateProd()
    ui.setupUi(DlgCalculateProd)
    DlgCalculateProd.show()
    sys.exit(app.exec_())