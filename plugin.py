# -*- coding: utf-8 -*-
"""
/***************************************************************************
 MISLAND - A QGIS plugin
 This plugin supports monitoring and reporting of land degradation to the UNCCD 
 and in support of the SDG Land Degradation Neutrality (LDN) target.
                              -------------------
        begin                : 2017-05-23
        git sha              : $Format:%H$
        copyright            : (C) 2017 by Conservation International
        email                : trends.earth@conservation.org
 ***************************************************************************/
"""

from builtins import object
import os
import subprocess
import sys
ext_libs_path = os.path.join(os.path.dirname(__file__), 'ext_libs')
if ext_libs_path not in sys.path:
    sys.path.append(ext_libs_path)
from qgis.PyQt.QtCore import QSettings, QTranslator, qVersion, QCoreApplication
from qgis.PyQt.QtWidgets import QAction, QMessageBox, QApplication, QMenu
from qgis.PyQt.QtGui import QIcon
from qgis.core import QgsApplication, QgsMessageLog, Qgis
from qgis.utils import showPluginHelp
try:
    from MISLAND import __version__
    from MISLAND.settings import DlgSettings
    from MISLAND.download_data import DlgDownload
    from MISLAND.calculate import DlgCalculate
    from MISLAND.jobs import DlgJobs
    from MISLAND.timeseries import DlgTimeseries
    from MISLAND.visualization import DlgVisualization
    from MISLAND.data_io import DlgDataIO
    from MISLAND.about import DlgAbout
    from MISLAND.processing_provider.provider import Provider
    # Initialize Qt resources from file resources.py
    import MISLAND.resources

except (ImportError, ModuleNotFoundError,FileNotFoundError):
    plugin_path = os.path.dirname(__file__)
    req_path = os.path.join(plugin_path,"requirements.txt").replace("\\","/")
    target_path = os.path.join(plugin_path,"ext-libs").replace("\\","/")
    print(req_path)
    print(target_path)
    print("Installing python paths")

    # Find the full path to the pip executable
    py_version = sys.version.split(".")
    py_folder = "Python"+py_version[0]+py_version[1]
    pip_executable =""
    if os.name == "nt":
        pip_executable = os.path.join(os.path.dirname(os.path.dirname(sys.executable)), 'apps',py_folder,"Scripts", 'pip3.exe')
    elif os.name == "posix":
        pip_executable = os.path.join(os.path.dirname(os.path.dirname(sys.executable)), 'bin','pip3')
    # pip_executable=pip_executable.replace("\\bin","")
    print(pip_executable)

    try:
        command = [pip_executable,"install","--upgrade","-r",req_path,"--target",target_path]
        print(command)
        subprocess.run(command)
    except subprocess.CalledProcessError as e:
        print(e)

    except FileNotFoundError as e:
        print(e)

    finally:
        from MISLAND import __version__
        from MISLAND.settings import DlgSettings
        from MISLAND.download_data import DlgDownload
        from MISLAND.calculate import DlgCalculate
        from MISLAND.jobs import DlgJobs
        from MISLAND.timeseries import DlgTimeseries
        from MISLAND.visualization import DlgVisualization
        from MISLAND.data_io import DlgDataIO
        from MISLAND.about import DlgAbout
        from MISLAND.processing_provider.provider import Provider
        # Initialize Qt resources from file resources.py
        import MISLAND.resources
    

def showHelp(file='index', section=None):
    locale = QSettings().value('locale/userLocale')[0:2]
    help_base_path = os.path.join(os.path.dirname(__file__), 'help', 'build', 'html')
    locale_path = os.path.join(help_base_path, locale)
    QgsMessageLog.logMessage(u'Checking for plugin help in {}'.format(locale_path), tag="MISLAND", level=Qgis.Info)
    if os.path.exists(locale_path):
        help_path = os.path.join(locale_path, file)
    else:
        help_path = os.path.join(help_base_path, 'en', file)
    QgsMessageLog.logMessage(u'Showing plugin help from {}'.format(help_path), tag="MISLAND", level=Qgis.Info)
    if section:
        showPluginHelp(filename=help_path, section=section)
    else:
        showPluginHelp(filename=help_path)


class LDMPPlugin(object):
    """QGIS Plugin Implementation."""

    def __init__(self, iface):
        """Constructor.

        :param iface: An interface instance that will be passed to this class
            which provides the hook by which you can manipulate the QGIS
            application at run time.
        :type iface: QgsInterface
        """
        # Save reference to the QGIS interface
        self.iface = iface
        # initialize plugin directory
        self.plugin_dir = os.path.dirname(__file__)

        self.provider = None

        # initialize locale and translation
        locale = QSettings().value('locale/userLocale')[0:2]
        locale_path = os.path.join(
            self.plugin_dir,
            'i18n',
            u'LDMP_{}.qm'.format(locale))
        QgsMessageLog.logMessage(u'Starting MISLAND version {} using locale "{}" in path {}.'.format(__version__, locale, locale_path),
                                 tag="MISLAND", level=Qgis.Info)

        if os.path.exists(locale_path):
            self.translator = QTranslator()
            self.translator.load(locale_path)
            if qVersion() > '4.3.3':
                QCoreApplication.installTranslator(self.translator)
                QgsMessageLog.logMessage("Translator installed.", tag="MISLAND",
                                         level=Qgis.Info)

        # Declare instance attributes
        self.actions = []
        self.menu = QMenu(QApplication.translate('MISLAND', u'&MISLAND'))
        self.menu.setIcon(QIcon(':/plugins/MISLAND/misland_logo.png'))
        self.raster_menu = self.iface.rasterMenu()
        self.raster_menu.addMenu(self.menu)
        self.toolbar = self.iface.addToolBar(u'MISLAND')
        try:
            self.dlg_settings = DlgSettings()
            self.dlg_calculate = DlgCalculate()
            self.dlg_jobs = DlgJobs()
            self.dlg_timeseries = DlgTimeseries()
            self.dlg_visualization = DlgVisualization()
            self.dlg_download = DlgDownload()
            self.dlg_data_io = DlgDataIO()
            self.dlg_about = DlgAbout()
        except (ImportError,ModuleNotFoundError,FileNotFoundError):
            plugin_path = os.path.dirname(__file__)
            req_path = os.path.join(plugin_path,"requirements.txt")
            target_path = os.path.join(plugin_path,"ext-libs")
            pip_executable = os.path.join(os.path.dirname(os.path.dirname(sys.executable)), 'apps',py_folder,"Scripts", 'pip3.exe') if os.name == 'nt' else 'pip'
            # print("Installing python paths")
            try:
                subprocess.run([pip_executable,"install","--upgrade","-r",req_path,"--target",target_path])
            except subprocess.CalledProcessError as e:
                print(e)

    # noinspection PyMethodMayBeStatic
    def tr(self, message):
        """Get the translation for a string using Qt translation API.

        We implement this ourselves since we do not inherit QObject.

        :param message: String for translation.
        :type message: str, QString

        :returns: Translated version of message.
        :rtype: QString
        """
        # noinspection PyTypeChecker,PyArgumentList,PyCallByClass
        return QCoreApplication.translate('MISLAND', message)

    def initProcessing(self):
        self.provider = Provider()
        QgsApplication.processingRegistry().addProvider(self.provider)
        

    def add_action(
            self,
            icon_path,
            text,
            callback,
            enabled_flag=True,
            add_to_menu=True,
            add_to_toolbar=True,
            status_tip=None,
            whats_this=None,
            parent=None):
        """Add a toolbar icon to the toolbar.

        :param icon_path: Path to the icon for this action. Can be a resource
            path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
        :type icon_path: str

        :param text: Text that should be shown in menu items for this action.
        :type text: str

        :param callback: Function to be called when the action is triggered.
        :type callback: function

        :param enabled_flag: A flag indicating if the action should be enabled
            by default. Defaults to True.
        :type enabled_flag: bool

        :param add_to_menu: Flag indicating whether the action should also
            be added to the menu. Defaults to True.
        :type add_to_menu: bool

        :param add_to_toolbar: Flag indicating whether the action should also
            be added to the toolbar. Defaults to True.
        :type add_to_toolbar: bool

        :param status_tip: Optional text to show in a popup when mouse pointer
            hovers over the action.
        :type status_tip: str

        :param parent: Parent widget for the new action. Defaults None.
        :type parent: QWidget

        :param whats_this: Optional text to show in the status bar when the
            mouse pointer hovers over the action.

        :returns: The action that was created. Note that the action is also
            added to self.actions list.
        :rtype: QAction
        """

        icon = QIcon(icon_path)
        action = QAction(icon, text, parent)
        action.triggered.connect(callback)
        action.setEnabled(enabled_flag)

        if status_tip is not None:
            action.setStatusTip(status_tip)

        if whats_this is not None:
            action.setWhatsThis(whats_this)

        if add_to_toolbar:
            self.toolbar.addAction(action)

        if add_to_menu:
            self.menu.addAction(action)

        self.actions.append(action)

        return action

    def initGui(self):
        self.initProcessing()

        """Create the menu entries and toolbar icons inside the QGIS GUI."""
        self.add_action(
            ':/plugins/MISLAND/icons/wrench2.svg',
            text=QApplication.translate('MISLAND', u'User Settings'),
            callback=self.run_settings,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'LDMT User Settings'))

        self.add_action(
            ':/plugins/MISLAND/icons/calculator1.svg',
            text=QApplication.translate('MISLAND', u'Calculate indicators'),
            callback=self.run_calculate,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'Calculate indicators'))

        # self.add_action(
        #     ':/plugins/MISLAND/icons/graph1.svg',
        #     text=QApplication.translate('MISLAND', u'Plot data'),
        #     callback=self.run_plot,
        #     parent=self.iface.mainWindow(),
        #     status_tip=QApplication.translate('MISLAND', 'Plot time series datasets'))

        self.add_action(
            ':/plugins/MISLAND/icons/folder1.svg',
            text=QApplication.translate('MISLAND', u'Load data'),
            callback=self.data_io,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'Load local data'))

        self.add_action(
            ':/plugins/MISLAND/icons/globe1.svg',
            text=QApplication.translate('MISLAND', u'Download raw data'),
            callback=self.run_download,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'Download raw datasets'))

        self.add_action(
            ':/plugins/MISLAND/icons/document1.svg',
            text=QApplication.translate('MISLAND', u'Visualization tool'),
            callback=self.run_visualization,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'Visualize and summarize data'))
            
                    
        self.add_action(
            ':/plugins/MISLAND/icons/cloud-download1.svg',
            text=QApplication.translate('MISLAND', u'View Google Earth Engine tasks'),
            callback=self.get_jobs,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'View cloud processing tasks'))


        self.add_action(
            ':/plugins/MISLAND/icons/info1.svg',
            text=QApplication.translate('MISLAND', u'About'),
            callback=self.run_about,
            parent=self.iface.mainWindow(),
            status_tip=QApplication.translate('MISLAND', 'About MISLAND'))

    def unload(self):
        self.provider = Provider()
        """Removes the plugin menu item and icon from QGIS GUI."""
        for action in self.actions:
            self.iface.removePluginRasterMenu(
                QApplication.translate('MISLAND', u'&MISLAND'),
                action)
            self.iface.removeToolBarIcon(action)
        # remove the menu
        self.raster_menu.removeAction(self.menu.menuAction())
        # remove the toolbar
        del self.toolbar

        QgsApplication.processingRegistry().removeProvider(self.provider)

    def run_settings(self):
        self.dlg_settings.show()
        result = self.dlg_settings.exec_()

    def run_download(self):
        self.dlg_download.show()
        result = self.dlg_download.exec_()

    def run_calculate(self):
        # show the dialog
        self.dlg_calculate.show()
        result = self.dlg_calculate.exec_()

    def get_jobs(self):
        # show the dialog
        self.dlg_jobs.show()
        result = self.dlg_jobs.exec_()

    def run_plot(self):
        self.dlg_timeseries.show()
        result = self.dlg_timeseries.exec_()

    def run_visualization(self):
        self.dlg_visualization.show()
        result = self.dlg_visualization.exec_()

    def data_io(self):
        self.dlg_data_io.show()
        result = self.dlg_data_io.exec_()

    def run_about(self):
        #showHelp()
        self.dlg_about.show()
        result = self.dlg_about.exec_()
