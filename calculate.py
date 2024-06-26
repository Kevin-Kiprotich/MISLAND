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
from pathlib import Path
import json
import tempfile

from osgeo import gdal, ogr, osr

from qgis.PyQt import QtWidgets
from qgis.PyQt.QtGui import QIcon, QPixmap, QDoubleValidator,QDesktopServices
from qgis.PyQt.QtCore import (QTextCodec, QSettings, pyqtSignal,
    QCoreApplication, QDate,QUrl,Qt)

from qgis.core import (QgsFeature, QgsPointXY, QgsGeometry, QgsJsonUtils,
    QgsVectorLayer, QgsCoordinateTransform, QgsCoordinateReferenceSystem,
    Qgis, QgsProject, QgsLayerTreeGroup, QgsLayerTreeLayer,
    QgsVectorFileWriter, QgsFields, QgsWkbTypes)
from qgis.utils import iface
from qgis.gui import QgsMapToolEmitPoint, QgsMapToolPan

from MISLAND import log, GetTempFilename
from MISLAND.api import run_script
# from MISLAND.timeseries import DlgTimeseries
from MISLAND.gui.DlgCalculate import Ui_DlgCalculate
from MISLAND.gui.DlgCalculateLD import Ui_DlgCalculateLD
from MISLAND.gui.DlgCalculateTC import Ui_DlgCalculateTC
from MISLAND.gui.DlgCalculateRestBiomass import Ui_DlgCalculateRestBiomass
from MISLAND.gui.DlgCalculateUrban import Ui_DlgCalculateUrban
from MISLAND.gui.DlgCalculateMedalus import Ui_DlgCalculateMedalus
from MISLAND.gui.DlgCalculateForest import Ui_DlgCalculateForest
from MISLAND.gui.DlgTimeseries import Ui_DlgTimeseries
from MISLAND.gui.DlgCalculateForestFire import Ui_DlgCalculateForestFire

#Add the coastal erosion ui
from MISLAND.gui.DlgCalculateCE import Ui_DlgCalculateCE

#Add the soil erosion ui
from MISLAND.gui.DlgCalculateSE import Ui_DlgCalculateSE


from MISLAND.gui.WidgetSelectArea import Ui_WidgetSelectArea
from MISLAND.gui.WidgetCalculationOptions import Ui_WidgetCalculationOptions
from MISLAND.gui.WidgetCalculationOutput import Ui_WidgetCalculationOutput
from MISLAND.download import read_json, get_admin_bounds, get_cities
from MISLAND.worker import AbstractWorker

from MISLAND.misland_africa.downloads import get_Adm0, get_Adm1, get_Adm2,downloadGeomorphology,downloadCVI
from MISLAND.misland_africa.water_erosion import getRUSLE
from MISLAND.misland_africa.wind_erosion import getILSWE

mb = iface.messageBar()


if QSettings().value("MISLAND/binaries_enabled", False) == 'True':
    try:
        from trends_earth_binaries.calculate_numba import *
        log("Using numba-compiled version of calculate_numba.")
    except (ModuleNotFoundError, ImportError) as e:
        from MISLAND.calculate_numba import *
        log("Failed import of numba-compiled code, falling back to python version of calculate_numba.")
else:
    from MISLAND.calculate_numba import *
    log("Using python version of calculate_numba.")


class tr_calculate(object):
    def tr(message):
        return QCoreApplication.translate("tr_calculate", message)


# Make a function to get a script slug from a script name, including the script 
# version string
with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                       'data', 'scripts.json')) as script_file:
    scripts = json.load(script_file)
def get_script_slug(script_name):
    # Note that dots and underscores can't be used in the slugs, so they are 
    # replaced with dashesk
    return script_name + '-' + scripts[script_name]['script version'].replace('.', '-')

# Transform CRS of a layer while optionally wrapping geometries
# across the 180th meridian
def transform_layer(l, crs_dst, datatype='polygon', wrap=False):
    log('Transforming layer from "{}" to "{}". Wrap is {}. Datatype is {}.'.format(l.crs().toProj(), crs_dst.toProj(), wrap, datatype))

    crs_src_string = l.crs().toProj()
    if wrap:
        if not l.crs().isGeographic():
            QtWidgets.QMessageBox.critical(None,tr_calculate.tr("Error"),
                   tr_calculate.tr("Error - layer is not in a geographic coordinate system. Cannot wrap layer across 180th meridian."))
            log('Can\'t wrap layer in non-geographic coordinate system: "{}"'.format(crs_src_string))
            return None
        crs_src_string = crs_src_string + ' +lon_wrap=180'
    crs_src = QgsCoordinateReferenceSystem()
    crs_src.createFromProj(crs_src_string)
    t = QgsCoordinateTransform(crs_src, crs_dst, QgsProject.instance())

    l_w = QgsVectorLayer("{datatype}?crs=proj4:{crs}".format(datatype=datatype, 
                         crs=crs_dst.toProj()), "calculation boundary (transformed)",  
                         "memory")
    feats = []
    for f in l.getFeatures():
        geom = f.geometry()
        if wrap:
            n = 0
            p = geom.vertexAt(n)
            # Note vertexAt returns QgsPointXY(0, 0) on error
            while p != QgsPointXY(0, 0):
                if p.x() < 0:
                    geom.moveVertex(p.x() + 360, p.y(), n)
                n += 1
                p = geom.vertexAt(n)
        geom.transform(t)
        f.setGeometry(geom)
        feats.append(f)
    l_w.dataProvider().addFeatures(feats)
    l_w.commitChanges()
    if not l_w.isValid():
        log('Error transforming layer from "{}" to "{}" (wrap is {})'.format(crs_src_string, crs_dst.toProj(), wrap))
        return None
    else:
        return l_w


def json_geom_to_geojson(txt):
    d = {'type': 'FeatureCollection',
         'features': [{'type': 'Feature',
                       'geometry': json.loads(txt)}]
         }
    return d


class AOI(object):
    def __init__(self, crs_dst):
        self.crs_dst = crs_dst

    def get_crs_dst_wkt(self):
        return self.crs_dst.toWkt()

    def update_from_file(self, f, wrap=False):
        log(u'Setting up AOI from file at {}"'.format(f))
        l = QgsVectorLayer(f, "calculation boundary", "ogr")
        if not l.isValid():
            QtWidgets.QMessageBox.critical(None,
                    tr_calculate.tr("Error"),
                    tr_calculate.tr(u"Unable to load area of interest from {}. There may be a problem with the file or coordinate system. Try manually loading this file into QGIS to verify that it displays properly. If you continue to have problems with this file, send us a message at trends.earth@conservation.org.".format(f)))
            log("Unable to load area of interest.")
            return
        if l.wkbType() == QgsWkbTypes.Polygon \
                or l.wkbType() == QgsWkbTypes.PolygonZ \
                or l.wkbType() == QgsWkbTypes.MultiPolygon \
                or l.wkbType() == QgsWkbTypes.MultiPolygonZ:
            self.datatype = "polygon"
        elif l.wkbType() == QgsWkbTypes.Point \
                or l.wkbType() == QgsWkbTypes.PointZ \
                or l.wkbType() == QgsWkbTypes.MultiPoint \
                or l.wkbType() == QgsWkbTypes.MultiPointZ:
            self.datatype = "point"
        else:
            QtWidgets.QMessageBox.critical(None,
                    tr_calculate.tr("Error"),
                    tr_calculate.tr("Failed to process area of interest - unknown geometry type: {}".format(l.wkbType())))
            log("Failed to process area of interest - unknown geometry type.")
            self.datatype = "unknown"
            return

        self.l = transform_layer(l, self.crs_dst, datatype=self.datatype, wrap=wrap)

    def update_from_geojson(self, geojson, crs_src='epsg:4326', datatype='polygon', wrap=False):
        log('Setting up AOI with geojson. Wrap is {}.'.format(wrap))
        self.datatype = datatype
        # Note geojson is assumed to be in 4326
        l = QgsVectorLayer("{datatype}?crs={crs}".format(datatype=self.datatype, crs=crs_src), "calculation boundary", "memory")
        ds = ogr.Open(json.dumps(geojson))
        layer_in = ds.GetLayer()
        feats_out = []
        for i in range(0, layer_in.GetFeatureCount()):
            feat_in = layer_in.GetFeature(i)
            feat = QgsFeature(l.fields())
            geom = QgsGeometry()
            geom.fromWkb(feat_in.geometry().ExportToWkb())
            feat.setGeometry(geom)
            feats_out.append(feat)
        l.dataProvider().addFeatures(feats_out)
        l.commitChanges()
        if not l.isValid():
            QtWidgets.QMessageBox.critical(None,tr_calculate.tr("Error"),
                                      tr_calculate.tr("Failed to add geojson to temporary layer."))
            log("Failed to add geojson to temporary layer.")
            return
        self.l = transform_layer(l, self.crs_dst, datatype=self.datatype, wrap=wrap)

    def meridian_split(self, out_type='extent', out_format='geojson', warn=True):
        """
        Return list of bounding boxes in WGS84 as geojson for GEE

        Returns multiple geometries as needed to avoid having an extent 
        crossing the 180th meridian
        """

        if out_type not in ['extent', 'layer']:
            raise ValueError('Unrecognized out_type "{}"'.format(out_type))
        if out_format not in ['geojson', 'wkt']:
            raise ValueError(u'Unrecognized out_format "{}"'.format(out_format))

        # Calculate a single feature that is the union of all the features in 
        # this layer - that way there is a single feature to intersect with 
        # each hemisphere.
        geometries = []
        n = 1
        for f in self.get_layer_wgs84().getFeatures():
            # Get an OGR geometry from the QGIS geometry
            geom = f.geometry()
            if not geom.isGeosValid():
                log(u'Invalid feature in row {}.'.format(n))
                QtWidgets.QMessageBox.critical(None,
                                               tr_calculate.tr("Error"),
                                               tr_calculate.tr('Invalid geometry in row {}. Check that all input geometries are valid before processing. Try using the check validity tool on the "Vector" menu on the toolbar for more information on which features are invalid (Under "Vector" - "Geometry Tools" - "Check Validity").'.format(n)))
                return None
            geometries.append(geom)
            n += 1
        union = QgsGeometry.unaryUnion(geometries)

        log(u'Calculating east and west intersections to test if AOI crosses 180th meridian.')
        hemi_e = QgsGeometry.fromWkt('POLYGON ((0 -90, 0 90, 180 90, 180 -90, 0 -90))')
        hemi_w = QgsGeometry.fromWkt('POLYGON ((-180 -90, -180 90, 0 90, 0 -90, -180 -90))')
        intersections = [hemi.intersection(union) for hemi in [hemi_e, hemi_w]]

        if out_type == 'extent':
            pieces = [QgsGeometry.fromRect(i.boundingBox()) for i in intersections if not i.isEmpty()]
        elif out_type == 'layer':
            pieces = [i for i in intersections if not i.isEmpty()]
        pieces_union = QgsGeometry.unaryUnion(pieces)

        if out_format == 'geojson':
            pieces_txt = [json.loads(piece.asJson()) for piece in pieces]
            pieces_union_txt = json.loads(pieces_union.asJson())
        elif out_format == 'wkt':
            pieces_txt = [piece.asWkt() for piece in pieces]
            pieces_union_txt = pieces_union.asWkt()

        if (len(pieces) == 0) or (sum([piece.area() for piece in pieces]) > (pieces_union.area() / 2)):
            # If there is no area in one of the hemispheres, return the
            # original layer, or extent of the original layer. Also return the 
            # original layer (or extent) if the area of the combined pieces
            # from both hemispheres is not significantly smaller than that of 
            # the original polygon.
            log("AOI being processed in one piece (does not cross 180th meridian)")
            return (False, [pieces_union_txt])
        else:
            log("AOI crosses 180th meridian - splitting AOI into two geojsons.")
            if warn:
                QtWidgets.QMessageBox.information(None,tr_calculate.tr("Warning"),
                       tr_calculate.tr('The chosen area crosses the 180th meridian. It is recommended that you set the project coordinate system to a local coordinate system (see the "CRS" tab of the "Project Properties" window from the "Project" menu.)'))
            return (True, pieces_txt)

    def get_aligned_output_bounds(self, f):
        wkts = self.meridian_split(out_format='wkt', warn=False)[1]
        if not wkts:
            return None
        out = []
        for wkt in wkts:
            # Compute the pixel-aligned bounding box (slightly larger than 
            # aoi).
            # Use this to set bounds in vrt files in order to keep the
            # pixels aligned with the chosen layer
            geom = ogr.CreateGeometryFromWkt(wkt)
            (minx, maxx, miny, maxy) = geom.GetEnvelope()
            gt = gdal.Open(f).GetGeoTransform()
            left = minx - (minx - gt[0]) % gt[1]
            right = maxx + (gt[1] - ((maxx - gt[0]) % gt[1]))
            bottom = miny + (gt[5] - ((miny - gt[3]) % gt[5]))
            top = maxy - (maxy - gt[3]) % gt[5]
            out.append([left, bottom, right, top])
        return out

    def get_aligned_output_bounds_deprecated(self, f):
        # Compute the pixel-aligned bounding box (slightly larger than aoi).
        # Use this to set bounds in vrt files in order to keep the
        # pixels aligned with the chosen layer
        bb = self.bounding_box_geom().boundingBox()
        minx = bb.xMinimum()
        miny = bb.yMinimum()
        maxx = bb.xMaximum()
        maxy = bb.yMaximum()
        gt = gdal.Open(f).GetGeoTransform()
        left = minx - (minx - gt[0]) % gt[1]
        right = maxx + (gt[1] - ((maxx - gt[0]) % gt[1]))
        bottom = miny + (gt[5] - ((miny - gt[3]) % gt[5]))
        top = maxy - (maxy - gt[3]) % gt[5]
        return [left, bottom, right, top]

    def get_area(self):
        wgs84_crs = QgsCoordinateReferenceSystem('EPSG:4326')
        
        # Returns area of aoi components in sq m
        wkts = self.meridian_split(out_format='wkt', warn=False)[1]
        if not wkts:
            return None
        area = 0.
        for wkt in wkts:
            geom = QgsGeometry.fromWkt(wkt)
            # Lambert azimuthal equal area centered on polygon centroid
            centroid = geom.centroid().asPoint()
            laea_crs = QgsCoordinateReferenceSystem.fromProj('+proj=laea +lat_0={} +lon_0={}'.format(centroid.y(), centroid.x()))
            to_laea = QgsCoordinateTransform(wgs84_crs, laea_crs, QgsProject.instance())

            try:
                ret = geom.transform(to_laea)
            except:
                log('Error buffering layer while transforming to laea')
                QtWidgets.QMessageBox.critical(None,tr_calculate.tr("Error"),
                                          tr_calculate.tr("Error transforming coordinates. Check that the input geometry is valid."))
                return None
            this_area = geom.area()
            area += this_area
        return area

    def get_layer(self):
        """
        Return layer
        """
        return self.l

    def get_layer_wgs84(self):
        """
        Return layer in WGS84 (WPGS:4326)
        """
        # Setup settings for AOI provided to GEE:
        wgs84_crs = QgsCoordinateReferenceSystem()
        wgs84_crs.createFromProj('+proj=longlat +datum=WGS84 +no_defs')
        return transform_layer(self.l, wgs84_crs, datatype=self.datatype, wrap=False)

    def bounding_box_geom(self):
        'Returns bounding box in chosen destination coordinate system'
        return QgsGeometry.fromRect(self.l.extent())

    def bounding_box_gee_geojson(self):
        '''
        Returns two values - first is an indicator of whether this geojson 
        includes two geometries due to crossing of the 180th meridian, and the 
        second is the list of bounding box geojsons.
        '''
        if self.datatype == 'polygon':
            return self.meridian_split()
        elif self.datatype == 'point':
            # If there is only on point, don't calculate an extent (extent of 
            # one point is a box with sides equal to zero)
            n = 0
            for f in self.l.getFeatures():
                n += 1
                if n == 1:
                    # Save the first geometry in case it is needed later 
                    # for a layer that only has one point in it
                    geom = f.geometry()
            if n == 1:
                log('Layer only has one point')
                return (False, [json.loads(geom.asJson())])
            else:
                log('Layer has many points ({})'.format(n))
                return self.meridian_split()
        else:
            QtWidgets.QMessageBox.critical(None,tr_calculate.tr("Error"),
                   tr_calculate.tr("Failed to process area of interest - unknown geometry type: {}".format(self.datatype)))
            log("Failed to process area of interest - unknown geometry type.")


    def buffer(self, d):
        log('Buffering layer by {} km.'.format(d))

        feats = []
        for f in self.l.getFeatures():
            geom = f.geometry()
            # Setup an azimuthal equidistant projection centered on the polygon 
            # centroid
            centroid = geom.centroid().asPoint()
            geom.centroid()
            wgs84_crs = QgsCoordinateReferenceSystem('EPSG:4326')
            aeqd_crs = QgsCoordinateReferenceSystem.fromProj('+proj=aeqd +lat_0={} +lon_0={}'.format(centroid.y(), centroid.x()))
            to_aeqd = QgsCoordinateTransform(wgs84_crs, aeqd_crs, QgsProject.instance())

            try:
                ret = geom.transform(to_aeqd)
            except:
                log('Error buffering layer while transforming to aeqd')
                QtWidgets.QMessageBox.critical(None,tr_calculate.tr("Error"),
                                          tr_calculate.tr("Error transforming coordinates. Check that the input geometry is valid."))
                return None
            # Need to convert from km to meters
            geom_buffered = geom.buffer(d * 1000, 100)
            log('Feature area in sq km after buffering (and in aeqd) is: {}'.format(geom_buffered.area()/(1000 * 1000)))
            geom_buffered.transform(to_aeqd, QgsCoordinateTransform.TransformDirection.ReverseTransform)
            f.setGeometry(geom_buffered)
            feats.append(f)
            log('Feature area after buffering (and in WGS84) is: {}'.format(geom_buffered.area()))

        l_buffered = QgsVectorLayer("polygon?crs=proj4:{crs}".format(crs=self.l.crs().toProj()),
                                    "calculation boundary (transformed)",  
                                    "memory")
        l_buffered.dataProvider().addFeatures(feats)
        l_buffered.commitChanges()


        if not l_buffered.isValid():
            log('Error buffering layer')
            raise
        else:
            self.l = l_buffered
            self.datatype = 'polygon'
        return True

    def isValid(self):
        return self.l.isValid()

    def calc_frac_overlap(self, geom):
        """
        Returns fraction of AOI that is overlapped by geom (where geom is a QgsGeometry)

        Used to calculate "within" with a tolerance
        """
        aoi_geom = ogr.CreateGeometryFromWkt(self.bounding_box_geom().asWkt())
        in_geom = ogr.CreateGeometryFromWkt(geom.asWkt())

        geom_area = aoi_geom.GetArea()
        if geom_area == 0:
            # Handle case of a point with zero area
            frac = aoi_geom.Within(in_geom)
        else:
            frac = aoi_geom.Intersection(in_geom).GetArea() / geom_area
            log('Fractional area of overlap: {}'.format(frac))
        return frac


class DlgCalculate(QtWidgets.QDialog, Ui_DlgCalculate):
    def __init__(self, parent=None):
        super(DlgCalculate, self).__init__(parent)

        self.setupUi(self)

        self.dlg_calculate_ld = DlgCalculateLD()
        self.dlg_calculate_tc = DlgCalculateTC()
        self.dlg_calculate_rest_biomass = DlgCalculateRestBiomass()
        # self.dlg_calculate_urban = DlgCalculateUrban()
        self.dlg_calculate_forest = DlgCalculateForest()
        self.dlg_calculate_medalus = DlgCalculateMedalus()

        # add the calculate cvi sequence
        self.dlg_calculate_ce = DlgCalculateCE()

        # add the soil erosion sequence
        self.dlg_calculate_se = DlgCalculateSE() 

        from MISLAND.timeseries import DlgTimeseries
        self.dlg_timeseries = DlgTimeseries()
        
        self.pushButton_ld.clicked.connect(self.btn_ld_clicked)
        self.pushButton_timeseries.clicked.connect(self.btn_timeseries_clicked)
        # self.pushButton_tc.clicked.connect(self.btn_tc_clicked)
        # self.pushButton_rest_biomass.clicked.connect(self.btn_rest_biomass_clicked)
        # self.pushButton_urban.clicked.connect(self.btn_urban_clicked)
        self.pushButton_forest.clicked.connect(self.btn_forest_clicked)
        self.pushButton_medalus.clicked.connect(self.btn_medalus_clicked)
        self.pushbutton_CE.clicked.connect(self.btn_ce_clicked)
        self.pushButton_SE.clicked.connect(self.btn_se_clicked)
    
    def btn_ld_clicked(self):
        self.close()
        result = self.dlg_calculate_ld.exec_()

    def btn_tc_clicked(self):
        self.close()
        result = self.dlg_calculate_tc.exec_()

    def btn_rest_biomass_clicked(self):
        self.close()
        result = self.dlg_calculate_rest_biomass.exec_()

    def btn_forest_clicked(self):
        self.close()
        result = self.dlg_calculate_forest.exec_()
        # self.close()
        # result = self.dlg_calculate_forest.exec_()
        # QtWidgets.QMessageBox.information(None, self.tr("Coming soon!"),
        #                         self.tr("Forest coming soon!"))

    def btn_medalus_clicked(self):
        self.close()
        result = self.dlg_calculate_medalus.exec_()
        # QtWidgets.QMessageBox.information(None, self.tr("Coming soon!"),
        #                         self.tr("MEDALUS coming soon!"))

    def btn_timeseries_clicked(self):
        self.close()
        result = self.dlg_timeseries.exec_()
        # QtWidgets.QMessageBox.information(None, self.tr("Coming soon!"),
        #                         self.tr("MEDALUS coming soon!"))
    
    #This listens for when the coastal erosion button is clicked
    def btn_ce_clicked(self):
        self.close()
        result=self.dlg_calculate_ce.exec_()

    def btn_se_clicked(self):
        self.close()
        result = self.dlg_calculate_se.exec_() 
    # def btn_urban_clicked(self):
    #     self.close()
    #     result = self.dlg_calculate_urban.exec_()

"""
    Add the ui for computing soil erosion indicators
"""  
class DlgCalculateSE(QtWidgets.QDialog, Ui_DlgCalculateSE):
    def __init__(self,parent=None):
        super(DlgCalculateSE,self).__init__(parent)
        self.setupUi(self)
        self.folder=""
        self.countries=[]
        
        for item in get_Adm0():
            self.countries.append(item['name_0'])
        
        self.se_country_combobox.addItems(sorted(self.countries))
        self.populate_region()

        #handle events
        self.pushButton_se_wind.clicked.connect(self.btn_wind_clicked)
        self.pushButton_se_water.clicked.connect(self.btn_water_clicked)
        self.se_country_combobox.currentTextChanged.connect(self.populate_region)
        self.pushButton_se_getFolder.clicked.connect(self.getFolder)

    def btn_wind_clicked(self):
        if self.folder:
            self.close()
            result = self.run_download_ILSWE()
        else:
            self.show_error_message("Please select an output folder.")

    def btn_water_clicked(self):
        if self.folder:
            self.close()
            result = self.run_download_RUSLE()
        else:
            self.show_error_message("Please select an output folder.")
        


    def getFolder(self):
        output_dir=QtWidgets.QFileDialog.getExistingDirectory(None,"Select folder")
        self.folder = output_dir
        self.se_lineEdit.setText(output_dir)

    def show_error_message(self,message):
        msg_box = QtWidgets.QMessageBox()
        msg_box.setIcon(QtWidgets.QMessageBox.Critical)
        msg_box.setWindowTitle("Error")
        msg_box.setText(message)
        msg_box.exec_()

    def run_download_ILSWE(self):
        progress_dialog = QtWidgets.QProgressDialog("Processing...", "Cancel", 0, 100, self)
        progress_dialog.setWindowTitle("Downloading ILSWE")
        progress_dialog.setWindowModality(Qt.WindowModal)

        progress_bar = QtWidgets.QProgressBar(progress_dialog)
        progress_dialog.setBar(progress_bar)

        progress_dialog.show()
        QtWidgets.QApplication.processEvents()  # Ensure the dialog is shown before starting
        selected_region=self.se_region_combobox.currentText()
        selected_country = self.se_country_combobox.currentText()
        year = self.se_dateEdit.date().year()
        getILSWE(self.folder, progress_dialog, progress_bar,selected_country,selected_region,year)

    def run_download_RUSLE(self):
        progress_dialog = QtWidgets.QProgressDialog("Processing...", "Cancel", 0, 100, self)
        progress_dialog.setWindowTitle("Downloading RUSLE")
        progress_dialog.setWindowModality(Qt.WindowModal)

        progress_bar = QtWidgets.QProgressBar(progress_dialog)
        progress_dialog.setBar(progress_bar)

        progress_dialog.show()
        QtWidgets.QApplication.processEvents()  # Ensure the dialog is shown before starting
        selected_region=self.se_region_combobox.currentText()
        selected_country = self.se_country_combobox.currentText()
        year = self.se_dateEdit.date().year()
        getRUSLE(self.folder, progress_dialog, progress_bar,selected_country,selected_region,year)
    def populate_region(self):
        self.se_region_combobox.clear()
        selected_country = self.se_country_combobox.currentText()
        
        regions = []
        for item in get_Adm0():
            if selected_country == item['name_0']:
                id = item['id']
                self.regions = get_Adm1(id)
                for region in get_Adm1(id):
                    regions.append(region['name_1'])

        self.se_region_combobox.addItems(sorted(regions))
"""
    Add the ui interactivity for computing the cvi sub-indicators
"""
class DlgCalculateCE(QtWidgets.QDialog, Ui_DlgCalculateCE):
    def __init__(self,parent=None):
        super(DlgCalculateCE,self).__init__(parent)

        self.setupUi(self)
        # self.downloadGeomorphology = downloadGeomorphology()
        self.get_Adm0 = get_Adm0()
        self.folder = ""
        self.countries = []
        self.regions=[]
        self.sub_regions=[]
        
        for item in self.get_Adm0:
            self.countries.append(item['name_0'])
        
        self.cvi_country_combobox.addItems(sorted(self.countries))
        self.populate_region()
        self.cvi_country_combobox.currentTextChanged.connect(self.populate_region)
        # self.cvi_region_combobox.currentTextChanged.connect(self.populate_sub_region)
        self.pushButton_cvi_sub.clicked.connect(self.btn_cvi_sub_clicked)
        self.pushButton_cvi_getFolder.clicked.connect(self.getFolder)
        self.pushButton_cvi.clicked.connect(self.btn_cvi_clicked)
    def show_error_message(self,message):
        msg_box = QtWidgets.QMessageBox()
        msg_box.setIcon(QtWidgets.QMessageBox.Critical)
        msg_box.setWindowTitle("Error")
        msg_box.setText(message)
        msg_box.exec_()
    
    def getFolder(self):
        output_dir=QtWidgets.QFileDialog.getExistingDirectory(None,"Select folder")
        self.folder = output_dir
        self.cvi_lineEdit.setText(output_dir)
    
    def btn_cvi_sub_clicked(self):
        if self.folder:
            self.close()
            result = self.run_download_sub()
        else:
            self.show_error_message("Please select an output folder.")

    def btn_cvi_clicked(self):
        if self.folder:
            self.close()
            result = self.run_download_cvi()
        else:
            self.show_error_message("Please select an output folder")

    def run_download_sub(self):
        progress_dialog = QtWidgets.QProgressDialog("Processing...", "Cancel", 0, 100, self)
        progress_dialog.setWindowTitle("Downloading Sub-indicators")
        progress_dialog.setWindowModality(Qt.WindowModal)

        progress_bar = QtWidgets.QProgressBar(progress_dialog)
        progress_dialog.setBar(progress_bar)

        progress_dialog.show()
        QtWidgets.QApplication.processEvents()  # Ensure the dialog is shown before starting
        selected_region=self.cvi_region_combobox.currentText()
        downloadGeomorphology(self.folder, progress_dialog, progress_bar,selected_region)
    
    def run_download_cvi(self):
        progress_dialog = QtWidgets.QProgressDialog("Processing...", "Cancel", 0, 100, self)
        progress_dialog.setWindowTitle("Downloading CVI")
        progress_dialog.setWindowModality(Qt.WindowModal)

        progress_bar = QtWidgets.QProgressBar(progress_dialog)
        progress_dialog.setBar(progress_bar)

        progress_dialog.show()
        QtWidgets.QApplication.processEvents()  # Ensure the dialog is shown before starting
        selected_region=self.cvi_region_combobox.currentText()
        downloadCVI(self.folder, progress_dialog, progress_bar,selected_region)
    def populate_region(self):
        self.cvi_region_combobox.clear()
        selected_country = self.cvi_country_combobox.currentText()
        
        regions = []
        for item in self.get_Adm0:
            if selected_country == item['name_0']:
                id = item['id']
                self.regions = get_Adm1(id)
                for region in get_Adm1(id):
                    regions.append(region['name_1'])

        self.cvi_region_combobox.addItems(sorted(regions))



class DlgCalculateLD(QtWidgets.QDialog, Ui_DlgCalculateLD):
    def __init__(self, parent=None):
        super(DlgCalculateLD, self).__init__(parent)

        self.setupUi(self)

        # TODO: Bad style - fix when refactoring
        from MISLAND.calculate_prod import DlgCalculateProd
        from MISLAND.calculate_lc import DlgCalculateLC
        from MISLAND.calculate_soc import DlgCalculateSOC
        from MISLAND.calculate_ldn import DlgCalculateOneStep, DlgCalculateLDNSummaryTableAdmin
        self.dlg_calculate_prod = DlgCalculateProd()
        self.dlg_calculate_lc = DlgCalculateLC()
        self.dlg_calculate_soc = DlgCalculateSOC()
        self.dlg_calculate_ldn_onestep = DlgCalculateOneStep()
        self.dlg_calculate_ldn_advanced = DlgCalculateLDNSummaryTableAdmin()

        self.btn_prod.clicked.connect(self.btn_prod_clicked)
        self.btn_lc.clicked.connect(self.btn_lc_clicked)
        self.btn_soc.clicked.connect(self.btn_soc_clicked)
        self.btn_sdg_onestep.clicked.connect(self.btn_sdg_onestep_clicked)
        self.btn_summary_single_polygon.clicked.connect(self.btn_summary_single_polygon_clicked)
        # self.btn_summary_multi_polygons.clicked.connect(self.btn_summary_multi_polygons_clicked)

    def btn_prod_clicked(self):
        self.close()
        result = self.dlg_calculate_prod.exec_()

    def btn_lc_clicked(self):
        self.close()
        result = self.dlg_calculate_lc.exec_()

    def btn_soc_clicked(self):
        self.close()
        result = self.dlg_calculate_soc.exec_()

    def btn_sdg_onestep_clicked(self):
        self.close()
        result = self.dlg_calculate_ldn_onestep.exec_()

    def btn_summary_single_polygon_clicked(self):
        self.close()
        result = self.dlg_calculate_ldn_advanced.exec_()

    # def btn_summary_multi_polygons_clicked(self):
    #     QtWidgets.QMessageBox.information(None, self.tr("Coming soon!"),
    #                                   self.tr("Multiple polygon summary table calculation coming soon!"))


class DlgCalculateTC(QtWidgets.QDialog, Ui_DlgCalculateTC):
    def __init__(self, parent=None):
        super(DlgCalculateTC, self).__init__(parent)

        self.setupUi(self)

        # TODO: Bad style - fix when refactoring
        from MISLAND.calculate_tc import DlgCalculateTCData
        from MISLAND.calculate_tc import DlgCalculateTCSummaryTable
        self.dlg_calculate_tc_data = DlgCalculateTCData()
        self.dlg_calculate_tc_summary = DlgCalculateTCSummaryTable()

        self.btn_calculate_carbon_change.clicked.connect(self.btn_calculate_carbon_change_clicked)
        self.btn_summary_single_polygon.clicked.connect(self.btn_summary_single_polygon_clicked)

    def btn_calculate_carbon_change_clicked(self):
        self.close()
        result = self.dlg_calculate_tc_data.exec_()

    def btn_summary_single_polygon_clicked(self):
        self.close()
        result = self.dlg_calculate_tc_summary.exec_()


class DlgCalculateRestBiomass(QtWidgets.QDialog, Ui_DlgCalculateRestBiomass):
    def __init__(self, parent=None):
        super(DlgCalculateRestBiomass, self).__init__(parent)

        self.setupUi(self)

        # TODO: Bad style - fix when refactoring
        from MISLAND.calculate_rest_biomass import DlgCalculateRestBiomassData
        from MISLAND.calculate_rest_biomass import DlgCalculateRestBiomassSummaryTable
        self.dlg_calculate_rest_biomass_data = DlgCalculateRestBiomassData()
        self.dlg_calculate_rest_biomass_summary = DlgCalculateRestBiomassSummaryTable()

        self.btn_calculate_rest_biomass_change.clicked.connect(self.btn_calculate_rest_biomass_change_clicked)
        self.btn_summary_single_polygon.clicked.connect(self.btn_summary_single_polygon_clicked)

    def btn_calculate_rest_biomass_change_clicked(self):
        self.close()
        result = self.dlg_calculate_rest_biomass_data.exec_()

    def btn_summary_single_polygon_clicked(self):
        self.close()
        result = self.dlg_calculate_rest_biomass_summary.exec_()

class DlgCalculateForest(QtWidgets.QDialog, Ui_DlgCalculateForest):
    def __init__(self, parent=None):
        super(DlgCalculateForest, self).__init__(parent)

        self.setupUi(self)
        from MISLAND.calculate_ff import DlgCalculateForestFire
                                  
        self.dlg_calculate_tc = DlgCalculateTC()
        self.dlg_calculate_ff = DlgCalculateForestFire()

        self.pushButton_tc.clicked.connect(self.btn_tc_clicked)
        self.pushButton_ff.clicked.connect(self.btn_ff_clicked)

    def btn_tc_clicked(self):
        self.close()
        result = self.dlg_calculate_tc.exec_()

    def btn_ff_clicked(self):
        self.close()
        result = self.dlg_calculate_ff.exec_()

class DlgCalculateMedalus(QtWidgets.QDialog, Ui_DlgCalculateMedalus):
    def __init__(self, parent=None):
        super(DlgCalculateMedalus, self).__init__(parent)

        self.setupUi(self)

        # TODO: Bad style - fix when refactoring
        from MISLAND.calculate_cqi import DlgCalculateCQI
        from MISLAND.calculate_mqi import DlgCalculateMQI
        from MISLAND.calculate_sqi import DlgCalculateSQI
        from MISLAND.calculate_vqi import DlgCalculateVQI
        from MISLAND.calculate_sdi import DlgCalculateSDISummaryTableAdmin

        self.dlg_calculate_cqi = DlgCalculateCQI() 
        self.dlg_calculate_mqi = DlgCalculateMQI()
        self.dlg_calculate_sqi = DlgCalculateSQI()
        self.dlg_calculate_vqi = DlgCalculateVQI()
        self.dlg_calculate_sdi = DlgCalculateSDISummaryTableAdmin()

        self.btn_calculate_cqi.clicked.connect(self.btn_calculate_cqi_clicked)
        self.btn_calculate_mqi.clicked.connect(self.btn_calculate_mqi_clicked)
        self.btn_calculate_sqi.clicked.connect(self.btn_calculate_sqi_clicked)
        self.btn_calculate_vqi.clicked.connect(self.btn_calculate_vqi_clicked)
        self.btn_calculate_sdi.clicked.connect(self.btn_calculate_sdi_clicked)

    def btn_calculate_cqi_clicked(self):
        self.close()
        result = self.dlg_calculate_cqi.exec_()
    def btn_calculate_mqi_clicked(self):
        self.close()
        result = self.dlg_calculate_mqi.exec_()

    def btn_calculate_sqi_clicked(self):
        self.close()
        result = self.dlg_calculate_sqi.exec_()

    def btn_calculate_vqi_clicked(self):
        self.close()
        result = self.dlg_calculate_vqi.exec_()
        
    def btn_calculate_sdi_clicked(self):
        self.close()
        result = self.dlg_calculate_sdi.exec_()
        # QtWidgets.QMessageBox.information(None, self.tr("Coming soon!"),
        #                         self.tr("Final Medalus Calculation coming soon!"))


class CalculationOptionsWidget(QtWidgets.QWidget, Ui_WidgetCalculationOptions):
    def __init__(self, parent=None):
        super(CalculationOptionsWidget, self).__init__(parent)

        self.setupUi(self)

        self.radioButton_run_in_cloud.toggled.connect(self.radioButton_run_in_cloud_changed)
        self.btn_local_data_folder_browse.clicked.connect(self.open_folder_browse)

        self.task_name.textChanged.connect(self.task_name_changed)
        self.task_notes.textChanged.connect(self.task_notes_changed)

    def showEvent(self, event):
        super(CalculationOptionsWidget, self).showEvent(event)

        local_data_folder = QSettings().value("MISLAND/localdata_dir", None)
        if local_data_folder and os.access(local_data_folder, os.R_OK):
            self.lineEdit_local_data_folder.setText(local_data_folder)
        else:
            self.lineEdit_local_data_folder.setText(None)
        self.task_name.setText(QSettings().value("MISLAND/CalculationOptionsWidget/task_name", None))
        self.task_notes.setHtml(QSettings().value("MISLAND/CalculationOptionsWidget/task_notes", None))

    def task_name_changed(self, value=None):
        QSettings().setValue("MISLAND/CalculationOptionsWidget/task_name", value)

    def task_notes_changed(self):
        QSettings().setValue("MISLAND/CalculationOptionsWidget/task_notes", self.task_notes.toHtml())

    def radioButton_run_in_cloud_changed(self):
        if self.radioButton_run_in_cloud.isChecked():
            self.lineEdit_local_data_folder.setEnabled(False)
            self.btn_local_data_folder_browse.setEnabled(False)
        else:
            self.lineEdit_local_data_folder.setEnabled(True)
            self.btn_local_data_folder_browse.setEnabled(True)

    def open_folder_browse(self):
        self.lineEdit_local_data_folder.clear()

        folder = QtWidgets.QFileDialog.getExistingDirectory(self,
                                                        self.tr('Select folder containing data'),
                                                        QSettings().value("MISLAND/localdata_dir", None))
        if folder:
            if os.access(folder, os.R_OK):
                QSettings().setValue("MISLAND/localdata_dir", os.path.dirname(folder))
                self.lineEdit_local_data_folder.setText(folder)
                return True
            else:
                QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                           self.tr(u"Cannot read {}. Choose a different folder.".format(folder)))
                return False
        else:
            return False

    def toggle_show_where_to_run(self, enable):
        if enable:
            self.where_to_run_enabled = True
            self.groupBox_where_to_run.show()
        else:
            self.where_to_run_enabled = False
            self.groupBox_where_to_run.hide()

class CalculationOutputWidget(QtWidgets.QWidget, Ui_WidgetCalculationOutput):
    def __init__(self, suffixes, subclass_name, parent=None):
        super(CalculationOutputWidget, self).__init__(parent)

        self.output_suffixes = suffixes
        self.subclass_name = subclass_name

        self.setupUi(self)

        self.browse_output_basename.clicked.connect(self.select_output_basename)

    def select_output_basename(self):
        local_name = QSettings().value("MISLAND/output_basename_{}".format(self.subclass_name), None)
        if local_name:
            initial_path = local_name
        else:
            initial_path = QSettings().value("MISLAND/output_dir", None)


        f, _ = QtWidgets.QFileDialog.getSaveFileName(self,
                self.tr('Choose a prefix to be used when naming output files'),
                initial_path,
                self.tr('Base name (*)'))

        if f:
            if os.access(os.path.dirname(f), os.W_OK):
                QSettings().setValue("MISLAND/output_dir", os.path.dirname(f))
                QSettings().setValue("MISLAND/output_basename_{}".format(self.subclass_name), f)
                self.output_basename.setText(f)
                self.set_output_summary(f)
            else:
                QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                           self.tr(u"Cannot write to {}. Choose a different file.".format(f)))

    def set_output_summary(self, f):
        out_files = [f + suffix for suffix in self.output_suffixes]
        self.output_summary.setText("\n".join(["{}"]*len(out_files)).format(*out_files))

    def check_overwrites(self):
        overwrites = []
        for suffix in self.output_suffixes: 
            if os.path.exists(self.output_basename.text() + suffix):
                overwrites.append(os.path.basename(self.output_basename.text() + suffix))

        if len(overwrites) > 0:
            resp = QtWidgets.QMessageBox.question(self,
                    self.tr('Overwrite file?'),
                    self.tr('Using the prefix "{}" would lead to overwriting existing file(s) {}. Do you want to overwrite these file(s)?'.format(
                        self.output_basename.text(),
                        ", ".join(["{}"]*len(overwrites)).format(*overwrites))),
                    QtWidgets.QMessageBox.Yes, QtWidgets.QMessageBox.No)
            if resp == QtWidgets.QMessageBox.No:
                QtWidgets.QMessageBox.information(None, self.tr("Information"),
                                           self.tr(u"Choose a different output prefix and try again."))
                return False

        return True


class AreaWidget(QtWidgets.QWidget, Ui_WidgetSelectArea):
    def __init__(self, parent=None):
        super(AreaWidget, self).__init__(parent)

        self.setupUi(self)

        self.canvas = iface.mapCanvas()

        self.admin_bounds_key = get_admin_bounds()
        if not self.admin_bounds_key:
            raise ValueError('Admin boundaries not available')

        self.cities = get_cities()
        if not self.cities:
            raise ValueError('Cities list not available')

        # Populate
        self.area_admin_0.addItems(sorted(self.admin_bounds_key.keys()))
        self.populate_admin_1()
        self.populate_cities()
        self.area_admin_0.currentIndexChanged.connect(self.populate_admin_1)
        self.area_admin_0.currentIndexChanged.connect(self.populate_cities)

        # Handle saving to qsettings
        self.area_admin_0.activated.connect(self.admin_0_changed)
        self.secondLevel_area_admin_1.activated.connect(self.admin_1_changed)
        self.secondLevel_city.activated.connect(self.admin_city_changed)
        self.area_frompoint_point_x.textChanged.connect(self.point_x_changed)
        self.area_frompoint_point_y.textChanged.connect(self.point_y_changed)
        self.area_fromfile_file.textChanged.connect(self.file_changed)
        self.groupBox_buffer.clicked.connect(self.buffer_changed)
        self.buffer_size_km.valueChanged.connect(self.buffer_changed)

        self.area_fromfile_browse.clicked.connect(self.open_vector_browse)
        self.area_fromadmin.clicked.connect(self.area_type_toggle)
        self.area_fromfile.clicked.connect(self.area_type_toggle)

        self.radioButton_secondLevel_region.clicked.connect(self.radioButton_secondLevel_toggle)
        self.radioButton_secondLevel_city.clicked.connect(self.radioButton_secondLevel_toggle)

        icon = QIcon(QPixmap(':/plugins/MISLAND/icons/map-marker.svg'))
        self.area_frompoint_choose_point.setIcon(icon)
        self.area_frompoint_choose_point.clicked.connect(self.point_chooser)
        #TODO: Set range to only accept valid coordinates for current map coordinate system
        self.area_frompoint_point_x.setValidator(QDoubleValidator())
        #TODO: Set range to only accept valid coordinates for current map coordinate system
        self.area_frompoint_point_y.setValidator(QDoubleValidator())
        self.area_frompoint.clicked.connect(self.area_type_toggle)

        # Setup point chooser
        self.choose_point_tool = QgsMapToolEmitPoint(self.canvas)
        self.choose_point_tool.canvasClicked.connect(self.set_point_coords)

        proj_crs = QgsCoordinateReferenceSystem(self.canvas.mapSettings().destinationCrs().authid())
        self.mQgsProjectionSelectionWidget.setCrs(QgsCoordinateReferenceSystem('epsg:4326'))

    def showEvent(self, event):
        super(AreaWidget, self).showEvent(event)

        buffer_checked = QSettings().value("MISLAND/AreaWidget/buffer_checked", False) == 'True'
        area_from_option = QSettings().value("MISLAND/AreaWidget/area_from_option", None)

        if area_from_option == 'admin':
            self.area_fromadmin.setChecked(True)
        elif area_from_option == 'point':
            self.area_frompoint.setChecked(True)
        elif area_from_option == 'file':
            self.area_fromfile.setChecked(True)
        self.area_frompoint_point_x.setText(QSettings().value("MISLAND/AreaWidget/area_frompoint_point_x", None))
        self.area_frompoint_point_y.setText(QSettings().value("MISLAND/AreaWidget/area_frompoint_point_y", None))
        self.area_fromfile_file.setText(QSettings().value("MISLAND/AreaWidget/area_fromfile_file", None))
        self.area_type_toggle(False)

        admin_0 = QSettings().value("MISLAND/AreaWidget/area_admin_0", None)
        if admin_0:
            self.area_admin_0.setCurrentIndex(self.area_admin_0.findText(admin_0))
            self.populate_admin_1()

        area_from_option_secondLevel = QSettings().value("MISLAND/AreaWidget/area_from_option_secondLevel", None)
        if area_from_option_secondLevel == 'admin':
            self.radioButton_secondLevel_region.setChecked(True)
        elif area_from_option_secondLevel == 'city':
            self.radioButton_secondLevel_city.setChecked(True)
        self.radioButton_secondLevel_toggle(False)

        secondLevel_area_admin_1 = QSettings().value("MISLAND/AreaWidget/secondLevel_area_admin_1", None)
        if secondLevel_area_admin_1:
            self.secondLevel_area_admin_1.setCurrentIndex(self.secondLevel_area_admin_1.findText(secondLevel_area_admin_1))
        secondLevel_city = QSettings().value("MISLAND/AreaWidget/secondLevel_city", None)
        if secondLevel_city:
            self.populate_cities()
            self.secondLevel_city.setCurrentIndex(self.secondLevel_city.findText(secondLevel_city))

        buffer_size = QSettings().value("MISLAND/AreaWidget/buffer_size", None)
        if buffer_size:
            self.buffer_size_km.setValue(int(buffer_size))
        self.groupBox_buffer.setChecked(buffer_checked)

    def admin_0_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/area_admin_0", self.area_admin_0.currentText())

    def admin_1_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/secondLevel_area_admin_1", self.secondLevel_area_admin_1.currentText())
        
    def admin_city_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/secondLevel_city", self.secondLevel_city.currentText())

    def file_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/area_fromfile_file", self.area_fromfile_file.text())

    def point_x_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/area_frompoint_point_x", self.area_frompoint_point_x.text())

    def point_y_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/area_frompoint_point_y", self.area_frompoint_point_y.text())

    def buffer_changed(self):
        QSettings().setValue("MISLAND/AreaWidget/buffer_checked", str(self.groupBox_buffer.isChecked()))
        QSettings().setValue("MISLAND/AreaWidget/buffer_size", self.buffer_size_km.value())

    def populate_cities(self):
        self.secondLevel_city.clear()
        adm0_a3 = self.admin_bounds_key[self.area_admin_0.currentText()]['code']
        self.current_cities_key = {value['name_en']: key for key, value in self.cities[adm0_a3].items()}
        self.secondLevel_city.addItems(sorted(self.current_cities_key.keys()))

    def populate_admin_1(self):
        self.secondLevel_area_admin_1.clear()
        # self.secondLevel_area_admin_1.addItems(['All regions'])
        self.secondLevel_area_admin_1.addItems(sorted(self.admin_bounds_key[self.area_admin_0.currentText()]['admin1'].keys()))

    def area_type_toggle(self, save=True):
        if self.area_frompoint.isChecked():
            if save: QSettings().setValue("MISLAND/AreaWidget/area_from_option", 'point')
            self.area_frompoint_point_x.setEnabled(True)
            self.area_frompoint_point_y.setEnabled(True)
            self.area_frompoint_choose_point.setEnabled(True)
        else:
            self.area_frompoint_point_x.setEnabled(False)
            self.area_frompoint_point_y.setEnabled(False)
            self.area_frompoint_choose_point.setEnabled(False)

        if self.area_fromadmin.isChecked():
            if save: QSettings().setValue("MISLAND/AreaWidget/area_from_option", 'admin')
            self.groupBox_first_level.setEnabled(True)
            self.groupBox_second_level.setEnabled(True)
        else:
            self.groupBox_first_level.setEnabled(False)
            self.groupBox_second_level.setEnabled(False)

        if self.area_fromfile.isChecked():
            if save: QSettings().setValue("MISLAND/AreaWidget/area_from_option", 'file')
            self.area_fromfile_file.setEnabled(True)
            self.area_fromfile_browse.setEnabled(True)
        else:
            self.area_fromfile_file.setEnabled(False)
            self.area_fromfile_browse.setEnabled(False)

    def radioButton_secondLevel_toggle(self, save=True):
        if self.radioButton_secondLevel_region.isChecked():
            if save: QSettings().setValue("MISLAND/AreaWidget/area_from_option_secondLevel", 'admin')
            self.secondLevel_area_admin_1.setEnabled(True)
            self.secondLevel_city.setEnabled(False)
        else:
            if save: QSettings().setValue("MISLAND/AreaWidget/area_from_option_secondLevel", 'city')
            self.secondLevel_area_admin_1.setEnabled(False)
            self.secondLevel_city.setEnabled(True)

    def point_chooser(self):
        log("Choosing point from canvas...")
        self.canvas.setMapTool(self.choose_point_tool)
        self.window().hide()
        QtWidgets.QMessageBox.critical(None, self.tr("Point chooser"), self.tr("Click the map to choose a point."))

    def set_point_coords(self, point, button):
        log("Set point coords")
        #TODO: Show a messagebar while tool is active, and then remove the bar when a point is chosen.
        self.point = point
        # Disable the choose point tool
        self.canvas.setMapTool(QgsMapToolPan(self.canvas))
        # Don't reset_tab_on_show as it would lead to return to first tab after
        # using the point chooser
        self.window().reset_tab_on_showEvent = False
        self.window().show()
        self.window().reset_tab_on_showEvent = True
        self.point = self.canvas.getCoordinateTransform().toMapCoordinates(self.canvas.mouseLastXY())
        log("Chose point: {}, {}.".format(self.point.x(), self.point.y()))
        self.area_frompoint_point_x.setText("{:.8f}".format(self.point.x()))
        self.area_frompoint_point_y.setText("{:.8f}".format(self.point.y()))

    def open_vector_browse(self):
        initial_path = QSettings().value("MISLAND/input_shapefile", None)
        if not initial_path:
            initial_path = QSettings().value("MISLAND/input_shapefile_dir", None)
        if not initial_path:
            initial_path = str(Path.home())

        vector_file, _ = QtWidgets.QFileDialog.getOpenFileName(self,
                                                        self.tr('Select a file defining the area of interest'),
                                                        initial_path,
                                                        self.tr('Vector file (*.shp *.kml *.kmz *.geojson)'))
        if vector_file:
            if os.access(vector_file, os.R_OK):
                QSettings().setValue("MISLAND/input_shapefile", vector_file)
                QSettings().setValue("MISLAND/input_shapefile_dir", os.path.dirname(vector_file))
                self.area_fromfile_file.setText(vector_file)
                return True
            else:
                QtWidgets.QMessageBox.critical(None,
                                               self.tr("Error"),
                                               self.tr(u"Cannot read {}. Choose a different file.".format(vector_file)))
                return False
        else:
            return False
                

class DlgCalculateBase(QtWidgets.QDialog):
    """Base class for individual indicator calculate dialogs"""
    firstShowEvent = pyqtSignal()

    @classmethod
    def get_subclass_name(cls):
        return cls.__name__

    def __init__(self, parent=None):
        super(DlgCalculateBase, self).__init__(parent)

        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'data', 'scripts.json')) as script_file:
            self.scripts = json.load(script_file)

        with open(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                               'data', 'gee_datasets.json')) as datasets_file:
            self.datasets = json.load(datasets_file)

        self._has_output = False
        self._firstShowEvent = True
        self.reset_tab_on_showEvent = True
        self._max_area = 0.63e6 # maximum size task the tool supports

        self.firstShowEvent.connect(self.firstShow)

    def add_output_tab(self, suffixes=['.json', '.tif']):
        self._has_output = True
        self.output_suffixes = suffixes

    def showEvent(self, event):
        super(DlgCalculateBase, self).showEvent(event)

        if self._firstShowEvent:
            self._firstShowEvent = False
            self.firstShowEvent.emit()

        if self.reset_tab_on_showEvent:
            self.TabBox.setCurrentIndex(0)
        
        # If this dialog has an output_basename widget then set it up with any 
        # saved values in QSettings
        if self._has_output:
            f = QSettings().value("MISLAND/output_basename_{}".format(self.get_subclass_name()), None)
            if f:
                self.output_tab.output_basename.setText(f)
                self.output_tab.set_output_summary(f)

    def firstShow(self):
        self.area_tab = AreaWidget()
        self.area_tab.setParent(self)
        self.TabBox.addTab(self.area_tab, self.tr('Area'))

        if self._has_output:
            self.output_tab = CalculationOutputWidget(self.output_suffixes, self.get_subclass_name())
            self.output_tab.setParent(self)
            self.TabBox.addTab(self.output_tab, self.tr('Output'))

        self.options_tab = CalculationOptionsWidget()
        self.options_tab.setParent(self)
        self.TabBox.addTab(self.options_tab, self.tr('Options'))

        # By default show the local or cloud option
        #self.options_tab.toggle_show_where_to_run(True)
        self.options_tab.toggle_show_where_to_run(False)

        # By default hide the custom crs box
        self.area_tab.groupBox_custom_crs.hide()
        
        self.button_calculate.clicked.connect(self.btn_calculate)
        self.button_prev.clicked.connect(self.tab_back)
        self.button_next.clicked.connect(self.tab_forward)

        # Start on first tab so button_prev and calculate should be disabled
        self.button_prev.setEnabled(False)
        self.button_calculate.setEnabled(False)
        self.TabBox.currentChanged.connect(self.tab_changed)

    def tab_back(self):
        if self.TabBox.currentIndex() - 1 >= 0:
            self.TabBox.setCurrentIndex(self.TabBox.currentIndex() - 1)

    def tab_forward(self):
        if self.TabBox.currentIndex() + 1 < self.TabBox.count():
            self.TabBox.setCurrentIndex(self.TabBox.currentIndex() + 1)

    def tab_changed(self):
        if self.TabBox.currentIndex() > 0:
            self.button_prev.setEnabled(True)
        else:
            self.button_prev.setEnabled(False)

        if self.TabBox.currentIndex() < (self.TabBox.count() - 1):
            self.button_next.setEnabled(True)
        else:
            self.button_next.setEnabled(False)

        if self.TabBox.currentIndex() == (self.TabBox.count() - 1):
            self.button_calculate.setEnabled(True)
        else:
            self.button_calculate.setEnabled(False)

    def btn_cancel(self):
        self.close()

    def get_city_geojson(self):
        adm0_a3 = self.area_tab.admin_bounds_key[self.area_tab.area_admin_0.currentText()]['code']
        wof_id = self.area_tab.current_cities_key[self.area_tab.secondLevel_city.currentText()]
        return (self.area_tab.cities[adm0_a3][wof_id]['geojson'])

    def get_admin_poly_geojson(self):
        adm0_a3 = self.area_tab.admin_bounds_key[self.area_tab.area_admin_0.currentText()]['code']
        admin_polys = read_json(u'admin_bounds_polys_{}.json.gz'.format(adm0_a3), verify=False)
        if not admin_polys:
            return None

        if not self.area_tab.secondLevel_area_admin_1.currentText() or self.area_tab.secondLevel_area_admin_1.currentText() == 'All regions':
            return (admin_polys['geojson'])
        else:
            admin_1_code = self.area_tab.admin_bounds_key[self.area_tab.area_admin_0.currentText()]['admin1'][self.area_tab.secondLevel_area_admin_1.currentText()]['code']
            return (admin_polys['admin1'][admin_1_code]['geojson'])

    def btn_calculate(self):
        if self.area_tab.groupBox_custom_crs.isChecked():
            crs_dst = self.area_tab.mQgsProjectionSelectionWidget.crs()
        else:
            crs_dst = QgsCoordinateReferenceSystem('epsg:4326')

        self.aoi = AOI(crs_dst)

        if self.area_tab.area_fromadmin.isChecked():
            if self.area_tab.radioButton_secondLevel_city.isChecked():
                if not self.area_tab.groupBox_buffer.isChecked():
                    QtWidgets.QMessageBox.critical(None,tr_calculate.tr("Error"),
                           tr_calculate.tr("You have chosen to run calculations for a city. You must select a buffer distance to define the calculation area when you are processing a city."))
                    return False
                geojson = self.get_city_geojson()
                self.aoi.update_from_geojson(geojson=geojson, 
                                             wrap=self.area_tab.checkBox_custom_crs_wrap.isChecked(),
                                             datatype='point')
            else:
                if not self.area_tab.area_admin_0.currentText():
                    QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                               self.tr("Choose a first level administrative boundary."))
                    return False
                self.button_calculate.setEnabled(False)
                geojson = self.get_admin_poly_geojson()
                self.button_calculate.setEnabled(True)
                if not geojson:
                    QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                               self.tr("Unable to load administrative boundaries."))
                    return False
                self.aoi.update_from_geojson(geojson=geojson, 
                                             wrap=self.area_tab.checkBox_custom_crs_wrap.isChecked())
        elif self.area_tab.area_fromfile.isChecked():
            if not self.area_tab.area_fromfile_file.text():
                QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                           self.tr("Choose a file to define the area of interest."))
                return False
            if not os.access(self.area_tab.area_fromfile_file.text(), os.R_OK):
                QtWidgets.QMessageBox.critical(None,
                                               self.tr("Error"),
                                               self.tr("Unable to read {}.".format(self.area_tab.area_fromfile_file.text())))
                return False
            self.aoi.update_from_file(f=self.area_tab.area_fromfile_file.text(),
                                      wrap=self.area_tab.checkBox_custom_crs_wrap.isChecked())
        elif self.area_tab.area_frompoint.isChecked():
            # Area from point
            if not self.area_tab.area_frompoint_point_x.text() or not self.area_tab.area_frompoint_point_y.text():
                QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                           self.tr("Choose a point to define the area of interest."))
                return False
            point = QgsPointXY(float(self.area_tab.area_frompoint_point_x.text()), float(self.area_tab.area_frompoint_point_y.text()))
            crs_src = QgsCoordinateReferenceSystem(self.area_tab.canvas.mapSettings().destinationCrs().authid())
            point = QgsCoordinateTransform(crs_src, crs_dst, QgsProject.instance()).transform(point)
            geojson = json.loads(QgsGeometry.fromPointXY(point).asJson())
            self.aoi.update_from_geojson(geojson=geojson, 
                                         wrap=self.area_tab.checkBox_custom_crs_wrap.isChecked(),
                                         datatype='point')
        else:
            QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                       self.tr("Choose an area of interest."))
            return False

        if self.aoi and (self.aoi.datatype == "unknown" or not self.aoi.isValid()):
            QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                                       self.tr("Unable to read area of interest."))
            return False

        if self.area_tab.groupBox_buffer.isChecked():
            ret = self.aoi.buffer(self.area_tab.buffer_size_km.value())
            if not ret:
                QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                        self.tr("Error buffering polygon"))
                return False

        # Check that a bounding box can be created successfully from this aoi
        ret = self.aoi.bounding_box_gee_geojson()
        if not ret:
            QtWidgets.QMessageBox.critical(None,
                                           self.tr("Error"),
                                           self.tr("Unable to calculate bounding box."))
            return False
        else:
            self.gee_bounding_box = ret

        # Limit processing area to be no greater than 10^7 sq km if using a 
        # custom shapefile
        # if not self.area_tab.area_fromadmin.isChecked():
        aoi_area = self.aoi.get_area() / (1000 * 1000)
        if aoi_area > self._max_area:
            QtWidgets.QMessageBox.critical(None, self.tr("Error"),
                    self.tr("The bounding box for the requested area (approximately {0}) sq km is too large. Choose a smaller area to process. Maximum area {1} sq km".format(aoi_area, self._max_area)))
            return False

        if self._has_output:
            if not self.output_tab.output_basename.text():
                QtWidgets.QMessageBox.information(None, self.tr("Error"),
                                              self.tr("Choose an output base name."))
                return False

            # Check if the chosen basename would lead to an  overwrite(s):
            ret = self.output_tab.check_overwrites()
            if not ret:
                return False

        return True


class ClipWorker(AbstractWorker):
    def __init__(self, in_file, out_file, geojson, output_bounds=None):
        AbstractWorker.__init__(self)

        self.in_file = in_file
        self.out_file = out_file
        self.output_bounds = output_bounds

        self.geojson = geojson

    def work(self):
        self.toggle_show_progress.emit(True)
        self.toggle_show_cancel.emit(True)

        json_file = GetTempFilename('.geojson')
        with open(json_file, 'w') as f:
            json.dump(self.geojson, f, separators=(',', ': '))

        gdal.UseExceptions()
        res = gdal.Warp(self.out_file, self.in_file, format='GTiff',
                        cutlineDSName=json_file, srcNodata=-32768, 
                        outputBounds=self.output_bounds,
                        dstNodata=-32767,
                        dstSRS="epsg:4326",
                        outputType=gdal.GDT_Int16,
                        resampleAlg=gdal.GRA_NearestNeighbour,
                        creationOptions=['COMPRESS=LZW'],
                        callback=self.progress_callback)
        os.remove(json_file)

        if res:
            return True
        else:
            return None

    def progress_callback(self, fraction, message, data):
        if self.killed:
            return False
        else:
            self.progress.emit(100 * fraction)
            return True


class MaskWorker(AbstractWorker):
    def __init__(self, out_file, geojson, model_file=None):
        AbstractWorker.__init__(self)

        self.out_file = out_file
        self.geojson = geojson
        self.model_file = model_file

    def work(self):
        self.toggle_show_progress.emit(True)
        self.toggle_show_cancel.emit(True)


        json_file = GetTempFilename('.geojson')
        with open(json_file, 'w') as f:
            json.dump(self.geojson, f, separators=(',', ': '))

        gdal.UseExceptions()

        if self.model_file:
            # Assumes an image with no rotation
            gt = gdal.Info(self.model_file, format='json')['geoTransform']
            x_size, y_size= gdal.Info(self.model_file, format='json')['size']
            x_min = min(gt[0], gt[0] + x_size * gt[1])
            x_max = max(gt[0], gt[0] + x_size * gt[1])
            y_min = min(gt[3], gt[3] + y_size * gt[5])
            y_max = max(gt[3], gt[3] + y_size * gt[5])
            output_bounds = [x_min, y_min, x_max, y_max]
            x_res = gt[1]
            y_res = gt[5]
        else:
            output_bounds = None
            x_res = None
            y_res = None

        res = gdal.Rasterize(self.out_file, json_file, format='GTiff',
                             outputBounds=output_bounds,
                             initValues=-32767, # Areas that are masked out
                             burnValues=1, # Areas that are NOT masked out
                             xRes=x_res,
                             yRes=y_res,
                             outputSRS="epsg:4326",
                             outputType=gdal.GDT_Int16,
                             creationOptions=['COMPRESS=LZW'],
                             callback=self.progress_callback)
        os.remove(json_file)

        if res:
            return True
        else:
            return None

    def progress_callback(self, fraction, message, data):
        if self.killed:
            return False
        else:
            self.progress.emit(100 * fraction)
            return True