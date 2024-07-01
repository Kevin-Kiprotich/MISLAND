import requests
import os
import ee
import geemap
from natsort import natsorted
import rasterio
import glob
from rasterio.merge import merge
from pyproj import CRS
import numpy as np
import shutil
from qgis.PyQt.QtWidgets import QApplication,QMessageBox
from qgis.PyQt.QtCore import Qt
from qgis.core import QgsProject, QgsRasterLayer
from qgis.utils import iface
from MISLAND import log
URL_PATH = "http://41.227.30.136:1337"

"""
    This function retrieves level 0 admin boundaries to be used in the country comboboxes.
"""

def get_Adm0():
    path = "/api/vect0/?include=all"
    #create a response object
    try:
        response = requests.get(URL_PATH + path)
        return response.json()
    except requests.exceptions.ConnectionError:
        log("Unable to access the MISLAND server")
        return []
    except requests.exceptions.Timeout:
        log('API unable to login - general error')
        return []
    


"""
    This function retrieves level 1 admin boundaries to be added region comboboxes.
"""
def get_Adm1(adm0_id):
    path= "/api/vect1/?include=all"
    try:
        response = requests.get(URL_PATH +path)
        all_adm1=response.json()
        adm1=[]
        for item in all_adm1:
            if item['admin_zero_id']==adm0_id:
                adm1.append(item)
        return adm1
    except requests.exceptions.ConnectionError:
        log("Unable to access the MISLAND server")
        return []
    except requests.exceptions.Timeout:
        log('API unable to login - general error')
        return []


"""
    This function retrieves level 2 admin boundaries to be added in the sub-region combo boxes
"""
def get_Adm2(adm1_id):
    path="/api/vect2/"
    try:
        response = requests.get(URL_PATH+path)
        all_adm2=response.json()
        adm2=[]
        for item in all_adm2:
            if item["admin_one_id"]==adm1_id:
                adm2.append(item)
        
        return adm2

    except requests.exceptions.ConnectionError:
        log("Unable to access the MISLAND server")
        return []
    
    except requests.exceptions.Timeout:
        log('API unable to login - general error')
        return []

def show_error_message(message):
        msg_box = QMessageBox()
        msg_box.setIcon(QMessageBox.Critical)
        msg_box.setWindowTitle("Error")
        msg_box.setText(message)
        msg_box.exec_()

def initialize():
    json_file="ee-kevinkiprotich0089-449197507985.json"
    script_directory = os.path.dirname(os.path.abspath(__file__))
    json_file_path = os.path.join(script_directory, json_file)
    service_account = 'kevinkiprotich@ee-kevinkiprotich0089.iam.gserviceaccount.com'
    credentials = ee.ServiceAccountCredentials(service_account, json_file_path)
    print(credentials)
    ee.Initialize(credentials)


#This is a function to add the layer to the map
def addMapLayer(path):
    log("Adding layer to QGIS")
    project = QgsProject.instance()
    filename=path.split('/')[-1][0:-4]
    layer=QgsRasterLayer(path,filename,'gdal')
    if layer.isValid():
        project.instance().addMapLayer(layer)
        canvas = iface.mapCanvas()
        extent = layer.extent()
        canvas.setExtent(extent)
        canvas.refresh()
    else:
        print("Failed to load layer")

def MergeGrids(in_folder,out_folder, bands,data_type,no_data=0):   
    q = in_folder
    # sorted files by name
    fp = natsorted(glob.glob(q)) 
    
    # List for storing the raster image
    src_files = []
    
    # Open each raster files by iterating and then append to our list
    for raster in fp:
        # open raster file
        files = rasterio.open(raster)
    
        # add each file to our list
        src_files.append(files)
    
    # Merge function returns a single mosaic array and the transformation info
    mosaic, out_trans = merge(src_files, dtype=data_type)
    
    # Set metadata
    out_meta = src_files[0].meta.copy()
    out_dtype = src_files[0].dtypes[0]
    print(f"Data type:\t{out_dtype}")
    out_meta.update({"driver": "GTiff",
                    "dtype": data_type,
                    "nodata": no_data,
                    "height": mosaic.shape[1],
                    "width": mosaic.shape[2],
                    "transform": out_trans,
                    "crs": src_files[0].crs
                })
                    
    # Write the mosaic raster
    # output = os.path.join('./output/'+folder+'/WaterMask/pre-processed-images'+'/WaterMask_'+file_name+'.tif')
    output = out_folder
    with rasterio.open(output, "w", tiled=True, compress='lzw', **out_meta) as dest:
        dest.write(mosaic.astype(data_type))

    for src in src_files:
        src.close()

def rasterInRegion(raster,region):
    mask = raster.mask().reduce(ee.Reducer.anyNonZero())
    masked = mask.reduceRegion(
        reducer= ee.Reducer.sum(),
        geometry= region,
        scale=100,
        maxPixels= 1e13
    )

    # Evaluate the result
    result = masked.getInfo()

    # Check if there's any non-null value
    contains_data = result['any'] > 0
    return contains_data

def downloadGeomorphology(output_folder,progress_dialog,progress_bar,region):
    initialize()
    
    log("Fetching data. Please wait!!!...")
    progress_bar.setValue(0)
    progress_dialog.setLabelText("Initializing...")
    africa_regions = ee.FeatureCollection('projects/ee-kevinkiprotich0089/assets/MISLAND/Africa_Regions')
    shapefile = africa_regions.filter(ee.Filter.eq('NAME_1',region))
    #These are the datasets
    progress_bar.setValue(10)
    progress_dialog.setLabelText("Loading geomorphology data...")
    geomorphology_CA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_CA')
    geomorphology_EA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_EA')
    geomorphology_NA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_NA')
    geomorphology_SA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_SA')
    geomorphology_WA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_WA')

    geomorphology = ee.ImageCollection([geomorphology_CA,geomorphology_EA,geomorphology_NA,geomorphology_SA,geomorphology_WA]).mosaic()
    progress_bar.setValue(20)
    referenceProjection=geomorphology.projection()
    geomorphology = geomorphology.resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    progress_bar.setValue(20)
    coastal_slope = ee.Image('projects/ee-kevinkiprotich0089/assets/coastal_slope').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    sea_level_rise = ee.Image('projects/ee-kevinkiprotich0089/assets/Sea_Level_Rise').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    shoreline_change = ee.Image('projects/ee-kevinkiprotich0089/assets/shoreline_change').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    tidal_range = ee.Image('projects/ee-kevinkiprotich0089/assets/tidal_range').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    wave_height = ee.Image('projects/ee-kevinkiprotich0089/assets/wave_height').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    
    #clip to test the dataset rasterinpolygon
    geomorphology_clipped=geomorphology.clip(shapefile)
    coastal_slope_clipped=coastal_slope.clip(shapefile)
    sea_level_rise_clipped=sea_level_rise.clip(shapefile)
    shoreline_change_clipped=shoreline_change.clip(shapefile)
    tidal_range_clipped=tidal_range.clip(shapefile)
    wave_height_clipped=wave_height.clip(shapefile)
    progress_bar.setValue(30)

    if rasterInRegion(geomorphology_clipped,shapefile):
        grid = geemap.fishnet(shapefile, h_interval=0.5, v_interval=0.5, delta=1)
        gridList = grid.toList(grid.size())
        grid_num = grid.toList(grid.size()).length()
        progress_bar.setValue(40)
        progress_dialog.setLabelText("Generating grid features...")
        # Create list of each grid feature

        ls_feature = []
        for i in range(grid_num.getInfo()):
            feature = ee.Feature(gridList.get(i)).geometry()
            ls_feature.append(feature)
            progress_bar.setValue(40 + int(10*(i/grid_num.getInfo())))
            QApplication.processEvents()

        # create a temporary folder to hold the grid datasets
        if not os.path.exists(output_folder+'/temp/'):
            os.makedirs(output_folder+'/temp/')
        
        progress_dialog.setLabelText("Downloading image...")
        for i in range(grid_num.getInfo()):
            log(f"Downloading image_{i}")
            download = ee.data.getDownloadId({
                'image':geomorphology_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/geomorphology_{i}.tif",'wb') as fd:
                fd.write(response.content)

            download = ee.data.getDownloadId({
                'image':coastal_slope_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/coastal_slope_{i}.tif",'wb') as fd:
                fd.write(response.content)
            

            download = ee.data.getDownloadId({
                'image':sea_level_rise_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/sea_level_rise_{i}.tif",'wb') as fd:
                fd.write(response.content)

            download = ee.data.getDownloadId({
                'image':shoreline_change_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/shoreline_change_{i}.tif",'wb') as fd:
                fd.write(response.content)

            download = ee.data.getDownloadId({
                'image':tidal_range_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/tidal_range_{i}.tif",'wb') as fd:
                fd.write(response.content)

            download = ee.data.getDownloadId({
                'image':wave_height_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/wave_height_{i}.tif",'wb') as fd:
                fd.write(response.content)
            
            progress_bar.setValue(60 + int(30*(i/grid_num.getInfo())))
            QApplication.processEvents()
        
        #merging raster data
        # create in folders
        geomorphology_infolder = f"{output_folder}/temp/geomorphology_*.tif"
        coastal_slope_infolder=f"{output_folder}/temp/coastal_slope_*.tif"
        sea_level_rise_infolder=f"{output_folder}/temp/sea_level_rise_*.tif"
        shoreline_change_infolder=f"{output_folder}/temp/shoreline_change_*.tif"
        tidal_range_infolder=f"{output_folder}/temp/tidal_range_*.tif"
        wave_height_infolder=f"{output_folder}/temp/tidal_range_*.tif"
        infolders=[geomorphology_infolder,coastal_slope_infolder,sea_level_rise_infolder,shoreline_change_infolder,tidal_range_infolder,wave_height_infolder]
        #create outfolder
        geomorphology_folder=f"{output_folder}/geomorphology.tif"
        coastal_slope_folder=f"{output_folder}/coastal slope.tif"
        sea_level_rise_folder=f"{output_folder}/Sea level rise.tif"
        shoreline_change_folder=f"{output_folder}/Shoreline change.tif"
        tidal_range_folder=f"{output_folder}/Tidal Range.tif"
        wave_height_folder=f"{output_folder}/Wave Height.tif"
        outfolders=[geomorphology_folder,coastal_slope_folder,sea_level_rise_folder,shoreline_change_folder,tidal_range_folder,wave_height_folder]

        progress_dialog.setLabelText("Merging...")
        for i in range(len(infolders)):
            MergeGrids(infolders[i],outfolders[i],1,"float32")
            addMapLayer(outfolders[i])
            progress_bar.setValue(90 + int(5*(i/len(infolders))))
            QApplication.processEvents()

        progress_dialog.setLabelText("Downloading image...")
        src_files = None
        shutil.rmtree(f"{output_folder}/temp")
        progress_bar.setValue(100)
    else:
        progress_dialog.close()
        show_error_message("Your region is not at the coast line!!")

        

def downloadCVI(output_folder,progress_dialog,progress_bar,region):
    initialize()
   
    progress_bar.setValue(0)
    progress_dialog.setLabelText("Initializing...")
    africa_regions = ee.FeatureCollection('projects/ee-kevinkiprotich0089/assets/MISLAND/Africa_Regions')
    shapefile = africa_regions.filter(ee.Filter.eq('NAME_1',region))

    progress_bar.setValue(10)
    progress_dialog.setLabelText("Loading coastal data...")
    geomorphology_CA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_CA')
    geomorphology_EA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_EA')
    geomorphology_NA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_NA')
    geomorphology_SA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_SA')
    geomorphology_WA = ee.Image('projects/ee-kevinkiprotich0089/assets/Geomorphology_WA')

    geomorphology = ee.ImageCollection([geomorphology_CA,geomorphology_EA,geomorphology_NA,geomorphology_SA,geomorphology_WA]).mosaic()
    referenceProjection=geomorphology.projection()
    geomorphology = geomorphology.resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    progress_bar.setValue(20)
    coastal_slope = ee.Image('projects/ee-kevinkiprotich0089/assets/coastal_slope').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    sea_level_rise = ee.Image('projects/ee-kevinkiprotich0089/assets/Sea_Level_Rise').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    shoreline_change = ee.Image('projects/ee-kevinkiprotich0089/assets/shoreline_change').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    tidal_range = ee.Image('projects/ee-kevinkiprotich0089/assets/tidal_range').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )
    wave_height = ee.Image('projects/ee-kevinkiprotich0089/assets/wave_height').resample('bicubic').reproject(
      crs= referenceProjection,
      scale= 100,
    )

    # cvi= ee.Image().expression(
    # 'sqrt((i1 * i2 * i3 * i4 * i5 * i6)/6)',
    #     {
    #         'i1': geomorphology,
    #         'i2': coastal_slope,
    #         'i3': sea_level_rise,
    #         'i4': shoreline_change,
    #         'i5': tidal_range,
    #         'i6': wave_height
    #     }
    # )

    cvi= ((geomorphology.multiply(coastal_slope).multiply(sea_level_rise).multiply(shoreline_change).multiply(tidal_range).multiply(wave_height)).divide(6)).sqrt()

    progress_bar.setValue(30)
    progress_dialog.setLabelText("Clipping data...")
    cvi_clipped = cvi.clip(shapefile.geometry()).selfMask()
    if rasterInRegion(cvi_clipped, shapefile):
        grid = geemap.fishnet(shapefile, h_interval=0.5, v_interval=0.5, delta=1)
        gridList = grid.toList(grid.size())
        grid_num = grid.toList(grid.size()).length()
        progress_bar.setValue(40)
        progress_dialog.setLabelText("Generating grid features...")
        # Create list of each grid feature

        ls_feature = []
        for i in range(grid_num.getInfo()):
            feature = ee.Feature(gridList.get(i)).geometry()
            ls_feature.append(feature)
            progress_bar.setValue(40 + int(20*(i/grid_num.getInfo())))
            QApplication.processEvents()

        # create a temporary folder to hold the grid datasets
        if not os.path.exists(output_folder+'/temp/'):
            os.makedirs(output_folder+'/temp/')
        
        progress_dialog.setLabelText("Downloading image...")
        for i in range(grid_num.getInfo()):
            log(f"Downloading image_{i}")
            download = ee.data.getDownloadId({
                'image':cvi_clipped,
                'bands':['b1'],
                'region':ls_feature[i],
                'scale':100,
                'format':'GEO_TIFF',
                'crs':'EPSG:4326'
            })

            response = requests.get(ee.data.makeDownloadUrl(download))
            with open (f"{output_folder}/temp/cvi_{i}.tif",'wb') as fd:
                fd.write(response.content)
            
            progress_bar.setValue(60 + int(30*(i/grid_num.getInfo())))
            QApplication.processEvents()

        cvi_infolder=f"{output_folder}/temp/cvi_*.tif"
        cvi_outfolder=f"{output_folder}/cvi.tif"
        MergeGrids(cvi_infolder,cvi_outfolder,1,"float32")
        progress_bar.setValue(95)
        addMapLayer(cvi_outfolder)

        src_files = None
        shutil.rmtree(f"{output_folder}/temp")
        progress_bar.setValue(100)
        progress_dialog.close()
    else:
        progress_dialog.close()
        show_error_message("Your region is not at the coast line!!")