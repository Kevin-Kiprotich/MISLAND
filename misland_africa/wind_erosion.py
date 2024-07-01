from MISLAND.misland_africa.downloads import (show_error_message,
                                              initialize,
                                              addMapLayer,
                                              MergeGrids,
                                              rasterInRegion)
from qgis.PyQt.QtWidgets import QApplication,QFileDialog
from qgis.core import QgsProject, QgsRasterLayer,QgsVectorLayer,QgsProcessing
from qgis.utils import iface
from MISLAND import log
import processing
import requests
import rasterio
import geemap
import shutil
import glob
import ee
import os


# Linear Fuzzification
def fuzzy_linear(img, lowBound, highBound) :  
    if lowBound < highBound:
        maskedImage = img.multiply(0) \
                         .where(img.lte(lowBound), 0) \
                         .where(img.gte(highBound),1) \
                         .where(img.gt(lowBound)and(img.lt(highBound)), img )
                            
        inverted_img = maskedImage.expression('( (img - lowBound )/(highBound - lowBound ) )', {
                            'img': maskedImage,
                            'lowBound': lowBound,
                            'highBound': highBound
                            })


    elif lowBound > highBound: 
        maskedImage = img.multiply(0)\
                                .where(img.lte(highBound), 1)\
                                .where(img.gte(lowBound), 0)\
                                .where(img.gt(highBound)and(img.lt(lowBound)), img )#(raster_array > high_bound) & (raster_array < low_bound))

        inverted_img = maskedImage.expression('( (img - lowBound )/(highBound - lowBound ) )', {
                            'img': maskedImage,
                            'lowBound': lowBound,
                            'highBound': highBound
                            })


    return inverted_img

def fuzzify_exponetial(img, lowBound, highBound):
  
    if(lowBound < highBound):
    # Zero_mask = img.multiply(0).where(img.lte(lowBound), img)
    # one_mask = img.multiply(1).where(img.gte(highBound), img)

        maskedImg = img.multiply(0)\
                                    .where(img.lte(lowBound), 0)\
                                    .where(img.gte(highBound),1)\
                                    .where(img.gt(lowBound)and(img.lt(highBound)), img)
                                    # .rename('img')

        f_inverted = maskedImg.expression('( (img - lowBound )/(highBound - lowBound ) ) **2', {
                            'img': maskedImg,#.select('img'),
                            'lowBound': lowBound,
                            'highBound': highBound
                            })

    elif(lowBound > highBound):

        maskedImg = img.multiply(0)\
                                    .where(img.lte(highBound), 1)\
                                    .where(img.gte(lowBound), 0)\
                                    .where(img.gt(highBound)and(img.lt(lowBound)), img)#(raster_array > high_bound) & (raster_array < low_bound)
                                    # .rename('img')

        f_inverted = maskedImg.expression('( (img - lowBound )/(highBound - lowBound ) ) **2', {
                            'img': maskedImg,#.select('img'),
                            'lowBound': lowBound,
                            'highBound': highBound
                            })


    return f_inverted

def getILSWE(output_folder,progress_dialog,progress_bar,country,region,year):
    initialize()
    log("Fetching data. Please wait!!!...")
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"
    filename = f"ILSWE__{country}_{region}_{year}"

    progress_bar.setValue(0)
    progress_dialog.setLabelText("Initializing ...")

    africa_regions = ee.FeatureCollection('projects/ee-kevinkiprotich0089/assets/MISLAND/Africa_Regions')
    shapefile = africa_regions.filter(ee.Filter.eq('NAME_1',region))

    ##############################################################################
    #===== CLIMATE EROSIVITY (CE) ================================================
    ##############################################################################
    log("Computing R-Factor ...")
    progress_dialog.setLabelText("Computing R-Factor ...")
    progress_bar.setValue(10)
    climate_dataset = (ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
                        .filter(ee.Filter.date(start_date, end_date)))
    windSpeed= climate_dataset.select('vs').mean()
    windSpeed=windSpeed.multiply(0.01).rename('wind')

    petRescaled = climate_dataset.select('pet').mean()
    petRescaled = petRescaled.multiply(0.1).rename('PET')

    rainfall = climate_dataset.select('pr').mean().rename('rainfall')

    ce_data = windSpeed.addBands(petRescaled).addBands(rainfall)
    ce_data = ce_data.clip(shapefile)

    CE = ce_data.expression('0.01*(U**3)*((PET - P)/PET)', {
        'U': ce_data.select('wind'),
        'PET': ce_data.select('PET'),
        'P': ce_data.select('rainfall')
    }).rename('CE')

    if not rasterInRegion(CE,shapefile):
        show_error_message(f"Could not find Climate Erosivity for {region}.")
        return {}

    ##############################################################################
    #===== Vegetation Cover (VC) =================================================
    ##############################################################################

    log("Computing R-Factor ...")
    progress_dialog.setLabelText("Computing Fractional Vegetation Cover () ...")
    progress_bar.setValue(20)

    ndvi = (ee.ImageCollection('MODIS/006/MOD13Q1')
                  .filter(ee.Filter.date(start_date, end_date))
                  .select(['NDVI'], ['ndvi'])
                  .mean().divide(10000)
                  .clip(shapefile))
    
    VC = ndvi.expression('(ndvi - 0.05)/( 0.86 - 0.05)',{
            'ndvi': ndvi.select('ndvi'),
        }).rename('F_cover')
    
    if not rasterInRegion(VC,shapefile):
        show_error_message(f"Could not find Vegetation Cover (VC) for {region}.")
        return {}
    
    ##############################################################################
    #===== Soil Erodible Factor (EF) =============================================
    ##############################################################################
    log("Computing Fractional Vegetation Cover (VC) ...")
    progress_dialog.setLabelText("Computing Fractional Vegetation Cover (VC) ...")
    progress_bar.setValue(30)
    clay_content = ee.Image("ISDASOIL/Africa/v1/clay_content").select('mean_0_20').clip(shapefile).rename('cl')
    sand_content = ee.Image("ISDASOIL/Africa/v1/sand_content").select('mean_0_20').clip(shapefile).rename('sa')
    silt_content = ee.Image("ISDASOIL/Africa/v1/silt_content").select('mean_0_20').clip(shapefile).rename('si')
    organic_matter = ee.Image("ISDASOIL/Africa/v1/carbon_organic").select('mean_0_20').clip(shapefile).rename('om')
    calcium_carbonate = ee.Image("ISDASOIL/Africa/v1/calcium_extractable").select('mean_0_20').clip(shapefile).rename('Ca')
    soilData = clay_content.addBands(sand_content).addBands(silt_content).addBands(organic_matter).addBands(calcium_carbonate)
    EF = soilData.expression("(29.09 + (0.31*Sa) + (0.17 * Si) + (0.33 * (Sa/Cl)) - (2.59 * OM) - (0.95 * CaCO3))/100", {
          'Sa': soilData.select('sa'),
          'Si': soilData.select('si'),
          'Cl': soilData.select('cl'),
          'OM': soilData.select('om'),
          'CaCO3': soilData.select('Ca')
        }).rename('EF')
    if not rasterInRegion(EF,shapefile):
        show_error_message(f"Could not find Soil Erodible Factor for {region}.")
        return {}
    ##############################################################################
    #===== Soil Crust Factor (SCF) ===============================================
    ##############################################################################
    log("Computing Soil Crust Factor (SCF) ...")
    progress_dialog.setLabelText("Computing Soil Crust Factor (SCF) ...")
    progress_bar.setValue(40)
    SC = soilData.expression("1/(1 + (0.0066) * (Cl**2) + (0.21) * (OM**2) )", {
            'Cl': soilData.select('cl'),
            'OM': soilData.select('om')
        }).rename('SC')
    
    if not rasterInRegion(SC,shapefile):
        show_error_message(f"Could not find Soil Crust Factor for {region}.")
        return {}
    ##############################################################################
    #===== Soil Roughness (SR) ===================================================
    ##############################################################################
    log("Computing Surface Roughness (SR) ...")
    progress_dialog.setLabelText("Computing Surface Roughness (SR) ...")
    progress_bar.setValue(45)
    SR = ee.Image('users/mutindarisper/gwa3_250_RIX').clip(shapefile).rename('SR')

    if not rasterInRegion(SR,shapefile):
        show_error_message(f"Could not find Surface roughness image for {region}.")
        return {}
        
    #Fuzzifying features

    log("Fuzzifying ...")
    progress_dialog.setLabelText("Fuzzifying ...")
    progress_bar.setValue(50)
    CE_Fuzzy = fuzzy_linear(CE, 1, 0).rename('CE')
    progress_bar.setValue(52)
    EF_Fuzzy = fuzzy_linear(EF, 1, 0).rename('EF')
    progress_bar.setValue(54)
    SC_Fuzzy = fuzzy_linear(SC, 1, 0).rename('SC')
    progress_bar.setValue(56)
    VC_Fuzzy = fuzzify_exponetial(VC, 1, 0).rename('VC')
    progress_bar.setValue(58)

    log("Computing ILSWE ...")
    progress_dialog.setLabelText("Computing ILSWE ...")
    progress_bar.setValue(60)
    ILSWE = CE_Fuzzy.multiply(EF_Fuzzy).multiply(SR).multiply(SC_Fuzzy).multiply(VC_Fuzzy).rename('ILSWE')

    progress_dialog.setLabelText("Generating grid ...")
    grid = geemap.fishnet(shapefile, h_interval=0.5, v_interval=0.5, delta=1)
    gridList = grid.toList(grid.size())
    grid_num = grid.toList(grid.size()).length()

    ls_feature = []
    for i in range(grid_num.getInfo()):
        feature = ee.Feature(gridList.get(i)).geometry()
        ls_feature.append(feature)
        progress_bar.setValue(60 + int(10*(i/grid_num.getInfo())))
        QApplication.processEvents()

    if not os.path.exists(output_folder+'/temp/'):
            os.makedirs(output_folder+'/temp/')
    
    log("Downloading ILSWE ...")
    progress_dialog.setLabelText("Downloading ILSWE ...")

    for i in range(grid_num.getInfo()):
        log(f"Downloading image_{i}")
        download = ee.data.getDownloadId({
            'image':ILSWE,
            'bands':['ILSWE'],
            'region':ls_feature[i],
            'scale':100,
            'format':'GEO_TIFF',
            'crs':'EPSG:4326'
        })

        response = requests.get(ee.data.makeDownloadUrl(download))
        with open (f"{output_folder}/temp/{filename}_{i}.tif",'wb') as fd:
            fd.write(response.content)
        
        progress_bar.setValue(70 + int(20*(i/grid_num.getInfo())))
        QApplication.processEvents()

    infolder=f"{output_folder}/temp/{filename}_*.tif"
    outfolder=f"{output_folder}/temp/{filename}.tif"
    MergeGrids(infolder,outfolder,1,"float32",-9999)
    progress_bar.setValue(95)
    Add_and_ClipLayer(outfolder,shapefile,country,region,output_folder,year)

    shutil.rmtree(f"{output_folder}/temp")
    progress_bar.setValue(100)
    progress_dialog.close()

def Add_and_ClipLayer(path,shapefile,country,region,folder,year):
    log("Adding layer to QGIS")
    project = QgsProject.instance()
    filename=path.split('/')[-1][0:-4]
    shapefile_path = f"{folder}/temp/{region}.shp"
    geemap.ee_to_shp(shapefile,filename=shapefile_path)

    layer=QgsRasterLayer(path,filename,'gdal')
    if not layer.isValid():
        show_error_message("Invalid raster file")
        return

    vector = QgsVectorLayer(shapefile_path, "vector_layer", "ogr")
    if not vector.isValid():
        show_error_message("Invalid clip shapefile..")
        return 
    
    params= {
        'INPUT':layer,
        'MASK':vector,
        'OUTPUT':f"{folder}/ILSWE_{country}_{region}_{year}.tif"
    }

    result = processing.run("gdal:cliprasterbymasklayer", params)
    print(result)
    clipped = QgsRasterLayer(result['OUTPUT'],filename,'gdal')
    if not clipped.isValid():
        show_error_message("Invalid clip output")
        return
    project.instance().addMapLayer(clipped)

    canvas = iface.mapCanvas()
    extent = clipped.extent()
    canvas.setExtent(extent)
    canvas.refresh()
    return result