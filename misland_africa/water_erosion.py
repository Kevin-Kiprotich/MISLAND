from MISLAND.misland_africa.downloads import (show_error_message,
                                              initialize,
                                              addMapLayer,
                                              MergeGrids,
                                              rasterInRegion)
from qgis.PyQt.QtWidgets import QApplication
from MISLAND import log
import requests
import geemap
import shutil
import ee
import os




"""
    Main function to compute the RUSLE model
"""
def getRUSLE(output_folder,progress_dialog,progress_bar,country,region,year):
    initialize()
    log("Fetching data. Please wait!!!...")
    start_date = f"{year}-01-01"
    end_date = f"{year}-12-31"
    log(f"RUSLE:\t{year}")

    filename = f"Soil_loss_{country}_{region}_{year}"
    print(filename)
    progress_bar.setValue(0)
    progress_dialog.setLabelText("Initializing ...")

    #Get the area of interest
    africa_regions = ee.FeatureCollection('projects/ee-kevinkiprotich0089/assets/MISLAND/Africa_Regions')
    shapefile = africa_regions.filter(ee.Filter.eq('NAME_1',region))

    def clipCollection(img):
        return img.clip(shapefile)

    def monthlySquare(img):
        return img.expression('Pi**2',{
            'Pi': img.select('pr')
        }).rename('Pi')

    

    # get the relevant datasets
    progress_bar.setValue(10)
    ##############################################################################################################
    #==================================== R-Factor (Rainfall Erosivity) ==========================================
    ##############################################################################################################
    log("Computing R-Factor ...")
    progress_dialog.setLabelText("Computing R-Factor ...")
    climateData = (ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filter(ee.Filter.date(start_date, end_date)))
    precipitation = climateData.select('pr').map(clipCollection)
    monthly_P = precipitation.map(monthlySquare)
    annualTotal = precipitation.reduce(ee.Reducer.sum()).rename('annualTotal')

    def computeMFI(img):
        return img.expression('Pi/P',{
            "Pi":img.select('Pi'),
            'P':annualTotal.select('annualTotal')
        })
    
    monthly_MFI = monthly_P.map(computeMFI)
    R_Factor = monthly_MFI.reduce(ee.Reducer.sum()).rename('R')
    print(f"R-Factor:\t{R_Factor.getInfo()}") 

    ############################################################################################################
    #==================================== S-Factor (Slope Steepeness) ==========================================
    ############################################################################################################
    progress_bar.setValue(20)
    log("Computing S-Factor ...")
    progress_dialog.setLabelText("Computing S-Factor ...")
    elevation = ee.Image('CGIAR/SRTM90_V4').clip(shapefile)
    slope = ee.Terrain.slope(elevation).rename('slope')

    S_Factor_lt = slope.expression('(10.8 * sin(A)) + 0.03',{
            'A': slope.select('slope')
        }).rename('S')
    
    S_Factor_gt = slope.expression('(16.8 * sin(A)) - 0.05',{
            'A': slope.select('slope')
        }).rename('S')
    
    S_Factor = slope.multiply(0).where(slope.lt(5.15), S_Factor_lt).where(slope.gte(5.15), S_Factor_gt).rename('S')
    print(f"S-Factor:\t{S_Factor.getInfo()}")
    #################################################################################################################
    #==================================== P-Factor (Conservation practice) ==========================================
    #################################################################################################################

    progress_bar.setValue(30)
    log("Computing P-Factor ...")
    progress_dialog.setLabelText("Computing P-Factor ...")
    landcover = ee.Image("COPERNICUS/Landcover/100m/Proba-V-C3/Global/2015").clip(shapefile)
    cropland = landcover.select('discrete_classification').eq(40)
    
    #convert slope to percentage
    slope_pc = slope.expression('(tan(s * 0.0174533)) * 100',{
            's': slope.select('slope')
        }).rename('slope_pc')
    
    # Reclassify the % slope
    slope_reclassed = ee.Image([0, 3, 9, 13, 17, 21, 25])
    slopeReclassed = slope_pc.gt(slope_reclassed).reduce(ee.Reducer.sum()).toInt().select('sum')
    print(f"Slope-Reclassed:\t{slopeReclassed.getInfo()}")
    P_Factor = cropland.multiply(0) \
                       .where(cropland.eq(0), 1) \
                       .where(cropland.eq(1).And(slopeReclassed.eq(1)), 0.6) \
                       .where(cropland.eq(1).And(slopeReclassed.eq(2)), 0.5) \
                       .where(cropland.eq(1).And(slopeReclassed.eq(3)), 0.6) \
                       .where(cropland.eq(1).And(slopeReclassed.eq(4)), 0.7) \
                       .where(cropland.eq(1).And(slopeReclassed.eq(5)), 0.8) \
                       .rename('P')
    
    print(f"P-Factor:\t{P_Factor.getInfo()}")
    log("P-Factor complete ...")

    
    ##############################################################
    #=====K-Factor (Soil Erodibility Factor) =====================
    ##############################################################
    progress_bar.setValue(40)
    log("Computing K-Factor ...")
    progress_dialog.setLabelText("Computing K-Factor ...")

    clay = ee.Image("OpenLandMap/SOL/SOL_CLAY-WFRACTION_USDA-3A1A1A_M/v02").select('b30').clip(shapefile)
    sand = ee.Image("OpenLandMap/SOL/SOL_SAND-WFRACTION_USDA-3A1A1A_M/v02").select('b30').clip(shapefile)
    org = ee.Image("OpenLandMap/SOL/SOL_ORGANIC-CARBON_USDA-6A1C_M/v02").select('b30').clip(shapefile).multiply(0.5)
    silt = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").select('b30').clip(shapefile)
    soil_texture = ee.Image("OpenLandMap/SOL/SOL_TEXTURE-CLASS_USDA-TT_M/v02").select('b30').clip(shapefile)

    #soil permiability classes
    permiability = soil_texture.multiply(0) \
                              .where(soil_texture.eq(12), 1) \
                              .where(soil_texture.eq(11).Or(soil_texture.eq(9)), 2) \
                              .where(soil_texture.eq(7).Or(soil_texture.eq(8)).Or(soil_texture.eq(10)), 3) \
                              .where(soil_texture.eq(6).Or(soil_texture.eq(4)), 4) \
                              .where(soil_texture.eq(5).Or(soil_texture.eq(3)), 5) \
                              .where(soil_texture.eq(2).Or(soil_texture.eq(1)), 6)
    
    # soil textural groups(Structure) form soil textural class reference https://passel2.unl.edu/view/lesson/0cff7943f577/2
    soil_structure = soil_texture.multiply(0) \
                              .where(soil_texture.eq(12).Or(soil_texture.eq(11)), 1) \
                              .where(soil_texture.eq(9), 2) \
                              .where(soil_texture.eq(7).Or(soil_texture.eq(8)).Or(soil_texture.eq(10)), 3) \
                              .where(soil_texture.eq(6).Or(soil_texture.eq(5)).Or(soil_texture.eq(4)).Or(soil_texture.eq(3)).Or(soil_texture.eq(2)).Or(soil_texture.eq(1)), 4)

    soil_bands = sand.addBands(clay).addBands(org).addBands(silt).addBands(permiability).addBands(soil_structure)
    bands_named = soil_bands.rename(['sand','clay','org','silt','permiability','structure'])

    particleSize = bands_named.expression('Silt * (100 - clay)',{
            'Silt': bands_named.select('silt'),
            'clay': bands_named.select('clay')
        }).rename('Ps')
    
    bands_named = bands_named.addBands(particleSize)

    K_Factor = bands_named.expression('((2.1 * 0.0001 * (M ** 1.14) * (12 - OM)) + (3.25 * (s - 2)) + (2.5 * (p -3))) * 0.1317', {
            'M':bands_named.select('Ps'),
            'OM':bands_named.select('org'),
            'p' : bands_named.select('permiability'),
            's' : bands_named.select('structure')
            }).rename ('K')
    
    log("K-Factor complete")
    print(f"K-Factor:\t{K_Factor.getInfo()}")
    ##############################################################
    #=====C-Factor (Soil Erodibility Factor) =====================
    ##############################################################
    progress_bar.setValue(50)
    log("Computing C-Factor ...")
    progress_dialog.setLabelText("Computing C-Factor ...")
    landcoverType = landcover.select('discrete_classification')
    non_crop_MAXc = landcoverType.multiply(0) \
                           .where(landcoverType.eq(40), 0.38) \
                           .where(landcoverType.eq(111).Or(landcoverType.eq(112)).Or(landcoverType.eq(113)).Or(landcoverType.eq(114)) \
                           .Or(landcoverType.eq(115)).Or(landcoverType.eq(116)).Or(landcoverType.eq(121)).Or(landcoverType.eq(122)) \
                           .Or(landcoverType.eq(123)).Or(landcoverType.eq(124)).Or(landcoverType.eq(125)).Or(landcoverType.eq(126)), 0.003) \
                           .where(landcoverType.eq(20).Or(landcoverType.eq(30)).Or(landcoverType.eq(90)), 0.15) \
                           .where(landcoverType.eq(60).Or(landcoverType.eq(50)), 0.5) \
                           .where(landcoverType.eq(200), 0) \
                           .rename('MAXc')
    
    non_crop_MINc = landcoverType.multiply(0) \
                               .where(landcoverType.eq(40), 0.15) \
                               .where(landcoverType.eq(111).Or(landcoverType.eq(112)).Or(landcoverType.eq(113)).Or(landcoverType.eq(114)) \
                               .Or(landcoverType.eq(115)).Or(landcoverType.eq(116)).Or(landcoverType.eq(121)).Or(landcoverType.eq(122)) \
                               .Or(landcoverType.eq(123)).Or(landcoverType.eq(124)).Or(landcoverType.eq(125)).Or(landcoverType.eq(126)), 0.0001) \
                               .where(landcoverType.eq(20).Or(landcoverType.eq(30)).Or(landcoverType.eq(90)), 0.01) \
                               .where(landcoverType.eq(60).Or(landcoverType.eq(50)), 0.1) \
                               .where(landcoverType.eq(200), 0) \
                               .rename('MINc')
    
    NDVI = (ee.ImageCollection('MODIS/006/MOD13Q1')
                  .filter(ee.Filter.date(start_date, end_date))
                  .select('NDVI')
                  .mean().divide(10000)
                  .clip(shapefile))
    
    F_cover =  NDVI.expression('(NDVI - 0.05)/( 0.86 - 0.05)',{
            'NDVI': NDVI.select('NDVI'),
        }).rename('F_cover')
    
    F_cover_factor = non_crop_MAXc.addBands(non_crop_MINc).addBands(F_cover)

    C_Factor = F_cover_factor.expression('MINc + (MAXc - MINc) * (1 - Fc)',{
            'MINc': F_cover_factor.select('MINc'),
            'MAXc': F_cover_factor.select('MAXc'),
            'Fc': F_cover_factor.select('F_cover')
        }).rename('C')
    
    log("C-Factor complete")
    print(f"C-Factor:\t{C_Factor.getInfo()}")
    ##############################################################
    #=====RUSLE (Soil Loss) ======================================
    ##############################################################
    progress_bar.setValue(60)
    log("Computing RUSLE soil loss ...")
    progress_dialog.setLabelText("Computing RUSLE soil loss ...")
    landcoverType = landcover.select('discrete_classification')
    RUSLE_bands = R_Factor.addBands(S_Factor).addBands(P_Factor).addBands(K_Factor).addBands(C_Factor)
    rusle = RUSLE_bands.expression('R * S * P * K * C', {
            'R': RUSLE_bands.select('R'),
            'S': RUSLE_bands.select('S'),
            'P': RUSLE_bands.select('P'),
            'K': RUSLE_bands.select('K'),
            'C': RUSLE_bands.select('C'),
        }).rename('RUSLE')
    print(rusle.getInfo())

    rusle_thresholds = ee.Image([2,10,20,50,1000])
    soil_loss = rusle.lt(rusle_thresholds).reduce(ee.Reducer.sum()).toInt().select('sum')

    print(f"Soil loss:\t{soil_loss.getInfo()}")
    progress_bar.setValue(64)
    progress_dialog.setLabelText("Generating grid ...")
    grid = geemap.fishnet(shapefile, h_interval=0.5, v_interval=0.5, delta=1)
    gridList = grid.toList(grid.size())
    grid_num = grid.toList(grid.size()).length()

    ls_feature = []
    for i in range(grid_num.getInfo()):
        feature = ee.Feature(gridList.get(i)).geometry()
        ls_feature.append(feature)
        progress_bar.setValue(64 + int(6*(i/grid_num.getInfo())))
        QApplication.processEvents()

    if not os.path.exists(output_folder+'/temp/'):
            os.makedirs(output_folder+'/temp/')

    progress_bar.setValue(71)
    log("Downloading RUSLE soil loss ...")
    progress_dialog.setLabelText("Downloading RUSLE soil loss ...")

    for i in range(grid_num.getInfo()):
        log(f"Downloading image_{i}")
        download = ee.data.getDownloadId({
            'image':soil_loss,
            'bands':['sum'],
            'region':ls_feature[i],
            'scale':100,
            'format':'GEO_TIFF',
            'crs':'EPSG:4326'
        })

        response = requests.get(ee.data.makeDownloadUrl(download))
        with open (f"{output_folder}/temp/{filename}_{i}.tif",'wb') as fd:
            fd.write(response.content)
        
        progress_bar.setValue(71 + int(20*(i/grid_num.getInfo())))
        QApplication.processEvents()

    infolder=f"{output_folder}/temp/{filename}_*.tif"
    outfolder=f"{output_folder}/{filename}.tif"
    MergeGrids(infolder,outfolder,1,"int32",-2147483648)
    progress_bar.setValue(95)
    addMapLayer(outfolder)

    shutil.rmtree(f"{output_folder}/temp")
    progress_bar.setValue(100)
    progress_dialog.close()
