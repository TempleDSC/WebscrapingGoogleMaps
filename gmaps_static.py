#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 15:13:57 2019

@author: tug76662
"""
import time
import math
#import urllib
#import pandas as pd
import shapely
#import earthpy
import geopandas as gpd
#import rasterio as rio
#from rasterio.plot import plotting_extent
#import earthpy.spatial as es
#from matplotlib import pyplot
#from obspy.geodetics import kilometers2degrees
import numpy as np
import gdal
import osr
#import geopy
import LatLon.LatLon as ll
from PIL import Image
import requests
#import geojson

GMAPS_BASE_URL = 'https://maps.googleapis.com/maps/api/staticmap?'
GMAPS_API_KEY = 'KEY'
GMAPS_SIZE_PX_X = 256
GMAPS_SIZE_PX_Y = GMAPS_SIZE_PX_X + 20
GMAPS_ZOOM = 20
GMAP_SCALE_20 = 1128.497220
DATA_DIR = '/Users/tug76662/Documents/census processing/'



#NOTES - 20m earth, 15m on qgis, 5 m off -- why?


def georeference_png_and_save_geotiff(npimg, centerpoint, tileszkm, filename, row):
    xmin = centerpoint.offset(270,tileszkm/2.0).lon.decimal_degree
    xmax = centerpoint.offset(90,tileszkm/2.0).lon.decimal_degree
    xres = math.fabs(xmin - xmax) / GMAPS_SIZE_PX_X
    ymin = centerpoint.offset(0,tileszkm/2.0).lat.decimal_degree
    ymax = centerpoint.offset(180,tileszkm/2.0).lat.decimal_degree
    yres = math.fabs(ymin - ymax) / GMAPS_SIZE_PX_X
    ulx = xmin - (xres / 2.0)
    uly = ymax - (yres / 2.0)
    
    gt = [ulx, xres, 0, uly, 0, yres ]
    
    driver = gdal.GetDriverByName('GTiff')
    options = ['PHOTOMETRIC=RGB', 'PROFILE=GeoTIFF']
    ds = driver.Create(DATA_DIR+filename, GMAPS_SIZE_PX_X, GMAPS_SIZE_PX_X, 3, gdal.GDT_Float32, options=options)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    #reference: https://stackoverflow.com/questions/16095699/python-gdal-georeference-array-using-other-file-for-projection
    #           https://stackoverflow.com/questions/27166739/description-of-parameters-of-gdal-setgeotransform
    ds.SetGeoTransform(gt)
    
    outband_r = ds.GetRasterBand(1)
    outband_r.WriteArray(npimg[:,:,0])
    outband_g = ds.GetRasterBand(2)
    outband_g.WriteArray(npimg[:,:,1])
    outband_b = ds.GetRasterBand(3)
    outband_b.WriteArray(npimg[:,:,2])
    
    ds = None
#    # NOW CROP AGAINST SHAPEFILE  
    # SKIP THIS STEP FOR NOW AS IT'S UNDESIRED FOR THE ML PROCESS
#  https://rasterio.readthedocs.io/en/stable/topics/masking-by-shapefile.html
#    crop = None
#    crop_meta = None
#    with rio.open(DATA_DIR+filename) as src:
#        df = pd.DataFrame()
#        gdf = gpd.GeoDataFrame(df,geometry=gpd.points_from_xy(row.geometry.exterior.coords.xy[0], row.geometry.exterior.coords.xy[1]))
#        crop, crop_meta = es.crop_image(src,gdf)
#        pyplot.imshow(crop)
#    with rio.open(DATA_DIR+filename, 'w', **crop_meta) as ff:
#        ff.write(crop, 3)
    
    #end georeference_png_and_save_geotiff
        

#crop_extent_phl = gpd.read_file(DATA_DIR+'/Census_Tracts_2010/Census_Tracts_2010.shp')
crop_extent = gpd.read_file(DATA_DIR+'/500 Cities/500Cities_Tracts_Clip.shp')
allentown = crop_extent['PlaceName'].str.match('Allentown')
crop_extent = crop_extent[allentown]
# MAKE SURE WE ARE IN WGS84 !
crop_extent = crop_extent.to_crs({'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})

#tile_lon = crop_extent.iloc[0]['geometry'].bounds[0]
#tile_lat = crop_extent.iloc[0]['geometry'].bounds[1]

# FOR EVERY CENSUS TRACT...
for index, row in crop_extent.iterrows():
    bounding_box = row['geometry'].bounds
    tile_lon = bounding_box[0]

    slice_num = 0
    while tile_lon < bounding_box[2]:
        tile_lat = bounding_box[1]
        while tile_lat < bounding_box[3]:
            centerpoint = ll.LatLon(tile_lat, tile_lon)
            # CALCULATE EXTENTS WE'LL NEED TO GEOREFERENCE AND FIND NEXT TILE POS
            # notes:  resolution in degrees per pixel ...  resolution = 156543.03 meters/pixel * cos(latitude) / (2 ** zoomlevel)
            # NOTE: LATITUDE IN RADIANS NOT DEGREES
            resmpp = math.fabs(156543.03  * math.cos(math.radians(tile_lat)) / (2.0 ** GMAPS_ZOOM))
            tileszmeters = resmpp * GMAPS_SIZE_PX_X
            tileszkm = tileszmeters/1000.0
            xmin = centerpoint.offset(270.0,tileszkm/2.0).lon.decimal_degree
            xmax = centerpoint.offset(90.0,tileszkm/2.0).lon.decimal_degree
            ymin = centerpoint.offset(0.0,tileszkm/2.0).lat.decimal_degree
            ymax = centerpoint.offset(180.0,tileszkm/2.0).lat.decimal_degree
            tile = shapely.geometry.polygon.Polygon([(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)])
            # IGNORE TILES THAT FALL OUTSIDE CENSUS TRACT
            if row['geometry'].intersects(tile):
                GMAPS_PARAMS = 'center={},{}&zoom={}&size={}x{}&maptype=satellite'.format(tile_lat,tile_lon,GMAPS_ZOOM,GMAPS_SIZE_PX_X,GMAPS_SIZE_PX_Y)
                url = GMAPS_BASE_URL + GMAPS_PARAMS + '&key=' + GMAPS_API_KEY
                # GET GMAPS TILE
                time.sleep(5)
                response = requests.get(url, stream=True)
                response.raw.decode_content = True
                image = Image.open(response.raw)
                rgbimg = image.convert('RGB')
                # CROP OUT WATERMARK AND FLIP PIXEL ARRAY
                npimg = np.flipud(np.array(rgbimg)[:GMAPS_SIZE_PX_X,:])
                # CREATE GEOTIFF
                # TRY TO WORK AROUND VARYING METADATA IN DIFFERENT SHAPEFILES
                prefix = ''
                if 'STATEFP10' in row:
                    prefix = "{}{}_{}".format(row['STATEFP10'],row['COUNTYFP10'],row['TRACTCE10'])
                elif 'PlaceName' in row:
                    prefix = "{}_{}".format(row['tract2010'],row['PlaceName'])
                filename = "{}_{}_{}_{}.tif".format(prefix,slice_num,tile_lat,tile_lon)
                georeference_png_and_save_geotiff(npimg, centerpoint, tileszkm, filename, row)
            else:
                print("Skipping tile {} as out of bounds".format(slice_num))
            # MOVE TO NEXT CROP POINT
            tile_lat = centerpoint.offset(0,tileszkm).lat.decimal_degree
            slice_num += 1
            # end lat loop
        tile_lon = centerpoint.offset(90,tileszkm).lon.decimal_degree
        # end lon loop
    # end tract loop
