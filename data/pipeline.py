import os
from os.path import join as pathjoin

import rasterio
from rasterio.plot import reshape_as_image
import rasterio.mask
from rasterio.features import rasterize

import pandas as pd
import geopandas as gpd
from shapely.geometry import mapping, Point, Polygon, MultiPolygon
from shapely.ops import cascaded_union

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

from tqdm import tqdm
from PIL import Image
import ee
import os
import time
from ee import batch



ee.Initialize()

def getNaipTask(coords,name):
    # takes as input coordinates of the boundary to export to 
    # google cloud bucket
        # takes as input coordinates of the boundary to export to 
    # google cloud bucket
    """
    Generate GEE Task for NAIP image
    :param coords: object coordinates
    :param name: object name
    :return: gee task to be executed
    """
    
    geom = ee.Geometry.Rectangle([coords[0],coords[1],coords[2],coords[3]]);
    collection = ee.ImageCollection("USDA/NAIP/DOQQ") \
                .filter(ee.Filter.date('2015-01-01', '2017-12-31'));
    trueColor = collection.select(['R', 'G', 'B','N'])
    trueColorVis = {
      min: 0.0,
      max: 255.0,
    }
    image = collection.sort('system:index', False).mosaic()
    image = image.clip(geom)
    image.projection()

    task = ee.batch.Export.image.toCloudStorage(image=image,
                                        region=image.geometry().bounds().\
                                        getInfo()['coordinates'],
                                        description='power_plant',
                                        outputBucket='earth_engine_data',
                                        fileNamePrefix=name,
                                        scale=10)
    return task

def getLandsatTask(coords,name):
    """
    Generate GEE Task for LANDSAT image
    :param coords: object coordinates
    :param name: object name
    :return: gee task to be executed
    """
    def maskL8sr(image):
        # Bits 3 and 5 are cloud shadow and cloud, respectively.
        cloudShadowBitMask = (1 << 3)
        cloudsBitMask = (1 << 5)
        # Get the pixel QA band.
        qa = image.select('pixel_qa')
        # Both flags should be set to zero, indicating clear conditions.
        a = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        b = qa.bitwiseAnd(cloudsBitMask).eq(0)
        mask = (a and b)
        return image.updateMask(mask).divide(10000)
    
    collection = ee.ImageCollection('LANDSAT/LC08/C01/T1_SR') \
                  .filterDate('2016-01-01', '2016-12-31').map(maskL8sr)
    
    bands = ['B4', 'B3', 'B2'] 
    collection = collection.select(bands)
    
    #geom = ee.Geometry.Rectangle([coords[0],coords[1],coords[2],coords[3]])
    geom = ee.Geometry.Rectangle([coords[0],coords[1],coords[2],coords[3]]);
    
    image = collection.sort('system:index', False).mosaic()
    image = image.clip(geom)
    image = image.select('B.+')
    
    imageRGB = image.visualize(bands=bands,min=0,max=3000,gamma=1.4)
    
    task = ee.batch.Export.image.toCloudStorage(scale=10,
                                                image=imageRGB,                                        
                                                description='power_plant',
                                                fileFormat = 'GeoTIFF',
                                                outputBucket='earth_engine_data',
                                                fileNamePrefix=name,
                                                region=image.geometry().bounds().getInfo()['coordinates']
                                                )
    return task

def getSarTask(coords,name):
    """
    Generate GEE Task for SAR image
    :param coords: object coordinates
    :param name: object name
    :return: gee task to be executed
    """
    imgVV = ee.ImageCollection('COPERNICUS/S1_GRD') \
            .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
            .filter(ee.Filter.eq('instrumentMode', 'IW')) \
            .select('VV') 
    
    desc = imgVV.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'));
    asc = imgVV.filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));
    
    spring = ee.Filter.date('2016-01-01', '2017-1-01');
    descChange = ee.Image.cat(desc.filter(spring).mean());
    ascChange = ee.Image.cat(asc.filter(spring).mean());
    
    geom = ee.Geometry.Rectangle([coords[0],coords[1],coords[2],coords[3]]);
    
    image = desc.mosaic()
    image = image.clip(geom)
    image = image.visualize(bands=['VV'],min=-25,max=5)
    
    task = ee.batch.Export.image.toCloudStorage(image=image,
                                        region=image.geometry().bounds().getInfo()['coordinates'],
                                        description='power_plant',
                                        outputBucket='earth_engine_data',
                                        fileNamePrefix=name,
                                        scale=10)
    return task


def getTask(coords, input_id_name, image_type):
    """
    Generate Task for general image
    :param coords: object coordinates
    :param name: object name
    :param image_type: type: NAIP, LANDSAT or SAR
    :return: gee task to be executed
    """
    if image_type == "NAIP":
        return getNaipTask(coords,input_id_name)
    elif image_type == "LANDSAT":
        return getLandsatTask(coords,input_id_name)
    elif image_type == "SAR":
        return getSarTask(coords,input_id_name)
    else:
        print("Incorrect input image type")
        return -1

# All the code for rasterizing image:
def crop_image(img, y, x, h, w):
    """
    Crop the image with given top-left anchor and corresponding width & height
    :param img: image to be cropped
    :param y: height of anchor
    :param x: width of anchor
    :param h: height of the patch
    :param w: width of the patch
    :return:
    """
    if len(img.shape) == 2:
        return img[y:y+w, x:x+h]
    else:
        return img[y:y+w, x:x+h, :]


def make_grid(tile_size, patch_size, overlap=0):
    """
    Extract patches at fixed locations. Output coordinates for Y,X as a list (not two lists)
    :param tile_size: size of the tile (input image)
    :param patch_size: size of the output patch
    :param overlap: #overlapping pixels
    :return:
    """
    max_h = int(tile_size[0] - patch_size[0])
    max_w = int(tile_size[1] - patch_size[1])

    if max_h > 0 and max_w > 0:
        h_step = int(np.ceil(tile_size[0] / (patch_size[0] - overlap)))
        w_step = int(np.ceil(tile_size[1] / (patch_size[1] - overlap)))
    else:
        h_step = 1
        w_step = 1
    patch_grid_h = np.floor(np.linspace(0, max_h, h_step)).astype(np.int32)
    patch_grid_w = np.floor(np.linspace(0, max_w, w_step)).astype(np.int32)

    y, x = np.meshgrid(patch_grid_h, patch_grid_w)

    return list(zip(y.flatten(), x.flatten()))


def patch_tile(rgb, gt, patch_size, pad=0, overlap=0):
    """
    Extract the given rgb and gt tiles into patches
    :param rgb:
    :param gt:
    :param patch_size: size of the patches, should be a tuple of (h, w)
    :param pad: #pixels to be padded around each tile, should be either 
    one element or four elements
    :param overlap: #overlapping pixels between two patches in both vertical
    and horizontal direction
    :return: rgb and gt patches as well as coordinates
    """
    # rgb = misc_utils.load_file(rgb_file)
    # gt = misc_utils.load_file(gt_file)[:, :, 0]
    np.testing.assert_array_equal(rgb.shape[:2], gt.shape)
    grid_list = make_grid(
        np.array(rgb.shape[:2]) + 2 * pad, patch_size, overlap)
    
    for y, x in grid_list:
        rgb_patch = crop_image(
            rgb, y, x, patch_size[0], patch_size[1])
        gt_patch = crop_image(
            gt, y, x, patch_size[0], patch_size[1])

        yield rgb_patch, gt_patch, y, x


def read_geotiff(geotiff_path):
    """Read geotiff, return reshaped image and metadata."""
    with rasterio.open(geotiff_path, 'r') as src:
        img = src.read()
        img_meta = src.meta
    return reshape_as_image(img), img_meta

def read_labels(labels_path, geotiff_crs, plant_id):
    """Read geojson labels and convert projection, return geopandas dataframe."""
    labels = gpd.read_file(labels_path)
    labels = labels[labels.geometry.notnull()]#[labels.building == 'yes']
    lb = labels.to_crs({'init': geotiff_crs['init']})
    lb_plant = lb.loc[lb['id'] == plant_id,:]
    return lb_plant

def make_dir_if_not_exists(path, return_path=False):
    if not os.path.exists(path):
        os.makedirs(path)
    if return_path:
        return path

def save_image(img, path, name):
    make_dir_if_not_exists(path)
    data = Image.fromarray(img.astype(np.uint8))
    data.save(pathjoin(path, name))

def rasterize_labels(labels, img_size,img_meta):
    """
    Draw rasterized labeled imagery based on corresponding geotiff image size.
    :param labels: geopandas dataframe, must have 'geometry' column with Polygon objects
    :img_size: corresponding geotiff image size
    """
    new_polygons = []

    for _, row in labels.iterrows():
        if isinstance(row['geometry'], Polygon):
            new_polygons.append(convert_polygon(
                row['geometry'], img_meta['transform']))
        elif isinstance(row['geometry'], MultiPolygon):
            for poly in list(row['geometry']):
                new_polygons.append(convert_polygon(
                    poly, img_meta['transform']))
        else:
            continue
    return rasterize(shapes=new_polygons, out_shape=img_size)

def convert_polygon(rowcol_polygon, transform):
    """
    Convert polygons from geojson rowcol coordinates to pixel positions
    :param rowcol_polygon: geojson polygon(s)
    :param transform: affine.Affine object, read from geotiff meta
    """
    polygon_points = []

    for point in np.array(rowcol_polygon.exterior.coords):
        # transform rowcol coords to geotiff crs, using reverse affine transformation
        polygon_points.append(~transform * point)

    return Polygon(polygon_points)

import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from OSMPythonTools.api import Api

'''This file contains the following functions:

1. get_ids(geojson_file_path)

Given a file_path, gets the ids for all the *ways* it contains
Each id will be an input to the next function:

2. get_coord(way_id)

Returns a list with the coordinates of the vertices of this way.
This function can be mapped on the output of the previous function
to return a list of the coordinates of each *way*.

The output of the get_coord function serves as an input to the next function:

3. bound_box(list_of_coordinates)

Constructs a box which bounds a *way* object and returns it as a Polygon object.

'''
def get_ids(geojson_file_path):
    '''Given a geojson file path extracted from OCM
    with different identified structures (*ways*),
    this function returns a list with their ids.'''

    #Read geojson file
    geojson = gpd.read_file(geojson_file_path)

    #Get the id colum and transform it in a list
    ids = list(geojson.loc[1:,'id'])
    names = list(geojson.loc[1:,'name'])
    dic = {}
    for i in range(len(ids)):
        if (ids[i].split('/')[0] == 'way'):
            dic[ids[i]] = names[i]
    #print(ids)
    #print(names)
    #Return list with ids
    return dic

def get_coord(uniq_id):
    '''Given the id of one identified structure (*way*),
    from OCM, this function returns a list with its coordinates'''

    #Define empty list to store coordinates
    coord = []

    #Execute query using OCM API
    #Requires command
    #from OSMPythonTools.api import Api
    api = Api()
    query_results = api.query(uniq_id)
    idType = uniq_id.split('/')[0]
    nodes = []
    if (idType == 'way'):
        nodes = query_results.nodes()
    elif(idType == 'relation'):
        nodes = query_results.members()
    #Get the coordinates for each node and store in list
    for node in nodes:
        coord.append(node.geometry()['coordinates'])

    return(coord)


def bound_box(list_of_coordinates):
    '''Given a list of coordinates
    (such as the output of function get_coords)
    returns a rectangle which contains the polygon
    defined by the coordinates.
    This rectangle is itself returned as a polygon object.
    Note that this is not the minimum bound box,
    but rather the minimum bound box which is parallel to the x and y axis (parallels and meridians).

    Requires
    from shapely.geometry import Polygon.'''

    #Collects x's and y's of polygons vertices
    x = []
    y = []
    for vertex in list_of_coordinates:
        x.append(vertex[0])
        y.append(vertex[1])

    #Define minimum and maximum values of x and y
    min_x = min(x)-0.01
    min_y = min(y)-0.01
    max_x = max(x)+0.01
    max_y = max(y)+0.01

    #Define vertices of box
    #Note that these are not point objects
    #as polygon objects in shapely are not created from point objects
    point_A = (min_x, min_y)
    point_B = (min_x, max_y)
    point_C = (max_x, max_y)
    point_D = (max_x, min_y)

    #Define bound box object as Polygon
    box = Polygon([point_A, point_B, point_C, point_D])

    #Return bound_box
    return(box)

def get_rast_image(filename, img_path, json_filepath, output_path, plant_id, image_type):
    """
    Generate rasterized image
    :param filename: input image name
    :param img_path: input image path
    :param json_filepath: geojson file path
    :param output_path: the path output to
    :param plant_id: id of the power plant
    :param image_type: type: NAIP, LANDSAT or SAR
    :return: 
    """
    if (image_type != "NAIP"):
        print("You cannot rasterize an image that is not type of NAIP")
        return -1
    patch_size = (500,500)
    # This will read the image and the meta data

    img,img_meta = read_geotiff(img_path)
    # get the image sizes
    img_size = (img_meta['height'], img_meta['width'])
    # readthe crs
    #labels = read_labels('jsons/{}.json'.format(filename),  img_meta['crs'])
    labels = read_labels(json_filepath, img_meta['crs'], plant_id)
    gt =rasterize_labels(labels, img_size,img_meta)
    img_patches_dir = output_path+'rgb'
    gt_patches_dir = output_path+'gt'
    loc='nia'
    idx=0
    for img_patch, gt_patch, y, x in patch_tile(img, gt, patch_size):
        idx+=1
        img_patchname='img-{}-{}.png'.format(filename,idx)
        gt_patchname='gt-{}-{}.png'.format(filename,idx)
        save_image(img_patch, img_patches_dir, img_patchname)
        save_image(gt_patch*255, gt_patches_dir, gt_patchname)
        break
        
        
def ggeToGoogleDr(allPlants, allPlantsGeojson, coords, input_id_name, image_type):
    """
    Export from Google Earth Engine to google storage
    :param allPlants: geojson file name
    :param allPlantsGeojson: full geojson file name
    :param coords: coordinates
    :param input_id_name: full id and name combination of the power plant
    :param image_type: type: NAIP, LANDSAT or SAR
    :return: 
    """
    #This geojson file contains the following *ways* in it:
    


    #The coordinates of the first *way* in the geojson file above is:
    
    #print(coords)

    #We can also get the coordinates of all the *ways* in the previous file:
    #all_coords = list(map(get_coord, ids))

    # These are the coordinates of the bounding box
    a = list(coords[0,:])
    b = list(coords[2,:])
    a.extend(b)
    task = getTask(a, input_id_name, image_type)
    if task == -1:
        return -1
    task.start()
    while (task.status()['state'] != 'COMPLETED'):
        time.sleep(1)
    return task.status()
    
# Step 2: download image from google cloud bucket to local directory
def googleDrToLocal(allPlants, gd_path, local_path):
    """
    Download from google storage to local machine
    :param allPlants: geojson file namee
    :param gd_path: google storage path
    :param local_path: local path
    :return: 
    """
    os.system('export GOOGLE_APPLICATION_CREDENTIALS=\"key.json\"')
    os.system('gsutil cp -r '+ gd_path + ' ' + local_path)