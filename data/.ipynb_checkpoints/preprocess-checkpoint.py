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
    :param pad: #pixels to be padded around each tile, should be either one element or four elements
    :param overlap: #overlapping pixels between two patches in both vertical and horizontal direction
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


def rasterize_labels(labels, img_size):
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


def read_geotiff(geotiff_path):
    """Read geotiff, return reshaped image and metadata."""
    with rasterio.open(geotiff_path, 'r') as src:
        img = src.read()
        img_meta = src.meta

    return reshape_as_image(img), img_meta


def read_labels(labels_path, geotiff_crs):
    """Read geojson labels and convert projection, return geopandas dataframe."""
    labels = gpd.read_file(labels_path)
    labels = labels[labels.geometry.notnull()][labels.building == 'yes']

    return labels.to_crs({'init': geotiff_crs['init']})


def make_dir_if_not_exists(path, return_path=False):
    if not os.path.exists(path):
        os.makedirs(path)
    if return_path:
        return path


def save_image(img, path, name):
    make_dir_if_not_exists(path)
    data = Image.fromarray(img.astype(np.uint8))
    data.save(pathjoin(path, name))


if __name__ == '__main__':
    patch_size = (5000, 5000)
    train_folder = 'train_tier_1'
    locations = [f for f in os.listdir(train_folder) if f != 'patches']

    for i, loc in enumerate(locations):

        print(f'Processing location {i}/{len(locations)}')
        prefix = pathjoin(train_folder, loc)
        tiles = [f for f in os.listdir(prefix) if os.path.isdir(
            pathjoin(prefix, f)) and not f.endswith('labels')]

        for tile in tqdm(tiles):
            img_patches_dir = make_dir_if_not_exists(
                pathjoin(train_folder, 'patches', 'rgb'), return_path=True)
            gt_patches_dir = make_dir_if_not_exists(
                pathjoin(train_folder, 'patches', 'gt'), return_path=True)

            img, img_meta = read_geotiff(pathjoin(
                train_folder, loc, tile, f'{tile}.tif'))
            img_size = (img_meta['height'], img_meta['width'])
            labels = read_labels(
                pathjoin(train_folder, loc,
                         f'{tile}-labels', f'{tile}.geojson'),
                img_meta['crs'])
            gt = rasterize_labels(labels, img_size)

            for img_patch, gt_patch, y, x in patch_tile(img, gt, patch_size):
                img_patchname = f'{tile}_y{int(y)}x{int(x)}.png'
                gt_patchname = f'{tile}_y{int(y)}x{int(x)}.png'

                if len(np.unique(img_patch)) == 1 and np.unique(img_patch) == 0:
                    continue

                save_image(img_patch, pathjoin(
                    img_patches_dir, loc), img_patchname)
                save_image(gt_patch*255, pathjoin(
                    gt_patches_dir, loc), gt_patchname)
