import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from OSMPythonTools.api import Api

'''
Author: Felipe Buchbinder

This file contains the following functions:

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
    ids = list(geojson.loc[:,'id'])

    #Return list with ids
    return(ids)


def get_coord(way_id):
    '''Given the id of one identified structure (*way*),
    from OCM, this function returns a list with its coordinates'''

    #Define empty list to store coordinates
    coord = []

    #Execute query using OCM API
    #Requires command
    #from OSMPythonTools.api import Api
    api = Api()
    query_results = api.query("way/617488238")
    nodes = query_results.nodes()

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
    min_x = min(x)
    min_y = min(y)
    max_x = max(x)
    max_y = max(y)

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

if __name__ == "__main__":
    #---------------------------------------------------------
    #An example
    example_file = ""

    #This geojson file contains the following *ways* in it:
    ids = get_ids(example_file)
    ids

    #The coordinates of the first *way* in the geojson file above is:
    coords = get_coord(ids[0])
    coords

    #We can also get the coordinates of all the *ways* in the previous file:
    all_coords = list(map(get_coord, ids))
    all_coords

    #Lets get the bound box which envolves the first *way* of the file above:
    bound_box(coords)

    #We can also get the bound boxes to evolve each one of the *ways* of the file above.
    #Note that the boxes are given as Polygon objects
    all_boxes = list(map(bound_box, all_coords))
    all_boxes
