#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
############################################################################
#
# MODULE:       v.ssn.sites
# AUTHOR(S):    Jason Lessels
# PURPOSE:      Finds the contributing area and statistics for multiple points
#
# COPYRIGHT:
#
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
#############################################################################


#%module
#% description: Determines contributing area and statistics for multiple points along a stream network.
#% keyword: vector
#%end

#%option G_OPT_V_MAP
#%  key: input
#%  description: Name of map with sampling sites.
#%  required: yes
#%  multiple: no
#%end
#%option G_OPT_R_INPUT
#% key: direction
#% label: Name of the direction map
#% description: The flow direction map.
#% required: yes
#%end
#%option G_OPT_R_MAPS
#% key: cont_maps
#% label: Map(s) of continuous covariates
#% description:
#% required: no
#% guisection: Settings
#%end
#%option G_OPT_R_MAPS
#% key: cat_maps
#% label: Map(s) of categorical covariates
#% description:
#% required: no
#% guisection: Settings
#%end
#%option G_OPT_V_FIELD
#% key: layer
#% description: Vector layer number to add catchment statistics
#% answer: 1
#% guisection: Settings
#%end

"""
Created on 12 March 2015

@author: Jason Lessels <jlessels gmail.com>
"""

import sys
import os


import grass.script as grass
from grass.pygrass import vector
from pygrass.modules import Module
from pygrass.modules import stdout2dict
from pygrass import raster
from grass.exceptions import CalledModuleError

if not os.environ.has_key("GISBASE"):
    grass.message( "You must be in GRASS GIS to run this program." )
    sys.exit(1)

def process_inputs(map_strings = ""):
    if(len(map_strings)>1):
        string_list = str.split(map_strings, ",")
        for i in range(len(string_list)):
            string_list[i] = string_list[i].strip()
        return string_list
    else:
        return None

def map_exists(map_name, type = "raster"):
    if(type == "raster"):
        result = grass.find_file(name = map_name, element = 'cell')
        if not result['file']:
            grass.fatal("Raster map <%s> not found" % map_name)
    else:
        result = grass.find_file(name = map_name, element = 'vector')
        if not result['file']:
            grass.fatal("Vector map <%s> not found" % map_name)

def add_column_names_cont(map_name, vect_name, layer):
    master_string = "map_name_min double precision,map_name_max double precision,map_name_mean double precision,map_name_sd double precision"
    col_master_name = grass.find_file(name = map_name, element = 'cell')['name']
    master_string = master_string.replace("map_name",col_master_name)
    Module('v.db.addcolumn', map=vect_name, columns=master_string,layer=layer)

def add_column_names_cat(map_name, vect_name, layer):
    tmp = raster.RasterRow(map_name)
    tmp.has_cats()
    if(tmp.num_cats()>0):
        col_names = []
        col_master_name = grass.find_file(name = map_name, element = 'cell')['name']
        for i in range(len(tmp.cats)):
            col_names.append(col_master_name + "_" + str(tmp.cats[i][1]))
        master_string = " double precision,".join(col_names) + " double precision"
        Module('v.db.addcolumn', map=vect_name, columns=master_string, layer = layer)
    else:
        grass.fatal("Map <%s> does not have layer number" % map_name)


def get_stats_cont(map_to_query, table_name, cat):
    covar_stats = grass.read_command('r.univar', map=map_to_query,quiet=True)
    min = float(covar_stats.split('\n')[6].split(':')[1])
    max = float(covar_stats.split('\n')[7].split(':')[1])
    mean = float(covar_stats.split('\n')[10].split(':')[1])
    sd = float(covar_stats.split('\n')[11].split(':')[1])
    col_master_name = grass.find_file(name = map_to_query, element = 'cell')['name']
    query = "UPDATE " + table_name + " SET " + col_master_name  + "_min" + "=" + str(min) + " WHERE cat = " + str(cat)
    grass.run_command('db.execute', sql=str(query))
    query = "UPDATE " + table_name + " SET " + col_master_name  + "_max" + "=" + str(max) + " WHERE cat = " + str(cat)
    grass.run_command('db.execute', sql=str(query))
    query = "UPDATE " + table_name + " SET " + col_master_name  + "_mean" + "=" + str(mean) + " WHERE cat = " + str(cat)
    grass.run_command('db.execute', sql=str(query))
    query = "UPDATE " + table_name + " SET " + col_master_name  + "_sd" + "=" + str(sd) + " WHERE cat = " + str(cat)
    grass.run_command('db.execute', sql=str(query))




def get_stats_cat(map_name, table_name, cat):
    ### Get the col names - sometimes not all values are found in the sub-catchment.
    tmp = raster.RasterRow(map_name)
    tmp.has_cats()
    col_names = []
    col_master_name = grass.find_file(name = map_name, element = 'cell')['name']
    for i in range(len(tmp.cats)):
        col_names.append(col_master_name + "_" + str(tmp.cats[i][1]))
    ## Get the stats for this layer and save it to a dict.
    stats  = grass.read_command("r.stats", input=map_name, flags='pn', quiet = True)
    stats = stats.split("\n")
    area_dict = {}
    for i in range(len(stats)):
        if(len(stats[i])>0):
            tmp = stats[i].split(" ")
            area_dict[col_master_name + "_" + tmp[0]] = float(tmp[1].replace("%",""))
    ## Save the results - look up the col name in the dict from above. If nothing is found than land cover area must equal 0.
    for i in range(len(col_names)):
        if col_names[i] in area_dict:
            #query = "UPDATE " + prediction_map + "_2 SET " + col_names[i] + "=" + str(area_dict[col_names[i]]) + " WHERE cat = " + str(cat)
            query = "UPDATE " + table_name + " SET " + col_names[i] + "=" + str(area_dict[col_names[i]]) + " WHERE cat = " + str(cat)
            grass.run_command('db.execute', sql=str(query))
        else:
            query = "UPDATE " + table_name + " SET " + col_names[i] + "=" + str(0.0) + " WHERE cat = " + str(cat)
            grass.run_command('db.execute', sql=str(query))


def get_catchment_area(map_name, table_name, area_column, cat):
    area = grass.read_command("r.surf.area",map=map_name, quiet = True, units="kilometers").split("\n")[5].split(":")[1]
    ## And added to the vector attr table.
    query = "UPDATE " + table_name + " SET " + area_column + "=" + area + " WHERE cat = " + str(cat)
    grass.run_command('db.execute', sql=str(query))


def check_for_mask():
    if "MASK" in grass.list_grouped('rast')[grass.gisenv()['MAPSET']]:
        grass.fatal(_('There is already a MASK in place.'
                        'Please remove it before using this module. '))

def get_coords(map_name,layer):
    ## Load the vector with the points in it.
    data = vector.VectorTopo(map_name) # Create a VectorTopo object
    data.open('r', layer = int(layer)) # Open this object for reading
    coords = []
    for i in range(len(data)):
        coords.append(data.read(i+1).coords()) # gives a tuple
    data.close()
    return coords

def get_table_name(map_name, layer):
    ### Test to make sure the vector has the correct layer - and get the table name from it
    try:
        table_name = grass.vector_layer_db(map_name, layer)['name']
    except:
        grass.fatal("Map <%s> does not have layer number %s" % (map_name,layer))
    return table_name


def main():
    ## Cehck to see if a mask has already been applied.
    check_for_mask()
    ## Create a name for the temp file
    tmp_map_name = "tmp_mask_%s" % os.getpid()


    ## Load the inputs from the user
    sites = options['input']
    direction = options['direction']
    layer = options['layer']
    cont_maps = process_inputs(options['cont_maps'])
    cat_maps = process_inputs(options['cat_maps'])
    ## Check all of the chosen maps to make sure they exist.
    map_exists(direction, "raster")
    ## Check to make sure the vector map exists
    map_exists(sites, "vector")
    ## Check to make sure the layer exists and return the table name of that layer
    table_name = get_table_name(sites, layer)
    ## Check if the covariate maps exist and add columns if they do.
    if(isinstance(cont_maps, list)):
        for i in range(len(cont_maps)):
            map_exists(cont_maps[i], "raster")
            add_column_names_cont(cont_maps[i], vect_name = sites, layer = layer)
    if(isinstance(cat_maps, list)):
        for i in range(len(cat_maps)):
            map_exists(cat_maps[i], "raster")
            add_column_names_cat(cat_maps[i], vect_name = sites, layer = layer)

    ## Add the area column if it does not already exist.
    addcolumn = Module('v.db.addcolumn', map=sites, columns="area_km2 double precision", layer =layer)




    coords = get_coords(sites, layer)

    cat = grass.read_command("v.db.select", map=sites, layer=layer, col='cat', flags='c')
    cat = cat.split("\n")



    for i in range(len(coords)):
        ## This first part creates the basin based on the location.
        grass.run_command("r.water.outlet", input = direction, output = tmp_map_name, overwrite=True, quiet = True, coordinates =map(str, coords[i]))
        ## Now the area is calculated.
        grass.run_command("r.mask", rast=tmp_map_name, quiet = True)
        get_catchment_area(tmp_map_name, table_name=table_name, area_column="area_km2", cat=cat[i])
        ## The next section has to loop through all the desired maps to get the stats for each.
        if(isinstance(cont_maps, list)):
            for j in range(len(cont_maps)):
                get_stats_cont(cont_maps[j], table_name= sites, cat=cat[i])
        if(isinstance(cat_maps, list)):
            for j in range(len(cat_maps)):
                get_stats_cat(cat_maps[j], table_name= sites, cat=cat[i])
        #Remove the mask.
        grass.run_command("r.mask", flags="r", quiet = True)

    ## Remove the catchment raster.
    grass.run_command("g.remove", name=tmp_map_name, flags="f", type="raster", quiet = True)




if __name__ == "__main__":
    options, flags = grass.parser()
    main()
