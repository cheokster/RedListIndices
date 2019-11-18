# .................................................................................................
# Script: GIS-component of calculating country-specific Red List Index for any groups of species
# country-specific-rli.py
# Trialling on: Wedgefishes and giant guitarfishes
# Description: Calculating and extraction for each country the proportional range (i.e. proportion of
# species' total range contained within the country's EEZ) of all species.
# Date: 8 May 2019
# .................................................................................................

## Modules
import os, arcpy, pandas, shutil
from arcpy import env  # this is for setting working environments to run certain arcpy functions
from arcpy.sa import *  # this module is for extracting raster by extent

## Check out necessary licenses
arcpy.CheckOutExtension("spatial")
arcpy.env.overwriteOutput = True

## Set local variables (TOCHANGE)
# Directories
# base directory (this needs to contain the good bathymetry raster and other necessary files, for e.g. EEZ shapefile)
basedir = r'Z:\vmshare\rli'
# geodatabase containing bathymetry related data
bathgdb = r'Z:\vmshare\bathymetry\bathymetry.gdb'
# input directory where original species distribution shapefiles are located
indir = r'Z:\vmshare\sea\shpfiles'
# directory with relevant tables of information, e.g. depth ranges, specist lists, etc.
tabdir = r'Z:\vmshare\sea\Tables'
# output directory where refined (clipped to high-res bathmetry depths) species distribution ranges will be saved to
# create the GDB first (TODO commented out once the GDB is created, to ensure it can't be overwritten!)
#arcpy.CreateFileGDB_management(out_folder_path=basedir, out_name="sea-refine", out_version="CURRENT")
outdir = basedir + '\\sea-refine.gdb'
# directory to hold any temporary files in the middle of analysis loops
# create the GDB first (TODO commented out once the GDB is created, to ensure it can't be overwritten!)
#arcpy.CreateFileGDB_management(out_folder_path=basedir, out_name="tmp-files", out_version="CURRENT")
tmpdir = basedir + '\\tmp-files.gdb'
# directory to hold distribution range proportional files for each species and respective countries
# create the GDB first (TODO commented out once the GDB is created, to ensure it can't be overwritten!)
#arcpy.CreateFileGDB_management(out_folder_path=basedir, out_name="proportional-files", out_version="CURRENT")
propout = basedir + '\\proportional-files.gdb'
# directory for final exports of proportional ranges of each species for each country (TODO: Make sure this dir exists)
propTables = indir[:-8] + 'Tables\\country_props'

# Files
# name of original bathymetry raster
bathyras = bathgdb + '\\mar_bathy'
# name of global EEZs
mareez = bathgdb + '\\mareez_final_prj'
# Tables
# all country codes
# countryCodes = list(pandas.read_csv(basedir + '\\tables\\complete_country_codes.csv'))
# species' countries of occurrence information TODO: Changing filename to match
spCountries = pandas.read_csv(basedir + '\\tables\\se-africa-countries.csv')


# create an object for the desired projected coordinate system
# (all files will need to be in the same PCS prior to any analyses)
# we need the PCS in "well-known text" (WKT) format
# (online database of reference systems in various formats: https://spatialreference.org/)
pcs = 'PROJCS["World_Behrmann",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]],PROJECTION["Behrmann"],PARAMETER["False_Easting",0],PARAMETER["False_Northing",0],PARAMETER["Central_Meridian",0],UNIT["Meter",1]]' # World Behrmann Equal Area Cylindrical PCS
# WKT of geographic coordinate system that unprojected data is in
gcs = 'GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137,298.257223563]],PRIMEM["Greenwich",0],UNIT["Degree",0.017453292519943295]]'


## Prep work for all distribution range shapefiles
# create geodatabase to contain all files (TODO commented out once the GDB is created, to ensure it can't be overwritten!)
#arcpy.CreateFileGDB_management(out_folder_path=basedir, out_name="sea_spp", out_version="CURRENT")
# create local variable for this geodatabase
spGDB = basedir + '\\sea_spp.gdb'
# make a list of all species whose distribution ranges we are working with
# first arrange directories so all shapefiles (and related other associated files) are all contained in the
# main directory (i.e. not split into subdirectories)

# Subfolders = os.listdir(indir)  # get the list of all species subfolders
# for Subfolder in Subfolders:  # iterate through each subfolder
#     sfiles = os.listdir(indir + '\\' + Subfolder)  # get list of file at each subfolder and move into main directory
#     for sfile in sfiles:
#         shutil.move(indir + '\\' + Subfolder + '\\' + sfile, indir)
# (there is now a separate script that does this properly: 'dir_shpfile_organisation.py')

# now create a list of all files in the main directory
#listofFiles = os.listdir(indir)
# # from this create a subset list to contain the shapefiles only
# listofShps = []
#
# for i in listofFiles:  # loop through the all file list and add the ones finishing with shp to the empty list
#     if str(i).endswith(".shp"):
#         listofShps.append(i)
#
# print('List of species maps: ')
# print(listofShps)
# Copy these shapefiles into geodatabase that contains all species ranges shapefiles
arcpy.env.workspace = indir
# List all shapefiles in the arcpy workspace for copying into the first geodatabase
listofShp = arcpy.ListFeatureClasses()
# if ever have unwanted species in this list, a way to get rid of them:
# unwanted = [0, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14]
# for ele in sorted(unwanted, reverse=True):
#     del listofShps[ele]

### In instances where the list of shapefiles we want to work with is SMALLER than the original list of species
### (for e.g. coastal endemics of the region only, as in Sub-equatorial Africa),
### reconcile the listofShps with the list of species names from the depth-ranges csv file,
### which should have the list of species names that we want to work with only
# filepath for depth-range csv file (TODO: need to change filepath to be correct)
depth_csv_path = tabdir + '\\SEA_spp_depths.csv'
# Read in all relevant information from a csv. The 'converters' call is to force read_csv to read it is as a string
# object; range(50) is number equal or greater than number of columns
# https://stackoverflow.com/questions/16988526/pandas-reading-csv-as-string-type
# read in depth-range csv
sp_depth = pandas.read_csv(depth_csv_path,
                           converters={i: str for i in range(50)})  # (TOCHANGE: local filepath)
print(sp_depth)
# Take first column of species name and turn into list
final_sp_names = sp_depth['Name'].tolist()
# Create a list of all shapefiles in the original/master directory
l=os.listdir(indir)
li=[x.split('.')[0] for x in l]  # making sure the list does not include file extensions (shapefiles have many)
# to remove duplicates in list, convert to a set (and, back to a list)
li = list(set(li))
final_namelist = [i for i in li if i in final_sp_names]
# now have to convert final_namelist into a format that will be read as actual shapefiles
final_namelist = map(unicode, final_namelist)  # need to convert to unicode since this is the format arcpy listFeatureClass is
final_namelist = [s + '.shp' for s in final_namelist]  # add .shp extension
len(final_namelist)  # check to see that list length is as it should be

# Execute CopyFeatures for each input shapefile
for shapefile in final_namelist:
    # Determine the new output feature class path and name
    outShp = os.path.join(spGDB, os.path.splitext(shapefile)[0])
    arcpy.CopyFeatures_management(shapefile, outShp)
    arcpy.Delete_management('in_memory')

### To figure out the maximum extent of bathymetry layer that we need to extract, so that we can minimise
### the extent of bathymetry information to handle (by clipping relevant section only)
# First set workspace to the right geodatabase
env.workspace = spGDB
# List all shapefiles in geodatabase to union
# FYI if this function takes too long to execute, consider using the merge function
listSp = arcpy.ListFeatureClasses("*", "POLYGON")
# check there is the right number of species
len(listSp)

arcpy.Union_analysis(listSp, "Allsp_union")
# Deriving extent
fullExt = spGDB + '\\Allsp_union'
desc = arcpy.Describe(fullExt)
# Turn this information into an extent objects
xmin = desc.extent.XMin
xmax = desc.extent.XMax
ymin = desc.extent.YMin
ymax = desc.extent.YMax
# Important !! - set the environmental extent so that the output raster is set to the correct extent of the area
# we are extracting
# arcpy.env.extent = arcpy.Extent(xmin, ymin, xmax, ymax)
arcpy.env.extent = spGDB + '\\Allsp_union'  # the shapefile of all union distribution ranges created earlier
# Set local variables
# from arcpy.sa import *
inRectangle = Extent(XMin=xmin, YMin=ymin, XMax=xmax, YMax=ymax)
# Extract bathymetry based on relevant extent
rasExport = ExtractByRectangle(bathyras, inRectangle, "INSIDE")
# Save the output
rasExport.save(bathgdb + '\\bath_exp')

## Doing the basic work of projecting all relevant files (bathymetry raster, species distribution shapefiles)
# bathymetry raster
arcpy.ProjectRaster_management(in_raster=bathgdb + '\\bath_exp', out_raster=bathgdb + '\\bath_exp_prj',
                               out_coor_system=pcs, resampling_type="BILINEAR", vertical="NO_VERTICAL")

# species distribution files (located in a geodatabase)
# Create new geodatabase that will hold the reprojected distribution maps and output files
# (TODO commented out once the GDB is created, to ensure it can't be overwritten!)
#arcpy.CreateFileGDB_management(out_folder_path=basedir, out_name="sea_prj", out_version="CURRENT")
workGDB = basedir + '\\sea_prj.gdb'  # (TODO)

# make a list of the distribution ranges that need projecting
arcpy.env.workspace = spGDB
listSp = arcpy.ListFeatureClasses("*", "POLYGON")
# Remove the union shapefile of all ranges that is also in the spGDB
listSp.remove('Allsp_union')
# if ever have unwanted species in this list, a way to get rid of them:
# unwanted = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
# for ele in sorted(unwanted, reverse=True):
#     del listSp[ele]
# Print results
print('List of species: ')
print(listSp)
# Begin counter for looping
counter = 0
# Loop through shapefiles in GDB and project into desired projection
for Sp in listSp:
    counter = counter + 1
    arcpy.Project_management(in_dataset=spGDB + '\\' + Sp,
                             out_dataset=workGDB + '\\' + Sp,
                             out_coor_system=pcs,
                             transform_method="",
                             preserve_shape="NO_PRESERVE_SHAPE", max_deviation="", vertical="NO_VERTICAL")
    print('Species ' + Sp + ' projected into geodatabase, ' + str(counter))
    arcpy.Delete_management('in_memory')

##### Clip all species shapefiles to correct bathymetry ranges
# Read in all relevant information from a csv. The 'converters' call is to force read_csv to read it is as a string
# object; range(50) is number equal or greater than number of columns
# https://stackoverflow.com/questions/16988526/pandas-reading-csv-as-string-type
#sp_depth = pandas.read_csv(r'Z:/vmshare/wedge-guitarfish/Tables/spp_depths_fin.csv',
#                           converters={i: str for i in range(50)})  # (TOCHANGE: local filepath)
print(sp_depth)
# if ever need to extract only some rows of the species depths list
# sp_depth = sp_depth.drop([0,1,2,3,4,5,6,7,9,11,12,13,14], axis=0)  # numbers represent row index values
# OR
# sp_depth.loc[ [42,56,60,74] , : ]  # numbers represent row index values
# Extract all species name into a list
sp_names = sp_depth["Name"].tolist()
# sp_names[:] = [name.replace(' ', '_') for name in sp_names]  # to replace unwanted symbols for naming
# List of upper (shallower) depth
upper_depths = sp_depth["Upper_depth"].tolist()
# List of lower (deeper) depth
lower_depths = sp_depth["Lower_depth"].tolist()
# Create tuple list of upper and lower depth values
depth_ranges = list(zip(sp_names, upper_depths, lower_depths))

## Customised 'refine to bathmetry' arcpy function
# https://gis.stackexchange.com/questions/246331/saving-python-script-as-a-tool-in-arcgis


def RefineToBathymetry(workspace, species, mindepth, maxdepth, outspace):
    # Import arcpy module
    import arcpy, os
    arcpy.env.overwriteOutput = True
    # Check out any necessary licenses
    arcpy.CheckOutExtension('spatial')
    # Local variables:
    bathymetry_ras = bathgdb + '\\' + 'bath_exp_prj'  # file path to the desired projected bathymetry raster
    # Process: Get extant of polygon to clip raster
    desc = arcpy.Describe(workspace + '\\' + species)
    frame = str(desc.extent)
    # Process: Clip raster to polygon extent
    arcpy.Clip_management(bathymetry_ras, frame, bathgdb + '/tmp_ras')
    print('... Clipping bathymetry raster to coarse distribution shapefile.')
    # Process: Raster Calculator
    arcpy.gp.ExtractByAttributes_sa(bathgdb + '/tmp_ras', 'Value >= ' + mindepth + ' AND Value <= ' + maxdepth, bathgdb + '/tmp_ras2')
    print('... Clipping bathymetry raster further to accurate depth ranges.')
    # Process: Reclassify raster so all cell values equal 1 (to create uniform mask when convert to polygon)
    arcpy.gp.Reclassify_sa(bathgdb + '/tmp_ras2', 'VALUE', mindepth + ' ' + maxdepth + ' 1', bathgdb + '/tmp_ras3', 'DATA')
    print('... Reclassified bathymetry raster to binary values.')
    # Process: Raster to Polygon
    arcpy.RasterToPolygon_conversion(bathgdb + '/tmp_ras3', 'in_memory/depth_shp', 'NO_SIMPLIFY', 'VALUE')
    print('... Accurate depth raster converted to shapefile.')
    # Process: Clip
    arcpy.Clip_analysis(workspace + '\\' + species, 'in_memory/depth_shp', outspace + '\\' + species + '_clip')
    print('... Coarse distribution shapefile clipped to accurate depth range.')
    # End of loop
    arcpy.Delete_management('in_memory')
    print('End of single loop iteration.')


# Now for the actual work
counter = 0  # counter to keep track of where we are in the loop

for spp, mind, maxd, in depth_ranges:
    counter = counter + 1
    # Print current species working on in this iteration and number tracking
    print('Working on species: ' + spp + ', ' + str(counter) + " out of " + str(len(depth_ranges)))
    # Apply IUCN clipping tool
    RefineToBathymetry(workGDB, spp, mind, maxd, outdir)
    print('... Coarse shapefile clipped to bathymetry.')
    arcpy.Delete_management('in_memory')

# running for loop successfully completed
print('Iterative refinement to species range depths complete.')

## Looping through all species range shapefiles and deriving total geographic range of species (in square kilometres)
# make a list of the distribution ranges that need projecting
arcpy.env.workspace = outdir
listSp_areas = arcpy.ListFeatureClasses("*", "POLYGON")
## if ever have unwanted species in this list, a way to get rid of them:
# unwanted = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
# for ele in sorted(unwanted, reverse=True):
#     del listSp_areas[ele]
counter = 0  # begin counter for looping
# create an empty pandas dataframe to append outputs rows to
totalSpRange = pandas.DataFrame(columns=['species', 'area_sqkm'])
# looping through refined ranges to calculate total area
for sparea in listSp_areas:
    counter = counter + 1
    # add new field to attribute table to calculate area in
    arcpy.AddField_management(sparea, "range_area", "DOUBLE")
    # expression to calculate area in square kilometres
    exp = "!SHAPE.AREA@SQUAREKILOMETERS!"
    arcpy.CalculateField_management(sparea, "range_area", exp, "PYTHON_9.3")
    print('... Total area calculated for ' + sparea)
    # to extract the newly calculated area in square kilometres
    totalSpArea = 0
    with arcpy.da.SearchCursor(sparea, ["range_area"]) as cursor:
        for row in cursor:
            # areas for each row are being summed but this should not be necessary as the output shapefiles should
            # only have one feature/polygon
            totalSpArea = totalSpArea + row[0]
    print('... Total range area for ' + sparea + ' saved.')
    # now append this information to empty dataframe that will contain all proportional areas calculated for
    # each species within each country
    totalSpRange = totalSpRange.append({'species': sparea[:-5], 'area_sqkm': totalSpArea}, ignore_index=True)
    print(totalSpRange.tail(1))  # print the last line that was appended to the dataframe
    # clear any space taken up in arcpy memory
    arcpy.Delete_management('in_memory')
# export all calculated proportional areas to a .csv file
totalSpRange.to_csv(tabdir + '\\outputs\\all_species_total_range.csv')
print('... Loop successfully run until completion.')

## Intersecting with each country's EEZ and calculating proportion of total range for each species, in each country
# note that the ISO-2-digit country codes are used for linking between the SIS exported information and the EEZ data
# extract a list of all unique countries that are included in SIS export of species' countries occurrence
# extracting the 3-digit ISO codes
listCountries = spCountries.iso_3code.unique()
counter = 0
spCounter = 0
# create an empty pandas dataframe to append outputs rows to
countrySp_prop = pandas.DataFrame(columns=['country', 'species', 'area_sqkm'])

# for country in list of countries, do:
for country in listCountries:
    counter = counter + 1
    print('... Iterating through proportional range calculations for country ' + str(counter) + ' out of ' + str(len(listCountries)))
    print(country)
    # export EEZ of country x to temporary directory that holds in-between files
    # first create a table view of the EEZ shapefile (basically an object that represents the attribute table)
    arcpy.MakeFeatureLayer_management(mareez, tmpdir + '\\mareez_table')
    # select country's EEZ by 3-digit ISO code
    arcpy.SelectLayerByAttribute_management(tmpdir + '\\mareez_table', 'NEW_SELECTION', ' "ISO_3digit" = \'' + country + '\' ')  # ridiculous syntax required with backslashes for arcpy
    print('... Country EEZ attribute selected.')
    # save selected attribute to temporary directory
    arcpy.CopyFeatures_management(tmpdir + '\\mareez_table', tmpdir + '\\' + country + '_eez')
    print('... Country EEZ exported to new shapefile in temporary geodatabase.')

    # filter list of all species that occur in country x using boolean expression
    speciesExtract = spCountries[spCountries['iso_3code'] == country]
    print(speciesExtract)
    # then extract all rows with country x in column country_name to derive subset list of species
    speciesList = speciesExtract['taxonid'].tolist()
    # sometimes there will be duplicate rows from the original SIS-exported spreadsheet which need to be removed here
    # otherwise the species loop will be repeating analyses for the same species and country combo
    # do this by turning into a python set (and back to a list)
    speciesList = list(set(speciesList))
    print(speciesList)

    # loop through this list of species names to calculate proportion of range in country x's EEZ
    for spec in speciesList:
        # set species counter
        spCounter = spCounter + 1
        print('... Entering species loop; working on species ' + str(spCounter) + ' out of ' + str(len(speciesList)) + ' for country ' + country)
        # have to replace the ' ' in the species name with an '_' for filename purposes
        spec = spec.replace(' ', '_')
        # little name check necessary: sometimes the species filenames have 'S' appended on the end in which case
        # it needs to happen here too
        # checking against one of the original lists of species filenames
        if listSp[0].endswith('S'):
            spec = spec + 'S'
        print(spec)
        # take species x distribution range and clip to country x EEZ
        arcpy.Clip_analysis(in_features=outdir + '\\' + spec + '_clip', clip_features=tmpdir + '\\' + country + '_eez',
                            out_feature_class=propout + '\\' + country + '_' + spec)
        print('... Distribution of ' + spec + ' clipped to country ' + country)
        # add new field to attribute table to calculate area in
        arcpy.env.workspace = propout
        infileName = propout + '\\' + country + '_' + spec
        arcpy.AddField_management(infileName, "range_area", "DOUBLE")
        # expression to calculate area in square kilometres
        exp = "!SHAPE.AREA@SQUAREKILOMETERS!"
        arcpy.CalculateField_management(infileName, "range_area", exp, "PYTHON_9.3")
        print('... Proportional area calculated in attribute table for ' + spec)
        # to extract the newly calculated area in square kilometres
        spAreainEEZ = 0
        with arcpy.da.SearchCursor(infileName, ["range_area"]) as cursor:
            for row in cursor:
                # areas for each row are being summed but this should not be necessary, as the output shapefiles
                # should only have one feature/polygon
                spAreainEEZ = spAreainEEZ + row[0]
        print('... Proportional area for ' + spec + ' in country ' + country + ', saved.')
        # now append this information to empty dataframe that will contain all proportional areas calculated
        # for each species, within each country
        countrySp_prop = countrySp_prop.append({'country': country, 'species': spec, 'area_sqkm': spAreainEEZ}, ignore_index=True)
        print(countrySp_prop.tail(1))  # print the last line that was appended to the dataframe
        # clear any space taken up in arcpy memory
        arcpy.Delete_management('in_memory')

    # reset species counter
    spCounter = 0
    print('... End of one iteration of listCountries loop.')
# export all calculated proportional areas to a .csv file
countrySp_prop.to_csv(propTables + '\\all_countries_species_range_proportions.csv')
print('... Loop successfully run until completion.')
# End script
print(' (-:  (-:  (-:  All country-specific RLI values calculated for ' + str(len(listCountries)) + ' countries!  :-)  :-)  :-) ')
