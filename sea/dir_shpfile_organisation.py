### Script to re-organise shapefiles from original IUCN assessment mapping
### Extracting all shapefile files from species subdirectories into just one directory

# Module
import os, shutil, glob
from os.path import join

indir = r'Z:\vmshare\sea\Sub_Equatorial_Africa_Maps'
outdir = r'Z:\vmshare\sea\shpfiles'


# Create list of all files in directory (including all subdirectories, i.e. recursive list)
Subfolders = os.listdir(indir)  # get the list of all species subfolders
for Subfolder in Subfolders:  # iterate through each subfolder
    sfile = glob.glob(indir + '\\' + Subfolder)  # filename for subfolder
    spnames = os.listdir(indir + '\\' + Subfolder)  # shapefiles within species subfolder
    for name in spnames:
        source = sfile[0] + '\\' + str(name)  # the shapefile-file we want to move
        destination = outdir + '\\' + name  # new folder that will contain all the shapfiles
        shutil.move(source, destination)
