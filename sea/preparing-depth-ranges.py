### Script to prep the csv file containing detailed depth ranges for species
### created: 23 July 2019

# Sub-equatorial Africa subregion
# import necessary modules
import os, csv

# workspace and inputs
seaDir = r'Z:\vmshare\sea'
inDir = r'Z:\vmshare\sea\shpfiles'
# create a new directory where relevant tabled information (in csv format) will go
tableDir = seaDir + '\\Tables'
os.mkdir(tableDir)

# creating the beginning of the csv file of detailed species depth ranges
# list all files of directory containing species shapefiles
listofFiles = os.listdir(inDir)
# from this create a subset list to contain the files with extension .shp only
spList = []  # empty list for appending to

for i in listofFiles:  # loop through the list and append names only to empty list (without extensions)
    if str(i).endswith(".shp"):
        spList.append(i[:-4])

# sort list alphabetically
spList = sorted(spList)

# # If want to pull species names without underscores and proper depth ranges
# # create a list with the folder names and replace '_' with ' ' in each name
# namelist = sorted(os.listdir(inDir))
# if '.DS_Store' in namelist:  # to remove hidden files in the directory
#     namelist.remove('.DS_Store')

## If need to get rid of any symbols for naming purposes
# namelist[:] = [name.replace('_', ' ') for name in namelist]

# set workdir
os.chdir(tableDir)

# save this as .csv file
with open('SEA_spp_depths2.csv', 'wb') as f:  # 'wb' opens file in binary mode, which disables universal newlines
    writer = csv.writer(f)
    writer.writerow(["Name", "Upper_depth", "Lower_depth"])
    for val in spList:
        writer.writerow([val])

