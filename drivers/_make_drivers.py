from tempfile import mkstemp
from shutil import move
from os import fdopen, remove
import os
import shutil
parent = os.getcwd()

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

DEM_folder = '' #MODIFY THIS: folder with DEM in tif format

for elevation_file in os.listdir(DEM_folder):
    if elevation_file.endswith('.tif'):
        fips = elevation_file[:-13]
        shutil.copyfile('_driver_template.py',fips +'.py')
        replace(fips +'.py','COUNTY_FIPS',fips)
