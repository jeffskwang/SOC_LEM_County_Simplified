===SOC_LEM_County_Simplified===

A model to simulates soil and SOC redistribution at county-scale.

This model is a simplification of SOC_LEM (https://zenodo.org/badge/latestdoi/218085995). The main differences between the models is that this model does not utilize memory to track the change in the 3-dimensional distribution of soil organic carbon (SOC). This simplification allows the tracking of accurate SOC redistribution at large scales. The model requires data that initializes topography and SOC stocks. Datasets can be found (https://scholarworks.umass.edu/data/148/). The model files need to be modified to point the code to the location of the input data. All data should be placed in a main folder named "input". The topographic data needs to be saved in a folder named "DEMs" (digital elevation models). Each topographic file for each county is a single .tif file. Another folder needs to be made for curvature data, named "curvature". The SOC data needs to be saved in a folder named "SOC". For each county, 7 files are required. The files hold data that quantify SOC concentrations at 0-5, 5-20, 20-50, 50-100, 100-150, and 150+ cm, as well as the total soil depth.

===Folder Structure===

+++Folders+++

++drivers
   This folder contains the .py files that drive the code for each county simulation. The "driver_template.py" file is a template for each county driver. The   "make_drivers.py" file will make driver files for each county that has a DEM (digital elevation model) .tif file. 

++utilities
   This folder contains a module that calculates model statistics each time step.

+++FIles+++

+county_list.txt - this file contains the FIPs codes for each county that you want the model to simulate soil and SOC redistribution. Needs moficiation to reflect counties that you want to simulate.

+run_non_parallel.sh and run_parallel.sh - The code can be run on a single processor or in parallel. Example script files have been provided (.sh files). Needs modication to run the counties listed in county_list.txt

+SOC_LEM_County.py - main code, does not require modification to run.

