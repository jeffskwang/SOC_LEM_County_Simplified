#folders
input_folder = '' #MODIFY THIS: location of input DEMs and SOC folder
results_folder = '' #MODIFY THIS: location folder to output data
#simulation name
county_FIPs = 'COUNTY_FIPS'
#inputs
dem_file = input_folder +'/DEMs/'+county_FIPs +'_buffered.tif' #dem 
soc_a_file = input_folder +'/SOC/'+county_FIPs +'_buffered_soc_000_005_cm.tif' #soc 0 to 5 cm
soc_b_file = input_folder +'/SOC/'+county_FIPs +'_buffered_soc_005_020_cm.tif' #soc 5 to 20 cm
soc_c_file = input_folder +'/SOC/'+county_FIPs +'_buffered_soc_020_050_cm.tif' #soc 20 to 50 cm
soc_d_file = input_folder +'/SOC/'+county_FIPs +'_buffered_soc_050_100_cm.tif' #soc 50 to 100 cm
soc_e_file = input_folder +'/SOC/'+county_FIPs +'_buffered_soc_100_150_cm.tif' #soc 100 to 150 cm
soc_f_file = input_folder +'/SOC/'+county_FIPs +'_buffered_soc_150_999_cm.tif' #soc 150 cm to max depth
soil_depth_file  = input_folder + '/SOC/'+county_FIPs +'_buffered_soil_depth.tif' #soil depth
curv_file  = input_folder + '/curvature/'+county_FIPs +'_buffered_gauss_curv.tif' #curvature

#NaN Value
NaN_value = -3.402823466385288598e+38

###physical parameters###
D = 0.2 #[m^2/yr] diffusion coefficient
La = 0.2 #[m] mixing layer

###numerical parameters###
T = 2000. # [yr] dimensioned simulation time
dt = 1.0# [yr] model timestep

###output parameters###
plots_on = 1 #0 no tiffs made, 1 tiff are made according to dt_plot
dt_plot = 200. # [yr] plot timestep
SOC_levels = [0.0025,0.005,0.0075,0.01,
              0.0125,0.015,0.0175,0.02,
              0.0225,0.025,0.0275,0.03,
              0.0325,0.035,0.0375,0.04,
              0.0425,0.045,0.0475,0.05,
              0.0525,0.055,0.0575,0.06,
              0.0625,0.065,0.0675,0.07,
              0.0725,0.075,0.0775,0.08]

