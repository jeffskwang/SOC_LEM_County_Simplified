#importing libaries
import numpy as np
import os
import importlib
import shutil
import time
import sys
from osgeo import gdal, osr, gdalnumeric
from numba import jit

#read county list
county_file = open('county_list.txt','r')
county_list = county_file.read()
county_list = county_list.split('\n')

#####HOUSEKEEPING#####
#import parameters
code_folder = os.getcwd()
sys.path.append(code_folder +'/drivers')
sys.path.append(code_folder +'/utilities')
parameters = importlib.import_module(county_list[int(sys.argv[1])-1])
globals().update(parameters.__dict__)

#make output folder
if os.path.isdir(results_folder+'/'+county_FIPs):
    shutil.rmtree(results_folder+'/'+county_FIPs)
    time.sleep(3)
os.makedirs(results_folder+'/'+county_FIPs)
os.makedirs(results_folder+'/'+county_FIPs+'/input')
shutil.copyfile(code_folder +'/drivers/'+str(county_list[int(sys.argv[1])-1])+'.py',results_folder+'/'+county_FIPs+'/input/'+str(county_list[int(sys.argv[1])-1])+'.py')

from stats_module import *
#make time series file with header
stats_file = results_folder+'/'+county_FIPs+'/_'+county_FIPs+'_time_series.csv'
with open(stats_file, 'w') as f:
    SOC_less_than_area_header = ''
    for k in range(0,len(SOC_levels)):
        SOC_less_than_area_header += ',area less than ' + str(SOC_levels[k]) + ' g/cm^3 [m^2]'
    f.write('time [yr],total initial soc stock [g],surface soc stock [g],total area [m^2],total area with NaN [m^2],negative curvature surface soc stock [g],positive curvature surface soc stock [g],negative curvature area [m^2],positive curvature area [m^2],eroded soc [g], additional buried soc [g],eroded area [m^2],eroded volume [m^3],deposited area [m^2], deposited volume [m^3]'+SOC_less_than_area_header)

def load_data_tif(dem_fname,
                  curv_fname,
                  soc_a_fname,
                  soc_b_fname,
                  soc_c_fname,
                  soc_d_fname,
                  soc_e_fname,
                  soc_f_fname,
                  soil_depth_fname):
    """
    Load GeoTIFF with gdal and import into landlab
    """
    t = gdal.Open(dem_fname)
    gt = t.GetGeoTransform()
    cs = t.GetProjection()
    cs_sr = osr.SpatialReference()
    cs_sr.ImportFromWkt(cs)

    #NOTE Values are exactly the same as Landlab when using dtype = float64. Our data is in 32-bit float, but landlab automatically converts it to 64-bit.
    #Results are slighlty different when using less precise data, which is necessary for large datasets.
    eta = gdalnumeric.LoadFile(dem_fname).astype(dtype='float32')
    curv = gdalnumeric.LoadFile(curv_fname).astype(dtype='float32')
    soc_a = gdalnumeric.LoadFile(soc_a_fname).astype(dtype='float32')
    soc_b = gdalnumeric.LoadFile(soc_b_fname).astype(dtype='float32')
    soc_c = gdalnumeric.LoadFile(soc_c_fname).astype(dtype='float32')
    soc_d = gdalnumeric.LoadFile(soc_d_fname).astype(dtype='float32')
    soc_e = gdalnumeric.LoadFile(soc_e_fname).astype(dtype='float32')
    soc_f = gdalnumeric.LoadFile(soc_f_fname).astype(dtype='float32')
    soil_depth = gdalnumeric.LoadFile(soil_depth_fname).astype(dtype='float32')
    cellsx =eta.shape[0]
    cellsy =eta.shape[1]

    #convert nodata to zeros
    soc_a[soc_a==NaN_value] = 0.0
    soc_b[soc_b==NaN_value] = 0.0
    soc_c[soc_c==NaN_value] = 0.0
    soc_d[soc_d==NaN_value] = 0.0
    soc_e[soc_e==NaN_value] = 0.0
    soc_f[soc_f==NaN_value] = 0.0
    soil_depth[soil_depth==NaN_value]=200.
    soil_depth[soil_depth<200.]=200.
    
    #convert g/m^2 to g/cm^3
    soc_a *= 1./0.05 / 100. / 100. / 100. #we could divide by bulk density to get gram per gram.
    soc_b *= 1./0.15 / 100. / 100. / 100.
    soc_c *= 1./ 0.30 / 100. / 100. / 100.
    soc_d *= 1./ 0.50 / 100. / 100. / 100.
    soc_e *= 1./0.50 / 100. / 100. / 100.
    soc_f *= 1. / (soil_depth - 150.) / 100. / 100.
    
    return eta,curv,soc_a,soc_b,soc_c,soc_d,soc_e,soc_f ,soil_depth,gt, cs, cellsx,cellsy

#load initial condition
eta,curv,soc_a,soc_b,soc_c,soc_d,soc_e,soc_f,soil_depth,gt, cs, cellsx,cellsy = load_data_tif(dem_file,curv_file,soc_a_file,soc_b_file,soc_c_file,soc_d_file,soc_e_file,soc_f_file,soil_depth_file)

#define arrays
SOC_La = np.zeros_like(eta)
SOC_transfer = np.zeros_like(eta)

#grid size and number of cells
dx = gt[1]
dy = -gt[5]
nt = int(T / dt)
nt_plot = int(dt_plot/dt)

#timer
start_time = time.time()
#########################

#####FUNCTIONS#####
@jit(nopython=True)
def neighbor(x,y,eta,SOC_La):
    #overwrite if cell is at the boundary, where there is no flow, i.e., slope = 0
    N_cell = y + 1
    S_cell = y - 1
    W_cell = x-1
    E_cell = x + 1
    
    if N_cell > cellsy-1:
        N_cell = cellsy - 1
    if S_cell < 0:
        S_cell = 0
    if W_cell < 0:
        W_cell = 0
    if E_cell > cellsx - 1:
        E_cell = cellsx - 1
        
    N_eta_cell = eta[x,N_cell]
    S_eta_cell = eta[x,S_cell]
    W_eta_cell = eta[W_cell,y]
    E_eta_cell = eta[E_cell,y]

    #overwrite if cell is NaN, no flow to NaN, i.e. slope = 0 
    if N_eta_cell == NaN_value:
        N_eta_cell =  eta[x,y]
    if S_eta_cell == NaN_value:
        S_eta_cell =  eta[x,y]
    if W_eta_cell == NaN_value:
        W_eta_cell =  eta[x,y]
    if E_eta_cell == NaN_value:
        E_eta_cell =  eta[x,y]
        
    #initialize as zero
    N_C_link = 0.0 
    S_C_link = 0.0
    W_C_link = 0.0
    E_C_link = 0.0
    #Flux of carbon determined by cell with higher elevation, only flux if slope != 0
    if N_eta_cell > eta[x,y]:
        N_C_link = SOC_La[x,y+1]
    elif N_eta_cell < eta[x,y]:
        N_C_link = SOC_La[x,y]
        
    if S_eta_cell > eta[x,y]:
        S_C_link = SOC_La[x,y-1]
    elif S_eta_cell < eta[x,y]:
        S_C_link = SOC_La[x,y]
        
    if W_eta_cell > eta[x,y]:
        W_C_link = SOC_La[x-1,y]
    elif W_eta_cell < eta[x,y]:
        W_C_link = SOC_La[x,y]
        
    if E_eta_cell > eta[x,y]:
        E_C_link = SOC_La[x+1,y]
    elif E_eta_cell < eta[x,y]:
        E_C_link = SOC_La[x,y]
        
    return N_eta_cell, S_eta_cell, W_eta_cell, E_eta_cell, N_C_link, S_C_link, W_C_link, E_C_link
    
@jit(nopython=True)
def soil_and_SOC_transport(eta,SOC_La):
    dqda = np.zeros_like(eta)
    dqcda = np.zeros_like(eta)
    for x in range (0,cellsx):
        for y in range (0,cellsy):
                if eta[x,y] != NaN_value: #proceed if not NaN
                    #neighbor cells
                    N_eta_cell, S_eta_cell, W_eta_cell, E_eta_cell, N_C_link, S_C_link, W_C_link, E_C_link = neighbor(x,y,eta,SOC_La)

                    #NOTE: S to N is positive y-direction, W to E is positive x-direction
                    N_soil_flux = - D * (N_eta_cell - eta[x,y]) / dy
                    S_soil_flux = - D * (eta[x,y] - S_eta_cell) / dy
                    W_soil_flux = - D * (eta[x,y] - W_eta_cell) / dx
                    E_soil_flux = - D * (E_eta_cell - eta[x,y]) / dx
                    
                    dqdy =  (N_soil_flux - S_soil_flux) / dy
                    dqdx =  (E_soil_flux - W_soil_flux) / dx
                    dqda[x,y] = dqdx + dqdy

                    #non-zero if (1) not on boundary, (2) not flowing to/from NaN, (3) elevations are different
                    N_carbon_flux = N_soil_flux * N_C_link
                    S_carbon_flux = S_soil_flux * S_C_link
                    W_carbon_flux = W_soil_flux * W_C_link
                    E_carbon_flux = E_soil_flux * E_C_link
                    
                    dqcdy = (N_carbon_flux - S_carbon_flux) / dy
                    dqcdx = (E_carbon_flux - W_carbon_flux) / dx
                    dqcda[x,y] = dqcdx + dqcdy
                    
    return dqda,dqcda

@jit(nopython=True)
def SOC_transfer_function_simple(eta_old,eta_ini,soc_a,soc_b,soc_c,soc_d,soc_e,soc_f,soil_depth,dzdt,SOC_La,SOC_transfer):
    interface = -1.0 * (eta_old - eta_ini - La)
    interface_new = -1.0 * (eta_old+ dzdt * dt  - eta_ini - La)
    for i in range(0,cellsx):
        for j in range(0,cellsy):
            if dzdt[i,j] < 0.0: #erosion unearths substrate
                if interface[i,j] <= 0.05:
                    if interface_new[i,j] > 0.05: #if interface crosses boundary, calculate proportion from soc_a and soc_b
                        SOC_transfer[i,j] = soc_a[i,j] * (0.05 - interface[i,j]) / (-dzdt[i,j] * dt)\
                                                                 + soc_b[i,j] * (interface_new[i,j] - 0.05) / (-dzdt[i,j] * dt)
                    else:
                        SOC_transfer[i,j] = soc_a[i,j]
                elif interface[i,j] > 0.05 and interface[i,j] <= 0.2:
                    if interface_new[i,j] > 0.2: #if interface crosses boundary
                        SOC_transfer[i,j] = soc_b[i,j] * (0.2 - interface[i,j]) / (-dzdt[i,j] * dt)\
                                                                 + soc_c[i,j] * (interface_new[i,j] - 0.2) / (-dzdt[i,j] * dt)
                    else:
                        SOC_transfer[i,j] = soc_b[i,j]
                elif interface[i,j] > 0.2 and interface[i,j] <= 0.5:
                    if interface_new[i,j] > 0.5: #if interface crosses boundary
                        SOC_transfer[i,j] = soc_c[i,j] * (0.5 - interface[i,j]) / (-dzdt[i,j] * dt)\
                                                                 + soc_d[i,j] * (interface_new[i,j] - 0.5) / (-dzdt[i,j] * dt)
                    else:
                        SOC_transfer[i,j] = soc_c[i,j]
                elif interface[i,j] > 0.5 and interface[i,j] <= 1.0:
                    if interface_new[i,j] > 1.0: #if interface crosses boundary
                        SOC_transfer[i,j] = soc_d[i,j] * (1.0 - interface[i,j]) / (-dzdt[i,j] * dt)\
                                                                 + soc_e[i,j] * (interface_new[i,j] - 1.0) / (-dzdt[i,j] * dt)
                    else:
                        SOC_transfer[i,j] = soc_d[i,j]
                elif interface[i,j] > 1.0 and interface[i,j] <= 1.5:
                    if interface_new[i,j] > 1.5: #if interface crosses boundary
                        SOC_transfer[i,j] = soc_e[i,j] * (1.5 - interface[i,j]) / (-dzdt[i,j] * dt)\
                                                                 + soc_f[i,j] * (interface_new[i,j] - 1.5) / (-dzdt[i,j] * dt)
                    else:
                        SOC_transfer[i,j] = soc_e[i,j]
                        
                elif interface[i,j] > 1.5 and interface[i,j] <= (soil_depth[i,j] / 100.):
                    if interface_new[i,j] > (soil_depth[i,j] / 100.): #if interface crosses boundary
                        SOC_transfer[i,j] = soc_f[i,j] * ((soil_depth[i,j] / 100.) - interface[i,j]) / (-dzdt[i,j] * dt)
                    else:
                        SOC_transfer[i,j] = soc_f[i,j]                        
                elif interface[i,j] > (soil_depth[i,j] / 100.):
                        SOC_transfer[i,j] = 0.0                  
            elif dzdt[i,j] > 0.0: #deposition deposits the active layer material
                SOC_transfer[i,j] = SOC_La[i,j]
    return SOC_transfer

def write_data_tif(dem_out_fname,eta,gt,cs,cellsx,cellsy):   
    driver = gdal.GetDriverByName('GTiff')
    dataset = driver.Create(dem_out_fname,cellsy,cellsx, 1,gdal.GDT_Float32)
    dataset.GetRasterBand(1).WriteArray(eta)
    dataset.GetRasterBand(1).SetNoDataValue(NaN_value)

    dataset.SetGeoTransform(gt)
    dataset.SetProjection(cs)
    dataset.FlushCache()
    dataset=None

##### LOOP START #####
if La <= 0.05:
    SOC_La = soc_a
elif La > 0.05 and La <= 0.2:
    SOC_La = (0.05 * soc_a + (La - 0.05) * soc_b) / La
elif La > 0.2: #assuming La is always less than 0.5m
    SOC_La = (0.05 * soc_a + 0.2 * soc_b + (La - 0.2) * soc_c) / La

SOC_La[eta==NaN_value]=NaN_value
eta_ini = eta.copy()
for t in range(0,nt + 1):
    stats = calc_stats(t, La,dx,dy,dt,eta,eta_ini,curv,SOC_La,soc_a,soc_b,soc_c,soc_d,soc_e,soc_f,soil_depth,stats_file,NaN_value,np.array(SOC_levels))
    with open(stats_file, 'a') as f:
        #time [yr],total initial soc stock [g],surface soc stock [g],total area [m^2],total area with NaN [m^2],negative curvature surface soc stock [g],positive curvature surface soc stock [g],negative curvature area [m^2],positive curvature area [m^2],eroded soc [g], additional buried soc [g],eroded area [m^2],eroded volume [m^3],deposited area [m^2], deposited volume [m^3]
        SOC_less_than_area_text = ''
        for k in range(0,len(SOC_levels)):
            SOC_less_than_area_text += ',' + str(stats[15][k])
        f.write('\n'+str(stats[0])+','+str(stats[1])+','+str(stats[2])+','+str(stats[3])+','+str(stats[4])+','+str(stats[5])+','+str(stats[6])+','+str(stats[7])+','+str(stats[8])+','+str(stats[9])+','+str(stats[10])+','+str(stats[11])+','+str(stats[12])+','+str(stats[13])+','+str(stats[14])+SOC_less_than_area_text)
        
    if t%nt_plot == 0 and plots_on == 1:
        print ('Time = ' + str(t * dt) + ' years; ' + str(int(float(t)*dt/T*1000.)/10.) + '% done')
        write_data_tif(results_folder+'/'+county_FIPs+'/'+'_eta_'+ '%06d' % t +'yrs.tif',eta,gt,cs,cellsx,cellsy)
        write_data_tif(results_folder+'/'+county_FIPs+'/'+'_soc_'+ '%06d' % t +'yrs.tif',SOC_La,gt,cs,cellsx,cellsy)
    dqda,dqcda = soil_and_SOC_transport(eta,SOC_La)
    SOC_transfer = SOC_transfer_function_simple(eta,eta_ini,soc_a,soc_b,soc_c,soc_d,soc_e,soc_f,soil_depth,-dqda,SOC_La,SOC_transfer) 
    eta += dt *( - dqda)
    SOC_La  += dt/La * (SOC_transfer * dqda  - dqcda)

#end time
stop_time = time.time()
print (str(round((stop_time -start_time )/60.,1))+' mins')
