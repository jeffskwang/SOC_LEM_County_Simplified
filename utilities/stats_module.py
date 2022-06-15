from numba import jit
import numpy as np

cubic_m_to_cubic_cm = 100. * 100. * 100.

@jit(nopython=True)
def calc_stats(t,La,dx,dy,dt,eta,eta_ini,curv,SOC_La,soc_a,soc_b,soc_c,soc_d,soc_e,soc_f,soil_depth,stats_file,NaN_value,SOC_levels):
    #time in years
    time = float(t * dt)

    #initialize the sums
    total_initial_stock = 0.0
    total_surface_stock = 0.0

    total_area = 0.0
    total_area_with_NaN = 0.0
    
    negative_surface_stock = 0.0
    positive_surface_stock = 0.0
    negative_curvature_area = 0.0
    positive_curvature_area = 0.0
    
    eroded_area = 0.0
    eroded_volume = 0.0
    deposited_area = 0.0
    deposited_volume =0.0
    
    eroded_carbon_stock = 0.0
    change_in_deposited_surface_stock = 0.0
    SOC_less_than_area = [0.0 for k in range(0,len(SOC_levels))]
    for i in range(0,SOC_La.shape[0]):
        for j in range(0,SOC_La.shape[1]):
            total_area_with_NaN += dx * dy
            if eta[i,j] != NaN_value:               
                #initial stock & current surface stock of the cell, grams
                cell_initial_stock = dx * dy * (0.05 * soc_a[i,j] + 0.15 * soc_b[i,j] + 0.30 * soc_c[i,j] + 0.50 * soc_d[i,j] + 0.50 * soc_e[i,j] + (soil_depth[i,j] - 150.) / 100. * soc_f[i,j]) * cubic_m_to_cubic_cm
                cell_surface_stock = dx * dy * SOC_La[i,j] * La * cubic_m_to_cubic_cm          

                #sum total intial stock
                total_initial_stock +=  cell_initial_stock
                #sum total surface stock
                total_surface_stock += cell_surface_stock
                #sum total area
                total_area += dx * dy

                #negative curvature, surface stock
                if curv[i,j] > 0.0: #ARCMAP calcualtes curvature differently. ARC: + curvature is hilltops, but it should be negative.
                    #see https://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/curvature.htm
                    negative_surface_stock += cell_surface_stock
                    negative_curvature_area += dx * dy
                #positive curvature, surface stock
                elif curv[i,j] < 0.0: #ARCMAP calcualtes curvature differently. ARC: - curvature is hollows, but it should be positive   
                    positive_surface_stock  += cell_surface_stock    
                    positive_curvature_area  += dx * dy    

                #erosional regions
                if eta[i,j] < eta_ini[i,j]:

                    #interface distance from the initial surface
                    interface = La + eta_ini[i,j] - eta[i,j]
                    
                    #calculate the stock left in the substrate in eroded areas
                    if interface <= 0.05:
                        cell_substrate_stock = dx * dy * ((0.05-interface) * soc_a[i,j] + 0.15 * soc_b[i,j] + 0.30 * soc_c[i,j] + 0.50 * soc_d[i,j] + 0.50 * soc_e[i,j] + (soil_depth[i,j] - 150.) / 100. * soc_f[i,j]) * cubic_m_to_cubic_cm
                    elif  interface > 0.05 and interface <= 0.2:
                        cell_substrate_stock = dx * dy * ((0.2 - interface) * soc_b[i,j] + 0.30 * soc_c[i,j] + 0.50 * soc_d[i,j] + 0.50 * soc_e[i,j] + (soil_depth[i,j] - 150.) / 100. * soc_f[i,j]) * cubic_m_to_cubic_cm
                    elif  interface > 0.2 and interface <= 0.5:
                        cell_substrate_stock = dx * dy * ((0.5 - interface) * soc_c[i,j] + 0.50 * soc_d[i,j] + 0.50 * soc_e[i,j] + (soil_depth[i,j] - 150.) / 100. * soc_f[i,j]) * cubic_m_to_cubic_cm
                    elif  interface > 0.5 and interface <= 1.0:
                        cell_substrate_stock = dx * dy * ((1.0 - interface) * soc_d[i,j] + 0.50 * soc_e[i,j] + (soil_depth[i,j] - 150.) / 100. * soc_f[i,j]) * cubic_m_to_cubic_cm
                    elif  interface > 1.0 and interface <= 1.5:
                        cell_substrate_stock = dx * dy * ((1.5 - interface) * soc_e[i,j] + (soil_depth[i,j] - 150.) / 100. * soc_f[i,j]) * cubic_m_to_cubic_cm
                    elif  interface > 1.5 and interface <= (soil_depth[i,j] / 100.):
                        cell_substrate_stock = dx * dy * (((soil_depth[i,j] / 100.)-interface) * soc_f[i,j]) * cubic_m_to_cubic_cm
                    else:
                        cell_substrate_stock = 0.0
                        
                    #soil erosion    
                    eroded_area += dx * dy
                    eroded_volume += dx * dy * (eta_ini[i,j] - eta[i,j])
                    
                    #carbon erosion, the initial stock - what's left (surface + substrate)
                    eroded_carbon_stock += (cell_initial_stock - cell_surface_stock - cell_substrate_stock)

                #depositional regions
                elif eta[i,j] > eta_ini[i,j]:
                    #soil deposition
                    deposited_area += dx * dy
                    deposited_volume += dx * dy * (eta[i,j] - eta_ini[i,j])
                    
                    #initial surface stock, grams
                    if La <= 0.05:
                        cell_initial_surface_stock = dx * dy * La * soc_a[i,j] * cubic_m_to_cubic_cm
                    elif La > 0.05 and La <= 0.2:
                        cell_initial_surface_stock = dx * dy * (0.05 * soc_a[i,j] + (La - 0.05) * soc_b[i,j]) * cubic_m_to_cubic_cm
                    elif La > 0.2: #assuming La is always less than 0.5m
                        cell_initial_surface_stock = dx * dy * (0.05 * soc_a[i,j] + 0.2 * soc_b[i,j] + (La - 0.2) * soc_c[i,j]) * cubic_m_to_cubic_cm

                    #change in the surface stock in the deposited area
                    change_in_deposited_surface_stock += (cell_surface_stock - cell_initial_surface_stock)

                #SOC Less Than
                for k in range(0,len(SOC_levels)):
                    if SOC_La[i,j] <= SOC_levels[k]:
                        SOC_less_than_area[k] += dx * dy
                    

    #how much material is buried from agriculture. This is equal to the mass eroded plus mass from the depositional surface layer if the
    # depositional surface decreases in stock. Control volume is in the depositional area from the interface upward.
    additional_buried_carbon_stock = eroded_carbon_stock - change_in_deposited_surface_stock

    return time, total_initial_stock, total_surface_stock, total_area, total_area_with_NaN, negative_surface_stock, positive_surface_stock,negative_curvature_area,positive_curvature_area,eroded_carbon_stock, additional_buried_carbon_stock, eroded_area, eroded_volume, deposited_area, deposited_volume, SOC_less_than_area

            
    


