# -*- coding: utf-8 -*-
"""
Created on February 12, 2024 

This code is to calculate viscoelastic parameters from an oscillatory experiment 
of the first harmonics. 
Three experimental data are stored in the csv files "LAOS...csv". The data 
are directly read from the rheometer and incorporate stress and strain data of
all oscillatory cycles. So, Lissajous-Bowditch curves can be read and higher harmonics 
can be analyzed by applying Fourier-Transformation or Chebyshev polynomials in 
another code. 

In this code, following steps are conducted: 
    1. raw data averaging 
    2. reading the esperimental data of the first harmonic, directly extracted from 
    experimental raw data (inter-cyclic experimental data)
    3. plotting intercyclic data 
    4. calculation of the viscoelastic parameters LVE, G' at LVE plateau, yield strain, 
    crossover point, yield index

@author: Mareike Thiedeitz 
"""

import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import os
import csv
import codecs #text encoding
import numpy as np
import math

# Import Fonts 
import matplotlib.font_manager
from IPython.core.display import HTML

def make_html(fontname):
    return "<p>{font}: <span style='font-family:{font}; font-size: 24px;'>{font}</p>".format(font=fontname)

code = "\n".join([make_html(font) for font in sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))])

HTML("<div style='column-count: 2;'>{}</div>".format(code))
###############################################################################
#plot properties
cm = 1/2.54  # centimeters in inches
plt.rcParams["figure.figsize"] = (10*cm,9*cm) # Größe der Grafiken
csfont = {'fontname':'Times New Roman'}

# =============================================================================
# Functions to calculate visco-elastic parameters
# =============================================================================
def find_plateau(arr, threshold):
    max_val_idx = np.argmax(arr)
    idx_len = len(storage_mean)
    for i in range(max_val_idx, idx_len, 1):
        if abs( arr[i] / arr[max_val_idx] ) < (1-threshold):
            break 
    G_prime = arr[i+1]
    G_prime_idx = np.where(arr==G_prime)[0][0]
    return G_prime, G_prime_idx 


def tau_Hooke(G1, G2, strain):
    return math.sqrt(G1**2+G2**2)*strain


def gamma_f(x, y):

    # Find where the two arrays cross
    crossover = np.where(np.diff(np.sign(y - x)))[0]   
    # Check if there are multiple crossings
    if len(crossover) > 1:
        diffs = np.diff(crossover )      
        last_pair = np.argmax(diffs)        
        # Get the index of the last crossing point
        last_cross = crossover [last_pair + 1]
    elif len(crossover ) == 1:
        last_cross = crossover [0]
    else:
        last_cross = None
    
    return last_cross

# =============================================================================
# Data import
# =============================================================================
#read all data in one folder 
current_dir = os.path.dirname(os.path.abspath(__file__))
file_path = os.path.join(current_dir)

data_arr = {} #finally, time, shear rate and shear stress are stored for all files in the folder in this dictionary
for file in os.listdir(file_path): 
    if file.startswith('LAOS') and file.endswith('.csv'): 
        INPUT_file = file_path + '/' + file 
        csvFile = open(INPUT_file,'r')
        csvrdr = csv.reader(csvFile, delimiter=',')
        # Read the data into python
        raw_data = list(csv.reader(codecs.open(INPUT_file, 'rb', 'utf-16'), 
                                   delimiter='\t', skipinitialspace=True,
                                    quoting=csv.QUOTE_NONE))

        data_arr[file] = []
        df_rawdata = pd.DataFrame(raw_data[72:], columns=raw_data[72]) 
        df_rawdata = df_rawdata.iloc[: , :-9] #cut obsolete columns
        df_rawdata = df_rawdata.apply(lambda x: x.replace({',':'.'}, regex=True))
        data_arr[file].append(df_rawdata)

# =============================================================================
# Calculate Average of raw data 
# =============================================================================
# initialize empty dictionaries for calculated parameters
raw_dict = {}
vel_dict = {}
  
#read intercyclic data from experimental tests   
for key in data_arr.keys():
    macro_list=[]
    df_rawdata_LAOS = data_arr[key][0]
    raw_dict[key]=[]
    vel_dict[key]=[]
    
    # there are 51 amplitudes with increasing strain.One oscillatory cycle 
    #has 513 data pojnts. 
    #so all first harmonic data are read by starting with 51 and reading through 
    #the rows by i*513 
    #The code is read by manually inserting the to-be-calculated data. If rheometric adjustments are changed, 
    #this needs to be changed, too. 
    for i in range(0,51):
        df_LAOS_i = df_rawdata_LAOS.iloc[3+i*513]
        macro_list.append(df_LAOS_i)
     
    df_LAOS_in = pd.concat(macro_list, axis = 1).T
    df_LAOS = df_LAOS_in.applymap(lambda x: float(x.replace(',', '.')) if isinstance(x, str) and x != '' else x)
    raw_dict[key].append(df_LAOS)
    
    amplitude = np.array(df_LAOS['Scherdeformation'],dtype='float64') #read strains 
    amplitude = amplitude[:,0]
    storage = np.array(df_LAOS['Speichermodul'],dtype='float64') #read storage modulus
    loss = np.array(df_LAOS['Verlustmodul'],dtype='float64') #read loss modulus
    test = np.array([amplitude, storage, loss])
    vel_dict[key].append(test)
    
 
#average the raw data 
df_list = list(vel_dict.values())
raw_array = np.array(df_list)
mean_array = np.mean(raw_array, axis=0)
std_array = np.std(raw_array, axis=0)

strain_mean =mean_array[0][0] 
storage_mean = mean_array[0][1]
loss_mean = mean_array[0][2]

std_storage = std_array[0][1] 
std_loss = std_array[0][2] 

#calculate the phase shift angle 
delta_rad = np.arctan2(np.array(loss_mean),np.array(storage_mean))
# Convert radians to degrees
delta = np.degrees(delta_rad)

# =============================================================================
# LAOS first harmonics intercyclic plot of storage modulus, loss modulus and phase shift angle 
# =============================================================================
fig, ax1 = plt.subplots(1, 1, figsize=(4,4),dpi=300)
ax1.grid(False)
ax2 = ax1.twinx()
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.001, 100)
ax1.tick_params(axis='both', which='minor', labelsize=8)
ax1.set_xlabel(r'$\gamma$' " [%]")
ax1.set_ylabel("G' / G'' [Pa]")
ax2.set_ylabel(r"$\delta$" "[°]")
ax1.plot(strain_mean, storage_mean, color = 'maroon', label = "G' ",marker = 'o', markersize = 3)
ax1.plot(strain_mean, loss_mean,  color = 'grey', label = "G'' ", marker = 'v', markersize = 3)
ax2. plot(strain_mean, delta, color = 'lightgrey', linestyle = '--', label = r"$\delta$")
ax1.fill_between(strain_mean, loss_mean - std_loss, loss_mean + std_loss,
                  color='gray', alpha=0.2)

ax1.fill_between(strain_mean, storage_mean - std_storage, storage_mean + std_storage,
                  color='maroon', alpha=0.2)

fig.legend(loc="lower left", bbox_to_anchor=(0,0), bbox_transform=ax1.transAxes)  
plt.show()

# =============================================================================
# Calculate all relevant visco-elastic parameters
# =============================================================================
plateau = find_plateau(storage_mean, 0.1)
G_storage = plateau[0]
G_loss = loss_mean[plateau[1]]
gamma_yield = strain_mean[plateau[1]]
tau_Hooke = tau_Hooke(G_storage, G_loss, gamma_yield)/100
gamma_f = gamma_f(storage_mean, loss_mean)
yield_index = gamma_f / gamma_yield


# index_labels = ["G storage plateau", "G storage plateau", "G loss plateau", "gamma yield", "gamma fracture", "tau yield", "Yield index"]
data = {"G storage plateau": [G_storage], 
        "G loss plateau": [G_loss], 
        "gamma yield": [gamma_yield], 
        "gamma fracture": [gamma_f], 
        "tau yield": [tau_Hooke], 
        "Yield index": [yield_index]}

vel_par = pd.DataFrame(data)

# =============================================================================
# Store averaged raw data
# =============================================================================    
df_avedata = pd.DataFrame({"strain_mean" : strain_mean, "Storage modulus" : storage_mean, 
                           "Loss modulus": loss_mean, "std storage": 
                               std_storage, "std loss": 
                                   std_loss})
