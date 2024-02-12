# -*- coding: utf-8 -*-
"""
Created on February 12, 2024

This code is to analye structural buildup. 
From a static small amplitude oscillatory sweep test (SAOS) with an imposed 
strain below the critical strain, the raw data are read. 
The increase of storage modulus is fitted to a phenomenological function that 
analyzes the structural increase in two parts: First, with an exponential function to 
calculate a flocculation parameter theta. Second, the linear static increase after flocculation is 
analyzed as rigidification parameter. 

To find the best fit for both functions, the flocculation fit is performed with 
increasingly more raw data as input data starting from time zero. 
The regression metrics are stored, and the best fit is finally analyzed. 

The same applies for the rigidificagtion fit. However, data are read from the 
back towards the front. 

Having analazed both regression functions, it becomes clear that there is a 
gap between both functions. That indicates that an optimium regression 
function to describe structural buildup has not been found yet. 

As cementitious pastes exhibit both physical flocculation and hydration related 
rigidification, various phenomenological regression functions cannot clearly 
define and distinguish the different effects. The analysis requires more research. 

However, as current approximation, the linear Grigid and the flocculation 
parameter theta for the exponential intial increase suffice for structural 
buildup description. 

@author: Mareike Thiedeitz with Maximilian Prakesch 
"""
# =============================================================================
# Packages to load
# =============================================================================
import numpy as np
import pandas as pd
import scipy.optimize as scopt
import matplotlib.pyplot as plt 
import matplotlib 
import matplotlib.font_manager
import os
import csv
import codecs #text encoding

# For nice plotting 
from IPython.core.display import HTML
from matplotlib.legend_handler import HandlerLine2D

def make_html(fontname):
    return "<p>{font}: <span style='font-family:{font}; font-size: 24px;'>{font}</p>".format(font=fontname)

code = "\n".join([make_html(font) for font in sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))])

HTML("<div style='column-count: 2;'>{}</div>".format(code))

class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle,xdescent, ydescent, width, height, fontsize, trans):
        xx= 0.6*height
        return super(SymHandler, self).create_artists(legend, orig_handle,xdescent, xx, width, height, fontsize, trans)


# =============================================================================
# Define functions for percolation and rigidification rate 
# =============================================================================
def percfit(t, c):
    return G_0 + c*(1-np.exp(-(t/theta)))

def Grigid_fit(t, a, Grigid):
    return a + Grigid*t

guess = [0.001]
def solver_percfit(Fit_curve, X, Y, p0, MSE = True):
    param_bounds=([0.001],[np.inf])
    popt, pcov = scopt.curve_fit(Fit_curve, X, Y, p0=p0, bounds=param_bounds) #fit function
    residuals = Y - Fit_curve(X, *popt) #calculate residulas
    ss_res = np.sum(residuals**2) #calculate square of residuals
    Y_analytical = Fit_curve(X,popt[0])
    ss_tot = np.sum((Y-np.mean(Y))**2)
    MSE = np.square(Y - Y_analytical).mean()
    r_squared = 1 - (ss_res / ss_tot)    
    return popt, Y_analytical, residuals, ss_res, ss_tot, MSE, r_squared 

def solver_Grigid_fit(Fit_curve, X, Y, p0, MSE = True):
    param_bounds = ([-np.inf, 0.001],[np.inf, np.inf])
    popt, pcov = scopt.curve_fit(Fit_curve, X, Y, p0=p0, bounds=param_bounds) #fit function
    residuals = Y - Fit_curve(X, *popt) #calculate residulas
    ss_res = np.sum(residuals**2) #calculate square of residuals
    Y_analytical = Fit_curve(X,popt[0], popt[1])
    ss_tot = np.sum((Y-np.mean(Y))**2)
    MSE = np.square(Y - Y_analytical).mean()
    r_squared = 1 - (ss_res / ss_tot)    
    return popt, Y_analytical, residuals, ss_res, ss_tot, MSE, r_squared 

###############################################################################
#plot properties
cm_s = 1/2.54  # centimeters in inches
style = 'default'
plt.style.use(style)
plt.rcParams["figure.figsize"] = (10*cm_s,9*cm_s) # Größe der Grafiken

# =============================================================================
# Data import and read csv data
# =============================================================================
current_dir = os.path.dirname(os.path.abspath(__file__))
INPUT_file = os.path.join(current_dir, 'SAOS_0.55_250.csv')
csvFile = open(INPUT_file,'r')
csvrdr = csv.reader(csvFile, delimiter=',')
# Read the data into python
raw_data = list(csv.reader(codecs.open(INPUT_file, 'rb', 'utf-16'), delimiter='\t', skipinitialspace=True,
                           quoting=csv.QUOTE_NONE))

#read data from input line 72
df_rawdata = pd.DataFrame(raw_data[72:], columns=raw_data[72])
df_rawdata = df_rawdata.iloc[: , :-6] #cut obsolete columns
df_rawdata = df_rawdata.apply(lambda x: x.replace({',':'.'}, regex=True))

#read time data 
x_time = np.array(df_rawdata['Zeit'][3:], dtype=float) #load time data
#the first time step is at 35.63 seconds. For the fit, the data have to be shofted to 0. 
#For the plot, they are shifted back to the correct measurement time 
x_time = x_time - 35.63 #shift to t0 = 0 for better fitting 
y_G = np.array(df_rawdata['Speichermodul'][3:], dtype=float) #load storage modulus
y_delta = np.array(df_rawdata['Phasenverschiebungswinkel'][3:], dtype=float) #load sthe phase shift angle 
G_0 = y_G[0] #take start condition as first measurfement value 

# =============================================================================
# Flocculation fit 
# =============================================================================
# initial flocculation fit: 
# Try all fit possibilities and store fitting parameters in vectors 

#Initialize matrices for values from curvefit
MSE1 = []     
R_square_1 = []     
theta1 = []         
c1 = []  

#start condition, manually defined 
P0_perc = [0.01]   
for i in range(2,200):
    x_i = x_time[0:i] #von vorn nach hinten durcharbeiten bei x und y 
    y_i = y_G[0:i]
    theta = x_time[i] #temporäres theta bestimmen
    g = solver_percfit(percfit,x_i,y_i,p0=P0_perc,MSE = True) #:: bedeutet bis Ende
    MSE1.append(g[5])   #Vektor für SSE von Perkolation
    R_square_1.append(g[6])  #Vektor für R^2 von Perkolation
    theta1.append(theta)        #Vektor für theta
    c1.append(g[0])            #Vektor für c
    y_anal = percfit(x_time, g[0])

R1 = np.array(R_square_1)
MSE_1 = np.array(MSE1)
#find optimum for flocculation
k_perc = max(R_square_1[9:]) #find best fit for percolation but not in the very first range, this is set manually
k_idx = R_square_1.index(max(R_square_1[9:])) #find index of best fit

theta = x_time[k_idx] + 35.63 #take theta index of best fit 
theta_val = theta 

g_new = solver_percfit(percfit,x_time[0:k_idx],y_G[0:k_idx],p0=P0_perc,MSE = True)
y_perc = percfit(x_time, g_new[0]) #calculate analytical curve for newly fitted c


# =============================================================================
# #Rigidification fit
# =============================================================================
#Initialize matrices for values from curvefit
MSE2 = []          
R_square_2 = []    
Grigid_par = []

P0_rigid = [0.01, 0.01]
for j in range(195, 20,-1):
    x_j = x_time[j:200] 
    y_j = y_G[j:200]
    h = solver_Grigid_fit(Grigid_fit, x_j, y_j, p0 = P0_rigid,MSE = True)
    MSE2.append(g[5])   
    R_square_2.append(h[6])  
    Grigid_par.append(h[0][1])
    h_anal = Grigid_fit(x_time, h[0][0], h [0][1])

k_rig = max(R_square_2[0:])
k_idx_rig = R_square_2.index(max(R_square_2[0:]))
k_rig = 20 + k_idx_rig
h_new = solver_Grigid_fit(Grigid_fit,x_time[(200-k_rig):],y_G[(200-k_rig):],
                          p0=P0_rigid,MSE = True)
y_rigid = Grigid_fit(x_time, h_new[0][0], h_new[0][1])
G_rig = round(h_new[0][1], 2)

# =============================================================================
# Plot raw data and fits for flocculation time and Grigid evolution
# =============================================================================
fig, ax1 = plt.subplots(1, 1,dpi=300, figsize=(4,4))
ax1.grid(False)
ax2 = ax1.twinx()
ax1.set_yscale('log')
#ax1.set_ylim(0,30000)
ax1.set_xlabel('t [s]')
ax1.set_ylabel("G' [Pa] ")
ax2.set_ylabel('$\delta$')

ax1.plot(x_time, y_G, markersize=2, linewidth='0.75',linestyle='--',
          color='black', label = 'Experiment')

ax1.plot(x_time, y_perc, markersize=1, linewidth='1',linestyle='-.',
          color='darkred', label = r"$G'(t), [0 - \Theta]$")

ax1.plot(x_time, y_rigid, markersize=1, linewidth='1', linestyle='--', 
          color='navy', label = r"$G'(t), G_{rig} = $, " + f"{G_rig}" + "Pa/s")
ax1.plot(theta_val - 35.63, y_G[k_idx], 'x', markersize=5, color = 'black', 
         label = r"$\Theta$ = " +  f"{theta_val}" + "s")                 # Additional point

ax2.plot(x_time, y_delta, markersize=1, linewidth='2',linestyle=':',
          color='darkgrey', label = '$\delta$(t)')

fig.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
        fontsize=8, ncol=1, handleheight=0.5, labelspacing=0.02, 
        loc = 'lower right', handlelength= 1, columnspacing=0.5, 
        bbox_to_anchor=(1,0), bbox_transform=ax1.transAxes)
plt.show()

df_fitdata = pd.DataFrame({"theta_delta" : theta, 
                            "theta_perc" : theta_val, 
                            "c": g_new[0], 
                            "Grigid": G_rig, 
                                "a": h_new[0][0]})
   

