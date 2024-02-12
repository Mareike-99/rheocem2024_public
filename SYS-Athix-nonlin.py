# -*- coding: utf-8 -*-
"""
Created on February 12, 2024

This file is to calculate structural buildup in different ways from a static yield stress test.
In this example, raw data are stored in the same folder as the code. The procedure in the code is: 
    
    1. read the raw data 
    2. average the raw data 
    3. read the shear-ups and maaximum torque after rest, store them in arrays
    4. from the increasing torque values, structural builup is calculated as 
        a) Athix early 
        b) Athix late
        c) nonlinear thixotropy increase 
    5. Data are plotted, parameters are stored in dataframe 

@author: Mareike Thiedeitz
"""

###############################################################################
# Notwendige Imports
import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
import os
import math 
import csv
import numpy as np
import codecs #text encoding
import seaborn as sns
sns.set_theme(color_codes=True) 

from scipy import optimize
import scipy.optimize as scopt
from scipy import interpolate

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


###############################################################################
#global plot properties 
cm_s = 1/2.54  # centimeters in inches
style = 'default'
plt.style.use(style)
plt.rcParams["figure.figsize"] = (10*cm_s,9*cm_s) # Größe der Grafiken
plt.rcParams.update({'font.family':'sans-serif'})
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['legend.fontsize'] = 10 #legend fontsize

color_list = ['indianred', 'cornflowerblue', 'darkblue', 'gray', 'black']

# =============================================================================
# Functions to calculate Athix
# =============================================================================
# linear Athix formulation and corresponding solver to fit experimental data. 
# take care: boundary parameters must be defined manually 
def Athix(t, tau_0, Athix):
     return tau_0 +( Athix * t )
       
        # Solver mit MSE
def solverAthix(Fit_curve, X, Y, p0, MSE):
    param_bounds = ([-100, 0.001], [np.inf, np.inf])
    params, _ = optimize.curve_fit(Fit_curve, X, Y, p0=p0, bounds=param_bounds)  # p0 Startwert
    
    if MSE:
        Y_analytical = Fit_curve(X, *params)
        MSE = np.square(Y - Y_analytical).mean()
        return params, MSE
    else:
        return params
    
# nonlinear solver for thixotropy and rigidification
def Athix_nonlin(t, tau_0, c, theta):
    return tau_0 + c*(1+(lambda_0-1)*np.exp(-t/theta)) + Athix * t

def solver_Athix_nonlin(Fit_curve, X, Y, p0, MSE = True):
    param_bounds=([-500, 0.001, 0.001],[np.inf, np.inf, np.inf]) #c und theta 
    popt, pcov = scopt.curve_fit(Fit_curve, X, Y, p0=p0, bounds=param_bounds) #fit function
    residuals = Y - Fit_curve(X, *popt) #calculate residulas
    ss_res = np.sum(residuals**2) #calculate square of residuals
    Y_analytical = Fit_curve(X,popt[0], popt[1] , popt[2])
    ss_tot = np.sum((Y-np.mean(Y))**2)
    MSE = np.square(Y - Y_analytical).mean()
    r_squared = 1 - (ss_res / ss_tot)    
    return popt, Y_analytical, residuals, ss_res, ss_tot, MSE, r_squared 


# =============================================================================
# Data import
# =============================================================================
#read all data in one folder 
# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# get the path of the excel file 
file_path = os.path.join(current_dir)
data_arr = {} #finally, time, shear rate and shear stress are stored for all files in the folder in this dictionary

#read the raw data into a data dictionary data_arr
for file in os.listdir(file_path): 
    if file.startswith('SYS') and file.endswith('.csv'): #read data 
        INPUT_file = file_path + '/' + file 
        # Read the data into python 
        csvFile = open(INPUT_file,'r')
        csvrdr = csv.reader(csvFile, delimiter=',')
        raw_data = list(csv.reader(codecs.open(INPUT_file, 'rb', 'utf-16'), delimiter='\t', skipinitialspace=True,
                                    quoting=csv.QUOTE_NONE))
        data_arr[file] = []
        df_rawdata = pd.DataFrame(raw_data[7:], columns=raw_data[7]) 
        df_rawdata = df_rawdata.apply(lambda x: x.replace({',':'.'}, regex=True))
        data_arr[file].append(df_rawdata)

# # =============================================================================
# # Calculate Average of raw data & max. torque per static step
# # =============================================================================
#initialize empty arrays and dictionaries to store the raw data and maximum torque data 
mean_dict = {}
drehzahl_mean = []
moment_mean = []

torque_dict = {}
tau_dict = {}
torque_mean = []
tau_mean = []

#TAKE CARE 
#this array is manually defined and declares exactly the seconds where the 
#shear up takes place after rest. This could be automated, but is not for the
#example file. 
shearup = [60.1, 96.1, 162.1, 288.1, 534.1, 1020] 

#read the raw data rotational speed ('Drehzahl'), towrue ('Moment') and time ('Zeit')
#find the shear ups depnding on the defined input second from line 145
#from the torque, calculate the shear stress 
#store torque and shear stress 
for key in data_arr.keys():
    torque_dict[key]= []
    tau_dict[key]= [] 
    drehzahl = np.array(data_arr[key][0]['Drehzahl'][3:], dtype = 'float64')
    moment = np.array(data_arr[key][0]['Moment'][3:], dtype = 'float64')
    time = np.array(data_arr[key][0]['Zeit'][3:], dtype = 'float64')

    drehzahl_mean.append(drehzahl)
    moment_mean.append(moment)
    for i in shearup: 
        start = np.where(time == i)
        start = start[0][0]
        end = start + 59
        max_torque = moment[start:end].max()
        max_tau = max_torque/( 2* math.pi*0.06*0.02**2)/1000
        torque_dict[key].append(max_torque)
        tau_dict[key].append(max_tau)
    torque_mean.append(torque_dict[key])
    tau_mean.append(tau_dict[key])
 
#average the raw data 
y_moment = np.mean(np.array([ i for i in torque_mean]), axis=0 ) #storage
y_tau = np.mean(np.array([ i for i in tau_mean]), axis=0 ) #storage
std_dev_moment = np.std((torque_mean[0],torque_mean[1]), axis=0, ddof=1)
std_dev_tau = np.std((tau_mean[0],tau_mean[1]), axis=0, ddof=1)

## initialize arrays to store data 
paramatrix_dict= {}
Athix_early = []
Athix_late = []
tau0_early = []
tau0_late = []

#Start condition manually defined 
#tau_0 = 0 might not always yield the best result. 
#P0=[tau_0. Athix_e]
p0 = [0, 0.01] 

#calculate Athix_eary and Athix_late from structural increase using linear regression function 
for key in tau_dict.keys():
    paramatrix_dict[key] = []
    y = tau_dict[key][0:4] # nimmt nur die Werte 1,2,3 
    x = np.array(shearup)[0:4] #nimmt nur die Werte 1,2,3 
    params,MSE = solverAthix(Athix,x,y,p0, MSE = True) #:: bedeutet bis Ende
    paramatrix = pd.DataFrame({'Tau_0':params[0],'Athix_early':params[1], 'MSE':MSE}, index = [0])
    Athix_early.append(params[1])
    tau0_early.append(params[0])
    paramatrix_dict[key].append(paramatrix)

for key in tau_dict.keys():
    y = tau_dict[key][3:6]
    x = np.array(shearup)[3:6] 
    params,MSE = solverAthix(Athix,x,y,p0, MSE = True) #:: bedeutet bis Ende
    paramatrix = pd.DataFrame({'Tau_0':params[0],'Athix_late':params[1], 'MSE':MSE}, index = [0])
    Athix_late.append(params[1])
    tau0_late.append(params[0])
    paramatrix_dict[key].append(paramatrix)


#calculate average and standard deviation 
Athix_early_mean = np.mean(Athix_early)
Athix_late_mean = np.mean(Athix_late)

tau0_early_mean = np.mean(tau0_early)
tau0_late_mean = np.mean(tau0_late)

X_analytical_Athix_e = np.arange(shearup[0],shearup[3],0.1) 
X_analytical_Athix_l = np.arange(shearup[3],shearup[5],0.1) 

Y_analytical_Athix_e = Athix(X_analytical_Athix_e,tau0_early_mean, Athix_early_mean)
Y_analytical_Athix_l = Athix(X_analytical_Athix_l,tau0_late_mean, Athix_late_mean)

# =============================================================================
# Nonlinear fit with boundary conditions lambda_0 and late Athix 
# =============================================================================
lambda_0 = 0 #start after pre-shear. Theoretically, the start condition of the structural parameter 
#could vary. In this calculation, it is defined to be 0 after pre-shear. 
#Athix is taken as late strucgtural buildup.  
Athix = Athix_late_mean #bekannt 
x_nonlin = np.linspace(shearup[1], shearup[5], 100)

# start conditions for nonlinear fit. 
#P0 [tau_0, c, theta]
p0_nonlin = [0, 1, 10]
#interpolate from experimental data to apply nonlinear fit to the interpolated array 
func_interpol = interpolate.interp1d(shearup[1:], y_tau[1:])
y_interpol = func_interpol(x_nonlin[:])
fit_params = solver_Athix_nonlin(Athix_nonlin, x_nonlin, y_interpol, p0 = p0_nonlin)
theta = round(fit_params[0][2], 2)
c = round(fit_params[0][1], 2)



# =============================================================================
# Plot averaged raw data
# =============================================================================  

fig, ax1 = plt.subplots(1, 1,dpi=300)
ax1.grid(False)
ax1.set_xticks([0, 250, 500, 750, 1000])
ax1.set_xlabel("t [s]")
ax1.set_ylabel(r"$\tau$" " [Pa]")
Athix_early_mean = round(Athix_early_mean, 2)
Athix_late_mean = round(Athix_late_mean, 2)

count = 0
for key in tau_dict.keys():
    color = color_list[count]
    y = tau_dict[key]
    ax1.plot(shearup, y,'x', color = 'black')
    count +=1


ax1.plot(X_analytical_Athix_e, Y_analytical_Athix_e, linewidth='0.75', 
         color = 'lightcoral', linestyle = '--', label =  r"$\tau (t)_{Athix_e}$")
ax1.plot(X_analytical_Athix_l, Y_analytical_Athix_l, linewidth='0.75', 
         color = 'maroon', linestyle = '-.', label =  r"$\tau (t)_{Athix_l}$")
ax1.fill_between(shearup, y_tau - std_dev_tau, y_tau + std_dev_tau,
                  color='gray', alpha=0.2)

ax1.plot(x_nonlin, fit_params[1], linewidth='1', 
         color = 'darkblue', linestyle = 'dotted', label = r"$\tau (t)_{nonlin}$")

fig.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
        fontsize=8, ncol=1, handleheight=0.5, labelspacing=0.02, 
        loc = 'lower right', handlelength= 1, columnspacing=0.5, 
        bbox_to_anchor=(1,0), bbox_transform=ax1.transAxes)
plt.show()

#store parameters in dataframes 

df_nonlin = pd.DataFrame({'Theta': theta, 'c': c}, index = [0])
df_avedata = pd.DataFrame({"time" : shearup, "Torque" : y_moment, "standarad deviation": std_dev_moment, "Stress": y_tau, "standarad deviation tau": std_dev_tau,})


