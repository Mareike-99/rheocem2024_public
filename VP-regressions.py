# -*- coding: utf-8 -*-
"""
Created on February 12, 2024

This code is to analyze flow curve with viscoplastic phenomenological models. 
There are yield stress related models that analyze the flow curves of the 
shear stress tau. There are also flow curves that analyze the viscosity curves. 
The viscosity curves are the derivative of the shear stress curves. 

For this code to be run, the input data have to be in the same folder as the code. 
If the input data are somewhere else, the path needs to be declared. 
The input data are an excel file with shear rate and shear stress values as 
already averaged data from experiments, The standard deviation is in the excel 
file and also plotted.

Regression functions are applied only for the shear rate range beyond gammad_crit. 

Boundary and start parameters for all regression functions need to be inserted manually. 
This makes the analysis susceptible to errors if the boundary parametes are chosen wrong. 

However, an automated script with testing the best boundary and start parameters
was not implemented, to be aware of the start parameters.   

By running the code, the flow curve is optimized using all phenomenological models. 
The code does decide on the best fit itself. 

@author: Mareike Thiedeitz 
"""

###############################################################################
# required Imports
import matplotlib 
import matplotlib.pyplot as plt
import pandas as pd
import os
import numpy as np
from math import e

from scipy import optimize
import seaborn as sns
sns.set_theme(color_codes=True) 

# For nice plotting 
from IPython.core.display import HTML
from matplotlib.legend_handler import HandlerLine2D


class SymHandler(HandlerLine2D):
    def create_artists(self, legend, orig_handle,xdescent, ydescent, width, height, fontsize, trans):
        xx= 0.6*height
        return super(SymHandler, self).create_artists(legend, orig_handle,xdescent, xx, width, height, fontsize, trans)

from matplotlib import ticker
# Create a ScalarFormatter for the y-axis
formatter = ticker.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)

def make_html(fontname):
    return "<p>{font}: <span style='font-family:{font}; font-size: 24px;'>{font}</p>".format(font=fontname)

code = "\n".join([make_html(font) for font in sorted(set([f.name for f in matplotlib.font_manager.fontManager.ttflist]))])

HTML("<div style='column-count: 2;'>{}</div>".format(code))

###############################################################################
#plot options 
cm_s = 1/2.54  # centimeters in inches
style = 'default'

plt.style.use(style)
plt.rcParams["figure.figsize"] = (10*cm_s,9*cm_s) # Größe der Grafiken

# set the font globally
plt.rcParams.update({'font.family':'sans-serif'})
#.rcParams['text.usetex'] = True #latex font
plt.rcParams['mathtext.fontset'] = 'custom'
plt.rcParams['legend.fontsize'] = 10 #legend fontsize
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'

colorlist = ['indianred', 'gray', 'darkblue', 'black']
marker_list = ["o", "v", "s","*", "+"]

# =============================================================================
# Define functions 
# =============================================================================
#generic solver for all functions 
def solver(Fit_curve, X, Y, p0, param_bounds):
    params, _ = optimize.curve_fit(Fit_curve, X, Y, p0=p0, bounds=param_bounds)
    Y_Data = Fit_curve(X, *params)
    MSE = np.square(Y - Y_Data).mean()
    return params, Y_Data, MSE

# =============================================================================
# Yield stress models
# =============================================================================

def Bingham(gammad, tau_0, mu):
    tau = tau_0+ gammad*mu
    return tau

def HerschelBulkley(gammad, tau_0, k, n):
     tau = tau_0+(k*gammad**n)
     #eta = tau_0 / gammad + k*gammad**(n-1)
     return tau
 

def modBing(gammad, tau_0,mu, mu2):
     tau = tau_0+(mu*gammad)+mu2*gammad**2
     #eta = tau_0 / gammad + mu + mu2 * gammad
     return tau
 
def Casson(gammad, tau_0,eta_inf):
     tau = tau_0+ eta_inf *gammad + 2*(tau_0 * eta_inf)**0.5 * gammad**0.5
     #eta = tau_0 / gammad + eta_inf  +  2*(tau_0 * eta_inf)**0.5 * gammad**0.5
     return tau
 

def Toussaint(gammad, tau_0, beta, k, n):
     tau = tau_0 + (beta*gammad)+ 2*(tau_0*beta)**0.5 * (gammad)**0.5 + (k*gammad**n)
     #eta = tau_0 / gammad + 2*(tau_0)**0.5* (beta*gammad)**-0.5 + beta + (k*gammad**(n-1))
     return tau
 
    
# =============================================================================
# Purely viscosity models
# =============================================================================
 
def Cross(gammad, eta_z, eta_inf, c, d):
    eta = (eta_z-eta_inf)/(1+(c*gammad)**d)+eta_inf
    return eta

def Sisko(gammad, eta_inf, k, n):
    eta = eta_inf + k*gammad**(n-1)
    return eta

def Williamson(gammad, eta_z, k, n):
    eta = eta_z / (1+(k*gammad)**n)
    return eta

def SouzaMendez(gammad, eta_z, tau_0s, tau_0d, k, n, eta_inf, gammad_c):
    eta = (1-(e**-(eta_z*gammad)/tau_0s))*(((tau_0s-tau_0d)/gammad)*e**
                                           -(gammad/gammad_c)+(tau_0d/gammad)
                                           +k*gammad**(n-1))+eta_inf
    return eta

# =============================================================================
# read the data to be analyzed
# =============================================================================
# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# get the path of the excel file 
file_path = os.path.join(current_dir, 'PPdyn_0.45-250.xlsx')

# Read the Excel file
flow_data = pd.read_excel(file_path, sheet_name='Flow data')

# =============================================================================
# Herschel-Bulkley analysis 
# =============================================================================
#Initialize empty parameter matrices
param_list_B = []
param_list_HB = []
param_list_mB = []
param_list_C = []
param_list_T = []

# read the flow data 
rate_arr = flow_data['shear rate']
stress_arr1 = flow_data['shear stress']
std_dev = flow_data['standard deviation']

#find the critical shear rate 
min_stress = np.amin(stress_arr1)
min_stress_ind1 = np.where(stress_arr1==np.amin(stress_arr1))
min_stress_ind1 = int(min_stress_ind1[0]) #calls the integer
    
crit_shear = rate_arr[min_stress_ind1]
#cut the flowcurve at the critical shear rate and read the arrays for x and y 
#beyond the critical sher rate 
x_rate=rate_arr[0:min_stress_ind1]
y_stress=stress_arr1[0:min_stress_ind1]
std_dev = std_dev[0:min_stress_ind1]
    
# =============================================================================
# Regression analysis 
# =============================================================================
#First guess for all models 

#P0=[tau_0, mu1]  # Bingham regression
P0_B = [0.001, 0.001] #Startwerte definieren! Wichtig, daimt die Formel Sinn ergibt     
param_bounds_B=([0.001, 0.001],[np.inf,np.inf])   


#P0=[ tau_0, k, n,]  # Herschel-Bulkley regression
P0_HB = [0.001, 0.001, 0.001] #Startwerte definieren! Wichtig, daimt die Formel Sinn ergibt     
param_bounds_HB=([0.001, 0.001, 0.001],[np.inf,np.inf,np.inf])   

#Input data for mB regression
#P0=[tau_0, mu1a, mu1b]
P0_mB = [0.001, 0.00001, 0.00001] #Startwerte definieren! Wichtig, daimt die Formel Sinn ergibt     
param_bounds_mB=([0.001, 0.00001, -np.inf],[np.inf,np.inf,np.inf])

#Input data for Casson regression
#P0=[tau_0,eta_inf]
P0_C = [0.001, 0.001] #Startwerte definieren! Wichtig, daimt die Formel Sinn ergibt     
param_bounds_C=([0.001, 0.001],[np.inf,np.inf])

#Input data for Toussaint regression
#P0=[tau_0, beta, k, n]
P0_T = [1, 0.01, 0.00001, 1] #Startwerte definieren! Wichtig, daimt die Formel Sinn ergibt     
param_bounds_T=([1, 0.01, 0.00001, 0.8],[np.inf,np.inf,np.inf, np.inf])


# =============================================================================
# Regression analysis for Bingham
# =============================================================================
paramsB, Y_B, MSE_B = solver(Bingham, x_rate,y_stress, P0_B, param_bounds_B) #:: bedeutet bis Ende
paramatrix_B = pd.DataFrame({'Tau_0':paramsB[0],'mu':paramsB[1], 'gamma_crit':crit_shear ,'MSE':MSE_B}, index=['Bingham'])
param_list_B.append(paramatrix_B)
MSE_B = round(MSE_B, 2)

# =============================================================================
# Regression analysis for Herschel Bulkley
# =============================================================================
paramsHB, Y_HB, MSE_HB = solver(HerschelBulkley,x_rate,y_stress, P0_HB, param_bounds_HB) #:: bedeutet bis Ende
paramatrix_HB = pd.DataFrame({'Tau_0':paramsHB[0],'k':paramsHB[1], 'n':paramsHB[2], 'gamma_crit':crit_shear ,'MSE':MSE_HB}, index=['Herschel-Bulkley'])
param_list_HB.append(paramatrix_HB)
MSE_HB = round(MSE_HB, 2)

# =============================================================================
# Regression analysis for modified Bingham
# =============================================================================
paramsmB, Y_mB, MSE_mB = solver(modBing,x_rate,y_stress, P0_mB, param_bounds_mB) #:: bedeutet bis Ende
paramatrix_mB = pd.DataFrame({'Tau_0':paramsmB[0],'mu1':paramsmB[1], 'mu2':paramsmB[2], 'gamma_crit':crit_shear ,'MSE':MSE_mB}, index=['mod. Bingham'])
param_list_mB.append(paramatrix_mB)
MSE_mB = round(MSE_mB, 2)


# =============================================================================
# Regression analysis for Casson
# =============================================================================
paramsC,Y_C, MSE_C = solver(Casson,x_rate,y_stress, P0_C, param_bounds_C) #:: bedeutet bis Ende
paramatrix_C = pd.DataFrame({'Tau_0':paramsC[0],'eta_inf':paramsC[1], 'gamma_crit':crit_shear ,'MSE':MSE_C}, index=['Casson'])
param_list_C.append(paramatrix_C)
MSE_C = round(MSE_C, 2)

# =============================================================================
# Regression analysis for Toussaint
# =============================================================================
paramsT, Y_T, MSE_T = solver(Toussaint,x_rate,y_stress, P0_T, param_bounds_T) #:: bedeutet bis Ende
paramatrix_T = pd.DataFrame({'Tau_0':paramsT[0],'beta':paramsT[1], 'k':paramsT[2], 'n':paramsT[3], 'gamma_crit':crit_shear ,'MSE':MSE_T}, index=['Toussaint'])
param_list_T.append(paramatrix_T)
MSE_T = round(MSE_T, 2)

y_data1 = [Y_B, Y_HB, Y_mB, Y_C, Y_T]

#=============================================================================
# Comparative yield flow curves shear thinning and shear thickening
# =============================================================================
color_list = ['black', 'grey', 'navy', 'maroon', 'indianred']
line_labels = [r'$\tau_{Bingham}$', r'$\tau_{H.-B.}$', 
                r'$\tau_{modB}$', r'$\tau_{C}$', r'$\tau_{T}$',]
linestyle_list = ['solid', '-', '--', ':', '-.']


fig, ax1 = plt.subplots(1, 1,dpi= 300, figsize=(4,4))


ax1.grid(False)
#ax.set_xscale('log')
# ax.set_yscale('log')
ax1.set_ylabel(r'Shear stress $\tau$' " [Pa]")
ax1.set_xlabel('Shear rate $\dot{\gamma}$ [1/s]')
ax1.set_xticks([0, 25, 50, 75, 100])
#ax1.set_xlim(0, phi_m)
#ax1.set_ylim(0, 200)
#ax1.set_yscale('log')

ax1.tick_params(axis='both', which='minor', labelsize=8)
ax1.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.0f'))

ax1.fill_between(x_rate, y_stress - std_dev, y_stress + std_dev,
                 color= 'grey', alpha=0.2)
ax1.plot(x_rate, y_stress, linewidth = '0.1', marker = 'x', 
         color = 'black', label = r'$\tau_{exp}$')
for i in range(0, len(y_data1)):
    color = color_list[i]
    label = line_labels[i]
    ls = linestyle_list[i]
    y = y_data1[i]
    ax1.plot(x_rate, y, linewidth = '1', linestyle = ls,  color = color, label = label)


leg = ax1.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
    fontsize=12, ncol=1, handleheight=1.2, labelspacing=0.02, 
    loc = 'upper left', handlelength= 0.8, columnspacing=1.0)

ax1.text(0.5, -0.2, "(a)", transform=ax1.transAxes,
          fontsize=10, va='center', ha='center')


plt.tight_layout()


# # =============================================================================
# # Storefitted data
# # =============================================================================    
df_regressionB = pd.concat(param_list_B)
df_regressionHB = pd.concat(param_list_HB)
df_regressionmB = pd.concat(param_list_mB)
df_regressionC = pd.concat(param_list_C)
df_regressionT = pd.concat(param_list_T)



    

         



