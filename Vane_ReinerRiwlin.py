# -*- coding: utf-8 -*-
"""
Created on February 12, 2024

This code is to analyze raw data from a decreasing dynamic rotatioanl step-rate 
profile in a large-gap Vane-in-Cup device. 

Raw data are read from the file 

"Vane_0.45-250" inside this folder. 
Several functions were implemented: 
    --> Reiner-Riwlin equation for the whole gap 
    --> Reiner-Riwlin equation with a partially sheared gap
    --> Second-Krieger solution (not applied in this code, but implemented) 
    --> function to loop through the raw data with varying input ranges and input range steps 

By running this code, the variation of calculated yield stress and viscosity values
 is plotted, together with the MSE. The data can be stored in pickle files.
 
 You can read, for comparison
 
 'Vane_0.45-250.xlsx' --rs_list 1
 'Vane_0.55-250.xlsx' -- rs list 2
 
 
@author: Mareike Thiedeitz
"""
# =============================================================================
# # import librabries
# =============================================================================
import numpy as np
import pandas as pd
import math
import matplotlib 
import matplotlib.pyplot as plt
import os 
import pickle
from scipy import optimize
from matplotlib import cm
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
# #functions: Reiner Rivlin Equation for partial and full shear, Second Krieger solution, looping fits 
# =============================================================================

def RR_full(Torque_Nm,tau0,mu):
    R2= np.sqrt(Torque_Nm/(2*math.pi*vane_height*tau0)) #effective R_2
    fitted_speed = Torque_Nm/(4*math.pi*vane_height*mu)* (1/vane_radius**2-1/R2**2)+tau0/mu*np.log(vane_radius/R2)
    return fitted_speed

#rr_partial
def RR_part(Torque_Nm,tau0,mu):
    fitted_speed = tau0 /(2* mu) * (np.log(2*math.pi*vane_height*tau0*vane_radius**2)-1) \
        + Torque_Nm / (4 * math.pi * vane_height * mu * vane_radius**2) \
            - tau0 / (2 * mu) * np.log(Torque_Nm)
    return fitted_speed

#solver to fit tau0 and mu with the RR equation 
def solverRR(Fit_curve, X, Y, p0, param_bounds, MSE):
    params,_ = optimize.curve_fit(Fit_curve, X, Y, p0=p0, bounds=param_bounds) #p0 Startwert
    if MSE == False:
          return params
      
    elif MSE == True:
        Y_analytical = Fit_curve(X,params[0],params[1])
        MSE = np.square(Y - Y_analytical).mean()
        return params,MSE 
    
# =============================================================================
# #second Krieger solution
# =============================================================================
def n_value(tau_b, Rotational_speed):
    ln_tau_b = np.log(tau_b)
    ln_Rotational_speed = np.log(Rotational_speed)
    n = np.gradient(ln_tau_b , ln_Rotational_speed)
    return n 

def gammad_secKrieg(Rotational_speed, r1,r2, n):
    gammad = (2 *Rotational_speed/n ) / (1-r2/r1**(-2/n))
    return gammad

def looping_fits(x,y, start, limit, range_size):

    rheol_params = [] #store here the dictionary of yield stress and viscosity for every step of calculation
    RR_rs = {}
    
    RR_rs['full'] = []
        # Perform the initial fit for the entire arrays x and y
    params, MSE = solverRR(RR_full, y, x, p0=[0.1, 0.1], param_bounds=param_bounds, MSE=True)
    RR_speed = RR_full(Torque_Nm, params[0], params[1])
    
    RR_rs['full'] = RR_speed

    # Store data for the initial fit
    paramatrix_RR = pd.DataFrame({'max_w': Rotational_speed[0],
                                  'min_w': Rotational_speed[-1],
                                  'tau_0': params[0],
                                  'mu': params[1],
                                  'MSE': MSE,
                                  }, index=['full'])
    rheol_params.append(paramatrix_RR)
    
    for i in range(start, limit):
        if i + range_size < limit:
            # define range to fit the data rotational speed and torque 
            X = x[i:i + range_size]
            Y = y[i:i + range_size]
            label = 'range = ' + str(i) + '-' + str(i + range_size) + ', o_max = ' + str(round(Rotational_speed[i],2))
            RR_rs[label] = []
            #fit Reiner-Riwlin equation to evaluate tau_0 and mu 
            params,MSE= solverRR(RR_full, Y, X, p0=[0.1,0.1], 
                                  param_bounds = param_bounds, MSE = True)
            
            RR_speed= RR_full(Torque_Nm,params[0],params[1]) #
            RR_rs[label] = RR_speed
            
            #store data 
            paramatrix_RR = pd.DataFrame({'max_w':Rotational_speed[i], 
                                          'min_w':Rotational_speed[i+range_size],
                                          'tau_0':params[0],'mu':params[1], 
                                          'MSE': MSE}, index = [label])
            
            rheol_params.append(paramatrix_RR) 
        rheol_values = pd.concat(rheol_params)
    return rheol_values, RR_rs

# =============================================================================
# #plot properties 
# =============================================================================
cm_s = 1/2.54  # centimeters in inches
style = 'default'
plt.style.use(style)
plt.rcParams["figure.figsize"] = (10*cm_s,9*cm_s) # Größe der Grafiken
plt.rcParams['legend.fontsize'] = 10 #legend fontsize

#colorlist = ['indianred', 'gray', 'darkblue', 'black']
marker_list = ["o", "v", "s","*", "+", "o", "v", "s","*", "+"]

linestyle_tuple = [
     ('loosely dotted',        (0, (1, 10))),
     ('dotted',                (0, (1, 1))),
     ('densely dotted',        (0, (1, 1))),
     ('long dash with offset', (5, (10, 3))),
     ('loosely dashed',        (0, (5, 10))),
     ('dashed',                (0, (5, 5))),
     ('densely dashed',        (0, (5, 1))),

     ('loosely dashdotted',    (0, (3, 10, 1, 10))),
     ('dashdotted',            (0, (3, 5, 1, 5))),
     ('densely dashdotted',    (0, (3, 1, 1, 1))),

     ('dashdotdotted',         (0, (3, 5, 1, 5, 1, 5))),
     ('loosely dashdotdotted', (0, (3, 10, 1, 10, 1, 10))),
     ('densely dashdotdotted', (0, (3, 1, 1, 1, 1, 1)))]


step = 12
gradient = np.linspace(0.5, 2, step)
colorlist = []
for i in range(0, len(gradient)):
    grad = gradient[i]
    color = cm.PuBuGn(grad)
    colorlist.append(color)
    
# =============================================================================
# Read data 
# =============================================================================
current_dir = os.path.dirname(os.path.abspath(__file__))
# get the path of the excel file 
file_path = os.path.join(current_dir, 'Vane_0.55-250.xlsx')
flow_data = pd.read_excel(file_path, sheet_name='Raw data')

data_arr = []
input = pd.read_excel(file_path , sheet_name='Raw data')
time = input['time'].to_numpy() #change all to array of objects
rs = input['rotational speed'].to_numpy()
torque = input['Torque'].to_numpy()
std = input['standard deviation'].to_numpy()
data_arr.append(flow_data) #read the calculated averaged raw data 

# =============================================================================
# Write the raw data in the corresponding rate-steps to calculate a flow curve from the raw data 
# =============================================================================
#manual list of rotational speed: The raw data are not round. Could be automated in a better way. 
#rs_list = [80, 60,  40, 20, 10,  5, 3, 2, 1, 0.5, 0.25, 0.125, 0.06] 
rs_list = [80, 70, 60,  50, 40, 30,  15, 10, 7.5,  5, 2.5, 1.25, 0.63, 0.32]
rs_arr = np.array(rs_list)
  
tbls = {}
tbls_std = {}
y_torque = []
y_std_full = []

#loop through raw data and find corresponding torque values to the required rotational speed steps       
for i in range(0,len(rs_arr)):
    index = rs_arr[i]
    tbls["rs_{0}".format(i)]=[]
    tbls_std["rs_{0}".format(i)]=[]
    val = rs_arr[i]
    per = val * 0.1
    lower = val - per 
    upper = val + per
    indices = np.where((rs > lower ) & (rs <= upper))
    torque_arr = np.take(torque, indices)
    std_arr = np.take(std, indices)
    tbls["rs_{0}".format(i)].append(torque_arr)
    tbls_std["rs_{0}".format(i)].append(std_arr)
    torque_fin = torque_arr[0][-3] #take last value of stress array
    std_fin = std_arr[0][-3]
    y_torque.append(torque_fin)
    y_std_full.append(std_fin)

# find critical rotational speed to cut the data 
min_torque = np.amin(y_torque)
min_torque_ind = np.where(y_torque==np.amin(y_torque))
min_torque_ind = int(min_torque_ind[0]) #calls the integer
crit_rs = rs_arr[min_torque_ind]

#x and y data for flowcurve
Drehzahl=rs_arr[0:min_torque_ind]
Torque_mNm=np.array(y_torque[0:min_torque_ind], dtype = 'float')
Standard_dev = y_std_full[0:min_torque_ind]

#convert data for units 
Rotational_speed=2*math.pi*Drehzahl/60 # convert rotational speed to [rad/s]
Torque_Nm= Torque_mNm/1000 #convert torque to [Nm]

# Vane cup geometry constants for Reiner Rivlin equation
vane_height=0.06
vane_radius=0.02
#cylinder_radius=0.035 
cylinder_radius=0.05
#calculate stress data at inner and outer radius of the gap 
tau_b = Torque_Nm / (2* math.pi * vane_radius**2 * vane_height) #inner shear stress
tau_c = Torque_Nm / (2* math.pi * cylinder_radius**2 * vane_height) #outer shear stress


# =============================================================================
# Optimum regression analysis
# =============================================================================
limit = len(Drehzahl)
param_bounds=([0.001, 0.001],[np.inf,np.inf])
Statistical_analysis = {}
Statistical_analysis['Torque_data'] = [Torque_Nm, Rotational_speed, np.array(Standard_dev)]
Statistical_analysis['Raw_data'] = [np.array(y_torque), np.array(rs_arr), np.array(y_std_full)]

# the value of start range can be changed 
# the range size varies in the loop function 
for j in range(3,limit):
    Statistical_analysis['range = ' + str(j)] = []
    rheol_values, RR_rs =  looping_fits(Rotational_speed, Torque_Nm,  0,limit, j)
    Statistical_analysis['range = ' + str(j)] = [rheol_values, RR_rs]

# =============================================================================
# Plot all MSE values from Statistical analysis
# =============================================================================
fig, ax1 = plt.subplots(1, 1,dpi=300)
plt.grid(False)  
#plt.title("Experiment data vs Reiner-Rivlin")
ax1.set_xlabel(r"$\omega$ steps for RR-iteration [-]")
ax1.set_ylabel("MSE [-]")

MSE_arr= []
count = 0
for key, df in Statistical_analysis.items():
    if key == 'Torque_data' or key == 'Raw_data':
        continue 
        
    color = colorlist[count]
    val = int(key.split("=")[-1])
    marker = marker_list[count]
    
    MSE = np.array(Statistical_analysis[key][0]['MSE'])
    val = [val] * len(MSE)
    #MSE = np.array(df[0]['MSE'])
    MSE_arr.append(MSE)
    ax1.scatter(val, MSE, marker = marker, color = color, label = key)
    count +=1


leg = ax1.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
    fontsize=8,loc='lower center', bbox_to_anchor=(0.5, -0.4),
    ncol=4, handleheight=1.2, labelspacing=0.02, 
    handlelength= 0.8, columnspacing=1.0)


# =============================================================================
# Plot same but for yield stress
# =============================================================================

fig, ax1 = plt.subplots(1, 1,dpi=300)
plt.grid(False)  
#plt.title("Experiment data vs Reiner-Rivlin")
ax1.set_xlabel(r"$\omega$ steps for RR-iteration [-]")
ax1.set_ylabel(r"Yield stress $\tau_0$ [Pa]")

tau_0_arr= []
count = 0
for key, df in Statistical_analysis.items():
    if key == 'Torque_data' or key == 'Raw_data':
        continue 
    color = colorlist[count]
    val = int(key.split("=")[-1])
    marker = marker_list[count]    
    tau_0 = np.array(Statistical_analysis[key][0]['tau_0'])
    val = [val] * len(tau_0)
    #MSE = np.array(df[0]['MSE'])
    tau_0_arr.append(tau_0)
    ax1.scatter(val, tau_0, marker = marker , color = color, label = key)
    count +=1


leg = ax1.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
    fontsize=8,loc='lower center', bbox_to_anchor=(0.5, -0.4),
    ncol=4, handleheight=1.2, labelspacing=0.02, 
    handlelength= 0.8, columnspacing=1.0)

# =============================================================================
# Save the whole dictionary with all values
# =============================================================================
# OUTPUT_file = filepath + '/' + 'Statistical_analysis.pkl'
# with open(OUTPUT_file, 'wb') as file:
#     pickle.dump(Statistical_analysis, file)

# =============================================================================
#  single plot with selected defined range    
# =============================================================================  

fig, ax1 = plt.subplots(1, 1,dpi=300, figsize=(4,4))
plt.grid(False)  
ax1.set_xlabel(r"$\omega$ [rad/s]")
ax1.set_ylabel("T [Nm]")

rheol_values, RR_rs =  looping_fits(Rotational_speed, Torque_Nm,  0,limit, 4)
ax1.plot(Rotational_speed,Torque_Nm,'o', color = 'black', markersize = 3, label = 'experimental data')

count = 0
for key in RR_rs.keys():
    color = colorlist[count]
    linestyle = linestyle_tuple[count][1]
    o_max = np.round(rheol_values['max_w'][count],1)
    tau_now = np.round(rheol_values['tau_0'][count],1)
    mu = np.round(rheol_values['mu'][count],1)
    label = str(key)
    x_vals = RR_rs[key]
    ax1.plot(x_vals,Torque_Nm,color = color, linewidth = 1, linestyle = linestyle, 
             label = 'RR-reg., '+ r"$\omega_{max}$ = " + str(o_max) + r"; $\tau_{0} = $" + str(tau_now) + r"; $\mu = $" + str(mu)) # add fitted curve
    count +=1

leg = ax1.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
    fontsize=8, loc='lower right', bbox_to_anchor=(1, 0),
    ncol=1, handleheight=1.2, labelspacing=0.02, 
    handlelength= 0.8, columnspacing=1.0)

# # Pickle the DataFrame to a file
# OUTPUT_file = filepath + '/' + 'rheol_vals.pkl'
# with open(OUTPUT_file, 'wb') as file:
#     pickle.dump(rheol_values, file)
    
# =============================================================================
#  single plot with selected defined range    
# =============================================================================  
    
fig, ax1 = plt.subplots(1, 1,dpi=300, figsize=(4,4))
plt.grid(False)  
#plt.title("Experiment data vs Reiner-Rivlin")
ax1.set_xlabel(r"$\omega$ [rad/s]")
ax1.set_ylabel("T [Nm]")
ax1.set_xlim(0.5,2.5)
ax1.set_ylim(0.002,0.006)

rheol_values, RR_rs =  looping_fits(Rotational_speed, Torque_Nm,  0,limit, 4)
ax1.plot(Rotational_speed,Torque_Nm,'s', label = 'exp')

count = 0
for key in RR_rs.keys():
    color = colorlist[count]
    linestyle = linestyle_tuple[count][1]
    o_max = np.round(rheol_values['max_w'][count],1)
    tau_now = np.round(rheol_values['tau_0'][count],1)
    mu = np.round(rheol_values['mu'][count],1)
    label = str(key)
    x_vals = RR_rs[key]
    ax1.plot(x_vals,Torque_Nm,color = color, linewidth = 1, linestyle = linestyle, 
             label = 'RR-reg., '+ r"$\omega_{max}$ = " + str(o_max) + r"; $\tau_{0} = $" + str(tau_now) + r"; $\mu = $" + str(mu)) # add fitted curve
    count +=1

leg = ax1.legend(handler_map={matplotlib.lines.Line2D: SymHandler()}, 
    fontsize=8, loc='lower right', bbox_to_anchor=(1, 0),
    ncol=1, handleheight=1.2, labelspacing=0.02, 
    handlelength= 0.8, columnspacing=1.0)

# # Pickle the DataFrame to a file
# OUTPUT_file = filepath + '/' + 'rheol_vals.pkl'
# with open(OUTPUT_file, 'wb') as file:
#     pickle.dump(rheol_values, file)
    
