# rheocem2024_public
This repository is to show minimal examples on how rheological data of cementitious pastes, originating from rheometric tests, were analyzed. 


The folder can be cloned and directly run locally - the codes solely refer to data that are
stored in the same folder. The code is open to use and share, licenced under the MIT Licence
(see https://opensource.org/license/mit/.)

————————————-

This repository is to show minimal examples on how rheological data of cementitious pastes,
originating from rheometric tests, were analyzed.

Following analysis methods are available:

————————————- 

"VP-regressions.py"

The file fits experimental flow curves, consisting of the raw data shear rate and shear stress,
to phenomenological models:

• Bingham model

• Herschel-Bulkley model

• modified Bingham model

• Casson model

• Toussaint model

• Cross model

• Sisko model

• Williamson model

Flow curve data are stored in the files "PPdyn_0.45-250.xlsx"

"PPdyn_0.55-250.xlsx"

————————————-

LAOS-intercyclic.py

The code analyzes raw data of a large amplitude oscillatory sweep test with 51 strain amplitudes.
Both first harmonic rheological data and oscillatory wave data are stored. The
code plots the data as strain-sorage modulus / loss modulus plot. The linear viscoelastic
range, critical strain, crossover point and yield index are calculated. The code reads the
data:

"LAOS_0.55_250_1.csv"

"LAOS_0.55_250_2.csv"

"LAOS_0.55_250_3.csv"

————————————-

"Vane_RenerRiwlin.py"

The code analyzes rheometric raw data from large gap Vane-in-Cup rheometry. Several
functions were implemented:

• Reiner-Riwlin equation for the whole gap

• Reiner-Riwlin equation with a partially sheared gap

• Second-Krieger solution (not applied in this code, but implemented)

• function to loop through the raw data with varying input ranges and input range steps

By running this code, the variation of calculated yield stress and viscosity values is plotted,
together with the MSE.

The code reads the data:

"Vane_0.45-250.xlsx"

"Vane_0.55-250.xlsx"

————————————-

"SYS-Athix-nonlin.py"

This file is to calculate structural buildup in different ways from a static yield stress test.
The procedure in the code is:

• Read the raw data

• Average the raw data

• Read the shear-ups and maximum torque after rest, and store them in arrays


• From the increasing torque values, calculate structural buildup as:

1. Athix,early
   
3. Athix,late
   
5. Nonlinear thixotropy increase
   
• Plot the data and store parameters in a dataframe

The code reads the data:

"Vane_0.45-250.xlsx"

"Vane_0.55-250.xlsx"

————————————-

"SYS-Athix-nonlin.py"

This code is to analyze structural buildup. The increase of storage modulus is fitted to a
phenomenological function that analyzes the structural increase in two parts: First, with an
exponential function to calculate a flocculation parameter theta. Second, the linear static
increase after flocculation is analyzed as rigidification parameter.

The procedure in the code is:

• read the raw data

• fit the flocculation function with different input arrays

• fit the rigidification function with different input arrays

• find the best fit

The code reads the data:

"SAOS_0.55_250.csv"
