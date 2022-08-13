# ADRC
Alpha Dose Rate Calculator

----THIS PROJECT HAS MOVED TO THE RAD-DOSE PROJECT (http://radiation-dose.com/) FOR ALPHA,BETA RADIATION THROUGH MULTIPLE MATERIALS----

##################################################################################

Version 2.1

Features on recent update:
- Bug Fixes
- Improvement on point picking regime in 2D providing a result matching the results from Hansson et al 2022
- Fixed radionulcide stacking issue in the Spherical model
- No longer available in macOS due to incombatability with macOS Big Sur v11

Possible development: (As develoments have moved to the Rad-Dose project (Radiation-dose.com) there are no planed developments currently for the ADRC project)
- Making available in macOS pre and post macOS Big Sur v11
- Run-time warning for large simulations
- Material input for geometry
- Surrounding atmosphere options


##################################################################################

About: 
The Alpha Dose Rate Calculator (ADRC) is a mathematical model utilising data from the Projected Range Algorithm (PRAL) from the Stopping Ranges of Ions in Matter (SRIM) software (Ziegler & Biersack 2003), to calculate dose rates from a variety of geometries, activities and energies through water. The current working version is for UO2 and water only but is currently being developed for a multitude of radionuclide contains materials and surrounding environments. For more information regarding the underlying physics please refer to the referenced article below.

Referencing:
A complete technical report is currently being written (a short one is attached). For the time being, If you are using or publishing results from this software please reference the following article

- A. Siberry, D. Hambley, A. Adamska, R. Springell, “A mathematical model to describe the alpha dose rate from a uo2 surface”, Radiation Physics and Chemistry (2020), Volume 182, 2021, 109359, ISSN 0969-806X

https://doi.org/10.1016/j.radphyschem.2021.109359. 

A second article laying out the geometrical cosiderations of a particle geometry can be found here

- A. Siberry, D. Hambley, A. Adamska, R. Springell, “A geometrical model to describe the alpha dose rates from particulates of UO2 in water”, Radiation Physics and Chemistry (2021), Volume 188, 2021, 109677, ISSN 0969-806X

https://doi.org/10.1016/j.radphyschem.2021.109677.

##################################################################################

Troubleshooting:

Note: GUI varies with operating system but back-end is the same for all versions.

Windows:
- Make sure figures are kept in the same directory as the executable otherwise it won’t start up.
- Can use shortcut to pin to Desktop.

Mac OS:
- All required data to run is fully enclosed into one app so no need to worry about directory issues.
- Easily transferred to the desktop without enclosing folder.

GNU/Linux: 
- Make sure figures are kept in the same directory as the executable otherwise it won’t start up.


#########################################################################################

USER GUIDE

Isotopes:
- Ensure you have logged every isotope before running, the model can currently take an array of 100 isotope-activity combinations for a single simulation. 
 - Note each logged input will be run by the simulation number, significantly increasing run time

Crack simulation:
- Ensure you tick the white tick-box next to the slider to tell the simulation to run for the desired crack width.

Decay simulation number:
- Run-times will vary
- For multiple radionuclide inputs reduce the simulation number as it is multiplied by the radionuclide number
- Planar and Crack options have a faster computation time than spherical due to increased complexity of calculation

#########################################################################################

Get in touch about any bugs or developments you would like to see.

Angus Siberry

Email: as14659@bristol.ac.uk
