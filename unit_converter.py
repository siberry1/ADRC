#######################------Copyright 2022, Angus Siberry, All rights reserved------#################################
#This file is part of the ADRC software
#ADRC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ADRC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ADRC. If not, see <https://www.gnu.org/licenses/>.
#######################################################################################################################

import math

def surface_area(R):
    A = 4 * math.pi * (R ** 2)
    return A

def Volume(R):
    V = 4 * math.pi * (R ** 3) / 3
    return V

def shell(Rout,thick):
    Rin = Rout - thick
    V_shell = 4 * math.pi * ((Rout ** 3) - (Rin ** 3)) / 3
    return V_shell

def Flux_to_dose(Flux,dEdx,R,thickness):
    Conversion_factor = 5.76e-7
    Water_layer_mass = shell(R, thickness)
    Mev_g = dEdx / Water_layer_mass #MeV per gram of water
    Dose_rate = Mev_g * Flux * Conversion_factor
    return Dose_rate
