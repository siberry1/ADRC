#######################------Copyright 2022, Angus Siberry, All rights reserved------#################################
#This file is part of the ADRC software
#ADRC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ADRC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ADRC. If not, see <https://www.gnu.org/licenses/>.
#######################################################################################################################

import math
import numpy as np

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

def range_float(start, stop, step):
    y = start
    while y <= stop:
        yield y
        y = y + step

def UO2_range(E):
    f_E = 0.15924410713457773 + (1.4025489663291928 * E) + (0.1692401330655026 * (E ** 2)) + (0.00030296058261558676 * (E ** 3))
    return f_E

def prob(R, delta):
    global p_alpha

    delt = delta * 10
    surface_array = np.zeros((int(delt)*10) + 2)
    p_array = np.zeros((int(delt) * 10) + 2)

    if R > delta:

        lowerlim1 = R - delta
        upperlim1 = math.sqrt(abs((delta ** 2) - (R ** 2)))
        lowerlim2 = upperlim1
        upperlim2 = R

        i = 0
        s = 0

        rf1 = range_float(lowerlim1, upperlim1, .1)

        for r in rf1:
            i = i + 1
            x = ((r ** 2) - (delta ** 2) + (R ** 2)) / (2 * r)
            c = abs(r - x)
            h = delta - c
            p = h / (2 * delta)
            surface_array[i] = 4*math.pi*(r**2)
            p_array[i] = p
            sum = s + p
            s = sum
            m = sum / i  # Here m represents the mean probability over the whole rage of r for a given R
            p_alpha = m

        rf2 = range_float(lowerlim2, upperlim2, .1)

        for r in rf2:
            i = i + 1
            x = ((r ** 2) - (delta ** 2) + (R ** 2)) / (2 * r)
            c = abs(r - x)
            h = delta - c
            p = 1 - (h / (2 * delta))
            surface_array[i] = 4*math.pi*(r**2)
            p_array[i] = p
            sum = s + p
            s = sum
            m = sum / i  # Here m represents the mean probability over the whole rage of r for a given R
            p_alpha = m

    if (delta / 2) < R <= delta:

        lowerlim1 = 0
        upperlim1 = abs(R - delta)
        lowerlim2 = upperlim1
        upperlim2 = R

        i = 0
        s = 0

        rf12 = range_float(lowerlim1, upperlim1, .1)

        for r in rf12:
            i = i + 1
            p = 1
            surface_array[i] = 4*math.pi*(r**2)
            p_array[i] = p
            sum = s + p
            s = sum
            m = sum / i  # Here m represents the mean probability over the whole rage of r for a given R
            p_alpha = m

        rf22 = range_float(lowerlim2, upperlim2, .1)

        for r in rf22:
            i = i + 1
            x = -((r ** 2) - (delta ** 2) + (R ** 2)) / (2 * r)
            c = r + x
            h = delta - c
            p = 1 - (h / (2 * delta))
            surface_array[i] = 4*math.pi*(r**2)
            p_array[i] = p
            sum = s + p
            s = sum
            m = sum / i  # Here m represents the mean probability over the whole rage of r for a given R
            p_alpha = m

    if R < (delta / 2):
        m = 1  # Here m represents the mean probability over the whole rage of r for a given R
        p_alpha = m

    surface_array = np.delete(surface_array,0,0)
    p_array = np.delete(p_array, 0, 0)

    if R > (delta / 2):

        Avg_surface = np.mean(surface_array)
        scaling = surface_array / Avg_surface
        scaled_p = p_array * scaling
        p_alpha = np.mean(scaled_p)

    return p_alpha

def palpha(R,E):

    delta = UO2_range(E)

    prob_alpha = prob(R,delta)

    return prob_alpha

def palpha_particle(R,E):
    delta = UO2_range(E)
    smaller_r = R - delta
    V_shell = shell(R,smaller_r)
    V_particle = Volume(R)
    ratio = 1 + V_shell/V_particle
    if R < delta:
        ratio = 1
    prob_shell = prob(R, delta)
    new_palpha = prob_shell/ratio
    return new_palpha
