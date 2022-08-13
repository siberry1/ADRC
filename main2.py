#######################------Copyright 2022, Angus Siberry, All rights reserved------#################################
#This file is part of the ADRC software
#ADRC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ADRC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ADRC. If not, see <https://www.gnu.org/licenses/>.
#######################################################################################################################

import math
import numpy as np
import random

def range_float(start, stop, step):
    y = start
    while y <= stop:
        yield y
        y = y + step

def sphere_point_picking(thea_limit):
    cycle = range_float(0,100000,1)
    for cyc in cycle:
        u = np.random.uniform(0,1)
        v = np.random.uniform(0,1)
        theta = u * 2 * math.pi
        phi = math.acos((2 * v) - 1)
        r = 1
        sinTheta = math.sin(theta)
        cosTheta = math.cos(theta)
        sinPhi = math.sin(phi)
        cosPhi = math.cos(phi)
        x = r * sinPhi * cosTheta
        z = r * cosPhi

        sign = 1
        if z < 0:
            sign = -sign

        y_ = sign * math.sqrt(1 - (x ** 2))
        theta_ = math.atan(y_ / x)

        if abs(theta_) <= abs(thea_limit):
            break

    return theta_

def uniform_dist(lower_lim, upper_lim):
    low_lim = lower_lim/upper_lim
    low = math.sqrt(5)*(low_lim**3)/3
    upper = math.sqrt(5)/3
    x_in = random.uniform(low,upper)
    x = math.pow(((3*x_in)/math.sqrt(5)),1/3)*upper_lim
    return x

def UO2_range(E):
    f_E = 0.15924410713457773 + (1.4025489663291928 * E) + (0.1692401330655026 * (E ** 2)) + (0.00030296058261558676 * (E ** 3))
    return f_E

def exit_E(E, x):
    f_x = E - ((-0.0748023532999904) + (0.6225037315332709 * x) - (0.0166026423895669 * (x ** 2)) - (0.000095207815644925 * (x ** 3)) + (0.00002310477716497566 * (x ** 4)) - (4.13159291327e-7 * (x ** 5)))
    return f_x

def H2O_range(Ex_E):
    f_range = (0.4910892067117468 + (5.721214198823148 * Ex_E) - (0.9563591666422817 * (Ex_E ** 2)) + (0.43854133180117644 * (Ex_E ** 3)) - (0.046634729935241065 * (Ex_E ** 4)) + (0.0017474358784242412 * (Ex_E ** 5)))
    return f_range

def Bragg(x_stop):
    x_prime = 75 - x_stop
    dedx = 5.4554233561908985e-002 + (1.9933820603675028e-003 * (x_prime)) + (-4.2195888893602207e-004 * (x_prime ** 2)) + (2.6337308106341468e-005 * (x_prime ** 3)) + (3.6869382022456037e-006 * (x_prime ** 4)) + (-6.7103661027537596e-007 * (x_prime ** 5)) + (4.4203929796380202e-008 * (x_prime ** 6)) + (-1.5462542470487985e-009 * (x_prime ** 7)) + (3.1397729731166389e-011 * (x_prime ** 8)) + (-3.2425227038800709e-013 * (x_prime ** 9)) + (-4.4057899127830233e-015 * (x_prime ** 10)) + (3.0899274705689867e-016 * (x_prime ** 11)) + (-5.7420910975847635e-018 * (x_prime ** 12)) + (2.6094497784818177e-020 * (x_prime ** 13)) + (3.3293893160012215e-023 * (x_prime ** 14)) + (3.5427862040635494e-024 * (x_prime ** 15)) + (9.8220260272660832e-026 * (x_prime ** 16)) + (-1.1720811884889775e-027 * (x_prime ** 17)) + (-4.6150898211917356e-029 * (x_prime ** 18)) + (6.0430680146204443e-031 * (x_prime ** 19)) + (-5.8181197611899666e-033 * (x_prime ** 20)) + (5.8345116251414505e-035 * (x_prime ** 21)) + (1.0066059129047870e-036 * (x_prime ** 22)) + (7.7422198035494042e-039 * (x_prime ** 23)) + (-2.8735490315447147e-040 * (x_prime ** 24)) + (-3.1493577849888386e-042 * (x_prime ** 25)) + (5.4859278588090131e-044 * (x_prime ** 26)) + (-5.6580441918456349e-046 * (x_prime ** 27)) + (4.6639869601726191e-048 * (x_prime ** 28)) + (-5.3839486549933306e-050 * (x_prime ** 29)) + (2.7323601619909468e-051 * (x_prime ** 30)) + (-4.8456810349738356e-054 * (x_prime ** 31)) + (-3.7470004217638227e-055 * (x_prime ** 32)) + (-3.1914363863625329e-057 * (x_prime ** 33)) + (8.5377709871860822e-059 * (x_prime ** 34)) + (1.7531365560182840e-061 * (x_prime ** 35)) + (-8.3716828165646151e-063 * (x_prime ** 36)) + (-9.6040653037546997e-066 * (x_prime ** 37)) + (3.7294190095962618e-068 * (x_prime ** 38)) + (1.2689456112273356e-068 * (x_prime ** 39)) + (-8.9920192488068123e-071 * (x_prime ** 40))
    bragg = dedx
    return bragg

def dose(E,Activity,decays,layers):
    density = 10.97
    dR = 0.1
    f_x = np.zeros(layers)
    rf = range_float(0, decays, 1)
    delta = UO2_range(E)  # Gets max range in fuel from energy
    print("Delta=",delta)

    for n in rf:

        x_1 = random.uniform(-delta, 0) #Uses range for depth limit
        y_1 = 0

        theta_lim = math.acos((-x_1) / delta) #Works out maximum angle limits to escape
        theta = sphere_point_picking(theta_lim)

        x_2 = x_1 + (delta * math.cos(theta)) #Creates particle trajectory
        y_2 = y_1 + (delta * math.sin(theta))

        P1 = np.array([x_1, y_1]) #Creates coordinates
        P2 = np.array([x_2, y_2])

        A = np.array([P1, P2]) #Solving for intersection coordinate on the surface
        B = np.array([[0, -delta], [0, delta]])
        s = np.linalg.solve(np.array([A[1] - A[0], B[0] - B[1]]).T, B[0] - A[0])
        inter = ((1 - s) * B[0] + s * B[1])
        d_uo2 = P1 - inter
        dx_uo2 = np.linalg.norm(d_uo2) #Distance travelled in the fuel

        Ex_E = exit_E(E,dx_uo2) #Energy left once exited the fuel

        water_range = H2O_range(Ex_E) #range of particle through water along its trajectory

        x3 = water_range * math.cos(theta) #stopping x coordinate as y is irrelavent in a 1D solution

        layer = range_float(0,((x3/dR)+1),1) #each path is divided into 0.1 increments
        angle = math.sqrt((theta ** 2))
        dx_prime = dR/math.cos(angle) #working out the path distance in each interval

        i = 0

        for i in layer:
            if i < layers:
                position = x3 - i*dR
                LET = f_x[i] + (Bragg(position) * (dx_prime/dR))
                f_x[i] = LET


    flux = density*Activity*(delta*1e-4)*0.25 #0.25 as this is the probanilty of alpha escape for a planar surface
    time = decays / flux #simulated time compared to 1s worth of decays
    dose = f_x / time #scaling the simulation number to decays in a second
    final_result_s = ((dose * 1e6 * 1.62e-19) / ((1e-4 * 1 * 1)/1000)) #scaling to J/kgs as 1kg water = 1000cm^3
    #here we have equated 1 micron in (cm) x 1cm x 1cm for the layer volume * 1 gcm^3 to give it in grams, then / 1000 to Kg
    #final_result_h = final_result_s*3600 #switch from Gy/s to Gy/h

    return final_result_s
