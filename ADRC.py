#######################------Copyright 2022, Angus Siberry, All rights reserved------#################################
#This file is part of the ADRC software
#ADRC is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#ADRC is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#You should have received a copy of the GNU General Public License along with ADRC. If not, see <https://www.gnu.org/licenses/>.
#######################################################################################################################

from tkinter import *
from tkinter import scrolledtext
from tkinter import messagebox
from tkinter import ttk
from tkinter.filedialog import asksaveasfile
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import math
import numpy as np
import main2
import maincrack
import mainsphere2
import palpha2
import unit_converter

# function to call when user press
# the save button, a filedialog will
# open and ask to save file

def save(func, title):
    files = [('All Files', '*.*')]
    savefile = asksaveasfile(initialfile=str(
        title) + ".csv", filetypes=files, defaultextension=files)
    if savefile:
        np.savetxt(savefile, func, delimiter=",")


def range_float(start, stop, step):
    y = start
    while y <= stop:
        yield y
        y = y + step


def UO2_range(E):
    f_E = 0.15924410713457773 + (1.4025489663291928 * E) + (
        0.1692401330655026 * (E ** 2)) + (0.00030296058261558676 * (E ** 3))
    return f_E


def range_float(start, stop, step):
    y = start
    while y <= stop:
        yield y
        y = y + step


def SuperScriptinate(number):
    return number.replace('0', '⁰').replace('1', '¹').replace('2', '²').replace('3', '³').replace('4', '⁴').replace('5', '⁵').replace('6', '⁶').replace('7', '⁷').replace('8', '⁸').replace('9', '⁹').replace('-', '⁻')


def sci_notation(number, sig_fig=2):
    ret_string = "{0:.{1:d}e}".format(number, sig_fig)
    a, b = ret_string.split("e")
    b = int(b)
    return a + " x10" + SuperScriptinate(str(b))


def spherical_results_window():
    global Ele, HEIGHT, WIDTH, w1, w3, variable4, variable5, variable6, var1, w2, Iso_num, data_collector, activity_collector, linear_results_window, background_image1, text_area, size_image

    if w2.get() == 0:
        messagebox.showerror(title="Invalid input",
                             message="   Error: Please select a Radius > 0 ")
    else:
        if str(Ele.get()) == ' Isotope ':
            messagebox.showerror(
                title="Invalid input", message=" Error: Make sure you 'log' a valid an input before running")
        else:
            if variable4.get() == 'Unit':
                messagebox.showerror(
                    title="Invalid input", message=" Error: Make sure you 'log' a valid an input before running")
            else:
                totE = np.sum(data_collector)
                if totE == 0:
                    messagebox.showerror(title="Invalid input",
                                         message="  Error: Average Energy = 0, please correct this")
                else:
                    totA = np.sum(activity_collector)
                    if totA == 0:
                        messagebox.showerror(title="Invalid input",
                                             message="  Error: Total Activity = 0, please correct this")
                    else:

                        window = Toplevel(master)
                        layers = 1000

                        def savepress():
                            save(dEdxsave, A_dE_dx)
                            save(useless, fixv)

                        Radius = w2.get()
                        Rad = variable5.get()
                        Activ = w3.get()
                        Pref = variable4.get()
                        density = 10.97

                        if str(Pref) == 'GBq':
                            Activity = Activ * 1e9
                            print(Activity)

                        if str(Pref) == 'MBq':
                            Activity = Activ * 1e6
                            print(Activity)

                        if str(Pref) == 'kBq':
                            Activity = Activ * 1e3
                            print(Activity)

                        if str(Pref) == 'Bq':
                            Activity = Activ
                            print(Activity)

                        if str(Rad) == '100\u03BCm':
                            Radius = 100
                            print(Radius)

                        if str(Rad) == '1,000\u03BCm':
                            Radius = 1000
                            print(Radius)

                        if str(Rad) == '10,000\u03BCm':
                            Radius = 10000
                            print(Radius)

                        number = (w1.get() * 1000)
                        print(Radius)
                        print(number)

                        canvas = Canvas(window, height=HEIGHT, width=WIDTH)
                        canvas.pack()

                        results_frame = Frame(window, bg='white', bd=0)
                        results_frame.place(
                            relwidth=2, relheight=2, anchor='n')

                        # bd = boarder bg = background and colour
                        frame = Frame(window, bg='white')
                        frame.place(relx=0.5, rely=0, relwidth=0.5,
                                    relheight=0.05, anchor='n')

                        title = Label(
                            frame, bg='white', text="Full Report", font=(None, 30, 'bold'))
                        title.place(relwidth=1, relheight=1)

                        # bd = boarder bg = background and colour
                        save_frame = Frame(frame, bg="white", bd=0)
                        save_frame.place(relx=0.8, rely=0.2,
                                         width=100, relheight=1, anchor='n')

                        save_button = Button(
                            save_frame, bg="white", text='Save', command=savepress)
                        save_button.pack()

                        lower_frame = Frame(window, bg='white', bd=0)
                        lower_frame.place(
                            relx=0.5, rely=0.1, relwidth=0.8, relheight=0.1, anchor='n')

                        l = Label(lower_frame, bg='white')
                        l.place(relwidth=1, relheight=1)

                        Pref = variable4.get()

                        print(str(Pref))

                        unit = variable6.get()

                        sph = range_float(0, (Iso_num-1), 1)
                        sphere = np.zeros(layers)
                        Dose = np.zeros(layers)
                        palph_ = 0

                        for q in sph:
                            sphere_topup, Dose_topup = mainsphere2.spherical_dose(
                                Radius, number, data_collector[q], activity_collector[q])
                            palph_ = palph_ + palpha2.palpha(Radius,data_collector[q])

                            sphere = sphere + sphere_topup
                            Dose = Dose + Dose_topup

                        lin = range_float(0, (Iso_num-1), 1)
                        linear = np.zeros(layers)

                        for f in lin:
                            linear_topup = main2.dose(
                                data_collector[f], activity_collector[f], number, layers)
                            linear = linear + linear_topup

                        Activity = np.sum(activity_collector)

                        sc = range_float(0, (Iso_num-1), 1)
                        Scaling = np.zeros(Iso_num)
                        for s in sc:
                            Scaling[s] = data_collector[s] * \
                                activity_collector[s]

                        En = np.sum(Scaling) / Activity
                        E_ = round(En, 3)
                        R = Radius
                        dR = 0.1
                        palph = palph_ / (Iso_num)
                        if palph > 1:
                            palph = 1
                        Avg_dE_dx = sphere


                        arr = data_collector
                        sorted_index_array = np.argsort(arr)
                        sorted_array = arr[sorted_index_array]
                        E_max = sorted_array[99]
                        delt = UO2_range(E_max)
                        delta = delt
                        print(E_max)

                        # Using cm to get cm^3 volumes against the gcm^-3 density
                        Particle_Radius = Radius * 1e-4
                        Volume = (4 * math.pi * (Particle_Radius ** 3)) / 3
                        Dose_volume = Volume - \
                            (4 * math.pi * (Particle_Radius - (delt * 1e-4)) ** 3) / 3
                        if R <= delt:
                            Dose_volume = Volume

                        real_flux = density * Activity * Dose_volume * palph
                        E_int = np.sum(Avg_dE_dx)*dR
                        E_per_g = real_flux * E_int
                        f_x_standard = np.zeros(layers)
                        for asd in range_float(0, 999, 1):
                            water_r = (Radius + 0.1*asd) * 1e-4
                            thickness = 1 * 1e-4
                            f_x_standard[asd] = unit_converter.Flux_to_dose(
                                real_flux, Avg_dE_dx[asd], water_r, thickness)
                        Dose_rate = f_x_standard / 3600
                        Tot = np.sum(Dose_rate)
                        Volume_ratio = Dose_volume / Volume
                        Prob_escape = Volume_ratio * palph
                        if R <= (delta/2):
                            Prob_escape = 1

                        particle_mass = density * Volume
                        f_x_particle = f_x_standard

                        tot_particle = np.sum(f_x_particle)
                        Avg = tot_particle / 300

                        if str(unit) == 'Gy/h':
                            Tot = Tot * 3600
                            Dose_rate = Dose_rate * 3600
                            Avg = Avg * 3600
                            linear = linear * 3600
                            tot_particle = tot_particle * 3600

                        fixv = str("FixedVol_R=" + str(Radius) +
                                   "-Activity" + str(Activity) + "Av_E=" + str(E_))

                        A_dE_dx = str("Average_LET_R=" + str(Radius) +
                                      "-Activity" + str(Activity) + "Av_E=" + str(E_))

                        Tot_ = sci_notation(Tot, 3)
                        Avg_ = sci_notation(Avg, 3)
                        tot_particle_ = sci_notation(tot_particle, 3)
                        particle_mass_ = sci_notation(particle_mass, 3)
                        Prob_escape_ = sci_notation(Prob_escape, 3)
                        Activity_ = sci_notation(Activity, 3)
                        E_int_ = sci_notation(E_int, 3)
                        E_per_g_ = sci_notation(E_per_g, 3)


                        results_frame1 = Frame(window, bg='white', bd=0)
                        results_frame1.place(
                            relx=0.25, rely=0.05, relwidth=0.34, relheight=0.4, anchor='n')

                        results_frame3 = Frame(window, bg='white', bd=0)
                        results_frame3.place(
                            relx=0.75, rely=0.05, relwidth=0.35, relheight=0.4, anchor='n')

                        #Figure plotting
                        fig1 = Figure(figsize=(3, 4), dpi=74)
                        t = np.arange(0, layers, 1)
                        x = t * 100 / layers
                        ax = fig1.add_subplot(111)
                        ax.grid(color='grey', linestyle='-', linewidth=0.5)
                        ax.plot(x, Avg_dE_dx, color='red', lw=1,
                                label='R = ' + str(Radius))
                        ax.legend(loc='upper left')
                        ax.set(xlabel='Distance from Interface (x) $\mu$m', ylabel='dE/dx (MeV)',
                               title='Average energy deposition per alpha')
                        dEdxsave = np.column_stack([x, Avg_dE_dx])

                        fig3 = Figure(figsize=(3, 4), dpi=75)
                        t = np.arange(0, layers, 1)
                        x = t * 100 / layers
                        ax = fig3.add_subplot(111)
                        ax.grid(color='grey', linestyle='-', linewidth=0.5)
                        ax.plot(x, Dose, color='blue',
                                lw=1, label='R = ' + str(Radius) + '' + ' $\mu$m')
                        ax.plot(x, linear, label='Planar surface model',
                                color='orange', lw=1, )
                        ax.legend(loc='upper right')
                        ax.set(xlabel='Distance from Interface $\mu$m', ylabel='Dose rate (' + str(unit) + ')',
                               title='Dose rate profile per particulate of UO\u2082')
                        useless = np.column_stack([x, Dose_rate])

                        # A tk.DrawingArea.
                        canvas1 = FigureCanvasTkAgg(
                            fig1, master=results_frame1)
                        canvas1.draw()
                        canvas1.get_tk_widget().place(relx=0.02, rely=0, relwidth=1, relheight=1)

                        # A tk.DrawingArea.
                        canvas3 = FigureCanvasTkAgg(
                            fig3, master=results_frame3)
                        canvas3.draw()
                        canvas3.get_tk_widget().place(relx=0, rely=0, relwidth=1, relheight=1)

                        # bd = boarder bg = background and colour
                        bottom_frame = Frame(window)
                        bottom_frame.place(
                            relx=0.166, rely=0.55, relwidth=0.33, relheight=0.4, anchor='n')

                        # bd = boarder bg = background and colour
                        bottom_frame2 = Frame(window)
                        bottom_frame2.place(
                            relx=0.834, rely=0.55, relwidth=0.33, relheight=0.4, anchor='n')

                        # bd = boarder bg = background and colour
                        bottom_frame3 = Frame(window)
                        bottom_frame3.place(
                            relx=0.395, rely=0.57, relwidth=0.12, relheight=0.4, anchor='n')

                        # bd = boarder bg = background and colour
                        bottom_frame31 = Frame(window)
                        bottom_frame31.place(
                            relx=0.5, rely=0.47, relwidth=0.33, relheight=0.1, anchor='n')

                        # bd = boarder bg = background and colour
                        bottom_frame32 = Frame(window)
                        bottom_frame32.place(
                            relx=0.5, rely=0.57, relwidth=0.11, relheight=0.4, anchor='n')

                        # bd = boarder bg = background and colour
                        bottom_frame33 = Frame(window)
                        bottom_frame33.place(
                            relx=0.61, rely=0.57, relwidth=0.11, relheight=0.4, anchor='n')

                        title_lab = Label(
                            bottom_frame, text="Properties", font=(None, 18, 'bold'))
                        title_lab.place(relx=0.5, rely=0.05, anchor='n')

                        r_lab = Label(bottom_frame, text="Particle radius = " +
                                      str(R) + "\u03BCm", font=(None, 13, 'bold'))
                        r_lab.place(relx=0.5, rely=0.21, anchor='n')

                        f_lab = Label(
                            bottom_frame, text="Fuel type = UO\u2082", font=(None, 13))
                        f_lab.place(relx=0.5, rely=0.32, anchor='n')

                        f_lab = Label(
                            bottom_frame, text="Average decay energy =" + str(E_) + " MeV", font=(None, 13))
                        f_lab.place(relx=0.5, rely=0.39, anchor='n')

                        f_lab = Label(bottom_frame, text="Total activity = " +
                                      str(Activity_) + " Bqg\u207B\u00B9", font=(None, 13))
                        f_lab.place(relx=0.5, rely=0.46, anchor='n')

                        d_lab = Label(bottom_frame, text="Density = " +
                                      str(density) + "gcm\u207B\u00B3", font=(None, 13))
                        d_lab.place(relx=0.5, rely=0.53, anchor='n')

                        d_lab = Label(bottom_frame, text="Particle mass = " +
                                      str(particle_mass_) + "g", font=(None, 13))
                        d_lab.place(relx=0.5, rely=0.60, anchor='n')

                        p = round(palph,2)

                        d_lab = Label(
                            bottom_frame, text="Probability of alpha escape = " + str(Prob_escape_), font=(None, 13))
                        d_lab.place(relx=0.5, rely=0.73, anchor='n')

                        d_lab = Label(bottom_frame, text="Probability of alpha escape =" + str(p) + " \n (within decay depth)   ",
                                      font=(None, 13))
                        d_lab.place(relx=0.5, rely=0.80, anchor='n')

                        title2_lab = Label(
                            bottom_frame2, text="Dosimetry Results", font=(None, 18, 'bold'))
                        title2_lab.place(relx=0.5, rely=0.05, anchor='n')

                        Dose_lab = Label(
                            bottom_frame2, text="Infinite body of water", font=(None, 13, 'bold'))
                        Dose_lab.place(relx=0.5, rely=0.25, anchor='n')

                        Dose_lab = Label(bottom_frame2, text="Total dose rate per particle = " + str(Tot_) + " " + str(unit),
                                         font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0.32, anchor='n')

                        Dose_lab = Label(bottom_frame2,
                                         text="Average escape energy = " +
                                         str(E_int_) +
                                         " " + "MeV",
                                         font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0.39, anchor='n')

                        Dose_lab = Label(bottom_frame2,
                                         text="Energy released per gram = " +
                                              str(E_per_g_) +
                                              " " + "MeV/s",
                                         font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0.46, anchor='n')

                        Dose_lab = Label(
                            bottom_frame2, text="Diffusion layer", font=(None, 13, 'bold'))
                        Dose_lab.place(relx=0.5, rely=0.59, anchor='n')

                        Dose_lab = Label(
                            bottom_frame2, text="(30\u03BCm thick spherical shell)", font=(None, 13,))
                        Dose_lab.place(relx=0.5, rely=0.66, anchor='n')

                        Dose_lab = Label(bottom_frame2, text="Particle dose rate  = " + str(tot_particle_) + " " + str(unit),
                                         font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0.73, anchor='n')

                        AvgDose_lab = Label(bottom_frame2,
                                            text="Average dose rate over diffusion layer = " +
                                            str(Avg_) + " " + str(unit),
                                            font=(None, 13))
                        AvgDose_lab.place(relx=0.5, rely=0.8, anchor='n')

                        # RADIOLYSIS DATA (micro mol per J)

                        H2 = 0.1284
                        H2O2 = 0.104
                        e_aq = 0.01
                        H_rad = 0.0104
                        OH_rad = 0.0364
                        HO2_rad = 0.0104
                        H_plus = 0.01872
                        OH_minus = 0.00312

                        # DIFFUSION LAYER

                        H2_p = 0.1284 * tot_particle
                        H2O2_p = 0.104 * tot_particle
                        e_aq_p = 0.01 * tot_particle
                        H_rad_p = 0.0104 * tot_particle
                        OH_rad_p = 0.0364 * tot_particle
                        HO2_rad_p = 0.0104 * tot_particle
                        H_plus_p = 0.01872 * tot_particle
                        OH_minus_p = 0.00312 * tot_particle

                        H2_p_ = sci_notation(H2_p, 3)
                        H2O2_p_ = sci_notation(H2O2_p, 3)
                        e_aq_p_ = sci_notation(e_aq_p, 3)
                        H_rad_p_ = sci_notation(H_rad_p, 3)
                        OH_rad_p_ = sci_notation(OH_rad_p, 3)
                        HO2_rad_p_ = sci_notation(HO2_rad_p, 3)
                        H_plus_p_ = sci_notation(H_plus_p, 3)
                        OH_minus_p_ = sci_notation(OH_minus_p, 3)

                        # PER GRAM

                        H2_g = 0.1284 * Tot
                        H2O2_g = 0.104 * Tot
                        e_aq_g = 0.01 * Tot
                        H_rad_g = 0.0104 * Tot
                        OH_rad_g = 0.0364 * Tot
                        HO2_rad_g = 0.0104 * Tot
                        H_plus_g = 0.01872 * Tot
                        OH_minus_g = 0.00312 * Tot

                        H2_g_ = sci_notation(H2_g, 3)
                        H2O2_g_ = sci_notation(H2O2_g, 3)
                        e_aq_g_ = sci_notation(e_aq_g, 3)
                        H_rad_g_ = sci_notation(H_rad_g, 3)
                        OH_rad_g_ = sci_notation(OH_rad_g, 3)
                        HO2_rad_g_ = sci_notation(HO2_rad_g, 3)
                        H_plus_g_ = sci_notation(H_plus_g, 3)
                        OH_minus_g_ = sci_notation(OH_minus_g, 3)

                        title2_lab = Label(
                            bottom_frame31, text="Radiolysis Rates", font=(None, 18, 'bold'))
                        title2_lab.place(relx=0.5, rely=0.1, anchor='n')

                        title22_lab = Label(bottom_frame31,
                                            text="G-values taken from Pastina et al (2001) \n at 25\u00B0C (\u03BCmolL\u207B\u00B9)",
                                            font=(None, 13,))
                        title22_lab.place(relx=0.5, rely=0.5, anchor='n')

                        Dose_lab = Label(bottom_frame3,
                                         text="H\u2082 \n\n H\u2082O\u2082 \n\n e(aq) \n\n H(rad) \n\n OH(rad) \n\n HO\u2082(rad) \n\n H+ \n\n OH- ",
                                         font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0.1, anchor='n')

                        Dose_lab = Label(bottom_frame32,
                                         text="Diffusion Layer \n\n" + str(H2_p_) + "\n\n" + str(H2O2_p_) + "\n\n" + str(
                                             e_aq_p_) + "\n\n" + str(H_rad_p_) + '\n\n' + str(OH_rad_p_) + '\n\n' + str(
                                             HO2_rad_p_) + '\n\n' + str(H_plus_p_) + '\n\n' + str(OH_minus_p_), font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0, anchor='n')

                        Dose_lab = Label(bottom_frame33, text="Per gram \n\n" + str(H2_g_) + "\n\n" + str(H2O2_g_) + "\n\n" + str(
                            e_aq_g_) + "\n\n" + str(H_rad_g_) + '\n\n' + str(OH_rad_g_) + '\n\n' + str(HO2_rad_g_) + '\n\n' + str(
                            H_plus_g_) + '\n\n' + str(OH_minus_g_), font=(None, 13))
                        Dose_lab.place(relx=0.5, rely=0, anchor='n')


def linear_results_window():
    global Ele, variable4, Iso_num

    if str(Ele.get()) == ' Isotope ':
        messagebox.showerror(
            title="Invalid input", message="Error:  Make sure you 'log' a valid an input before running")
    else:
        if variable4.get() == 'Unit':
            messagebox.showerror(
                title="Invalid input", message="Error:  Make sure you 'log' a valid an input before running")
        else:
            totE = np.sum(data_collector)
            if totE == 0:
                messagebox.showerror(title="Invalid input",
                                     message="  Error: Average Energy = 0, please correct this")
            else:
                totA = np.sum(activity_collector)
                if totA == 0:
                    messagebox.showerror(title="Invalid input",
                                         message="  Error: Total Activity = 0, please correct this")
                else:
                    def buttonpress():
                        save(planar, strin)

                    window = Toplevel(master)
                    layers = 1000
                    number = w1.get() * 1000
                    unit = variable6.get()
                    density = 10.97

                    canvas = Canvas(window, height=800, width=1400)
                    canvas.pack()

                    results_frame = Frame(window, bg='white', bd=0)
                    results_frame.place(relwidth=2, relheight=2, anchor='n')

                    # bd = boarder bg = background and colour
                    frame = Frame(window)
                    frame.place(relx=0.5, rely=0.05, relwidth=0.5,
                                relheight=0.05, anchor='n')

                    title = Label(frame, bg='white',
                                  text="Planar Report", font=(None, 30, 'bold'))
                    title.place(relwidth=1, relheight=1)

                    crack = var1.get()
                    scale = 1000
                    width = w2.get()
                    fx = np.zeros((scale))
                    e = range_float(0, (Iso_num-1), 1)

                    for l in e:
                        fx_topup = main2.dose(
                            data_collector[l], activity_collector[l], number, scale)
                        fx = fx + fx_topup

                    if crack == 1:

                        c = range_float(0, (Iso_num-1), 1)
                        d = range_float(0, (Iso_num-1), 1)
                        layers2 = width*10
                        fx_ = np.zeros((layers2+1))

                        for m in c:
                            fx_topup = maincrack.dose(
                                data_collector[m], activity_collector[m], number, width)
                            fx_ = fx_ + fx_topup

                        for a in d:
                            fx_topup = main2.dose(
                                data_collector[a], activity_collector[a], number, (layers2+1))
                            fx_ = fx_ + fx_topup

                    Activity = np.sum(activity_collector)

                    sc = range_float(0, Iso_num, 1)
                    Scaling = np.zeros(Iso_num+1)
                    for s in sc:
                        Scaling[s] = data_collector[s] * activity_collector[s]

                    Tot = np.sum(fx)
                    En = np.sum(Scaling)/Activity

                    if str(unit) == 'Gy/h':
                        Tot = Tot*3600
                        fx = fx*3600

                    Avg = Tot / 300
                    Tot_ = sci_notation(Tot, 3)
                    Avg_ = sci_notation(Avg, 3)
                    Activity_ = sci_notation(Activity, 3)
                    E_ = round(En, 2)
                    results_frame = Frame(window, bg='white', bd=0)
                    results_frame.place(
                        relx=0.5, rely=0.1, relwidth=0.5, relheight=0.6, anchor='n')

                    fig = Figure(figsize=(7, 7), dpi=100)
                    t = np.arange(0, layers, 1)
                    x = t*scale / (10 * layers)
                    ax = fig.add_subplot(111)
                    ax.grid(which='both', color="white",
                            linestyle='-', linewidth=0.5)
                    ax.plot(x, fx, label='E=' + str(E_) +
                            'MeV', color='Blue', lw=1, )
                    ax.legend(loc='upper right')
                    ax.set(xlabel='Distance from Interface $\mu$m', ylabel='Dose rate (' + str(unit) + ')',
                           title='Dose rate profile')

                    planar = np.column_stack([x, fx])
                    strin = str("Planar_Activity=" +
                                str(Activity) + "AvgEnergy=" + str(E_))

                    if crack == 1:
                        fig = Figure(figsize=(7, 7), dpi=100)
                        t = np.arange(0, (layers2+1), 1)
                        x = t * width / (layers2+1)
                        ax = fig.add_subplot(111)
                        ax.grid(which='both', color="white",
                                linestyle='-', linewidth=0.5)
                        ax.plot(x, fx_, label='E=' + str(E_) +
                                'MeV', color='Blue', lw=1, )
                        ax.legend(loc='upper right')
                        ax.set(xlabel='Distance from Interface $\mu$m', ylabel='Dose rate (' + str(unit) + ')',
                               title='Dose rate profile')
                        ax.set_ylim(ymin=0)
                        planar = np.column_stack([x, fx_])
                        strin = str("Crack_Activity=" +
                                    str(Activity) + "AvgEnergy=" + str(E_))

                    # A tk.DrawingArea.
                    canvas1 = FigureCanvasTkAgg(fig, master=results_frame)
                    toolbar = NavigationToolbar2Tk(canvas1, results_frame)
                    toolbar.config(background='white')
                    toolbar._message_label.config(background='white')
                    toolbar.update()
                    canvas1.draw()
                    canvas1.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)

                    # bd = boarder bg = background and colour
                    bottom_frame = Frame(window, bg='white')
                    bottom_frame.place(relx=0.15, rely=0.15,
                                       relwidth=0.2, relheight=0.55, anchor='n')

                    # bd = boarder bg = background and colour
                    bottom_frame2 = Frame(window, bg='white')
                    bottom_frame2.place(
                        relx=0.85, rely=0.15, relwidth=0.2, relheight=0.55, anchor='n')

                    # bd = boarder bg = background and colour
                    bottom_frame31 = Frame(window, bg='white')
                    bottom_frame31.place(
                        relx=0.5, rely=0.7, relwidth=0.9, relheight=0.25, anchor='n')

                    # bd = boarder bg = background and colout
                    sav_frame = Frame(window, bg="white", bd=0)
                    sav_frame.place(relx=0.7, rely=0.057,
                                    width=100, height=30, anchor='n')

                    sbutton = Button(sav_frame, bg="white",
                                     text='Save', command=buttonpress)
                    sbutton.pack()

                    title_lab = Label(bottom_frame, bg='white',
                                      text="Properties", font=(None, 16, 'bold'))
                    title_lab.place(relx=0.5, rely=0.05, anchor='n')

                    f_lab = Label(
                        bottom_frame, bg='white', text="Average decay energy = " + str(E_) + " MeV", font=(None, 13))
                    f_lab.place(relx=0.5, rely=0.2, anchor='n')

                    f_lab = Label(bottom_frame, bg='white', text="Activity = " +
                                  str(Activity_) + " Bqg\u207B\u00B9", font=(None, 13))
                    f_lab.place(relx=0.5, rely=0.3, anchor='n')

                    d_lab = Label(bottom_frame, bg='white', text="Density = " +
                                  str(density) + " gcm\u207B\u00B3", font=(None, 13))
                    d_lab.place(relx=0.5, rely=0.4, anchor='n')

                    d_lab = Label(
                        bottom_frame, bg='white', text="Probability of \n alpha escape = " + str(0.25), font=(None, 13))
                    d_lab.place(relx=0.5, rely=0.55, anchor='n')

                    d_lab = Label(bottom_frame, bg='white', text="Total no. of alpha \n decays simulated = " + str(number) + "\n (before scaling)",
                                  font=(None, 13))
                    d_lab.place(relx=0.5, rely=0.7, anchor='n')

                    title2_lab = Label(
                        bottom_frame31, bg='white', text="Dosimetry Results", font=(None, 16, 'bold'))
                    title2_lab.place(relx=0.5, rely=0.35, anchor='n')

                    title2_lab = Label(bottom_frame31, bg='white', text="Total dose rate from the interface = " + str(Tot_) + " " + str(unit),
                                       font=(None, 13))
                    title2_lab.place(relx=0.5, rely=0.55, anchor='n')

                    title2_lab = Label(bottom_frame31, bg='white',
                                       text="Average dose rate within the 30 \u03BCm from the interface = " +
                                       str(Avg_) + " " + str(unit),
                                       font=(None, 13))
                    title2_lab.place(relx=0.5, rely=0.75, anchor='n')

                    Pref = variable4.get()

                    print(str(Pref))

                    # bd = boarder bg = background and colour
                    getridofgrey_frame = Frame(window, bg='white')
                    getridofgrey_frame.place(
                        relx=0.747, rely=0.648, relwidth=0.01, relheight=0.05, anchor='n')

                    # RADIOLYSIS DATA (micro mol per J)

                    H2 = 0.1284
                    H2O2 = 0.104
                    e_aq = 0.01
                    H_rad = 0.0104
                    OH_rad = 0.0364
                    HO2_rad = 0.0104
                    H_plus = 0.01872
                    OH_minus = 0.00312

                    # PER GRAM

                    H2_g = 0.1284 * Tot
                    H2O2_g = 0.104 * Tot
                    e_aq_g = 0.01 * Tot
                    H_rad_g = 0.0104 * Tot
                    OH_rad_g = 0.0364 * Tot
                    HO2_rad_g = 0.0104 * Tot
                    H_plus_g = 0.01872 * Tot
                    OH_minus_g = 0.00312 * Tot

                    H2_g_ = sci_notation(H2_g, 3)
                    H2O2_g_ = sci_notation(H2O2_g, 3)
                    e_aq_g_ = sci_notation(e_aq_g, 3)
                    H_rad_g_ = sci_notation(H_rad_g, 3)
                    OH_rad_g_ = sci_notation(OH_rad_g, 3)
                    HO2_rad_g_ = sci_notation(HO2_rad_g, 3)
                    H_plus_g_ = sci_notation(H_plus_g, 3)
                    OH_minus_g_ = sci_notation(OH_minus_g, 3)

                    title2_lab = Label(
                        bottom_frame2, bg='white', text="Radiolysis Rates", font=(None, 16, 'bold'))
                    title2_lab.place(relx=0.5, rely=0, anchor='n')

                    title22_lab = Label(bottom_frame2, bg='white',
                                        text="G-values taken from \n Pastina et al (2001) \n at 25\u00B0C (\u03BCmolL\u207B\u00B9)",
                                        font=(None, 13,))
                    title22_lab.place(relx=0.5, rely=0.1, anchor='n')

                    Dose_lab = Label(bottom_frame2, bg='white',
                                     text="H\u2082 \n\n H\u2082O\u2082 \n\n e(aq) \n\n H(rad) \n\n OH(rad) \n\n HO\u2082(rad) \n\n H+ \n\n OH- ",
                                     font=(None, 13,))
                    Dose_lab.place(relx=0.3, rely=0.3, anchor='n')

                    Dose_lab = Label(bottom_frame2, bg='white',
                                     text=str(H2_g_) + "\n\n" + str(H2O2_g_) + "\n\n" + str(e_aq_g_) + "\n\n" + str(
                                         H_rad_g_) + '\n\n' + str(OH_rad_g_) + '\n\n' + str(HO2_rad_g_) + '\n\n' + str(
                                         H_plus_g_) + '\n\n' + str(OH_minus_g_),font=(None, 13,))
                    Dose_lab.place(relx=0.7, rely=0.3, anchor='n')

                    def on_key_press(event):
                        print("you pressed {}".format(event.key))
                        key_press_handler(event, canvas, toolbar)

                    canvas1.mpl_connect("key_press_event", on_key_press)


def _quit():
    master.quit()  # stops mainloop
    master.destroy()  # this is necessary on Windows to prevent
    # Fatal Python Error: PyEval_RestoreThread: NULL tstate


def combine_funcs(*funcs):
    def combined_func(*args, **kwargs):
        for f in funcs:
            f(*args, **kwargs)
    return combined_func


def linear_window():
    global w1, w3, variable4, variable6, var1, w2, Iso_num, data_collector, activity_collector, linear_results_window, background_image1, text_area, Ele

    linwindow = Toplevel(master)
    Ele = ' Isotope '
    variable4 = 'Unit'
    Iso_num = 0
    data_collector = np.zeros(100)
    activity_collector = np.zeros(100)

    def show_values(): #List of isotopes and assosiated energy
        global Iso_num, text_area, data_collector

        if str(Ele.get()) == ' Ac-225  ':
            data_collector[Iso_num] = 5.78676497633266
        if str(Ele.get()) == ' Ac-227  ':
            data_collector[Iso_num] = 4.93030222712238
        if str(Ele.get()) == ' Am-241  ':
            data_collector[Iso_num] = 5.47905653355526
        if str(Ele.get()) == ' Am-242m ':
            data_collector[Iso_num] = 5.20757734113712
        if str(Ele.get()) == ' Am-243  ':
            data_collector[Iso_num] = 5.2700509632566
        if str(Ele.get()) == ' At-211  ':
            data_collector[Iso_num] = 5.86692757826832
        if str(Ele.get()) == ' At-217  ':
            data_collector[Iso_num] = 7.065707158
        if str(Ele.get()) == ' Bi-211  ':
            data_collector[Iso_num] = 6.56708094999348
        if str(Ele.get()) == ' Bi-212  ':
            data_collector[Iso_num] = 6.05084850736504
        if str(Ele.get()) == ' Bi-213  ':
            data_collector[Iso_num] = 5.84625194444445
        if str(Ele.get()) == ' Cf-248  ':
            data_collector[Iso_num] = 6.25319986399728
        if str(Ele.get()) == ' Cf-249  ':
            data_collector[Iso_num] = 5.83168542760296
        if str(Ele.get()) == ' Cf-250  ':
            data_collector[Iso_num] = 6.02359330476776
        if str(Ele.get()) == ' Cf-251  ':
            data_collector[Iso_num] = 5.78330637254902
        if str(Ele.get()) == ' Cf-252  ':
            data_collector[Iso_num] = 6.11127538893686
        if str(Ele.get()) == ' Cf-253  ':
            data_collector[Iso_num] = 5.97593161290323
        if str(Ele.get()) == ' Cf-254  ':
            data_collector[Iso_num] = 5.82686
        if str(Ele.get()) == ' Cm-242  ':
            data_collector[Iso_num] = 6.10162420344283
        if str(Ele.get()) == ' Cm-243  ':
            data_collector[Iso_num] = 5.8086536484323
        if str(Ele.get()) == ' Cm-244  ':
            data_collector[Iso_num] = 5.79499868034311
        if str(Ele.get()) == ' Cm-245  ':
            data_collector[Iso_num] = 5.36066984769848
        if str(Ele.get()) == ' Cm-246  ':
            data_collector[Iso_num] = 5.37696980214856
        if str(Ele.get()) == ' Cm-247  ':
            data_collector[Iso_num] = 4.946722
        if str(Ele.get()) == ' Cm-248  ':
            data_collector[Iso_num] = 5.07064524844216
        if str(Ele.get()) == ' Cm-250  ':
            data_collector[Iso_num] = 5.19
        if str(Ele.get()) == ' Es-253  ':
            data_collector[Iso_num] = 6.62720945951081
        if str(Ele.get()) == ' Es-254  ':
            data_collector[Iso_num] = 6.42347304845296
        if str(Ele.get()) == ' Es-254m ':
            data_collector[Iso_num] = 6.39970888147923
        if str(Ele.get()) == ' Es-255  ':
            data_collector[Iso_num] = 6.2934665
        if str(Ele.get()) == ' Fm-254  ':
            data_collector[Iso_num] = 7.18186174085926
        if str(Ele.get()) == ' Fm-255  ':
            data_collector[Iso_num] = 7.01783615847222
        if str(Ele.get()) == ' Fm-256  ':
            data_collector[Iso_num] = 6.915
        if str(Ele.get()) == ' Fr-221  ':
            data_collector[Iso_num] = 6.30282713881863
        if str(Ele.get()) == ' Gd-152  ':
            data_collector[Iso_num] = 2.1496
        if str(Ele.get()) == ' Np-237  ':
            data_collector[Iso_num] = 4.76842282937581
        if str(Ele.get()) == ' Os-186  ':
            data_collector[Iso_num] = 2.7564
        if str(Ele.get()) == ' Pa-231  ':
            data_collector[Iso_num] = 4.93788777037308
        if str(Ele.get()) == ' Po-209  ':
            data_collector[Iso_num] = 4.88049493776203
        if str(Ele.get()) == ' Po-210  ':
            data_collector[Iso_num] = 5.30449141450859
        if str(Ele.get()) == ' Po-211  ':
            data_collector[Iso_num] = 7.44262767827678
        if str(Ele.get()) == ' Po-212  ':
            data_collector[Iso_num] = 8.7849
        if str(Ele.get()) == ' Po-213  ':
            data_collector[Iso_num] = 8.3769694
        if str(Ele.get()) == ' Po-214  ':
            data_collector[Iso_num] = 7.68701576146305
        if str(Ele.get()) == ' Po-215  ':
            data_collector[Iso_num] = 7.386157912
        if str(Ele.get()) == ' Po-216  ':
            data_collector[Iso_num] = 6.77848571697143
        if str(Ele.get()) == ' Po-218  ':
            data_collector[Iso_num] = 6.00249096170138
        if str(Ele.get()) == ' Pu-236  ':
            data_collector[Iso_num] = 5.75449315725211
        if str(Ele.get()) == ' Pu-238  ':
            data_collector[Iso_num] = 5.48696511658139
        if str(Ele.get()) == ' Pu-239  ':
            data_collector[Iso_num] = 5.14750943911272
        if str(Ele.get()) == ' Pu-240  ':
            data_collector[Iso_num] = 5.1563360143253
        if str(Ele.get()) == ' Pu-242  ':
            data_collector[Iso_num] = 4.89058114167105
        if str(Ele.get()) == ' Pu-244  ':
            data_collector[Iso_num] = 4.58053341647129
        if str(Ele.get()) == ' Ra-222  ':
            data_collector[Iso_num] = 6.54611374428159
        if str(Ele.get()) == ' Ra-223  ':
            data_collector[Iso_num] = 5.66395904817902
        if str(Ele.get()) == ' Ra-224  ':
            data_collector[Iso_num] = 5.67390446641391
        if str(Ele.get()) == ' Ra-226  ':
            data_collector[Iso_num] = 4.7743400002797
        if str(Ele.get()) == ' Rn-218  ':
            data_collector[Iso_num] = 7.13224054
        if str(Ele.get()) == ' Rn-219  ':
            data_collector[Iso_num] = 6.75459109927581
        if str(Ele.get()) == ' Rn-220  ':
            data_collector[Iso_num] = 6.287774939
        if str(Ele.get()) == ' Rn-222  ':
            data_collector[Iso_num] = 5.48930458956884
        if str(Ele.get()) == ' Sm-147  ':
            data_collector[Iso_num] = 2.2476
        if str(Ele.get()) == ' Th-226  ':
            data_collector[Iso_num] = 6.310050565526
        if str(Ele.get()) == ' Th-227  ':
            data_collector[Iso_num] = 5.88334404178184
        if str(Ele.get()) == ' Th-228  ':
            data_collector[Iso_num] = 5.3998401340134
        if str(Ele.get()) == ' Th-229  ':
            data_collector[Iso_num] = 4.87212365602172
        if str(Ele.get()) == ' Th-230  ':
            data_collector[Iso_num] = 4.67071594926595
        if str(Ele.get()) == ' Th-232  ':
            data_collector[Iso_num] = 3.99655688622755
        if str(Ele.get()) == ' U-230  ':
            data_collector[Iso_num] = 5.86439142396811
        if str(Ele.get()) == ' U-232  ':
            data_collector[Iso_num] = 5.30207449527094
        if str(Ele.get()) == ' U-233  ':
            data_collector[Iso_num] = 4.81667520161068
        if str(Ele.get()) == ' U-234  ':
            data_collector[Iso_num] = 4.76111943443504
        if str(Ele.get()) == ' U-235  ':
            data_collector[Iso_num] = 4.39162387161484
        if str(Ele.get()) == ' U-236  ':
            data_collector[Iso_num] = 4.48087293038101
        if str(Ele.get()) == ' U-238  ':
            data_collector[Iso_num] = 4.18439559014267
        if str(Ele.get()) == ' Avg E ':
            data_collector[Iso_num] = (AE.get()/100)

        Pref = variable4.get()

        if str(Ele.get()) == ' Isotope ':
            messagebox.showerror(title="Invalid input",
                                 message="  Error:  Please select an option for Isotope")

        else:
            if str(Pref) == 'Unit':
                messagebox.showerror(title="Invalid input",
                                     message="  Error:  Please select an option for Unit")
            else:

                if str(Pref) == 'GBqg\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 1e9)

                if str(Pref) == 'MBqg\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 1e6)

                if str(Pref) == 'kBqg\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 1e3)

                if str(Pref) == 'Bqg\u207B\u00B9':
                    activity_collector[Iso_num] = w3.get()

                if str(Pref) == 'Cig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37 * 1e9)

                if str(Pref) == 'mCig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37 * 1e6)

                if str(Pref) == '\u03BCCig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37 * 1e3)

                if str(Pref) == 'nCig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37)


                Iso_num = Iso_num + 1

                if str(Ele.get()) == ' Avg E ':
                    text_area.insert(END, "    Avg E   =    " + str(
                        (AE.get() / 100)) + " MeV" + "          Activity      " + str(Iso_num) + "   =   " + str(
                        w3.get()) + " " + str(
                        variable4.get()) + " \n")
                    if AE.get() == 0:
                        messagebox.showerror(
                            title="Invalid input", message="  Error:    Please select an Energy > 0    \n (Correct and press 'log' to overwrite)")
                        Iso_num = Iso_num - 1

                else:
                    text_area.insert(END, "      Isotope    " + str(Iso_num) + "  =    " + str(
                        Ele.get()) + "          Activity    " + str(Iso_num) + "   =   " + str(w3.get()) + " " + str(
                        variable4.get()) + " \n")

                if w3.get() == 0:
                    messagebox.showerror(
                        title="Invalid input", message="  Error:    Please select an Activity > 0 \n    (Correct and press 'log' to overwrite)")
                    Iso_num = Iso_num - 1

        print("Energies = ", data_collector)
        print("Activities = ", activity_collector)

    canvas_ = Canvas(linwindow, height=860, width=830)
    canvas_.pack()

    background_label1 = Label(linwindow, image=background_image1)
    background_label1.place(relwidth=1, relheight=1,
                            relx=0.5, rely=0, anchor='n')

    lower_frame = Frame(linwindow, bg="white", bd=1)
    lower_frame.place(relx=0.5, rely=0.22, relwidth=0.65,
                      relheight=0.72, anchor='n')

    log_frame = Frame(linwindow, bg="white", bd=0)
    log_frame.place(relx=0.52, rely=0.45, relwidth=0.6,
                    relheight=0.15, anchor='n')

    text_area = scrolledtext.ScrolledText(log_frame,
                                          width=50,
                                          height=5,
                                          font=(None,
                                                12))

    text_area.grid(column=0, pady=0, padx=0)

    label = Label(lower_frame, bg="white")
    label.place(relwidth=1, relheight=1, anchor='n')

    ELEMENT_OPTIONS = [
        " Avg E ",
        " Ac-225  ",
        " Ac-227  ",
        " Am-241  ",
        " Am-242m ",
        " Am-243  ",
        " At-211  ",
        " At-217  ",
        " Bi-211  ",
        " Bi-212  ",
        " Bi-213  ",
        " Cf-248  ",
        " Cf-249  ",
        " Cf-250  ",
        " Cf-251  ",
        " Cf-252  ",
        " Cf-253  ",
        " Cf-254  ",
        " Cm-242  ",
        " Cm-243  ",
        " Cm-244  ",
        " Cm-245  ",
        " Cm-246  ",
        " Cm-247  ",
        " Cm-248  ",
        " Cm-250  ",
        " Es-253  ",
        " Es-254  ",
        " Es-254m ",
        " Es-255  ",
        " Fm-254  ",
        " Fm-255  ",
        " Fm-256  ",
        " Fr-221  ",
        " Gd-152  ",
        " Np-237  ",
        " Os-186  ",
        " Pa-231  ",
        " Po-209  ",
        " Po-210  ",
        " Po-211  ",
        " Po-212  ",
        " Po-213  ",
        " Po-214  ",
        " Po-215  ",
        " Po-216  ",
        " Po-218  ",
        " Pu-236  ",
        " Pu-238  ",
        " Pu-239  ",
        " Pu-240  ",
        " Pu-242  ",
        " Pu-244  ",
        " Ra-222  ",
        " Ra-223  ",
        " Ra-224  ",
        " Ra-226  ",
        " Rn-218  ",
        " Rn-219  ",
        " Rn-220  ",
        " Rn-222  ",
        " Sm-147  ",
        " Th-226  ",
        " Th-227  ",
        " Th-228  ",
        " Th-229  ",
        " Th-230  ",
        " Th-232  ",
        " U-230  ",
        " U-232  ",
        " U-233  ",
        " U-234  ",
        " U-235  ",
        " U-236  ",
        " U-238  "
    ]

    Ele = StringVar(lower_frame)
    Ele.set(" Isotope ")

    w = ttk.Combobox(lower_frame, textvariable=Ele, values=[
                   ELEMENT_OPTIONS[0],
                   ELEMENT_OPTIONS[1],
                   ELEMENT_OPTIONS[2],
                   ELEMENT_OPTIONS[3],
                   ELEMENT_OPTIONS[4],
                   ELEMENT_OPTIONS[5],
                   ELEMENT_OPTIONS[6],
                   ELEMENT_OPTIONS[7],
                   ELEMENT_OPTIONS[8],
                   ELEMENT_OPTIONS[9],
                   ELEMENT_OPTIONS[10],
                   ELEMENT_OPTIONS[11],
                   ELEMENT_OPTIONS[12],
                   ELEMENT_OPTIONS[13],
                   ELEMENT_OPTIONS[14],
                   ELEMENT_OPTIONS[15],
                   ELEMENT_OPTIONS[16],
                   ELEMENT_OPTIONS[17],
                   ELEMENT_OPTIONS[18],
                   ELEMENT_OPTIONS[19],
                   ELEMENT_OPTIONS[20],
                   ELEMENT_OPTIONS[21],
                   ELEMENT_OPTIONS[22],
                   ELEMENT_OPTIONS[23],
                   ELEMENT_OPTIONS[24],
                   ELEMENT_OPTIONS[25],
                   ELEMENT_OPTIONS[26],
                   ELEMENT_OPTIONS[27],
                   ELEMENT_OPTIONS[28],
                   ELEMENT_OPTIONS[29],
                   ELEMENT_OPTIONS[30],
                   ELEMENT_OPTIONS[31],
                   ELEMENT_OPTIONS[32],
                   ELEMENT_OPTIONS[33],
                   ELEMENT_OPTIONS[34],
                   ELEMENT_OPTIONS[35],
                   ELEMENT_OPTIONS[36],
                   ELEMENT_OPTIONS[37],
                   ELEMENT_OPTIONS[38],
                   ELEMENT_OPTIONS[39],
                   ELEMENT_OPTIONS[40],
                   ELEMENT_OPTIONS[41],
                   ELEMENT_OPTIONS[42],
                   ELEMENT_OPTIONS[43],
                   ELEMENT_OPTIONS[44],
                   ELEMENT_OPTIONS[45],
                   ELEMENT_OPTIONS[46],
                   ELEMENT_OPTIONS[47],
                   ELEMENT_OPTIONS[48],
                   ELEMENT_OPTIONS[49],
                   ELEMENT_OPTIONS[50],
                   ELEMENT_OPTIONS[51],
                   ELEMENT_OPTIONS[52],
                   ELEMENT_OPTIONS[53],
                   ELEMENT_OPTIONS[54],
                   ELEMENT_OPTIONS[55],
                   ELEMENT_OPTIONS[56],
                   ELEMENT_OPTIONS[57],
                   ELEMENT_OPTIONS[58],
                   ELEMENT_OPTIONS[59],
                   ELEMENT_OPTIONS[60],
                   ELEMENT_OPTIONS[61],
                   ELEMENT_OPTIONS[62],
                   ELEMENT_OPTIONS[63],
                   ELEMENT_OPTIONS[64],
                   ELEMENT_OPTIONS[65],
                   ELEMENT_OPTIONS[66],
                   ELEMENT_OPTIONS[67],
                   ELEMENT_OPTIONS[68],
                   ELEMENT_OPTIONS[69],
                   ELEMENT_OPTIONS[70],
                   ELEMENT_OPTIONS[71],
                   ELEMENT_OPTIONS[72],
                   ELEMENT_OPTIONS[73],
                   ELEMENT_OPTIONS[74]])
    w.place(relx=0.1, rely=0.08, width=100, anchor='n')

    w1 = Scale(lower_frame, bg="white", from_=1, to=100,
               tickinterval=99, orient=HORIZONTAL)
    w1.set(1)
    w1.place(relx=0.5, rely=0.78, relwidth=0.7, anchor='n')

    text = Label(lower_frame, bg="white",
                 text="Move the Slider below to choose how many decays to simulate \n (1=1000 decays)", font=(None, 12,))
    text.place(relx=0.5, rely=0.7, relwidth=1, anchor='n')

    entry_frame = Frame(linwindow, bg="white")
    entry_frame.place(relx=0.5, rely=0.57, relwidth=0.7,
                      relheight=0.12, anchor='n')

    text = Label(lower_frame, bg="white",
                 text="If Isotope = Avg E use slider to choose \n average decay energy (1 = 10 GeV)", font=(None, 12))
    text.place(relx=0.25, rely=0.17, relwidth=1, anchor='n')

    AE = Scale(lower_frame, bg="white", from_=0, to=1000,
               tickinterval=500, orient=HORIZONTAL)
    AE.set(1)
    AE.place(relx=0.75, rely=0.15, relwidth=0.4, anchor='n')

    text2 = Label(entry_frame, bg="white",
                  text="Crack width (\u03BCm): ", font=(None, 12,))
    text2.place(relx=0.32, rely=0.18, width=150, anchor='n')

    var1 = IntVar()
    Crack = Checkbutton(entry_frame, bg="gray",bd=0.1, variable=var1)
    Crack.place(relx=0.15, rely=0.18, width=20, anchor='n')

    w2 = Scale(entry_frame, bg="white", from_=0, to=100,
               tickinterval=20, orient=HORIZONTAL)
    w2.set(1)
    w2.place(relx=0.7, rely=0.02, relwidth=0.5, anchor='n')

    PREFIX_OPTIONS = [
        "GBqg\u207B\u00B9",
        "MBqg\u207B\u00B9",
        "kBqg\u207B\u00B9",
        "Bqg\u207B\u00B9",
        "Cig\u207B\u00B9",
        "mCig\u207B\u00B9",
        "\u03BCCig\u207B\u00B9",
        "nCig\u207B\u00B9",
    ]

    variable4 = StringVar(linwindow)
    variable4.set("Unit")  # default value

    w4 = ttk.Combobox(lower_frame, textvariable=variable4, values=[PREFIX_OPTIONS[0], PREFIX_OPTIONS[1],  PREFIX_OPTIONS[2],
                    PREFIX_OPTIONS[3],  PREFIX_OPTIONS[4], PREFIX_OPTIONS[5],  PREFIX_OPTIONS[6], PREFIX_OPTIONS[7]])
    w4.place(relx=0.9, rely=0.08, width=100, anchor='n')

    w3 = Scale(lower_frame, bg="white", from_=0, to=1000,
               tickinterval=200, orient=HORIZONTAL)
    w3.set(500)
    w3.place(relx=0.5, rely=0.05, relwidth=0.6, anchor='n')

    text3 = Label(lower_frame, bg="white",
                  text="Enter Isotope, Activity and Unit then press 'Log' to save entry:", font=(None, 12))
    text3.place(relx=0.5, rely=0, width=600, anchor='n')

    UNIT_OPTIONS = [
        "Gy/s",
        "Gy/h"
    ]

    variable6 = StringVar(linwindow)
    variable6.set(UNIT_OPTIONS[0])  # default value

    text2 = Label(entry_frame, bg="white", text="Unit:", font=(None, 12))
    text2.place(relx=0.4, rely=0.7, width=100, anchor='n')

    w6 = ttk.Combobox(entry_frame, textvariable=variable6, values=[UNIT_OPTIONS[0], UNIT_OPTIONS[1]])
    w6.place(relx=0.6, rely=0.7, width=100, anchor='n')

    button2 = Button(lower_frame, bg="white", text='Log', command=show_values)
    button2.place(relx=0.5, rely=0.26, anchor='n')

    button = Button(lower_frame, bg="white", text='Solve',
                    command=linear_results_window)
    button.place(relx=0.5, rely=0.88, anchor='n')

    text2 = Label(linwindow, bg="white",
                  text="Version 2.1, \u00A9 2022, Angus Siberry. All rights reserved.", font=(None, 8))
    text2.place(relx=0.7, rely=0.9, width=300, anchor='n')


def spherical_window():
    global w1, w3, variable4, variable5, variable6, var1, w2, Iso_num, data_collector, activity_collector, linear_results_window, background_image1, text_area, size_image, Ele

    swindow = Toplevel(master)
    Ele = ' Isotope '
    variable4 = 'Unit'
    Iso_num = 0
    data_collector = np.zeros(100)
    activity_collector = np.zeros(100)

    def show_values():
        global Iso_num, text_area, data_collector

        if str(Ele.get()) == ' Ac-225  ':
            data_collector[Iso_num] = 5.78676497633266
        if str(Ele.get()) == ' Ac-227  ':
            data_collector[Iso_num] = 4.93030222712238
        if str(Ele.get()) == ' Am-241  ':
            data_collector[Iso_num] = 5.47905653355526
        if str(Ele.get()) == ' Am-242m ':
            data_collector[Iso_num] = 5.20757734113712
        if str(Ele.get()) == ' Am-243  ':
            data_collector[Iso_num] = 5.2700509632566
        if str(Ele.get()) == ' At-211  ':
            data_collector[Iso_num] = 5.86692757826832
        if str(Ele.get()) == ' At-217  ':
            data_collector[Iso_num] = 7.065707158
        if str(Ele.get()) == ' Bi-211  ':
            data_collector[Iso_num] = 6.56708094999348
        if str(Ele.get()) == ' Bi-212  ':
            data_collector[Iso_num] = 6.05084850736504
        if str(Ele.get()) == ' Bi-213  ':
            data_collector[Iso_num] = 5.84625194444445
        if str(Ele.get()) == ' Cf-248  ':
            data_collector[Iso_num] = 6.25319986399728
        if str(Ele.get()) == ' Cf-249  ':
            data_collector[Iso_num] = 5.83168542760296
        if str(Ele.get()) == ' Cf-250  ':
            data_collector[Iso_num] = 6.02359330476776
        if str(Ele.get()) == ' Cf-251  ':
            data_collector[Iso_num] = 5.78330637254902
        if str(Ele.get()) == ' Cf-252  ':
            data_collector[Iso_num] = 6.11127538893686
        if str(Ele.get()) == ' Cf-253  ':
            data_collector[Iso_num] = 5.97593161290323
        if str(Ele.get()) == ' Cf-254  ':
            data_collector[Iso_num] = 5.82686
        if str(Ele.get()) == ' Cm-242  ':
            data_collector[Iso_num] = 6.10162420344283
        if str(Ele.get()) == ' Cm-243  ':
            data_collector[Iso_num] = 5.8086536484323
        if str(Ele.get()) == ' Cm-244  ':
            data_collector[Iso_num] = 5.79499868034311
        if str(Ele.get()) == ' Cm-245  ':
            data_collector[Iso_num] = 5.36066984769848
        if str(Ele.get()) == ' Cm-246  ':
            data_collector[Iso_num] = 5.37696980214856
        if str(Ele.get()) == ' Cm-247  ':
            data_collector[Iso_num] = 4.946722
        if str(Ele.get()) == ' Cm-248  ':
            data_collector[Iso_num] = 5.07064524844216
        if str(Ele.get()) == ' Cm-250  ':
            data_collector[Iso_num] = 5.19
        if str(Ele.get()) == ' Es-253  ':
            data_collector[Iso_num] = 6.62720945951081
        if str(Ele.get()) == ' Es-254  ':
            data_collector[Iso_num] = 6.42347304845296
        if str(Ele.get()) == ' Es-254m ':
            data_collector[Iso_num] = 6.39970888147923
        if str(Ele.get()) == ' Es-255  ':
            data_collector[Iso_num] = 6.2934665
        if str(Ele.get()) == ' Fm-254  ':
            data_collector[Iso_num] = 7.18186174085926
        if str(Ele.get()) == ' Fm-255  ':
            data_collector[Iso_num] = 7.01783615847222
        if str(Ele.get()) == ' Fm-256  ':
            data_collector[Iso_num] = 6.915
        if str(Ele.get()) == ' Fr-221  ':
            data_collector[Iso_num] = 6.30282713881863
        if str(Ele.get()) == ' Gd-152  ':
            data_collector[Iso_num] = 2.1496
        if str(Ele.get()) == ' Np-237  ':
            data_collector[Iso_num] = 4.76842282937581
        if str(Ele.get()) == ' Os-186  ':
            data_collector[Iso_num] = 2.7564
        if str(Ele.get()) == ' Pa-231  ':
            data_collector[Iso_num] = 4.93788777037308
        if str(Ele.get()) == ' Po-209  ':
            data_collector[Iso_num] = 4.88049493776203
        if str(Ele.get()) == ' Po-210  ':
            data_collector[Iso_num] = 5.30449141450859
        if str(Ele.get()) == ' Po-211  ':
            data_collector[Iso_num] = 7.44262767827678
        if str(Ele.get()) == ' Po-212  ':
            data_collector[Iso_num] = 8.7849
        if str(Ele.get()) == ' Po-213  ':
            data_collector[Iso_num] = 8.3769694
        if str(Ele.get()) == ' Po-214  ':
            data_collector[Iso_num] = 7.68701576146305
        if str(Ele.get()) == ' Po-215  ':
            data_collector[Iso_num] = 7.386157912
        if str(Ele.get()) == ' Po-216  ':
            data_collector[Iso_num] = 6.77848571697143
        if str(Ele.get()) == ' Po-218  ':
            data_collector[Iso_num] = 6.00249096170138
        if str(Ele.get()) == ' Pu-236  ':
            data_collector[Iso_num] = 5.75449315725211
        if str(Ele.get()) == ' Pu-238  ':
            data_collector[Iso_num] = 5.48696511658139
        if str(Ele.get()) == ' Pu-239  ':
            data_collector[Iso_num] = 5.14750943911272
        if str(Ele.get()) == ' Pu-240  ':
            data_collector[Iso_num] = 5.1563360143253
        if str(Ele.get()) == ' Pu-242  ':
            data_collector[Iso_num] = 4.89058114167105
        if str(Ele.get()) == ' Pu-244  ':
            data_collector[Iso_num] = 4.58053341647129
        if str(Ele.get()) == ' Ra-222  ':
            data_collector[Iso_num] = 6.54611374428159
        if str(Ele.get()) == ' Ra-223  ':
            data_collector[Iso_num] = 5.66395904817902
        if str(Ele.get()) == ' Ra-224  ':
            data_collector[Iso_num] = 5.67390446641391
        if str(Ele.get()) == ' Ra-226  ':
            data_collector[Iso_num] = 4.7743400002797
        if str(Ele.get()) == ' Rn-218  ':
            data_collector[Iso_num] = 7.13224054
        if str(Ele.get()) == ' Rn-219  ':
            data_collector[Iso_num] = 6.75459109927581
        if str(Ele.get()) == ' Rn-220  ':
            data_collector[Iso_num] = 6.287774939
        if str(Ele.get()) == ' Rn-222  ':
            data_collector[Iso_num] = 5.48930458956884
        if str(Ele.get()) == ' Sm-147  ':
            data_collector[Iso_num] = 2.2476
        if str(Ele.get()) == ' Th-226  ':
            data_collector[Iso_num] = 6.310050565526
        if str(Ele.get()) == ' Th-227  ':
            data_collector[Iso_num] = 5.88334404178184
        if str(Ele.get()) == ' Th-228  ':
            data_collector[Iso_num] = 5.3998401340134
        if str(Ele.get()) == ' Th-229  ':
            data_collector[Iso_num] = 4.87212365602172
        if str(Ele.get()) == ' Th-230  ':
            data_collector[Iso_num] = 4.67071594926595
        if str(Ele.get()) == ' Th-232  ':
            data_collector[Iso_num] = 3.99655688622755
        if str(Ele.get()) == ' U-230   ':
            data_collector[Iso_num] = 5.86439142396811
        if str(Ele.get()) == ' U-232   ':
            data_collector[Iso_num] = 5.30207449527094
        if str(Ele.get()) == ' U-233   ':
            data_collector[Iso_num] = 4.81667520161068
        if str(Ele.get()) == ' U-234   ':
            data_collector[Iso_num] = 4.76111943443504
        if str(Ele.get()) == ' U-235   ':
            data_collector[Iso_num] = 4.39162387161484
        if str(Ele.get()) == ' U-236   ':
            data_collector[Iso_num] = 4.48087293038101
        if str(Ele.get()) == ' U-238   ':
            data_collector[Iso_num] = 4.18439559014267
        if str(Ele.get()) == ' Avg E ':
            data_collector[Iso_num] = (AE.get()/100)

        Pref = variable4.get()

        if str(Ele.get()) == ' Isotope ':
            messagebox.showerror(title="Invalid input",
                                 message="  Error: Please select an option for Isotope")

        else:
            if str(Pref) == 'Unit':
                messagebox.showerror(title="Invalid input",
                                     message="   Error: Please select an option for Unit")
            else:

                if str(Pref) == 'GBqg\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 1e9)

                if str(Pref) == 'MBqg\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 1e6)

                if str(Pref) == 'kBqg\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 1e3)

                if str(Pref) == 'Bqg\u207B\u00B9':
                    activity_collector[Iso_num] = w3.get()

                if str(Pref) == 'Cig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37 * 1e9)

                if str(Pref) == 'mCig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37 * 1e6)

                if str(Pref) == '\u03BCCig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37 * 1e3)

                if str(Pref) == 'nCig\u207B\u00B9':
                    activity_collector[Iso_num] = (w3.get() * 37)

                Iso_num = Iso_num + 1

                if str(Ele.get()) == ' Avg E ':
                    text_area.insert(END, "               Avg E   =    " + str(
                        (AE.get() / 100)) + " MeV" + "                Activity      " + str(Iso_num) + "   =   " + str(
                        w3.get()) + " " + str(
                        variable4.get()) + " \n")
                    if AE.get() == 0:
                        messagebox.showerror(
                            title="Invalid input", message="  Error:       Please select an Energy > 0 \n    (Correct and press 'log' to overwrite)")
                        Iso_num = Iso_num - 1

                else:
                    text_area.insert(END, "               Isotope    " + str(Iso_num) + "  =    " + str(Ele.get()) +
                                     "                Activity      " + str(Iso_num) + "   =   " + str(w3.get()) + " " + str(variable4.get()) + " \n")

                if w3.get() == 0:
                    messagebox.showerror(
                        title="Invalid input", message="  Error:      Please select an Activity > 0 \n    (Correct and press 'log' to overwrite)")
                    if AE.get() == 0:
                        Iso_num = Iso_num
                    else:
                        Iso_num = Iso_num - 1

        print("Energies = ", data_collector)
        print("Activities = ", activity_collector)

    canvas_ = Canvas(swindow, height=860, width=1080)
    canvas_.pack()

    background_label2 = Label(swindow, image=background_image2)
    background_label2.place(relwidth=1, relheight=1,
                            relx=0.5, rely=0, anchor='n')

    lower_frame = Frame(swindow, bg="white", bd=1)
    lower_frame.place(relx=0.5, rely=0.22, relwidth=0.65,
                      relheight=0.88, anchor='n')

    #log_frame = Frame(swindow,bg="white", bd=0)
    #log_frame.place(relx=0.55, rely=0.4, relwidth=0.65, relheight=0.1, anchor='n')

    text_area = scrolledtext.ScrolledText(lower_frame,
                                          width=60,
                                          height=3,
                                          font=(None,
                                                12))

    text_area.grid(column=0, pady=0, padx=0)
    text_area.place(relx=0.52, rely=0.21, relwidth=0.8,
                    relheight=0.09, anchor='n')

    size_label = Label(swindow, bg='white', image=size_image)
    size_label.place(relx=0.5, rely=0.485, width=500, height=160, anchor='n')

    ELEMENT_OPTIONS = [
        " Avg E ",
        " Ac-225  ",
        " Ac-227  ",
        " Am-241  ",
        " Am-242m ",
        " Am-243  ",
        " At-211  ",
        " At-217  ",
        " Bi-211  ",
        " Bi-212  ",
        " Bi-213  ",
        " Cf-248  ",
        " Cf-249  ",
        " Cf-250  ",
        " Cf-251  ",
        " Cf-252  ",
        " Cf-253  ",
        " Cf-254  ",
        " Cm-242  ",
        " Cm-243  ",
        " Cm-244  ",
        " Cm-245  ",
        " Cm-246  ",
        " Cm-247  ",
        " Cm-248  ",
        " Cm-250  ",
        " Es-253  ",
        " Es-254  ",
        " Es-254m ",
        " Es-255  ",
        " Fm-254  ",
        " Fm-255  ",
        " Fm-256  ",
        " Fr-221  ",
        " Gd-152  ",
        " Np-237  ",
        " Os-186  ",
        " Pa-231  ",
        " Po-209  ",
        " Po-210  ",
        " Po-211  ",
        " Po-212  ",
        " Po-213  ",
        " Po-214  ",
        " Po-215  ",
        " Po-216  ",
        " Po-218  ",
        " Pu-236  ",
        " Pu-238  ",
        " Pu-239  ",
        " Pu-240  ",
        " Pu-242  ",
        " Pu-244  ",
        " Ra-222  ",
        " Ra-223  ",
        " Ra-224  ",
        " Ra-226  ",
        " Rn-218  ",
        " Rn-219  ",
        " Rn-220  ",
        " Rn-222  ",
        " Sm-147  ",
        " Th-226  ",
        " Th-227  ",
        " Th-228  ",
        " Th-229  ",
        " Th-230  ",
        " Th-232  ",
        " U-230  ",
        " U-232  ",
        " U-233  ",
        " U-234  ",
        " U-235  ",
        " U-236  ",
        " U-238  "
    ]

    Ele = StringVar(lower_frame)
    Ele.set(" Isotope ")

    w = ttk.Combobox(lower_frame, textvariable=Ele, values=[
                   ELEMENT_OPTIONS[0],
                   ELEMENT_OPTIONS[1],
                   ELEMENT_OPTIONS[2],
                   ELEMENT_OPTIONS[3],
                   ELEMENT_OPTIONS[4],
                   ELEMENT_OPTIONS[5],
                   ELEMENT_OPTIONS[6],
                   ELEMENT_OPTIONS[7],
                   ELEMENT_OPTIONS[8],
                   ELEMENT_OPTIONS[9],
                   ELEMENT_OPTIONS[10],
                   ELEMENT_OPTIONS[11],
                   ELEMENT_OPTIONS[12],
                   ELEMENT_OPTIONS[13],
                   ELEMENT_OPTIONS[14],
                   ELEMENT_OPTIONS[15],
                   ELEMENT_OPTIONS[16],
                   ELEMENT_OPTIONS[17],
                   ELEMENT_OPTIONS[18],
                   ELEMENT_OPTIONS[19],
                   ELEMENT_OPTIONS[20],
                   ELEMENT_OPTIONS[21],
                   ELEMENT_OPTIONS[22],
                   ELEMENT_OPTIONS[23],
                   ELEMENT_OPTIONS[24],
                   ELEMENT_OPTIONS[25],
                   ELEMENT_OPTIONS[26],
                   ELEMENT_OPTIONS[27],
                   ELEMENT_OPTIONS[28],
                   ELEMENT_OPTIONS[29],
                   ELEMENT_OPTIONS[30],
                   ELEMENT_OPTIONS[31],
                   ELEMENT_OPTIONS[32],
                   ELEMENT_OPTIONS[33],
                   ELEMENT_OPTIONS[34],
                   ELEMENT_OPTIONS[35],
                   ELEMENT_OPTIONS[36],
                   ELEMENT_OPTIONS[37],
                   ELEMENT_OPTIONS[38],
                   ELEMENT_OPTIONS[39],
                   ELEMENT_OPTIONS[40],
                   ELEMENT_OPTIONS[41],
                   ELEMENT_OPTIONS[42],
                   ELEMENT_OPTIONS[43],
                   ELEMENT_OPTIONS[44],
                   ELEMENT_OPTIONS[45],
                   ELEMENT_OPTIONS[46],
                   ELEMENT_OPTIONS[47],
                   ELEMENT_OPTIONS[48],
                   ELEMENT_OPTIONS[49],
                   ELEMENT_OPTIONS[50],
                   ELEMENT_OPTIONS[51],
                   ELEMENT_OPTIONS[52],
                   ELEMENT_OPTIONS[53],
                   ELEMENT_OPTIONS[54],
                   ELEMENT_OPTIONS[55],
                   ELEMENT_OPTIONS[56],
                   ELEMENT_OPTIONS[57],
                   ELEMENT_OPTIONS[58],
                   ELEMENT_OPTIONS[59],
                   ELEMENT_OPTIONS[60],
                   ELEMENT_OPTIONS[61],
                   ELEMENT_OPTIONS[62],
                   ELEMENT_OPTIONS[63],
                   ELEMENT_OPTIONS[64],
                   ELEMENT_OPTIONS[65],
                   ELEMENT_OPTIONS[66],
                   ELEMENT_OPTIONS[67],
                   ELEMENT_OPTIONS[68],
                   ELEMENT_OPTIONS[69],
                   ELEMENT_OPTIONS[70],
                   ELEMENT_OPTIONS[71],
                   ELEMENT_OPTIONS[72],
                   ELEMENT_OPTIONS[73],
                   ELEMENT_OPTIONS[74]])
    w.place(relx=0.1, rely=0.025, width=100, anchor='n')

    w1 = Scale(lower_frame, bg="white", from_=1, to=100,
               tickinterval=99, orient=HORIZONTAL)
    w1.set(1)
    w1.place(relx=0.5, rely=0.75, relwidth=0.7, anchor='n')

    text = Label(lower_frame, bg="white",
                 text="Move the slider below to choose how many decays to simulate \n (1=1000 decays)", font=(None, 12,))
    text.place(relx=0.5, rely=0.7, relwidth=1, anchor='n')

    text_ = Label(lower_frame, bg="white",
                  text="Move the slider below to the desired particle Radius (\u03BCm)", font=(None, 12,))
    text_.place(relx=0.5, rely=0.515, relwidth=1, anchor='n')

    w2 = Scale(lower_frame, bg="white", from_=0, to=50,
               tickinterval=10, orient=HORIZONTAL)
    w2.set(1)
    w2.place(relx=0.5, rely=0.55, relwidth=1, anchor='n')

    variable5 = StringVar(master)
    variable5.set("R > 50 \u03BCm")  # default value

    R_OPTIONS = [
        "R < 50\u03BCm",
        "100\u03BCm",
        "1,000\u03BCm",
        "10,000\u03BCm"
    ]

    w5 = ttk.Combobox(lower_frame, textvariable=variable5, values=[
                    R_OPTIONS[0], R_OPTIONS[1], R_OPTIONS[2], R_OPTIONS[3]])
    w5.place(relx=0.65, rely=0.65, width=100, anchor='n')

    PREFIX_OPTIONS = [
        "GBqg\u207B\u00B9",
        "MBqg\u207B\u00B9",
        "kBqg\u207B\u00B9",
        "Bqg\u207B\u00B9",
        "Cig\u207B\u00B9",
        "mCig\u207B\u00B9",
        "\u03BCCig\u207B\u00B9",
        "nCig\u207B\u00B9",
    ]

    variable4 = StringVar(swindow)
    variable4.set("Unit")  # default value

    w4 = ttk.Combobox(lower_frame, textvariable=variable4, values=[PREFIX_OPTIONS[0], PREFIX_OPTIONS[1],  PREFIX_OPTIONS[2],
                    PREFIX_OPTIONS[3],  PREFIX_OPTIONS[4], PREFIX_OPTIONS[5],  PREFIX_OPTIONS[6], PREFIX_OPTIONS[7]])
    w4.place(relx=0.9, rely=0.025, width=100, anchor='n')

    w3 = Scale(lower_frame, bg="white", from_=0, to=1000,
               tickinterval=200, orient=HORIZONTAL)
    w3.set(500)
    w3.place(relx=0.5, rely=0, relwidth=0.6, anchor='n')

    text3 = Label(swindow, bg="white",
                  text="Enter Isotope and Activity then press 'log' to save entry:", font=(None, 12))
    text3.place(relx=0.5, rely=0.17, width=400, anchor='n')

    text = Label(lower_frame, bg="white",
                 text="If Isotope = Avg E use slider to choose \n average decay energy (1 = 10 GeV)", font=(None, 12))
    text.place(relx=0.25, rely=0.09, relwidth=1, anchor='n')

    AE = Scale(lower_frame, bg="white", from_=0, to=1000,
               tickinterval=500, orient=HORIZONTAL)
    AE.set(1)
    AE.place(relx=0.75, rely=0.08, relwidth=0.5, anchor='n')

    UNIT_OPTIONS = [
        "Gy/s",
        "Gy/h"
    ]

    variable6 = StringVar(swindow)
    variable6.set(UNIT_OPTIONS[0])  # default value

    text2 = Label(lower_frame, bg="white", text="Unit:", font=(None, 12))
    text2.place(relx=0.25, rely=0.65, width=100, anchor='n')

    w6 = ttk.Combobox(lower_frame, textvariable=variable6, values=[UNIT_OPTIONS[0], UNIT_OPTIONS[1]])
    w6.place(relx=0.35, rely=0.65, width=100, anchor='n')

    button2 = Button(lower_frame, bg="white", text='Log', command=show_values)
    button2.place(relx=0.5, rely=0.16, anchor='n')

    button = Button(lower_frame, bg="white", text='Solve',
                    command=spherical_results_window)
    button.place(relx=0.5, rely=0.81, anchor='n')

    text2 = Label(swindow, bg="white",
                  text="Version 2.1, \u00A9 2022, Angus Siberry. All rights reserved.", font=(None, 8))
    text2.place(relx=0.68, rely=0.97, width=300, anchor='n')


HEIGHT = 875
WIDTH = 1450

master = Tk()
master.title("Alpha Dose Rate Calculator v2.1")

canvas0 = Canvas(master, height=HEIGHT, width=WIDTH, bg='black')
canvas0.pack()

background_image = PhotoImage(file='images/BCG.png')
background_image1 = PhotoImage(file='images/ps51.png')
background_image2 = PhotoImage(file='images/sp3.png')
size_image = PhotoImage(file='images/Particle_size1.png')


background_label = Label(master, image=background_image)
background_label.place(relwidth=1, relheight=1, relx=0.5, rely=0, anchor='n')

# bd = boarder bg = background and colout
p_frame = Frame(master, bg="white", bd=0)
p_frame.place(relx=0.1, rely=0.495, width=68, height=20, anchor='n')

button = Button(p_frame, bg="white", text='Click here', command=linear_window)
button.pack()

# bd = boarder bg = background and colout
s_frame = Frame(master, bg="white", bd=0)
s_frame.place(relx=0.86, rely=0.495, width=68, height=20, anchor='n')

s_button = Button(s_frame, bg="white", text='Click here',
                  command=spherical_window)
s_button.pack()

mainloop()
