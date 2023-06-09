# -*- coding: utf-8 -*-
# This is a script file for Bistatic Radar Doppler Shift Simulation
# Author: Wolfgang Kaufmann, 2019
# contact@ars-electromagnetica.de
# Algorithm of bistatic radar Doppler shift: William P. Thomson (1985):
# “Airborne Bistatic Radar Limitations and Sample Calculations”,
# Air Force Institute of Technology, Ohio.

import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import tkinter
import tkinter.filedialog, tkinter.messagebox

# ----------------------------------------------------------------------------
# Speed of light, km/s
C = 299792.458
# ----------------------------------------------------------------------------
# Receiver
# receiver altitude, km
RALT = 0
# ----------------------------------------------------------------------------
#Transmitter:
# transmitter-altitude, km
XALT = 0
# ----------------------------------------------------------------------------
# Simulation parameters:
tStep = 0.1    # Time intervall of Simulation in s
# ----------------------------------------------------------------------------

def calcVelocity(TEL, TAZ):
    # calculate target velocity vectors
    VTZ = VT * math.sin(TEL)
    VTY = VT * math.cos(TEL) * math.sin(TAZ)
    VTX = VT * math.cos(TEL) * math.cos(TAZ)
    return VTZ, VTY, VTX

def calcDS(X, Y, Z, VTZ,VTY, VTX):
    # vector-normalization
    A1 = math.sqrt(X**2 + Y**2 + (Z-RALT)**2)
    B1 = math.sqrt((SEP-X)**2 + Y**2 + (XALT-Z)**2)
    # doppler shift when arriving at the target, Hz
    FD1 = ((SEP-X)*VTX + Y*(-VTY) + (XALT-Z)*VTZ) / (B1*L)
    # new wavelength when leaving the target, km
    LN = C / (FO + FD1)
    # doppler shift when arriving at the receiver, Hz:
    FF = ((SEP-X)*VTX + Y*(-VTY) + (XALT-Z)*VTZ) / (B1*L)\
         - ((VTX*X) + (VTY*Y) + VTZ*(Z-RALT)) / (A1*LN)
    return FF

def calcMov(X, Y, Z,VTZ, VTY, VTX, taz, tel, dur):
    # Calculate movement of target in terms of x,y,z-coordinates
    # for flight time in steps of tStep
    moveD = []
    tNr = int(dur/tStep)
    for t in range(0,tNr):
        Zt = Z + VTZ * t * tStep
        Yt = Y + VTY * t * tStep
        Xt = X + VTX * t * tStep
        DS = calcDS(Xt, Yt, Zt, VTZ, VTY, VTX)
        moveD.append([t*tStep, DS, Xt, Yt, Zt, taz, tel])
    return moveD

def plotData(moveD, runs):
    try:
        # length of data file of one run
        dInterval = int(len(moveD)/runs)
        # height limit
        hLimit = int(enter_lh.get())
        # frequency shift limit
        fuLimit = int(enter_lfu.get())
        flLimit = int(enter_lfl.get())
    except:
        tkinter.messagebox.showerror\
        ("Input error",  "Check your Simulation Parameter (integers only)")
        return
    fig = plt.figure()
    # prepare data
    partSize = []
    tm = []
    dm = []
    xm = []
    ym = []
    zm = []
    i = 0
    j = 0
    for line in moveD:
        j = j + 1
        # limit height and frequency shift by storing only adequate
        # parts of the trajectories
        if (line[4] > hLimit) and (line[1] > flLimit) and (line[1] < fuLimit):
            tm.append(line[0])
            dm.append(round(line[1],2))
            xm.append(round(line[2],1))
            ym.append(round(line[3],1))
            zm.append(round(line[4],1))
            i = i + 1
        if j >= dInterval:
            # save length of filtered data chunk of the finished run
            partSize.append(i)
            # reset counters
            i = 0
            j = 0
    # plot frequency-data
    ax = fig.add_subplot(1, 2, 1)
    dLen = 0
    for k in range(runs):
        if partSize[k] > 0:
            # filter for zero crossing Doppler-Shifts
            if enter_cZ.get():
                if max(dm[dLen:(dLen+partSize[k])]) >= 0 and\
                       min(dm[dLen:(dLen+partSize[k])]) <= 0:
                    ax.plot(tm[dLen:(dLen+partSize[k])],\
                            dm[dLen:(dLen+partSize[k])])
            else:
                ax.plot(tm[dLen:(dLen+partSize[k])],dm[dLen:(dLen+partSize[k])])
        # sum lengths of data chunks
        dLen = dLen + partSize[k]
    ax.set_xlabel('Flight-Time [s]')
    ax.set_ylabel('Doppler Shift [Hz]')
    ax.grid(True)
    # plot move-Data
    ax = fig.add_subplot(1, 2, 2, projection='3d')
    dLen = 0
    for k in range(runs):
        if partSize[k] > 0:
            # filter for zero crossing Doppler-Shifts
            if enter_cZ.get():
                if max(dm[dLen:(dLen+partSize[k])]) >= 0 and\
                       min(dm[dLen:(dLen+partSize[k])]) <= 0:
                           ax.plot(xm[dLen:(dLen+partSize[k])],\
                                   ym[dLen:(dLen+partSize[k])],\
                                   zm[dLen:(dLen+partSize[k])])
            else:
                ax.plot(xm[dLen:(dLen+partSize[k])],\
                        ym[dLen:(dLen+partSize[k])],\
                        zm[dLen:(dLen+partSize[k])])
            # Markierung Startpunkte Target
            # filter for zero crossing Doppler-Shifts
            if enter_cZ.get():
                if max(dm[dLen:(dLen+partSize[k])]) >= 0 and\
                       min(dm[dLen:(dLen+partSize[k])]) <= 0:
                           ax.scatter(xm[dLen], ym[dLen], zm[dLen], marker="o")
            else:
                ax.scatter(xm[dLen], ym[dLen], zm[dLen], marker="o")
        # sum lengths of data chunks
        dLen = dLen + partSize[k]
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # plot baseline
    ax.plot([SEP,0,0], [0,0,0])
    ax.text(SEP, 0, 0, "T")
    ax.text(0, 0, 0, "R")
    fig.show()

def saveSim(moveD):
    filename=tkinter.filedialog.asksaveasfilename(defaultextension="csv")
    fName = str(filename)
    if len(fName) == 0:
        return
    try:
        d = open(fName, "w")
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "Saving file failed")
        return
    # Schreibe Header
    d.write("Time[s], Doppler[Hz], x[km], y[km], z[km], Bearing[°], Inclination[°]\r")
    # Schreibe Liste mit den Spalten:
    # Flight-Time, Doppler-Shift, x(t), y(t), z(t), bearing, inclination
    for line in moveD:
        d.write(str(line[0]) + "," + str(line[1]) + "," + str(line[2]) + ","\
                + str(line[3]) + "," + str(line[4]) + "," + str(line[5]) + ","\
                + str(line[6]) + "\r")
    # Schreibe Parameter
    d.write("SEP, F0, VT, , , , \r")
    d.write(str(SEP) + "," + str(FO) + "," + str(VT) + ",,,,\r")
    d.close()

def runSimulation():
    global SEP,FO,L,X,Y,Z,VT,TAz,TEl
    allMoves = []
    try:
        SEP = int(enter_TXP.get())
        FO = int(enter_TXF.get())
        L = C / FO
        X = int(enter_tx.get())
        Y = int(enter_ty.get())
        Z = float(enter_tz.get())
        #### to avoid division by zero set Z always > 0 ( => 1 cm) ####
        if Z == 0:
            Z = 0.00001
        VT = float(enter_tv.get())
        TAz = int(enter_ta.get())
        TEl = int(enter_te.get())
        tXNr = int(enter_xnr.get())
        tXDist = int(enter_xd.get())
        tYNr = int(enter_ynr.get())
        tYDist = int(enter_yd.get())
        tELNr = int(enter_enr.get())
        tELDist = int(enter_ed.get())
        tAZNr = int(enter_anr.get())
        tAZDist = int(enter_ad.get())
        tDur = int(enter_lt.get())
    except:
        tkinter.messagebox.showerror\
        ("Input error",  "Check your Transmitter and Target Entries (integers only)")
        return
    # tAZNr Startpunkte mit unterschiedlichem Azimuth
    for tazdist in range(TAz, TAz+tAZNr*tAZDist, tAZDist):
        # tELNr Startpunkte mit unterschiedlichen Steigungswinkeln
        for teldist in range(TEl, TEl+tELNr*tELDist, tELDist):
            VZ, VY, VX = calcVelocity(math.radians(teldist), math.radians(tazdist))
            # tYNr Startpunkte in y-Richtung im Abstand tYDist
            for y in range(Y, Y+tYNr*tYDist, tYDist):
                # tXNr Startpunkte in x-Richtung im Abstand tXDist
                for x in range(X, X+tXNr*tXDist, tXDist):
                    moveDS = calcMov(x, y, Z, VZ, VY, VX, tazdist, teldist, tDur)
                    allMoves = allMoves + moveDS
    plotData(allMoves,tXNr*tYNr*tELNr*tAZNr)
    if enter_s.get():
        saveSim(allMoves)

def quit():
    plt.close("all")
    main.destroy()

# Main -----------------------------------------------------------------------
main = tkinter.Tk()
main.title("MDopplerShift v0.32")
# Titel, Programm-Info
titel = tkinter.Label(main, text="Dopplershift Simulation of a Bistatic Radar")
titel["font"] = "Arial 14"
titel.grid(row=0, columnspan=2)

itxt = "This program uses the 3D cartesian coordinate system (x,y,z). It as-\r\n\
sumes the receiver is coincident with the coordinate origin (0,0,0).\r\n\
The transmitter is located on the x-axis in distance X to the recei-\r\n\
ver (X,0,0). The target can be anywhere in the space (x,y,z). Its\r\n\
bearing must be indicated as angle between its forward direction\r\n\
and the receiver>>transmitter direction. Its inclination must be\r\n\
specified towards the horizon."
info = tkinter.Text(main, width=70, height = 7)
info.insert("end", itxt)
info["borderwidth"] = 2
info["relief"] = "ridge"
info.grid(row=1, columnspan=2, padx=10, pady=10)

# Transmitter
label_TXP = tkinter.Label(main, text="Transmitter Distance [km]:")
label_TXP.grid(row=2, column=0, sticky=tkinter.W, padx=10)
preset_TXP = tkinter.StringVar()
preset_TXP.set("800")
enter_TXP = tkinter.Entry(main, textvariable=preset_TXP)
enter_TXP.grid(row=2, column=1)

label_TXF = tkinter.Label(main, text="Transmitter Frequency [Hz]:")
label_TXF.grid(row=3, column=0, sticky=tkinter.W, padx=10)
preset_TXF = tkinter.StringVar()
preset_TXF.set("143050000")
enter_TXF = tkinter.Entry(main, textvariable=preset_TXF)
enter_TXF.grid(row=3, column=1)

# Target:
label_tx = tkinter.Label(main, text="x - Start Position of the Target [km]:")
label_tx.grid(row=4, column=0, sticky = tkinter.W, padx = 10)
preset_tx = tkinter.StringVar()
preset_tx.set("400")
enter_tx = tkinter.Entry(main, textvariable=preset_tx)
enter_tx.grid(row=4, column=1)

label_ty = tkinter.Label(main, text="y - Start Position of the Target, \
- = right side [km]:")
label_ty.grid(row=5, column=0, sticky = tkinter.W, padx = 10)
preset_ty = tkinter.StringVar()
preset_ty.set("0")
enter_ty = tkinter.Entry(main, textvariable=preset_ty)
enter_ty.grid(row=5, column=1)

label_tz = tkinter.Label(main, text="z - Start Position of the Target [km]:")
label_tz.grid(row=6, column=0, sticky = tkinter.W, padx = 10)
preset_tz = tkinter.StringVar()
preset_tz.set("100")
enter_tz = tkinter.Entry(main, textvariable=preset_tz)
enter_tz.grid(row=6, column=1)

label_tv = tkinter.Label(main, text="Velocity of the Target [km/s]:")
label_tv.grid(row=7, column=0, sticky = tkinter.W, padx = 10)
preset_tv = tkinter.StringVar()
preset_tv.set("10")
enter_tv = tkinter.Entry(main, textvariable=preset_tv)
enter_tv.grid(row=7, column=1)

label_ta = tkinter.Label(main, text="Bearing of the Target, counterclockwise [°]:")
label_ta.grid(row=8, column=0, sticky = tkinter.W, padx = 10)
preset_ta = tkinter.StringVar()
preset_ta.set("0")
enter_ta = tkinter.Entry(main, textvariable=preset_ta)
enter_ta.grid(row=8, column=1)

label_te = tkinter.Label(main, text="Angle of Inclination of the Target, \
- = falling [°]:")
label_te.grid(row=9, column=0, sticky = tkinter.W, padx = 10)
preset_te = tkinter.StringVar()
preset_te.set("0")
enter_te = tkinter.Entry(main, textvariable=preset_te)
enter_te.grid(row=9, column=1)

label_lt = tkinter.Label(main, text="Flight Time of the Target [s]:")
label_lt.grid(row=10, column=0, sticky = tkinter.W, padx = 10)
preset_lt = tkinter.StringVar()
preset_lt.set("10")
enter_lt = tkinter.Entry(main, textvariable=preset_lt)
enter_lt.grid(row=10, column=1)

# Multiple Targets:

label_mt = tkinter.Label(main, text="------------------ Multiple Targets \
---------------------")
label_mt.grid(row=11, column=0, sticky = tkinter.W, padx = 10)

label_xnr = tkinter.Label(main, text="Number of Targets in x-Direction:")
label_xnr.grid(row=12, column=0, sticky = tkinter.W, padx = 10)
preset_xnr = tkinter.StringVar()
preset_xnr.set("1")
enter_xnr = tkinter.Entry(main, textvariable=preset_xnr)
enter_xnr.grid(row=12, column=1)

label_xd = tkinter.Label(main, text="Distance of the Targets in x-Direction [km]:")
label_xd.grid(row=13, column=0, sticky = tkinter.W, padx = 10)
preset_xd = tkinter.StringVar()
preset_xd.set("100")
enter_xd = tkinter.Entry(main, textvariable=preset_xd)
enter_xd.grid(row=13, column=1)

label_ynr = tkinter.Label(main, text="Number of Targets in y-Direction:")
label_ynr.grid(row=14, column=0, sticky = tkinter.W, padx = 10)
preset_ynr = tkinter.StringVar()
preset_ynr.set("1")
enter_ynr = tkinter.Entry(main, textvariable=preset_ynr)
enter_ynr.grid(row=14, column=1)

label_yd = tkinter.Label(main, text="Distance of the Targets in y-Direction [km]:")
label_yd.grid(row=15, column=0, sticky = tkinter.W, padx = 10)
preset_yd = tkinter.StringVar()
preset_yd.set("50")
enter_yd = tkinter.Entry(main, textvariable=preset_yd)
enter_yd.grid(row=15, column=1)

label_enr = tkinter.Label(main, text="Number of Targets with increasing Inclinations per Loc.:")
label_enr.grid(row=16, column=0, sticky = tkinter.W, padx = 10)
preset_enr = tkinter.StringVar()
preset_enr.set("1")
enter_enr = tkinter.Entry(main, textvariable=preset_enr)
enter_enr.grid(row=16, column=1)

label_ed = tkinter.Label(main, text="Increase of the Angle of Inclination of the \
Targets [°]:")
label_ed.grid(row=17, column=0, sticky = tkinter.W, padx = 10)
preset_ed = tkinter.StringVar()
preset_ed.set("30")
enter_ed = tkinter.Entry(main, textvariable=preset_ed)
enter_ed.grid(row=17, column=1)

label_anr = tkinter.Label(main, text="Number of Targets with increasing Bearing per Loc.:")
label_anr.grid(row=18, column=0, sticky = tkinter.W, padx = 10)
preset_anr = tkinter.StringVar()
preset_anr.set("1")
enter_anr = tkinter.Entry(main, textvariable=preset_anr)
enter_anr.grid(row=18, column=1)

label_ad = tkinter.Label(main, text="Increase of the Angle of Bearing of the \
Targets [°]:")
label_ad.grid(row=19, column=0, sticky = tkinter.W, padx = 10)
preset_ad = tkinter.StringVar()
preset_ad.set("30")
enter_ad = tkinter.Entry(main, textvariable=preset_ad)
enter_ad.grid(row=19, column=1)

label_sp = tkinter.Label(main, text="---------------- Simulation Parameter \
-----------------")
label_sp.grid(row=20, column=0, sticky = tkinter.W, padx = 10)

label_lh = tkinter.Label(main, text="Lower Height Limit [km]:")
label_lh.grid(row=21, column=0, sticky = tkinter.W, padx = 10)
preset_lh = tkinter.StringVar()
preset_lh.set("0")
enter_lh = tkinter.Entry(main, textvariable=preset_lh)
enter_lh.grid(row=21, column=1)

label_lfu = tkinter.Label(main, text="Upper Doppler-Shift Limit [Hz]:")
label_lfu.grid(row=22, column=0, sticky = tkinter.W, padx = 10)
preset_lfu = tkinter.StringVar()
preset_lfu.set("100000")
enter_lfu = tkinter.Entry(main, textvariable=preset_lfu)
enter_lfu.grid(row=22, column=1)

label_lfl = tkinter.Label(main, text="Lower Doppler-Shift Limit [Hz]:")
label_lfl.grid(row=23, column=0, sticky = tkinter.W, padx = 10)
preset_lfl = tkinter.StringVar()
preset_lfl.set("-100000")
enter_lfl = tkinter.Entry(main, textvariable=preset_lfl)
enter_lfl.grid(row=23, column=1)

enter_cZ = tkinter.BooleanVar()
label_cz = tkinter.Checkbutton(main, text="Only Zero Crossing Doppler-Shifts",
                               variable=enter_cZ, onvalue=True, offvalue=False)
label_cz.grid(row=24, column=0, sticky = tkinter.W, padx = 10)

enter_s = tkinter.BooleanVar()
label_s = tkinter.Checkbutton(main, text="Save Simulation", variable=enter_s,
                              onvalue=True, offvalue=False)
label_s.grid(row=25, column=0, sticky = tkinter.W, padx = 10)

# Buttons
bRun = tkinter.Button(main, text="Run Simulation", command=runSimulation)
bRun.grid(row=26, column=0, sticky = tkinter.W, padx=10, pady=10)

bQuit = tkinter.Button(main, text="Quit", command=quit)
bQuit.grid(row=26, column=1, pady=10)
main.mainloop()