# -*- coding: utf-8 -*-
# Process Data-Files from Meteor Logger
# copyright 2017-2020, Wolfgang Kaufmann
# under the terms of the GNU Affero General Public License
# as published by the Free Software Foundation, version 3 of the License
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.


import tkinter
import tkinter.filedialog, tkinter.messagebox, tkinter.scrolledtext
import sys
from matplotlib import pyplot as plt
from datetime import datetime
from datetime import timedelta
import matplotlib.dates as mdates
import scipy as sc
from scipy import stats
from scipy import interpolate
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy import ndimage
import numpy as np
import ephem
import math
import os

DataList = []
fDist = 48000 / 2048  # FFT-Bin-Breite MeteorData
spM_F = 2*np.pi/24    # Frequenz Tagesgang sporadische Meteore


def convertDtoHMS(degree):
    frac = degree/360
    hours = frac * 24
    minutes = (hours - int(hours)) * 60
    seconds = (minutes - int(minutes)) * 60
    tString = str(int(hours))+":"+ str(int(minutes))+":"+ str(int(seconds))
    return tString

def convertDtoDMS(degree):
    minutes = (degree - int(degree)) * 60
    seconds = (minutes - int(minutes)) * 60
    tString = str(int(degree))+":"+ str(int(minutes))+":"+ str(int(seconds))
    return tString

def convertDMStoD(DMS):
    dms = str(DMS)
    index1 = dms.find(":")
    index2 = dms.find(":", index1+1)
    degree = int(dms[:index1])
    minutes = int(dms[index1+1:index2])
    seconds = float(dms[index2+1:])
    deg = degree + minutes/60 + seconds/3600
    return deg

def askDaylightSaving():
    answer = tkinter.messagebox.askyesno("Process Time",
                            "Correct for Daylight Saving (Local Time + 1 h)?")
    if answer == 1:
        return True
    else:
        return False

def displayFrame():
    # Zeitrahmen der Datei anzeigen (Start-/Enddatum und -zeit)
    endD = len(DataList)-1
    fformat = "{:.8}"
    fstring = DataList[0][1] + "  " + fformat.format(DataList[0][2]) + "  ->  "
    fstring = fstring + DataList[endD][1] + "  " + fformat.format(DataList[endD][2])
    fstring = fstring + "  UTC"
    label4["text"] = fstring
    return fstring

def readConfig():
    # Konfigurationsdatei auslesen:
    # zeilenweise Amplitude, Offset, Zeitverschiebung, geo. Länge, geo. Breite,
    # Höhe ü. Meeresspiegel, Mindestdauer eines overdense Meteors, Max. Zeitspanne
    # zwischen zwei Signalen für ihre Zusammenführung
    param = []
    try:
        cf = open("config.txt", "r")
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "Opening config-file failed")
        sys.exit()
    try:
        configData = cf.read()
    except:
        cf.close()
        tkinter.messagebox.showerror\
        ("I/O error",  "Reading config-file failed")
        sys.exit()
    cf.close()
    # Aus ConfigData die einzelnen Parameter extrahieren
    try:
        DataLines = configData.split(chr(10))
        for line in DataLines:
            if line.find("#") == -1:
                param.append(line)
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "config-file corrupted")
        sys.exit()
    # Test, ob alle 9 Parameter gefunden wurden
    if len(param) < 9:
        tkinter.messagebox.showerror\
        ("I/O error",  "config-file corrupted")
        sys.exit()
    try:
        return float(param[0]), float(param[1]), int(param[2]),\
                convertDtoDMS(float(param[3])), convertDtoDMS(float(param[4])),\
                float(param[5]), float(param[6]), float(param[7]), float(param[8])
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "config-file contains invalid numbers")
        sys.exit()

def colButton(c):
    bIDel["bg"] = c
    bConf["bg"] = c
    bThre["bg"] = c
    bSNR["bg"] = c

def renumber(DataList):
    # Meteor-Nummern aufsteigend?
    neuNum = False
    for line in range(len(DataList)-1):
        if DataList[line][0] > DataList[line+1][0]:
            neuNum = True
            break
    if neuNum:
        # Neu-Nummerierung erforderlich
        sigNumber = 1
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        index = 0
        while index < len(DataList):
            meteorNr = DataList[index][0]
            istart = index
            while DataList[index][0] == meteorNr:
                index = index + 1
                if index >= len(DataList):
                    break
            istop = index
            # Neu-Nummerierung
            for i in range(istart, istop):
                DataList[i][0] = sigNumber
            sigNumber = sigNumber + 1
    else:
        return

def openfb():
    global DataList, SNR, spM_A, spM_O, dUTC, GeoLon, dTime, fFlat
    global GeoLat, GeoElevation, tOverdense, tSpan, thresholdAlt
    filename=tkinter.filedialog.askopenfilename()
    fName = str(filename)
    if fName == "":
        return
    else:
        label1["text"]=fName
    # Datei öffnen
    try:
        d = open(fName, "r")
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "Opening file failed")
        return
    # Datei komplett einlesen
    try:
        Dataset = d.read()
    except:
        d.close()
        tkinter.messagebox.showerror\
        ("I/O error",  "Reading file failed")
        return
    d.close()
    # Alte Daten aus DataList löschen
    del DataList[:]
    # GUI zurücksetzen
    dataTable.delete("0.0","end")
    carrierF.delete("0", "end")
    threshold.delete("0", "end")
    colButton("SystemButtonFace")
    # Busy-Cursor
    main.config(cursor="watch")
    main.update()
    # Neue Daten einlesen
    try:
        # Zerlege in Zeilen
        DataLines = Dataset.split(chr(10))
        # Erzeuge 2D-Liste mit den Spalten:
        # Nr, Date, Time, f [Hz], Signal, Noise\r
        # 0   1     2     3       4       5
        for line in DataLines:
            dTemp = line.split(",")
            # Beschriftungszeilen, Leerzeilen und f = -1.0 - Zeilen entf.
            if (dTemp[0] != "Nr") and (dTemp[0] != "") and (dTemp[3] != "-1.0"):
                # Zeiteintrag auf fehlende Nachkommastellen prüfen
                if dTemp[2].find(".") < 0:
                    dTemp[2] = dTemp[2] + ".000"
                DataList.append([int(dTemp[0]), dTemp[1], dTemp[2],
                                float(dTemp[3]), float(dTemp[4]), float(dTemp[5])])
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "This is not a MeteorData-File")
        return
    # Meteor-Nummern auf Einmaligkeit prüfen, gfs. neu nummerieren
    renumber(DataList)
    # Konfigurationsparameter laden
    spM_A,spM_O,dUTC,GeoLon,GeoLat,GeoElevation,tOverdense,tSpan,thresholdAlt=readConfig()
    # Falls Roh-Datensatz, Zeitkorrektur vornehmen
    if not dImport.get():
        # UTC- und Daylight Saving Zeitkorrektur
        DS = askDaylightSaving()
        for line in DataList:
            timestr = line[1] + " " + line[2]
            # Zeitstring in Zeitobjekt umwandeln
            t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
            # in UTC umrechnen
            if dUTC < 0:
                t_object = t_object - timedelta(hours = abs(dUTC))
            else:
                t_object = t_object + timedelta(hours = abs(dUTC))
            if DS:
                # Daylight Saving berücksichtigen
                t_object = t_object - timedelta(hours = 1)
            line[1] = str(t_object)[:10]
            line[2] = str(t_object)[11:]
            if line[2].find(".") < 0:
                    line[2] = line[2] + ".000"
    #Start-/Endezeit der MeteorData-Messungen ausgeben
    displayFrame()
    # Verlauf Trägerfrequenz schätzen:
    dTime, fFlat = fReference()
    # Vorbereitungen:
    SNR = False
    colButton("lightgrey")
    main.config(cursor="")

def savefb():
    if len(DataList) == 0:
        return
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
    # Schreibe Liste mit den Spalten:
    # Nr, Date, Time, f [Hz], Signal/SNR, Noise
    # 0   1     2     3       4           5
    line = "Nr,Date,UTC,f [Hz],"
    if SNR:
        line = line + "SNR,Noise\r"
    else:
        line = line + "Signal,Noise\r"
    d.write(line)
    for i in DataList:
        line = str(i[0]) + "," + i[1] + "," + i[2] + "," + str(i[3])
        line = line + "," + str(i[4]) + "," + str(i[5]) + "\r"
        d.write(line)
    d.close()

def calcDistSlope(n, x):
    # Berechne die Steigung der Verteilung
    # 1. Logarithmiere die Klassenhäufigkeit
    y = []
    for i in n:
        y.append(np.lib.scimath.log10(i))
    # 2. Führe lieare Regression sukzessiv über je 75 Werte durch
    step = 75
    sl = []  # slope
    cl = []  # class
    for i in range(0, (len(x)-step), step):
        slope, intercept, r_value, p_value, std_err = \
                        sc.stats.linregress(x[i:i+step],y[i:i+step])
        sl.append(slope)
        index = int(i+round(step/2,0))
        cl.append(x[index])
    # 3. Index des negativsten Steigungswertes ermitteln
    isort = np.argsort(sl)
    # 4. nur die Steigungen bis zum negativsten Steigungswertes zurückgeben
    return cl[:isort[0]+1], sl[:isort[0]+1]

def my_kde_bandwidth(obj, fac=0.45):
    # We use Scott's Rule, multiplied by a constant factor.
    return np.power(obj.n, -1./(obj.d+4)) * fac

def plotData():
    if len(DataList) == 0:
        return
    # Extrahiere Date-, Time-, f- und Power-Spalte alternierend
    # nach Meteor-Nr in je zwei Variable
    f1 = []
    f2 = []
    p1 = []
    p2 = []
    t1 = []
    t2 = []
    index = 0
    while index < len(DataList):
        meteorI = DataList[index][0]
        while meteorI == DataList[index][0]:
            f1.append(DataList[index][3])
            p1.append(DataList[index][4])
            timestr = DataList[index][1] + " " + DataList[index][2]
            # Zeitstring in Zeitobjekt umwandeln
            t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
            t1.append(t_object)
            index = index + 1
            if index >= len(DataList):
                break
        if index >= len(DataList):
                break
        meteorII = DataList[index][0]
        while meteorII == DataList[index][0]:
            f2.append(DataList[index][3])
            p2.append(DataList[index][4])
            timestr = DataList[index][1] + " " + DataList[index][2]
            # Zeitstring in Zeitobjekt umwandeln
            t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
            t2.append(t_object)
            index = index + 1
            if index >= len(DataList):
                break
    # Grafikausgabe
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S.%f'))
    fig.autofmt_xdate(rotation=45)
    plt.title("Detected Signals")
    plt.ylabel("Frequency [Hz]")
    ax1.plot(t1, f1, ".", markersize=5, color="black")
    ax1.plot(t2, f2, ".", markersize=5, color="red")
    plt.grid(visible=True, which="major", axis="y", linestyle="dotted")

    ax2 = fig.add_subplot(2, 1, 2, sharex = ax1)
    ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M:%S.%f'))
    fig.autofmt_xdate(rotation=45)
    if SNR:
        plt.ylabel("SNR")
    else:
        plt.ylabel("Power [arbitrary units]")
    plt.yscale('log')
    plt.xlabel("UTC")
    ax2.plot(t1, p1, ".", markersize=5, color="black")
    ax2.plot(t2, p2, ".", markersize=5, color="red")
    plt.grid(visible=True, which="major", axis="y", linestyle="dotted")
    fig.show()

def plotDist(aMax, aSum, tDur, tDist):
    plt.rcParams.update({'font.size': 12})
    # Zeichne Verteilung der log. kumul. Amplitude
    fig = plt.figure()
    ax1 = fig.add_subplot(1,2,1)
    plt.title("Logarithmic Cumulative Amplitude Distribution")
    if SNR:
        plt.xlabel("log SNR")
    else:
        plt.xlabel("log Amplitude")
    # die aMax-Verteilung hat einen glatteren Verlauf, die Auflösung gegenüber
    # aSum ist im overdense-Bereich jedoch schlechter
    n1, m1, p1 = ax1.hist(aMax, bins=1000, cumulative=-1, histtype="step", log=True)
    ax1.set_ylabel("Number", color='b')
    ax1.tick_params('y', colors='b')
    plt.grid(visible=True, which="major", axis="x", linestyle="dotted")
    # Zeichne die Steigung der Verteilung als Punkte
    c1, s1 = calcDistSlope(n1, m1)
    ax3 = ax1.twinx()
    ax3.plot(c1, s1, "^", markersize=6, color="r")
    if SNR:
        ax3.set_ylabel("Slope [log Number/log SNR]", color="r")
    else:
        ax3.set_ylabel("Slope [log Number/log Amplitude]", color="r")
    ax3.tick_params('y', colors='r')

    # Zeichne Verteilung der Signalabstände
    # tDist logarithmieren für Detail-Darstellung der kleineren Werte
    ax2 = fig.add_subplot(1,2,2)
    plt.title("Distribution of the Time Lag between the Signals")
    plt.xlabel("log Time Lag [s]")
    # Ermitteln der Anzahl Klassen
    if len(tDist) < 2:
        binNr = 4
    else:
        binNr = 4 * int(np.lib.scimath.log10(len(tDist)) / np.lib.scimath.log10(2))
    # Histogramm zeichnen
    n2, m2, p2 = ax2.hist(tDist, bins=binNr, log=False)
    ax2.set_ylabel("Number", color='b')
    ax2.tick_params('y', colors='b')
    plt.grid(visible=True, which="major", axis="x", linestyle="dotted")
    fig.show()
    return n1, m1, n2, m2

def plotCounts(indication, cH, yFit, alt, tTitle, yErr):
    plt.rcParams.update({'font.size': 12})
    # zwei Variablen mit x-Beschriftung und y-Werten
    x = []
    y = []
    firstItem = True
    for line in cH:
        if (line[1] == "00") or firstItem:
            x.append(line[0] + "  " + line[1])
            firstItem = False
        else:
            # Nur jeden geraden x-Wert ausgeben
            if int(line[1])%2==0:
                x.append(line[1])
            else:
                x.append(" ")
        y.append(line[2])
    xpos = np.arange(len(cH))
    # Bar-Graph
    fig = plt.figure()
    plt.title(tTitle)
    plt.ylabel("Counts/h")
    plt.bar(xpos, y, yerr = yErr)
    plt.xticks(xpos, x, rotation=90)
    plt.grid(visible=True, which="major", axis="y", linestyle="dotted")
    plt.xlabel("UTC [h]")
    # Zeige angepasste Sinusfunktion als Trendkurve
    if spM_C.get() > 0:
        plt.plot(xpos, yFit, "k--")
    # Zeige Verlauf der Höhe ü.d. Horizont des Meteorstrom-Radianten
    if indication:
        plt.twinx()
        plt.plot(xpos, alt, "-.", color="r")
        plt.ylabel("Altitude [°]", color="r")
        plt.tick_params('y', colors='r')
        plt.ylim([0, 90])
    fig.show()

def plotHeadEchoes(slopeList):
    plt.rcParams.update({'font.size': 12})
    # Graphische Darstellung auf 0-12 kHz beschränken, da durch
    # Zeit-Jitter rechnerisch Steigungen bis 16 kHz entstehen

    # Slope in Array kopieren
    y = []
    for line in slopeList:
        y.append(line[4])
    # Zeitachse in Dezimal-Tagen erstellen
    t = []
    timestr = slopeList[0][1] + " 00:00:00.000" # Startzeit: 00.00 h erster Tag
    t_object1 = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
    for line in slopeList:
        timestr = line[1] + " " + line[2]
        # Zeitstring in Zeitobjekt umwandeln
        t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
        # Zeitdifferenz zum Startzeitpunkt festellen
        deltaT = t_object - t_object1
        tseconds = deltaT.total_seconds()
        t.append(tseconds/86400)

    # Min/Max-Werte
    ymin = min(y)
    ymax = max(y)
    tmin = min(t)
    tmax = max(t)

    # Peform the kernel density estimate
    xx, yy = np.mgrid[tmin:tmax:200j, ymin:ymax:200j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([t, y])
    kernel = sc.stats.gaussian_kde(values, bw_method=my_kde_bandwidth)
    f = np.reshape(kernel(positions).T, xx.shape)
    # Kernel density map
    fig = plt.figure()
    ax = fig.gca()
    ax.set_xlim(tmin, tmax)
    if ymin < -12000:
        ymin = -12000
    if ymax > 0:
        ymax = 0
    ax.set_ylim(ymin, ymax)
    # Contourf plot
    ax.contourf(xx, yy, f, cmap='Blues')
    # Contour plot: cset = ax.contour(xx, yy, f, colors='k')
    # Scatter plot
    ax.scatter(t, y, s=10, marker=".")
    # Label plot: ax.clabel(cset, inline=1, fontsize=6)
    ax.set_xlabel('Decimal Day')
    ax.set_ylabel('Frequency Gradient [Hz/s]')
    plt.title("Frequency Slopes of Head Echoes")
    fig.show()

def plotDur():
    plt.rcParams.update({'font.size': 12})
    if len(DataList) == 0:
        return
    # Hole Liste mit Nummern der Signale mit exp. Decay (überwiegend underdense M.)
    nrList = compileDecay()
    if len(nrList) < 10:
        tkinter.messagebox.showwarning("Underdense Trails",
                                       "Insufficient number of underdense trails")
        return
    # Berechne Dauer dieser Signale
    t = []
    for nr in nrList:
        di = 0
        for (i,line) in enumerate(DataList):
            if line[0] == nr:
                istop = i
                di = di+1
        istart = istop - (di - 1)
        # 1. Zeit-Strings in Zeitvariablen umwandeln
        timestr1 = DataList[istart][1] + " " + DataList[istart][2]
        t_object1 = datetime.strptime(timestr1, '%Y-%m-%d %H:%M:%S.%f')
        timestr2 = DataList[istop][1] + " " + DataList[istop][2]
        t_object2 = datetime.strptime(timestr2, '%Y-%m-%d %H:%M:%S.%f')
        # 2. Zeitdifferenz festellen
        sigDur = t_object2 - t_object1
        tseconds = sigDur.total_seconds()
        t.append(tseconds)
    # Histogramm anzeigen
    fig = plt.figure()
    plt.hist(t, bins=100, normed=1, cumulative=True, histtype="step")
    plt.grid(visible=True, which="major", axis="x", linestyle="dotted")
    plt.grid(visible=True, which="major", axis="y", linestyle="dotted")
    plt.title("Generic Cumulative Distribution")
    plt.xlabel("Duration of Underdense Meteors [s]")
    plt.ylabel("Relative Number")
    fig.show()

def getFC():
    fc_ = carrierF.get()
    try:
        fc = float(fc_)
        if (fc < 0) or (fc > 20000):
            # Falls fNotch außerhalb Analysegrenzen mit -1 besetzten
            fc = -1.0
            carrierF.delete("0", "end")
        else:
            answer = tkinter.messagebox.askokcancel("f-Reference",
                            "You are going to override f-reference!")
            if not answer:
                fc = -1.0
                carrierF.delete("0", "end")
    except:
        fc = -1.0
        carrierF.delete("0", "end")
    return fc

def fReference():
    # TX Trägerfrequenz schätzen (Statistische Kennzahl: Median)
    f = [row[3] for row in DataList]
    TXF = np.median(f)
    # nur den Bereich +- 3.05*fDist um die geschätzte Trägerfrequenz (+- 71 Hz)
    # berücksichtigen
    dTime = []
    fF = []
    fc1 = TXF + fDist * 3.05
    fc2 = TXF - fDist * 3.05
    for zeile in DataList:
        if (zeile[3] < fc1) and (zeile[3] > fc2):
            fF.append(zeile[3])
            # Datum und Uhrzeit in einen Zeitstring umwandeln
            timestr = zeile[1] + " " + zeile[2]
            # Zeitstring in Zeitobjekt umwandeln
            t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
            # Sekunden seit 1.1.1970
            seconds = t_object.timestamp()
            dTime.append(seconds)
    # # Savitzky-Golay-Filter anwenden, um den Frequenzverlauf zu glätten
    intervall = 5001
    if len(fF) < intervall:
        modus = "constant"
    else:
        modus = "interp"
    fFl = savgol_filter(fF, window_length=intervall, polyorder=1, mode=modus, cval=TXF)
    fFlat = fFl.tolist()
    # Anzeige Referenzfrequenz
    dataTable.insert("end", "\r\n" + "Estimated f-Reference [Hz] per 3h, Min = " +\
                     str(round(min(fFlat),1)) + ", Max = " +\
                     str(round(max(fFlat),1)) + "\r\n-- ")
    # Median-Funktion definieren
    medianF = sc.interpolate.interp1d(dTime, fFlat, fill_value='extrapolate')
    # f-Referenz im 3h-Takt ausgeben (3*60*60=10800s)
    for i in range (int(dTime[0]), int(dTime[-1]), 10800):
        dataTable.insert("end", "{:.1f}".format(medianF(i)) + " -- ")
    # Rückgabe zwei Vektoren mit den Zeit- und Frequenz-Werten für eine
    # Interpolation.
    return dTime, fFlat

def reduceInterference():
    global DataList, dTime, fFlat
    if len(DataList) == 0:
        return
    bIDel["bg"] = "lightgreen"
    # Prüfen, ob Trägerfrequenz eingegeben wurde
    fnew = getFC()
    if fnew > -1:
        # Zeit- und Frequenzvektor mit Anfangs- und Endwert besetzen
        t1 = dTime[0]
        t2 = dTime[-1]
        del dTime[:]
        dTime.append(t1)
        dTime.append(t2)
        del fFlat[:]
        fFlat.append(fnew)
        fFlat.append(fnew)
    # Median-Funktion definieren
    medianF = sc.interpolate.interp1d(dTime, fFlat, fill_value='extrapolate')
    # Anzeige eliminierter Signale
    dataTable.insert("end", "\r\nThe following signals are removed:\r\n")
    dataTable.insert("end", "Number:\r\n")
    # DataList Signal für Signal durchlaufen
    index = 0
    number = 0
    delCounts = []
    while index < len(DataList):
        flag = False
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        istart = index
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList):
                break
        istop = index
        # Bestimme zugehörigen fMedian
        timestr = DataList[istart][1] + " " + DataList[istart][2]
        # Zeitstring in Zeitobjekt umwandeln
        t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
        # Sekunden seit 1.1.1970
        seconds = t_object.timestamp()
        ynew = medianF(seconds)
        fMedian1 = ynew + 2.05 * fDist
        fMedian2 = ynew - 2.05 * fDist
        # Signalreihe muss mind. ein f im Bereich fMedian +- 2*fDist haben
        for i in range(istart,istop):
            if (DataList[i][3] < fMedian1) and (DataList[i][3] > fMedian2):
                flag = True
        if not flag:
            # zu eliminierende Counts anzeigen
            dataTable.insert("end", "{:>8}".format(str(DataList[istart][0])))
            number = number + 1
            if number > 8:
                number = 0
                dataTable.insert("end", "\r\n")
            # Counts notieren für späteres Löschen
            delCounts.append([istart, istop])
    # Notierte Counts löschen
    for line in reversed(delCounts):
        del DataList[line[0]:line[1]]

def conflateSignals():
    global DataList
    if len(DataList) == 0:
        return
    # Hole Liste mit Meteornummern, die ausgeschlossen werden sollen
    # 1. Signale, die exponentiell abfallen
    nrListD = compileDecay()
    # 2. Head-Echoes
    hL, nrListH = compileHeadEcho()
    # Hole Zeitspanne
    tN = timedelta(seconds=tSpan)
    dataTable.insert("end", "\r\n" + "Signals conflated!\r\n")
    dataTable.insert("end", "The following signals are conflated:\r\n")
    dataTable.insert("end", "1st,2nd Number:\r\n")
    bConf["bg"] = "lightgreen"
    # Gehe DataList schrittweise durch
    number = 0
    index = 1
    while index < len(DataList):
        metNr1 = DataList[index-1][0]
        # Wechsel der Meteor-Nr feststellen
        while DataList[index][0] == metNr1:
            index = index+1
            if index >= len(DataList):
                break
        if index >= len(DataList):
                break
        metNr2 = DataList[index][0]
        # Zeit-Strings in Zeitvariablen umwandeln
        timestr1 = DataList[index-1][1] + " " + DataList[index-1][2]
        t_object1 = datetime.strptime(timestr1, '%Y-%m-%d %H:%M:%S.%f')
        timestr2 = DataList[index][1] + " " + DataList[index][2]
        t_object2 = datetime.strptime(timestr2, '%Y-%m-%d %H:%M:%S.%f')
        # Zeitdifferenz festellen
        tD = t_object2 - t_object1
        # Zeitdifferenz mit max. Zeitspanne vergleichen
        if (tD < tN) and (metNr1 not in nrListD) and (metNr2 not in nrListH):
            # Zusammengeführte Counts anzeigen
            dataTable.insert("end", "{:>17}".format(str(metNr1) \
                             + "," + str(metNr2)))
            number = number + 1
            if number > 3:
                number = 0
                dataTable.insert("end", "\r\n")
            # Count-Nrn zusammenführen
            while DataList[index][0] == metNr2:
                DataList[index][0] = metNr1
                index = index+1
                if index >= len(DataList):
                    break
            # metNr2 geht in metNr1 auf, daher für den nächsten Durchlauf:
            if metNr2 in nrListD:
                nrListD.append(metNr1)
        else:
            index = index + 1

def applyThreshold():
    global DataList
    if len(DataList) == 0:
        return
    # Hole Maximale Zeitspanne und prüfe auf valide Zahleneingabe
    ths = threshold.get()
    try:
        t = float(ths)
        if t < 0:
            t = 0
    except:
        tkinter.messagebox.showerror("Input error",  "This is not a valid figure")
        threshold.delete(0,"end")
        return
    dataTable.insert("end", "\r\n" +  "Threshold applied!\r\n")
    dataTable.insert("end", "The following signals are removed:\r\n")
    dataTable.insert("end", "Number:\r\n")
    bThre["bg"] = "lightgreen"
    index = 0
    number = 0
    while index < len(DataList):
        flag = False
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        istart = index
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList):
                break
        istop = index
        # Signalreihe muss mind. 1 Wert >= Threshold haben
        for i in range(istart,istop):
            if DataList[i][4] >= t:
                flag = True
        if not flag:
            # zu eliminierende Counts anzeigen
            dataTable.insert("end", "{:>8}".format(str(DataList[istart][0])))
            number = number + 1
            if number > 8:
                number = 0
                dataTable.insert("end", "\r\n")
            # Count eliminieren
            del DataList[istart:istop]
            index = istart

def averageN():
    global DataList
    if len(DataList) == 0:
        return
    dataTable.insert("end", "\r\n" + "Noise averaged\r\n")
    index = 0
    while index < len(DataList):
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        istart = index
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList):
                break
        istop = index
        # Bilde Mittelwert des Meteor-Counts
        noiseSum = 0
        for i in range(istart, istop):
            noiseSum = noiseSum + DataList[i][5]
        # berechne Noise-Mittelwert des Meteor-Counts
        noiseMean = noiseSum / (istop-istart)
        # Schreibe Noise-Mittelwert
        for i in range(istart, istop):
            DataList[i][5] = noiseMean

def calcSNR():
    global DataList, SNR
    if (len(DataList) == 0) or SNR:
        return
    bSNR["bg"] = "lightgreen"
    SNR = True
    # Noise innerhalb eines jeden Counts mitteln
    averageN()
    # Ersetzte Power gegen SNR
    for line in DataList:
        line[4] = line[4] / line[5]
    dataTable.insert("end", "SNR calculated\r\n")
    dataTable.insert("end", "All following operations are done now on base of SNR!")

def assignCount():
    # Erzeuge Liste mit den Angaben lokales Datum \ Stunde \ Counts/h
    count = 1
    countH = []
    metNr = DataList[0][0]
    dateHour = DataList[0][1] + " " + DataList[0][2][0:2]
    for i in range(len(DataList)):
        # Meteor-Nr gewechselt?
        if DataList[i][0] != metNr:
            # noch die gleiche Stunde?
            if (DataList[i][1] + " " + DataList[i][2][0:2]) == dateHour:
                count = count+1
            else:
                # Stundenwechsel: Counts in CountH schreiben
                countH.append([dateHour[:10], dateHour[11:], count])
                count = 1
                # nächste Mess-Stunde in dateHour schreiben
                dateHour = DataList[i][1] + " " + DataList[i][2][0:2]
                # Stunden ohne counts berücksichtigen
                timestr = DataList[i-1][1] + " " + DataList[i-1][2][0:2]
                t_object1 = datetime.strptime(timestr, '%Y-%m-%d %H')
                timestr = DataList[i][1] + " " + DataList[i][2][0:2]
                t_object2 = datetime.strptime(timestr, '%Y-%m-%d %H')
                # Differenz zur nächsten Mess-Stunde in Anzahl Stunden
                deltaH = t_object2 - t_object1
                tseconds = deltaH.total_seconds()
                dH = int(tseconds/3600)
                if dH < 1:
                    dH = 1
                # Stundenlücken mit 0 counts füllen
                for h in range(1,dH,1):
                    t = t_object1 + timedelta(hours = h)
                    d = t.strftime('%Y-%m-%d')
                    h = t.strftime('%H')
                    countH.append([d, h, 0])
            metNr = DataList[i][0]
    # auch letzte unvollendete Stunde ausgeben
    countH.append([dateHour[:10], dateHour[11:], count])
    return countH

def applySmoother(hC):
    if len(hC) < 3:
        return hC
    # Counts extrahieren
    hc = [line[2] for line in hC]
    if Smoo.get() == 1:
        # Gauß-Filter
        hcs = sc.ndimage.gaussian_filter1d(hc, 1, order=0, mode="nearest")
    else:
        # Gleitender Mittelwert, mode = nearest
        a = hc[0]
        b = hc[len(hc)-1]
        hc.insert(0, a)
        hc.append(b)
        hcs = []
        for i in range(1, len(hc)-1, 1):
           hcs.append((hc[i-1]+hc[i]+hc[i+1])//3)
    # Ergebnis in hC zurückschreiben
    for (i,line) in enumerate(hC):
        line[2] = hcs[i]
    return hC

def sinusFunc(x, freq, amplitude, phase, offset):
    return np.sin(x * freq + phase) * amplitude + offset

def fitSinus(selection, cH):
    # Stunden-Array mit Länge cH erstellen
    t = np.arange(0, len(cH), 1)
    # Counts in einfache Liste schreiben
    y = []
    for line in cH:
        y.append(line[2])
    # Startwerte für die Iteration festlegen
    guessFreq = 2*np.pi/24
    guessAmpl = 3*np.std(y)/(2**0.5)
    guessPhase = 0.0
    guessOffset = min(y)
    p0=[guessFreq, guessAmpl, guessPhase, guessOffset]

    # Wahl: Sinus-fit
    if selection == 2:
        # Führe Näherung durch
        try:
            fit = curve_fit(sinusFunc, t, y, p0 = p0)
        except:
            yFit = np.empty(len(cH))
            yFit.fill(0)
            return yFit, False
        dataTable.insert("end", "\r\n" + "Sinus-Fit:\r\n")
        dataTable.insert("end", "amplitude = " + str(fit[0][1]) + \
                         "  offset = " + str(fit[0][3]) + "\r\n")
        dataTable.insert("end", "frequency = " + str(fit[0][0]) + \
                         "  phase = " + str(fit[0][2]))
        return sinusFunc(t, *fit[0]), True

    # Wahl: Parametrisierte Sinus-Fktn oder off
    fit = []
    # 1. Bestimmung der Phase: Bezugspunkt 0h
    for i in range(len(cH)):
        if cH[i][1] == "00":
            break
    # 2. Falls Bezugspunkt gefunden, setze Parameter:
    if i < (len(cH)-1):
        fit.append(spM_F)
        fit.append(spM_A)
        i = i + dUTC
        fit.append(-i * (2*np.pi/24))
        fit.append(spM_O)
        return sinusFunc(t, *fit), True
    else:
        yFit = np.empty(len(cH))
        yFit.fill(0)
        return yFit, False

def writeCount(indication, countH, spM, alt):
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
    # Schreibe Liste mit den Spalten:
    # UTC-Date, UTC-Hour, Counts/h,
    # 0        1          2
    # sporadic Meteor Counts/h or fitted sinus-funct., altitude
    # 3                                                4
    line = "Date(UTC),Hour(UTC),Counts/h"
    if spM_C.get() == 2:
        line = line + ",fitted sinus"
    if spM_C.get() == 1:
        line = line + ",sporadicMC/h"
    if indication:
        line = line + ",Altitude[°]\r"
    else:
        line = line + "\r"
    d.write(line)
    for i in range(len(countH)):
        line = countH[i][0] + "," + countH[i][1] + "," + str(countH[i][2])
        if spM_C.get() > 0:
            line = line + "," + str(round(spM[i]))
        if indication:
            line = line + "," + str(round(alt[i], 1)) + "\r"
        else:
            line = line + "\r"
        d.write(line)
    d.close()

def runCount():
    t = []
    alt = []
    az = []
    result = False
    if len(DataList) == 0:
        return
    # Stelle Counts per h zusammen:
    hCount = assignCount()
    # Glätte die Counts per h:
    if Smoo.get() > 0:
        hsCount = applySmoother(hCount)
        countH = hsCount
    else:
        countH = hCount
    # Berechne Höhe und Azimuth des Meteorstrom-Radianten
    if radiant.get():
        result, ra, dec = getRadiantPos()
        if result:
           for line in countH:
               # countH[0]: YYYY-mm-dd, countH[1]: %H
               t.append(line[0].replace("-", "/") + " " + line[1] + ":30:00")
           alt, az = calcRadiantPos(ra, dec, t)
    # Berechne Verlauf der sporadischen Meteore / passe eine Sinuskurve in den
    # Verlauf der count rate
    yFit, success = fitSinus(spM_C.get(), countH)
    # Für Grafikanzeige ein 0-Fehler-Array erstellen
    cErr = np.empty(len(countH))
    cErr.fill(0)
    # count rate anzeigen
    plotCounts(result, countH, yFit, alt, "Hourly Count Rates", cErr)
    # count rate in file schreiben
    writeCount(result, countH, yFit, alt)

def writeZHR(countH, zErr, alt):
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
    # Schreibe Liste mit den Spalten:
    # UTC-Date, UTC-Hour, Counts/h, Fehler, Höhe
    # 0        1          2         3       4
    line = "Date(UTC),Hour(UTC),ZHR_r,Error,Altitude[°]\r"
    d.write(line)
    for i in range(len(countH)):
        line = countH[i][0] + "," + countH[i][1] + "," + str(countH[i][2])
        line = line + "," + str(round(zErr[i]))
        line = line + "," + str(round(alt[i], 1)) + "\r"
        d.write(line)
    d.close()

def runRZHR():
    if len(DataList) == 0:
        return
    result, ra, dec = getRadiantPos()
    if not result:
        return
    # Stelle Counts per h zusammen:
    hCount = assignCount()
    # Glätte die Counts per h:
    if Smoo.get() > 0:
        hsCount = applySmoother(hCount)
        countH = hsCount
    else:
        countH = hCount
    # Bilde Zeitvektor
    t = []
    for line in countH:
        # countH[0]: YYYY-mm-dd, countH[1]: %H
        t.append(line[0].replace("-", "/") + " " + line[1] + ":30:00")
    # Berechne den zeitl. Verlauf von Azimuth und Höhe des Meteorstromes
    alt, az = calcRadiantPos(ra, dec, t)
    # Berechne Verlauf der sporadischen Meteore
    yFit, success = fitSinus(1, countH)
    if not success:
        tkinter.messagebox.showwarning("r_ZHR", "Radio-ZHR could not be calculated")
        return
    # Subtrahiere sporadische Meteorzahl von Gesamtzahl und führe Zenith-
    # korrektur durch
    Z_Err = []
    for i in range(len(countH)):
        if alt[i] > thresholdAlt:
            zCorr = math.sin(math.radians(alt[i]))
            cDiff = countH[i][2] - round(yFit[i])
            if cDiff < 0:
                cDiff = 0
            ZHR_R = round(cDiff / zCorr)
            ZHR_R_Err = ZHR_R / math.sqrt(1 + cDiff)
            countH[i][2] = ZHR_R
        else:
            countH[i][2] = 0
            ZHR_R_Err = 0
        Z_Err.append(ZHR_R_Err)
    plotCounts(radiant.get(), countH, yFit, alt, "Radio-ZHR", Z_Err)
    writeZHR(countH, Z_Err, alt)

def compileData():
    ShortList = []
    index = 0
    t = [" "]
    # Meteor-Schauer: Prüfe, ob Anzeige und Radiant-Daten eingegeben
    result = False
    if radiant.get():
        result, ra, dec = getRadiantPos()
    while index < len(DataList):
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        istart = index
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList):
                break
        istop = index
        # Finde Datenzeile mit dem höchsten Signalwert (sigIndex) und
        # summiere alle Signal-/Noisewerte eines Counts (sigSum/noiseSum)
        noiseSum =0
        sigSum = 0
        sig = 0
        sigIndex = istart
        for i in range(istart, istop):
            sigSum = sigSum + DataList[i][4]
            noiseSum = noiseSum + DataList[i][5]
            if DataList[i][4] > sig:
                sig = DataList[i][4]
                sigIndex = i
        # Dauer des Meteor-Counts berechnen (sigDur)
        # 1. Zeit-Strings in Zeitvariablen umwandeln
        timestr1 = DataList[istart][1] + " " + DataList[istart][2]
        t_object1 = datetime.strptime(timestr1, '%Y-%m-%d %H:%M:%S.%f')
        timestr2 = DataList[istop-1][1] + " " + DataList[istop-1][2]
        t_object2 = datetime.strptime(timestr2, '%Y-%m-%d %H:%M:%S.%f')
        # 2. Zeitdifferenz festellen
        sigDur = t_object2 - t_object1
        tseconds = sigDur.total_seconds()
        # 3. Berechne Höhe und Azimuth des Meteorstrom-Radianten
        if result:
            # DataList[1]: YYYY-mm-dd, DataList[2]: 00:00:00.000
            t[0] = DataList[sigIndex][1].replace("-", "/") + " " + DataList[sigIndex][2]
            alt, az = calcRadiantPos(ra, dec, t)
        # Resultate in ShortList schreiben
        if result:
            ShortList.append([DataList[sigIndex],sigSum,noiseSum,tseconds,alt[0]])
        else:
            ShortList.append([DataList[sigIndex],sigSum,noiseSum,tseconds])
    return ShortList, result

def summariseData():
    if len(DataList) == 0:
        return
    ShortList, result = compileData()
    tDist = tDistance()
    # ShortList und tDist in File schreiben
    # 1. Dateiname und Speicherpfad holen
    filename=tkinter.filedialog.asksaveasfilename(defaultextension="csv")
    fName = str(filename)
    if len(fName) == 0:
        return
    try:
        d = open(fName, "w")
    except:
        answer=tkinter.messagebox.askretrycancel\
        ("I/O error",  "Saving file failed")
        if answer == 1:
            return
        else:
            sys.exit()
    # 2. Schreibe Liste mit den Spalten:
    line = "Nr,Date,UTC,f [Hz],"
    if SNR:
        line = line + "Max-SNR,SNR-Sum,Duration [s], t-Distance [s]"
    else:
        line = line + "Max-Signal,Noise,Signal-Sum,Noise-Sum,Duration [s],t-Distance [s]"
    if result:
        line = line + ",Altitude [°]\r"
    else:
        line = line + "\r"
    d.write(line)
    # 3. Schreibe Daten
    for (j, dline) in enumerate(ShortList):
        # Nr, Datum, Uhrzeit, f, Max-Sig/SNR
        line = str(dline[0][0]) + "," + dline[0][1] + "," + dline[0][2] + ","
        line = line + str(dline[0][3]) + "," + str(dline[0][4]) + ","
        # Noise
        if not SNR:
            line = line + str(dline[0][5]) + ","
        # Sig-/SNR-Sum
        line = line + str(dline[1]) + ","
        # Noise-Sum
        if not SNR:
            line = line + str(dline[2]) + ","
        # Duration
        line = line + str(dline[3]) + ","
        # t-Distance
        if j < len(tDist):
            line = line + str(tDist[j])
        # Altitude
        if result:
            line = line + "," + str(dline[4]) + "\r"
        else:
            line = line + "\r"
        d.write(line)
    d.close()

def compileHeadEcho():
    fHigh = 5 * fDist     # Mindeststartfrequenz für ein Headecho
    fLow = 2 * fDist      # Untere Frequenzgrenze für ein Headecho
    # Prüfen, ob Trägerfrequenz eingegeben wurde
    fnew = getFC()
    if fnew > -1:
        # Zeit- und Frequenzvektor mit Anfangs- und Endwert besetzen
        t1 = dTime[0]
        t2 = dTime[-1]
        del dTime[:]
        dTime.append(t1)
        dTime.append(t2)
        del fFlat[:]
        fFlat.append(fnew)
        fFlat.append(fnew)
    # Median-Funktion definieren
    medianF = sc.interpolate.interp1d(dTime, fFlat, fill_value='extrapolate')
    # Finde Headechos
    HeadList = []
    NrHeadList = []  # Ajout de l'initialisation de la variable NrHeadList
    index = 0
    while index < len(DataList):
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        istart = index
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList):
                break
        istop = index-1
        # Jeden Meteor-Count prüfen, ob er ein Headecho enthält
        # 1. Bestimme zugehörigen fMedian
        timestr = DataList[istart][1] + " " + DataList[istart][2]
        # Zeitstring in Zeitobjekt umwandeln
        t_object = datetime.strptime(timestr, '%Y-%m-%d %H:%M:%S.%f')
        # Sekunden seit 1.1.1970
        seconds = t_object.timestamp()
        fMedian = medianF(seconds)
        # 2. f(istart) > fHigh-Schwelle
        if DataList[istart][3] >= (fHigh + fMedian):
            ihead = istart
            D0 = False
            fDesc = False
            mem = 0
            # 3. f absteigend bis fLow-Schwelle
            while ihead < istop:
                if DataList[ihead][3] < (fLow + fMedian):
                    D0 = True
                    break
                fDesc = DataList[ihead][3] >= DataList[ihead+1][3]
                if not fDesc:
                    break
                # Max. zwei gleiche Frequenzwerte nacheinander akzeptieren:
                if DataList[ihead][3] == DataList[ihead+1][3]:
                    if mem == 1:
                        break
                    else:
                        mem = 1
                else:
                    mem = 0
                ihead = ihead + 1
            # 4. head echo muss mindestens drei Messwerte umfassen
            fLength = ihead >= (istart + 3)
            # head echo-Daten in Liste schreiben
            if fDesc and fLength and D0:
                for i in range(istart,ihead):
                    HeadList.append(DataList[i])
            # zusätzlich eine Liste nur mit den head echo-Nummern erzeugen
            heNr = [row[0] for row in HeadList]
            NrHeadList = list(set(heNr))
            NrHeadList.sort()
    return HeadList, NrHeadList

def writeHeadEcho(NrHeadList):
    # 1. Dateiname und Speicherpfad holen
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
    # 2. Schreibe Liste mit den Spalten:
    line = "Nr,Date,UTC,f [Hz],"
    if SNR:
        line = line + "SNR\r"
    else:
        line = line + "Signal,Noise\r"
    d.write(line)
    # 3. Nur die head echo-Datensätze aus DataList schreiben
    index = 0
    while index < len(DataList):
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        if meteorNr in NrHeadList:
            istart = index
            while DataList[index][0] == meteorNr:
                index = index + 1
                if index >= len(DataList):
                    break
            istop = index
            # Daten schreiben
            for i in range(istart, istop):
                line = str(DataList[i][0]) + "," + DataList[i][1] + ","\
                + DataList[i][2] + "," + str(DataList[i][3]) + ","\
                + str(DataList[i][4]) + "," + str(DataList[i][5]) + "\r"
                d.write(line)
        else:
            index = index + 1
    d.close()

def exportHeadEchoes():
    if len(DataList) == 0:
        return
    # Liste mit Headechos erstellen
    HeadList, NrHeadList = compileHeadEcho()
    if len(HeadList) > 1:
        # Liste der Headechos auf Disk schreiben
        writeHeadEcho(NrHeadList)
    else:
        tkinter.messagebox.showwarning("Export Head Echoes",
                                       "No Head Echoes found")

def compileDecay():
    # Nach exponentiell fallendem Signalverlauf suchen:
    index = 0
    nrList = []
    y = []
    while index < len(DataList):
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        istart = index
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList):
                break
        istop = index
        # Finde Datenzeile mit dem höchsten Signalwert (sigIndex)
        sig = 0
        sigIndex = istart
        for i in range(istart, istop):
            if DataList[i][4] > sig:
                sig = DataList[i][4]
                sigIndex = i
        # Kurvenverlauf ab Max-Sig prüfen: x, y-Vektoren erstellen
        x = np.arange(1, istop-sigIndex+1, 1)
        # y linearisieren durch Logarithmieren
        del y[:]
        for j in range(sigIndex, istop, 1):
            y.append(np.lib.scimath.log(DataList[j][4]))
        # lineare Regression
        if len(y) > 2:
            slope, intercept, r_value, p_value, std_err = \
                                sc.stats.linregress(x, y)
            if abs(r_value) > 0.975:
                nrList.append(meteorNr)
    return nrList

def exportOverdense():
    if len(DataList) == 0:
        return
    # Struktur ShortList: DataList[signalpeak-index],sigSum,noiseSum,tseconds,alt[0]
    ShortList, result = compileData()
    # Liste mit Nummern der overdense-trails erstellen
    nrList = []
    for line in ShortList:
        if line[3] >= tOverdense:
            nrList.append(line[0][0])
    # Overdense trails exportieren
    writeTrails(nrList)

def exportUnderdense():
    if len(DataList) == 0:
        return
    # Struktur ShortList: DataList[signalpeak-index],sigSum,noiseSum,tseconds,alt[0]
    ShortList, result = compileData()
    # Liste mit Nummern der underdense-trails erstellen
    nrList = []
    for line in ShortList:
        if line[3] < tOverdense:
            nrList.append(line[0][0])
    # Underdense trails exportieren
    writeTrails(nrList)

def writeTrails(nrList):
     # 1. Dateiname und Speicherpfad holen
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
    # 2. Schreibe Liste mit den Spalten:
    line = "Nr,Date,UTC,f [Hz],"
    if SNR:
        line = line + "SNR\r"
    else:
        line = line + "Signal,Noise\r"
    d.write(line)
    # Schreibe Meteordaten entsprechend nrList:
    for dline in DataList:
        if dline[0] in nrList:
            line = str(dline[0]) + "," + dline[1] + "," + dline[2] + ","
            line = line + str(dline[3]) + "," + str(dline[4])
            if SNR:
                line = line + "\r"
            else:
                line = line + "," + str(dline[5]) + "\r"
            d.write(line)
    d.close()

def calcSlope(HeadList):
    slopeList = []
    f = []
    t = []
    index = 0
    while index < len(HeadList):
        # Ermittle die Indexpositionen der jeweiligen Headecho's
        meteorNr = HeadList[index][0]
        istart = index
        while HeadList[index][0] == meteorNr:
            index = index + 1
            if index >= len(HeadList):
                break
        istop = index
        # Bilde zwei Arrays für f und t
        for i in range(istart,istop):
            f.append(HeadList[i][3])
            # Zeit-Strings in Zeitvariablen umwandeln
            dtimes = HeadList[i][2].split(":")
            tseconds = int(dtimes[0])*3600 + int(dtimes[1])*60 + float(dtimes[2])
            t.append(tseconds)
        slope, intercept, r_value, p_value, std_err = sc.stats.linregress(t,f)
        deltaF = f[0] - f[len(f)-1]
        del t[:]
        del f[:]
        # Berechne die Signalsumme des Headechos
        hSum = 0
        for i in range(istart,istop):
            hSum = hSum + HeadList[i][4]
        # Schreibe Daten in Liste
        # Nr, Datum, Uhrzeit, Signalsumme/SNR-Summe, Slope, deltaF, Messpunkte
        slopeList.append([HeadList[istart][0], HeadList[istart][1],
                         HeadList[istart][2], hSum, slope, deltaF, (istop-istart)])
    return slopeList

def writeSlopeList(slopeList, result, alt, az):
    # 1. Dateiname und Speicherpfad holen
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
    # 2. Schreibe Liste mit den Spalten:
    # Nr, Datum, Uhrzeit, deltaF, Dauer, Signalsumme/SNR-Summe, Slope, Alt, Az
    line = "Nr,Date,UTC,deltaF [Hz],MPoints,"
    if SNR:
        line = line + "SNR-Sum,Slope [Hz/s]"
    else:
        line = line + "Signal-Sum,Slope [Hz/s]"
    if result:
        line = line + ",Altitude [°],Azimuth [°]\r"
    else:
        line = line + "\r"
    d.write(line)
    for i in range(len(slopeList)):
        line = str(slopeList[i][0]) + "," + slopeList[i][1] + "," + slopeList[i][2]
        line = line + "," + str(slopeList[i][5]) + "," + str(slopeList[i][6])+ ","
        line = line + str(slopeList[i][3]) + "," + str(slopeList[i][4])
        if result:
            line = line + "," + str(alt[i]) + "," + str(az[i]) + "\r"
        else:
            line = line + "\r"
        d.write(line)
    d.close()

def analyseHeadEchoes():
    if len(DataList) == 0:
        return
    # Liste mit Headechos erstellen
    HeadList, NrHeadList = compileHeadEcho()
    # Steigung der Headechos berechnen
    slopeList = calcSlope(HeadList)
    # Verteilung der Steigungen darstellen
    if len(slopeList) > 0:
        plotHeadEchoes(slopeList)
        # Berechne Höhe und Azimuth des Meteorstrom-Radianten
        t = []
        alt = []
        az = []
        result = False
        if radiant.get():
            result, ra, dec = getRadiantPos()
            if result:
               for line in slopeList:
                   # slopeList[1]: YYYY-mm-dd, slopeList[2]: 00:00:00.000
                   t.append(line[1].replace("-", "/") + " " + line[2])
               alt, az = calcRadiantPos(ra, dec, t)
        # slopeList auf Disc schreiben
        writeSlopeList(slopeList, result, alt, az)
    else:
        tkinter.messagebox.showwarning("Head Echoes", "No Head Echoes found")

def tDistance():
    # Ermittle zeitl. Abstände zwischen den Signalen
    tDistList = []
    index = 0
    while index < len(DataList)-1:
        # Ermittle die Indexpositionen der einzelnen Meteor-Counts
        meteorNr = DataList[index][0]
        while DataList[index][0] == meteorNr:
            index = index + 1
            if index >= len(DataList)-1:
                break
        istop = index
        # zeitl. Diff. für istop (= Start neues Sig.) und istop-1 (= Ende Vorgänger)
        # berechnen
        timestr1 = DataList[istop-1][1] + " " + DataList[istop-1][2]
        t_object1 = datetime.strptime(timestr1, '%Y-%m-%d %H:%M:%S.%f')
        timestr2 = DataList[istop][1] + " " + DataList[istop][2]
        t_object2 = datetime.strptime(timestr2, '%Y-%m-%d %H:%M:%S.%f')
        deltaH = t_object2 - t_object1
        tseconds = deltaH.total_seconds()
        tDistList.append(tseconds)
    # letzten Wert verwerfen
    return tDistList[:len(tDistList)-1]

def compileDist():
    if len(DataList) == 0:
        return
    # Struktur ShortList: DataList[signalpeak-index],sigSum,noiseSum,tseconds,alt[0]
    ShortList, result = compileData()
    # Extrahiere Power-Spalte und wandle Power in Amplitude um
    aMax = []
    for line in ShortList:
        aMax.append(np.lib.scimath.log10(np.lib.scimath.sqrt(line[0][4])))
    # Extrahiere Powersummen-Spalte und wandle in Power in Amplitude um
    aSum = []
    for line in ShortList:
        aSum.append(np.lib.scimath.log10(np.lib.scimath.sqrt(line[1])))
    # Extrahiere Duration-Spalte
    tDur = []
    for line in ShortList:
        tDur.append(line[3])
    # Ermittle zeitl. Abstände zwischen den Signalen
    tDist = tDistance()
    # Logarithmieren für eine bessere Aulösung der kleinen Werte
    for i in range(len(tDist)):
        tDist[i] = np.lib.scimath.log10(tDist[i])
    # Grafikausgabe
    n1, m1, n2, m2 = plotDist(aMax, aSum, tDur, tDist)
    # Kumulierte Verteilungen speichern
    writeDist(n1, m1, n2, m2)

def writeDist(n1, m1, n2, m2):
    # 1. Dateiname und Speicherpfad holen
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
    # 2. Schreibe Datumsrahmen
    d.write(displayFrame() + ",\r")
    # 3. Schreibe Liste mit den Spalten:
    # Klassenmitte, Klassenhäufigkeit
    if SNR:
        d.write("Log Cumulative Amplitude/Noise Distribution\r")
        d.write("Class [log SNR]" + "," + "Class Frequency\r")
    else:
        d.write("Log Cumulative Amplitude Distribution\r")
        d.write("Class [log Amplitude]" + "," + "Class Frequency\r")
    for i in range(len(n1)):
        classM = m1[i] + ((m1[i+1] - m1[i]) / 2)
        d.write(str(classM) + "," + str(np.lib.scimath.log10(n1[i])) + "\r")
    d.write("Distribution of Time Laps between the Signals\r")
    d.write("Class [log s]" + "," + "Class Frequency\r")
    for i in range(len(n2)):
        classM = m2[i] + ((m2[i+1] - m2[i]) / 2)
        d.write(str(classM) + "," + str(n2[i]) + "\r")
    d.close()

def getRadiantPos():
    # Hole RA, DEC und prüfe auf valide Zahleneingabe
    r = rascension.get()
    try:
        RA = float(r)
        if RA < 0:
            RA = 0
        if RA > 360:
            RA = 360
    except:
        tkinter.messagebox.showwarning("RA", "RA: This is not a valid figure")
        #r.delete(0,"end")
        return False, 0, 0
    d = declination.get()
    try:
        DEC = float(d)
        if DEC < -90:
            DEC = -90
        if DEC > 90:
            DEC = 90
    except:
        tkinter.messagebox.showwarning("DEC", "DEC: This is not a valid figure")
        #d.delete(0,"end")
        return False, 0, 0
    return True, RA, DEC

def calcRadiantPos(r, d, t):
    alt = []
    az = []
    # Define a new fixed body:
    RA = convertDtoHMS(r)   # Dezimalgrad in Stunde:Minute:Sekunde
    DEC = convertDtoDMS(d)  # Dezimalgrad in Grad:Minute:Sekunde
    line = "Radiant,f," + RA + "," + DEC + ",5,2015"
    NFB = ephem.readdb(line)
    # Observers position
    GRAVES = ephem.Observer()
    GRAVES.lon = GeoLon
    GRAVES.lat = GeoLat
    GRAVES.elevation = GeoElevation
    # Berechne alt and az für den übergebenen timestring
    # Date/Time im Format "2017/08/13 11:25:00"
    for i in t:
        GRAVES.date = i
        NFB.compute(GRAVES)
        alt.append(convertDMStoD(NFB.alt))  # Höhe (D:M:S)
        az.append(convertDMStoD(NFB.az))    # Azimuth (D:M:S)
    return alt, az

def exportRMOB():
    if len(DataList) == 0:
        return
    # Hourly count rates holen: YYYY-MM-DD / Stunde / Counts
    countH = assignCount()
    # Datei öffnen
    filename=tkinter.filedialog.asksaveasfilename(defaultextension="txt")
    fName = str(filename)
    if len(fName) == 0:
        return
    try:
        d = open(fName, "w")
    except:
        tkinter.messagebox.showerror\
        ("I/O error",  "Saving file failed")
        return
    # Header schreiben
    line = " |0|1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|20|21|22|23|\r\n"
    d.write(line)
    # Alle Tage von 1 bis 31 und Stunden von 0 bis 23 durchgehen
    index = 0
    for dd in range(1, 32):
        line = str(dd) + "|"
        for hh in range(24):
            tTemp = countH[index][0].split("-")
            # Tag prüfen
            if int(tTemp[2]) == dd:
                # Stunde prüfen
                if int(countH[index][1]) == hh:
                    # 0 counts mit nicht gemessen gleichsetzen: -1
                    if countH[index][2] == 0:
                        line = line + "-1|"
                    else:
                        line = line + str(countH[index][2]) + "|"
                    if hh == 23:
                        line = line + "\r\n"
                    if index < (len(countH)-1):
                        index = index + 1
                else:
                    line = line + "-1|"
                    if hh == 23:
                        line = line + "\r\n"
            else:
                line = line + "-1|"
                if hh == 23:
                    line = line + "\r\n"
        d.write(line)
    d.close()

def openPdf():
    # pdf Datei zur Anzeige öffnen
    if os.path.isfile("Manual_ProcessData.pdf"):
        os.system("Manual_ProcessData.pdf")
    else:
        tkinter.messagebox.showerror\
        ("I/O error",  "Manual not found")

def printAbout():
    global aboutWindow
    aboutWindow = tkinter.Toplevel(master=main)
    aboutWindow.title("About")
    t = tkinter.Text(aboutWindow, width=65, height=23, padx=5, pady=5)
    t.pack()
    aboutText = "Process MeteorData v1.67\r\n\n\
    This software processes and analyses the raw data\r\n\
    of Meteor Logger v1.2 and newer versions.\r\n\n\
    copyright 2017-2020 Wolfgang Kaufmann\r\n\
    contact@ars-electromagnetica.de\r\n\n\
    This program is free software: you can redistribute it\r\n\
    and/or modify it under the terms of the GNU Affero General\r\n\
    Public License as published by the Free Software Foundation,\r\n\
    version 3 of the License.\r\n\n\
    This program is distributed in the hope that it will be\r\n\
    useful, but WITHOUT ANY WARRANTY; without even the\r\n\
    implied warranty of MERCHANTABILITY or FITNESS FOR A \r\n\
    PARTICULAR PURPOSE. See the GNU Affero General Public\r\n\
    License for more details.\r\n\n\
    You should have received a copy of the GNU Affero General\r\n\
    Public License along with this program.  If not, see\r\n\
    <http://www.gnu.org/licenses/>.\r\n"
    t.insert("end", aboutText)
    t.tag_add("HIGH", "1.0", "1.24")
    t.tag_add("MARK", "6.0", "7.40")
    t.tag_config("HIGH", font=('Arial', 12, 'bold', 'italic'))
    t.tag_config("MARK", foreground = "blue")

def beenden():
    plt.close()
    try:
        aboutWindow.destroy()
    except:
        pass
    main.destroy()

# Hauptfenster erzeugen
main=tkinter.Tk()
main.title("ProcessData v1.67")
main.geometry("650x780")

# Überschrift
label0 = tkinter.Label(main, text="Process Data from \'Meteor Logger\'")
label0["font"] = "Arial 16"
label0.grid(row=0, column=0, sticky=tkinter.W, padx=10, pady=10)

# Menue-Eintrag erzeugen
mBar = tkinter.Menu(main)

mPath = tkinter.Menu(mBar)
mPath["tearoff"] = 0
mPath.add_command(label = "Open", command = openfb)
mPath.add_command(label = "Save as", command = savefb)
mPath.add_separator()
mPath.add_command(label = "Quit", command = beenden)
mBar.add_cascade(label = "File", menu = mPath)

mFunc = tkinter.Menu(mBar)
mFunc["tearoff"] = 0
mFunc.add_command(label = "Summarise Data", command = summariseData)
mFunc.add_command(label = "Export Zero Cross HE", command = exportHeadEchoes)
mFunc.add_command(label = "Export Overdense Trails", command = exportOverdense)
mFunc.add_command(label = "Export Underdense Trails", command = exportUnderdense)
mFunc.add_command(label = "Export RMOB-file", command = exportRMOB)
mFunc.add_command(label = "Hourly Count Rates", command = runCount)
mFunc.add_command(label = "Radio-ZHR", command = runRZHR)
mFunc.add_command(label = "Distributions", command = compileDist)
mFunc.add_command(label = "Analyse Zero Cross HE", command = analyseHeadEchoes)
mBar.add_cascade(label = "Functions", menu = mFunc)

mConf = tkinter.Menu(mBar)
mConf["tearoff"] = 0

mSpor = tkinter.Menu(mConf)
mSpor["tearoff"] = 0
spM_C = tkinter.IntVar()
spM_C.set(0)
mSpor.add_radiobutton(label = "Sporadic Meteor Count Rate", variable = spM_C,
                      value = 1)
mSpor.add_radiobutton(label = "Sinus-fit of Count Rate", variable = spM_C, value = 2)
mSpor.add_radiobutton(label = "Off", variable = spM_C, value = 0)
mConf.add_cascade(label = "Show Trend", menu = mSpor)

radiant = tkinter.BooleanVar()
radiant.set(False)
mConf.add_checkbutton(label = "Show Radiant", variable = radiant,
                       onvalue = True, offvalue = False)

mSmooth = tkinter.Menu(mConf)
mSmooth["tearoff"] = 0
Smoo = tkinter.IntVar()
Smoo.set(0)
mSmooth.add_radiobutton(label = "Gaussian Filter", variable = Smoo, value = 1)
mSmooth.add_radiobutton(label = "Running Mean", variable = Smoo, value = 2)
mSmooth.add_radiobutton(label = "Off", variable = Smoo, value = 0)
mConf.add_cascade(label = "Apply a Smoother", menu = mSmooth)

dImport = tkinter.BooleanVar()
dImport.set(False)
mConf.add_checkbutton(label = "Load Processed Data", variable = dImport,
                       onvalue = True, offvalue = False)

mBar.add_cascade(label = "Configuration", menu = mConf)

mInf = tkinter.Menu(mBar)
mInf["tearoff"] = 0
mInf.add_command(label = "Approximate Limit", command = plotDur)
mInf.add_command(label = "About", command = printAbout)
mInf.add_command(label = "Manual", command = openPdf)
mBar.add_cascade(label = "Info", menu = mInf)

main["menu"] = mBar

# Filenamen anzeigen
label1=tkinter.Label(main)
label1.grid(row=1, column=0, sticky=tkinter.W, padx=10, pady=5)

# Zeitrahmen der Messung anzeigen
label4=tkinter.Label(main)
label4.grid(row=2, column=0, sticky=tkinter.W, padx=10, pady=5)

# Labelframe einfügen
lF1 = tkinter.LabelFrame(main, relief="groove", text="Cleaning Up Records")
lF1.grid(row=3, column=0, sticky=tkinter.W, padx=10, pady=5)

# Button für das Reduzieren von Interferenzen
bIDel = tkinter.Button(lF1, text = "Reduce Interference",
                       command = reduceInterference)
bIDel.grid(row=0, column=0, sticky=tkinter.W, padx=10, pady=5)

# Eingabe der Audiofrequenz des Trägersignals
label2 = tkinter.Label(lF1, text = "f-Reference [Hz]:")
label2["font"] = "Arial 8 italic"
label2["fg"] = "gray"
label2.grid(row=0, column=0, sticky=tkinter.W, padx=150, pady=5)
carrierF = tkinter.Entry(lF1)
carrierF["bg"] = "#EEEEEE"
carrierF.grid(row=0, column=0, sticky=tkinter.W, padx=250, pady=5)

# Button für das Zusammenführen von unterbrochenen Signalreihen
bConf = tkinter.Button(lF1, text = "Conflate", command = conflateSignals)
bConf.grid(row=1, column=0, sticky=tkinter.W, padx=10, pady=5)

# Labelframe einfügen
lF2 = tkinter.LabelFrame(main, relief="groove", text="Optional Signal Processing")
lF2.grid(row=4, column=0, sticky=tkinter.W, padx=10, pady=20)

# Button für Threshold
bThre = tkinter.Button(lF2, text = "Apply Threshold", command = applyThreshold)
bThre.grid(row=0, column=0, sticky=tkinter.W, padx=10, pady=5)
# Threshold eingeben
label3 = tkinter.Label(lF2, text = "Threshold:")
label3.grid(row=0, column=0, sticky=tkinter.W, padx=180)
threshold = tkinter.Entry(lF2)
threshold.grid(row=0, column=0, sticky=tkinter.W, padx=250)

# Button für die Errechnung des SNR
bSNR = tkinter.Button(lF2, text = "Convert to SNR", command = calcSNR)
bSNR.grid(row=1, column=0, sticky=tkinter.W, padx=10, pady=8)

# Anzeigetafel für vorgenommenes Processing
dataTable = tkinter.scrolledtext.ScrolledText(main, width=75, height=23)
dataTable.grid(row=5, column=0, sticky=tkinter.W, padx=10, pady=5)

# Buttons für Grafikdarstellung
bPlot = tkinter.Button(main, text = "Plot Signal", command = plotData)
bPlot.grid(row=6, column=0, sticky=tkinter.W, padx=10, pady=10)

# Eingabe Radiant Meteorschauer
# Rektascension eingeben:
labelR = tkinter.Label(main, text="Meteorshower:   Right ascension [°]:")
labelR.grid(row=6, column=0, sticky=tkinter.W, padx=140)
rascension = tkinter.Entry(main)
rascension["width"] = 10
rascension.grid(row=6, column=0, sticky=tkinter.W, padx=340)
# Declination eingeben:
labelD = tkinter.Label(main, text="Declination [°]:")
labelD.grid(row=6, column=0, sticky=tkinter.W, padx=410)
declination = tkinter.Entry(main)
declination["width"] = 10
declination.grid(row=6, column=0, sticky=tkinter.E, padx=500)

# Endlosschleife
main.mainloop()
