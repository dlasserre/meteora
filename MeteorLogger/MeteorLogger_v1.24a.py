# -*- coding: utf-8 -*-
# Meteor detection and logging
# copyright 2019, Wolfgang Kaufmann
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

import pyaudio
import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from scipy import signal
from scipy import fftpack
from scipy import integrate
from scipy.optimize import minimize_scalar
import scipy.stats
import tkinter, tkinter.filedialog, tkinter.scrolledtext, tkinter.messagebox
import datetime
import struct
import threading
import copy
from numpy import nanmean

# Steuervariable:
sampFreq = 48000
chunkSize = 2048  # => fDist = 23,4 Hz
shift = 512       # => tDist = 10,6 ms
oLap = int(chunkSize / shift)   # Überlappung
tDist = (chunkSize/sampFreq) * (shift/chunkSize)  # zeitl. Auflösung in s
fDist = sampFreq / chunkSize    # Frequenzauflösung in Hz
fSlope = 11000 * tDist          # Max. beob. f-Änderung bei Meteor: 11000 Hz/s
soundformat = pyaudio.paInt16
channelnr = 1
# Globale Variable
fileName = "MeteorData.csv"
RUN = False
graphics1 = []
graphics2 = []
IBuffer = []


class GUI(object):
    def __init__(self, master):
        global fmin, fmax, dt, lStatus, lCounts, auDev, dMode, supp
        global lActive, resp, notch, fNotch
        master.title("Meteor Logger v1.24a")
        # Audio-Device Liste aufrufen
        devList = listAudio()
        # Menue-Leiste
        mBar = tkinter.Menu(master)
        # Dateipfad setzen
        mPath = tkinter.Menu(mBar)
        mPath["tearoff"] = 0
        mPath.add_command(label = "Set Directory", command = setPath)
        mPath.add_separator()
        mPath.add_command(label = "Quit", command = beenden)
        mBar.add_cascade(label = "File", menu = mPath)
        # Select Audio Device
        mDevice = tkinter.Menu(mBar)
        mDevice["tearoff"] = 0
        auDev = tkinter.IntVar()
        try:
            auDev.set(devList[0][0])
        except:
            tkinter.messagebox.showerror \
                ("Error", "No accessible audio input-device found")
            quit()
        for devNr in range(0, len(devList)):
            mDevice.add_radiobutton(label = str(devList[devNr][1]), variable = auDev,
                                               value = devList[devNr][0])
        mBar.add_cascade(label = "Audio Device", menu = mDevice)
        # Select Detector-Sensitivity
        mSens = tkinter.Menu(mBar)
        mSens["tearoff"] = 0
        dMode = tkinter.IntVar()
        dMode.set(6)
        mSens.add_radiobutton(label = "sensitive", variable = dMode,
                                value = 5)
        mSens.add_radiobutton(label = "robust", variable = dMode,
                                value = 6)
        mBar.add_cascade(label = "Detection", menu = mSens)
        # Select Auto-Notch
        mANotch = tkinter.Menu(mBar)
        mANotch["tearoff"] = 0

        mSupp = tkinter.Menu(mANotch)
        mSupp["tearoff"] = 0
        supp = tkinter.IntVar()
        supp.set = 0
        mSupp.add_radiobutton(label = "off", variable = supp, value = 0)
        mSupp.add_radiobutton(label = "fast", variable = supp, value = 50)
        mSupp.add_radiobutton(label = "middle", variable = supp, value = 150)
        mSupp.add_radiobutton(label = "slow", variable = supp, value = 300)
        mSupp.add_radiobutton(label = "very slow", variable = supp, value = 600)
        mANotch.add_cascade(label = "Response Time", menu = mSupp)

        mResp = tkinter.Menu(mANotch)
        mResp["tearoff"] = 0
        resp = tkinter.IntVar()
        resp.set = 0
        mResp.add_radiobutton(label = "low", variable = resp, value = 0)
        mResp.add_radiobutton(label = "medium", variable = resp, value = 1)
        mResp.add_radiobutton(label = "high", variable = resp, value = 2)
        mANotch.add_cascade(label = "Sensitivity", menu = mResp)

        mBar.add_cascade(label = "Auto-Notch", menu = mANotch)
        # Hilfe
        mHelp = tkinter.Menu(mBar)
        mHelp["tearoff"] = 0
        mHelp.add_command(label = "Instruction", command = printHelp)
        mHelp.add_command(label = "About Meteor Logger", command = printAbout)
        mBar.add_cascade(label = "Help", menu = mHelp)
        # Menue-Leiste sichtbar machen
        master["menu"] = mBar
        # Überschrift
        label0 = tkinter.Label(master, text="Meteor Logger")
        label0["font"] = "Arial 16"
        label0["width"] = 40
        label0.grid(row=0, column=0, columnspan=2, sticky="we", pady=5)
        # Festlegen Frequenzgrenzen für die Spektrumsdarstellung
        # und die Suche der max. Magnitude
        label7 = tkinter.Label(master, text="<= Set lower/upper frequency limit in Hz =>")
        label7["font"] = "Arial 10"
        label7.grid(row=11, column=0, columnspan=2, pady=5)
        fmin = tkinter.Entry(master)
        fmax = tkinter.Entry(master)
        fmin.grid(row=11, column=0, sticky=tkinter.W, padx=50)
        fmax.grid(row=11, column=1, sticky=tkinter.E, padx=50)
        # Starte Analyse
        b3 = tkinter.Button(master,  text = "Start Monitoring", command=startThread)
        b3.grid(row=13, column=0, sticky=tkinter.W, padx=50)
        # Analyse stoppen
        b4 = tkinter.Button(master,  text = "Stop Monitoring", command=stopAnalyse)
        b4.grid(row=14, column=0, sticky=tkinter.W, padx=50)
        # Aktivität Auto-Suppression anzeigen
        lActive = tkinter.Label(master)
        lActive["fg"] = "#0000FF"
        lActive.grid(row=13, column=0, columnspan=2)
        # Grafik ausgeben
        b5 = tkinter.Button(master,  text = "Show Diagram", command=plotDiagram)
        b5.grid(row=13, column=1, sticky=tkinter.E, padx=50, pady=10)
        # Notch-Frequenzeingabe
        label9 = tkinter.Label(master, text="Notch-Frequ. [Hz]:")
        label9["font"] = "Arial 10"
        label9.grid(row=14, column=0, sticky="E", pady=10)
        fNotch = tkinter.Entry(master)
        fNotch.grid(row=14, column=1, sticky="W", pady=10)
        # Status anzeigen
        lStatus = tkinter.Label(master)
        lStatus.grid(row=14, column=1, sticky=tkinter.E, padx=50)
        # Datenausgabe auf dem Bildschirm
        dt = tkinter.scrolledtext.ScrolledText(master, width=70, height=30)
        dt.grid(row=15, column=0, columnspan=2, padx=50, pady=10)
        # Meteor-Counts/h anzeigen
        lCounts = tkinter.Label(master)
        lCounts.grid(row=16, column=0, sticky=tkinter.W, padx=50)
        # Beenden Button
        b2 = tkinter.Button(master, text="Quit", command=beenden)
        b2.grid(row=17, column=0, sticky=tkinter.W, padx=50, pady=5)

def setPath():
    global fileName
    dirName = tkinter.filedialog.askdirectory()
    # Cancel-Button pressed?
    if len(dirName) == 0:
        return
    else:
        fileName = dirName + "\\" + "MeteorData.csv"


def listAudio():
    p = pyaudio.PyAudio()
    devList = []
    info = p.get_host_api_info_by_index(0)
    numdevices = info.get("deviceCount")
    for i in range(0, numdevices):
        if(p.get_device_info_by_host_api_device_index(0,i).get("maxInputChannels"))>0:
            devList.append([i, p.get_device_info_by_host_api_device_index(0,i).get("name")])
    return devList

def openSound():
    global SOUND, stream
    SOUND = pyaudio.PyAudio()  # Sound-Objekt erzeugen
    stream = SOUND.open(format=soundformat,
                    channels=channelnr,
                    rate=sampFreq,
                    input=True,
                    input_device_index=auDev.get(),
                    frames_per_buffer=chunkSize)

def closeSound():
    stream.stop_stream()
    stream.close()
    SOUND.terminate()

def record():
    data = stream.read(chunkSize)
    y = np.array(struct.unpack("%dh" % (chunkSize * channelnr), data))
    return y

def openFile():
    try:
        dfile = open(fileName, "a")
    except:
        tkinter.messagebox.showerror \
                ("Error", "Datafile could not be opened")
        return
    return dfile

def getMNotch(Fmin, Fmax):
    fn_ = fNotch.get()
    try:
        fn = float(fn_)
        if (fn < Fmin) or (fn > Fmax):
            # Falls fNotch außerhalb Analysegrenzen mit -1 besetzten
            fn = -1.0
    except:
        fn = -1.0
    return fn
    

def autoNotch(window, kmin, kmax, Fmin, Fmax, y):
    global lActive
    p = doFFT(window, y)
    # Finde die drei Frequenzen mit der höchsten Magnitude
    fM = findMax(p[kmin:kmax], kmin)
    # Stelle fest, ob ein Signal oder Rauschen vorliegt
    fD = detectPeak(fM)
    # Manual-Notch an?
    fn = getMNotch(Fmin, Fmax)
    if fn > 0:
        IFrequenz = fn
    else:
        # Prüfe auf Vorhandensein von Interferenz
        IFrequenz = detectI(fD[0], Fmin, Fmax)
    if IFrequenz > 0:
        lActive["text"] = "Notch active: " + str(round(IFrequenz, 1)) + " Hz"
        pos = int(round(IFrequenz / fDist))
        p[pos-1] = p[pos-1]/10000
        p[pos] = p[pos]/10000
        p[pos+1] = p[pos+1]/10000
    else:
        lActive["text"] = " "
    return p

def getFlimit():
    global fmin, fmax
    # Untere f-Grenze einlesen
    f1 = fmin.get()
    try:
        Fmin = int(f1)
        if (Fmin < 0) or (Fmin > int(sampFreq / 4)):
            Fmin = 0
    except:
        Fmin = 0
    # obere f-Grenze einlesen
    f2 = fmax.get()
    try:
        Fmax = int(f2)
        if (Fmax > int(sampFreq / 4)) or (Fmax < 0):
            Fmax = int(sampFreq / 4)
    except:
        Fmax = 4000
    # Prüfen auf richtige Reihenfolge
    if Fmin > Fmax:
        f = Fmax
        Fmax = Fmin
        Fmin = f
    # Mindestbreite erzwingen
    fIst = Fmax - Fmin
    fSoll = 65 * fDist
    fDiff = fSoll - fIst
    if fDiff > 0:
        Fmax = int(Fmax + fDiff)
    # Korrigierte Ausgabe im GUI
    fmin.delete(0, last=50)
    fmax.delete(0, last=50)
    fmin.insert("end", Fmin)
    fmax.insert("end", Fmax)
    # f in Array-Index umwandeln
    kmin = int(Fmin * (chunkSize / sampFreq))
    kmax = int(Fmax * (chunkSize / sampFreq))
    return kmin, kmax, Fmin, Fmax

def doFFT(window, y):
    y = window * y
    p = fftpack.fft(y)
    nUniquePts = int(np.ceil((chunkSize+1)/2.0))
    p = p[0:nUniquePts]
    p = abs(p)          # Umwandlung komplexes fft-Ergebnis in die Amplitude
    p = p / chunkSize   # scale by the number of points
    p = p**2            # square it to get the power
    if chunkSize % 2 > 0:
        p[1:len(p)] = p[1:len(p)] * 2           # we've got odd number of points fft
    else:
        p[1:len(p) - 1] = p[1:len(p) - 1] * 2   # we've got even number of points fft
    return p

def findMax(p, kmin):
    pos = []
    fM = []
    fPoints = 3
    # Indices der drei größten Werte finden
    isort = np.argsort(p)
    pos.append(int(isort[-1:]))
    for i in range(-1,-fPoints,-1):
        pos.append(int(isort[i-1:i]))
    # Rauschen errechnen als Median in p
    noise = scipy.stats.median_abs_deviation(p)
    # Frequenz [Hz] zu Index bestimmen
    # und paarweise abspeichern: [f,power,noise]
    for i in range(fPoints):
        fM.append([((kmin + pos[i]) * fDist), p[pos[i]], noise])
    # Listenpaare aufsteigend nach Frequenz sortieren
    fM.sort()
    return fM

def listMeteor(dB, ap, dFile, nr):
    # Formatvorgaben definieren
    ft4 = "{: >12.2E}"                  # Noise
    ft3 = "{:0>5d}"                     # Number
    ft2 = "{: >11.2E}"                  # Power
    ft1 = "{: >9.1f}"                   # Frequenz
    if not ap:
        header = "Nr-----Date---------Time--------------f/Hz------Signal"
        header = f"{header}-------Noise\r\n"
        dt.insert("end", header)
    signal = f"{ft3.format(nr)}  {dB[0]}   {dB[1]} "
    signal = f"{signal}{ft1.format(dB[2][0])} {ft2.format(dB[2][1])}"
    signal = f"{signal}{ft4.format(dB[2][2])}\r\n"
    dt.insert("end", signal)
    dt.see("end")
    # Wenn Listeneinträge 1 Mio Zeichen überschreiten, erste Einträge löschen
    n = list(dt.count("0.0", "end"))
    if n[0] > 1000000:
        dt.delete("0.0", "5.70")
    # Daten in file schreiben
    signal = f"{nr},{dB[0]},{dB[1]},{dB[2][0]},{dB[2][1]},{dB[2][2]}\r"
    dFile.write(signal)

def getTime(index, DT1, DT2):
    deltaT = (DT2 - DT1) / oLap
    lTime = str(DT1 + (deltaT * index))
    YYMMDD = lTime[:10]
    HHMMSS = lTime[11:23]
    return YYMMDD, HHMMSS

def detectPeak(fM):
    # x, y-vectors for curve-fit
    x = [fM[0][0], fM[1][0], fM[2][0]]
    y = [fM[0][1], fM[1][1], fM[2][1]]
    # Bestimme, ob f1, f2, f3 innerhalb 5*fDist liegen
    nb1 = abs(fM[1][0] - fM[0][0])
    nb2 = abs(fM[1][0] - fM[2][0])
    nb = 2.1 > round(((nb1 + nb2) / (2 * fDist)),1)
    if nb:
        # polynomial curve fit (2 degrees)
        z = np.polyfit(x, y, 2)
        ff = np.poly1d(z)
        # find maximum within the data boundaries
        sol = minimize_scalar(lambda x: -ff(x), bounds=[x[0],x[2]], method='bounded')
        # fmax interpoliert, Power und Noise wie gemessen zurückgeben
        return [sol["x"], fM[1][1], fM[1][2]]
        # return fM[1]
    else:
        # f=-1, p=0, noise wie gemessen zurückgeben
        return [-1.0,0,fM[1][2]]

def detectI (f, Fmin, Fmax):
    global IBuffer
    iFac = [[4.61E-06, 0.19, 50],[5.6E-06, 0.235, 60.8],[6.59E-06, 0.28, 71.6]]
    iDur = supp.get()
    iResp = resp.get()
    iLimit = -iFac[iResp][0] * (Fmax-Fmin)**2 + iFac[iResp][1] * (Fmax-Fmin)\
             - iFac[iResp][2]
    # Auto-Suppression aus
    if iDur == 0:
        return -1.0
    # Auto-Suppression an
    if f > 0:
        IBuffer.append(f)
    if len(IBuffer) > iDur:
        del IBuffer[:1]
    if len(IBuffer) == iDur:
        fI = scipy.stats.median_abs_deviation(IBuffer)
        stdI = sc.std(IBuffer)
        if stdI < iLimit:
            return fI
        else:
            return -1.0
    else:
        return -1.0

def markOutlier(dB, nFollow, r, shift, index):
    fArray = []
    for m in range(nFollow):
        if m != index:
            fArray.append(dB[m][2][0])
    fMedian = scipy.stats.median_abs_deviation(fArray)
    # Outlier markieren
    if abs(dB[r-shift][2][0] - fMedian) < abs(dB[r+1][2][0] - fMedian):
        dB[r+1][2][0] = -1.0
        dB[r+1][2][1] = 0
        return r+1
    else:
        dB[r-shift][2][0] = -1.0
        dB[r-shift][2][1] = 0
        return r-shift

def checkTimeline(dB, nFollow):
    # Führe die Analyse für die =>nFollow Datenreihen durch:
    # 1. Bestimme, ob mind. nFollow-1 Signale in den Datenreihen vorliegen.
    #    Signal: f >= 0, kein Signal (gap): f = -1
    sCrit = 0
    for j in range(nFollow):
        Crit = np.sign(dB[j][2][0])
        if Crit == 0:
            Crit = 1
        sCrit = sCrit + Crit
    if sCrit < (nFollow-2):
        return False
    # 2. Bestimme den Grad der f-Variabilität der peak-f bezüglich fSlope
    #    Wenn ein Datensatz mit f = -1 besteht, dessen Index feststellen:
    index = -1
    for j in range(nFollow):
        if dB[j][2][0] < 0.0:
            index = j
    #    Prüfen von fSlope unter Umgehung eines gaps:
    for j in range(nFollow-1):
        if (j != index) and (j+1 != index):
            if (abs(dB[j][2][0] - dB[j+1][2][0])) > fSlope:
                pos = markOutlier(dB, nFollow, j, 0, index)
                if index > -1:
                    return False
                else:
                    index = pos
        if (j == index) and (j != 0):
            if (abs(dB[j-1][2][0] - dB[j+1][2][0])) > (2*fSlope):
                pos = markOutlier(dB, nFollow, j, 1, index)
                return False
    return True

def analyseFmax(p,kmin,kmax,Fmin,Fmax,allPrinted,dBuffer,dFile,index,DT1,DT2,
                nr,nFollow, switch, cycle):
    global lActive
    # Hole Uhrzeit
    YYMMDD, HHMMSS = getTime(index, DT1, DT2)
    # Finde die drei Frequenzen mit der höchsten Magnitude
    fM = findMax(p[kmin:kmax], kmin)
    # Stelle fest, ob ein Signal oder Rauschen vorliegt
    fD = detectPeak(fM)
    # Speichere die Datenreihe im Schieberegister
    del dBuffer[0]
    dBuffer.append([YYMMDD,HHMMSS,fD])
    # Daten für die graphische Darstellung bereitstellen
    pFFT = [fM[0][0],fM[1][0],fM[2][0]]
    # Prüfe Zeitreihe, ob ein Meteorsignal vorliegt
    if not switch:
        if checkTimeline(dBuffer[0:nFollow], nFollow):
            if not allPrinted:
                nr = nr + 1
                for k in range(nFollow):
                    listMeteor(dBuffer[k], allPrinted, dFile, nr)
                    allPrinted = True
            else:
                listMeteor(dBuffer[nFollow-1], allPrinted, dFile, nr)
        else:
            # Prüfen, ob nur eine einfache Signalunterbrechung vorliegt
            if allPrinted:
                if checkTimeline(dBuffer[nFollow:], nFollow):
                    listMeteor(dBuffer[nFollow-1], allPrinted, dFile, nr)
                else:
                    allPrinted = False
                    switch = True
    # nach Signalabruch nFollow-mal leer durchlaufen, um Puffer mit neuen
    # Daten zu füllen
    else:
        cycle = cycle + 1
        if cycle == nFollow:
            switch = False
            cycle = 0
    # Funktion mit Rückgaben verlassen
    return pFFT, allPrinted, nr, switch, cycle

def countM(aP, nr, DT2, nrc, counts, h):
    # nr-Wechsel bedeutet neue Meteor-Detektion
    if nr > nrc:
        nrc = nr
        counts = counts + 1
    # Count-Ausgabe aktualisieren
    lCounts["text"] = "Hour: " + h + "  =>  " + str(counts) + " counts"
    # Zum Stundenwechsel Zählerstand abspeichern und Zähler zurücksetzen
    if (str(DT2)[11:13] != h) and not aP:
        h = str(DT2)[11:13]
        counts = 0
    return nrc, counts, h

def runAnalyse(window, kmin, kmax, Fmin, Fmax, dFile):
    global RUN, graphics1, graphics2, IBuffer
    del graphics1[:]
    del graphics2[:]
    del IBuffer[:]
    # Hole Detektions-Empfindlichkeit (Anzahl aufeinander folgender slices
    # für die Detektionsanalyse):
    nFollow = dMode.get()
    dBuffer = []
    aP = False
    nr = 0
    counts = 0
    nrc = 0
    switch = False
    cyc = 0
    # dBuffer (fungiert als Schieberegister) vorbesetzen
    for j in range(nFollow*2):
        dBuffer.append(["x","x",[-1,0,0]])
    yA = []
    yB = []
    yC = []
    # get first chunk
    yA = record()
    DT1 = datetime.datetime.now()
    h = str(DT1)[11:13]
    while RUN:
        # get next chunk
        yC = record()
        DT2 = datetime.datetime.now()
        # FFT of sliced chunks out of yA and yC
        for i in range(oLap):
            yB = yA[shift*i:]
            yB = np.append(yB, yC[:shift*i])
            p = autoNotch(window, kmin, kmax, Fmin, Fmax, yB)
            pF,aP,nr,switch,cyc = analyseFmax(p,kmin,kmax,Fmin,Fmax,aP,dBuffer,
                                         dFile,i,DT1,DT2,nr,nFollow,switch,cyc)
            graphics1.append(pF)
            graphics2.append(0)
            # Eintrag bei Signaldetektion
            if aP:
                x1 = len(graphics2) - 2*nFollow
                for m in range(nFollow):
                    graphics2[x1 + m] = 1
            # count meteors/h
            nrc, counts, h = countM(aP, nr, DT2, nrc, counts, h)
        # next run
        yA = yC
        DT1 = DT2
        # remain Graphic size at 500 points
        if len(graphics1) > 500:
            del graphics1[:oLap]
            del graphics2[:oLap]

def analyseAudio():
    global RUN
    RUN = True   # Echtzeitanalyse
    # Audio Stream öffenen
    openSound()
    # Data file öffnen
    dataFile = openFile()
    #Frequenzgrenzen der Meteorsuche holen
    kmin, kmax, Fmin, Fmax = getFlimit()
    # Vorbereitung für Rohdaten-Speicherung als .csv-file
    ## hText = "date,time,f1 [Hz],magn1 [dB],f2 [Hz],magn2 [dB],f3 [Hz],magn3 [dB]\r"
    ## dataFile.write(hText)
    # Vorbereitung für die Meteordaten-Speicherung als csv-file
    hText = "Nr,Date,Time,f [Hz],Signal,Noise\r"
    dataFile.write(hText)
    # Window-Funktion berechnen
    window = signal.blackmanharris(chunkSize)
    # FFT-Analyse
    runAnalyse(window, kmin, kmax, Fmin, Fmax, dataFile)
    # Audio Stream beenden
    closeSound()
    # Data file schließen
    dataFile.close()

def plotDiagram():
    # statisches Abbild von graphics erstellen
    d1 = copy.deepcopy(graphics1)
    d2 = copy.deepcopy(graphics2)
    plt.close()
    plt.figure(1)
    # Grafik der drei Frequenzen mit den höchsten Magnituden pro Slice
    ax = plt.subplot(211)
    plt.title("Peak-Frequencies per Slice")
    plt.ylabel("Frequency[Hz]")
    plt.plot(d1,"o")
    # Grafik der entdeckten Meteorsignale
    plt.subplot(212, sharex = ax)
    plt.title("Meteor-Detection")
    plt.xlabel("Slice-Number")
    plt.ylabel("Detection")
    plt.axis([0,500,0,1.05])
    plt.plot(d2)
    plt.show()

def startThread():
    global lStatus, Process
    lStatus["text"] = "Monitoring is running ... "
    lStatus["bg"] = "#80D057"  # green
    if RUN == False:
        Process = threading.Thread(target=analyseAudio)
        Process.start()

def stopAnalyse():
    global RUN, lStatus
    RUN = False
    lStatus["text"] = "Monitoring is stopped ... "
    lStatus["bg"] = "#FFCE50"  # yellow

def printHelp():
    global helpWindow
    helpWindow = tkinter.Tk()
    helpWindow.title("Instruction")
    t = tkinter.scrolledtext.ScrolledText(helpWindow, width = 65, height = 40,
                                          padx = 5, pady = 5)
    t.pack()
    helptext =  "Instruction\r\n\n\
    IMPORTANT: Exit Meteor Logger only by its QUIT-Button or\r\n\
    QUIT-Menue entry! Otherwise you may risk to loose all data!\r\n\n\
    General Information: Meteor Logger has a time resolution of\r\n\
    10.7 ms, a fft-frequency resolution of 23,4 Hz and runs with\r\n\
    a sampling rate of 48 kHz. To avoid lossy resampling it\r\n\
    is recommended to adapt the output sampling rate of a\r\n\
    SDR-software to 48 kHz also.\r\n\n\
    A signal must last at least 40-50 ms to be registered.\r\n\
    To a small amount the detection algorithm is fooled by\r\n\
    noise. Noise produces from time to time a random event\r\n\
    that fulfills the criteria of a real signal and therefore\r\n\
    will be registered. These erroneous registrations stand\r\n\
    out by divergent frequencies and low SNR. The amount of\r\n\
    such fake-signals depends on the selected bandwidth (upper\r\n\
    - lower frequency) and the detection mode. The smaller\r\n\
    the bandwidth the higher the error rate. White noise\r\n\
    gives rise to an error rate of about 0.7/hour in\r\n\
    robust mode and 1.5/hour in sensitive mode at a\r\n\
    bandwidth of 2500 Hz.\r\n\n\
    Data: Meteor Logger automatically writes a datafile (.csv)\r\n\
    to the directory where it is executed from. You can\r\n\
    choose any other directory via menue.\r\n\n\
    Operation: Selecting frequency limits, detection mode and\r\n\
    directory should be done only before monitoring is (re)\r\n\
    started.\r\n\n\
    Selecting lower/upper frequency determines the\r\n\
    frequency range in which a signal detection is performed.\r\n\
    Because the detection algorithm depends on the variability\r\n\
    of frequency to distinguish between noise and signal the\r\n\
    minimal bandwidth is limited to 1523 Hz. The lower/upper\r\n\
    frequency limits of audio input must correspond or exceed\r\n\
    the selected limits within Meteor Logger.\r\n\n\
    Detector mode can be set to \'sensitive\' or \'robust\'.\r\n\
    The robust mode has a more rigorous signal checking.\r\n\
    This leads to a reduced detection level of signals com-\r\n\
    pared to the sensitive mode. In return robust mode pro-\r\n\
    duces less fake-signals and is less prone to inter-\r\n\
    ference.\r\n\n\
    Meteor Logger has a built in auto-notch function. Within\r\n\
    the selected bandwidth it eliminates the strongest\r\n\
    continuous signal. The algorithm can not distinguish\r\n\
    between a meteor signal and an unwanted interference\r\n\
    signal. Its only criterion is the duration of the present\r\n\
    strongest signal. You can select between different time\r\n\
    delays before notching starts. So if you want to study\r\n\
    longer lasting meteor signals without getting them trun-\r\n\
    cated by Auto-Notch select \'slow\' or \'very slow\'. Auto-\r\n\
    Notch allows the unaffected registration of meteor signals\r\n\
    at the presence of interference as long as the frequency\r\n\
    of the interference is +/- 50 Hz off the frequency of\r\n\
    the meteor signal.\r\n\n\
    Also you manually can enter a notch frequency. It overrides\r\n\
    the Auto-Notch-function.\r\n\n\
    Diagnostics: <Show diagram> displays a plot of the last\r\n\
    5 s: the upper window shows the distribution of the three\r\n\
    frequencies with the strongest magnitudes within each\r\n\
    timeslot. The lower window indicates the position of\r\n\
    detected signals. You can zoom into the graph to get a\r\n\
    detailed insight.\r\n"
    t.insert("end", helptext)
    t.tag_add("HIGH", "1.0", "1.20")
    t.tag_add("IMPORTANT", "3.0", "3.15")
    t.tag_add("MARK", "6.0", "6.25")
    t.tag_add("MARK", "25.0", "25.10")
    t.tag_add("MARK", "29.0", "29.15")
    t.tag_add("MARK", "65.0", "65.17")
    t.tag_config("HIGH", font=('Arial', 12, 'bold', 'italic'))
    t.tag_config("MARK", foreground = "blue")
    t.tag_config("IMPORTANT", foreground="red")

def printAbout():
    global aboutWindow
    aboutWindow = tkinter.Tk()
    aboutWindow.title("About")
    t = tkinter.Text(aboutWindow, width=65, height=23, padx=5, pady=5)
    t.pack()
    aboutText = "Meteor Logger v1.24a\r\n\n\
    This software detects and loggs radio meteor signals\r\n\
    from an audio-stream.\r\n\n\
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
    t.tag_add("HIGH", "1.0", "1.20")
    t.tag_add("MARK", "6.0", "7.40")
    t.tag_config("HIGH", font=('Arial', 12, 'bold', 'italic'))
    t.tag_config("MARK", foreground = "blue")

def beenden():
    global RUN, lStatus
    if RUN:
        RUN = False
        lStatus["text"] = "Monitoring is stopped ... "
        lStatus["bg"] = "#FFCE50"  # yellow
    answer = tkinter.messagebox.askokcancel \
                ("Warning", "Close Meteor Logger?")
    if answer == 1:
        plt.close()
        try:
            helpWindow.destroy()
        except:
            pass
        try:
            aboutWindow.destroy()
        except:
            pass
        mainWindow.destroy()


# Hauptfenster erzeugen
mainWindow=tkinter.Tk()
app = GUI(mainWindow)
#Endlosschleife
mainWindow.mainloop()
