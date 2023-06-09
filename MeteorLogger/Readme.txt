Python Script of Meteor Logger v.1.24a

This script runs under Windows (tested with Windows7/10). Under Linux there are not yet identified
issues with the audio handling.

Requirements:
Python 3.4 or higher must be installed on your computer. It is not tested whether Meteor Logger
also runs under a former version of python 3.4. At least the following modules are required:
- pyaudio
- numpy
- scipy
- matplotlib
- tkinter
- datetime
- struct
- threading
- copy

It is recommended to install a complete python-package. I developed Meteor Logger
under Windows 7/10 with WinPython (http://winpython.sourceforge.net/) and used the 32bit version
3.4.4.6Qt5 (https://sourceforge.net/projects/winpython/files/WinPython_3.4/3.4.4.6/) and
WinPython32-3.6.8.0Qt5 (https://sourceforge.net/projects/winpython/files/
WinPython_3.6/3.6.8.0/) respectively.


Windows:
How to install:
Copy <MeteorLogger_v1.24a.py> to any folder if the path to python.exe is
known to your OS (WinPython offers in the <Control paneel.exe> a simple way to register python
to Windows). Otherwise copy it into the folder that contains python.exe.

How to start:
Try one of the following methods:
-  If your python-package is fully registered to your OS double click <MeteorLogger_v1.24a.py>
   otherwise right click and select python.exe as program for opening .py-files.
-  Open IDLE(X) from the Windows start menue (or <IDLE(X).exe> in the python folder). In IDLE(X)
   open <MeteorLogger_v1.24a.py> and run it by pressing <F5>.
-  Under Windows open the Windows command line, change to the directory containing <MeteorLogger_v1.24a.py>
   and enter the following command: python MeteorLogger_v1.24a.py

If you encounter the error "pyaudio-OSError: [Errno -9999] Unanticipated host error" under Windows 10
you must enable the access to the microphone in Settings/System/Sound/Microphone.

