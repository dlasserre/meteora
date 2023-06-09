Python Script <ProcessData>

Requirements:
Python 3.4 or higher must be installed on your computer. It is not tested whether
<ProcessData> also runs under a former version of python 3.4. At least the following
modules are required:
- numpy
- scipy
- matplotlib
- tkinter
- datetime
- sys
- os
- math
- ephem

It is recommended to install a complete python-package. Meteor Logger was developped 
with WinPython (http://winpython.sourceforge.net/) and used the 32bit version
WinPython32-3.6.8.0Qt5 (https://sourceforge.net/projects/winpython/files/
WinPython_3.6/3.6.8.0/).

Pyephem is not part of the above mentioned package. It can be downloaded here:
http://rhodesmill.org/pyephem/ Installation is easily done with the winpython control
panel which is part of the WinPython-package.

How to install:
Copy <ProcessData_v*.**.py>, <config.txt> and <Manual_ProcessData> to any common folder 
if the path to python.exe is known to your OS. Otherwise copy them into the folder that 
contains python.exe. If your computer has a low screen resolution use
<ProcessData_v*.**_smallwindow.py instead.

How to start:
Try one of the following methods:
-  If your python-package is fully registered to your OS double click
   <ProcessData_v*.**.py> otherwise right click and select python.exe as program
   for opening .py-files.
-  Open IDLE(X) from the Windows start menue (or double click IDLE(X).exe in the
   python folder). In IDLE(X) open <ProcessData_v*.**.py> and run it by pressing <F5>.
-  Under Windows open the Windows command line, change to the directory containing
   <ProcessData_v*.**.py> and enter the following command: python ProcessData_v*.**.py


*.** = version number of ProcessData
