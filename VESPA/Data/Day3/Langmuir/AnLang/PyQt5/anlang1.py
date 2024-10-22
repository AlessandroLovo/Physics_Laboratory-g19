"""
ANLANG

Program to perform Langmuir probe analysis.

author: Emilio Martines
last edited: April 2017
"""

import sys, shelve, os, traceback
import QtGui
import QtCore
import QtWidgets
#from PyQt4 import QtGui, QtCore
import pyqtgraph as pg
import numpy as np
from scipy.optimize import curve_fit
from math import sqrt, log, pi
from emiliolib import readYoko
#from subprocess import Popen

class anLang(QtWidgets.QWidget):
#class anLang(QtWidgets.QMainWindow):
    
    def __init__(self):

#        self.hal = Popen('Theremino_HAL')
#        sleep(1)

        super(anLang, self).__init__()

        self.data_sources =['ASCII','Yokogawa','Synthetic']       
#        self.data_sources =['ASCII','Yokogawa','Theremino']       
        self.atomic_mass = {'Hydrogen': '1', 'Helium': '4', 'Nitrogen': '14', 
                            'Argon': '39.9', 'Xenon':  '131.3'}
        self.readConfig()
        self.initUI()
        self.analyzeData()
        self.plotGeneric()

        
    def initUI(self):

        QtWidgets.QToolTip.setFont(QtGui.QFont('SansSerif', 10))        
#        self.setToolTip('ANLANG program for Langmuir probe measurements')

        hbox_source = QtWidgets.QHBoxLayout()
        hbox_source.addWidget(QtWidgets.QLabel("Data source: "))
        self.combo_source = QtWidgets.QComboBox()
        for source in self.data_sources: 
            self.combo_source.addItem(source)
        index = self.combo_source.findText(self.source, QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.combo_source.setCurrentIndex(index)   
        self.combo_source.setToolTip('Source of Langmuir probe data')
        self.combo_source.currentIndexChanged.connect(self.setSource)
        hbox_source.addWidget(self.combo_source)

#        runButton = QtGui.QPushButton('Measure')
#        runButton.setToolTip('Perform measurement')
#        runButton.clicked.connect(self.sweepVoltage)
#        self.pbar = QtGui.QProgressBar()
#        self.pbar.setToolTip('Measurement progress indicator')
        readButton = QtWidgets.QPushButton('Read data')
        readButton.setToolTip('Read data from file')
        readButton.clicked.connect(self.readData)
#        writeButton = QtGui.QPushButton('Write to file')
#        writeButton.setToolTip('Write data to file')
#        writeButton.clicked.connect(self.writeData)
#        voltageButton = QtGui.QPushButton('View V')
#        voltageButton.setToolTip('View measured voltage values (for debug)')
#        voltageButton.clicked.connect(self.viewVoltage)
#        currentButton = QtGui.QPushButton('View I')
#        currentButton.setToolTip('View measured current values (for debug)')
#        currentButton.clicked.connect(self.viewCurrent)

        hbox_file = QtWidgets.QHBoxLayout()
        hbox_file.addWidget(QtWidgets.QLabel("Filename: "))
        self.filenameLabel = QtWidgets.QLabel(self.filename)
        hbox_file.addWidget(self.filenameLabel)

        self.TeLabel = QtWidgets.QLabel("Te = 0 eV")
        self.TeLabel.setToolTip('Electron temperature')
        self.IsLabel = QtWidgets.QLabel("Is = 0 mA")
        self.IsLabel.setToolTip('Ion saturation current')
        self.VfLabel = QtWidgets.QLabel("Vf = 0 V")
        self.VfLabel.setToolTip('Floating potential')
        self.RLabel = QtWidgets.QLabel("R = 0")
        self.RLabel.setToolTip('R parameter')
        self.nLabel = QtWidgets.QLabel("n = 0 m<sup>-3</sup>")
        self.nLabel.setToolTip('Plasma density')
        self.VpLabel = QtWidgets.QLabel("Vp = 0 V")
        self.VpLabel.setToolTip('Plasma potential')

        self.plotToolBar = QtWidgets.QToolBar()
        self.plotToolBar.addWidget(QtWidgets.QLabel('Graph: '))
        self.V = self.plotToolBar.addAction('V',self.plotVoltage)
        self.V.setCheckable(True)
        self.I = self.plotToolBar.addAction('I',self.plotCurrent)
        self.I.setCheckable(True)
        self.VI = self.plotToolBar.addAction('V-I',self.plotCharacteristic)
        self.VI.setCheckable(True)
        self.VI.setChecked(True)
        self.VIele = self.plotToolBar.addAction('V-Iel',self.plotCharacteristicEle)
        self.VIele.setCheckable(True)
        group = QtWidgets.QActionGroup(self)
        self.V.setActionGroup(group)
        self.I.setActionGroup(group)
        self.VI.setActionGroup(group)
        self.VIele.setActionGroup(group)

        self.fitToolBar = QtWidgets.QToolBar()
#        self.fitToolBar.addWidget(QtWidgets.QLabel('Fit type: '))
        self.threepar = self.fitToolBar.addAction('3 params',self.analyzeData)
        self.threepar.setCheckable(True)
        if self.method == 0: self.threepar.setChecked(True)
        self.linexp = self.fitToolBar.addAction('Lin-exp',self.analyzeData)
        self.linexp.setCheckable(True)
        if self.method == 1: self.linexp.setChecked(True)
        self.fourpar = self.fitToolBar.addAction('4 params',self.analyzeData)
        self.fourpar.setCheckable(True)
        if self.method == 2: self.fourpar.setChecked(True)
        group2 = QtWidgets.QActionGroup(self)
        self.threepar.setActionGroup(group2)
        self.linexp.setActionGroup(group2)
        self.fourpar.setActionGroup(group2)
    
#        hbox_slot_out = QtGui.QHBoxLayout()
#        hbox_slot_out.addWidget(QtGui.QLabel("Slot Vout: "))
#        line_slot_out = QtGui.QLineEdit('{0}'.format(self.slot_out))
#        line_slot_out.setInputMask('000')
#        line_slot_out.setToolTip('Theremino slot for sweep voltage')
#        line_slot_out.editingFinished.connect(
#            lambda: setattr(self, 'slot_out', line_slot_out.text()))
#        hbox_slot_out.addWidget(line_slot_out)
#
#        hbox_slot_V = QtGui.QHBoxLayout()
#        hbox_slot_V.addWidget(QtGui.QLabel("Slot V: "))
#        line_slot_V = QtGui.QLineEdit('{0}'.format(self.slot_V))
#        line_slot_V.setInputMask('000')
#        line_slot_V.setToolTip('Theremino slot for voltage measurement')
#        line_slot_V.editingFinished.connect(
#            lambda: setattr(self, 'slot_V', line_slot_V.text()))
#        hbox_slot_V.addWidget(line_slot_V)
#
#        hbox_slot_I = QtGui.QHBoxLayout()
#        hbox_slot_I.addWidget(QtGui.QLabel("Slot I: "))
#        line_slot_I = QtGui.QLineEdit('{0}'.format(self.slot_I))
#        line_slot_I.setInputMask('000')
#        line_slot_I.setToolTip('Theremino slot for current measurement')
#        line_slot_I.editingFinished.connect(
#            lambda: setattr(self, 'slot_I', line_slot_I.text()))
#        hbox_slot_I.addWidget(line_slot_I)

#         hbox_np = QtGui.QHBoxLayout()
#         hbox_np.addWidget(QtGui.QLabel("N. points: "))
#         self.line_np = QtGui.QLineEdit('{0}'.format(100))
#         self.line_np.setInputMask('000')
#         self.line_np.setToolTip('Number of points in I-V characteristic ' +
#             '(if changed a new measurement is required)')
#         self.line_np.editingFinished.connect(self.setNp)
#         hbox_np.addWidget(self.line_np)
# 
        hbox_gas = QtWidgets.QHBoxLayout()
        hbox_gas.addWidget(QtWidgets.QLabel("Gas: "))
        self.combo_gas = QtWidgets.QComboBox()
        for gas in self.atomic_mass: 
            self.combo_gas.addItem(gas)
        index = self.combo_gas.findText(self.gas, QtCore.Qt.MatchFixedString)
        if index >= 0:
            self.combo_gas.setCurrentIndex(index)   
        self.combo_gas.setToolTip('Gas used to produce the plasma')
        self.combo_gas.currentIndexChanged.connect(self.setGas)
        hbox_gas.addWidget(self.combo_gas)

        hbox_area = QtWidgets.QHBoxLayout()
        hbox_area.addWidget(QtWidgets.QLabel("Area (mm<sup>2</sup>): "))
        self.line_area = QtWidgets.QLineEdit('{0}'.format(self.area))
#        self.line_area.setInputMask('0000')
        self.line_area.setToolTip(
            'Collecting area of the Langmuir probe (in mm<sup>2</sup>)')
        self.line_area.editingFinished.connect(self.setArea)
        hbox_area.addWidget(self.line_area)

        hbox_shunt = QtWidgets.QHBoxLayout()
        hbox_shunt.addWidget(QtWidgets.QLabel("Rshunt (Ohm): "))
        self.line_shunt = QtWidgets.QLineEdit('{0}'.format(self.Rshunt))
        self.line_shunt.setToolTip(
            'Shunt resistor used for current measurement (in Ohm)')
        self.line_shunt.editingFinished.connect(self.setRshunt)
        hbox_shunt.addWidget(self.line_shunt)

        quitButton = QtWidgets.QPushButton('Quit')
#        quitButton.clicked.connect(QtCore.QCoreApplication.instance().quit)
        quitButton.clicked.connect(self.close)
        quitButton.setToolTip('Close program and exit')
#        quitButton.resize(quitButton.sizeHint())
#        quitButton.move(50, 50)       

        ## Switch to using white background and black foreground
        pg.setConfigOption('background', 'w')
        pg.setConfigOption('foreground', 'k')
        self.graph = pg.PlotWidget()

        hbox = QtWidgets.QHBoxLayout()
        vbox1 = QtWidgets.QVBoxLayout()
        vbox1.addLayout(hbox_source)
#        vbox1.addWidget(runButton)
#        vbox1.addWidget(self.pbar)
        vbox1.addWidget(readButton)
#        vbox1.addWidget(writeButton)
#        vbox1.addLayout(hbox_slot_out)
#        vbox1.addLayout(hbox_slot_V)
#        vbox1.addLayout(hbox_slot_I)
#        vbox1.addLayout(hbox_np)
        vbox1.addLayout(hbox_file)
        vbox1.addWidget(self.plotToolBar)
        vbox1.addWidget(self.fitToolBar)
        vbox1.addLayout(hbox_gas)
        vbox1.addLayout(hbox_area)
        vbox1.addLayout(hbox_shunt)
        vbox1.addWidget(self.TeLabel)
        vbox1.addWidget(self.IsLabel)
        vbox1.addWidget(self.VfLabel)
        vbox1.addWidget(self.RLabel)
        vbox1.addWidget(self.nLabel)
        vbox1.addWidget(self.VpLabel)
        vbox1.addWidget(quitButton)
        vbox1.addStretch(1)
        hbox.addLayout(vbox1)
        vbox2 = QtWidgets.QVBoxLayout()        
        vbox2.addWidget(self.graph)        
        hbox.addLayout(vbox2)
        self.setLayout(hbox)

#        self.resize(800, 450)
        self.resize(900, 650)
        self.centerWindow()
        self.setWindowTitle('anLang')
#        self.setWindowIcon(QtGui.QIcon('web.png'))        
    
        self.show()


    def centerWindow(self):
        qr = self.frameGeometry()
        cp = QtGui.QDesktopWidget().availableGeometry().center()
        qr.moveCenter(cp)
        self.move(qr.topLeft())


    def readConfig(self):
        self.config = shelve.open('.anlang',protocol=2) #protocol serve per Python2 -> Python3
        if 'area' in self.config: self.area = self.config['area']
        else: self.area = 7.07 #sq. mm
        if 'rshunt' in self.config: self.Rshunt = self.config['rshunt']
        else: self.Rshunt = 1. #Ohm
        if 'source' in self.config: self.source = self.config['source']
        else: self.source = 'Yokogawa'
        if 'gas' in self.config: self.gas = self.config['gas']
        else: self.gas = 'Argon'
        if 'filename' in self.config: self.filename = self.config['filename']
        else: self.filename = ' '
        if 'method' in self.config: self.method = self.config['method']
        else: self.method = 0
#        if self.config.has_key('slot_out'): self.slot_out = self.config['slot_out']
#        else: self.slot_out = 550
#        if self.config.has_key('slot_V'): self.slot_V = self.config['slot_V']
#        else: self.slot_V = 551
#        if self.config.has_key('slot_I'): self.slot_I = self.config['slot_I']
#        else: self.slot_I= 552
        if 'trange' in self.config: self.trange = self.config['trange']
        else: self.trange = [0.,100.]
        if 'vrangeIon' in self.config: self.vrangeIon = self.config['vrangeIon']
        else: self.vrangeIon = [-50.,0.]
        if 'vrangeEle' in self.config: self.vrangeEle = self.config['vrangeEle']
        else: self.vrangeEle = [0.,100.]
        if 'time' in self.config: self.time = self.config['time']
        else: self.time = np.linspace(0.,1.,num=100)
        if 'voltage' in self.config: 
            self.voltageRaw = self.config['voltage']
        else: 
            self.voltageRaw = np.linspace(-50.,50.,num=100)
        if 'current' in self.config: 
            self.currentRaw = self.config['current']
        else: 
            self.currentRaw = np.zeros(self.voltageRaw.size)

        
    def writeConfig(self):
        self.config['area'] = self.area
        self.config['rshunt'] = self.Rshunt
#        self.config['slot_out'] = self.slot_out
#        self.config['slot_V'] = self.slot_V
#        self.config['slot_I'] = self.slot_I
        self.config['trange'] = self.trange
        self.config['vrangeIon'] = self.vrangeIon
        self.config['vrangeEle'] = self.vrangeEle
        self.config['time'] = self.time
        self.config['voltage'] = self.voltageRaw
        self.config['current'] = self.currentRaw
        self.config['source'] = self.source
        self.config['gas'] = self.gas
        self.config['filename'] = self.filename
        self.config['method'] = self.method
        self.config.close()


    def readData(self):
        if self.source == 'Synthetic':
            self.filename = ' '
            self.fakeData()
        else:
            fileName = QtGui.QFileDialog.getOpenFileName(self, 
                           "Open input file")[0]
            try:
                self.filename = os.path.basename(os.path.splitext(fileName)[0])
                if self.source == 'ASCII':
                    data = np.genfromtxt(fileName, skip_header=1)
                    self.time = data[:,0]
                    self.voltageRaw = data[:,2]
                    self.currentRaw = data[:,1]
                elif self.source == 'Yokogawa':
                    fileName, fileExtension = os.path.splitext(fileName)
                    if fileExtension == '.WVF' or fileExtension == '.HDR':
                        data, time = readYoko(fileName)
                        self.time = time
                        self.voltageRaw = data[:,0]
                        self.currentRaw = data[:,1]
                    else:
                        QtGui.QMessageBox.critical(self, "Message", 
                            "Wrong file type: " + fileExtension)
                else:
                    print('Error in source specification')
            except Exception as e:
                QtGui.QMessageBox.critical(self, "Message", str(e))
                traceback.print_exc()
        self.filenameLabel.setText(self.filename)
        # Adjust time range selection to new data
        self.trange = [max(self.trange[0],self.time.min()),
                       min(self.trange[1],self.time.max())]
        self.analyzeData()
        self.plotGeneric()


    def setSource(self):
        self.source = str(self.combo_source.currentText())
        

    def setGas(self):
        self.gas = str(self.combo_gas.currentText())
        self.computeDerived()


    def setArea(self):
        self.area = float(self.line_area.text())
        self.computeDerived()

    def setRshunt(self):
        self.Rshunt = float(self.line_shunt.text())
        self.analyzeData()
        self.plotGeneric()


#    def sweepVoltage(self):
#        np = 200
#        for index in range(0,np):
#            Theremino.write_slot(self.slot_out, index)
#            self.voltageRaw[index] = Theremino.read_slot_no_nan(self.slot_V)
#            self.currentRaw[index] = Theremino.read_slot_no_nan(self.slot_I)
#            sleep(0.002)
#            self.pbar.setValue(index/float(np)*100) # Update progress bar
#        self.pbar.setValue(0)
#        self.fakeData()  # DA ELIMINARE
#        self.analyzeData()
#        self.plotGeneric()


    def fakeData(self):
        npoints = 100
        self.time = np.linspace(0, 0.1, num=npoints)
        self.voltageRaw = np.linspace(-50, 50, num=npoints)
        self.currentRaw = self.langmuirCharacteristic(self.voltageRaw, 
                                                      2., 0.005, 10.)
        np.place(self.currentRaw, self.currentRaw > 0.1, 0.1)
        self.currentRaw = self.currentRaw * np.random.normal(1., 0.1, npoints)


    def langmuirCharacteristic(self, V, Te, Is, Vf):
        return Is*(np.exp((V-Vf)/Te)-1.)


    def langmuirCharacteristic4par(self, V, Te, Is, Vf, R):
        return Is*(1.+R*(V-Vf))*(np.exp((V-Vf)/Te)-1.)


    def plotGeneric(self):
        if self.VI.isChecked(): self.plotCharacteristic()
        elif self.V.isChecked(): self.plotVoltage()
        elif self.I.isChecked(): self.plotCurrent()
        elif self.VIele.isChecked(): self.plotCharacteristicEle()
        else: print('Error in checking action from plot type menu')


    def plotVoltage(self):
        self.graph.clear()
        self.graph.plot(self.time, self.voltageRaw, symbol='o', symbolSize=4)
        self.graph.setTitle('Voltage')
        self.graph.setLabel('bottom', 'Time', 's')
        self.graph.setLabel('left', 'Voltage', 'V')
        self.graph.setLogMode(y=False)
        self.region = pg.LinearRegionItem(self.trange)
        self.graph.addItem(self.region, ignoreBounds=True)
        self.region.sigRegionChangeFinished.connect(self.regionChangedTrace)

    
    def plotCurrent(self):
        self.graph.clear()
        self.graph.plot(self.time, self.currentRaw/self.Rshunt, 
                        symbol='o', symbolSize=4)
        self.graph.setTitle('Current')
        self.graph.setLabel('bottom', 'Time', 's')
        self.graph.setLabel('left', 'Current', 'A')
        self.graph.setLogMode(y=False)
        self.region = pg.LinearRegionItem(self.trange)
        self.graph.addItem(self.region, ignoreBounds=True)
        self.region.sigRegionChangeFinished.connect(self.regionChangedTrace)


    def plotCharacteristic(self):
        self.graph.clear()
        self.graph.plot(self.voltage, self.current, symbol='o', symbolSize=4,
                        pen=None)
        self.graph.setTitle('I-V characteristic')
        self.graph.setLabel('bottom', 'Voltage', 'V')
        self.graph.setLabel('left', 'Current', 'A')
        self.graph.setLogMode(y=False)
        self.region = pg.LinearRegionItem(self.vrangeIon)
#        self.region.setZValue(10)
        self.graph.addItem(self.region, ignoreBounds=True)
        self.region.sigRegionChangeFinished.connect(self.regionChangedIon)
        if self.threepar.isChecked(): 
            sel = ((self.voltage >= self.vrangeIon[0]) & 
                   (self.voltage <= self.vrangeIon[1]))
            self.graph.plot(self.voltage[sel], self.langmuirCharacteristic(
                self.voltage[sel], self.Te, self.Is ,self.Vf), 
                symbol=None, pen=pg.mkPen(color='r', width=2))            
        elif self.linexp.isChecked(): 
            self.graph.plot(self.voltage, self.currentIon, symbol=None, 
                        pen=pg.mkPen(color='r', width=2))
            self.graph.addLine(x=self.Vf, pen=pg.mkPen(color='g', 
                        width=1, style=QtCore.Qt.DashLine))
            self.graph.addLine(y=0, pen=pg.mkPen(color='g', 
                        width=1, style=QtCore.Qt.DashLine))
        elif self.fourpar.isChecked():
            sel = ((self.voltage >= self.vrangeIon[0]) & 
                   (self.voltage <= self.vrangeIon[1]))
            self.graph.plot(self.voltage[sel], self.langmuirCharacteristic4par(
                self.voltage[sel], self.Te, self.Is ,self.Vf, self.R), 
                symbol=None, pen=pg.mkPen(color='r', width=2))            
        else: print('Error in checking action from fit type menu')


    def plotCharacteristicEle(self):
        sel = np.where((self.currentEle > 0.))
        self.graph.clear()
        self.graph.plot(self.voltage[sel], self.currentEle[sel], symbol='o', 
                        symbolSize=4, pen=None)
        self.graph.setTitle('I-V characteristic - electron component')
        self.graph.setLabel('bottom', 'Voltage', 'V')
        self.graph.setLabel('left', 'Electron Current', 'A')
        self.graph.setLogMode(y=True)
        self.region = pg.LinearRegionItem(self.vrangeEle)
#        self.region.setZValue(10)
        self.graph.addItem(self.region, ignoreBounds=True)
        self.region.sigRegionChangeFinished.connect(self.regionChangedEle)
        if self.linexp.isChecked(): 
            self.graph.plot(self.voltageEleFit, self.currentEleFit, symbol=None, 
                            pen=pg.mkPen(color='r', width=2))
            self.graph.addLine(x=self.Vf, pen=pg.mkPen(color='g', 
                            width=1, style=QtCore.Qt.DashLine))


    def regionChangedTrace(self):
#        self.region.setZValue(10) # replace region au premier plan 
        self.trange = self.region.getRegion() # recupera il range selezionato
        self.analyzeData()


    def regionChangedIon(self):
        self.vrangeIon = self.region.getRegion() # recupera il range selezionato
        self.analyzeData()
        self.plotGeneric()


    def regionChangedEle(self):
        self.vrangeEle = self.region.getRegion() # recupera il range selezionato
        self.analyzeData()
        self.plotGeneric()


    def analyzeData(self):
        # Select time window
        self.selection = (self.time >= self.trange[0]) & (self.time <= self.trange[1])
        self.voltage = self.voltageRaw[self.selection]
        self.current = self.currentRaw[self.selection] / self.Rshunt
        d = np.argsort(self.voltage)
        self.voltage = self.voltage[d]
        self.current = self.current[d]
        # Adjust voltage range selections to new data
        self.vrangeIon = [max(self.vrangeIon[0],self.voltage.min()),
                       min(self.vrangeIon[1],self.voltage.max())]
        self.vrangeEle = [max(self.vrangeEle[0],self.voltage.min()),
                       min(self.vrangeEle[1],self.voltage.max())]
        # Find floating potential
        try:
            Vf1 = self.voltage[np.where(self.current >= 0.)[0][0]]
            Vf2 = self.voltage[np.where(self.current <= 0.)[0][-1]]
            self.Vf = (Vf1 + Vf2) / 2.
#           print Vf1, Vf2, self.Vf
        except:
            self.Vf = 0.
        sel = ((self.voltage >= self.vrangeIon[0]) & 
               (self.voltage <= self.vrangeIon[1]))
        self.R = 0.
        if self.threepar.isChecked():
            self.method = 0
            p0 = [2., -self.current[0], self.Vf]
            popt, pcov = curve_fit(self.langmuirCharacteristic, 
                self.voltage[sel], self.current[sel], p0=p0, method='lm')
            perr = np.sqrt(np.diag(pcov))
            self.Te = popt[0]
            self.Is = popt[1]
            self.Vf = popt[2]
            self.currentIon = np.ones(self.voltage.size)*self.Is
            self.currentEle = self.current - self.currentIon
        elif self.linexp.isChecked(): 
            self.method = 1
            # Fit ion current and compute electron current
            try:
                sel = ((self.voltage >= self.vrangeIon[0]) & 
                       (self.voltage <= self.vrangeIon[1]))
                voltageIonFit = self.voltage[sel]
                currentIonFit = self.current[sel]
                coefficients = np.polyfit(voltageIonFit,currentIonFit,1)
                polynomial = np.poly1d(coefficients)
                self.currentIon = polynomial(self.voltage)
            except:
                self.currentIon = np.zeros(self.voltage.size)
            self.currentEle = self.current - self.currentIon
            # Fit electron current
            try:
                sel = ((self.voltage >= self.vrangeEle[0]) &
                       (self.voltage <= self.vrangeEle[1]) &
                       (self.currentEle > 0.))
                self.voltageEleFit = self.voltage[sel]
                self.currentEleFit = self.currentEle[sel]
                coefficients = np.polyfit(self.voltageEleFit,
                                          np.log(self.currentEleFit),1)
                self.Te = 1./coefficients[0]
                polynomial = np.poly1d(coefficients)
                self.Is = np.exp(polynomial(self.Vf))
                self.currentEleFit = np.exp(polynomial(self.voltageEleFit))
            except:
                self.Te = 0.
                self.Is = 0.
        elif self.fourpar.isChecked():
            self.method = 2
            p0 = [2., -self.current[0], self.Vf, 0.]
            popt, pcov = curve_fit(self.langmuirCharacteristic4par, 
                self.voltage[sel], self.current[sel], p0=p0, method='lm')
            perr = np.sqrt(np.diag(pcov))
            self.Te = popt[0]
            self.Is = popt[1]
            self.Vf = popt[2]
            self.R = popt[3]
        else: print('Error in checking action from fit type menu')           
        self.TeLabel.setText('Te = {0:.4f} eV'.format(self.Te))
        self.IsLabel.setText('Is = {0:.6f} mA'.format(self.Is*1000.))
        self.VfLabel.setText('Vf = {0:.4f} V'.format(self.Vf))
        self.RLabel.setText('R = {0:.4f}'.format(self.R))
        self.computeDerived()
        self.plotGeneric()


    def computeDerived(self):
        e = 1.6022e-19
        mp = 1.67E-27
        me = 9.1e-31
        mi = mp * float(self.atomic_mass[self.gas])
        alpha = 0.5*(log(mi/(2*pi*me))+1.)
        try:
            cs = sqrt(e*self.Te/mi)
            self.n = 2*self.Is/(e*cs*self.area*1E-6)
        except:
            self.n = 0.
        self.Vp = self.Vf + alpha*self.Te
        self.nLabel.setText('n = {0:.2e} m<sup>-3</sup>'.format(self.n))
        self.VpLabel.setText('Vp = {0:.2f} V'.format(self.Vp))
                

    def close(self):
        self.writeConfig()
#        self.hal.terminate()
        QtCore.QCoreApplication.instance().quit()


def main():    
#   Create QApplication if it doesnt exist 
#    app = QtGui.QApplication.instance() # checks if QApplication already exists 
#    if not app: app = QtGui.QApplication(sys.argv)
#    app.aboutToQuit.connect(app.deleteLater)
#    an = anLang()
#    sys.exit(app.exec_())
    app = QtWidgets.QApplication(sys.argv)
    an = anLang()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
