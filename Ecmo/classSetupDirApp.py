'''
Created on Oct 5, 2016

@author: fredrik
'''
from scipy.integrate import simps

import h5py
import numpy as np
import os
from os.path import dirname, join
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, Range1d, Plot, Span
from bokeh.models.widgets import Slider, Select, MultiSelect
from bokeh.models.glyphs import ImageURL
import time as timeModule

import pickle
#from matplotlib.pyplot import xlabel, ylabel

class SetupApp:
    
    def __init__(self):
        
        self.dataNumber = '100'

        self.referenceNetwork = 'Ecmo97_ref'
        self.ecmoNetwork = 'Ecmo97_aorta'
        
        
        self.cur = dirname(os.path.realpath('__file__'))
        print self.cur
        self.period = 0.8
        self.Pv = 5.*133.32
        self.vesselIdList =  []
        self.vesselId = 1
        self.node = 0
        
        print os.listdir(os.getcwd())
        #optAllVesselsWeighted_tot_constant_average_R1_preOpt.p

        
        referenceVesseldata = self.loadData(self.referenceNetwork, self.dataNumber, self.vesselId)
        
        time = referenceVesseldata['time']
        P = referenceVesseldata['P']
        Q = referenceVesseldata['Q']
        vesselName = referenceVesseldata['name']

        ecmoVesseldata = self.loadData(self.ecmoNetwork, self.dataNumber, self.vesselId)
        
        time_ecmo = ecmoVesseldata['time']
        P_ecmo = ecmoVesseldata['P']
        Q_ecmo = ecmoVesseldata['Q']
        
        self.sourceP = ColumnDataSource(data=dict(x=time, y=P))
        self.sourcePEcmoAorta = ColumnDataSource(data=dict(x=[time_ecmo], y=[P_ecmo]))
        self.sourcePEcmoAxillaris = ColumnDataSource(data=dict(x=[], y=[]))
        
        # Set up plotP
        self.plotP = Figure(plot_height=400, plot_width=600, title="Pressure",
                            x_axis_label="t [s]", y_axis_label="P [mmHg]",
                            tools="crosshair,pan,reset,resize,save,wheel_zoom"
                            )
        
        self.plotP.line('x', 'y', source=self.sourceP, color='black', line_width=2, legend='Reference (young adult)')
        self.plotP.line('x', 'y', source=self.sourcePEcmoAorta, color='green', line_width=2, legend='ECMO (aorta)')
        self.plotP.line('x', 'y', source=self.sourcePEcmoAxillaris, color='black', line_width=2, legend='ECMO (axillaris)')
        
        
        

        self.sourceQ = ColumnDataSource(data=dict(x=time, y=Q))
        self.sourceQEcmoAorta = ColumnDataSource(data=dict(x=[time_ecmo], y=[Q_ecmo]))
        self.sourceQEcmoAxillaris = ColumnDataSource(data=dict(x=[], y=[]))
        
        
        # Set up plotQ
        self.plotQ = Figure(plot_height=400, plot_width=600, title="flow",
                            x_axis_label="t [s]", y_axis_label="Q [ml/s]",
                            tools="crosshair,pan,reset,resize,save,wheel_zoom",
                            )
        
        self.plotQ.line('x', 'y', source=self.sourceQ, color='black', line_width=2)
        self.plotQ.line('x', 'y', source=self.sourceQEcmoAorta, color='green', line_width=2)
        self.plotQ.line('x', 'y', source=self.sourceQEcmoAxillaris, color='red', line_width=2)
        
        
        # Set up network plot loaded from image
        localUrl = "Ecmo/static/96model_first.svg"
        
        self.countImages = 0
        
        self.imageUrl = localUrl
        
        self.sourceImg = ColumnDataSource(data=dict(url = [localUrl]))
         
        self.plot_Img = Figure(plot_width=600, plot_height=800, title=vesselName, 
                               tools="crosshair, pan, reset, resize, save, wheel_zoom, box_zoom")
        self.plot_Img.x_range = Range1d(start=0, end=1)
        self.plot_Img.y_range = Range1d(start=0, end=1)
        self.plot_Img.image_url(url='url', x=0, y=1, h=1, w=1, source=self.sourceImg)
        
        # Set up widgets
        vesselList = []
        for Id in range(1, 98):
            vesselList.append(str(Id))
        
        self.internalVessels = [1, 2, 3, 4, 5, 7, 9, 14, 15, 18, 19, 21, 23, 27, 28, 29, 30, 35, 37, 39, 41, 42, 43, 44, 46, 50, 52]
        
        #self.sectionSelect = Select(title='section', value="upperLeft", options=sectionList)
        self.vesselIdSelect = Select(title='VesselID', value="1", options=vesselList)
        self.showVesselIdSelect = Select(title='Show vesselID on networks', value="No", options=["No", "Yes"])
        self.nodeSlider = Slider(title="node", value=0, start=0, end=2, step=1)
        self.canuleSelect = Select(title='Select site oof cannulation', value="aorta", options=['aorta', 'axillaris'])
        
        
        self.Widgetlist = [self.vesselIdSelect, self.showVesselIdSelect, self.nodeSlider, self.canuleSelect]
    
        self.initialize()
        
        listOfPlots = [self.plotP, self.plotQ]
        
        self.ajustFontSize(listOfPlots, '15pt', '10pt')
        self.plot_Img.title.text_font_size = '15pt'

        
    def initialize(self):
        
        self.vesselIdList =  []
        self.pathDict = None
        self.startLayerLines = None
        self.endLayerLines = None
        
        self.removeOldSVG()
        
        self.loadSVG()
    
    def ajustFontSize(self, listOfObjects, font_size, axis_size):
    
        
        for ob in listOfObjects:
            ob.xaxis.major_label_text_font_size = axis_size
            ob.yaxis.major_label_text_font_size = axis_size
            
            ob.xaxis.axis_label_text_font_size = font_size
            ob.yaxis.axis_label_text_font_size = font_size
    
    def update_data(self, attrname, old, new):
        

        #additionalPlot = self.additionalSelect.value
        self.node = self.nodeSlider.value
        self.assignDataNumber()
        
        
        vesselId = int(self.vesselIdSelect.value)
        self.vesselId = vesselId
        
        
        referenceVesseldata = self.loadData(self.referenceNetwork, self.dataNumber, vesselId, relativeFilePath=True)
            
        networkName = self.getNetworkName()
        vesselData = self.loadData(networkName, self.dataNumber, vesselId, relativeFilePath=True)
        print networkName
        
        showVesselId = self.showVesselIdSelect.value
        if showVesselId == "Yes":
            showNumbers = True
        else:
            showNumbers = False
        self.changeSVG(showNumbers=showNumbers)
        
        time = referenceVesseldata['time']
        P = referenceVesseldata['P']
        Pm = referenceVesseldata['Pm']
        R = referenceVesseldata['R']
        Pf = referenceVesseldata['Pf']
        Pb = referenceVesseldata['Pb']
        
        P_ecmo = vesselData['P']
        Pm_ecmo = vesselData['Pm']
        R_ecmo = vesselData['R']

        

        Q = referenceVesseldata['Q']
        Qm = referenceVesseldata['Qm']
        Qf = referenceVesseldata['Qf']
        
        Q_ecmo = vesselData['Q']
        Qm_ecmo = vesselData['Qm']

        
        RMS_P = round(self.calcRMS(P, P_ecmo, data_type="P")*100, 2)
        RMS_Q = round(self.calcRMS(Q, Q_ecmo, data_type="Q")*100, 2)
        
        
        self.plotP.title.text = "Pm = ({0}, {1}), R = ({2}, {3}); format = (reference, ecmo)".format(Pm, Pm_ecmo, R, R_ecmo)
        self.plot_Img.title.text = vesselData['name']
        self.plotQ.title.text = "Qm = ({0}, {1}), epsilon = {2}".format(Qm, Qm_ecmo, RMS_Q)
        
        self.sourceP.data = dict(x=time, y=P)

        self.sourcePEcmoAorta.data = dict(x=time, y=P_ecmo)

        self.sourceQ.data = dict(x=time, y=Q)
        self.sourceQEcmoAorta.data = dict(x=time, y=Q_ecmo)
        
        
        localUrl = self.imageUrl
        self.sourceImg.data = dict(url=[localUrl])
        self.plot_Img.image_url(url='url', x=0, y=1, h=1, w=1, source=self.sourceImg)
    
     
    def assignDataNumber(self):

        self.dataNumber = '100'
        
        
    def getNetworkName(self):
        canuleSite = self.canuleSelect.value

        if canuleSite == 'aorta':
            networkName = 'Ecmo97_aorta'
        elif canuleSite == 'axillaris':
            networkName = 'Ecmo97_axillaris'
        
        return networkName
            
            
        
    def loadData(self, networkName, dataNumber, vesselId, medical=True, relativeFilePath=False):
        
        
        file_name = networkName + "_SolutionData_" + dataNumber + ".hdf5"
        relFilePath = "Ecmo/data/" + file_name
        if relativeFilePath:
            relFilePath = relFilePath#"dirAppDev/" +  "data/" + networkName + "/SolutionData_" + dataNumber+ "/" + networkName + "_SolutionData_" + dataNumber + ".hdf5"

        hdf5File = h5py.File(relFilePath,'r')
        time = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
        
        node = int(self.node)
        dt = time[1] - time[0]
        period = self.period
        N = int(period/dt)
        N = len(time) - N - 1
        time = time[N:] - time[N]
        vesseldict = {}
        
        self.vesselIdList =  []
        for vesselName in hdf5File['vessels'].keys():
            
            tmpvesselId = vesselName.split(' - ')[-1]
            tmpvesselId = int(tmpvesselId)
            self.vesselIdList.append(tmpvesselId)
            if tmpvesselId == vesselId:
                
                P  = hdf5File['vessels'][vesselName]['Psol'][N:, node]
                Q  = hdf5File['vessels'][vesselName]['Qsol'][N:, node]
                A  = hdf5File['vessels'][vesselName]['Asol'][N:, node]
                
                Pf = hdf5File['vessels'][vesselName]['Psol_f'][N:, node]
                Pb = hdf5File['vessels'][vesselName]['Psol_b'][N:, node]
                Qf = hdf5File['vessels'][vesselName]['Qsol_f'][N:, node]
                Qb = hdf5File['vessels'][vesselName]['Qsol_b'][N:, node]
                
                Pm = self.findMean(time, P)
                Qm = self.findMean(time, Q)
                
                vName = vesselName.split(' - ')[0]
        if medical:
            #NB! no conversion of A
            P = P/133.32
            Pf = Pf/133.32
            Pb = Pb/133.32
            Qf = Qf*1e6
            Qb = Qb*1e6
            Pm = Pm/133.32
            Q = Q*1e6
            Qm = Qm*1e6
            Pv = self.Pv/133.32
            deltaP = Pm - Pv
            R = deltaP/Qm
        
        vesseldict['name'] = vName
        vesseldict['time'] = time
        vesseldict['P'] = P
        vesseldict['Pf'] = Pf
        vesseldict['Pb'] = Pb
        vesseldict['Q'] = Q
        vesseldict['Qf'] = Qf
        vesseldict['Qb'] = Qb
        vesseldict['Pm'] = round(Pm, 2)
        vesseldict['Qm'] = round(Qm, 2)
        vesseldict['R'] = round(R, 2)
        vesseldict['A'] = A
        
        return vesseldict
        

    def findMean(self, x, f):
        
        F = simps(f, x)
        
        f_mean = F/(x[-1] - x[0])
        
        return f_mean
    
    def calcRMS(self, dataFull, dataReduced, data_type="P"):
        
        if data_type == "P":
            RMS = np.sum(np.abs((dataReduced - dataFull)/dataFull))/len(dataFull)
        elif data_type == "Q":
            RMS = np.sum(np.abs((dataReduced - dataFull)/np.amax(dataFull)))/len(dataFull)
            
        return RMS
    
    def removeOldSVG(self, staticDirRel="/Ecmo/static"):

        staticDirAbs = ''.join([self.cur, staticDirRel])

        for fileName in os.listdir(staticDirAbs):
            if "_reduced" in fileName:
                print "removing: ", fileName
                fileNamePath = ''.join([staticDirAbs, '/', fileName])
                os.remove(fileNamePath)
    
    
    def loadSVG(self):
        
        """ Method thats load the """

        old_file = "Ecmo/static/96model.svg"
        f = open(old_file, "r")
        
        pathDict = {}
        
        path = 0
        pathLines = None
        addStartLayerLines = False
        addPaths = False
        for line in f:
            if "<g" in line:
                addStartLayerLines = True
                
                startLayerLines = []
            if "<path" in line:
                addPaths = True
                addStartLayerLines = False
                if path > 0:
                    pathDict[str(path)] = pathLines
                path += 1
                pathLines = []
        
            if "</g>" in line:
                addPaths = False
                endLayerLines = []
                endLayerLines.append(line)
            
            if addStartLayerLines:
                startLayerLines.append(line)
            if addPaths:
                pathLines.append(line)
            
            
            if ')">' in line:
                #addPaths = False
                addStartLayerLines = False
        
        f.close()
        pathDict[str(path)] = pathLines
        
        self.pathDict = pathDict
        self.startLayerLines = startLayerLines
        self.endLayerLines = endLayerLines
        
        
    def changeSVG(self, showNumbers=True):
        """ Method that copy the lines in old_file which contain svg paths for each vessels.
            Every vessel which is in the reduced model is then represented with a path in a second layer, 
            with a different color.
        """

        old_file = "Ecmo/static/96model.svg"
        vesselNumberFile = "Ecmo/static/vesselNumbers.svg"
        self.countImages += 1
        new_file = "Ecmo/static/96model_tmp" + str(timeModule.time()) + ".svg"
        previousFile = self.imageUrl
        print previousFile, "Ecmo/static/96model_first.svg", previousFile == "Ecmo/static/96model_first.svg"
        if "_tmp" in previousFile:
            os.remove(previousFile)
            
        self.imageUrl = new_file

        f = open(old_file, "r")
        f1 = open(vesselNumberFile, "r")
        f2 = open(new_file, "w")
        
        for oldLine in f:
            if "#Fill in Layer 2 Here" in oldLine:
        
                for line in self.startLayerLines:
                    if "Layer 1" in line:
                        line = line.replace("Layer 1", "Layer 2")
                        
                    elif "layer1" in line:
                        line = line.replace("layer1", "layer2")
                    f2.write(line)
                
                for n in range(len(self.pathDict)):
                    
                    if n + 1 in self.vesselIdList:
        
                        for line in self.pathDict[str(n + 1)]:
                            if "id=" in line:
                                line = '        id="path' + str(n+1) + '"\n'
                            if "#000000" in line:
                                if n + 1 == self.vesselId:
                                    line = line.replace("#000000", "#0000ff")
                                else:
                                    line = line.replace("#000000", "#ff0000")
                            
                            f2.write(line)
                f2.write(self.endLayerLines[0])
                
                if showNumbers:
                    for numberLines in f1:
                        f2.write(numberLines)
                
            else:
                f2.write(oldLine)
                
        f.close()
        f1.close()
        f2.close()
