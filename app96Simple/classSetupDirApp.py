'''
Created on Oct 5, 2016

@author: fredrik
'''
from scipy.integrate import simps
from scipy import interpolate

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

        self.reductionNetwork = 'Full96Model_baseline'
        
        
        self.cur = dirname(os.path.realpath('__file__'))
        print self.cur
        self.period = 0.8
        self.Pv = 5.*133.32
        self.vesselIdList =  []
        self.vesselId = 1
        self.node = 0
        
        print os.listdir(os.getcwd())
        #optAllVesselsWeighted_tot_constant_average_R1_preOpt.p

        self.aortaRMSPressureDict = pickle.load(open("app96Simple/optDicts/aortaRMSPressureDict.p", 'rb'))
        self.brachialPulseDict = pickle.load(open("app96Simple/optDicts/brachialPulsePressureDict.p", 'rb'))
        self.carotidRMSFlowDict = pickle.load(open("app96Simple/optDicts/carotidRMSFlowDict.p", 'rb'))
        self.internalCarotidRMSFlowDict = pickle.load(open("app96Simple/optDicts/internalCarotidRMSPressureDict.p", 'rb'))
        self.femoralRMSPressureDict = pickle.load(open("app96Simple/optDicts/femoralRMSPressureDict.p", 'rb'))
        self.anteriorCommunicatingFlowDict = pickle.load(open("app96Simple/optDicts/anteriorCommunicatingFlowDict.p", 'rb'))
        self.anteriorCommunicatingPressureDict = pickle.load(open("app96Simple/optDicts/anteriorCommunicatingPressureDict.p", 'rb'))
        self.rightMiddleCerebralFlowDict = pickle.load(open("app96Simple/optDicts/rightMiddleCerebralFlowDict.p", 'rb'))
        self.rightMiddleCerebralPressureDict = pickle.load(open("app96Simple/optDicts/rightMiddleCerebralPressureDict.p", 'rb'))

        self.leftPosteriorCerebralPressureDict = pickle.load(open("app96Simple/optDicts/leftPosteriorCerebralPressureDict.p", 'rb'))
        self.leftPosteriorCerebralFlowDict = pickle.load(open("app96Simple/optDicts/leftPosteriorCerebralFlowDict.p", 'rb'))

        self.rightAnteriorCerebralFlowDict = pickle.load(open("app96Simple/optDicts/rightAnteriorCerebralFlowDict.p", 'rb'))
        self.rightAnteriorCerebralPressureDict = pickle.load(open("app96Simple/optDicts/rightAnteriorCerebralPressureDict.p", 'rb'))
        
        baseVesseldata = self.loadData(self.reductionNetwork, self.dataNumber, self.vesselId)
        
        time = baseVesseldata['time']
        P = baseVesseldata['P']
        Q = baseVesseldata['Q']
        vesselName = baseVesseldata['name']
        self.sourceP = ColumnDataSource(data=dict(x=time, y=P))
        self.sourcePf = ColumnDataSource(data=dict(x=[], y=[]))
        self.sourcePb = ColumnDataSource(data=dict(x=[], y=[]))
        
        self.sourceP_reduced = ColumnDataSource(data=dict(x=time, y=P))

        self.sourcePf_reduced = ColumnDataSource(data=dict(x=[], y=[]))
        self.sourcePb_reduced = ColumnDataSource(data=dict(x=[], y=[]))
        
        
        # Set up plotP
        self.plotP = Figure(plot_height=250, plot_width=300, title="Pressure",
                            x_axis_label="t [s]", y_axis_label="P [mmHg]",
                            tools="crosshair,pan,reset,resize,save,wheel_zoom"
                            )
        
        self.plotP.line('x', 'y', source=self.sourceP, color='black', line_width=2, legend='Full model')
        self.plotP.line('x', 'y', source=self.sourcePf, color='black', line_width=2, line_dash='dashed')
        self.plotP.line('x', 'y', source=self.sourcePb, color='black', line_width=2, line_dash='dotted')
        
        
        self.plotP.line('x', 'y', source=self.sourceP_reduced, color='red', line_width=2, legend='Reduced model')
        self.plotP.line('x', 'y', source=self.sourcePf_reduced, color='red', line_width=2, line_dash='dashed')
        self.plotP.line('x', 'y', source=self.sourcePb_reduced, color='red', line_width=2, line_dash='dotted')
        

        self.sourceQ = ColumnDataSource(data=dict(x=time, y=Q))
        self.sourceQf = ColumnDataSource(data=dict(x=[], y=[]))
        self.sourceQb = ColumnDataSource(data=dict(x=[], y=[]))
        
        self.sourceQ_reduced = ColumnDataSource(data=dict(x=time, y=Q))
        self.sourceQf_reduced = ColumnDataSource(data=dict(x=[], y=[]))
        self.sourceQb_reduced = ColumnDataSource(data=dict(x=[], y=[]))
        
        # Set up plotQ
        self.plotQ = Figure(plot_height=250, plot_width=300, title="flow",
                            x_axis_label="t [s]", y_axis_label="Q [ml/s]",
                            tools="crosshair,pan,reset,resize,save,wheel_zoom",
                            )
        
        self.plotQ.line('x', 'y', source=self.sourceQ, color='black', line_width=2)
        self.plotQ.line('x', 'y', source=self.sourceQf, color='black', line_width=2, line_dash='dashed')
        self.plotQ.line('x', 'y', source=self.sourceQb, color='black', line_width=2, line_dash='dotted')
               
        self.plotQ.line('x', 'y', source=self.sourceQ_reduced, color='red', line_width=2)
        self.plotQ.line('x', 'y', source=self.sourceQf_reduced, color='red', line_width=2, line_dash='dashed', legend='forward')
        self.plotQ.line('x', 'y', source=self.sourceQb_reduced, color='red', line_width=2, line_dash='dotted', legend='backward')
        
        
        # Set up network plot loaded from image
        localUrl = "app96Simple/static/96model_first.svg"
        
        self.countImages = 0
        
        self.imageUrl = localUrl
        
        self.sourceImg = ColumnDataSource(data=dict(url = [localUrl]))
         
        self.plot_Img = Figure(plot_width=300, plot_height=600, title=vesselName, 
                               tools="crosshair, pan, reset, resize, save, wheel_zoom, box_zoom")
        self.plot_Img.x_range = Range1d(start=0, end=1)
        self.plot_Img.y_range = Range1d(start=0, end=1)
        self.plot_Img.image_url(url='url', x=0, y=1, h=1, w=1, source=self.sourceImg)
        
        # Set up widgets
        vesselList = []
        for Id in range(1, 97):
            vesselList.append(str(Id))
        constraintValues_array = np.arange(0.1, 3.5, 0.1)
        constrantValues = []
        for constrainValue in constraintValues_array:
            constrantValues.append(str(round(constrainValue, 1)))
        
        
        self.internalVessels = [1, 2, 3, 4, 5, 7, 9, 14, 15, 18, 19, 21, 23, 27, 28, 29, 30, 35, 37, 39, 41, 42, 43, 44, 46, 50, 52]
        
        self.vesselIdSelect = Select(title='VesselId', value="1", options=vesselList)
        

        self.waveSplitSelect = Select(title='Show waveSplit', value="False", options=["False", "True"])
        
        self.minimizeSelect = Select(title='minimization case:', value="Aorta (Pressure epsilon avg)", 
                                     options=['Aorta (Pressure epsilon avg)', 
                                              'Carotid (Flow epsilon avg)', 
                                              'InternalCarotid (Pressure epsilon avg)', 
                                              'Brachial (Pressure Pulse)', 
                                              'Femoral (Pressure epsilon avg)', 
                                              'ant. Com. (Pressure epsilon avg)', 
                                              'ant. Com. (Flow epsilon avg)', 
                                              'r. M. Cerebral (Pressure epsilon avg)', 
                                              'r. M. Cerebral (Flow epsilon avg)', 
                                              'l. P. Cerebral (Pressure epsilon avg)', 
                                              'l. P. Cerebral (Flow epsilon avg)', 
                                              'r. ant. Cerebral (Pressure epsilon avg)',
                                              'r. ant. Cerebral (Flow epsilon avg)'])
        
        self.minimizationDicts = {'Aorta (Pressure epsilon avg)': self.aortaRMSPressureDict, 
                                  'Carotid (Flow epsilon avg)': self.carotidRMSFlowDict, 
                                  'InternalCarotid (Pressure epsilon avg)': self.internalCarotidRMSFlowDict, 
                                  'Brachial (Pressure Pulse)': self.brachialPulseDict,
                                  'Femoral (Pressure epsilon avg)': self.femoralRMSPressureDict, 
                                  'ant. Com. (Pressure epsilon avg)': self.anteriorCommunicatingPressureDict,
                                  'ant. Com. (Flow epsilon avg)': self.anteriorCommunicatingFlowDict, 
                                  'r. M. Cerebral (Pressure epsilon avg)': self.rightMiddleCerebralPressureDict,
                                  'r. M. Cerebral (Flow epsilon avg)': self.rightMiddleCerebralFlowDict, 
                                  'l. P. Cerebral (Pressure epsilon avg)': self.leftPosteriorCerebralPressureDict, 
                                  'l. P. Cerebral (Flow epsilon avg)': self.leftPosteriorCerebralFlowDict,
                                  'r. ant. Cerebral (Pressure epsilon avg)':self.rightAnteriorCerebralPressureDict,
                                  'r. ant. Cerebral (Flow epsilon avg)':self.rightAnteriorCerebralFlowDict}
        
        #self.constraintSlider = Slider(title='Choose constraint', value=0.1, start=0.1, end=3.5, step=0.1)
        self.constraintSlider = Select(title='Choose constraint', value='0.1', options=constrantValues)
        
        
        self.Widgetlist = [self.waveSplitSelect, self.minimizeSelect, self.constraintSlider]
    

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
        
        
        networkName = self.getNetworkName()
        vesselId = self.vesselId
        baseVesseldata = self.loadData(self.reductionNetwork, self.dataNumber, vesselId, relativeFilePath=True)
            
        

        vesselData = self.loadData(networkName, self.dataNumber, vesselId, relativeFilePath=True)
        print networkName
        self.changeSVG()
        
        time = baseVesseldata['time']
        P = baseVesseldata['P']
        Pm = baseVesseldata['Pm']
        R = baseVesseldata['R']
        Pf = baseVesseldata['Pf']
        Pb = baseVesseldata['Pb']
        
        P_reduced = vesselData['P']
        Pm_reduced = vesselData['Pm']
        R_reduced = vesselData['R']
        Pf_reduced = vesselData['Pf']
        Pb_reduced = vesselData['Pb']
        

        Q = baseVesseldata['Q']
        Qm = baseVesseldata['Qm']
        Qf = baseVesseldata['Qf']
        Qb = baseVesseldata['Qb']
        
        Q_reduced = vesselData['Q']
        Qm_reduced = vesselData['Qm']
        Qf_reduced = vesselData['Qf']
        Qb_reduced = vesselData['Qb']
        
        RMS_P = round(self.calcRMS(P, P_reduced, data_type="P")*100, 2)
        RMS_Q = round(self.calcRMS(Q, Q_reduced, data_type="Q")*100, 2)
        
        Nv_baseline = baseVesseldata['Nv']
        Nv_reduced = vesselData['Nv']
        
        self.plotP.title.text = "Pm = ({0}, {1}), R = ({2}, {3}); format = (Full, Reduced), epsilon = {4}".format(Pm, Pm_reduced, R, R_reduced, RMS_P)
        self.plot_Img.title.text = vesselData['name']
        self.plotQ.title.text = "Qm = ({0}, {1}), epsilon = {2}. N = ({3}, {4})".format(Qm, Qm_reduced, RMS_Q, Nv_baseline, Nv_reduced)
        
        self.sourceP.data = dict(x=time, y=P)

        self.sourceP_reduced.data = dict(x=time, y=P_reduced)

        self.sourceQ.data = dict(x=time, y=Q)
        self.sourceQ_reduced.data = dict(x=time, y=Q_reduced)
        
        
        if self.waveSplitSelect.value == "True":
            
            
            
            self.sourcePf.data = dict(x=time, y=Pf)
            self.sourcePb.data = dict(x=time, y=Pb)

            self.sourcePf_reduced.data = dict(x=time, y=Pf_reduced)
            self.sourcePb_reduced.data = dict(x=time, y=Pb_reduced)
            
            
            self.sourceQf.data = dict(x=time, y=Qf)
            self.sourceQb.data = dict(x=time, y=Qb)

            self.sourceQf_reduced.data = dict(x=time, y=Qf_reduced)
            self.sourceQb_reduced.data = dict(x=time, y=Qb_reduced)
            
            
            #Pb = Pb - P[0]
            
            #Pb_reduced = Pb_reduced - P_reduced[0]
            
            A = baseVesseldata['A']
            dPf = 133.32*(Pf[1:] - Pf[0:-1])
            dQf = 1e-6*(Qf[1:] - Qf[0:-1])
            dUf = dQf/A[0:-1]
            
            dPb = 133.32*(Pb[1:] - Pb[0:-1])
            dQb = 1e-6*(Qb[1:] - Qb[0:-1])
            dUb = dQb/A[0:-1]

            A_reduced = vesselData['A']
            dPf_reduced = 133.32*(Pf_reduced[1:] - Pf_reduced[0:-1])
            dQf_reduced = 1e-6*(Qf_reduced[1:] - Qf_reduced[0:-1])
            dUf_reduced = dQf_reduced/A_reduced[0:-1]
            
            dPb_reduced = 133.32*(Pb_reduced[1:] - Pb_reduced[0:-1])
            dQb_reduced = 1e-6*(Qb_reduced[1:] - Qb_reduced[0:-1])
            dUb_reduced = dQb_reduced/A_reduced[0:-1]
            
            dt = time[1] - time[0]
            
            WI_f = 1e-6*dPf*dUf/(dt**2)
            WI_b = 1e-6*dPb*dUb/(dt**2)
            
            WI_f_reduced = 1e-6*dPf_reduced*dUf_reduced/(dt**2)
            WI_b_reduced = 1e-6*dPb_reduced*dUb_reduced/(dt**2)
            
        
        else:
            
            self.sourcePf.data = dict(x=[], y=[])
            self.sourcePb.data = dict(x=[], y=[])

            self.sourcePf_reduced.data = dict(x=[], y=[])
            self.sourcePb_reduced.data = dict(x=[], y=[])

            self.sourceQf.data = dict(x=[], y=[])
            self.sourceQb.data = dict(x=[], y=[])

            self.sourceQf_reduced.data = dict(x=[], y=[])
            self.sourceQb_reduced.data = dict(x=[], y=[])


        
        
        localUrl = self.imageUrl
        self.sourceImg.data = dict(url=[localUrl])
        self.plot_Img.image_url(url='url', x=0, y=1, h=1, w=1, source=self.sourceImg)
    
    def assignAddionalData(self):
        
        method = self.methodSelect.value
        additionalPlot = self.additionalSelect.value
        minSize = 5
        maxSize = 30
        relradiusMin = 0.12925170068
        Ntmp = 31
        relRadius = np.linspace(relradiusMin, 1, Ntmp)
        markerSizeList = np.linspace(minSize, maxSize, Ntmp)
        if method == 'Total':
            case = self.optParamsDictLumpedMean
            key_C = 'C'
        else:
            case = self.optParamsDictAlgebraicMean
            key_C = 'Cw'
        
        radius_rel_all = []
        
        R1_div_R2_0_all = []
        R1_div_R2_all = []
        R1_div_R2_init_div_opt_all = []
        
        R1_rel_all = []
        C_rel_all = []
        marker_all = []
        
        radius_rel_vessel = []
        
        R1_rel_vessel = []
        C_rel_vessel = []
        marker_vessel = []
        
        R1_rel_vessel_phys = []
        C_rel_vessel_phys = []

        R1_div_R2_0_vessel = []
        R1_div_R2_vessel = []
        R1_div_R2_init_div_opt_vessel = []
        
        node = -1

        r_end_aorta = self.lumpedValues[1]['radius'][node]
        
        if additionalPlot in ['ParamScatter', 'R1_div_R2_opt', 'R1_div_R2_init', 'R1_div_R2_init_div_opt']:
            vesselColor = self.vesselId
            
            for vessel in self.internalVessels:
                
                r_end = self.lumpedValues[vessel]['radius'][node]
                radius_rel_all.append(r_end/r_end_aorta)
                #if vessel == 1:
                #print r_end, r_end/r_end_aorta
                R1_0 = case[vessel][node]["Wk3"]['R1LC']['R1_0']
                C_0 = case[vessel][node]["Wk3"]['R1LC']['C_0']
                R2_0 = case[vessel][node]["Wk3"]['R1LC']['R2_0']
                
                R1 = case[vessel][node]["Wk3"]['R1LC']['R1']
                C = case[vessel][node]["Wk3"]['R1LC']['C']
                R2 = case[vessel][node]["Wk3"]['R1LC']['R2']
                
                for n, relRadiusTmp in enumerate(relRadius):
                    if relRadiusTmp >= r_end/r_end_aorta:
                        n_use = n
                        break
                R1_rel_all.append(R1_0/R1)
                C_rel_all.append(C_0/C)
                marker_all.append(markerSizeList[n])
                
                R1_div_R2_0_all.append(R1_0/R2_0)
                R1_div_R2_all.append(R1/R2)
                R1_div_R2_init_div_opt_all.append((R1_0/R2_0)/(R1/R2))
                
                if vessel == vesselColor:
                    radius_rel_vessel.append(r_end/r_end_aorta)
                    
                    R1_rel_vessel.append(R1_0/R1)
                    C_rel_vessel.append(C_0/C)
                    marker_vessel.append(markerSizeList[n])
                    
                    R1_div_R2_0_vessel.append(R1_0/R2_0)
                    R1_div_R2_vessel.append(R1/R2)
                    R1_div_R2_init_div_opt_vessel.append((R1_0/R2_0)/(R1/R2))
                    
                    if additionalPlot == 'ParamScatter':
                        
                        R1_rel_Diastolic =  (10**-6*self.lumpedValues['Diastolic'][vessel]['R1'][node]/133.32)/R1
                        C_rel_Diastolic =  133.32*10**6*self.lumpedValues['Diastolic'][vessel][key_C][node]/C
                        R1_rel_Systolic =  (10**-6*self.lumpedValues['Systolic'][vessel]['R1'][node]/133.32)/R1
                        C_rel_Systolic =  133.32*10**6*self.lumpedValues['Systolic'][vessel][key_C][node]/C
                        
                        R1_rel_vessel_phys = [R1_rel_Diastolic, R1_rel_Diastolic, R1_rel_Systolic, R1_rel_Systolic]
                        C_rel_vessel_phys = [C_rel_Diastolic, C_rel_Systolic, C_rel_Systolic, C_rel_Diastolic]
                    
                    
            
            if additionalPlot == 'ParamScatter':
                sourceDict_all = dict(x=R1_rel_all, y=C_rel_all, size=marker_all)
                sourceDict_vessel = dict(x=R1_rel_vessel, y=C_rel_vessel, size=marker_vessel)
                self.plotAdditional.xaxis.axis_label= 'R1/R1_opt'
                self.plotAdditional.yaxis.axis_label= 'C/C_opt'
                
                sourceDict_vessel_phys = dict(x=R1_rel_vessel_phys, y=C_rel_vessel_phys)
                
            elif additionalPlot == 'R1_div_R2_opt':
                sourceDict_all = dict(x=radius_rel_all, y=R1_div_R2_all, size=marker_all)
                sourceDict_vessel = dict(x=radius_rel_vessel, y=R1_div_R2_vessel, size=marker_vessel)
                self.plotAdditional.xaxis.axis_label= 'r/r_aorta'
                self.plotAdditional.yaxis.axis_label= 'R1/R2'
                
                sourceDict_vessel_phys = dict(x=[], y=[])
                
            elif additionalPlot == 'R1_div_R2_init':
                sourceDict_all = dict(x=radius_rel_all, y=R1_div_R2_0_all, size=marker_all)
                sourceDict_vessel = dict(x=radius_rel_vessel, y=R1_div_R2_0_vessel, size=marker_vessel)
                self.plotAdditional.xaxis.axis_label= 'r/r_aorta'
                self.plotAdditional.yaxis.axis_label= 'R1/R2'
                
                sourceDict_vessel_phys = dict(x=[], y=[])
                
            elif additionalPlot == 'R1_div_R2_init_div_opt':
                sourceDict_all = dict(x=radius_rel_all, y=R1_div_R2_init_div_opt_all, size=marker_all)
                sourceDict_vessel = dict(x=radius_rel_vessel, y=R1_div_R2_init_div_opt_vessel, size=marker_vessel)
                self.plotAdditional.xaxis.axis_label= 'r/r_aorta'
                self.plotAdditional.yaxis.axis_label= ''
                self.plotAdditional.title.text = ''
                
                sourceDict_vessel_phys = dict(x=[], y=[])
            
            
        
        elif additionalPlot == 'MeanValues':

            vesselColor = self.vesselId
            P_rel_all = []
            P_rel_vessel = []
            meanValueSum = 0
            varianceSum = 0
            for vessel in self.vesselIdList:
                
                r_end = self.lumpedValues[vessel]['radius'][node]
                radius_rel_all.append(r_end/r_end_aorta)
                #if vessel == 1:
                #print r_end, r_end/r_end_aorta
                if self.averageValueSelect.value == 'False':
                    P_lumped = self.lumpedValues[vessel]['Pressure'][0]
                elif self.averageValueSelect.value == 'True (Total Pressure)':
                    P_lumped = self.lumpedValues['MeanValues'][vessel]['Pressure'][0]
                else:
                    P_lumped = self.lumpedValues['MeanValuesStatic'][vessel]['Pressure'][0]
                P_1D = self.lumpedValues['1DValues'][vessel]['Pressure'][0]
                
                
                for n, relRadiusTmp in enumerate(relRadius):
                    if relRadiusTmp >= r_end/r_end_aorta:
                        n_use = n
                        break
                
                marker_all.append(markerSizeList[n])
                P_rel_all.append(P_lumped/P_1D)
                
                meanValueSum += 100*(P_lumped/P_1D)
                varianceSum += (100*(P_lumped/P_1D) - 100)**2
                
                if vessel == vesselColor:
                    radius_rel_vessel.append(r_end/r_end_aorta)
                    
                    P_rel_vessel.append(P_lumped/P_1D)
                    marker_vessel.append(markerSizeList[n])
            
            meanValue =  round(meanValueSum/len(self.vesselIdList), 3)
            standardDev = round(np.sqrt(varianceSum/len(self.vesselIdList)), 3)
            sourceDict_all = dict(x=radius_rel_all, y=P_rel_all, size=marker_all)
            sourceDict_vessel = dict(x=radius_rel_vessel, y=P_rel_vessel, size=marker_vessel)
            self.plotAdditional.xaxis.axis_label= 'r/r_aorta'
            self.plotAdditional.yaxis.axis_label= 'P/P_1D'
            self.plotAdditional.title.text = '({0}, {1})'.format(meanValue, standardDev)
            
            sourceDict_vessel_phys = dict(x=[], y=[])
                
        
        else:

            sourceDict_all = dict(x=[], y=[], size=[])
            sourceDict_vessel = dict(x=[], y=[], size=[])
            sourceDict_vessel_phys = dict(x=[], y=[])
            self.plotAdditional.xaxis.axis_label= ''
            self.plotAdditional.yaxis.axis_label= ''
            self.plotAdditional.title.text = ''
    
        return sourceDict_all, sourceDict_vessel, sourceDict_vessel_phys
    
        
        
    def getNetworkName(self):
        minimizationCase = self.minimizeSelect.value
        constraint = str(str(float(self.constraintSlider.value)))
        networkName = self.minimizationDicts[minimizationCase][constraint]['networkName']
        self.vesselId = self.minimizationDicts[minimizationCase]['vesselId']
        self.node = self.minimizationDicts[minimizationCase]['node']
            
        return networkName
    
            
            
        
    def loadData(self, networkName, dataNumber, vesselId, medical=True, relativeFilePath=False):
        
        
        file_name = networkName + "_SolutionData_" + dataNumber + ".hdf5"
        relFilePath = "app96Simple/data/Full96Model/ReducedNetworks/Final/" + networkName + "/SolutionData_" + dataNumber+ "/" + file_name
        if relativeFilePath:
            relFilePath = relFilePath#"dirAppDev/" +  "data/" + networkName + "/SolutionData_" + dataNumber+ "/" + networkName + "_SolutionData_" + dataNumber + ".hdf5"

        hdf5File = h5py.File(relFilePath,'r')
        time = hdf5File['VascularNetwork']['simulationTime'][:] - hdf5File['VascularNetwork']['simulationTime'][0]
        
        node = int(self.node)
        dt = 0.005 #time[1] - time[0]
        period = self.period
        vesseldict = {}

        N = int(round(period/dt))
        
        nCycles = int(round(time[-1]/period))
        
        
        t_start = period*(nCycles - 1)
        t_end = period*nCycles
        print t_start, t_end
        
        time_compare = np.linspace(t_start, t_end, N + 1)
        
        self.vesselIdList =  []
        Nv = 0
        for vesselName in hdf5File['vessels'].keys():
            
            tmpvesselId = vesselName.split(' - ')[-1]
            tmpvesselId = int(tmpvesselId)
            self.vesselIdList.append(tmpvesselId)
            if tmpvesselId == vesselId:
                
                P  = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Psol'][:, node], time_compare)
                Q  = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Qsol'][:, node], time_compare)
                A  = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Asol'][:, node], time_compare)
                
                Pf = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Psol_f'][:, node], time_compare)
                Pb = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Psol_b'][:, node], time_compare)
                Qf = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Qsol_f'][:, node], time_compare)
                Qb = self.mySplineInterpolater(time, hdf5File['vessels'][vesselName]['Qsol_b'][:, node], time_compare)
                
                Pm = self.findMean(time_compare, P)
                Qm = self.findMean(time_compare, Q)
                
                vName = vesselName.split(' - ')[0]
            Nv += 1
                
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
        vesseldict['time'] = time_compare
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
        vesseldict['Nv'] = Nv
        
        return vesseldict
    

    def mySplineInterpolater(self, x, y, x_new):
        tck = interpolate.splrep(x, y, s=0)

        y_new = interpolate.splev(x_new, tck, der=0)
        
        return y_new
        

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
    
    def removeOldSVG(self, staticDirRel="/app96Simple/static"):

        staticDirAbs = ''.join([self.cur, staticDirRel])

        for fileName in os.listdir(staticDirAbs):
            if "_reduced" in fileName:
                print "removing: ", fileName
                fileNamePath = ''.join([staticDirAbs, '/', fileName])
                os.remove(fileNamePath)
    
    
    def loadSVG(self):
        
        """ Method thats load the """

        old_file = "app96Simple/static/96model.svg"
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
        
        
    def changeSVG(self):
        """ Method that copy the lines in old_file which contain svg paths for each vessels.
            Every vessel which is in the reduced model is then represented with a path in a second layer, 
            with a different color.
        """

        old_file = "app96Simple/static/96model.svg"
        self.countImages += 1
        new_file = "app96Simple/static/96model_tmp" + str(timeModule.time()) + ".svg"
        previousFile = self.imageUrl

        if "_tmp" in previousFile:
            os.remove(previousFile)
        
        self.imageUrl = new_file
        f = open(old_file, "r")
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
                
            else:
                f2.write(oldLine)
                
        f.close()
        f2.close()
