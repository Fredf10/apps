'''
Created on Oct 5, 2016

@author: fredrik
'''
'''
Created on Oct 3, 2016

@author: fredrik
'''

import numpy as np
import pickle
from os.path import dirname, join
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.widgets import Select, Slider



class SetupApp:
    
    def __init__(self):

        #===================================================
        # conversion factors
        #===================================================
        m_to_mm = 1e3
        Pa_to_mmHg = 1./133.32
        #===================================================
        # load dictionaries with data format:
        # solutionDict_time[sample] = {'r':np.array(r_x),    [m]
        #                              'uz':np.array(U[2]),  [m/s]
        #                              'ur':np.array(U[0]),  [m/s]
        #                              'P':np.array(P)}      [Pa]
        #     
        # solutionDict_space = {'x':np.array(X),                                     [m]
        #                       'r':np.array(r_x),                                   [m]
        #                       'U':np.array(U_mean),                                [m/s]
        #                       'P':np.array(P_mean),                                [Pa]
        #                       'P_tot':np.array(P_total)}                           [Pa]
        #                       'deltaP_convective':deltaP_convective,               [Pa]
        #                       'deltaP_convective_flat':deltaP_convective_flat,     [Pa]
        #                       'deltaP_viscous':deltaP_viscous,                     [Pa]
        #                        'deltaP_viscous_pouseille':deltaP_viscous_pouseille [Pa]
        #===================================================

        self.solutionDict, self.solutionDict_space = self.load_data()
        #===================================================
        # create sample-points along longitudinal (z) axis
        #===================================================
        #self.samplepointsLongitudinal = np.arange(-20, 71, 0.1)#np.arange(-20, 70, 1.)
        #self.samplepointsLongitudinal = np.arange(-100, 150.1, 0.1)
        self.samplepointsLongitudinal = np.arange(-20, 120.1, 0.1)
        #===================================================
        # retrieve datoa along the z-axis
        #===================================================
        X = self.solutionDict_space['x']*m_to_mm
        R = self.solutionDict_space['r']*m_to_mm
        P = self.solutionDict_space['P']*Pa_to_mmHg
        P_tot = self.solutionDict_space['P_tot']*Pa_to_mmHg
        U = self.solutionDict_space['U']

        #===================================================
        # retrieve data at a certain position along the z axis
        #===================================================
        xPos_init = -20.
        r = self.solutionDict[xPos_init]['r']*m_to_mm
        uz = self.solutionDict[xPos_init]['uz']
        ur = self.solutionDict[xPos_init]['ur']
        p = self.solutionDict[xPos_init]['P']*Pa_to_mmHg
        
        indx = np.where(self.samplepointsLongitudinal == xPos_init)
        xPos = X[indx]
        rPos = R[indx]
        pPos = P[indx]
        
        self.source_space = ColumnDataSource(data=dict(x=X, y=R, P=P, U=U, P_tot=P_tot,))
        self.source_space_conv = ColumnDataSource(data=dict(x=[], P_convective=[], P_convective_flat=[],))
        self.source_space_visc = ColumnDataSource(data=dict(x=[], P_visc=[], P_visc_pouseille=[],))
        self.source_space_oneD = ColumnDataSource(data=dict(x=[], P_1D=[], P_1D_conventional=[],))
        self.source = ColumnDataSource(data=dict(x=r, y=uz, P=p,))
        self.source_pos = ColumnDataSource(data=dict(x=xPos, r=rPos, P=pPos,))
        # Set up plot_line y = a*x + b
        self.plot_radius = Figure(plot_height=450, plot_width=600, title="radius",
                                  x_axis_label="z [mm]", y_axis_label="r [mm]",
                                  tools="crosshair,pan,reset,resize,save,wheel_zoom",
                                  y_range=[0, 1.6])
        self.plot_radius.line('x', 'y', source=self.source_space, color='green', line_alpha=0.6, line_width=2)
        self.plot_radius.circle(x="x", y="r", source=self.source_pos, size=12, color='green', fill_alpha=1)
        
        self.plot_pressure = Figure(plot_height=450, plot_width=600, title="pressure",
                                    x_axis_label="z [mm]", y_axis_label="P [mmHg]",
                                    tools="crosshair,pan,reset,resize,save,wheel_zoom")
        
        self.plot_pressure.line('x', 'P', source=self.source_space, color='black', line_alpha=1, line_width=2, legend="Pressure")
        self.plot_pressure.line('x', 'P_tot', source=self.source_space, color='grey', line_alpha=0.5, line_width=2, legend="total Pressure")
        self.plot_pressure.line('x', 'P_convective', source=self.source_space_conv, color='blue', line_alpha=1, line_width=2, legend="1D conv.")
        self.plot_pressure.line('x', 'P_convective_flat', source=self.source_space_conv, color='green', line_alpha=1, line_width=2, legend="1D conv. f.")
        self.plot_pressure.line('x', 'P_visc', source=self.source_space_visc, color='cyan', line_alpha=1, line_width=2, legend="1D visc")
        self.plot_pressure.line('x', 'P_visc_pouseille', source=self.source_space_visc, color='red', line_alpha=1, line_width=2, legend="1D visc p.")
        self.plot_pressure.line('x', 'P_1D', source=self.source_space_oneD, color='darkorange', line_alpha=1, line_width=2, legend="1D")
        self.plot_pressure.line('x', 'P_1D_conventional', source=self.source_space_oneD, color='chartreuse', line_alpha=1, line_width=2, legend="1D f. p.")

        self.plot_pressure.circle(x="x", y="P", source=self.source_pos, size=12, color='black', fill_alpha=1)
        self.plot_pressure.legend.background_fill_alpha = 0.
        self.plot_velocity_pos = Figure(plot_height=450, plot_width=600, title="radial velocity profile at given position",
                                        x_axis_label="Uz [m/s]", y_axis_label="r [mm]",
                                        tools="crosshair,pan,reset,resize,save,wheel_zoom",
                                        x_range=[-0.2, 3.0], y_range=[0, 1.6])
        
        self.plot_velocity_pos.line('y', 'x', source=self.source, color='black', line_alpha=1, line_width=2)

        self.plot_pressure_pos = Figure(plot_height=450, plot_width=600, title="radial pressure profile at given position",
                                        x_axis_label="P [mmHg]", y_axis_label="r [mm]" ,
                                        tools="crosshair,pan,reset,resize,save,wheel_zoom",
                                        x_range=[np.min(P)*1.9, np.max(P)*1.1], y_range=[0, 1.6])
        
        self.plot_pressure_pos.line('P', 'x', source=self.source, color='black', line_alpha=1, line_width=2)
        
         

        self.positionSlider = Slider(title="position", value=-20, start=-20, end=120, step=0.1)
        self.viscousSelect = Select(title="show viscous part", value="False", options=["True", "False"])
        #self.viscousPouseilleSelect = Select(title="show viscous part (pouseille)", value="False", options=["True", "False"])
        self.convectiveSelect = Select(title="show convective part", value="False", options=["True", "False"])
        #self.convectiveFlatSelect = Select(title="show convective part (flat)", value="False", options=["True", "False"])
        self.oneDSelect = Select(title="show 1D", value="False", options=["True", "False"])
        #self.oneDconventionalSelect = Select(title="show 1D", value="False", options=["True", "False"])
        
        #self.lineSelect = Select(title="selcet line", value="linear", options=["linear", "power"])
        #self.imageSelect = Select(title="select image", value="red", options=["red", "black"])


        self.Widgetlist = [self.positionSlider, self.viscousSelect, self.convectiveSelect, self.oneDSelect]
        
    def update_data(self, attrname, old, new):

        #===================================================
        # conversion factors
        #===================================================
        m_to_mm = 1e3
        Pa_to_mmHg = 1./133.32
        
        xPos = self.positionSlider.value
        
        viscous = self.viscousSelect.value
        convective = self.convectiveSelect.value
        oneD = self.oneDSelect.value

        #===================================================
        # retrieve datoa along the z-axis
        #===================================================
        X = self.solutionDict_space['x']*m_to_mm
        R = self.solutionDict_space['r']*m_to_mm
        P = self.solutionDict_space['P']*Pa_to_mmHg
        P_tot = self.solutionDict_space['P_tot']*Pa_to_mmHg
        U = self.solutionDict_space['U']
        
        P_convective = P[0] - self.solutionDict_space['deltaP_convective']*Pa_to_mmHg
        P_convective_flat = P[0] - self.solutionDict_space['deltaP_convective_flat']*Pa_to_mmHg
        P_viscous = P[0] - self.solutionDict_space['deltaP_viscous']*Pa_to_mmHg
        P_viscous_pouseille = P[0] - self.solutionDict_space['deltaP_viscous_pouseille']*Pa_to_mmHg
        
        P_oneD = P_convective + P_viscous - P[0]
        P_oneD_conventional = P_convective_flat + P_viscous_pouseille - P[0]

        #===================================================
        # retrieve data at a certain position along the z axis
        #===================================================
        r = self.solutionDict[xPos]['r']*m_to_mm
        uz = self.solutionDict[xPos]['uz']
        ur = self.solutionDict[xPos]['ur']
        p = self.solutionDict[xPos]['P']*Pa_to_mmHg
        
        indx = np.where(np.abs(self.samplepointsLongitudinal - xPos) < 0.0001)
        
        xPos = X[indx]
        rPos = R[indx]
        pPos = P[indx]

        #self.source = ColumnDataSource(data=dict(x=r, y=uz, P=p,))
        #self.source_pos = ColumnDataSource(data=dict(x=xPos, r=rPos, P=pPos,))
        
        self.source.data = dict(x=r, y=uz, P=p)
        self.source_pos.data = dict(x=xPos, r=rPos, P=pPos,)
        
        if viscous == "True":
            self.source_space_visc.data = dict(x=X, P_visc=P_viscous, P_visc_pouseille=P_viscous_pouseille)
        else:
            self.source_space_visc.data = dict(x=[], P_visc=[], P_visc_pouseille=[])
        if convective == "True":
            self.source_space_conv.data = dict(x=X, P_convective=P_convective, P_convective_flat=P_convective_flat)

        else:
            self.source_space_conv.data = dict(x=[], P_convective=[], P_convective_flat=[])
        if oneD == "True":
            self.source_space_oneD.data = dict(x=X, P_1D=P_oneD, P_1D_conventional=P_oneD_conventional)
        else:
            self.source_space_oneD.data = dict(x=[], P_1D=[], P_1D_conventional=[])
        
        
        
    def load_data(self, relFilePath=False):

        if relFilePath:
            solutionDict_path = join('stenosis', 'data', 'solutionDict.p')
            solutionDict_space_path = join('stenosis', 'data', 'solutionDict_space.p')
            
        else:
            solutionDict_path = join(dirname(__file__), 'data', 'solutionDict.p')
            solutionDict_space_path = join(dirname(__file__), 'data', 'solutionDict_space.p')
        solutionDict = pickle.load(open(solutionDict_path, "rb"))
        solutionDict_space = pickle.load(open(solutionDict_space_path, "rb"))
        
        return solutionDict, solutionDict_space
