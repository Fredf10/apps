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
        mm_to_m = 1e-3
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
        self.solutionDict, self.solutionDict_space, self.solutionDict_mu, self.solutionDict_space_mu, self.solutionDict_geom, self.solutionDict_space_geom, self.solutionDict_geom_mu, self.solutionDict_space_geom_mu = self.load_data()
        #===================================================
        # create sample-points along longitudinal (z) axis
        #===================================================
        #self.samplepointsLongitudinal = np.arange(-20, 71, 0.1)#np.arange(-20, 70, 1.)
        #self.samplepointsLongitudinal = np.arange(-100, 150.1, 0.1)
        self.samplepointsLongitudinal = np.arange(-20, 140.1, 0.1)
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
        
        U0 = max(uz)/2.
        R0 = max(r)*mm_to_m
        self.Q = U0*np.pi*R0**2
        
        indx = np.where(self.samplepointsLongitudinal == xPos_init)
        xPos = X[indx]
        rPos = R[indx]
        pPos = P[indx]
        
        self.source_space = ColumnDataSource(data=dict(x=X, y=R, P=P, U=U, P_tot=P_tot,))
        self.source_space_conv = ColumnDataSource(data=dict(x=[], P_convective=[], P_convective_flat=[],))
        self.source_space_visc = ColumnDataSource(data=dict(x=[], P_visc=[], P_visc_pouseille=[],))
        self.source_space_oneD = ColumnDataSource(data=dict(x=[], P_1D=[], P_1D_conventional=[],))
        self.source = ColumnDataSource(data=dict(x=r, y=uz, P=p,))
        self.source_pouseille = ColumnDataSource(data=dict(x=[], y=[],))
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
        self.plot_velocity_pos = Figure(plot_height=450, plot_width=600, title="velocity profile at given position",
                                        x_axis_label="Uz [m/s]", y_axis_label="r [mm]",
                                        tools="crosshair,pan,reset,resize,save,wheel_zoom",
                                        x_range=[-0.2, 5.0], y_range=[0, 1.6])
        
        self.plot_velocity_pos.line('y', 'x', source=self.source, color='black', line_alpha=1, line_width=2)
        self.plot_velocity_pos.line('y', 'x', source=self.source_pouseille, color='chartreuse', line_alpha=1, line_width=2)

        self.plot_pressure_pos = Figure(plot_height=450, plot_width=600, title="pressure profile at given position",
                                        x_axis_label="P [mmHg]", y_axis_label="r [mm]" ,
                                        tools="crosshair,pan,reset,resize,save,wheel_zoom",
                                        x_range=[np.min(P)*1.9, np.max(P)*1.1], y_range=[0, 1.6])
        
        self.plot_pressure_pos.line('P', 'x', source=self.source, color='black', line_alpha=1, line_width=2)
        
         

        self.positionSlider = Slider(title="position", value=-20, start=-20, end=140, step=0.1)
        self.viscositySelect = Select(title="select viscosity mP", value="3.5", options=["3.5", "1.75"])
        self.stenosisLengthSelect = Select(title="select stenosis length [cm]", value="1", options=["1", "2"])
        self.viscousSelect = Select(title="show viscous part", value="False", options=["True", "False"])
        self.velocityComponentSelect = Select(title="select velocity component", value="Uz", options=["Uz", "Ur"])
        #self.viscousPouseilleSelect = Select(title="show viscous part (pouseille)", value="False", options=["True", "False"])
        self.convectiveSelect = Select(title="show convective part", value="False", options=["True", "False"])
        #self.convectiveFlatSelect = Select(title="show convective part (flat)", value="False", options=["True", "False"])
        self.oneDSelect = Select(title="show 1D", value="False", options=["True", "False"])
        #self.oneDconventionalSelect = Select(title="show 1D", value="False", options=["True", "False"])
        
        #self.lineSelect = Select(title="selcet line", value="linear", options=["linear", "power"])
        #self.imageSelect = Select(title="select image", value="red", options=["red", "black"])


        self.Widgetlist = [self.positionSlider, self.viscositySelect, self.stenosisLengthSelect, self.velocityComponentSelect, self.viscousSelect, self.convectiveSelect, self.oneDSelect]
        
    def update_data(self, attrname, old, new):

        #===================================================
        # conversion factors
        #===================================================
        m_to_mm = 1e3
        mm_to_m = 1e-3
        Pa_to_mmHg = 1./133.32
        
        
        xPos = self.positionSlider.value
        viscosity = self.viscositySelect.value
        stenosisLength = self.stenosisLengthSelect.value
        velocityComponent = self.velocityComponentSelect.value
        
        viscous = self.viscousSelect.value
        convective = self.convectiveSelect.value
        oneD = self.oneDSelect.value
        
        if viscosity == "3.5":
            if stenosisLength == "1":
                solutionDict = self.solutionDict
                solutionDict_space = self.solutionDict_space
            else:
                solutionDict = self.solutionDict_geom
                solutionDict_space = self.solutionDict_space_geom
            mu = 3.5*1e-3
        else:
            if stenosisLength == "1":
                solutionDict = self.solutionDict_mu
                solutionDict_space = self.solutionDict_space_mu
            else:
                solutionDict = self.solutionDict_geom_mu
                solutionDict_space = self.solutionDict_space_geom_mu
            mu = 0.5*3.5*1e-3
        #===================================================
        # retrieve datoa along the z-axis
        #===================================================
        X = solutionDict_space['x']*m_to_mm
        R = solutionDict_space['r']*m_to_mm
        P = solutionDict_space['P']*Pa_to_mmHg
        P_tot = solutionDict_space['P_tot']*Pa_to_mmHg
        U = solutionDict_space['U']
        
        P_convective = P[0] - solutionDict_space['deltaP_convective']*Pa_to_mmHg
        P_convective_flat = P[0] - solutionDict_space['deltaP_convective_flat']*Pa_to_mmHg
        P_viscous = P[0] - solutionDict_space['deltaP_viscous']*Pa_to_mmHg
        P_viscous_pouseille = P[0] - solutionDict_space['deltaP_viscous_pouseille']*Pa_to_mmHg
        
        P_oneD = P_convective + P_viscous - P[0]
        P_oneD_conventional = P_convective_flat + P_viscous_pouseille - P[0]

        #===================================================
        # retrieve data at a certain position along the z axis
        #===================================================
        r = solutionDict[xPos]['r']*m_to_mm
        uz = solutionDict[xPos]['uz']
        ur = solutionDict[xPos]['ur']
        p = solutionDict[xPos]['P']*Pa_to_mmHg
        
        indx = np.where(np.abs(self.samplepointsLongitudinal - xPos) < 0.0001)
        
        xPos = X[indx]
        rPos = R[indx]
        pPos = P[indx]
        uPos = U[indx]
        rho = 1050

        Re = rho*uPos*2*(rPos*mm_to_m)/mu

        #self.source = ColumnDataSource(data=dict(x=r, y=uz, P=p,))
        #self.source_pos = ColumnDataSource(data=dict(x=xPos, r=rPos, P=pPos,))
        if velocityComponent == "Uz":
            self.source.data = dict(x=r, y=uz, P=p)
            self.plot_velocity_pos.xaxis.axis_label = "Uz [m/s]"
            self.plot_velocity_pos.title.text = "velocity profile at given position, Re: {0}".format(round(Re, 2))

        else:
            self.source.data = dict(x=r, y=ur, P=p)
            self.plot_velocity_pos.xaxis.axis_label = "Ur [m/s]"
            self.plot_velocity_pos.title.text = "velocity profile at given position, Re: {0}".format(round(Re, 2))

        self.source_space.data = dict(x=X, y=R, P=P, U=U, P_tot=P_tot,)
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
            
            if velocityComponent == "Uz":
                
                r0 = max(r)/1000.
                U0 = self.Q/(np.pi*r0**2)
                U_pouseille = self.velocity(U0, r0*1000, r, 2)
                
                self.source_pouseille.data = dict(x=r, y=U_pouseille)
            else:
                self.source_pouseille.data = dict(x=[], y=[])
            
            
        else:
            self.source_space_oneD.data = dict(x=[], P_1D=[], P_1D_conventional=[])
            self.source_pouseille.data = dict(x=[], y=[])
        
        
        
    def load_data(self, relFilePath=False):

        if relFilePath:
            solutionDict_path = join('stenosis', 'data', 'solutionDict.p')
            solutionDict_space_path = join('stenosis', 'data', 'solutionDict_space.p')
            solutionDict_path_mu = join('stenosis', 'data', 'solutionDict_mu.p')
            solutionDict_space_path_mu = join('stenosis', 'data', 'solutionDict_space_mu.p')
            solutionDict_path_geom_mu = join('stenosis', 'data', 'solutionDict_geom_mu.p')
            solutionDict_space_path_geom_mu = join('stenosis', 'data', 'solutionDict_space_geom_mu.p')
            solutionDict_path_geom = join('stenosis', 'data', 'solutionDict_geom.p')
            solutionDict_space_path_geom = join('stenosis', 'data', 'solutionDict_space_geom.p')
            
        else:
            solutionDict_path = join(dirname(__file__), 'data', 'solutionDict.p')
            solutionDict_space_path = join(dirname(__file__), 'data', 'solutionDict_space.p')
            solutionDict_path_mu = join(dirname(__file__), 'data', 'solutionDict_mu.p')
            solutionDict_space_path_mu = join(dirname(__file__), 'data', 'solutionDict_space_mu.p')
            solutionDict_path_geom_mu = join(dirname(__file__), 'data', 'solutionDict_geom_mu.p')
            solutionDict_space_path_geom_mu = join(dirname(__file__), 'data', 'solutionDict_space_geom_mu.p')
            solutionDict_path_geom = join(dirname(__file__), 'data', 'solutionDict_geom.p')
            solutionDict_space_path_geom = join(dirname(__file__), 'data', 'solutionDict_space_geom.p')
        solutionDict = pickle.load(open(solutionDict_path, "rb"))
        solutionDict_space = pickle.load(open(solutionDict_space_path, "rb"))
        solutionDict_mu = pickle.load(open(solutionDict_path_mu, "rb"))
        solutionDict_space_mu = pickle.load(open(solutionDict_space_path_mu, "rb"))
        solutionDict_geom_mu = pickle.load(open(solutionDict_path_geom_mu, "rb"))
        solutionDict_space_geom_mu = pickle.load(open(solutionDict_space_path_geom_mu, "rb"))
        solutionDict_geom = pickle.load(open(solutionDict_path_geom, "rb"))
        solutionDict_space_geom = pickle.load(open(solutionDict_space_path_geom, "rb"))
        
        return solutionDict, solutionDict_space, solutionDict_mu, solutionDict_space_mu, solutionDict_geom, solutionDict_space_geom, solutionDict_geom_mu, solutionDict_space_geom_mu
    
    def velocity(self, U, R, r, zetha):
        zetha=float(zetha)
        
        u=U*((zetha+2)/zetha)*(1-(r/R)**zetha)
        return u