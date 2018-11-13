'''
Created on Oct 5, 2016

@author: fredrik
'''
'''
Created on Oct 3, 2016

@author: fredrik
'''

import numpy as np
import os
from os.path import dirname, join
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, Range1d
from bokeh.models.widgets import Select
from classVeWK3 import varyingElastance

import h5py


class SetupApp:
    
    def __init__(self):
        
        
        self.T = 1.
        N = 1000*5
        self.cardiacCycles = 5
        self.t = np.linspace(0, self.T*self.cardiacCycles, N*self.cardiacCycles + 1)

        self.veWK3 = varyingElastance(self.t, self.T)
        
        t, P, P_LV, Q, E, V = self.veWK3.solveNCycle(Ncycles=self.cardiacCycles)
        
        self.source_p = ColumnDataSource(data=dict(x=t, y=P))
        self.source_p_LV = ColumnDataSource(data=dict(x=t, y=P_LV))
        # Set up plot_line y = a*x + b
        self.plot_line = Figure(plot_height=600, plot_width=800, title="aortic Pressure",
                                x_axis_label="t", y_axis_label="P [mmHg]",
                                tools="crosshair,pan,reset,save,wheel_zoom",
                                )
        
        self.plot_line.line('x', 'y', source=self.source_p, color='blue', line_alpha=0.6, line_width=2)
        self.plot_line.line('x', 'y', source=self.source_p_LV, color='green', line_alpha=0.6, line_width=2)

        self.resistanceSelect = Select(title="selcet total resistance", value="1.25", options=["1", "1.25", "1.5"])
        self.resistanceFactorSelect = Select(title="selcet factor for proximal resistance", value="0.1", options=["0.05", "0.1", "0.15", "0.2"])
        self.complianceSelect = Select(title="selcet total compliance", value="2.0", options=["1", "1.5", "2.0", "2.5"])
        self.eMaxSelect = Select(title="selcet E max", value="2.0", options=["1", "1.5", "2.0", "2.5"])
        


        self.Widgetlist = [self.resistanceSelect, self.resistanceFactorSelect, self.complianceSelect,
                           self.eMaxSelect]
        
    def update_data(self, attrname, old, new):
    
        
        R_tot = float(self.resistanceSelect.value)

        R_factor = float(self.resistanceFactorSelect.value)
        C_tot = float(self.complianceSelect.value)
        Emax = float(self.eMaxSelect.value)
        
        self.veWK3.initializeWKParams(R_tot=R_tot, C_tot=C_tot, R1_frac=R_factor)
        self.veWK3.Emax = Emax
        
        t, P, P_LV, Q, E, V = self.veWK3.solveNCycle(Ncycles=self.cardiacCycles)
        
        self.source_p.data = dict(x=t, y=P)
        self.source_p_LV.data = dict(x=t, y=P_LV)
        

