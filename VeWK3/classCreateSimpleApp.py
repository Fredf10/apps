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
        N = 250
        self.cardiacCycles = 5
        self.t = np.linspace(0, self.T*self.cardiacCycles, N*self.cardiacCycles + 1)

        self.veWK3 = varyingElastance(self.t, self.T)
        
        t, P, P_LV, Q, E, V = self.veWK3.solveNCycle(Ncycles=self.cardiacCycles)
        t_last_cycle = np.linspace(self.T*(self.cardiacCycles - 1), self.T*(self.cardiacCycles - 0), 1001)
        P_last, P_LV_last, Q_last, E_last, V_last = self.getLastCycle(t_last_cycle, t, P, P_LV, Q, E, V)
        
        self.source_p = ColumnDataSource(data=dict(x=t, y=P))
        self.source_p_LV = ColumnDataSource(data=dict(x=t, y=P_LV))
        self.source_LV_loop = ColumnDataSource(data=dict(x=V, y=P_LV))
        
        self.source_p_last = ColumnDataSource(data=dict(x=t_last_cycle, y=P_last))
        self.source_p_LV_last = ColumnDataSource(data=dict(x=t_last_cycle, y=P_LV_last))
        self.source_q_last = ColumnDataSource(data=dict(x=t_last_cycle, y=Q_last))
        self.source_E_last = ColumnDataSource(data=dict(x=[], y=[]))
        # Set up plot_line y = a*x + b
        self.plot_P = Figure(plot_height=500, plot_width=650, title="aortic and ventricular pressure",
                             x_axis_label="t", y_axis_label="P [mmHg]",
                             tools="crosshair,pan,reset,save,wheel_zoom",
                             )
        self.plot_LV_loop = Figure(plot_height=500, plot_width=650, title="PV-loop",
                                   x_axis_label="V [ml]", y_axis_label="P [mmHg]",
                                   tools="crosshair,pan,reset,save,wheel_zoom",
                                   )

        self.plot_P_last = Figure(plot_height=500, plot_width=650, title="aortic and ventricular pressure (last cardiac cycle)",
                                  x_axis_label="t", y_axis_label="P [mmHg]",
                                  tools="crosshair,pan,reset,save,wheel_zoom",
                                  )

        self.plot_flow_or_elastance_last = Figure(plot_height=500, plot_width=650, title="flow",
                                                  x_axis_label="t", y_axis_label="flow [ml/s]",
                                                  tools="crosshair,pan,reset,save,wheel_zoom",
                                                  )
        
        self.plot_P.line('x', 'y', source=self.source_p, color='blue', line_alpha=0.6, line_width=2)
        self.plot_P.line('x', 'y', source=self.source_p_LV, color='green', line_alpha=0.6, line_width=2)

        self.plot_LV_loop.line('x', 'y', source=self.source_LV_loop, color='blue', line_alpha=0.6, line_width=2)
        
        self.plot_P_last.line('x', 'y', source=self.source_p_last, color='blue', line_alpha=0.6, line_width=2, legend="aorta")
        self.plot_P_last.line('x', 'y', source=self.source_p_LV_last, color='green', line_alpha=0.6, line_width=2, legend="left ventricle")
        
        self.plot_flow_or_elastance_last.line('x', 'y', source=self.source_q_last, color='blue', line_alpha=0.6, line_width=2)
        self.plot_flow_or_elastance_last.line('x', 'y', source=self.source_E_last, color='blue', line_alpha=0.6, line_width=2)
        
        self.resistanceSelect = Select(title="selcet total resistance", value="1.25", options=["1", "1.25", "1.5"])
        self.resistanceFactorSelect = Select(title="selcet factor for proximal resistance", value="0.1", options=["0.05", "0.1", "0.15", "0.2"])
        self.complianceSelect = Select(title="selcet total compliance", value="2.0", options=["1", "1.5", "2.0", "2.5"])
        self.eMaxSelect = Select(title="selcet E max", value="2.0", options=["1", "1.5", "2.0", "2.5"])
        self.eMinSelect = Select(title="selcet E min", value="0.06", options=["0.03", "0.06", "0.09", "0.12"])
        self.tPeakSelect = Select(title="selcet time to peak", value="0.32", options=["0.25", "0.28", "0.30", "0.32"])
        self.RvSelect = Select(title="selcet mitral resistance", value="0.005", options=["0.0025", "0.005", "0.01", "0.05"])
        self.n1Select = Select(title="select elastance shape-function n1", value="1.32", options=["1.1", "1.2", "1.32", "1.4"])
        self.n2Select = Select(title="select elastance shape-function n2", value="21.9", options=["15", "21.9", "25", "30"])
        self.flowOrElastanceSelect = Select(title="show flow or elastance", value="flow", options=["flow", "elastance"])
        self.nCyclesSelect = Select(title="select number of cycles", value="5", options=["2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15"])
        self.nTimePointsSelect = Select(title="select time-points per cycle", value="250", options=["100", "150", "200", "250", "500", "1000", "2000", "3000", "4000", "5000", "6000", "7000", "8000", "9000", "10000"])

        self.symbolicSelect = Select(title="use symbolic differentiation", value="True", options=["True", "False"])
        self.isoSelect = Select(title="use integrated iso eq", value="True", options=["True", "False"])

        


        self.Widgetlist = [self.resistanceSelect, self.resistanceFactorSelect, self.complianceSelect,
                           self.eMaxSelect, self.eMinSelect, self.tPeakSelect, self.RvSelect,
                           self.n1Select, self.n2Select, self.flowOrElastanceSelect, self.nCyclesSelect, 
                           self.nTimePointsSelect, self.symbolicSelect, self.isoSelect]
        
    def update_data(self, attrname, old, new):
    
        
        R_tot = float(self.resistanceSelect.value)

        R_factor = float(self.resistanceFactorSelect.value)
        C_tot = float(self.complianceSelect.value)
        Emax = float(self.eMaxSelect.value)

        Emin = float(self.eMinSelect.value)
        TPeak = float(self.tPeakSelect.value)
        Rv = float(self.RvSelect.value)
        n1 = float(self.n1Select.value)
        n2 = float(self.n2Select.value)
        flowOrElastance = self.flowOrElastanceSelect.value
        self.cardiacCycles = int(self.nCyclesSelect.value)
        symbolic_differentiation = (self.symbolicSelect.value)
        integrated_iso_eq = (self.isoSelect.value)
        N = int(self.nTimePointsSelect.value)
        t = np.linspace(0, self.T*self.cardiacCycles, N*self.cardiacCycles + 1)
        
        
        
        self.veWK3.initializeWKParams(R_tot=R_tot, C_tot=C_tot, R1_frac=R_factor)
        self.veWK3.Emax = Emax
        self.veWK3.Emin = Emin
        self.veWK3.TPeak = TPeak
        self.veWK3.Rv = Rv
        self.veWK3.n1 = n1
        self.veWK3.n2 = n2
        
        if symbolic_differentiation == "True":
            self.veWK3.symbolic_differentiation = True
        else:
            self.veWK3.symbolic_differentiation = False

        if integrated_iso_eq == "True":
            self.veWK3.integrated_iso_eq = True
        else:
            self.veWK3.integrated_iso_eq = False
        
        self.veWK3.t = t
        t, P, P_LV, Q, E, V = self.veWK3.solveNCycle(Ncycles=self.cardiacCycles)
        
        t_last_cycle = np.linspace(self.T*(self.cardiacCycles - 1), self.T*(self.cardiacCycles - 0), 1001)
        
        P_last, P_LV_last, Q_last, E_last, V_last = self.getLastCycle(t_last_cycle, t, P, P_LV, Q, E, V)
        
        self.source_p.data = dict(x=t, y=P)
        self.source_p_LV.data = dict(x=t, y=P_LV)
        self.source_LV_loop.data = dict(x=V_last, y=P_LV_last)
        self.source_p_last.data = dict(x=t_last_cycle, y=P_last)
        self.source_p_LV_last.data = dict(x=t_last_cycle, y=P_LV_last)
        
        if flowOrElastance == "flow":
            self.source_q_last.data = dict(x=t_last_cycle, y=Q_last)
            self.source_E_last.data = dict(x=[], y=[])
            self.plot_flow_or_elastance_last.title.text = "flow"
            self.plot_flow_or_elastance_last.yaxis.axis_label = "flow [ml/s]"
        else:
            self.source_q_last.data = dict(x=[], y=[])
            self.source_E_last.data = dict(x=t_last_cycle, y=E_last)
            self.plot_flow_or_elastance_last.title.text = "elastance"
            self.plot_flow_or_elastance_last.yaxis.axis_label = "E [mmHg/ml]"
            
        
        #self.plotQ.title.text = "Qm = ({0}, {1}); (reference, ecmo)".format(Qm, Qm_ecmo)
    
    def getLastCycle(self, t_last, t, P, P_LV, Q, E, V):
        
        P_last = np.interp(t_last, t, P)
        P_LV_last = np.interp(t_last, t, P_LV)
        Q_last = np.interp(t_last, t, Q)
        E_last = np.interp(t_last, t, E)
        V_last = np.interp(t_last, t, V)
        
        return P_last, P_LV_last, Q_last, E_last, V_last
        

