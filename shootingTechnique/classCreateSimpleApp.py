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


def func(y, x):
    y0 = y[0]
    y1 = y[1]
    
    return np.array([y1, -2*y0*y1])
# define Euler solver
def euler(func, z0, time):
    """The Euler scheme for solution of systems of ODEs.
    z0 is a vector for the initial conditions,
    the right hand side of the system is represented by func which returns
    a vector with the same size as z0 ."""

    z = np.zeros((np.size(time), np.size(z0)))
    z[0,:] = z0

    for i in range(len(time)-1):
        dt = time[i+1] - time[i]
        z[i+1,:]=z[i,:] + np.asarray(func(z[i,:], time[i]))*dt

    return z

class SetupApp:
    
    def __init__(self):


        x_0, x_end = 1, 2
        y_0, y_end = 1, 0.5
        N = 3
        
        X = np.linspace(x_0, x_end, N +1)
        Y_analytic = 1/X
        
        dx = X[1] - X[0]
        s0 = 0
        
        x_s = np.linspace(x_0, 1.2, N +1)
        s_line = y_0 + s0*(x_s - x_0)
        
        Y_0 =[y_0, s0]
        self.source_bvp = ColumnDataSource(data=dict(x=[x_0, x_end], y=[y_0, y_end]))

        self.source_sol = ColumnDataSource(data=dict(X=X, Y_analytic=Y_analytic))
        self.source_sol_euler = ColumnDataSource(data=dict(x=[], y=[]))
        self.source_s = ColumnDataSource(data=dict(x=[], y=[]))

        self.source_phi = ColumnDataSource(data=dict(x=[], y=[]))
        self.source_phi_line = ColumnDataSource(data=dict(x=[], y=[]))
        self.source_phi_line_sol = ColumnDataSource(data=dict(x=[], y=[]))
        
        self.plot_sol = Figure(plot_height=450, plot_width=600,
                               x_axis_label="x ", y_axis_label="y",
                               tools="crosshair,pan,reset,resize,save,wheel_zoom")
    
        self.plot_phi = Figure(plot_height=450, plot_width=600,
                               x_axis_label="s ", y_axis_label="phi",
                               tools="crosshair,pan,reset,resize,save,wheel_zoom",
                               x_range=[-2, 1], y_range=[-1, 1])
        
        
        #self.plot_sol.line(x='X', y='Y_analytic', source=self.source_sol, color='black', line_alpha=1, line_width=2, line_dash="dashed",legend="analytical")
        self.plot_sol.circle(x="x", y="y", source=self.source_bvp, size=12, color='black', fill_alpha=1)
        self.plot_phi.circle(x="x", y="y", source=self.source_phi, size=12, color='blue', fill_alpha=0.5)
        
        self.plot_phi.line(x="x", y="y", source=self.source_phi_line, color='grey', line_alpha=1, line_width=1)
       
        
        self.plot_sol.line(x="x", y="y", source=self.source_sol_euler, color='blue', line_alpha=1, line_width=2, legend="y(x;s)")
        self.plot_sol.line(x='x', y='y', source=self.source_s, color='red', line_alpha=1, line_width=2, legend="s")
        self.plot_sol.line(x="x", y="y", source=self.source_phi_line_sol, color='grey', line_alpha=1, line_width=1)
        self.plot_sol.legend.location = "bottom_left"
        self.stepSlider = Slider(title="step", value=0, start=0, end=17, step=1)
        
        self.Widgetlist = [self.stepSlider]
        
        self.createLists()

        listOfPlots = [self.plot_sol, self.plot_phi]
        
        self.ajustFontSize(listOfPlots, '20pt', '15pt')
        
    
    def ajustFontSize(self, listOfObjects, font_size, axis_size):
    
        
        for ob in listOfObjects:
            ob.xaxis.major_label_text_font_size = axis_size
            ob.yaxis.major_label_text_font_size = axis_size
            
            ob.xaxis.axis_label_text_font_size = font_size
            ob.yaxis.axis_label_text_font_size = font_size
        
        
    def update_data(self, attrname, old, new):

        
        step = self.stepSlider.value
        print step
        if step == 1:
            self.source_s.data = dict(x=self.x_s, y=self.s_vector_list[0])
        if step == 2:
            self.source_sol_euler.data = dict(x=self.X, y=self.y_list[0])
            
        elif step == 3:
            self.source_phi.data = dict(x=self.s_list[0:1], y=self.phi_list[0:1])
            self.source_phi_line_sol.data = dict(x=self.x_end_array, y=self.phi_vector_list[0])
        elif step == 4:
            self.source_phi_line_sol.data = dict(x=[], y=[])
            self.source_s.data = dict(x=self.x_s, y=self.s_vector_list[1])
        elif step == 5:
            self.source_sol_euler.data = dict(x=self.X, y=self.y_list[1])

        elif step == 6:
            self.source_phi.data = dict(x=self.s_list[0:2], y=self.phi_list[0:2])
            self.source_phi_line_sol.data = dict(x=self.x_end_array, y=self.phi_vector_list[1])
        elif step == 7:
            self.source_phi_line.data = dict(x=self.s_list[0:3], y=self.phi_list[0:3])
            self.source_phi_line_sol.data = dict(x=[], y=[])
        elif step == 8:
            self.source_s.data = dict(x=self.x_s, y=self.s_vector_list[2])

        elif step == 9:
            self.source_sol_euler.data = dict(x=self.X, y=self.y_list[2])

        elif step == 10:
            self.source_phi.data = dict(x=self.s_list[0:3], y=self.phi_list[0:3])
            self.source_phi_line_sol.data = dict(x=self.x_end_array, y=self.phi_vector_list[2])
        elif step == 11:
            self.source_phi_line.data = dict(x=self.s_list[1:4], y=self.phi_list[1:4])
            self.source_phi_line_sol.data = dict(x=[], y=[])
        elif step == 12:
            self.source_s.data = dict(x=self.x_s, y=self.s_vector_list[3])

        elif step == 13:
            self.source_sol_euler.data = dict(x=self.X, y=self.y_list[3])

        elif step == 14:
            self.source_phi.data = dict(x=self.s_list[0:4], y=self.phi_list[0:4])
            self.source_phi_line_sol.data = dict(x=self.x_end_array, y=self.phi_vector_list[3])
        elif step == 15:
            self.source_phi_line.data = dict(x=self.s_list[2:5], y=self.phi_list[2:5])
            self.source_phi_line_sol.data = dict(x=[], y=[])
        elif step == 16:
            self.source_s.data = dict(x=self.x_s, y=self.s_vector_list[4])

        elif step == 17:
            self.source_sol_euler.data = dict(x=self.X, y=self.y_list[4])


            
    
    def createLists(self):
        
        s_list = []
        s_vector_list = []
        y_list = []
        phi_list = []
        phi_vector_list = []
        x_0, x_end = 1, 2
        y_0, y_end = 1, 0.5
        N = 3
        
        X = np.linspace(x_0, x_end, N +1)
        
        dx = X[1] - X[0]
        s0 = 0
        x_s = np.linspace(x_0, 1.2, N +1)
        s0_line = y_0 + s0*(x_s - x_0)
        Y_0 =[y_0, s0]
        
        Y_num = euler(func, Y_0, X)
        Y0 = Y_num[:, 0]
        Y1 = Y_num[:, 1]
        
        phi0 = Y0[-1] - y_end
        
        s_list.append(s0)
        s_vector_list.append(s0_line)
        y_list.append(Y0)
        phi_list.append(phi0)
        phi_vector_list.append([Y0[-1], y_end])
        
        print "n: ", 0, "s: ", s0, "phi: ", phi0, "y(2): ", Y0[-1], "y'(2): ", Y1[-1]
        s1 = -0.5
        
        itMax = 10
        tol = 10**(-10)
        for n in range(itMax):
            Y_0 =[y_0, s1]
        
            Y_num = euler(func, Y_0, X)
            Y0 = Y_num[:, 0]
            Y1 = Y_num[:, 1]
        
            phi1 = Y0[-1] - y_end

            
            s1_line = y_0 + s1*(x_s - x_0)
            s_list.append(s1)
            s_vector_list.append(s1_line)
            y_list.append(Y0)
            phi_list.append(phi1)
            phi_vector_list.append([Y0[-1], y_end])
            #print n + 1, s1, phi1, Y0[-1], Y1[-1]
            print "n: ", n + 1, "s: ", s1, "phi: ", phi1, "y(2): ", Y0[-1], "y'(2): ", Y1[-1]
            s2 = s1 - phi1*(s1 - s0)/(phi1 - phi0)
            
            s0, s1, phi0 = s1, s2, phi1
        
            if abs(s1 - s0) < tol:
                break
            
    

        self.s_list = s_list
        self.s_vector_list = s_vector_list
        self.y_list = y_list
        self.phi_list = phi_list
        self.phi_vector_list = phi_vector_list
        
        self.X = X
        self.x_end_array = [x_end, x_end]
        self.x_s = x_s




        
