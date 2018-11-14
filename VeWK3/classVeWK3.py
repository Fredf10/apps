'''
Created on Apr 12, 2018

@author: fredrik
'''

import numpy as np
#import matplotlib.pylab as plt
#import scipy.optimize
from ODEschemes import rk4


class varyingElastance:
    
    def __init__(self, t, T):
        
        self.t = t
        self.T = T
        self.TPeak = 0.28
        self.Emax = 1.6
        self.Emin = 0.06
        self.Patrium = 7.5
        self.V0 = 20.
        self.Vinit = 140.
        self.Rv = 0.005
        self.Pao_0 = 80
        ## Shape parameters
        self.alpha = 1.672
        self.n1 = 1.32
        self.n2 = 21.9
        self.T_shift = 0.0
        
        self.a1 = 0.708
        self.a2Factor = 1.677
        
        self.it = 0
        
        self.initializeWKParams()
        
    def initializeWKParams(self, R_tot=1.26, C_tot=2.0, Pv=7.5, R1_frac=0.1):
        self.Rtot = R_tot
        self.R1 = R1_frac*R_tot
        self.R2 = R_tot - self.R1
        self.C = C_tot
        self.Pv = Pv

    def Elastance(self, tin):
        """Computes the value of the elastance at time t, according to the shape parameters given by Stergiopolus and scaled
        according to Tpeak, T, Emax and Emin. """
        t = tin - int(tin/self.T)*self.T
        t = t - self.T_shift
        if t < 0:
            t = self.T + t
            
            
        alpha_inv = ((1./self.a1)**self.n1)/(1 + (1./self.a1)**self.n1)*(1./(1 + (1./(self.a1*self.a2Factor))**self.n2))
        alpha = 1/alpha_inv
        a1 = self.a1 * self.TPeak
        a2 = self.a2Factor * a1
        n1, n2 = self.n1, self.n2
        shapeFunction1 = (t/a1)**n1 / (1.0 + (t / (a1)) ** n1)
        shapeFunction2 = (1.0 + (t/a2)**n2)
        E = (self.Emax - self.Emin) * alpha * shapeFunction1/shapeFunction2 + self.Emin
        alpha = self.alpha
        Emax = self.Emax
        Emin = self.Emin
        #dEdt = (Emax-Emin)*alpha*((t/a1)**(2*n1)*(t/a2)**n2*n2-(t/a2)**n2*(t/a1)**n1*n1+(t/a1)**n1*(t/a2)**n2*n2-(t/a1)**n1*n1)/((1+(t/a1)**n1)**2*(1+(t/a2)**n2)**2*t)
        return E
    
    def readDatFile(self, fileName, scaleFactor=1):
    
        f = open(fileName, 'r')
        t_array = []
        P_array = []
        for line in f:
            
            line_split = line.split('    ')
            t, P = float(line_split[0]), float(line_split[1])
            t_array.append(t)
            P_array.append(P*scaleFactor)
        
        f.close()
        
        return np.array(t_array) - t_array[0], np.array(P_array)

    def makeDataPeriodic(self, x, nCycles):
        x = x - x[0]
        x_new = np.array([])
        
        for n in range(nCycles):
            if n > 0:
                x_new = np.append(x_new, x + x_new[-1])
            
            else:
                x_new = np.append(x_new, x)
        return x_new
    
    def solveNCycle(self, Ncycles=10):
        
        t = self.t.copy()
        
        N = len(t)
        Q = np.zeros(N)
        P = np.zeros(N)
        P_LV = np.zeros(N)
        V = np.zeros(N)
        E = np.zeros(N)
        P[0] = self.Pao_0
        P_LV[0] = self.Elastance(t[0])*(self.Vinit - self.V0)
        E[0] = self.Elastance(t[0])
        V[0] = self.Vinit
        
        
        for n in range(len(t) - 1):
            dt = t[n+1] - t[n]
            if Q[n]>=0 and P_LV[n] >= P[n]:
                Vnext, P_ao_next, P_LV_next, Qnext = self.solveNextSystoleV2(t[n], dt, V[n], Q[n], P[n])
            elif self.Pv < P_LV[n] < P[n]:
                Vnext, P_ao_next, P_LV_next, Qnext = self.solveNextIsoV2(t[n], dt, V[n], P[n], P_LV[n])
            else:
                Vnext, P_ao_next, P_LV_next, Qnext = self.solveNextDiastoleV2(t[n], dt, V[n], P[n], P_LV[n])
            
            if Qnext < 0:
                P_LV_next -= 1e-14 
            V[n + 1] = Vnext
            Q[n + 1] = Qnext
            P_LV[n + 1] = P_LV_next
            P[n + 1] = P_ao_next
            E[n + 1] = self.Elastance(t[n + 1])
        
        return t, P, P_LV, Q, E, V
        
    

    def solveNextSystoleV2(self, t, dt, V, Q, P):
        
        [Vnext, Qnext, P_ao_next] = rk4(self.fsystole, [V, Q, P], [t, t + dt])[1]
        P_LV_next = P_ao_next
        return Vnext, P_ao_next, P_LV_next, Qnext
    
    def solveNextIsoV2(self, t, dt, V, P_ao, P_LV, dtt=1e-6):


        [Vnext, P_LV_next] = rk4(self.fisocomp, [V, P_LV], [t, t + dt])[1]
        P_ao_next = P_ao*np.exp(-dt/(self.R2*self.C))
        
        return V, P_ao_next, P_LV_next, 0
    
    def solveNextDiastoleV2(self, t, dt, V, P_ao, P_LV, dtt=1e-6):
        
        [Vnext, P_LV_next] = rk4(self.fdiastole, [V, P_LV], [t, t + dt])[1]
        P_a0_next = P_ao*np.exp(-dt/(self.R2*self.C))
        
        return Vnext, P_a0_next, P_LV_next, 0

    def fsystole(self, u, t, dtt=1e-3):
        
        """The coupled differential equation of windkessel and Varying Elastance during systole"""
        
        E = self.Elastance(t)
        dEdt = (self.Elastance(t + dtt) - self.Elastance(t))/dtt
        A = self.C*((dEdt)*(u[0]-self.V0)- E*u[1])
        B = (u[2]/self.R2 -u[1]*(1+ self.R1/self.R2))
        D = dEdt*(u[0]-self.V0)-E*u[1]
        
        
        return[-u[1],(A + B)/(self.R1*self.C),D]
        
    def fdiastole(self, u, t, dtt=1e-3):
        
        """The coupled differential equations during diastole, Varying Elastance and Volume flow from left Atrium"""
            
        
        E = self.Elastance(t)
        dEdt = (self.Elastance(t + dtt) - self.Elastance(t))/dtt
        A = -(u[1]-self.Pv)/self.Rv
        B = dEdt*(u[0] - self.V0)+E*A
        
        return [A,B]
    
    def fisocomp(self, u, t, dtt=1e-3):
        
        """The Varying Elastance during isoVolumetric Contraction and Relaxation
           in differential form"""
        
        E = self.Elastance(t)
        dEdt = (self.Elastance(t + dtt) - self.Elastance(t))/dtt
        
        
        return [0, dEdt*(u[0] - self.V0)]
    
        

# if __name__ == "__main__":
#     
#     T = 1.
#     N = 1000*5
#     cardiacCycles = 5
#     t = np.linspace(0, T*cardiacCycles, N*cardiacCycles + 1)
#     
#     veWK3 = varyingElastance(t, T)
#     
#     t, P, P_LV, Q, E, V = veWK3.solveNCycle(Ncycles=cardiacCycles)
#     
#     plt.figure()
#     plt.plot(t, P)
#     plt.plot(t, P_LV)
#     
#     plt.figure()
#     plt.plot(t, Q)
# 
#     plt.figure()
#     plt.plot(t, V)
# 
#     plt.figure()
#     plt.plot(P_LV, V)
# 
#     plt.figure()
#     plt.plot(t[:-1], (E[1:]-E[:-1])/(t[1:]-t[:-1]))
#     plt.figure()
#     plt.plot(t, E)
# 
#     plt.show()

#     fig, ax1 = plt.subplots()
#     
#     color = 'tab:red'
#     ax1.set_xlabel('time (s)')
#     ax1.set_ylabel('exp', color=color)
#     ax1.plot(t, E, color=color)
#     ax1.tick_params(axis='E', labelcolor=color)
#     
#     ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
#     
#     color = 'tab:blue'
#     ax2.set_ylabel('dEdt', color=color)  # we already handled the x-label with ax1
#     ax2.plot(t[:-1], (E[1:]-E[:-1])/(t[1:]-t[:-1]), color=color)
#     ax2.tick_params(axis='y', labelcolor=color)
#     
#     fig.tight_layout()  # otherwise the right y-label is slightly clipped
#     plt.show()
    
    