import matplotlib.pyplot as plt
import pandas as pd
import numpy as np    

    

class CR():

    m_p = 0.9383                      # proton mass in GeV/c^2
    m_u = 0.9315                      # atomic mass in GeV/c^2
    m_e = 0.000510998918            # electron mass in GeV/c^2
    
    def __init__(self):

        return print("you can use mod, rig_to_En and En_to_Rig methods")

    
    def func(x, a, b):
        return a * pow(x, -b)#x**(-b)


    def fitfunc(x, a1, a2, A, xc):
        if len(x)>1:
            for x_ in x:
                if x_ < xc:
                    y.append(A * x_**(-a1))
                elif x_ > xc:
                    y.append(A * xc**(a1-a2) * x_**(-a2))
                else:
                    y.append(0.)
        else:
            if x < xc:
                y.append(A * x**(-a1))
            elif x > xc:
                y.append(A * xc**(a1-a2)*x**(-a2))
            else:
                    y.append(0.)
        return np.array(y)


    def fit2func(x, a1, a2, a3, A, xc, xc2):
        y = []
        Ec = 1e6
        if len(x)>1:
            for x_ in x:
                if x_ < xc:
                    y.append(A * x_**(-a1) * np.exp(-x_/Ec))
                elif x_ > xc and x_ <xc2:
                    y.append(A * xc**(a1-a2 +1.95)*x_**(-a2)* np.exp(-x_/Ec))
                elif x_ > xc2:
                    y.append( A * xc**(a1-a2+1.2)* xc2**(a2-a3+1.9)* x_**(-a3)* np.exp(-x_/Ec))
                else:
                    y.append(0.)
        else:
            if x < xc:
                y.append(A * x**(-a1))
            elif x > xc and x <xc2:
                y.append(A * xc**(a1-a2 +1.95)*x**(-a2))
            elif x > xc2:
                y.append( A * xc**(a1-a2+1.2)* xc2**(a2-a3+1.9)* x**(-a3))
            else:
                y.append(0.)
        return np.array(y)


    def fit3func(x, a1, a2, a3, a4, A, xc, xc2, xc3):
        y = []
        Ec = 1e6
        if len(x)>1:
            for x_ in x:
                if x_ < xc:
                    y.append(A * x_**(-a1) * np.exp(-x_/Ec))
                elif x_ > xc and x_ <xc2:
                    y.append(A * xc**(a1-a2)*x_**(-a2)* np.exp(-x_/Ec))
                elif x_ > xc2 and x_ <xc3:
                    y.append(A * xc**(a1-a2)*xc2**(a2-a3)*x_**(-a3)* np.exp(-x_/Ec))
                elif x_ > xc3:
                    y.append(A * xc**(a1-a2) * xc2**(a2-a3) * xc3**(a3-a4) * x_**(-a4)
                             * np.exp(-x_/Ec))
                else:
                    y.append(0.)
        else:
            if x < xc:
                y.append(A * x_**(-a1) * np.exp(-x/Ec))
            elif x > xc and x_ <xc2:
                y.append(A * xc**(a1-a2)*x_**(-a2)* np.exp(-x/Ec))
            elif x > xc2 and x_ <xc3:
                y.append(A * xc**(a1-a2)*xc2**(a2-a3)*x**(-a3)* np.exp(-x/Ec))
            elif x > xc3:
                y.append(A * xc**(a1-a2) * xc2**(a2-a3) * xc3**(a3-a4) * x**(-a4)
                         * np.exp(-x_/Ec))
        return np.array(y)


    
    def modul(E, flux, Z, A, phi_nuc, Charge_sign=False, phi_n= 0.7, fluxerr = 0):  ##IMPLEMENTED FOR THE FORCE-FIELD APPROX WITH THE MODIFICATION OF 
                                                                                    ##CHOLIS-HOOPER-LINDEN (arXiv:1511.01507) FOR CHARGE-SIGN DEPENDENCE
        if Z==1 and A==1:
            m = CR.m_p
        else:
            m = CR.m_u
        
        E = np.array(E)
          
        if Charge_sign == True:         
            R0 = 1.
            R = (A/Z) * np.sqrt(E**2 + 2*m*E)
            phi_ = phi_n
            phi_n = phi_nuc + phi_*(R0/R)
            En = E + (Z/A)*phi_n
            mod = E * (E + 2*m) / ((E + (Z/A)*phi_n + m)**2 - m**2)
            flux_mod = np.interp(En, E, flux)*mod
            
        else:
            En = E + (Z/A)*phi_nuc
            mod = E * (E + 2*m) / ((E + (Z/A)*phi_nuc + m)**2 - m**2)
            flux_mod = np.interp(En, E, flux)*mod
        
        if fluxerr != 0: 
            fluxerr_mod = np.interp(En, E, fluxerr)*mod
            return np.array(flux_mod), np.array(fluxerr_mod)
          
        return np.array(flux_mod)#, np.array(fluxerr_mod)
    

    def Rig_to_En(R, A, Z, flux, errflux):
        if Z==1 and A==1:
            m = CR.m_p
        else:
            m = CR.m_u

        Rmean = np.array(R)
        
        #M = A*m
        #T = (np.sqrt((Z*Rmean)**2 + M**2) - M )/A       # T is the Emean!!
        
        T = (np.sqrt((Z*Rmean/A)**2 + m**2) - m )
        dRdE = (A/Z)*(np.array(T)+m)/np.sqrt((np.array(T)**2) + 2. * np.array(T)*m)
        Flux = dRdE * np.array(flux)
        ErrFlux = dRdE * np.array(errflux)
        

        return T, Flux, ErrFlux
        

    def En_to_Rig(E, A, Z, flux, errflux):

        if Z==1 and A==1:
            m = CR.m_p
        else:
            m = CR.m_u
            
        Emean = np.array(E)
        R = (A/Z) * np.sqrt(Emean**2 + 2*m*Emean)
        dEdR = (Z/A)**2 * R/((R*Z/A)**2 + m**2)
        Flux = dEdR * np.array(flux)
        ErrFlux = dEdR * np.array(errflux)
        
        return R, Flux, ErrFlux