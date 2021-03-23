# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 19:21:59 2021

@author: anshulsa
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from oral_ode_file import ORAL_ODE



def Oral_TMZ(th, div, TMZ_Dose,TimeOfDoses):

    Vd = 30.; #% L, apparent volume of distribution

    #%TMZ Cell death parameters
    kDeadTMZ = 2.0e-4; #% hr-1, death rate constant due to TMZ
      
    #%EQUATION:
    aTgi_0 = 0;  #%aTgi = TMZ dose into gut and is equal to 1300
    
    #% ICs for ODE
    aTbl_0 = 0;
    Tif_0 = 0;
    Tic_0 = 0;
    Mic_0 = 0;
    meCic_0 = 0;
    DNAadd_0 = 0;
        
    u0 = [aTgi_0, aTbl_0, Tif_0, Tic_0, Mic_0, meCic_0, DNAadd_0];
    
    hr = 3600;
    
    c_meth_list = [];  
    DNA_adduct_list =  [];
    dose_reg = np.size(TMZ_Dose);
    
    if TimeOfDoses[0] == 0:
       
        C0 = TMZ_Dose; #%in microM
        u0[0] = C0;
    
        #%Time in hours
        delta=(TimeOfDoses[1]*hr)/div;
        t= np.linspace(0, TimeOfDoses[1],delta);
                
        u_average = odeint(ORAL_ODE, u0, t);
        
        for i in range(0, len(t)):
            c_meth_list.append(u_average[i, 5])
            DNA_adduct_list.append(u_average[i, 6])
        
        for i in range(0, len(TimeOfDoses)-1):
            
            time_dif = TimeOfDoses[i+1]-TimeOfDoses[i];
            delta=(time_dif*hr)/div;
            t = np.linspace(0,time_dif,delta);
            
            u0 = u_average[-1, :];
            
            if dose_reg > 1:
                C0 = TMZ_Dose;
                u0[0] = u0+C0;
            else:
                C0 = TMZ_Dose;
                u0[0] = u0[0]+C0;
        
            u_average = odeint(ORAL_ODE, u0, t);
            
            for i in range(0, len(t)):
                c_meth_list.append(u_average[i, 5])
                DNA_adduct_list.append(u_average[i, 6])
           
        
        time_dif = th-TimeOfDoses[-1];
        delta = int((time_dif*hr)/div);
        t = np.linspace(0,time_dif,delta);
        
        
        u0 = u_average[-1, :];
        C0 = TMZ_Dose;
        u0[0] = u0[0]+C0;
                
        u_average = odeint(ORAL_ODE, u0, t);
        for i in range(0, len(t)):
            c_meth_list.append(u_average[i, 5])
            DNA_adduct_list.append(u_average[i, 6])
  
    
    else:
       
        C0 = 0; # %in microM
        u0[0] = C0;
        print(u0)
        delta = int((TimeOfDoses[0]*hr)/div);
        t = np.linspace(0,TimeOfDoses[0],delta);
   #     print(len(t))
        
        u_average = odeint(ORAL_ODE, u0, t);
        
        for i in range(0, len(t)):
            c_meth_list.append(u_average[i, 5])
            DNA_adduct_list.append(u_average[i, 6])
   #     print(len(c_meth_list))
            #print(i)
        
        for i in range(0, len(TimeOfDoses)-1):
            
            time_dif = TimeOfDoses[i+1]-TimeOfDoses[i];
            delta = int((time_dif*hr)/div);
            t = np.linspace(0,time_dif,delta);
            #print(i)
    #       print(len(t))
            u0 = u_average[-1, :];
            #print(u0)
            if dose_reg > 1:
                C0 = TMZ_Dose;
                u0[0] = C0 + u0[0];
            else:
                C0 = TMZ_Dose;
                u0[0] = C0 + u0[0];
            #print(u0)       
            u_average = odeint(ORAL_ODE, u0, t);
            #print(u_average[-1]);
            for i in range(0, len(t)):
                c_meth_list.append(u_average[i, 5])
                DNA_adduct_list.append(u_average[i, 6])
  #          print(len(c_meth_list))
        
        
        time_dif = th - TimeOfDoses[-1];
        delta = int((time_dif*hr)/div);
        t = np.linspace(0, time_dif, delta);
 #      print(len(t))
        
        u0 = u_average[-1, :];
        C0 = TMZ_Dose;
        u0[0] = C0 + u0[0];
        
        u_average = odeint(ORAL_ODE, u0, t);   

        for i in range(0, len(t)):
            c_meth_list.append(u_average[i, 5])
            DNA_adduct_list.append(u_average[i, 6])
        
   #     print(len(c_meth_list))
    
    
    c_meth_list_cali = []
    DNA_adduct_list_cali = []
    
    c_meth_list_cali = [i*(7*(10**-12)/(2*(10**-12)))*1000 for i in c_meth_list] #%Convert to nM
    DNA_adduct_list_cali = [i*(7*(10**-12)/(2*(10**-12)))*1000 for i in DNA_adduct_list] #%Convert to nM

#    amt_bl = plasma/Vd;
#     
    tplot = np.linspace(0, th, int(th*hr/div));
#    print(len(tplot))
    plt.plot(tplot, c_meth_list_cali);
    plt.ylabel('Methylating Cation (nM)')
    plt.xlabel('Time Hours')
    plt.show()
    
    plt.plot(tplot, DNA_adduct_list_cali);
    plt.ylabel('DNA_adduct (nM)')
    plt.xlabel('Time Hours')
    plt.show()
    
    return DNA_adduct_list_cali
