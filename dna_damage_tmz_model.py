# -*- coding: utf-8 -*-
"""
Created on Thu Dec 17 21:07:24 2020

@author: Anshul Saini
"""
import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
from oral_tmz_file import Oral_TMZ
from dna_damage_ode_file import derivative
from dna_damage_repair_ode_file import derivative_mgmt
from numpy import genfromtxt




def dna_damage_tmz_model(t0, tf, kD_brca_input, kD_msh6_input, kD_mgmt_input):


    TMZ_Dose = 10000;   #%in microM  %334.6;
    TimeOfDoses = [24,48,72,96,120];    #%Time in hours of doses of TMZ
    th = 168;        #%time of simulation(s) in hours
    div = 1800;
    hr = 3600
    
    TMZ_input = Oral_TMZ(th, div, TMZ_Dose, TimeOfDoses);
    simulation_time_hr = np.linspace(0, th, int(th*hr/div))
    
    
    p53in_total= [];
    p53_total  = [];
    Mdm2_total = [];
    wip1_total = [];
    ATMP_total = [];
    ATRac_total= [];
    BRCA2_total= [];
    MSH6_total = [];
    MGMT_total = [];
    DSB_total = [];
    SSB_total = [];
    t_total = [];
    
    #Species initial values without intialization
    
    s0 = np.zeros(38)
    s0[0] = 2.966386e+02; s0[1] = 6.245528e+00; s0[2] = 2.056385e+02; s0[3] = 2.230736e+00; s0[4] = 0.000000e+00; s0[5] = 0.000000e+00; s0[6] = 6.245528e+00; 
    s0[7] = 6.245528e+00; s0[8] = 6.245528e+00; s0[9] = 6.245528e+00; s0[10] = 6.245528e+00; s0[11] = 6.245528e+00; s0[12] = 6.245528e+00; s0[13] = 6.245528e+00; 
    s0[14] = 6.245528e+00; s0[15] = 6.245528e+00; s0[16] = 6.245528e+00; s0[17] = 6.245528e+00; s0[18] = 6.245528e+00; s0[19] = 6.245528e+00; s0[20] = 6.245528e+00; 
    s0[21] = 6.245528e+00; s0[22] = 6.245528e+00; s0[23] = 6.245528e+00; s0[24] = 6.245528e+00; s0[25] = 6.245528e+00; s0[26] = 8.218084e-01; s0[27] = 4.550432e+01; 
    s0[28] = 1.202569e+02; s0[29] = 2.431635e-04; s0[30] = 6.027760e-06; s0[31] = 0.000000e+00; s0[32] = 0.000000e+00; s0[33] = 0.000000e+00; s0[34] = 2.721609e-02; 
    s0[35] = 0.000000e+00; s0[36] = 1.503781e+00; s0[37] = 2.400357e+00;
    
    
    t_max = div;
    y0 = s0;
    
    for i in range(t0, tf):
                
        t = np.linspace(0, t_max, 100);
        
                
        if i< 335:
            tmz_stim = TMZ_input[i];
            print(i);
        else:
            tmz_stim = 0;
        
        y = spi.odeint(derivative, y0, t, args=(tmz_stim, kD_brca_input, kD_msh6_input, kD_mgmt_input, ));
        
        
        for j in range(0, len(t)):
            p53in_total.append(y[j, 0])
            p53_total.append(y[j, 1])
            Mdm2_total.append(y[j, 2])
            wip1_total.append(y[j, 3])
            ATMP_total.append(y[j, 4])
            ATRac_total.append(y[j, 5])
            BRCA2_total.append(y[j, 26])
            MSH6_total.append(y[j, 27])
            MGMT_total.append(y[j, 28])
            DSB_total.append(y[j, 29])
            SSB_total.append(y[j, 30])
            t_total.append((i)*t_max + t[j])
            
        y0 = y[-1, :];    
    
        
    p53_combine = [];
    for j in range(0, len(p53_total)):
            p53_combine.append(p53_total[j] + p53in_total[j])
            
    #stoichometric = genfromtxt('SM_DNA_damage_mgmt.csv', delimiter=',')
       
    #plt.plot(t_total, p53in_total, label="p53 inactive")
    #plt.plot(t_total, p53_total, label="p53 active")
    plt.plot(t_total, p53_combine, label="p53 total")
    #plt.plot(t_total, Mdm2_total, label="Mdm2")
    #plt.plot(t_total, wip1_total, label="Wip1")
    #plt.plot(t_total, ATMP_total, label="ATM-P")
    #plt.plot(t_total, ATRac_total, label="ATRac")
    #plt.plot(t_total, BRCA2_total, label="BRCA2")
    #plt.plot(t_total, MSH6_total, label="MSH6")
    #plt.plot(t_total, MGMT_total, label="MGMT")
    #plt.plot(t_total, DSB_total, label="DSB")
    #plt.plot(t_total, SSB_total, label="SSB")
    #plt.legend(loc="upper right")
    #
    #plt.xlim(0, 55);
    plt.show()    
    #plt.plot(t_total, MGMT_total, label="MGMT_total")  
    #plt.show()    
    return p53_combine