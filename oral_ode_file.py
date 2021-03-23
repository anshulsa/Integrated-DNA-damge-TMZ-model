# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 18:07:42 2021

@author: anshulsa
"""

import numpy as np


def ORAL_ODE(u, t):
 
    fif = 0.3;# % unitless, volume fraction of  interstitial fluid (IF) in brain tumor
    fbl = 0.05;# % unitless, volume fraction of blood in brain tumor
    fic = 0.65; #% unitless, volume fraction of intracellular (IC) compartment in brain tumor
    Vbl = 0.0025; #% L, blood volume of brain tumor, Vtu*fbl
    Vtu = 0.05; #% L, total volume of the brain tumor
    Vif = 0.015; #%L, volume of the tumor IF compartment
    Vic = 0.0325;# % L, volume of the IC compartment
    Nic0 = 4.063e9;# % number of toal cells in IC compartment at time=0
    
    fGP = 0.10;# % unitless, initial fraction of IC cells occupied by GP cells
    fQ = 0.80; #% unitless, initial fraction of IC cells occupied by Q cells
    fGS0 = 0.05; #% unitless, initial fraction of GS0 cells
    fGS1 = 0.025; #% unitless, initial fraction of GS1 cells
    fGS2 = 0.025; #% unitless, initial fraction of GS2 cells 
    
    conNV = 1.25e11; #% conversion of cell number to L
    
    #% TMZ pH-dependent metabolic degradation
    pHbl = 7.4; #% blood pH
    pHif = 7.3; #% IF pH
    pHic = 7.2; #% IC pH
    kT0 = 1.1e-7;# % hr-1, baseline TMZ degradation constant
    kM0 = 292; #% hr-1, baseline MTIC degradation constant
    lambdaT = 2.08; #% TMZ empirical
    lambdaM = 0.31; #% MTIC empirical
    
    
    #% TMZ Transport Parameters
    qT = 0.019; #% L/hr, TMZ blood to IF clearance
    qT2 = 0.069; #% L/hr, TMZ IF to blood clearance
    pT = 17.5; #% L/hr, TMZ IF to IC clearance
    pT2 = 36.3; #% L/hr, TMZ IC to IF clearance
    
    #% TMZ forcing function/systemic
    ka = 3.0; #% hr-1, TMZ absorption rate constant
    kel = 0.35; #% hr-1, TMZ elimination rate constant
    Vd = 30.; #% L, apparent volume of distribution
    
    #%TMZ PD
    kadd = 1.81; #% hr-1, DNA adduct formation rate constant
    kdeg = 6000.; #% hr-1, methylating cation degradation constant
    kDNAoff = 0.0041; #% hr-1, DNA adduct degradation rate constant
    
    #%Cell-Type Parameters
    DNAaddMAX = 7.4e-3; #% uM, Maximum DNA adduct concentration without kDNAoff
    sigDNAaddGP = 1.0; #% unitless, empirial amplification exponent for DNA adduct effect
    sigDNAaddGS0 = 1.0; #% unitless, empirial amplification exponent for DNA adduct effect
    proGP =1.8e-3; #% 1/hr, proliferation rate constant of GP cells
    sigGP = 0.18;  #% unitless, exponent in logistic equation for GP 
    
    proGS0 =1.8e-3; #% 1/hr, proliferation rate constant of glioma stem cell state 0 (GS0) cells
    sigGS0 = 0.18;  #% unitless, exponent in logistic equation for GS0
    
    kQGP = 2.0e-5; #% hr-1, transfer rate constant of Q cells to proliferating glioma (GP) cells 
    kGS2GP = 2.0e-5;  #% hr-1, transfer rate constant of glioma stem cell state 2 (GS2) to GP
    kGPQ = 4.0e-4; #% hr-1, transfer rate constant of GP cells to Q cells
    kGPGS0 = 2.0e-5; #% hr-1, transfer rate constant of GP cells to GS0 cells
    kGS0Q = 4.0e-4;  #% hr-1, transfer rate constnat of glioma stem cell state 0 (GS0) to Q cells
    kQGS0 = 2.0e-5;  #% hr-1, transfer rate constnat of  Q cells to GS0 cells
    GStau = 2.0e-5;  #% hr-1, transfer rate constant of GS between states
    
#   rename input u_0
    aTgi = u[0];
    aTbl = u[1];
    Tif = u[2];
    Tic = u[3];
    Mic = u[4];
    meCic = u[5];
    DNAadd = u[6];

#   pH-dependent rate constants
    kTif = kT0*np.exp(lambdaT*pHif);
    kTic = kT0*np.exp(lambdaT*pHic);
    kMic = kM0*np.exp(lambdaM*pHic);
    
#   ODEs
    ddt_aTgi = -ka*aTgi;
    ddt_aTbl = ka*aTgi - kel*aTbl;
    Tbl = aTbl/Vd;
    
    ddt_Tif = (qT*Tbl - qT2*Tif - pT*Tif + pT2*Tic)/Vif - kTif*Tif;
    ddt_Tic = ( pT*Tif - pT2*Tic)/Vic - kTic*Tic;
    
    ddt_Mic = kTic*Tic - kMic*Mic;
    ddt_meCic = kMic*Mic - kadd*meCic - kdeg*meCic;
    
    ddt_DNAadd = kadd*meCic - kDNAoff*DNAadd;
    
    
    
#    ddt_u1 = -ka*u[1];
#    ddt_aTbl = ka*u[1] - kel*u[2];
#    Tbl = u[2]/Vd;
#    
#    ddt_Tif = (qT*Tbl - qT2*u[3] - pT*u[3] + pT2*u[4])/Vif - kTif*u[3];
#    ddt_Tic = ( pT*u[3] - pT2*u[4])/Vic - kTic*u[4];   
#    ddt_Mic = kTic*u[4] - kMic*u[5];
#    ddt_meCic = kMic*u[5] - kadd*u[6] - kdeg*u[6];
#    ddt_DNAadd = kadd*u[6] - kDNAoff*u[7];
       
    w = np.zeros(7);
    w[0] = ddt_aTgi;
    w[1] = ddt_aTbl;
    w[2] = ddt_Tif;
    w[3] = ddt_Tic;
    w[4] = ddt_Mic;
    w[5] = ddt_meCic;
    w[6] = ddt_DNAadd;
    
    dudt = np.transpose(w);
    return dudt   
    