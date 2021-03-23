# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 19:34:49 2021

@author: anshulsa
"""

import numpy as np
from numpy import genfromtxt

stoichometric_mat = genfromtxt('SM_DNA_damage.csv', delimiter=',')

def derivative(s, t, tmz_stim):
    
    #dS = np.zeros(38)
    vD = np.zeros((66, 1)) 
    kD = np.zeros(60)
    dS = np.zeros(38)
    
   # kD **
    kD[1]=0.9*(1000/3600);
    kD[2]=5/(1000*3600);
    kD[3]=2/3600;
    kD[4]=10/3600;
    kD[5]=4; 
    kD[6]=1*1000;
    kD[7]=1.4/(1000*3600) ;
    kD[8]=0.14/(1000*3600);
    kD[9]=0.9/3600;
    kD[10]=0.2*(1000/3600);
    kD[11]=1/3600;
    kD[12]=0.5/(1000*3600);
    kD[13]=0.25/3600;
    kD[14]=0.7/3600;
    kD[15]=0.1/(1000*3600); 
    kD[16]=50/3600; 
    kD[17]=4; 
    kD[18]=0.2*1000;
    kD[19]=7.5/3600;
    kD[20]=10*(1000/3600);
    kD[21]=4.5*(1000/3600);#0.60754*(1000/3600);
    kD[22]=0.033948*3600;
    kD[23]=0.058614*3600;
    kD[24]=0.0035778*3600;
    kD[25]=0.33705*3600;
    kD[26]=0.013496*3600;
    kD[27]=0.030939*3600;
    kD[28]=0.0042462*3600;
    kD[29]=0.081109*3600;
    kD[30]=0.10093*3600;
    kD[31]=0.13295*3600;
    kD[32]=0.080377*3600;
    kD[33]=0.12247*3600;
    kD[34]=0.16522*3600;
    kD[35]=0.099235*3600;
    kD[36]=0.13629*3600;
    kD[37]=0.28505*3600;
    kD[38]=0.10575*3600;
    kD[39]=0.066299*3600;
    kD[40]=0.11092*3600;
    kD[41]=0.042434*3600;
    kD[42]=0;
    kD[43]=0;
    kD[44]=0;
    kD[45]=0;

    kD[46]=1E-5*5;
    kD[47]=1E-5;
    kD[48]=1E-5;

    kD[49]=0.0005;
    kD[50]=0;
    kD[51]=(1E-4)*500000;
    kD[52]=5000;
    kD[53]=tmz_stim;
    kD[54]=5; 
    kD[55]=20; 
    kD[56]=20; 
    kD[57]=20; 
    kD[58]= 2;
    kD[59]=2;
    Ma = 5.230808e+00;
    Me = 2.746016e+00;
    
   
    damage = 1;
    c= 0;
    
    p53inac =s[0];
    p53ac =s[1];
    Mdm2 =s[2];
    Wip1 =s[3];
    ATMP =s[4];
    ATRac =s[5];
    Mdm2product1=s[6];
    Mdm2product2=s[7];
    Mdm2product3=s[8];
    Mdm2product4=s[9];
    Mdm2product5=s[10];
    Mdm2product6=s[11];
    Mdm2product7=s[12];
    Mdm2product8=s[13];
    Mdm2product9=s[14];
    Mdm2pro =s[15];
    Wip1product1=s[16];
    Wip1product2=s[17];
    Wip1product3=s[18];
    Wip1product4=s[19];
    Wip1product5=s[20];
    Wip1product6=s[21];
    Wip1product7=s[22];
    Wip1product8=s[23];
    Wip1product9=s[24];
    Wip1pro =s[25];
    BRCA2=s[26];
    MSH6=s[27];
    MGMT=s[28];
    damageDSB=s[29];
    damageSSB=s[30];
    ppAKT_Mdm2=s[31];
    pMdm2=s[32];
    ARF=s[33];
    MDM4=s[34];
    p53ac_MDM4=s[35];
    ATMinac=s[36];
    ATRinac=s[37];
    

  
    bp=kD[1];
    ampi=kD[2];
    api=kD[3];
    bsp=kD[4];
    ns=kD[5];
    Ts=kD[6];
    ampa=kD[7];
    awpa=kD[8];
    bm=kD[9];
    bmi=kD[10];
    am=kD[11];
    asm=kD[12];
    bw=kD[13];
    aw=kD[14];
    asm2=kD[15];
    aws=kD[16];
    nw=kD[17];
    Tw=kD[18];
    as1=kD[19]; # as is replaced by as1 due to python command
    bs=kD[20];
    bs2=kD[21];
    tau1=kD[22];
    tau2=kD[23];
    tau3=kD[24];
    tau4=kD[25];
    tau5=kD[26];
    tau6=kD[27];
    tau7=kD[28];
    tau8=kD[29];
    tau9=kD[30];
    tau10=kD[31];
    tau11=kD[32];
    tau12=kD[33];
    tau13=kD[34];
    tau14=kD[35];
    tau15=kD[36];
    tau16=kD[37];
    tau17=kD[38];
    tau18=kD[39];
    tau19=kD[40];
    tau20=kD[41];
    kD[42]=kD[42];
    kD[43]=kD[43];
    kD[44]=kD[44];
    kD[45]=kD[45];
    fixdsb1 =kD[46];
    fixmsh =kD[47];
    fixmgmt =kD[48];
    basalp53act=kD[49];
    kDDbasal=kD[50];
    kDDE=kD[51];
    kDEtop=kD[52];
    Etop=kD[53];
    kDnSP=kD[54];
    kDkmSP=kD[55];
    kDnSS=kD[56]; 
    kDnDS=kD[57]; 
    kDkmSS=kD[58]; 
    kDkmDS=kD[59];
    
    
    

    vD[0] = bp;
    vD[1] = ampi*Mdm2*p53inac; 
    vD[2] = basalp53act + bsp*p53inac*(ATMP**ns/(ATMP**ns+Ts**ns)+ATRac**ns/(ATRac**ns+Ts**ns)); 
    vD[3] = awpa*Wip1*p53ac;
    vD[4] = api*p53inac; 
    vD[5] = ampa*Mdm2*p53ac; 
    vD[6] = bm*Mdm2pro;
    vD[7] = bmi;
    vD[8] = asm*Mdm2*ATMP; 
    vD[9] = asm2*Mdm2*ATRac; 
    vD[10] = am*Mdm2;
    vD[11] = bw*Wip1pro;
    vD[12] = aw*Wip1; 
    vD[13] = bs*((damageDSB**kDnDS)/((kDkmDS**kDnDS)+(damageDSB**kDnDS))); #(damage - c)#
    vD[14] = aws*ATMP*Wip1**nw/(Wip1**nw+Tw**nw); 
    vD[15] = as1*ATMP; 
    vD[16] = bs2*((damageSSB**kDnSS)/((kDkmSS**kDnSS)+(damageSSB**kDnSS)));#(damage -c)
    vD[17] = as1*ATRac; 
    vD[18] = Mdm2product1/tau1;
    vD[19] = p53ac/tau1;
    vD[20] = Mdm2product2/tau2;
    vD[21] = Mdm2product1/tau2;
    vD[22] = Mdm2product3/tau3;
    vD[23] = Mdm2product2/tau3;
    vD[24] = Mdm2product4/tau4;
    vD[25] = Mdm2product3/tau4;
    vD[26] = Mdm2product5/tau5;
    vD[27] = Mdm2product4/tau5;
    vD[28] = Mdm2product6/tau6;
    vD[29] = Mdm2product5/tau6;
    vD[30] = Mdm2product7/tau7;
    vD[31] = Mdm2product6/tau7;
    vD[32] = Mdm2product8/tau8;
    vD[33] = Mdm2product7/tau8;
    vD[34] = Mdm2product9/tau9;
    vD[35] = Mdm2product8/tau9;
    vD[36] = Mdm2pro/tau10;
    vD[37] = Mdm2product9/tau10;
    vD[38] = Wip1product1/tau11;
    vD[39] = p53ac/tau11;
    vD[40] = Wip1product2/tau12;
    vD[41] = Wip1product1/tau12;
    vD[42] = Wip1product3/tau13;
    vD[43] = Wip1product2/tau13;
    vD[44] = Wip1product4/tau14;
    vD[45] = Wip1product3/tau14;
    vD[46] = Wip1product5/tau15;
    vD[47] = Wip1product4/tau15;
    vD[48] = Wip1product6/tau16;
    vD[49] = Wip1product5/tau16;
    vD[50] = Wip1product7/tau17;
    vD[51] = Wip1product6/tau17;
    vD[52] = Wip1product8/tau18;
    vD[53] = Wip1product7/tau18;
    vD[54] = Wip1product9/tau19;
    vD[55] = Wip1product8/tau19;
    vD[56] = Wip1pro/tau20;
    vD[57] = Wip1product9/tau20;
    vD[58] = kD[42]*ARF*Mdm2;
    vD[59] = kD[43]*ARF*pMdm2;
    vD[60] = kD[44]*MDM4*p53ac;
    vD[61] = kD[45]*p53ac_MDM4;
    vD[62] = fixdsb1*BRCA2*damageDSB;
    vD[63] = fixmsh*MSH6*damageSSB;
    vD[64] = fixmgmt*MGMT*damageSSB;
    vD[65] = (kDDbasal+kDDE*(Etop/(Etop+kDEtop)))*(((Ma+Me)**kDnSP)/(((Ma+Me)**kDnSP)+(kDkmSP**kDnSP))); #(damage -c)

    
    dS1 = np.dot(stoichometric_mat, vD) 
#    print(dS1)
 #   print(stoichometric_mat)
    for i in range(0, 38):
         dS[i] = dS1[i]
    return dS

