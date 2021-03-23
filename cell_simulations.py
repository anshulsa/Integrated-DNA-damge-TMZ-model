"""
Created on Thu Dec 17 21:07:24 2020

@author: Anshul Saini
"""
import numpy as np
import random as rand
import scipy.integrate as spi
import matplotlib.pyplot as plt
from oral_tmz_file import Oral_TMZ
from dna_damage_ode_file import derivative
from scipy.interpolate import InterpolatedUnivariateSpline




TMZ_Dose = 1000;   # in microM  334.6;
TimeOfDoses = [24,48,72,96,120];    # Time in hours of doses of TMZ
th = 200;        # time of simulation(s) in hours
div = 1800;      # Number of seconds after which TMZ output is updated
hr = 3600;

TMZ_input = Oral_TMZ(th, div, TMZ_Dose, TimeOfDoses);  # TMZ output function

simulation_time_hr = np.linspace(0, th, int(th*hr/div)); # Array of simulation time in hours.

# The cell vector containing all the details regarding the cells in the tumor. 

cell_vector =  [[1, 1, 1, 1, 1],          # Cell array each entry represnt a single cell
                [0, 0, 0, 0, 0],          # Time of creation of corresponding cells
                [1, 1, 1, 1, 1],          # BRCA2 concentration for each cell
                [1, 1, 1, 1, 1],          # MSH6 concentration for each cell
                [1, 1, 1, 1, 1],          # MGMT concentration for each cell
                [1, 2, 3, 4, 5]]          # 1 - Q cells  2- GS_Tr 3- GS_Tp 4- GP_Tr 5- GP_Tp

cell_count = [[1], [1], [1], [1], [1]]
time_count = [[0], [0], [0], [0], [0]]

brca2 = 1; msh6 = 1; mgmt = 1;

DNA_damage_integral_list = [];

#Species initial values without intialization

s0 = np.zeros(38);

s0[0] = 2.966386e+02; s0[1] = 6.245528e+00; s0[2] = 2.056385e+02; s0[3] = 2.230736e+00; s0[4] = 0.000000e+00; s0[5] = 0.000000e+00; s0[6] = 6.245528e+00; 
s0[7] = 6.245528e+00; s0[8] = 6.245528e+00; s0[9] = 6.245528e+00; s0[10] = 6.245528e+00; s0[11] = 6.245528e+00; s0[12] = 6.245528e+00; s0[13] = 6.245528e+00; 
s0[14] = 6.245528e+00; s0[15] = 6.245528e+00; s0[16] = 6.245528e+00; s0[17] = 6.245528e+00; s0[18] = 6.245528e+00; s0[19] = 6.245528e+00; s0[20] = 6.245528e+00; 
s0[21] = 6.245528e+00; s0[22] = 6.245528e+00; s0[23] = 6.245528e+00; s0[24] = 6.245528e+00; s0[25] = 6.245528e+00; 

s0[26] = 8.218084e-01;           # Initial BRCA2 concentration
s0[27] = 4.550432e+01;           # Initial MSH6 concentration 
s0[28] = 1.202569e+02;           # Initial MGMTconcentration 

s0[29] = 2.431635e-04; s0[30] = 6.027760e-06; s0[31] = 0.000000e+00; s0[32] = 0.000000e+00; s0[33] = 0.000000e+00; s0[34] = 2.721609e-02; 
s0[35] = 0.000000e+00; s0[36] = 1.503781e+00; s0[37] = 2.400357e+00;

k = 0;

            
while k in range(len(cell_vector[0])):
    
    
    
    brca2 = cell_vector[2][k]; msh6 =  cell_vector[3][k]; mgmt =  cell_vector[4][k]; # repair protein concentration for Kth cell
    
    
    #Stochastic fluctuation in repair protein expressions. For each cell state, the reapir protein value 
    #fluctuate based on Normal function. For transcriptionally permissive cell state the fluctuations are set
    #to higher values. For transcriptionally recessive cells states fluctuations are set to be lower than Q cells.
    
    if cell_vector[5][k] == 1:   # Qcells 
        brca2 = np.random.normal(brca2, brca2/10, 1); msh6 = np.random.normal(msh6, msh6/10, 1); mgmt =  np.random.normal(mgmt, mgmt/10, 1); 
    elif cell_vector[5][k] == 2: # GS_Tr
        brca2 = np.random.normal(brca2, brca2/20, 1); msh6 = np.random.normal(msh6, msh6/20, 1); mgmt =  np.random.normal(mgmt, mgmt/20, 1);   
    elif cell_vector[5][k] == 3: # GS_Tp
        brca2 = np.random.normal(brca2, brca2/5, 1); msh6 = np.random.normal(msh6, msh6/5, 1); mgmt =  np.random.normal(mgmt, mgmt/5, 1); 
    elif cell_vector[5][k] == 4: # GP_Tr
        brca2 = np.random.normal(brca2, brca2/20, 1); msh6 = np.random.normal(msh6, msh6/20, 1); mgmt =  np.random.normal(mgmt, mgmt/20, 1); 
    elif cell_vector[5][k] == 5: # GP_Tp
        brca2 = np.random.normal(brca2, brca2/5, 1); msh6 = np.random.normal(msh6, msh6/5, 1); mgmt =  np.random.normal(mgmt, mgmt/5, 1); 
    

    t = np.linspace(0, div, 100);    
    t_cell_cycle = 0;
    initial_time = cell_vector[1][k];
    tmz_stim = TMZ_input[initial_time];
    p53in_total= [];
    p53_total  = [];
    p53_combine = [];
    DSB_total = [];
    SSB_total = [];
    t_total = [];
    
    #Updated value of repair proteins after the fluctuations
    
    s0[26] = brca2; 
    s0[27] = msh6; 
    s0[28] = mgmt; 
    
        
    y0 = s0;
    
    for i in range(initial_time, len(simulation_time_hr)):
                
        tmz_stim = TMZ_input[i];         #TMZ input into the DNA damage repair model
        
        y = spi.odeint(derivative, y0, t, args=(tmz_stim, ));
            
        for j in range(1, len(t)):
            
            DSB_total.append(y[j, 29])  # DSB value stored
            SSB_total.append(y[j, 30])  # SSB value stored
            
            t_total.append((i)*div + t[j]);
        
        t_total_hr = [i/(3600) for i in t_total];
        
        DSB_fun = InterpolatedUnivariateSpline(t_total_hr, DSB_total, k = 1);
        SSB_fun = InterpolatedUnivariateSpline(t_total_hr, SSB_total, k = 1);
        DSB_int = DSB_fun.integral(0, t_total[-1]);             # Integrated value of total Double strand breaks 
        SSB_int = SSB_fun.integral(0, t_total[-1]);             # Integrated value of total Single strand breaks 
        DNA_damage_integral = DSB_int + SSB_int;                # Integrated total DNA damage
        apoptosis_limit = 10000;                   # Apoptosis limit. If DNA damage exceed this apoptosis occurs.
        
        #print(mgmt, msh6, brca2)
        #print(s0[26],s0[27], s0[28])
        
        if  DNA_damage_integral > apoptosis_limit:   #Apoptosis loop
            
            cell_vector[0][k] = 0;
            t_cell_cycle = 0;
            
            if cell_vector[5][k] == 1:
                cell_count[0].append(-1)
                time_count[0].append(i*div/hr)
            elif cell_vector[5][k] == 2:
                cell_count[1].append(-1)
                time_count[1].append(i*div/hr) 
            elif cell_vector[5][k] == 3:
                cell_count[2].append(-1)
                time_count[2].append(i*div/hr)    
            elif cell_vector[5][k] == 4:
                cell_count[3].append(-1)
                time_count[3].append(i*div/hr)    
            elif cell_vector[5][k] == 5:
                cell_count[4].append(-1)
                time_count[4].append(i*div/hr)      
            
            print('Apoptosis happened')
            break
        
        brca2_ave = sum(y[:,26])/len(y[:,26]);
        msh6_ave  = sum(y[:,27])/len(y[:,27]);
        mgmt_ave  = sum(y[:,28])/len(y[:,28]);
        DSB_ave   = sum(y[:,29])/len(y[:,29]);
        SSB_ave   = sum(y[:,30])/len(y[:,30]);
        
        if  DSB_ave < 2 and SSB_ave < 2:   # Cell cycle arrest condition
            
            t_cell_cycle  = t_cell_cycle + div/hr;
       
        else:
            
            t_cell_cycle  = t_cell_cycle;       
            
      #  print(t_cell_cycle,'    ', DSB_ave, '    ', SSB_ave)  
        
        if t_cell_cycle > 24:      # Cell reproduction cycle
            
            if cell_vector[5][k] > 1: #Birth process all cells can reproduce except Q cells
            
                cell_vector[0].append(1);
                cell_vector[1].append(i);
                cell_vector[2].append(brca2_ave);
                cell_vector[3].append(msh6_ave);
                cell_vector[4].append(mgmt_ave);
                
            
                print('cell created', DSB_ave, SSB_ave, brca2)
                       

            if cell_vector[5][k] == 2:       #for stem cells there are 50% chance duagter cell is a stem cell or proliferating cell
                            
               u = rand.uniform(0, 1)
               a = np.zeros(2)
               b = [0.5, 0.5]
               a[0] = b[0];
               a[1] = a[0] + b[1];
               tot1 = a[1]
               if 0 < u < a[0]/tot1:                  # Daughter cell is GS-tr
                   cell_vector[5].append(2);          
                   cell_count[1].append(1)
                   time_count[1].append(i*div/hr)
               elif a[0]/tot1 < u < a[1]/tot1:        # Daughter cell si GP-tr
                   cell_vector[5].append(4);       
                   cell_count[3].append(1)
                   time_count[3].append(i*div/hr)
                   print('GStr --> GPtr')
               else:
                   print('error')
        
            elif cell_vector[5][k] == 3:
                
                u = rand.uniform(0, 1)
                a = np.zeros(2)
                b = [0.5, 0.5]
                a[0] = b[0];
                a[1] = a[0] + b[1];
                tot1 = a[1]
                
                if 0 < u < a[0]/tot1:                # Daughter cell is GS-tp
                   cell_vector[5].append(3);
                   cell_count[2].append(1)
                   time_count[2].append(i*div/hr)
                elif a[0]/tot1 < u < a[1]/tot1:      # Daughter cell is GP-tp
                   cell_vector[5].append(5);
                   cell_count[4].append(1)
                   time_count[4].append(i*div/hr)
                   print('GStp --> GPtp')
                else:
                   print('error')
               
                
            
            elif cell_vector[5][k] == 4:
                
                 cell_vector[5].append(4);
                 cell_count[3].append(1)
                 time_count[3].append(i*div/hr)    
                
                
            elif cell_vector[5][k] == 5:
                cell_vector[5].append(5);
                cell_count[4].append(1)
                time_count[4].append(i*div/hr)      
           
            t_cell_cycle = 0
            
            
        arr = [];
        for i in range(5):
            if i != cell_vector[5][k]:
                arr.append(i)
        
        
        ut = rand.uniform(0, 1)
        at = np.zeros(5)
        bt = [0.0001, 0.0001, 0.0001, 0.0001, 0.9996]; # Cell state transitions with 0.01% probability of transition
        at[0] = bt[0];
        at[1] = at[0] + bt[1];
        at[2] = at[1] + bt[2];
        at[3] = at[2] + bt[3];
        at[4] = at[3] + bt[4];
        tot_t = at[4]
        
        if 0 < ut < at[0]/tot_t: 
            cell_vector[5][k] == arr[0];  print('Cell Transistion')
        elif at[0]/tot_t < ut < at[1]/tot_t:
            cell_vector[5][k] == arr[1];  print('Cell Transistion')  
        elif at[1]/tot_t < ut < at[2]/tot_t:
            cell_vector[5][k] == arr[2];  print('Cell Transistion')  
        elif at[2]/tot_t < ut < at[3]/tot_t:
            cell_vector[5][k] == arr[3];  print('Cell Transistion')
        elif at[3]/tot_t < ut < at[4]/tot_t:
            cell_vector[5][k] == cell_vector[5][k];
        else:
            print('error')
            
        y0 = y[-1, :];   
        
    DNA_damage_integral_list.append(DNA_damage_integral) 
    #print(len(cell_vector));
    #print(k);    
    k = k + 1;    

# Re-count of cells for plotting 


cell_count_org = [[], [], [], [], []];
time_count_org = [[], [], [], [], []];
cell_number = [0, 0, 0, 0, 0]

for t in np.linspace(0, simulation_time_hr[-1], 1700):
    
    for i in range(5):
    
        for j in range(len(time_count[i])):
            if t <= time_count[i][j]<= t + 0.1:
                cell_number[i] = cell_number[i] + cell_count[i][j]
                cell_count_org[i].append(cell_number[i])    
                time_count_org[i].append(time_count[i][j])  
            else:
                cell_number[i] = cell_number[i]
                cell_count_org[i].append(cell_number[i]) 
                time_count_org[i].append(t)
       
plt.plot(time_count_org[0], cell_count_org[0], label="Q cells")  
plt.plot(time_count_org[1], cell_count_org[1], label="GStr cells")  
plt.plot(time_count_org[2], cell_count_org[2], label="GStp cells")  
plt.plot(time_count_org[3], cell_count_org[3], label="GPtr cells")  
plt.plot(time_count_org[4], cell_count_org[4], label="GPtp cells")  

plt.ylabel('Cell number')
plt.xlabel('Time (Hours)')
plt.legend(loc="upper right")
#plt.ylim(0, 350);
plt.xlim(0, simulation_time_hr[-1]);
plt.show()    
