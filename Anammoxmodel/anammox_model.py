# This code solves the anammox model based on a selection of equations obtained from the
# following sources:
# FILL 

import math
import numpy as np
import matplotlib.pyplot as plt

#parameters
Y_H = 0.52
Y_NH = 0.15
Y_NO = 0.041
Y_AN = 0.159
Y_HNO2 = 0.44
Y_HNO3 = 0.44
i_NBM = 0.07
f_p = 0.08
i_nxi = 0.02
f_i = 0.1
k_H = 3.0
mu_H_max = 6.00#8.72
mu_NH_max = 0.26#0.8 # of 2.02
mu_NO_max = 0.79#0.79 # of 1.36
mu_AN_max = 0.1#0.019
K_X = 0.03
K_O_H = 0.2
K_S_H = 50 # of 20
K_NO3_H = 1.0
K_NO2_H = 1.0
K_pH_NH = 8.21
K_O_NH = 0.6 # of 0.235
K_NH3_NH = 0.75 # of 0.85
K_pH_NO = 8.66
K_O_NO = 1.5
K_HNO2_NO = 0.0008723
K_pH_AN = 12.41
K_NH_AN = 0.07
K_TNO2_AN = 0.05
K_O_AN = 0.01
b_H = 0.62#2.32
b_NO = 0.033 # of 0.092
b_NH = 0.05 # of 0.19
b_AN = 0.0025
theta_H = 0.069
theta_NO = 0.061
theta_NH = 0.094
theta_AN = 0.096
eta_NO2 = 0.6
eta_NO3 = 0.6
pH_opt_NH = 7.23
pH_opt_NO = 7.85
pH_opt_AN = 7.65
Temp_r = 20

def K_e_NH(Temp):
    K_e_NH = math.exp(-6344.0/(Temp+273.15))
    return K_e_NH

def K_e_NO(Temp):
    K_e_NO = math.exp(-2300.0/(Temp+273.15))
    return K_e_NO

def S_O(p_2,p_6,p_8):
    s_o = p_2*(Y_H-1)/Y_H+p_6*(Y_NH-3.43)/Y_NH+p_8*(Y_NO-1.14)/Y_NO
    return s_o

def S_S(p_1,p_2,p_4,p_5):
    s_s = p_1-p_2/Y_H-p_4/Y_HNO3-p_5/Y_HNO2
    return s_s

def S_NH(p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11):
    s_nh = -i_NBM*(p_2+p_4+p_5+p_8)+(i_NBM-f_p*i_nxi)*(p_3+p_7+p_9+p_11)+(-1/Y_NH-i_NBM)*p_6+(-1/Y_AN-i_NBM)*p_10
    return s_nh

def S_NH3(s_nh,pH,K_e_NH):
    s_nh3 = s_nh/(1+(10**(-pH))/K_e_NH)
    return s_nh3

def S_TNO2(p_4,p_5,p_6,p_8,p_10):
    s_tno2 = p_4*(1-Y_HNO3)/(1.14*Y_HNO3)+p_5*(Y_HNO2-1)/(1.71*Y_HNO2)+p_6/Y_NH-p_8/Y_NO+p_10*(-1.52-1/Y_AN)
    return s_tno2

def S_HNO2(s_tno2,pH,K_e_NO):
    s_hno2 = s_tno2/(1+K_e_NO/(10**(-pH)))
    return s_hno2

def S_TNO3(p_4,p_8,p_10):
    s_tno3 = p_4*(Y_HNO3-1)/(1.14*Y_HNO3)+p_8/Y_NO+1.52*p_10
    return s_tno3

def S_N2(p_5,p_10):
    s_n2 = p_5*(1-Y_HNO2)/(1.71*Y_HNO2)+p_10*2/Y_AN
    return s_n2

def X_H(p_2,p_3,p_4,p_5):
    x_h = p_2-p_3+p_4+p_5
    return x_h

def X_NH(p_6,p_7):
    x_nh = p_6-p_7
    return x_nh

def X_NO(p_8,p_9):
    x_no = p_8-p_9
    return x_no

def X_S(p_1,p_3,p_7,p_9,p_11):
    x_s = -p_1 + (1-f_i)*(p_3+p_7+p_9+p_11)
    return x_s

def X_I(p_3,p_7,p_9,p_11):
    x_i = f_i*(p_3+p_7+p_9+p_11)
    return x_i

def X_AN(p_10,p_11):
    x_an = p_10-p_11
    return x_an

def P_1(x_s,x_h):
    if x_h>0:
        p_1 = k_H*(x_s/x_h)/(K_X+x_s/x_h)*x_h
    else:
        p_1 = 0.0
    return p_1

def P_2(s_o,s_s,x_h,Temp):
    p_2 = mu_H_max*math.exp(theta_H*(Temp-Temp_r))*s_o/(K_O_H+s_o)*s_s/(K_S_H+s_s)*x_h
    return p_2

def P_3(x_h,Temp):
    p_3 = b_H*math.exp(theta_H*(Temp-Temp_r))*x_h
    return p_3

def P_4(s_o,s_tno3,s_s,x_h,Temp):
    p_4 = mu_H_max*math.exp(theta_H*(Temp-Temp_r))*eta_NO3*K_O_H/(K_O_H+s_o)*s_tno3/(K_NO3_H+s_tno3)*s_s/(K_S_H+s_s)*x_h
    return p_4

# TODO Wil deze functie TNO2 of NO2 (zonder HNO2)?
def P_5(s_o,s_tno2,s_s,x_h,Temp):
    p_5 = mu_H_max*math.exp(theta_H*(Temp-Temp_r))*eta_NO2*K_O_H/(K_O_H+s_o)*s_tno2/(K_NO2_H+s_tno2)*s_s/(K_S_H+s_s)*x_h
    return p_5

def P_6(s_o,s_nh3,x_nh,Temp,pH):
    p_6 = mu_NH_max*math.exp(theta_NH*(Temp-Temp_r))*K_pH_NH/(K_pH_NH-1+10**abs(pH_opt_NH-pH))*s_o/(K_O_NH+s_o)*s_nh3/(K_NH3_NH+s_nh3)*x_nh
    return p_6

def P_7(x_nh):
    p_7 = b_NH*math.exp(theta_NH*(Temp-Temp_r))*x_nh
    return p_7

def P_8(s_o,s_hno2,x_no):
    p_8 = mu_NO_max*math.exp(theta_NO*(Temp-Temp_r))*K_pH_NO/(K_pH_NO-1+10**abs(pH_opt_NO-pH))*s_o/(K_O_NO+s_o)*s_hno2/(K_HNO2_NO+s_hno2)*x_no
    return p_8    

def P_9(x_no,Temp):
    p_9 = b_NO*math.exp(theta_NO*(Temp-Temp_r))*x_no
    return p_9

def P_10(s_nh,s_tno2,s_o,x_an,Temp,pH):
    p_10 = mu_AN_max*math.exp(theta_AN*(Temp-Temp_r))*K_pH_AN/(K_pH_AN-1+10**abs(pH_opt_AN-pH))*s_nh/(K_NH_AN+s_nh)*s_tno2/(K_TNO2_AN+s_tno2)*K_O_AN/(K_O_AN+s_o)*x_an    
    return p_10

def P_11(x_an,Temp):
    p_11 = b_AN*math.exp(theta_AN*(Temp-Temp_r))*x_an
    return p_11

if __name__ == "__main__":
    # Influent information
    Temp = 32.3
    pH = 7.0
    NH_aanvoer_c = 1266.0
    TNO2_aanvoer_c = 0
    SS_aanvoer_c = 246
    CZV_c = 1107
    Dagaanvoer = 212
    Batches = 4.25
    Perc_X_S = 0.05
    Q_aanvoer = Dagaanvoer/(24/Batches)

    # Fractionation feed
    K_e_NH = K_e_NH(Temp); K_e_NO = K_e_NO(Temp)
    NH_aanvoer_v = NH_aanvoer_c*Q_aanvoer/1000
    NH3_aanvoer_c = S_NH3(NH_aanvoer_c,pH,K_e_NH); NH3_aanvoer_v = NH3_aanvoer_c*Q_aanvoer/1000
    NH4_aanvoer_c = NH_aanvoer_c-NH3_aanvoer_c; NH4_aanvoer_v = NH4_aanvoer_c*Q_aanvoer/1000
    TNO2_aanvoer_v = TNO2_aanvoer_c*Q_aanvoer/1000
    HNO2_aanvoer_c = S_HNO2(TNO2_aanvoer_c,pH,K_e_NO); HNO2_aanvoer_v = HNO2_aanvoer_c*Q_aanvoer/1000
    NO2_aanvoer_c = TNO2_aanvoer_c-HNO2_aanvoer_c; NO2_aanvoer_v = NO2_aanvoer_c*Q_aanvoer/1000
    X_S_aanvoer_c = Perc_X_S*SS_aanvoer_c; X_S_aanvoer_v = X_S_aanvoer_c*Q_aanvoer/1000
    X_I_aanvoer_c = SS_aanvoer_c-X_S_aanvoer_c; X_I_aanvoer_v = X_I_aanvoer_c*Q_aanvoer/1000

    # Reactor information
    reactor_volume = 412
    nitrogen_load_reactor = reactor_volume/((NH_aanvoer_v+TNO2_aanvoer_v)*Batches)
    verblijftijd = reactor_volume/Dagaanvoer
    zuurstof = 0.2
    dsg_max = 4.5
    uitspoeling_slib = 0.0
    invangen = 0.0

    # Additional values
    steps = 100
    Td = Batches/(24*(steps-6))
    Dh = 1/Td

    ##############################################################################
    """Step process"""

    #initialize solution matrix
    rows = ['S_O','S_S','S_NH','S_NH3','S_NH4','S_TNO2',"S_HNO2",'S_NO2','S_TNO3','S_N2','X_H','X_NH','X_NO','X_S','X_I','X_AN']
    solution_c = np.zeros((len(rows),steps))
    solution_v = np.zeros((len(rows),steps))
    volumes = np.zeros(steps)

    feed_c = np.array([0.0,SS_aanvoer_c,NH_aanvoer_c,NH3_aanvoer_c,NH4_aanvoer_c,TNO2_aanvoer_c,HNO2_aanvoer_c,NO2_aanvoer_c,0.0,0.0,0.0,0.0,0.0,X_S_aanvoer_c,X_I_aanvoer_c,0.0])
    feed_v = feed_c*Q_aanvoer/1000
    volumes[0] = reactor_volume
    volumes[1] = volumes[0]+Q_aanvoer
    start = np.array([zuurstof,100.0,400.0,0.0,400.0,0.0,0.0,0.0,0.0,0.0,1000.0,100.0,10.0,1000.0,100.0,200.0])

    solution_c[:,0] = start

    iterations = 1000
    result_matrix = np.zeros((len(rows),iterations))

    spuislib_matrix = np.zeros((3,iterations))

    for i in range(iterations):

        if i >0:
            solution_c[:,0] = solution_c[:,-1]

        solution_v[:,0] = solution_c[:,0]*volumes[0]/1000
        solution_v[:,1] = solution_v[:,0]+feed_v
        solution_c[:,1] = solution_v[:,1]*1000/volumes[1]

        # Allready did two steps, and the final 4 are for effluent and sludge calculations
        for step in range(steps-6):
            # Input
            input_vector = solution_c[:,1+step]
        
            s_o = input_vector[0]; s_s = input_vector[1]; s_nh = input_vector[2]; s_nh3 = input_vector[3]; s_nh4 = input_vector[4]; s_tno2 = input_vector[5];
            s_hno2 = input_vector[6]; s_no2 =input_vector[7]; s_tno3 = input_vector[8]; s_n2 = input_vector[9]; x_h = input_vector[10]; x_nh = input_vector[11];
            x_no = input_vector[12]; x_s = input_vector[13]; x_i = input_vector[14]; x_an = input_vector[15]

            # New kinetics
            p_1 = P_1(x_s,x_h); p_2 = P_2(zuurstof,s_s,x_h,Temp); p_3 = P_3(x_h,Temp); p_4 = P_4(zuurstof,s_tno3,s_s,x_h,Temp); p_5 = P_5(zuurstof,s_tno2,s_s,x_h,Temp);
            p_6 = P_6(zuurstof,s_nh3,x_nh,Temp,pH); p_7 = P_7(x_nh); p_8 = P_8(zuurstof,s_hno2,x_no); p_9 = P_9(x_no,Temp); p_10 = P_10(s_nh,s_tno2,zuurstof,x_an,Temp,pH); 
            p_11 = P_11(x_an, Temp)
            
            # New concentrations
            solution_c[0,2+step] = zuurstof#s_o + Td*S_O(p_2,p_6,p_8);
            solution_c[1,2+step] = s_s + Td*S_S(p_1,p_2,p_4,p_5); solution_c[2,2+step] = s_nh + Td*S_NH(p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11);
            solution_c[3,2+step] = s_nh3 + Td*S_NH3(s_nh,pH,K_e_NH); solution_c[4,2+step] = solution_c[2,2+step]-solution_c[3,2+step]; solution_c[5,2+step] = s_tno2 + Td*S_TNO2(p_4,p_5,p_6,p_8,p_10);
            solution_c[6,2+step] = s_hno2 + Td*S_HNO2(s_tno2,pH,K_e_NO); solution_c[7,2+step] = solution_c[5,2+step] - solution_c[6,2+step]; solution_c[8,2+step] = s_tno3 + Td*S_TNO3(p_4,p_8,p_10)
            solution_c[9,2+step] = s_n2 + Td*S_N2(p_5,p_10); solution_c[10,2+step] = x_h + Td*X_H(p_2,p_3,p_4,p_5); solution_c[11,2+step] = x_nh + Td*X_NH(p_6,p_7); solution_c[12,2+step] = x_no + Td*X_NO(p_8,p_9); 
            solution_c[13,2+step] = x_s + Td*X_S(p_1,p_3,p_7,p_9,p_11); solution_c[14,2+step] = x_i + Td*X_I(p_3,p_7,p_9,p_11); solution_c[15,2+step] = x_an + Td*X_AN(p_10,p_11)
            # if i == iterations-1:
            #     print(i, step, "and",p_2,p_3,p_4,p_5,"and",Td*X_H(p_2,p_3,p_4,p_5))

            # New volume (I will change this when the fill time will become finite, now it's infinitessimal)
            volumes[step+2] = volumes[1]

            # New loads 
            solution_v[:,2+step] = solution_c[:,2+step]*volumes[step+2]/1000.0

            # if i==iterations-1:
            #     if step==steps-8:
            #         plt.plot(solution_c[-6])
            #         plt.show()
        # Test number of steps
        # print("Number of steps gone right?: ",solution_v[:,-5],solution_v[:,-4])

        # Sludge information
        x_h = solution_c[-6,-5]; x_nh = solution_c[-5,-5]; x_no = solution_c[-4,-5]; x_s = solution_c[-3,-5]; x_i = solution_c[-2,-5]; x_an = solution_c[-1,-5]
        dsg_reactor = ((x_h+x_no+x_nh+x_s+x_an)/1.6+x_i)/1000
        spuislib = 0
        if dsg_reactor>dsg_max:
            spuislib = (dsg_reactor-dsg_max)*reactor_volume

        spuislib_m3 = spuislib/dsg_reactor
        slibleeftijd = dsg_reactor*reactor_volume/spuislib/Batches

        spuislib_matrix[0,i] = dsg_reactor
        spuislib_matrix[1,i] = spuislib
        spuislib_matrix[2,i] = spuislib_m3

        # Effluent
        solution_v[:,-4] = solution_v[:,-5]*Q_aanvoer/volumes[-5]
        solution_v[-6:,-4] = solution_v[-6:,-4]*uitspoeling_slib
        volumes[-4] = Q_aanvoer
        solution_c[:,-4] = solution_v[:,-4]*1000/Q_aanvoer
        # Stap 5 na aflaat
        solution_v[:,-3] = solution_v[:,-5]-solution_v[:,-4]
        volumes[-3] = reactor_volume
        solution_c[:,-3] = solution_v[:,-3]*1000/volumes[-3]
        # Stap 6 spuislib
        solution_c[:,-2] = solution_c[:,-3]
        volumes[-2] = spuislib_m3
        solution_v[:,-2] = solution_c[:,-2]*volumes[-2]/1000
        # Stap 7 na spuislib
        solution_v[:,-1] = solution_v[:,-3] - solution_v[:,-2]
        volumes[-1] = volumes[-3]
        solution_c[:,-1] = solution_v[:,-1]*1000/volumes[-1]

        result_matrix[:,i] = solution_c[:,-1]
    
    plt.plot(spuislib_matrix[0],label = "dsg reactor")
    plt.plot(spuislib_matrix[1], label = "spuislib")
    plt.plot(spuislib_matrix[2], label = "spuislib_m3")
    plt.legend()
    plt.show()
    

    # Bij plotten: laat aanvoer, effluent en spuislib weg.
    fig, axes = plt.subplots(4,4,figsize=(20, 20))
    plt.subplots_adjust(hspace=0.3,wspace = 0.3)
    for index in range(len(rows)): #len(rows)
        #tank_plot = np.concatenate(([solution_c[index,0]],solution_c[index,2:-4],[solution_c[index,-3]],[solution_c[index,-1]]),axis=0)
        axes[math.floor(index/4),index%4].plot(result_matrix[index])
        axes[math.floor(index/4),index%4].set_title(rows[index])
    plt.show()

    

    #print(solution_c[:,0:5])
    #print(solution_c[:,-5:-1])
    # # Aeration
    # elevation = 0
    # stat_pressure =((29.92-elevation*0.00105)/29.92*14.7)
    # sat_vapor_pressure  =(10**(-2238/(273.15+Temp)+8.896)*0.01934)
    # diff_depth = 14
    # sidewater = 16
    # C_20 = 9.092*(stat_pressure+0.007*62.4*diff_depth*0.25-sat_vapor_pressure)/(14.7-sat_vapor_pressure)
    # Omega =(stat_pressure+0.007*62.4*sidewater*0.25-sat_vapor_pressure)/(14.7+0.007*62.4*sidewater*0.25-sat_vapor_pressure) 
    # tau = math.exp((-3.9/0.00198*(1/293-1/(273+Temp))))*(stat_pressure-sat_vapor_pressure)/(stat_pressure-0.339))
    # beta = 0.95
    # C_inf = C_20*Omega*tau*beta
    # print("Check aeration: stat pres = %f, sat vapor pres = %f, C_20 = %f, Omega = %f, tau = %f, C_inf = %f"%(stat_pressure,sat_vapor_pressure,C_20,Omega,tau,C_inf))

    # kLa_dirty = S_O/24 # S_O of solution vector, Ronnie chooses the maximum value of the 4 steps
    # alfa = 0.7
    # O2_max = 483.7/(32.71+Temp)
    # kLa_clean = kLa_vuil/(1.025**(Temp-Temp_r)*alfa)
    # OC_clean = kLa_clean*O2_max*(reactor_volume+Q_aanvoer)/1000
    # O2_capacity = OC_clean*24/((NH_aanvoer_v+TNO2_aanvoer_v)*Batches)

