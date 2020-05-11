"""" Deze code berekent de concentraties stoffen en biomassa in de DEMON-reactor in Echten.
    Dit betekent dat er een vultijd en behandeltijd in zit, en dat het slib alleen via
    uitspoeling en niet via spuien de tank verlaat. Bij een ander project zal dit spuien
    waarschijnlijk toegevoegd moeten worden.
    TODO spuien toevoegen
    TODO invoer maken van startconcentraties
""""


import math
import numpy as np
import matplotlib.pyplot as plt
import easygui
import tkinter
import csv

#parameters


def K_e_NH(Temp):
    K_e_NH = math.exp(-6344.0/(Temp+273.15))
    return K_e_NH

def K_e_NO(Temp):
    K_e_NO = math.exp(-2300.0/(Temp+273.15))
    return K_e_NO

def KLA(O_sat,zuurstof,p_2,p_6,p_8):
    kla =  -(p_2*(Y_H-1)/Y_H+p_6*(Y_NH-3.43)/Y_NH+p_8*(Y_NO-1.14)/Y_NO)/(O_sat-(zuurstof+0.1))
    return kla

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
    p_8 = mu_NO_max*math.exp(theta_NO*(Temp-Temp_r))*K_pH_NO/(K_pH_NO-1+10.0**abs(pH_opt_NO-pH))*s_o/(K_O_NO+s_o)*s_hno2/(K_HNO2_NO+s_hno2)*x_no
    return p_8    

def P_9(x_no,Temp):
    p_9 = b_NO*math.exp(theta_NO*(Temp-Temp_r))*x_no
    return p_9

def P_10(s_nh,s_tno2,s_o,x_an,Temp,pH):
    p_10 = mu_AN_max*math.exp(theta_AN*(Temp-Temp_r))*K_pH_AN/(K_pH_AN-1+10.0**abs(pH_opt_AN-pH))*s_nh/(K_NH_AN+s_nh)*s_tno2/(K_TNO2_AN+s_tno2)*K_O_AN/(K_O_AN+s_o)*x_an    
    return p_10

def P_11(x_an,Temp):
    p_11 = b_AN*math.exp(theta_AN*(Temp-Temp_r))*x_an
    return p_11

if __name__ == "__main__":
    print("Geef invoerbestand op.")
    input_file = easygui.fileopenbox(title="Please select an inputfile")

    print("Geef uitvoermap op.")
    root = tkinter.Tk()
    root.withdraw()
    output_path = tkinter.filedialog.askdirectory(parent=root,initialdir="/",title='Please select an output directory')

    try:
        # Read input and store in dictionary
        params = {}
        read_input = open(input_file)
        for line in read_input:
            line = line.strip()
            if not line.startswith("#"):
                key_value = line.split("=")
                if len(key_value) == 2:
                    params[key_value[0].strip()] = key_value[1].strip()

        print(params)
        # parameters
        Y_H = float(params["Y_H"])
        Y_NH = float(params["Y_NH"])
        Y_NO = float(params["Y_NO"])
        Y_AN = float(params["Y_AN"])
        Y_HNO2 = float(params["Y_HNO2"])
        Y_HNO3 = float(params["Y_HNO3"])
        i_NBM = float(params["i_NBM"])
        f_p = float(params["f_p"])
        i_nxi = float(params["i_nxi"])
        f_i = float(params["f_i"])
        k_H = float(params["k_H"])
        mu_H_max = float(params["mu_H_max"])
        mu_NH_max = float(params["mu_NH_max"])
        mu_NO_max = float(params["mu_NO_max"])
        mu_AN_max = float(params["mu_AN_max"])
        K_X = float(params["K_X"])
        K_O_H = float(params["K_O_H"])
        K_S_H = float(params["K_S_H"])
        K_NO3_H = float(params["K_NO3_H"])
        K_NO2_H = float(params["K_NO2_H"])
        K_pH_NH = float(params["K_pH_NH"])
        K_O_NH = float(params["K_O_NH"])
        K_NH3_NH = float(params["K_NH3_NH"])
        K_pH_NO = float(params["K_pH_NO"])
        K_O_NO = float(params["K_O_NO"])
        K_HNO2_NO = float(params["K_HNO2_NO"])
        K_pH_AN = float(params["K_pH_AN"])
        K_NH_AN = float(params["K_NH_AN"])
        K_TNO2_AN = float(params["K_TNO2_AN"])
        K_O_AN = float(params["K_O_AN"])
        b_H = float(params["b_H"])
        b_NO = float(params["b_NO"])
        b_NH = float(params["b_NH"])
        b_AN = float(params["b_AN"])
        theta_H = float(params["theta_H"])
        theta_NO = float(params["theta_NO"])
        theta_NH = float(params["theta_NH"])
        theta_AN = float(params["theta_AN"])
        eta_NO2 = float(params["eta_NO2"])
        eta_NO3 = float(params["eta_NO3"])
        pH_opt_NH = float(params["pH_opt_NH"])
        pH_opt_NO = float(params["pH_opt_NO"])
        pH_opt_AN = float(params["pH_opt_AN"])
        Temp_r = float(params["Temp_r"])

        print(K_X,K_O_AN,Temp_r)

        # Influent information
        Temp = float(params["Temp"])
        O_sat = 483.7/(32.71+Temp)
        pH = float(params["pH"])
        NH_aanvoer_c = float(params["NH_aanvoer_c"])
        TNO2_aanvoer_c = float(params["TNO2_aanvoer_c"])
        OB_aanvoer_c = float(params["OB_aanvoer_c"])
        CZV_aanvoer_c = float(params["CZV_aanvoer_c"])
        Dagaanvoer = float(params["Dagaanvoer"])
        Batches = float(params["Batches"])
        Perc_X_S = float(params["Perc_X_S"])
        Q_aanvoer = Dagaanvoer/(24/(Batches/60))
        Tijd_Q_aanvoer = float(params["Tijd_Q_aanvoer"])
        Tijd_bezinking = float(params["Tijd_bezinking"])
        Tijd_overig = Batches-Tijd_Q_aanvoer-Tijd_bezinking
        perc_CSV_S_s = float(params["perc_CSV_S_s"])

        # Fractionation feed
        K_e_NH = K_e_NH(Temp); K_e_NO = K_e_NO(Temp)
        NH_aanvoer_v = NH_aanvoer_c*Q_aanvoer/1000
        NH3_aanvoer_c = S_NH3(NH_aanvoer_c,pH,K_e_NH); NH3_aanvoer_v = NH3_aanvoer_c*Q_aanvoer/1000
        NH4_aanvoer_c = NH_aanvoer_c-NH3_aanvoer_c; NH4_aanvoer_v = NH4_aanvoer_c*Q_aanvoer/1000
        TNO2_aanvoer_v = TNO2_aanvoer_c*Q_aanvoer/1000
        HNO2_aanvoer_c = S_HNO2(TNO2_aanvoer_c,pH,K_e_NO); HNO2_aanvoer_v = HNO2_aanvoer_c*Q_aanvoer/1000
        NO2_aanvoer_c = TNO2_aanvoer_c-HNO2_aanvoer_c; NO2_aanvoer_v = NO2_aanvoer_c*Q_aanvoer/1000
        X_S_aanvoer_c = Perc_X_S*OB_aanvoer_c; X_S_aanvoer_v = X_S_aanvoer_c*Q_aanvoer/1000
        X_I_aanvoer_c = OB_aanvoer_c-X_S_aanvoer_c; X_I_aanvoer_v = X_I_aanvoer_c*Q_aanvoer/1000
        S_s_aanvoer_c = perc_CSV_S_s*CZV_aanvoer_c; S_s_aanvoer_v = S_s_aanvoer_c*Q_aanvoer/1000

        # Reactor information
        # PAS AAN TODO
        reactor_volume_totaal = float(params["reactor_volume"])
        reactor_volume = reactor_volume_totaal - Q_aanvoer
        nitrogen_load_reactor = reactor_volume/((NH_aanvoer_v+TNO2_aanvoer_v)*Batches)
        verblijftijd = reactor_volume/Dagaanvoer
        zuurstof = float(params["zuurstof"])
        alfa = float(params["alfa"])
        dsg_max = float(params["dsg_max"])
        uitspoeling_slib = float(params["uitspoeling_slib"])
        invangen = float(params["invangen"])

        # Additional values
        steps_per_minute = int(params["steps_per_minute"])
        Td = 1.0/(24*60*float(steps_per_minute))
        steps_filling = int(Tijd_Q_aanvoer*steps_per_minute)
        steps_overig = int(Tijd_overig*steps_per_minute)
        steps = steps_filling + steps_overig
        if Q_aanvoer>0:
            Q_aanvoer_per_step = Q_aanvoer/steps_filling
        else:
            Q_aanvoer_per_step = 0
        iterations = int(params["iterations"])

        # Effluent information       
        NH4_effluent_c = float(params["NH4_effluent_c"])
        NH_effluent_c = NH4_effluent_c*(1+K_e_NH/10**(-pH))
        NH3_effluent_c = NH_effluent_c - NH4_effluent_c
        CZV_effluent_c = float(params["CZV_effluent_c"])
        S_s_effluent_c = CZV_effluent_c*perc_CSV_S_s
        OB_effluent_c = float(params["OB_effluent_c"])
        X_S_effluent_c = Perc_X_S*OB_effluent_c
        X_I_effluent_c = OB_effluent_c - X_S_effluent_c
        TNO2_effluent_c = float(params["TNO2_effluent_c"]); HNO2_effluent_c = S_HNO2(TNO2_effluent_c,pH,K_e_NO); NO2_effluent_c = TNO2_effluent_c - HNO2_effluent_c
        TNO3_effluent_c = float(params["TNO3_effluent_c"])

        effluent_c = [zuurstof,S_s_effluent_c,NH_effluent_c,NH3_effluent_c,NH4_effluent_c,TNO2_effluent_c,HNO2_effluent_c,NO2_effluent_c,TNO3_effluent_c,0,0,0,0,0,X_S_effluent_c,X_I_effluent_c*invangen,0]
        
        print("Parameters zijn ingelezen.")

        ##############################################################################
        """Step process"""

        #initialize solution matrix
        rows = ['S_O','S_S','S_NH','S_NH3','S_NH4','S_TNO2',"S_HNO2",'S_NO2','S_TNO3','S_N2','X_H','X_NH','X_NO','X_S','X_I','X_AN']
        solution_c = np.zeros((len(rows),steps+3)) #klopt 4?
        solution_v = np.zeros((len(rows),steps+3))
        volumes = np.zeros(steps+3)

        feed_c = np.array([0.0,S_s_aanvoer_c,NH_aanvoer_c,NH3_aanvoer_c,NH4_aanvoer_c,TNO2_aanvoer_c,HNO2_aanvoer_c,NO2_aanvoer_c,0.0,0.0,0.0,0.0,0.0,X_S_aanvoer_c,X_I_aanvoer_c*invangen,0.0])
        feed_v = feed_c*Q_aanvoer_per_step/1000
        volumes[0] = reactor_volume

        start = np.array([100.0,28.0,357.0,3.3,354.0,16.0,0.02,16.0,51.0,800.0,35.0,220.0,1.0,500.0,4000.0,2000.0])

        solution_c[:,0] = start

        result_matrix = np.zeros((len(rows),iterations))
        result_p = np.zeros((11,iterations))
        result_vracht = np.zeros(iterations)

        print("Iteraties worden doorlopen ...")

        for i in range(iterations):
            if i >0:
                solution_c[:,0] = solution_c[:,-1]

            solution_v[:,0] = solution_c[:,0]*volumes[0]/1000

            for step in range(steps_filling):
                # Input
                input_vector = solution_c[:,step]
            
                s_o = input_vector[0]; s_s = input_vector[1]; s_nh = input_vector[2]; s_nh3 = input_vector[3]; s_nh4 = input_vector[4]; s_tno2 = input_vector[5];
                s_hno2 = input_vector[6]; s_no2 =input_vector[7]; s_tno3 = input_vector[8]; s_n2 = input_vector[9]; x_h = input_vector[10]; x_nh = input_vector[11];
                x_no = input_vector[12]; x_s = input_vector[13]; x_i = input_vector[14]; x_an = input_vector[15]

                # New kinetics
                p_1 = P_1(x_s,x_h); p_2 = P_2(zuurstof,s_s,x_h,Temp); p_3 = P_3(x_h,Temp); p_4 = P_4(zuurstof,s_tno3,s_s,x_h,Temp); p_5 = P_5(zuurstof,s_tno2,s_s,x_h,Temp);
                p_6 = P_6(zuurstof,s_nh3,x_nh,Temp,pH); p_7 = P_7(x_nh); p_8 = P_8(zuurstof,s_hno2,x_no); p_9 = P_9(x_no,Temp); p_10 = P_10(s_nh,s_tno2,zuurstof,x_an,Temp,pH); 
                p_11 = P_11(x_an, Temp)

                # New concentrations
                kla = KLA(O_sat,zuurstof,p_2,p_6,p_8); kla_vies = kla/24.0;
                kla_schoon = (kla_vies/1.025**(Temp-Temp_r))/alfa
                OC_schoon = kla_schoon*O_sat
                solution_c[0,1+step] = OC_schoon
                solution_c[1,1+step] = max(s_s + Td*S_S(p_1,p_2,p_4,p_5),0.0); solution_c[2,1+step] = max(0.0,s_nh + Td*S_NH(p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11));
                solution_c[3,1+step] = max(0.0,S_NH3(s_nh,pH,K_e_NH)); solution_c[4,1+step] = max(0.0,solution_c[2,1+step]-solution_c[3,1+step]); solution_c[5,1+step] = max(0.0,s_tno2 + Td*S_TNO2(p_4,p_5,p_6,p_8,p_10));
                solution_c[6,1+step] = max(0.0,s_hno2 + Td*S_HNO2(s_tno2,pH,K_e_NO)); solution_c[7,1+step] = max(0.0,solution_c[5,1+step] - solution_c[6,1+step]); solution_c[8,1+step] = max(0.0,s_tno3 + Td*S_TNO3(p_4,p_8,p_10))
                solution_c[9,1+step] = max(0.0,s_n2 + Td*S_N2(p_5,p_10)); solution_c[10,1+step] = max(0.0,x_h + Td*X_H(p_2,p_3,p_4,p_5)); solution_c[11,1+step] = max(0.0,x_nh + Td*X_NH(p_6,p_7)); solution_c[12,1+step] = max(0.0,x_no + Td*X_NO(p_8,p_9)); 
                solution_c[13,1+step] = max(0.0,x_s + Td*X_S(p_1,p_3,p_7,p_9,p_11)); solution_c[14,1+step] = max(0.0,x_i + Td*X_I(p_3,p_7,p_9,p_11)); solution_c[15,1+step] = max(0.0,x_an + Td*X_AN(p_10,p_11))

                solution_v[:,1+step] = solution_c[:,1+step]*volumes[step]/1000.0 + feed_v
                volumes[step+1] = volumes[step] + Q_aanvoer_per_step
                solution_c[:,1+step] = solution_v[:,1+step]*1000.0/volumes[step+1]

            for step in range(steps_filling,steps):
                # Input
                input_vector = solution_c[:,step]
            
                s_o = input_vector[0]; s_s = input_vector[1]; s_nh = input_vector[2]; s_nh3 = input_vector[3]; s_nh4 = input_vector[4]; s_tno2 = input_vector[5];
                s_hno2 = input_vector[6]; s_no2 = input_vector[7]; s_tno3 = input_vector[8]; s_n2 = input_vector[9]; x_h = input_vector[10]; x_nh = input_vector[11];
                x_no = input_vector[12]; x_s = input_vector[13]; x_i = input_vector[14]; x_an = input_vector[15]

                # New kinetics
                p_1 = P_1(x_s,x_h); p_2 = P_2(zuurstof,s_s,x_h,Temp); p_3 = P_3(x_h,Temp); p_4 = P_4(zuurstof,s_tno3,s_s,x_h,Temp); p_5 = P_5(zuurstof,s_tno2,s_s,x_h,Temp);
                p_6 = P_6(zuurstof,s_nh3,x_nh,Temp,pH); p_7 = P_7(x_nh); p_8 = P_8(zuurstof,s_hno2,x_no); p_9 = P_9(x_no,Temp); p_10 = P_10(s_nh,s_tno2,zuurstof,x_an,Temp,pH); 
                p_11 = P_11(x_an, Temp)
                
                # New concentrations
                solution_c[0,1+step] = KLA(O_sat,zuurstof,p_2,p_6,p_8);
                solution_c[1,1+step] = max(0.0,s_s + Td*S_S(p_1,p_2,p_4,p_5)); solution_c[2,1+step] = max(0.0,s_nh + Td*S_NH(p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11));
                solution_c[3,1+step] = max(0.0,S_NH3(s_nh,pH,K_e_NH)); solution_c[4,1+step] = max(0.0,solution_c[2,1+step]-solution_c[3,1+step]); solution_c[5,1+step] = max(0.0,s_tno2 + Td*S_TNO2(p_4,p_5,p_6,p_8,p_10));
                solution_c[6,1+step] = max(0.0,s_hno2 + Td*S_HNO2(s_tno2,pH,K_e_NO)); solution_c[7,1+step] = max(0.0,solution_c[5,1+step] - solution_c[6,1+step]); solution_c[8,1+step] = max(0.0,s_tno3 + Td*S_TNO3(p_4,p_8,p_10));
                solution_c[9,1+step] = max(0.0,s_n2 + Td*S_N2(p_5,p_10)); solution_c[10,1+step] = max(0.0,x_h + Td*X_H(p_2,p_3,p_4,p_5)); solution_c[11,1+step] = max(0.0,x_nh + Td*X_NH(p_6,p_7)); solution_c[12,1+step] = max(0.0,x_no + Td*X_NO(p_8,p_9)); 
                solution_c[13,1+step] = max(0.0,x_s + Td*X_S(p_1,p_3,p_7,p_9,p_11)); solution_c[14,1+step] = max(0.0,x_i + Td*X_I(p_3,p_7,p_9,p_11)); solution_c[15,1+step] = max(0.0,x_an + Td*X_AN(p_10,p_11))

                volumes[step+1] = volumes[step]
                solution_v[:,1+step] = solution_c[:,1+step]*volumes[step+1]/1000.0

            result_p[:,i] = [p_1,p_2,p_3,p_4,p_5,p_6,p_7,p_8,p_9,p_10,p_11]
            
            # Effluent
            x_h = solution_c[-6,-3]; x_nh = solution_c[-5,-3]; x_no = solution_c[-4,-3]; x_s = solution_c[-3,-3]; x_i = solution_c[-2,-3]; x_an = solution_c[-1,-3]
            dsg_reactor = ((x_h+x_no+x_nh+x_s+x_an)/1.6+x_i)/1000
            dsg_naeffluent = dsg_reactor*reactor_volume_totaal/reactor_volume
            # laat zoveel uitspoelen dat er 5,5 g/l overblijft
            dsg_uitspoelen_c = max(dsg_naeffluent-dsg_max,0.0)
            factor_uitspoelen = (dsg_uitspoelen_c*reactor_volume)/(dsg_reactor*reactor_volume_totaal)

            slibleeftijd = ((dsg_reactor*reactor_volume_totaal)/(dsg_uitspoelen_c*reactor_volume))/(24*60/Batches)

            volumes[-2] = Q_aanvoer
            solution_v[:,-2] = solution_v[:,-3]*volumes[-2]/volumes[-3]
            solution_v[-6:,-2] = solution_v[-6:,-3]*factor_uitspoelen
            solution_c[:,-2] = solution_v[:,-2]*1000/volumes[-2]
            
            # Stap na aflaat
            volumes[-1] = volumes[-3]-volumes[-2]
            solution_v[:,-1] = solution_v[:,-3]-solution_v[:,-2]
            solution_c[:,-1] = solution_v[:,-1]*1000/volumes[-1]
            result_matrix[:,i] = solution_c[:,-1]

            #check dsg
            x_h = solution_c[-6,-1]; x_nh = solution_c[-5,-1]; x_no = solution_c[-4,-1]; x_s = solution_c[-3,-1]; x_i = solution_c[-2,-1]; x_an = solution_c[-1,-1]
            dsg_reactor_na = ((x_h+x_no+x_nh+x_s+x_an)/1.6+x_i)/1000

            result_vracht[i] = solution_c[0,-1]*(reactor_volume+Q_aanvoer)/1000
     
        #print("Het dgs, gewicht aan spuislib, volume en de slibleeftijd zijn: %f, %f, %f, %f"%(dsg_reactor, spuislib, volumes[-3],slibleeftijd))
        print("Slibleeftijd, dsg_na_effluent, dsg_effl, dsg_na: %f, %f, %f, %f"%(slibleeftijd,dsg_reactor_na,(dsg_uitspoelen_c*reactor_volume/Q_aanvoer), dsg_naeffluent))
        print("Alle iteraties zijn doorlopen.")
        print("Plots wegschrijven ...")
        plott = True
        if plott == True:
            # plt.plot(spuislib_matrix[0],label = "dsg reactor")
            # plt.plot(spuislib_matrix[1], label = "spuislib")
            # plt.plot(spuislib_matrix[2], label = "spuislib_m3")
            # plt.legend()
            # output_dgs = str(output_path)+"\drogestofgehaltes.png"
            # plt.savefig(output_dgs)
            # plt.close()

            figp, axesp = plt.subplots(4,3,figsize=(20,20))
            plt.subplots_adjust(hspace=0.3,wspace = 0.3)
            for indexp in range(11): #len(rows)
                axesp[math.floor(indexp/3),indexp%3].plot(result_p[indexp])
                axesp[math.floor(indexp/3),indexp%3].set_title("p_%i"%(indexp+1))
            output_kineticts = output_path + "\kinetische_factoren.png"
            plt.savefig(output_kineticts)
            plt.close()

            fig, axes = plt.subplots(4,4,figsize=(20, 20))
            plt.subplots_adjust(hspace=0.3,wspace = 0.3)
            for index in range(len(rows)): #len(rows)
                if effluent_c[index] != 0:
                    axes[math.floor(index/4),index%4].hlines(y=effluent_c[index],xmin=0,xmax=iterations, color='r', linestyle='-')
                    #axes[math.floor(index/4),index%4].plot(effluent_c[index])
                axes[math.floor(index/4),index%4].plot(result_matrix[index])
                axes[math.floor(index/4),index%4].set_title(rows[index])
            output_concentrations = output_path + "\concentraties.png"
            plt.savefig(output_concentrations)
            plt.close()

            plt.plot(result_vracht)
            output_zuurstof = output_path + "\zuurstofvracht.png"
            plt.savefig(output_zuurstof)
            plt.close()
        print("Alle plots zijn weggeschreven.")

        
        #zuurstof_inbreng = OC_schoon*24/(NH_aanvoer_v+TNO2_aanvoer_v)

        print("Tabel wegschrijven ...")
        my_table = {}
        for tab in range(len(rows)):
            my_table[rows[tab]] = result_matrix[tab,-1]
        my_table["OC schoon"] = OC_schoon
        my_table["slibleeftijd"] = slibleeftijd
        output_table = str(output_path) + "\\tabel.csv"
        with open(output_table, 'w') as f:
            for key in my_table.keys():
                f.write("%s;%s\n"%(key,my_table[key]))
        print("Tabel weggeschreven")        



    except Exception:
        import traceback
        traceback.print_exc()
        input("Program crashed; press Enter to exit")
        
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


