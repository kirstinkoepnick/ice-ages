import numpy as np
from milankovitch import milankovitch
from noise_term import noise_term#, create_nus

asi=0
def mod_sis(t, V, params, noise_amp=0, nu_time=0, nus=0):
    '''
    Function that 
    Input: 
        t: time (kyr)
        V: land ice volume nondimensional
        Other parameters (params): 
            precipitation rate (p0)
            growth rate (k)
            solar constanst (S0 and SM)
            ice volume threshold (V_max and V_min)
            sea ice switch (asi_on)
        Including stochastic forcing: 
            Use the function "create_nus" to create variables "nu_time" and "nus"
            Set desired noise amplitude with noise_amp
    Output: 
        dVdt: RHS of equation 3 from Koepnick & Tziperman 2023 (modified SIS with a threshold as a function of insolation)
    '''
    
    global asi
    p0, k, S0, SM, V_max, V_min, asi_on = params    
        
    alpha = 0.4 #threshold value for turning on the sea ice switch--for unmodified SIS: set to 0
    
    if V - alpha*milankovitch(t) > V_max:
        asi = asi_on - alpha*asi_on*milankovitch(t)

    elif V - alpha*milankovitch(t) < V_min:
        asi = 0

    if noise_amp != 0:
        #including stochastic/noise forcing
        dVdt = (p0-k*V)*(1-asi) - (S0) - SM*milankovitch(t, True) + noise_amp*noise_term(t, nu_time, nus) ##equation A2 from paper
    else:
        #not including stochastic/noise forcing
        dVdt = (p0-k*V)*(1-asi) - (S0) - SM*milankovitch(t, True)

    epsilon = V_min/10
    
    if V < epsilon and dVdt < 0:
        return dVdt*0.0 #makes sure that the ice volume does not go below zero
    else:
        return dVdt