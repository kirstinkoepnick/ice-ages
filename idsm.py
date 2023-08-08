import numpy as np
from milankovitch import milankovitch
from scipy import *
from scipy.integrate import solve_ivp
from noise_term import noise_term#, create_nus

def get_eq_line(V, pt1, pt2):
    '''
    Input: ice volume (V) and points (pt1 and pt2) which are of the form
        pt1 = (x1, y1)
    Output: equation of the line between the two points pt1 and pt2
    '''
    x1, y1 = pt1
    x2, y2 = pt2
    slope = ((y2-y1)/(x2-x1))
    intercept = -slope*(x2)+y2
    return slope*V + intercept
    
def f(V, s, want_shape = False):
    '''
    Function that returns the piece-wise linear RHS of the IDSM without noise and milankovitch forcing
    Inputs: 
        V: land ice volume
        s: constant solar forcing (this changes location of the RHS--used to solve for the hystersis)
        want_shape: set to True if you want simply the shape of RHS with s = 0 (i.e. no solar constant)
    '''
    #returns the RHS without noise and milankovitch
    YEAR = 31536000
    kyr = YEAR*1000
    v0 = 0 
    v1 = 0.02 
    v2 = 0.175
    v3 = 0.25
    v4 = 0.7 
    v5 = 0.92
    v6 = 1 
    vdot0 = 0.005*kyr
    vdot1 = -0.0005*kyr
    vdot2 = 0.0031*kyr
    vdot3 = 0.0025*kyr 
    vdot4 = 0.00125*kyr
    vdot5 = -0.00125*kyr
    vdot6 = -0.004970001*kyr

    #creating points that establish the shape of the RHS plot
    
    point1 = [v0, vdot0] 
    point2 = [v1, vdot1]
    point3 = [v2, vdot2]
    point4 = [v3, vdot3]
    point5 = [v4, vdot4]
    point6 = [v5, vdot5]
    point7 = [v6, vdot6]
    
    if V <= point2[0]:
        line = get_eq_line(V, point1, point2)
    elif V <= point3[0] and V >= point2[0]:
        line = get_eq_line(V, point2, point3)
    elif V > point3[0] and V <= point4[0]:
        line = get_eq_line(V, point3, point4)
    elif V > point4[0] and V <= point5[0]:
        line = get_eq_line(V, point4, point5)
    elif V > point5[0] and V <= point6[0]:
        line = get_eq_line(V, point5, point6)
    else:
        line = get_eq_line(V, point6, point7)
    
    if want_shape == True:
        return line/10e8
    
    else:
        return line + (1e-02)*(0.5 - s)*kyr


def rhs_ivp(t, V, noise_amp=0, nu_time=0, nus=0):
    '''
    Function that returns the RHS to be solved using solve_ivp
    Inputs:
        t: time (kyr)
        V: land ice volume (nondimensional)
        noise_amp: amplitude of noise term if including (if noise_amp != 0, then must include non-zero nu_time and nus)
        nu_time: time values of noise term if including
        nus: y values of noise term if including
    Outputs:
        dVdt: RHS of the ODE used for the IDSM (nondimensional)
    '''

    YEAR = 31536000 #seconds in a year

    tau = 10e8 / 10e8
    alpha = 6 
    beta = 5 
    gamma = 0.9 / 10e8
    
    scaled_solar = (1e-02)*(0.5 -  milankovitch(t, paillard_truncation = False))

    alpha_inv = 1/(alpha)
    tau_inv = 1/(tau)
    
    dVdt_test = f(V, 0, want_shape = True) + (gamma)*(scaled_solar)*YEAR*1000 

    if noise_amp != 0:
    
        if dVdt_test < 0:
            dVdt = (beta)*(dVdt_test) + noise_amp*noise_term(t, nu_time, nus)
        else:
            dVdt = alpha_inv*dVdt_test + noise_amp*noise_term(t, nu_time, nus)
    else:
        if dVdt_test < 0:
            dVdt = (beta)*(dVdt_test)
        else:
            dVdt = alpha_inv*dVdt_test
        
    return dVdt*tau_inv
    

def integrate_model(KYRS_TO_RUN, V0, noise_amp=0, nus=0, nu_time=0):
    '''
    Function that integrates the IDSM
    Inputs:
        KYRS_TO_RUN: positive run time in kyr
        V0: initial condition
        noise_amp: amplitude of noise term if including (if noise_amp != 0, then must include non zero nu_time and nus)
        nu_time: time values of noise term if including
        nus: y values of noise term if including
    Outputs:
        ice_vals: nondimensional ice volume spaced every 1 kyr
    '''
    time_start =  -KYRS_TO_RUN
    time_end = 0.0
    time_step = 1
    times_save = np.arange(time_start, time_end, time_step)
    tspan = (times_save[0], times_save[-1])
    
    if noise_amp != 0:
        sol = solve_ivp(rhs_ivp, t_span=tspan, y0=V0, args = (noise_amp, nus, nu_time), method='BDF', t_eval = times_save)
    else:
        sol = solve_ivp(rhs_ivp, t_span=tspan, y0=V0, method='BDF', t_eval = times_save)
    ice_vals = sol.y[0]
    return ice_vals
