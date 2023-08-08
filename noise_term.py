import numpy as np

def create_nus():
    '''
    Function that outputs nus from equation 2 from Koepnick & Tziperman 2023
    The output of this function is interpolated to create a noise function
    '''
    YEAR = 31536000 ##year in seconds
    dt = 100/1000 #time steps

    n = int((1500000/1000)/dt) 
    nus = np.zeros(n)
    nu_time = np.zeros(n)
    R = np.exp(-dt/3) 

    sigma = np.sqrt(0.042)
    nu_n = 0.0 #start nu

    for i in range(n):
        theta = np.random.normal(0, sigma)
        nu_n_old = 1.0*nu_n
        nu_n = R*nu_n_old + np.sqrt(1-R**2)*theta
        nus[i] = nu_n 
        nu_time[i] = i*dt - 1500
    
    return nu_time, nus

def noise_term(t, nu_time, nus):
    '''
    Input: time (t) in units of kyrs; nu time series (nu_time and nus created from create_nus())
    Output: noise for that time t
    Note: the function create_nus() will need to be called before running this function
    to create the nu_time and nus as input variables for this function. 
    That function (create_nus()) creates the randomized stochastic forcing
    which is then interpolated for any desired time in this function (noise_term(t))
    '''
    nu_interpolated = np.interp(t, nu_time, nus)  
    
    return nu_interpolated