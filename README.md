# Are the Pleistocene 100 kyr ice ages driven or phase-locked by Milankovitch forcing?

An investigation of glacial-interglacial oscillations and the role of insolation variations.

## Description

Glacial-interglacial oscillations exhibit a periodicity of approximately 100 kyr during the late Pleistocene. Insolation variations are generally understood to play a vital role in these ice ages, yet their exact role is still unknown; the 100 kyr ice ages may fall into either of two categories. First, they could be purely insolation-driven, such that ice ages are a consequence of insolation variations. Or second, ice ages are self-sustained oscillations, where they would have existed even without insolation variations. We develop several observable measures that are used to differentiate between the two scenarios and attempt to determine which one is more likely based on the observed proxy record. For this purpose, we use example mechanisms and simple models representing each of these two scenarios. For questions regarding the individual analyses/diagonistics included in Koepnick & Tziperman (2023; submitted), please contact Kirstin Koepnick at kirstinkoepnick@g.harvard.edu.

## Getting Started

### Data

* Berger (1978)'s calculation of insolation forcing for June 21st at 65 $`^\circ`$N is found in `Milankovitch.dat` file and interpolated for any desired time in `milankovitch.py`.
    * Data file has columns are time, latitude, insolation, eccentricity, precession, longitude, obliquity
* We use the untuned oxygen-18 compilation of Lisiecki & Raymo (2005, LR04) found in `LR04untuned.txt`. 

### Running the simple models

* The python functions are written to use `scipy.integrate.solve_ivp` but can be adapted for `scipy.integrate.ode_int`, if desired.
* To run the modified sea-ice switch, for example:
```
from scipy.integrate import solve_ivp
import numpy as np
from mod_sis import mod_sis
from sis_params import *

YEARS_TO_RUN_SIS = 1100000

V0, time_start = np.array([0.0]), -YEARS_TO_RUN_SIS/1000 #defining initial condition and time to start run in kyrs

time_end = 0 
dt = 1 #evaluation steps
times_save = np.arange(time_start, time_end, dt)
tspan = (times_save[0], times_save[-1])

method = 'BDF' #method best used for stiff odes
asi = 0 #sea-ice switch initial 

print('Solving Modified SIS...')
sol_sis = solve_ivp(mod_sis, args = (p,), t_span=tspan, y0=V0, t_eval = times_save, method = method,max_step = 0.1, rtol=1e-9, atol=1e-9)

time_vals_sis, ice_vals_sis = sol_sis.t, sol_sis.y[0]
print('Done.')
```
* To run the insolation-driven simple model (IDSM), simply use the function "integrate_model()" given in "idsm.py". For example:
```
from scipy.integrate import solve_ivp
import numpy as np
from idsm import * 
YEARS_TO_RUN = 800000
V0 = np.array([1])
print('Solving ISTM...')
ice_vals_idsm = integrate_model(YEARS_TO_RUN/1000, V0)
print('Done.')
```

## Authors

Kirstin Koepnick (kirstinkoepnick@g.harvard.edu)
Eli Tziperman
