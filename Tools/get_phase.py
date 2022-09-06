#retrieves phase in optics convention from omega and Eomega
import numpy as np
c=299792458

def get_phase(omega,Eomega, lambda0):
    #retrieves phase and subtracts the arbitrary phase at omega0
    om0 = 2*np.pi*c/lambda0
    om0_idx = np.argmin(np.abs(omega-om0))
 
    domega = omega[2] - omega[1]
    tau = np.pi/domega
    phase_raw = np.angle(Eomega)
    phase = np.unwrap(phase_raw - omega*tau)
    #phase -= phase[np.argmin(np.abs(omega - 2.4e-15))]
    phase=phase*(-1)
    phase -= phase[om0_idx]
    return phase