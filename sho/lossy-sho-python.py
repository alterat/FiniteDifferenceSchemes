"""
Simulation of a simple harmonic oscillator using the finite difference method. 
INPUT: f, sigma
	f: the frequency of the oscillator (in Hz), bigger than 0.
	sigma: the loss coefficient, a small positive number (e.g., sigma=0.1)
OUTPUT: out_u
	a list with the values of the oscillator's position at each time step.
"""
import numpy as np

def sho(f,sigma):
    # Set sampling rate
    SR = 44100
    # Set duration of simulation (in s)
    T = 1

    # Time step
    k = 1/SR
    # Number of time steps
    NF = int(np.floor(SR*T))
   
    # Frequency and loss parameter
    omega = 2*np.pi*np.abs(f)
    a0 = sigma * k
    
     # Check the stability condition    
    if k>(2/omega):
        print("Stability condition violated!")
        return
    
    # Initial conditions
    x0 = 1
    v0 = 0
    
    # Auxiliary coefficients
    B = (2-(omega**2 * k**2)) / (1 + a0)
    C = (1 - a0) / (1 + a0)
    
    # Initialise the system
    u0=x0
    u1=x0+k*v0

    # Initialise output
    out_u = []

    # Main Loop - Update
    for ii in range(1, NF):
        u = B*u1 - C*u0
        out_u.append(u)
        u0, u1 = u1, u
    
    return out_u

