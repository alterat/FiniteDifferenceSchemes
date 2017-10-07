# Simple harmonic oscillator simulation in Julia
# @INPUT:
#     f = frequency of the oscillator
#     σ = loss coefficients
# @OUTPUT:
#     out = time history of the oscillator's displacement
# AUTHOR: Alberto Torin
# =================================

function sho(f, σ)

# Initial conditions
x0=1
v0=0

SR=44100      # Sampling rate
T=1           # Duration of the simulation (s)

# Calculate the time step
k=1/SR
# Calculate the number of iterations
NF=floor(SR*T)

# Calculate ω0 and a0. They will be useful later
ω0=2*pi*f
a0=σ * k

# Calculate the recursion coefficients
B = (2-ω0^2*k^2)/(1+a0)
C = (1 - a0)/ (1 + a0)

# Check the stability condition
if k>(2/ω0)
    println("Stability condition violated!")
    return
end

# Initialise and allocate
u0=x0
u1=x0+k*v0
out=zeros(NF)
out[1]=u0
out[2]=u1

# Start the Main Loop
for nn=3:NF
  # Update the recursion
  u=B*u1-C*u0

  # Store the new value for the position
  out[nn]=u

  # Reset the known values of the recursion
  u0, u1 = u1, u;

end

  # At the end of the loop, return output
  return out
end
