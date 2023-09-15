# Nickelocene as Spin Sensor

Download these two files and run the code in the same directory to use the code.

``` python3 Spin_sensing_Nc.py parameters.txt ```

You will get an image.

In the parameter file, you can edit the input parameters. 
All the input parameters are: 

```
M_spin = 3/2           # Spin of Metal.
D_aniso = [0,0,5]      # Magneitc anisotropy of Metal.
Temp = 0.3             # Temperature.
B_field = 0            # Magnetic Field.
J_Nc = 3               # Magnetic Exchange between metal and Nc.
Nsites = 1             # Number of spin sites.
J1 = 6                 # Nearest Neighbour exchange when Nsites are greater than 1.
J2  = 0                # second Neighbour exchange when Nsites are greater than 2.
States_DOS = 3         # Number of energy states considered in plotting.
energy = 15            # Energy range for plotting.
contrast = 0.9         # Contrast for plotting should be between 0 and 1. The ideal value is between 0.7 to 0.9
coeff_Nc = 3           # coupling of nickelocene spin to the metallic tip apex.
coeff_M =  1           # coupling of surface spin to the underlying substrate.


```

Dependencies
1. Python 3.5 or later
2. Matplotlib and Numpy


