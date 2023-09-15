# Nickelocene as Spin Sensor

Download ```Spin_Sensing_Nc.py ``` and ``` parameters.txt``` these two files and run the code in the same directory to use the code.

``` python3 Spin_sensing_Nc.py parameters.txt ```

You will get an image.

In the parameter file, you can edit the input parameters. 
All the input parameters are: 

```
M_spin = 1/2             # Spin of Metal.
D_aniso = [0,0,5]      # Magneitc anisotropy of Metal [D_x, D_y, D_z].
Temp = 0.2             # Temperature.
B_field = 0            # Magnetic Field.
J_Nc =  3               # Magnetic Exchange between metal and Nc.
Nsites = 2             # Number of spin sites.
contrast = 0.7         # Contrast for plotting should be between 0 and 1. The ideal value is between 0.7 to 0.9.
### Only change the above parameters if you don't know what parameters are given below.
J1  = [6]              # Nearest Neighbour exchange when Nsites are greater than 1. Can be a list like [1,2,3]
J2  = 0                # second Neighbour exchange when Nsites are greater than 2.
States_DOS = 3         # Number of energy states considered in plotting.
energy = 15            # Energy range for plotting.
coeff_Nc = 3           # coupling of nickelocene spin to the metallic tip apex.
coeff_M =  1           # coupling of surface spin to the underlying substrate.
ciclo = 0              # For many spin sites, if you want to connect the last site to the first site, put 1.
Nc_M_connect = 1       # Which Nsite you want to scan with Nc.


```

Dependencies
1. Python 3.5 or later
2. Matplotlib and Numpy


