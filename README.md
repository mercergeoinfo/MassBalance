# Mass Balance Calculator

This programme goes step by step through the calculation of surface mass balance using accumulation and ablation data
The source files are set in the ''settings.txt'' file.  
The required source files are:

- density by depth
- snow depth probings (with coordinates)
- ablation stake readings (with coordinates)
- a DEM of the glacier
- a mask of the DEM



Run the massbalance.py script and follow the sequential instructions.
Note that the data is split into elevation bands after kriging/extrapolation.
The files used for input should be prepared so that the density profile is either a Toikka snow fork file

```
DEPTH,ATTENUATION,FREQUENCY,BAND,PERMITTIVITY,PERMITTIVITY,WETNESS,DENSITY,WETNESS,Sample Length,Cumulative Snow,Sample Mass,Cumulative Mass
5,810,708.9,13.5,1.57,0,0,0.296,0,5,0.05,0.015,0.015
10,837,664,15,1.79,0.005,0.68,0.376,1.7,5,0.10,0.019,0.034
15,861,659,14,1.82,0.003,0.41,0.398,1,5,0.15,0.020,0.054
```
or a density some other cumulative density data file.

The snow probing data should have snow depth plus coordinates for each point, e.g.

```
P,E,N,D
02C,651325,7536365,3.83
02N1,651331,7536463,4.68
03C,651225,7536370,2.77
03N1,651231,7536468,2.35
```
In this above example the depth is given in metres.  

The ablation stake data should give coordinates, stake height and snow depth in winter and stake height in summer.

```
Stake,Easting,Northing,Elevation,Hw,Dw,Hs,Ds
22S2,649311.28,7536285.52,1399.77,1.56,2.48,6.27,-2.23
20S5,649484.12,7535981.59,1385.87,0.9,2.19,6.17,-3.08
20S4,649505.86,7536084.37,1379.31,2.36,1.78,6.49,-2.35
```
