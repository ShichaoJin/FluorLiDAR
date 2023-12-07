# FluorLiDAR

## Step 1. PAR simulation from LiDAR points
Here, we provided an exe file for simulating PAR from LiDAR data. The detailed operation process is provided in a readme.txt in the _PAR_exe_share_ directory.
The source code of this step can be acquired by contacting Prof. Guang Zheng as declaimed in previous publications (Ma et al., 2021).  
_Ma, L., Zheng, G., Ying, Q., Hancock, S., Ju, W., & Yu, D. (2021). Characterizing the three-dimensional spatiotemporal variation of forest photosynthetically active radiation using terrestrial laser scanning data. Agricultural and Forest Meteorology, 301-302, 108346_  


## Setp2. Point-level SIF simulation from PAR simlutaion
This step includes two parts, one is the "_(1) SIF simulation at the leaf scale_" and the other is "_(2) Multiple scattering SIF simulation_".
The main source code is _SIF_daytime_Baima_distribution.m_.


## Step3. Directional SIF simulation
The main source code is _SIF_daytime_canopy_biDirectional_cylinder.m_.


Other scripts are predefined variables or key functions used by the main source codes.  
_RZT_ms2.txt_ and _RT_ms3.txt_ are generated by the _4-scale model_. 


If there are any confusion, please contact jschaon@njau.edu.cn
