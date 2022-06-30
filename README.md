# Compilation of Subduction Parameters on Global Trenches
Code to compile and analyse subduction parameters along the trenches in the GPlates models published by Müller et al. (2016) and Müller et al. (2019). Written and tested by Lior Suchoy, Imperial College London.

Please refer to Chapter 2 in the PhD thesis of Lior Suchoy at Imperial College London ***(change this with DOI after publication)*** for further details and proper referencing.

## Overview
The codes contained in this repository are used to create a compliation of parameters of subduction zones along the global subduction system back through time. It uses, as inputs, two plate reconstruction models by *Müller et al. (2016)*[^1] and *Müller et al. (2019)*[^2]. For each GPlates model it extracts age, velocities and IDs at the latitudes and lontitudes of trenches globaly. these parameters are extracted for the last 120 Myr (in 1 Myr steps). These codes also create a second compilation of subduction zone parameters along the main trenches in the Pacific and Indian Oceans during the last 60 Myr (in 10 Myr steps). Mean values and standard deviations for the age and velocities are calculated along these main trenches. This data is saved in .csv format for use in different applications.

The parameters used in the compilation for all trenches are:
1. Reconstruction time \[Myr\]
2. Location (latitude and longitude) - Latitude ranges from ±90°, with 0° at the equator. Longitude is set to range from ±180°, with 0° at the prime meridian.
3. Age of the subducting plate \[Ma\]
4. Trench azimuth \[°\] - Trench azimuth is between 0° and 360°, with 0° set to North.
5. Convergence velocity
6. Trench velocity
7. Subducting plate velocity
8. GPlates sub. plate ID
9. GPlates sub. zone ID

For each type of velocity, the compilation includes the angle to the trench \[°\] with 0° indicating trench-normal motion and ±90° indicating trench-parallel motion, magnitude \[cm/yr\], sign \[Boolean\], trench-normal velocity magnitude \[cm/yr\] and trench-parallel velocity magnitude \[cm/yr\]. Trench and subducting-plate velocities are in the absolute mantle reference frame of the reconstruction. For trench and subducting plate velocities, negative velocity is in the direction of the upper plate. For convergence velocity, negative velocity indicates divergence

Additional parameters added for the compilation of main trenches are:
1. Main plate ID
2. Main trench ID
3. Subducting plate surface area \[km$^2$\]
4. Age and velocity mean values
5. Age and velocity standard deviation

Additional compilation of parameters along spreading ridges is also generated using the attached codes. This compilation contains data for all spreading ridges during the last 120 Myr (in 1 Myr steps). It includes the following parameters: 
1. Reconstruction time \[Myr\]
2. Location (latitude and longitude) - Latitude ranges from ±90°, with 0° at the equator. Longitude is set to range from ±180°, with 0° at the prime meridian.
3. Spreading velocity \[cm/yr\]
4. Spreading obliquity \[°\]
5. Arc (segment) length \[m\]
6. Spreading arc normal azimuth \[°\]

[^1]: https://doi.org/10.1146/annurev-earth-060115-012211, GPlates model files can be found in https://www.earthbyte.org/

[^2]: https://doi.org/10.1029/2018TC005462, GPlates model files can be found in https://www.earthbyte.org/

## Files

### Python scripts
1. TrSubPropCompilation.py - Module which includes all required code to compile the data compilation of 'all trenches' and the 'main trenches'. Also includes the code to compile the spreading ridges compilation. Can also be used as script to generate all the compilations by running > python TrSubPropCompilation.py.
2. PlotsOfCompilation.py - Module with functions to generate examples of plots based on the comilations. Can be used as script to generate all plots by running > python PlotsOfCompilation.py.
3. subduction_convergence.py - Function from the .ptt module (https://github.com/EarthByte/PlateTectonicTools), used in TrSubPropCompilation.py.
4. ridge_spreading_rate.py - Function from the .ptt module (https://github.com/EarthByte/PlateTectonicTools), used in TrSubPropCompilation.py.
5. separate_ridge_transform_segments.py - Function from the .ptt module (https://github.com/EarthByte/PlateTectonicTools), used in ridge_spreading_rate.py.

### Compilation files
1. 201X/SubductionData_AllTr.csv - Compilation of the data for all trenches. 1st line includes headers.
2. 201X/SubductionData_MainTr.csv - Compilation of the data for the main trenches. 1st line includes headers.
3. 201X/RidgeData.csv - Compilation of the data for the spreading ridges. 1st line includes headers.
4. 201X/MainTrenches.csv - Table which associates reconstruction time (1st column), sub. plate ID (3rd column) and subduction zone ID (4th column) with main trench (2nd column). 1st-3rd lines include headers.
5. 201X/MainPlates.csv - Table which associates reconstruction time (1st column) and sub. plate ID (3rd column)with main plate (2nd column). 1st-4th lines include headers.
6. 201X/PointIDs.csv - Table of all combinations between plate ID (2nd column) and subduction zone ID (3rd column) through reconstruction times (1st column). 1st line includes headers.
7. 201X/IDsNamesLegend.csv - Table of all subduction zone names (as they appear in the GPlates data, 1st column) and associated subduction zone IDs (2nd column). 1st line includes headers.

### Folders
1. 201X/Muller201X_GPlatesModel - All datafiles of the GPlates models required for the codes.
2. 201X/Plots - Folder where plots are saved.

## Usage
The .csv data type is very common and can be used with any table reader (e.g. MS Excel). In the attached codes, the file is imported using Numpy loadtxt command (with delimiter=',' and skiprows=#of_header_rows). Once imported, the compilation can be manipulated as any other numpy array. More complex manipulation of the data (such as importing the headers) can be done using Pandas.

The Python scripts were built and executed using the following modules and versions:

- **Python 3.8.13**, **Cartopy 0.18.0**, **Matplotlib 3.5.1**, **Netcdf4 1.5.7** and **Numpy 1.22.3**.
- Some files from the **ptt** (plate tectonics tools) module are also included in this depository (subduction_convergence.py, ridge_spreading_rate.py and separate_ridge_transform_segments.py). The ptt module is a custom-made module for plate tectonics calculations of GPlates plate reconstruction models that was created and maintained by John Cannon and can be found at https://github.com/EarthByte/PlateTectonicTools. 
- **pyGPlates 2.8** module was used for these codes and is included in the repository.

## Work process

The compilation is generated using the code with the following workflow. Note that stage 7 is the only stage that is done manually by the user:

1. Load a GPlates plate reconstruction model (.gpml and .rot files and netcdf4 age grid files).
2. For each time step (i.e. 1 Myr), locate all the trench segments globally and sample them at intervals of 50 km to extract their location. Sample the corresponding raster age file for the nearest value to each location and assign subducting-plate age.
3. Extract the subducting plate and subduction zone IDs at each trench point. 
4. Calulate plates velocities.
5. Calculate trench obliquity and trench velocity.
6. Save the parameters as compilation of all trenches in .csv format. This file does not include the association to main trenches and main plates.
7. The user manually assigns trench segments (subduction zone ID and plate ID) to each of the main trenches (in the file MainTrenches.csv) and assigns plates (plate ID) to each of the main subducting plates (in the file MainPlates.csv).
8. Load the compilation of all trenches and the supplementary files, and add mian plate ID and main trench ID accordingly.
9. Calculate for the main trenches the mean values and standard deviation of the age and velocities and add these values to the compilation.
10. Save the extended set of parameters as the main subduction parameter compilation in .csv format.


