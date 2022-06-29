
"""
    Copyright (C) 2022 Lior Suchoy, Imperial College London, UK
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

def GenerateCompilationAllPoints(rotation_model, topology_filename, raster_filename_base,
               raster_filename_ext, WorkDir, YEAR=2016, min_time=0, max_time=120, time_delta=1, Sampling_Distance=50):
    import subduction_convergence, ridge_spreading_rate # edited from the ptt module
    import pygplates # included in the GPlates software
    import numpy as np

    print('  Extracting data from the GPlate model')
    ### Intiate arrays and variables
    tessellation_threshold_radians = (Sampling_Distance/pygplates.Earth.mean_radius_in_kms)
    SubductionData_AllTr, RidgeData = [],[]
    reconstruction_times = range(min_time,max_time+1,time_delta)
    ### Iterate over time
    for reconstruction_time in reconstruction_times:
        print('reconstruction_time: {0}\n************************'.format(reconstruction_time))
        # Determine raster filename
        raster_filename = '{0}-{1}.{2}'.format(raster_filename_base, reconstruction_time, raster_filename_ext)
        ####### Extract subduction parametrs #######
        # Use subduction convergence script to generate sample points along subduction zones at 'time'
        subduction_convergence_data = subduction_convergence.subduction_convergence(
                                      rotation_model, topology_filename, tessellation_threshold_radians, 
                                      reconstruction_time, output_distance_to_nearest_edge_of_trench = True,
                                      output_distance_to_start_edge_of_trench = True,
                                      output_convergence_velocity_components = True,
                                      output_trench_absolute_velocity_components = True,
                                      output_subducting_absolute_velocity = True,
                                      output_subducting_absolute_velocity_components = True)
        # Sample raster/grid at subduction points
        subduction_lons = [data[0] for data in subduction_convergence_data] # all lons
        subduction_lats = [data[1] for data in subduction_convergence_data] # all lats
        raster_values = sample_grid_using_scipy(subduction_lons, subduction_lats, raster_filename,YEAR) # get all ages
        # iterate over subduction points
        for point_index, age in enumerate(raster_values): 
            ### uncomment to exclude points with NaN age. Note that some inland subduction zones (e.g. Alps) may be removed
            # if np.isnan(age): # skip points with NaN age
            #     continue
            subduction_convergence_item = subduction_convergence_data[point_index]
            lon                                               = subduction_convergence_item[0]
            lat                                               = subduction_convergence_item[1]
            trench_normal_azimuth                             = subduction_convergence_item[7]
            convergence_obliquity_degrees                     = subduction_convergence_item[3]
            convergence_velocity_magnitude_cm_per_yr          = subduction_convergence_item[2]
            convergence_velocity_magnitude_sign               = np.sign(subduction_convergence_item[2])
            convergence_velocity_normal_cm_per_year           = subduction_convergence_item[12]
            convergence_velocity_parallel_cm_per_year         = subduction_convergence_item[13]
            trench_absolute_obliquity_degrees                 = subduction_convergence_item[5]
            trench_absolute_velocity_magnitude                = subduction_convergence_item[4]
            trench_absolute_velocity_magnitude_sign           = np.sign(subduction_convergence_item[4])
            trench_absolute_velocity_normal_cm_per_year       = subduction_convergence_item[14]
            trench_absolute_velocity_parallel_cm_per_year     = subduction_convergence_item[15]
            subducting_absolute_velocity_degrees              = subduction_convergence_item[17]
            subducting_absolute_velocity_magnitude_cm_per_yr  = subduction_convergence_item[16]      
            subducting_absolute_velocity_magnitude_sign       = np.sign(subduction_convergence_item[16])
            subducting_absolute_velocity_normal_cm_per_year   = subduction_convergence_item[18]
            subducting_absolute_velocity_parallel_cm_per_year = subduction_convergence_item[19]
            subducting_plate_id                               = subduction_convergence_item[8]
            subduction_zone_plate_id                          = subduction_convergence_item[9]
            ### subducting_length_km, distance_to_start_edge_of_trench_km and distance_to_nearest_edge_of_trench_km are not 
            ### working properly at the moment
            # subducting_length_km                              = (subduction_convergence_item[6] * 
            #                                                      pygplates.Earth.mean_radius_in_kms)
            # distance_to_start_edge_of_trench_km               = (subduction_convergence_item[11] * 
            #                                                      pygplates.Earth.mean_radius_in_kms)
            # distance_to_nearest_edge_of_trench_km             = (subduction_convergence_item[10] *
            #                                                      pygplates.Earth.mean_radius_in_kms)
            SubductionData_AllTr.append([
                reconstruction_time,
                lon,
                lat,
                age,
                trench_normal_azimuth,
                convergence_obliquity_degrees,
                convergence_velocity_magnitude_cm_per_yr,
                convergence_velocity_magnitude_sign,
                convergence_velocity_normal_cm_per_year,
                convergence_velocity_parallel_cm_per_year,
                trench_absolute_obliquity_degrees,
                trench_absolute_velocity_magnitude,
                trench_absolute_velocity_magnitude_sign,
                trench_absolute_velocity_normal_cm_per_year,
                trench_absolute_velocity_parallel_cm_per_year,
                subducting_absolute_velocity_degrees,
                subducting_absolute_velocity_magnitude_cm_per_yr,        
                subducting_absolute_velocity_magnitude_sign,
                subducting_absolute_velocity_normal_cm_per_year,
                subducting_absolute_velocity_parallel_cm_per_year,
                subducting_plate_id,
                subduction_zone_plate_id])#,
                ### subducting_length_km, distance_to_start_edge_of_trench_km and distance_to_nearest_edge_of_trench_km are not 
                ### working properly at the moment
                # subducting_length_km,
                # distance_to_start_edge_of_trench_km,
                # distance_to_nearest_edge_of_trench_km])
        
        ####### Extract spreading ridges parametrs #######
        # Use ridge spreading script to generate sample points along subduction zones at 'time'
        ridge_spreading_data = ridge_spreading_rate.spreading_rates_dense(
                               rotation_model, topology_filename, reconstruction_time, tessellation_threshold_radians)
        # Sample raster/grid at subduction points
        ridge_lons = [data[0] for data in ridge_spreading_data] # all lons
        ridge_lats = [data[1] for data in ridge_spreading_data] # all lats
        raster_values = sample_grid_using_scipy(ridge_lons, ridge_lats, raster_filename,YEAR) # get all ages
        # iterate over ridge points
        for point_index, age in enumerate(raster_values): 
            if np.isnan(age): # check for NaN age
                age = 0 # because this is a spreading ridge...
                ### uncomment to exclude points with NaN age. May remove some inland rifts (e.g. East Africa)
                # continue 
            ridge_spreading_item = ridge_spreading_data[point_index]
            lon                                    = ridge_spreading_item[0]
            lat                                    = ridge_spreading_item[1]
            spreading_velocity_magnitude_cm_per_yr = ridge_spreading_item[2]
            spreading_obliquity_degrees            = ridge_spreading_item[3]
            arc_length_metres                      = (ridge_spreading_item[4]*1e3*
                                                      pygplates.Earth.mean_radius_in_kms)
            spreading_arc_normal_azimuth           = ridge_spreading_item[5]
            RidgeData.append([reconstruction_time,
                               lon,
                               lat,
                               age,
                               spreading_velocity_magnitude_cm_per_yr,
                               spreading_obliquity_degrees,
                               arc_length_metres,
                               spreading_arc_normal_azimuth])
    
    print('  Saving compilation data')
    ### Save data for all trenches globally to SubductionData_AllTr.csv file
    print('    Saving subduction data')
    Titles = ["reconstruction_time",
              "lon",
              "lat",
              "age",
              "trench_normal_azimuth",
              "convergence_obliquity_degrees",
              "convergence_velocity_magnitude_cm_per_yr",
              "convergence_velocity_magnitude_sign",
              "convergence_velocity_normal_cm_per_year",
              "convergence_velocity_parallel_cm_per_year",
              "trench_absolute_obliquity_degrees",
              "trench_absolute_velocity_magnitude",
              "trench_absolute_velocity_magnitude_sign",
              "trench_absolute_velocity_normal_cm_per_year",
              "trench_absolute_velocity_parallel_cm_per_year",
              "subducting_absolute_velocity_degrees",
              "subducting_absolute_velocity_magnitude_cm_per_yr",        
              "subducting_absolute_velocity_magnitude_sign",
              "subducting_absolute_velocity_normal_cm_per_year",
              "subducting_absolute_velocity_parallel_cm_per_year",
              "subducting_plate_id",
              "subduction_zone_plate_id"]#,
              ### subducting_length_km, distance_to_start_edge_of_trench_km and distance_to_nearest_edge_of_trench_km are not 
              ### working properly at the moment
              #"subducting_length_km",
              #"distance_to_start_edge_of_trench_km",
              #"distance_to_nearest_edge_of_trench_km"]
    Save_data_as_CSV(Titles,SubductionData_AllTr,WorkDir,'SubductionData_AllTr')
    
    ### Saves the ridge data to RidgeData.csv file
    print('    Saving ridge data')
    Titles= ["reconstruction_time",
             "lon",
             "lat",
             "age",
             "spreading_velocity_magnitude_cm_per_yr",
             "spreading_obliquity_degrees",
             "arc_length_metres",
             "spreading_arc_normal_azimuth"]
    Save_data_as_CSV(Titles,RidgeData,WorkDir,'RidgeData')

def SaveIDsAndDescriptionsAllPoints(topology_filename,WorkDir):    
    import pygplates, os
    import numpy as np
    print('  Saving IDs legend and PointID data')
    ### Saves list of registered names of subduction zone IDs to IDsNamesLegend.csv file
    ### This is used to identify IDs of specific segments for the main trenches
    print('    Saving IDs legend')
    Features = []
    domain_features = pygplates.FeatureCollection(topology_filename)
    for F in domain_features:
        ### uncomment to get only the IDs of subduction zones. Note this sometimes miss important IDs
        #if F.get_feature_type() == pygplates.FeatureType.gpml_subduction_zone:
            Features.append([F.get_name(),F.get_reconstruction_plate_id()])
    Titles = ["Name","SubZoneID"]
    Save_data_as_CSV(Titles,Features,WorkDir,'IDsNamesLegend')
    
    ### Saves list of subduction zone IDs and subductiong plate IDs at each reconstruction 
    ### time to PointIDs.csv file
    print('    Saving PointID data')
    SubductionData_AllTr = np.loadtxt(os.path.join(WorkDir,'SubductionData_AllTr.csv'),delimiter=',', skiprows=1) # Load 'all trecnhes' data
    Titles = ["reconstruction_time",
              "subducting_plate_id",
              "subduction_zone_plate_id"]
    Save_data_as_CSV(Titles,np.unique(SubductionData_AllTr[:,[0,20,21]], axis=0).tolist(),WorkDir,'PointIDs')

def GenerateCompilationMainTrenches(WorkDir):
    import numpy as np
    import os
    print('  Creating compilation of main trenches and plates')
    print('    Compiling main trenches')
    SubductionData_AllTr = np.loadtxt(os.path.join(WorkDir,'SubductionData_AllTr.csv'),delimiter=',', skiprows=1) # Load 'all trecnhes' data
    temp = np.hstack((SubductionData_AllTr, np.zeros((SubductionData_AllTr.shape[0],2)), np.empty((SubductionData_AllTr.shape[0],26))*np.nan)) # add two columns of zeroes for the main trench and main plate IDs and 26 columns of NaNs for the mean and std values
    ### The MainTrenches.csv file is manually created by the user and includes which trench segments at which times be merged into any main trench
    Trenches = np.loadtxt(os.path.join(WorkDir,'MainTrenches.csv'), delimiter = ',', skiprows = 3)
    for i in np.arange(Trenches.shape[0]):
        # find all the poiunts with specific rec. time and IDs
        Trench = temp[(temp[:,0] == Trenches[i,0]) &
                      (temp[:,20] == Trenches[i,2]) &
                      (temp[:,21] == Trenches[i,3]),:]
        if Trench.any(): # check if there are points for the specific rec. time and IDs
            # Check if segments appear in more than one main trench
            temp_app = temp[(temp[:,0] == Trenches[i,0]) &
                            (temp[:,20] == Trenches[i,2]) &
                            (temp[:,21] == Trenches[i,3]) &
                            (temp[:,22] > 0.) ,:]
            if (temp_app.any()):
                temp_app[:,23] = Trenches[i,1]
                # add new trench line with the additional ID
                if len(temp_app[:,23])>1.:
                    temp = np.append(temp,temp_app,axis=0)
                else:
                    temp = np.append(temp,[temp_app],axis=0)
            else:
                ## Set trench value
                temp[(temp[:,0] == Trenches[i,0]) &
                     (temp[:,20] == Trenches[i,2]) &
                     (temp[:,21] == Trenches[i,3]) ,22] = Trenches[i,1]
    # SubductionData_MainTr = np.array(temp[temp[:,22]>0,:]) ### Uncomment to have the 'main trenches' compilation exclude all non-main trenches
    SubductionData_MainTr = np.array(temp) ### Add comment to have the 'main trenches' compilation exclude all non-main trenches

    print('    Compiling main plates')
    ### The MainPlates.csv file is manually created by the user and includes which plate IDs at which times be merged into any main plate
    Plates = np.loadtxt(os.path.join(WorkDir,'MainPlates.csv'), delimiter = ',', skiprows = 4)
    for i in np.arange(Plates.shape[0]):
        # find all the poiunts with specific rec. time and plate IDs
        Plate = SubductionData_MainTr[(SubductionData_MainTr[:,0] == Plates[i,0]) &
                                      (SubductionData_MainTr[:,20] == Plates[i,2]),:]
        if Plate.any(): # check if there are points for the specific rec. time and plate ID
            ## Set main plate  value
            SubductionData_MainTr[(SubductionData_MainTr[:,0] == Plates[i,0]) &
                                  (SubductionData_MainTr[:,20] == Plates[i,2]) ,23] = Plates[i,1]

    print('    Calculating mean and standard deviation values for main trenches')
    for time in np.unique(SubductionData_MainTr[:,0]):
        temp = SubductionData_MainTr[(SubductionData_MainTr[:,0] == time),:]
        for tr in np.unique(temp[:,22]):
            if tr>0:
                temp_tr = temp[temp[:,22] == tr, :]
                Ages_mean, Ages_std = np.average(temp_tr[(np.invert(np.isnan(temp_tr[:,3]))),3]), np.std(temp_tr[(np.invert(np.isnan(temp_tr[:,3]))),3])
                Vels_mean, Vels_std = np.average(temp_tr[:,[5,6,8,9,10,11,13,14,15,16,18,19]],axis=0), np.std(temp_tr[:,[5,6,8,9,10,11,13,14,15,16,18,19]],axis=0)
                SubductionData_MainTr[(SubductionData_MainTr[:,0] == time) & (SubductionData_MainTr[:,22] == tr),24] = Ages_mean
                SubductionData_MainTr[(SubductionData_MainTr[:,0] == time) & (SubductionData_MainTr[:,22] == tr),25:37] = Vels_mean
                SubductionData_MainTr[(SubductionData_MainTr[:,0] == time) & (SubductionData_MainTr[:,22] == tr),37] = Ages_std
                SubductionData_MainTr[(SubductionData_MainTr[:,0] == time) & (SubductionData_MainTr[:,22] == tr),38:50]  = Vels_std
    
    SubductionData_MainTr = SubductionData_MainTr.tolist()
    print('  Saving compilation of main trenches')
    ### Save data for the main trenches to SubductionData_MainTr.csv file
    Titles = ["time",
              "lon",
              "lat",
              "age",
              "trench_normal_azimuth",
              "convergence_obliquity_degrees",
              "convergence_velocity_magnitude_cm_per_yr",
              "convergence_velocity_magnitude_sign",
              "convergence_velocity_normal_cm_per_year",
              "convergence_velocity_parallel_cm_per_year",
              "trench_absolute_obliquity_degrees",
              "trench_absolute_velocity_magnitude",
              "trench_absolute_velocity_magnitude_sign",
              "trench_absolute_velocity_normal_cm_per_year",
              "trench_absolute_velocity_parallel_cm_per_year",
              "subducting_absolute_velocity_degrees",
              "subducting_absolute_velocity_magnitude_cm_per_yr",        
              "subducting_absolute_velocity_magnitude_sign",
              "subducting_absolute_velocity_normal_cm_per_year",
              "subducting_absolute_velocity_parallel_cm_per_year",
              "subducting_plate_id",
              "subduction_zone_plate_id",
              ### subducting_length_km, distance_to_start_edge_of_trench_km and distance_to_nearest_edge_of_trench_km are not 
              ### working properly at the moment
              # "distance_to_start_edge_of_trench_km",
              # "distance_to_nearest_edge_of_trench_km",
              # "subducting_length_km",
              "main_trench_ID",
              "mian_plate_ID",
              "mean_trench_age",
              "mean_convergence_azimuth",
              "mean_convergence_abs_vel",
              "mean_convergence_nor_vel",
              "mean_convergence_par_vel",
              "mean_trench_azimuth",
              "mean_trench_abs_vel",
              "mean_trench_nor_vel",
              "mean_trench_par_vel",
              "mean_subduction_azimuth",
              "mean_subduction_abs_vel",
              "mean_subduction_nor_vel",
              "mean_subduction_par_vel",
              "std_dev_trench_age",
              "std_dev_convergence_azimuth",
              "std_dev_convergence_abs_vel",
              "std_dev_convergence_nor_vel",
              "std_dev_convergence_par_vel",
              "std_dev_trench_azimuth",
              "std_dev_trench_abs_vel",
              "std_dev_trench_nor_vel",
              "std_dev_trench_par_vel",
              "std_dev_subduction_azimuth",
              "std_dev_subduction_abs_vel",
              "std_dev_subduction_nor_vel",
              "std_dev_subduction_par_vel"]#,
              # "plate_surface_area",
              # "calculated_trench_length"]
    Save_data_as_CSV(Titles,SubductionData_MainTr,WorkDir,'SubductionData_MainTr')
            
def sample_grid_using_scipy(lon_in,lat_in,agefile,YEAR):
    import numpy as np
    from netCDF4 import Dataset as netcdf
    from itertools import product
    Distance = 10 # x0.1 of degree, maximum distance to search raster (age) data
    data=netcdf(agefile,'r') # get the full data from the raster file
    Zg = data.variables['z'][:] # get the age data from the 'data' variable
    result = []
    for xi,yi in zip(lon_in,lat_in):
        # Change Lon/Lat format to start from 0 and not from -180/-90. Also round to the nearest 10th of a degree 
        if YEAR==2016:
            Lat, Lon = int(10*(yi+90)),int(10*xi)#int(np.round(10.*(yi+90))), int(np.round(10.*(xi+180)))
        elif YEAR==2019:
            Lat, Lon = int(np.round(10.*(yi+90))), int(np.round(10.*(xi+180)))
            Lat, Lon = min(max(Lat,0),1800), min(max(Lon,0),3600) # Limit the range (0-1800/0-3600)
        Age = float(Zg[Lat,Lon]) # get age at point
        if np.isnan(Age): # Check if age is NaN
            for R in range(1,Distance+1): # If NaN, check neighbours at increased distance up to 'Distance'
                temp = np.array([i for i in product((range(-R,R+1)), (range(-R,R+1)))]) # coordinates combinations
                AdjacentCells = np.unique(np.append(temp[np.abs(temp[:,0])==R], temp[np.abs(temp[:,1])==R], axis=0), axis=0) # unique combinations
                temp_ages = []
                for i in AdjacentCells: # check all neighbouring cells
                    Lat_Neigh = Lat+i[0]
                    Lon_Neigh = Lon+i[1]
                    Lat_Neigh = min(max(Lat+i[0],0),1800)
                    if YEAR==2016:
                        if Lon_Neigh<-1799: Lon_Neigh=Lon_Neigh+3600
                        if Lon_Neigh>1799: Lon_Neigh=Lon_Neigh-3600
                    elif YEAR==2019:
                        if Lon_Neigh<1: Lon_Neigh=Lon_Neigh+3600
                        if Lon_Neigh>3599: Lon_Neigh=Lon_Neigh-3600
                    temp_ages.append(float(Zg[Lat_Neigh,Lon_Neigh]))
                temp_ages = np.array(temp_ages)
                L = np.count_nonzero(~np.isnan(temp_ages))
                # Adjust age if any neighbour at distance R is non-NaN, increase distance if all NaN
                if L>0:
                    Age = np.nansum(temp_ages)/L # average non-nan values
                    break
                else:
                    Age = np.nan
        result.append(Age)
    return result
    
def Save_data_as_CSV(Titles,Data,Dir,Name):
    import os
    CleanData = str(Titles)+'\n'+str(Data)
    CleanData = CleanData.replace("], [", "\n")
    CleanData = CleanData.replace("]", "")
    CleanData = CleanData.replace("[", "")
    CleanData = CleanData.replace("(", "")
    CleanData = CleanData.replace(")", "")
    CleanData = CleanData.replace('"', "")
    CleanData = CleanData.replace("'", "")
    with open(os.path.join(Dir,Name+'.csv'), 'w') as O_file:
        O_file.write(CleanData)

def Merge_rot_files(Data_Folder, ROTFile):
    import os
    print("Merging .rot files:")
    temp = ""
    for name in os.listdir(Data_Folder): #Find all .rot files in Data_Folder
        if (name.endswith(".rot")): #Skip non .rot files 
            print('  adding {0} to {1} file '.format(name, ROTFile))
            fp = open(os.path.join(Data_Folder, name), encoding="utf8") #Encoding is important!
            temp += fp.read()
    with open (os.path.join(Data_Folder,ROTFile), 'w', encoding="utf8") as fp:
        Data = "".join([s for s in temp.strip().splitlines(True) if s.strip()]) #Remove empty lines
        fp.write(Data) #Write combined file

    
def Merge_gpml_files(Data_Folder, GPMLFile):
    import os, pygplates
    all_features = []
    for name in os.listdir(Data_Folder): #Find all .gpml files in Data_Folder
        if (name.endswith(".gpml")): #Skip non .gpml files 
            print('  adding {0} to {1} file '.format(name, GPMLFile))
            file = os.path.join(Data_Folder,name)
            features_in_file = pygplates.FeatureCollection(file)
            all_features.extend(features_in_file)
    pygplates.FeatureCollection(all_features).write(os.path.join(Data_Folder,GPMLFile))
    
if __name__ == '__main__':
    import gc, os, pygplates
    gc.collect()
    def main():
        ### Compile the 2019 compilation
        YEAR = 2019
        WorkDir = os.path.join(os.getcwd(),'2019_Compilation')
        Data_Folder = os.path.join(WorkDir,'Muller2019_GPlatesModel')
        raster_folder = os.path.join(Data_Folder,'Age_Raster_Files')
        raster_filename_base = os.path.join(raster_folder,'Muller_etal_2019_Tectonics_v2.0_AgeGrid')
        print('\nCreating 2019 compilation')
        set_variables(WorkDir,Data_Folder,raster_folder,raster_filename_base,YEAR)
        ### Compile the 2016 compilation
        YEAR = 2016
        WorkDir = os.path.join(os.getcwd(),'2016_Compilation')
        Data_Folder = os.path.join(WorkDir,'Muller2016_GPlatesModel')
        raster_folder = os.path.join(Data_Folder,'Age_Raster_Files')
        raster_filename_base = os.path.join(raster_folder,'EarthByte_AREPS_v1.15_Muller_etal_2016_AgeGrid')
        print('\nCreating 2016 compilation')
        set_variables(WorkDir,Data_Folder,raster_folder,raster_filename_base,YEAR)

    def set_variables(WorkDir,Data_Folder,raster_folder,raster_filename_base,YEAR):
        raster_filename_ext  = 'nc'
        ROTFile  = "MergedRot.rot"
        if ROTFile not in os.listdir(Data_Folder): #check if .rot files are merged
            Merge_rot_files(Data_Folder, ROTFile) #merge .rot files
        rotation_model = pygplates.RotationModel(os.path.join(Data_Folder, ROTFile))
        GPMLFile  = "MergedGPML.gpml"
        if GPMLFile not in os.listdir(Data_Folder): #check if .gpml files are merged
            Merge_gpml_files(Data_Folder, GPMLFile) #merge .gpml files
        topology_filename    = os.path.join(Data_Folder, GPMLFile)
        
        GenerateCompilationAllPoints(rotation_model, topology_filename, raster_filename_base,
                                     raster_filename_ext, WorkDir, YEAR, 0, 120, 1, 50)
        SaveIDsAndDescriptionsAllPoints(topology_filename, WorkDir)
        GenerateCompilationMainTrenches(WorkDir)

    main()