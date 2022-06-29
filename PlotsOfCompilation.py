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

def PlotTrenchThroughTime(SubductionData_MainTr,SaveDir):
    import os
    import numpy as np
    import cartopy as cpy
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    print('\n  Creating plots of main trenches locations through time')
    # extracting all trenches and reconstruction times
    ListReconstructionTimes = np.unique(SubductionData_MainTr[(SubductionData_MainTr[:,22]>0),0])
    ListTrenches = np.unique(SubductionData_MainTr[:,22])
    VMin = ListReconstructionTimes.min()
    VMax = ListReconstructionTimes.max()+10

    for Tr in ListTrenches: # looping through trenches
        if Tr==0: # the 0 trench plots all main trenches globally
            print('    All trenches')
            temp = SubductionData_MainTr[(SubductionData_MainTr[:,22]>0),:]
        else:
            print('    Trench no.: {0}'.format(Tr))
            temp = SubductionData_MainTr[(SubductionData_MainTr[:,22]==Tr),:]
        (ReconstructionTimes,Lons,Lats) = (temp[:,0],temp[:,1],temp[:,2]) # extrcating the parameters for specific trench
        Projection = ccrs.Orthographic(central_longitude=(0.5*(min(Lons)+max(Lons)))) #Set projection
        if Tr==0: # Correction for projection of the Global plot
            Projection = ccrs.Robinson(central_longitude=180.+(0.5*(min(Lons)+max(Lons))))
        if Tr in [3, 7, 9]: # Correction for projection of the Alleutian, Tonga and South Pacific plots
            Projection = ccrs.Orthographic(central_longitude=180.+(0.5*(min(Lons)+max(Lons))))
        Cmap = plt.get_cmap('jet', ListReconstructionTimes.size) # JET Color map

        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1, projection=Projection)
        ax.add_feature(cpy.feature.COASTLINE, color=[0,0,0,.1]) # COASTLINE (present day) feature from cartopy
        ax.add_feature(cpy.feature.LAND, color=[0,0,0,.1]) # LAND (landmass at present day) feature from cartopy
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color=[.9,.7,.7], linestyle='--',zorder=1) # Plot gridlines
        sctr = ax.scatter(Lons,Lats,c=ReconstructionTimes,vmin=VMin,vmax=VMax,cmap=Cmap,transform=ccrs.PlateCarree(),zorder=2) # Plot data
        cb = fig.colorbar(sctr,orientation='horizontal') # Colorbar (horizontal)
        cb.set_label('Reconstruction time [Myr]') # Colorbar label

        fig.tight_layout()
        Path = os.path.join(SaveDir,'TrenchesThroughTime') # Directory under 'plots'
        if not os.path.exists(Path): os.makedirs(Path) # check if directory exists and if not create it
        fig.savefig(os.path.join(Path,'Trench_{0}.png'.format(int(Tr))), facecolor='white', dpi=300) # Save .png
        fig.savefig(os.path.join(Path,'Trench_{0}.svg'.format(int(Tr))), facecolor='white', dpi=300) # Save .svg (for postprocessing)
        plt.close()

def GlobalAgePlots(SubductionData_AllTr, topology_filename, rotation_model, raster_filename_base, raster_filename_ext, coastline_features, YEAR, SaveDir, Vel_arrows_grid_spacing_degrees=10.0, AgeGridRes=5.):
    import pygplates
    import numpy as np
    from itertools import product
    from netCDF4 import Dataset as netcdf
    import matplotlib.pyplot as plt
    import cartopy as cpy
    import cartopy.crs as ccrs

    def Main_GlobalAgePlots(SubductionData_AllTr, topology_filename, rotation_model, raster_filename_base, raster_filename_ext, coastline_features, YEAR, SaveDir, Vel_arrows_grid_spacing_degrees, AgeGridRes):
        print('\n  Creating global plots of ages and velocities')
        Vel_arrows_loc_lats = np.array([i for i in range(int(np.floor(180./Vel_arrows_grid_spacing_degrees)))])
        Vel_arrows_loc_lons = np.array([i for i in range(int(np.floor(360./Vel_arrows_grid_spacing_degrees)))])
        Vel_arrows_locations = np.array([i for i in product((-90+(Vel_arrows_loc_lats+0.5)*Vel_arrows_grid_spacing_degrees),
                                                            (-180+(Vel_arrows_loc_lons+0.5)*Vel_arrows_grid_spacing_degrees))]) # coordinates combinations
        for time in range(0,120+1,10):
            print('    Reconstruction time: {0}'.format(time))
            temp = SubductionData_AllTr[(SubductionData_AllTr[:,0] == time),:] # extract global trenches at specific time
            if np.any(temp): # check if there are available trenches to plot (if not, skip time)
                Lons, Lats = temp[:,1],temp[:,2] # get lons and lats of trenches at reconstruction time
                ### extract plate boundaries at reconstructed time from the gpml files
                resolved_topologies, PlateBoundaries = [], [] # initialise
                topology_features = pygplates.FeatureCollection(topology_filename)
                pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time)
                points_in_resolved_topologies = {resolved_topology : [] for resolved_topology in resolved_topologies}
                for topology in resolved_topologies:
                    LatsLons = topology.get_resolved_geometry().to_lat_lon_array()
                    PlateBoundaries.append(LatsLons)
                ### calculate velocities at the velocity arrows locations at reconstructed time from the gpml files
                temp_latlon,temp_vel,Vel_arrows_points, Vel_arrows_vectors = [],[],[],[] # initialise
                for i in range(Vel_arrows_locations.shape[0]): # loop through vel arrows mesh
                    point = pygplates.PointOnSphere(Vel_arrows_locations[i,0],Vel_arrows_locations[i,1]) # point on mesh
                    for resolved_topology in resolved_topologies: # find to which polygon the point is associated
                        if resolved_topology.get_resolved_boundary().is_point_in_polygon(point):
                            points_in_resolved_topologies[resolved_topology].append(point)
                            break  
                for resolved_topology, points_in_resolved_topology in points_in_resolved_topologies.items(): # loop through mesh points
                    if not points_in_resolved_topology: continue # ignore points which are not associated with polygons (should not happen)
                    plate_id = resolved_topology.get_feature().get_reconstruction_plate_id()
                    plate_stage_rotation = rotation_model.get_rotation(time, plate_id, time+1,fixed_plate_id=0, anchor_plate_id=0)
                    plate_velocity_vectors = pygplates.calculate_velocities(points_in_resolved_topology,plate_stage_rotation,
                                                                            1,pygplates.VelocityUnits.cms_per_yr) # extract associated polygon velocity
                    for i in range(len(points_in_resolved_topology)):
                        temp_vel = pygplates.LocalCartesian.convert_from_geocentric_to_magnitude_azimuth_inclination(points_in_resolved_topology[i],
                                                                                                                     plate_velocity_vectors[i])
                        if temp_vel[0]>0.: # check that velocity magnitude is >0
                            temp_latlon = points_in_resolved_topology[i].to_lat_lon()
                            Vel_arrows_points.append([temp_latlon[0], temp_latlon[1]])
                            Vel_arrows_vectors.append([temp_vel[0]*np.cos(temp_vel[1]),temp_vel[0]*np.sin(temp_vel[1])])
                Vel_arrows_points, Vel_arrows_vectors = np.array(Vel_arrows_points), np.array(Vel_arrows_vectors)
                ### extract coastlines at reconstructed time from the gpml files
                Coastlines, resolved_coastlines = [], []
                pygplates.reconstruct(coastline_features, rotation_model, resolved_coastlines, time)
                Lim = 178
                if (YEAR==2019) and (time==30): Lim = 170
                elif (YEAR==2019) and (time==110): Lim = 130
                elif (YEAR==2019) and (time==100): Lim = 160
                elif (YEAR==2016) and (time==90): Lim = 170
                elif (YEAR==2016) and (time==110): Lim = 170
                elif (YEAR==2016) and (time==120): Lim = 170
                for feature in resolved_coastlines:
                    LatsLons = feature.get_reconstructed_geometry().to_lat_lon_array()
                    if any(LatsLons[:,1]<-Lim) and any(LatsLons[:,1]>Lim):
                        if any(LatsLons[:,1]>-10) and any(LatsLons[:,1]<10):
                            continue
                        else:
                            LatsLons[LatsLons[:,1]<0,1]=LatsLons[LatsLons[:,1]<0,1]+360
                    Coastlines.append(LatsLons)
                ### extract age grid at reconstructed time from the raster files
                AgeGridRes = 5.
                Lons_AgeGrid, Lats_AgeGrid = np.arange(0,3601,AgeGridRes), np.arange(0,1801,AgeGridRes)
                Lons_AgeGrid, Lats_AgeGrid = np.meshgrid(Lons_AgeGrid, Lats_AgeGrid)
                Lons_AgeGrid, Lats_AgeGrid = Lons_AgeGrid.flatten(), Lats_AgeGrid.flatten()
                raster_filename = '{0}-{1}.{2}'.format(raster_filename_base, time, raster_filename_ext)
                data=netcdf(raster_filename,'r')
                Zg = data.variables['z'][:]
                Age_AgeGrid = []
                for i in range(len(Lons_AgeGrid)):
                    Age_AgeGrid.append(float(Zg[int(Lats_AgeGrid[i]),int(Lons_AgeGrid[i])]))
                
                ### Plot and save
                fig = plt.figure()
                PlotMaps(Lons,Lats,PlateBoundaries,Coastlines,Lons_AgeGrid,Lats_AgeGrid,Age_AgeGrid,AgeGridRes,Vel_arrows_points,Vel_arrows_vectors,YEAR)
                fig.tight_layout()
                Path = os.path.join(SaveDir,"GlobalAgeMaps")
                if not os.path.exists(Path):
                    os.makedirs(Path)
                fig.savefig(os.path.join(Path,'Map_{0}.png'.format(int(time))), dpi=300, facecolor='white', transparent=False)
                fig.savefig(os.path.join(Path,'Map_{0}.svg'.format(int(time))), dpi=300, facecolor='white', transparent=False)
                plt.close()

    def PlotMaps(Lons,Lats,PlateBoundaries,Coastlines,Lons_AgeGrid,Lats_AgeGrid,Age_AgeGrid,AgeGridRes,Vel_arrows_points,Vel_arrows_vectors,YEAR):
        C_trenches, C_plates, C_gridlines, C_landmass, C_arrows = 'lime', 'black', '#e377c2', [.8,.8,.8,], '#d62728'
        ax = plt.subplot(111, projection=ccrs.Robinson(central_longitude=165))
        ax.set_global()
        # Plot plate boundaries
        for LatsLonsPlates in PlateBoundaries:
            B_Lons, B_Lats = LatsLonsPlates[:,1], LatsLonsPlates[:,0]
            if any(B_Lons<-150) & any(B_Lons>150):
                B_Lons_aftr = np.copy(B_Lons)
                B_Lats_aftr = np.copy(B_Lats)
                while any(B_Lons_aftr<-150) & any(B_Lons_aftr>150):
                    if any((B_Lons_aftr>-10) & (B_Lons_aftr<10)):
                        B_Lons_bfr, B_Lons_aftr, B_Lats_bfr, B_Lats_aftr = BreakLonsLats(B_Lons_aftr,B_Lats_aftr)
                        ax.plot(B_Lons_bfr, B_Lats_bfr, c=C_plates, transform=ccrs.PlateCarree(), zorder=5, rasterized=True)
                    else:
                        B_Lons_aftr[B_Lons_aftr<0] = B_Lons_aftr[B_Lons_aftr<0]+360.
                        ax.plot(B_Lons_aftr, B_Lats_aftr, c=C_plates, transform=ccrs.PlateCarree(), zorder=5, rasterized=True)
            else:
                ax.plot(B_Lons, B_Lats, c=C_plates, transform=ccrs.PlateCarree(), zorder=5, rasterized=True)
        # Plot landmass
        for LatsLonsCoastlines in Coastlines:
            ax.fill(LatsLonsCoastlines[:,1], LatsLonsCoastlines[:,0], transform=ccrs.PlateCarree(), facecolor=C_landmass, zorder=0, rasterized=True)      
        # Plot age grid
        if YEAR==2019: ax.scatter(0.1*Lons_AgeGrid-180, 0.1*Lats_AgeGrid-90, c=Age_AgeGrid, cmap='gist_earth_r', vmin=0, vmax=120, transform=ccrs.PlateCarree(), zorder=3, rasterized=True)
        elif YEAR==2016: ax.scatter(0.1*Lons_AgeGrid, 0.1*Lats_AgeGrid-90, c=Age_AgeGrid, cmap='gist_earth_r', vmin=0, vmax=120, transform=ccrs.PlateCarree(), zorder=3, rasterized=True)

        # Plot velocity arrows
        q = ax.quiver(Vel_arrows_points[:,1], Vel_arrows_points[:,0], Vel_arrows_vectors[:,1], 
                      Vel_arrows_vectors[:,0], color=C_arrows, transform=ccrs.PlateCarree(), zorder=9)
        ax.quiverkey(q, 0.1, 0.1, 5, r'$5 \frac{cm}{yr}$', labelpos='N',coordinates='figure')
        # Plot trenches
        sc = ax.scatter(Lons, Lats, c=C_trenches, transform=ccrs.PlateCarree(), zorder=7, rasterized=True)
        # Plot gridlines
        gl = ax.gridlines(crs=ccrs.PlateCarree(), color=C_gridlines, linestyle='--',zorder=3)

    def BreakLonsLats(Lons,Lats):
        Lons_Diff_Idx = np.flatnonzero(np.abs(np.diff(Lons))>300)
        if any(Lons_Diff_Idx):
            Lons_Diff_Idx = Lons_Diff_Idx[0] # index of first high gradient
            if Lons[Lons_Diff_Idx+1]>Lons[Lons_Diff_Idx]:
                Lons_MinIdx = np.copy(Lons_Diff_Idx)
                Lons_MaxIdx = Lons_Diff_Idx+1
                Add_Lon_bfr, Add_Lon_aftr = [-180], [180]
                Add_Lat_bfr, Add_Lat_aftr = Lons[Lons_Diff_Idx]+180., 180.-Lons[Lons_Diff_Idx+1]
            else:
                Lons_MaxIdx = np.copy(Lons_Diff_Idx)
                Lons_MinIdx = Lons_Diff_Idx+1                
                Add_Lon_bfr, Add_Lon_aftr = [180], [-180]
                Add_Lat_bfr, Add_Lat_aftr = 180.-Lons[Lons_Diff_Idx], Lons[Lons_Diff_Idx+1]+180.
            Slope = (Lats[Lons_MaxIdx]-Lats[Lons_MinIdx])/(360-(Lons[Lons_MaxIdx]-Lons[Lons_MinIdx]))
            Lons_bfr = np.hstack([Lons[0:Lons_Diff_Idx+1],Add_Lon_bfr])
            Lats_bfr = Lats[0:Lons_Diff_Idx+1]
            Lats_bfr = np.hstack([Lats_bfr,[Lats_bfr[-1]-Slope*Add_Lat_bfr]])
            Lons_aftr = np.hstack([Add_Lon_aftr,Lons[Lons_Diff_Idx+1:-1]])
            Lats_aftr = Lats[Lons_Diff_Idx+1:-1]
            Lats_aftr = np.hstack([[Lats_aftr[0]-Slope*Add_Lat_aftr],Lats_aftr])
        elif Lons_Diff_Idx.size>0:
            Lons_bfr = Lons
            Lons_bfr[0] = Lons_bfr[0]+360.
            Lons_aftr = np.array([0,0])
            Lats_bfr = Lats
            Lats_aftr = np.array([0,0])
        else:
            Lons_bfr = Lons
            Lons_aftr = np.array([0,0])
            Lats_bfr = Lats
            Lats_aftr = np.array([0,0])
        return(Lons_bfr, Lons_aftr, Lats_bfr, Lats_aftr)

    Main_GlobalAgePlots(SubductionData_AllTr, topology_filename, rotation_model, raster_filename_base, raster_filename_ext, coastline_features, YEAR, SaveDir, Vel_arrows_grid_spacing_degrees, AgeGridRes)
        
        
def GlobalTrencheSegments(SubductionData_AllTr,SaveDir):
    import os
    import numpy as np
    import cartopy as cpy
    import cartopy.crs as ccrs
    import matplotlib.pyplot as plt

    print('\n  Creating global plots trenches IDs and properties')
    def Main_GlobalTrencheSegments(SubductionData_AllTr,SaveDir):
        SubPlateID_List, SubZoneID_List   = np.unique(SubductionData_AllTr[:,20]), np.unique(SubductionData_AllTr[:,21]) # initialise
        SubZoneID_index, SubPlateID_index = np.arange(SubZoneID_List.size), np.arange(SubPlateID_List.size)
        for Time in np.arange(0,121+1,5): # loop through reconstruction times (at 5 Myr steps)
            print('    Reconstruction time: {0}'.format(Time))
            temp = SubductionData_AllTr[(SubductionData_AllTr[:,0]==Time),:] # extract only relevant reconstructed time data
            Lons, Lats, Ages, ConVel, TreVel, SubVel = temp[:,1], temp[:,2], temp[:,3], temp[:,6], temp[:,11], temp[:,16]
            (SubPlateID, SubPlateIdx, SubPlateID_Txt) = processVars(temp[:,20].astype(int), Lats, Lons) # get mask and masked lists for idx of sub. plates and texts and locations
            (SubZoneID, SubZoneIdx, SubZoneID_Txt) = processVars(temp[:,21].astype(int), Lats, Lons) # get mask and masked lists for idx of sub. zone and texts and locations
            PlotMaps(Lons, Lats, SubPlateIdx, SubZoneIdx, Ages, SubVel, TreVel, ConVel, Time, 
                       SubPlateID, SubPlateID_Txt, SubZoneID, SubZoneID_Txt, SaveDir) # plot
        
    def processVars(IDs, Lats, Lons):
        ID_Txt, Idx = [], np.zeros(Lons.size)  # initialise
        for i, ID in enumerate(np.unique(IDs)): # loop through unique IDs
            Idx[(IDs == ID)] = i # find all indices of ID in IDs
            MinLats = Lats[(IDs == ID)].min() # find minimum latitude of ID
            MinLons = Lons[(IDs == ID) & (Lats == MinLats)].min()# find the corresponding longtitude
            ID_Txt.append([MinLons,MinLats,str(ID)]) # add to list of texts and locations
        return(np.unique(IDs), Idx, ID_Txt)

    def PlotMaps(Lons, Lats, SubPlateIdx, SubZoneIdx, Ages, SubVel, TreVel, ConVel, Time, 
                 SubPlateID, SubPlateID_Txt, SubZoneID, SubZoneID_Txt, SaveDir):
        fig = plt.figure(figsize=(30, 10))
        ### Create plots
        CreateAxes(1, Lons, Lats, SubPlateIdx, 'jet'         , []      , 'PlateID'  ,SubPlateID, SubPlateID_Txt)
        CreateAxes(2, Lons, Lats, SubZoneIdx , 'jet'         , []      , 'TrenchID' ,SubZoneID , SubZoneID_Txt )
        CreateAxes(3, Lons, Lats, Ages       , 'gist_earth_r', [0,200] , 'Age'      ,[]        , []            )
        CreateAxes(4, Lons, Lats, SubVel     , 'gist_ncar'   , [-10,10], 'Sub. Vel.',[]        , []            )
        CreateAxes(5, Lons, Lats, TreVel     , 'gist_ncar'   , [-10,10], 'Tre. Vel.',[]        , []            )
        CreateAxes(6, Lons, Lats, ConVel     , 'gist_ncar'   , [-10,10], 'Con. Vel.',[]        , []            )
        ### figure parameters
        fig.suptitle('Global trenches at {0}'.format(int(Time)))
        fig.tight_layout()
        ### save
        Path = os.path.join(SaveDir,'GlobalTrenchSegments')
        if not os.path.exists(Path):
            os.makedirs(Path)
        fig.savefig(os.path.join(Path,'TrenchSegments_{0}Myr.png'.format(int(Time))), dpi=300, facecolor='white', transparent=False) # save as .png
        fig.savefig(os.path.join(Path,'TrenchSegments_{0}Myr.svg'.format(int(Time))), dpi=300, facecolor='white', transparent=False) # save as .svg
        plt.close()

    def CreateAxes(loc, X, Y, C, CM, Clim, Cname, Clabels, Texts):
        fig = plt.gcf()
        if Cname=='PlateID' or Cname=='TrenchID': # set parameters for IDs plots
            Cmap, Vmin, Vmax = plt.get_cmap(CM, np.unique(C).size), 0, np.unique(C).size
        else: # set parameters for other plots
            Cmap, Vmin, Vmax = plt.get_cmap(CM), min(Clim), max(Clim)
        ax = fig.add_subplot(2, 3, loc, projection=ccrs.Robinson(central_longitude=180)) # plot subfigure
        ax.add_feature(cpy.feature.COASTLINE, edgecolor=[0,0,0,.2]) # add coastlines
        ax.set_global() # global view plot
        ScatterPlot = ax.scatter(X, Y, c=C, s=1, vmin=Vmin, vmax=Vmax, cmap=Cmap, transform=ccrs.PlateCarree()) # plot the data
        if Cname=='PlateID' or Cname=='TrenchID': # additional plot parameters for IDs plots
            cb = fig.colorbar(ScatterPlot, ticks = np.unique(C)) # mask the colorbar
            cb.ax.set_yticklabels(Clabels) # mask the colorbar labels
            for i in range(len(Texts)): # add segment numbers on the figure
                plt.text(Texts[i][0], Texts[i][1], Texts[i][2], transform=ccrs.PlateCarree())          
        else:
            cb = fig.colorbar(ScatterPlot) # plot the colorbar in other plots
        cb.set_label(Cname)
        
    Main_GlobalTrencheSegments(SubductionData_AllTr,SaveDir)

if __name__ == '__main__':
    import gc, os, pygplates
    import numpy as np
    gc.collect()
    def main():
        ### Plots for 2019 compilation
        YEAR = 2019
        WorkDir = os.path.join(os.getcwd(),'2019_Compilation')
        Data_Folder = os.path.join(WorkDir,'Muller2019_GPlatesModel')
        raster_filename_base = os.path.join(Data_Folder,'Age_Raster_Files','Muller_etal_2019_Tectonics_v2.0_AgeGrid')
        coastline_features = pygplates.FeatureCollection(os.path.join(Data_Folder,"Coastlines","Global_coastlines_2019_v1_low_res.shp"))
        print('\nPlotting data for 2019 compilation')
        set_variables(WorkDir,Data_Folder,raster_filename_base,coastline_features,YEAR)
        ### Plots for 2019 compilation
        YEAR = 2016
        WorkDir = os.path.join(os.getcwd(),'2016_Compilation')
        Data_Folder = os.path.join(WorkDir,'Muller2016_GPlatesModel')
        raster_filename_base = os.path.join(Data_Folder,'Age_Raster_Files','EarthByte_AREPS_v1.15_Muller_etal_2016_AgeGrid')
        coastline_features = pygplates.FeatureCollection(os.path.join(Data_Folder, "Global_EarthByte_230-0Ma_GK07_AREPS_Coastlines.gpml"))
        print('\n\nPlotting data for 2016 compilation')
        set_variables(WorkDir,Data_Folder,raster_filename_base,coastline_features,YEAR)

    def set_variables(WorkDir,Data_Folder,raster_filename_base,coastline_features,YEAR):
        SaveDir = os.path.join(WorkDir,'Plots')
        raster_filename_ext = 'nc'
        rotation_model = pygplates.RotationModel(os.path.join(Data_Folder, "MergedRot.rot"))
        topology_filename = os.path.join(Data_Folder, "MergedGPML.gpml")
        SubductionData_AllTr  = np.loadtxt(os.path.join(WorkDir,'SubductionData_AllTr.csv') ,delimiter=',', skiprows=1) # Load 'all trenches' data
        SubductionData_MainTr = np.loadtxt(os.path.join(WorkDir,'SubductionData_MainTr.csv'),delimiter=',', skiprows=1) # Load 'main trenches' data

        
        GlobalTrencheSegments(SubductionData_AllTr,SaveDir)
        GlobalAgePlots(SubductionData_AllTr, topology_filename, rotation_model, raster_filename_base, raster_filename_ext, coastline_features, YEAR, SaveDir, 10., 5.)
        PlotTrenchThroughTime(SubductionData_MainTr,SaveDir)
        
    main()