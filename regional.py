## This section is to load in the modules that we need
from netCDF4 import Dataset # Our dataset is in a netCDF4 file so we will pull it here
import matplotlib # matplotlib is a library with many plotting utilities in it
matplotlib.use('agg') # This lets matplot draw graphics
import matplotlib.pyplot as plt # This lets us use matplotlib functions by typing plt.whatever 
import matplotlib.colors as colors # Lets us use matplotlib's library of colorbars
import glob # This helps us find files in a specific folder/directory (Dr. Hu's files)
import os # Lets us do things inside file
import numpy as np # Lets us make arrays
from matplotlib.cm import get_cmap # cmap helps us make colormaps
import cartopy.crs as crs # cartopy is a plotting feature that allows us to import maps like the US map
from cartopy.feature import NaturalEarthFeature # This is to allows to bring in state borders and other things
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords) # This allows us to pull in variables so that we can make 
#calculations and other things for our plots

################### Function to truncate color map ###################
def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100)'''    
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    return new_cmap

cmap_mod = truncate_colormap(minval=.2, maxval=1.0)  # calls function to truncate colormap 

#We are now starting with the first iteration (0=1 in code).
#We are digging into the certain directory with that date at any time and having it
#equal to filelist.
#We then begin to iterate each hour.

i = 0
filelist = glob.glob(os.path.join('/scratch/01178/tg804586/Run/CO2_and_otherGHG/WRFV4.3.3/CONUS/wrfchem4.3.3LES3d_Hu2021JGR_CH4NEI2017_Wetchart131_agwasteOce.2022021018/', 'wrfout_d01_2022-*:00:00'))
print(filelist)
print("retrieved all wrfout files")
for filename in sorted(filelist): 
    if i > -1:
        print(filename)
        # Open the NetCDF file
        ncfile = Dataset(filename)
        fileid = Dataset(filename, mode = 'r', format='cdf')
        
# We opened the NetCDF above and now we are digging into the file to find the variables needed.
# getvar is finding us the variable under its "name." We then define the array with [0,:,:].
# We then print the array shape of the variable to make for sure that the arrays are similar.
                
# U_Wind is the wind in the x direction.
# V_Wind is the wind in the y direction.

        URaw = getvar(ncfile, "U10")
        U_Wind = fileid.variables["U10"][0,:,:]
        
        VRaw = getvar(ncfile, "V10")
        V_Wind = fileid.variables["V10"][0,:,:]

#XLon is the longitude data.

#XLat is the latitude data.
        
        LonRaw = getvar(ncfile, "XLONG")
        XLon = fileid.variables["XLONG"][0,:,:]
        
        LatRaw = getvar(ncfile, "XLAT")
        XLat = fileid.variables["XLAT"][0,:,:]
   
        CH4_Ant_Raw = getvar(ncfile, "CH4_ANT")
        CH4_ANT = fileid.variables['CH4_ANT'][0,:,:]
        
        CH4_Bio_Raw = getvar(ncfile, "CH4_BIO")
        CH4_BIO = fileid.variables['CH4_BIO'][0,:,:]
        
        CH4_Bck_Raw = getvar(ncfile, "CH4_BCK")
        CH4_BCK = fileid.variables['CH4_BCK'][0,:,:]
        
        CH4_Tst_Raw = getvar(ncfile, "CH4_TST")
        CH4_TST = fileid.variables['CH4_TST'][0,:,:]
        
        CO2_Tst_Raw = getvar(ncfile, "CO2_TST")
        CO2_TST = fileid.variables['CO2_TST'][0,:,:]
   
# Doing the math
     
        CH4 = (CH4_BIO[0,:,:] + CH4_ANT[0,:,:] - CH4_BCK[0,:,:] + (CH4_TST[0,:,:] - CH4_BCK[0,:,:]) + (CO2_TST[0,:,:] - CH4_BCK[0,:,:]))
        
        # Get the latitude and longitude points
        lats, lons = latlon_coords(CH4_Ant_Raw)              
       
        # Get the cartopy mapping object
        cart_proj = get_cartopy(CH4_Ant_Raw)
        
        # Create a figure
        fig = plt.figure(figsize=(5,4))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
        
# Running iterations to plot out the data. We have the cropped line so that we only plot a small section of
# the CH4 calculations. It is plotted automatically so to have the best most accurate plot on display we only
# plot/caluclate the section we need.
        
        ispeed = 1

        if ispeed: 
            ret = ax.projection.transform_points(crs.PlateCarree(), np.array(lons),
                np.array(lats)) # This method only accept ndarray!
            xx = ret[..., 0]
            yy = ret[..., 1]

            cropped = (slice(40, 120, None), slice(160, 260, None))
            
            m = plt.contourf(xx[cropped], yy[cropped], to_np(CH4[cropped]),levels = 50,cmap=cmap_mod)

        else:
 
             m = plt.contourf(to_np(XLon[cropped]), to_np(XLat[cropped]), to_np(CH4[cropped]),
                     levels = 50, transform=crs.PlateCarrree(), cmap=cmap_mod) 
        
# This is making a blue star on the map where El Reno and Pampa Flogistix is so that
# we can see it on the map.
        
        name =['El Reno']
        lat1 = [35.54122]
        lon1 = [-97.95494]
        plt.plot(lon1, lat1, 'b+-', markersize=2, transform=crs.PlateCarree())
        
        name =['Pampa']
        lat2 = [35.544743]
        lon2 = [-100.963455]
        plt.plot(lon2, lat2, 'b+-', markersize=2, transform=crs.PlateCarree())
            
        str_xhu1 = b''.join(fileid.variables['Times'][0,:]).decode()
        str_xhu2 = "CH4 @"
        str_xhu = str_xhu2+str_xhu1
        plt.title(str_xhu,fontsize=8)
        
        # This is adding a color bar and determining its shape and location on the plot.
        cbaxes = fig.add_axes([0.1, 0.08, 0.84, 0.03]) 

        cb=plt.colorbar(cax= cbaxes, shrink=.68, orientation='horizontal', drawedges=True)
        cb.ax.tick_params(labelsize=6,direction="in")
        
        cb.set_label("CH4 Concentration (ppm)",fontsize=6, y=-100, rotation=0)
                
        # Download and add the states and coastlines
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                     facecolor="none",
                                     name="admin_1_states_provinces_lines")
        ax.add_feature(states, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
  
# This section here will add country borders so that our map looks complete.      
  
        country_borders = NaturalEarthFeature(category="cultural", scale="50m",
                                     facecolor="none",
                                     name="admin_0_boundary_lines_land")
        ax.add_feature(country_borders, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
        
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(CH4_Ant_Raw))
        ax.set_ylim(cartopy_ylim(CH4_Ant_Raw))
  
# This line here is for us to zoom in on the plot so that we only see the area that we need.      
  
        ax.set_extent([-104, -94, 32, 38], crs=crs.PlateCarree())        
        
        skip = (slice(None, None, 8), slice(None, None, 8))
        
# This line makes the arrows that are going to plot our wind.
        
        ax.quiver(XLon[skip], XLat[skip], U_Wind[skip], V_Wind[skip], headwidth=2, transform = crs.PlateCarree())        

# This here gives the figures we made a name then it prints it out so that we know it was made.      

        figname='wrfout_d01_CH4_BIO+ANT+Wet+agwaste_python_'+str(i)+'.png'
        print(figname)
        
# It then saves the figures in the specific directory we chose for it and to save it as a png.  
     
        plt.savefig(os.path.join('/scratch/08605/tg879045/Regional',figname),dpi=300,format='png', bbox_inches = 'tight',pad_inches = 0)
        plt.clf()#
        plt.close()
    i = i+1
    
# We then close everything and end the iteration to so that it can run again for the next hour scan.