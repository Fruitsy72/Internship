##This section is to load in the modules that we need
from netCDF4 import Dataset
import matplotlib 
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import glob
import os 
import numpy as np
from matplotlib.cm import get_cmap
import cartopy.crs as crs
from cartopy.feature import NaturalEarthFeature
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords)

################### Function to truncate color map ###################
def truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='jet', minval=0.0, maxval=1.0, n=100)'''    
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    return new_cmap

cmap_mod = truncate_colormap(minval=.2, maxval=1.0)  # calls function to truncate colormap 

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
                
        URaw = getvar(ncfile, "U10")
        U_Wind = fileid.variables["U10"][0,:,:]
#        print(np.shape(U_Wind))
        
        VRaw = getvar(ncfile, "V10")
        V_Wind = fileid.variables["V10"][0,:,:]
#        print(np.shape(V_Wind))
        
        LonRaw = getvar(ncfile, "XLONG")
        XLon = fileid.variables["XLONG"][0,:,:]
#        print(np.shape(XLon))
        
        LatRaw = getvar(ncfile, "XLAT")
        XLat = fileid.variables["XLAT"][0,:,:]
#        print(np.shape(XLat))
        
        CH4_Ant_Raw = getvar(ncfile, "CH4_ANT")
        CH4_ANT = fileid.variables['CH4_ANT'][0,:,:]
#        print(np.shape(CH4_ANT))
        
        CH4_Bio_Raw = getvar(ncfile, "CH4_BIO")
        CH4_BIO = fileid.variables['CH4_BIO'][0,:,:]
#        print(np.shape(CH4_BIO))
        
        CH4_Bck_Raw = getvar(ncfile, "CH4_BCK")
        CH4_BCK = fileid.variables['CH4_BCK'][0,:,:]
#        print(np.shape(CH4_BCK))
        
        CH4_Tst_Raw = getvar(ncfile, "CH4_TST")
        CH4_TST = fileid.variables['CH4_TST'][0,:,:]
#        print(np.shape(CH4_TST))
        
        CO2_Tst_Raw = getvar(ncfile, "CO2_TST")
        CO2_TST = fileid.variables['CO2_TST'][0,:,:]
#        print(np.shape(CO2_TST))
        
        CH4 = (CH4_BIO[0,:,:] + CH4_ANT[0,:,:] - CH4_BCK[0,:,:] + (CH4_TST[0,:,:] - CH4_BCK[0,:,:]) + (CO2_TST[0,:,:] - CH4_BCK[0,:,:]))
        
        # Get the latitude and longitude points#
        lats, lons = latlon_coords(CH4_Ant_Raw)
                
#        print(lats)
#        print(lons)
        
        # Get the cartopy mapping object
        cart_proj = get_cartopy(CH4_Ant_Raw)
        
        # Create a figure
        fig = plt.figure(figsize=(4,3))
        
        # Set the GeoAxes to the projection used by WRF
        ax = plt.axes(projection=cart_proj)
        
        ispeed = 1

        if ispeed: 
            ret = ax.projection.transform_points(crs.PlateCarree(), np.array(lons),
                np.array(lats)) # This method only accept ndarray!
            xx = ret[..., 0]
            yy = ret[..., 1]
 
# if you want full color plot comment out the lines with cropped in it            
 
#            cropped = (slice(40, 120, None), slice(160, 260, None))
            
#            m = plt.contourf(xx(cropped), yy(cropped), to_np(CH4),levels = 100,cmap=cmap_mod)

            m = plt.contourf(xx, yy, to_np(CH4),levels = 100,cmap=cmap_mod)

        else:
#            m = plt.contourf(to_np(XLon[cropped]), to_np(XLat[cropped]), to_np(CH4[cropped]),
#                             levels = 50, transform=crs.PlateCarrree(), cmap=cmap_mod)

             m = plt.contourf(to_np(XLon), to_np(XLat), to_np(CH4), levels = 100,
                              transform=crs.PlateCarree(),
                              cmap=cmap_mod)
        
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
        plt.title(str_xhu,fontsize=8,y=.979900091)
        plt.title(str_xhu,fontsize=8,pad=0.091)
        
        # This is adding a color bar and determining its shape and location on the plot.
        cbaxes = fig.add_axes([0.2, 0.15, 0.6, 0.021]) 

        cb=plt.colorbar(cax= cbaxes, shrink=.68, orientation='horizontal')
        cb.ax.tick_params(labelsize=6,direction="in")
        
        cb.set_label("CH4",labelpad=-8,fontsize=8, x=1.05, rotation=0)
        
         
        # Download and add the states and coastlines
        states = NaturalEarthFeature(category="cultural", scale="50m",
                                     facecolor="none",
                                     name="admin_1_states_provinces_lines")
        ax.add_feature(states, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
        
        country_borders = NaturalEarthFeature(category="cultural", scale="50m",
                                     facecolor="none",
                                     name="admin_0_boundary_lines_land")
        ax.add_feature(country_borders, linewidth=.5, edgecolor="black")
        ax.coastlines('50m', linewidth=0.8)
        
        # Set the map bounds
        ax.set_xlim(cartopy_xlim(CH4_Ant_Raw))
        ax.set_ylim(cartopy_ylim(CH4_Ant_Raw))
        
        #ax.set_extent([-112, -92, 26, 42], crs=crs.PlateCarree())        
        
        skip = (slice(None, None, 10), slice(None, None, 10))
        
        ax.quiver(XLon[skip], XLat[skip], U_Wind[skip], V_Wind[skip], headwidth=2, transform = crs.PlateCarree())        
        
        figname='wrfout_d01_CH4_BIO+ANT+Wet+agwaste_python_'+str(i)+'.png'
        print(figname)
        plt.savefig(os.path.join('/scratch/08605/tg879045/figures',figname),dpi=300,format='png', bbox_inches = 'tight',pad_inches = 0)
        plt.clf()#
        plt.close()
    i = i+1