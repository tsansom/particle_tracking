# -*- coding: utf-8 -*-
"""
Created on Mon Dec  4 11:17:21 2017

@author: tsansom
"""

import pandas as pd
import os
import matplotlib.pyplot as plt
from shapely.geometry.polygon import LinearRing
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import matplotlib.dates as mdates
from matplotlib.ticker import MultipleLocator
from matplotlib import gridspec
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
import matplotlib.gridspec as gridspec
import sys
import tbtools as tbt
import numpy as np
from shapely.geometry import Point
from matplotlib.mlab import griddata
from calendar import month_name
import seaborn as sns
from subprocess import run, PIPE
import errno

#this function will get the absolute path to ffmpeg
def get_ffmpeg_path():
    if os.sys.platform[:3] == 'win':
        fil = 'ffmpeg.exe'
        ffmpeg = run(['where', fil], stdout=PIPE)
    else: #<- need to test 'which' on mac
        fil = 'ffmpeg'
        ffmpeg = run(['which', fil], stdout=PIPE)
    if ffmpeg.returncode == 0:
        return ffmpeg.stdout.decode('ascii').rstrip()
    else:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), fil)

#this function is used to set up the salinity mask
def inpolygon(polygon, xp, yp):
    return np.array([Point(x,y).intersects(polygon) for x,y in zip(xp,yp)],dtype=np.bool)

def make_plot():
    #create the figure
    fig = plt.figure(figsize=(14, 8))
    #use gridspec for subplots
    gs = gridspec.GridSpec(2, 3, height_ratios=[3, 2])
    #create the three maps with miller projection
    ax1 = fig.add_subplot(gs[0], projection=ccrs.Miller())
    ax2 = fig.add_subplot(gs[1], projection=ccrs.Miller())
    ax3 = fig.add_subplot(gs[2], projection=ccrs.Miller())
    #create the three inflow plots and share the y axis so they are all on the same scale
    ax4 = fig.add_subplot(gs[3])
    ax5 = fig.add_subplot(gs[4], sharey=ax4)
    ax6 = fig.add_subplot(gs[5], sharey=ax4)
    return fig, ax1, ax2, ax3, ax4, ax5, ax6

#define the wet, dry, and average directories
ptrac_dirs = {'wet': ['1990-07', '2004-05', '2007-07'],
              'dry': ['1987-04', '2014-06', '2014-04'],
              'avg': ['2008-08', '1995-04', '1994-05']}

#this sets the x and y limits of each subplot
extent = (-97.3, -96.8, 27.8, 28.3)

#get current working directory
cwd = os.getcwd()

#point to the shapefile
shp_path = os.path.join(cwd, 'data', 'shapefile', 'CBclosed.shp')

#read the node coordinates
if 'coords' not in vars():
    print('Reading node coordinates\n')
    coords = tbt.read.coords(os.path.join(cwd, 'data', 'model_output', 'input'), 14)
else:
    print('Node coordinates already read in\n')

#clip the coordinates to only contain points inside the extent
#this is used for contouring salinity
avesalD_nodes = coords[(coords['lon'] > extent[0]) &
                       (coords['lon'] < extent[1]) &
                       (coords['lat'] > extent[2]) &
                       (coords['lat'] < extent[3])].index.tolist()

#create the salinity mask
if 'mask' not in vars():
    print('Preparing salinity mask\n')
    #get all nodes within the extent, and get their lons and lats
    avesalD_coords = coords.loc[avesalD_nodes]
    lon = np.array(avesalD_coords['lon'])
    lat = np.array(avesalD_coords['lat'])
    #create an array of lat and lon ranges for gridding
    loni = np.linspace(min(lon), max(lon), 150)
    lati = np.linspace(min(lat), max(lat), 150)
    #create mesh grid
    glon, glat = np.meshgrid(loni, lati)
    #read the shapefile
    shp = Reader(shp_path)
    #get the geometries from the shapefile
    geoms = list(shp.geometries())
    max_area = 0
    max_area_id = 0
    #find the largest polygon in the shapefile, the inside of which will be the bay
    for i in range(len(geoms)):
        if geoms[i].area > max_area:
            max_area = geoms[i].area
            max_area_id = i
            polygon = geoms[i]
    #create the mask based on which points are inside the largest polygon
    mask = inpolygon(polygon, glon.ravel(), glat.ravel())
    #reshape the mask to the same size as the meshgrid
    mask = mask.reshape(glon.shape)
else:
    print('Salinity mask already created\n')

if 'data' not in vars():
    data = {}

#loop through the three scenarios (wet, dry, avg)
    for key in ptrac_dirs.keys():
        #create a blank subdictionary for each scenario
        data[key] = {}
        #loop through the months for each scenario
        for p_dir in ptrac_dirs[key]:
            print('Reading data for {}'.format(p_dir))
            #set up file paths to the different data sources
            ptrac_path = os.path.join(cwd, 'data', 'model_output', p_dir)
            lon_path = os.path.join(ptrac_path, '{}_lon.csv'.format(p_dir))
            lat_path = os.path.join(ptrac_path, '{}_lat.csv'.format(p_dir))
            outflw1_dir = ptrac_path
            input_path = os.path.join(ptrac_path, 'input')
            avesalD_path = os.path.join(ptrac_path, 'avesalD.w')
            #now read in the data
            print('\tReading particle longitude')
            lon_data = pd.read_csv(lon_path, parse_dates=True, index_col=0)
            print('\tReading particle latitude')
            lat_data = pd.read_csv(lat_path, parse_dates=True, index_col=0)
            print('\tReading outflw1')
            outflw1_data = tbt.read.outflw1(outflw1_dir)
            for k in outflw1_data:
                #calculate the u and v components of velocity from the speed and direction
                outflw1_data[k]['u'] = -1 * outflw1_data[k]['velocity'] * np.sin(np.deg2rad(outflw1_data[k]['direction']))
                outflw1_data[k]['v'] = -1 * outflw1_data[k]['velocity'] * np.cos(np.deg2rad(outflw1_data[k]['direction']))
            print('\tReading daily average salinity\n')
            avesalD_data = tbt.read.avesalD(avesalD_path)
            #need to clip the daily average salinity to within the extent (saves compute time)
            avesalD_data = avesalD_data[avesalD_nodes]
            #store all the datas
            data[key][p_dir] = {'lat': lat_data,
                                'lon': lon_data,
                                'outflw1': outflw1_data,
                                'avesalD': avesalD_data}
else:
    print('All particle data already read in\n')

if 'inflow' not in vars():
    print('Reading inflow')
    inflow_path = os.path.join(cwd, 'data', 'model_input')
    files = ['inflw-aransas-7715', 'inflw-mission-7715', 'inflw-copano-7715', 'inflw-cavasso-7715']
    ctr = 0
    for f in files:
        #need to set up the dataframe on the first iteration
        if ctr == 0:
            inflow = tbt.read.inflow(os.path.join(inflow_path, f))
            inflow.columns = [f.split('-')[1]]
        #subsequent iterations append the data to the dataframe
        else:
            tmp = tbt.read.inflow(os.path.join(inflow_path, f))
            inflow[f.split('-')[1]] = tmp
        ctr += 1
    #scale the inflows units from cfs to 1000 cfs
    inflow['sum'] = inflow.sum(axis=1) / 1000

else:
    print('Inflow data already read in\n')

#this creates the ring to highlight the project area
lats = [28.142815, 28.142622, 28.128572, 28.128588]
lons = [-97.057931, -97.041549, -97.041542, -97.057923]
ring = LinearRing(list(zip(lons,lats)))

#create the big loop for scenarios
#for key, Key in zip(['wet'], ['Wet']):
#for key, Key in zip(['dry', 'avg'], ['Dry', 'Average']):
for key, Key in zip(['wet', 'dry', 'avg'], ['Wet', 'Dry', 'Average']):
    #create the figure and axes
    fig, ax1, ax2, ax3, ax4, ax5, ax6 = make_plot()
    #loop through the map axes
    for ax in [ax1, ax2, ax3]:
        #set the spatial bounds of the map
        ax.set_extent(extent, ccrs.Miller())
        #add the shapefile
        ax.add_geometries(shp.geometries(), ccrs.Miller(), facecolor='none', edgecolor='black')
        #add the project area ring
        ax.add_geometries([ring], ccrs.Miller(), facecolor='none', edgecolor='black')
    fig.tight_layout()
    #adjust the figure to make room for labels and titles
    fig.subplots_adjust(left=0.07, bottom=0.07, top=0.95)
    #create the title
    fig.suptitle(r'{} Months'.format(Key), fontsize=24, x=0.01, y=0.99,
                 horizontalalignment='left', verticalalignment='top')
    #get the check nodes and coordinates (same for all)
    cknodes = coords.loc[list(map(int, list(data['wet']['1990-07']['outflw1'].keys())))]
    #create empty dictionaries to store parts of the plots that will be updated
    index = {}
    time_text = {}
    sctr = {}
    quiv = {}
    contf = {}
    axes = {}
    vline = {}
    #create the levels and ticks for salinity contours and colorbar
    levels = np.arange(0, 30.5, 0.5)
    ticks = np.arange(0, 30.5, 2)
    #create separate lists for the two types of axes (map vs plot)
    ptrac_axes = [ax1, ax2, ax3]
    inflow_axes = [ax4, ax5, ax6]
    #get the list of months for the current scenario
    months = list(data[key].keys())
    #need this counter to set the ylabels of the plots on the left
    ctr = 0
    #loop through the axes to create initial figure
    for pax, iax, month in zip(ptrac_axes, inflow_axes, months):
        #get the first index
        index[month] = data[key][month]['lon'].index[0]
        #create the time box in the upper left of the maps
        time_text[month] = pax.text(0.03, 0.93, index[month], bbox=dict(facecolor='lightgray', alpha=1),
                 fontsize=14, weight='medium', transform=pax.transAxes)
        #set up the scatter by scattering the initial particle position
        #note: there is no particle position at time zero but this is a placeholder
        #so it can be updated in the animation loop
        sctr[month] = pax.scatter(data[key][month]['lon'].iloc[0], data[key][month]['lat'].iloc[0], marker='.',
            facecolor='yellow', edgecolor='k', linewidth=0.2, s=25)
        #calculate the u and v velocity components from speed and direction at check nodes
        u = [data[key][month]['outflw1'][k]['u'][0] for k in data[key][month]['outflw1']]
        v = [data[key][month]['outflw1'][k]['v'][0] for k in data[key][month]['outflw1']]
        #quiver the velocity vectors
        quiv[month] = pax.quiver(cknodes['lon'], cknodes['lat'], u, v, scale=10, facecolor='white',
            edgecolors='black', linewidth=0.3)
        #get the salinity for the initial time step
        sal = np.array(data[key][month]['avesalD'].iloc[0])
        #grid the salinity for the map
        sali = griddata(lon, lat, sal, loni, lati, interp='linear')
        #apply the mask
        sali.mask = ~mask
        #plot the filled contours on the maps
        contf[month] = pax.contourf(loni, lati, sali, 50, zorder=0, cmap=plt.cm.seismic,
             vmin=0, vmax=30, levels=levels, extend='max', alpha=0.7)
        #need to get the data range for the current month to plot inflow
        yr = int(month.split('-')[0])
        mth = int(month.split('-')[1])
        #slice the inflow for the specified month
        inflow_tmp = inflow[(inflow.index.year==yr) & (inflow.index.month==mth) & (inflow.index.day<=29)]
        #convert the date index to pydatetime (because plot_date doesn't work with numpy datetimes)
        pydt = inflow_tmp.index.to_pydatetime()
        #plot the inflow
        iax.plot_date(pydt, inflow_tmp['sum'], '-')
        #format the xticks to be display day number only
        iax.xaxis.set_major_formatter(mdates.DateFormatter('%d'))
        #set minor ticks to days
        iax.xaxis.set_minor_locator(MultipleLocator(1))
        #set the plot limits and labels
        iax.set(xlim=[inflow_tmp.index.min(), inflow_tmp.index.max()],
                      xlabel='{} {}'.format(month_name[mth], str(yr)))
        #if working on the leftmost plot (ctr=0) set the ylabel
        #if working on the two right panels, set the visibility of yticklabels to False
        if ctr != 0:
            plt.setp(iax.get_yticklabels(), visible=False)
        else:
            iax.set(ylabel='Total Daily Inflow (1,000 cfs)')
        #iterate the counter
        ctr += 1
        #make a red vertical line for the current time step
        vline[month] = iax.axvline(index[month], c='r')
        #store the axes for updating
        axes[month] = [pax, iax]
    #add a new axes to the plot for the colorbar
    cbar_ax = fig.add_axes([0.05, 0.455, 0.01, 0.465])
    #make the colorbar with ticks and labels on the left hand side
    fig.colorbar(contf[month], cax=cbar_ax, ticks=ticks, label='Salinity (ppt)')
    cbar_ax.yaxis.set_ticks_position('left')
    cbar_ax.yaxis.set_label_position('left')
    #now begin the animation updater
    def update(i):
        #loop through the months in the current scenario (wet, dry, average)
        for month in months:
            #get the current time index
            index[month] = data[key][month]['lon'].index[i]
            #update the time in the display box
            time_text[month].set_text(index[month])
            #get the position of the particles at the previous time step
            offsets = sctr[month].get_offsets()
            #update the positions for the current time step
            offsets[:, 0] = data[key][month]['lon'].loc[index[month]]
            offsets[:, 1] = data[key][month]['lat'].loc[index[month]]
            #update the particle positions
            sctr[month].set_offsets(offsets)
            #velocity is updated every hour, so need to check if the current time index
            #(index[month]) is in the outflw1 index, if so, update the velocity vectors
            if index[month] in data[key][month]['outflw1'][str(cknodes.index[0])].index:
                #get the velocity components for the current time step (if updating)
                u = [data[key][month]['outflw1'][k]['u'].loc[index[month]] for k in data[key][month]['outflw1']]
                v = [data[key][month]['outflw1'][k]['v'].loc[index[month]] for k in data[key][month]['outflw1']]
                #update the velocity vectors
                quiv[month].set_UVC(u, v)
            #salinity contours are updated once per day (at 00:00), so need to check if the
            #current time index is in the avesalD index, if so, update the salinity contours
            if index[month] in data[key][month]['avesalD'].index:
                #get the salinity for the current time step (if updating)
                sal = np.array(data[key][month]['avesalD'].loc[index[month]])
                #grid the salinity
                sali = griddata(lon, lat, sal, loni, lati, interp='linear')
                #apply the mask
                sali.mask = ~mask
                #remove all of the old salinity contours (haven't found a better way yet)
                for c in contf[month].collections: c.remove()
                #update the contours on the map
                contf[month] = axes[month][0].contourf(loni, lati, sali, 50, zorder=0, cmap=plt.cm.seismic,
                     vmin=0, vmax=30, alpha=0.7)
            #update the red vertical line's position to the current time step
            vline[month].set_xdata(index[month])
        #print the frame number in increments of 10
        if (i+1) % 10 == 0:
            print('Frame {} of 1345'.format(i+1))
    #create a quiverkey to show units of the velocity vectors
    qk = plt.quiverkey(quiv[month], 0.9, 0.96, 1, 'Velocity Vectors: ' + r'$1 \frac{ft}{s}$',
                       coordinates='figure', fontproperties={'size': 14}, labelpos='W')
    #set up the animator
    anim = FuncAnimation(fig, update, frames=1345, repeat=False)
    #adjust the frame rate here (7 fps -> 3:12 minutes, 20 fps -> 1:07 minutes)

    ##########################################################################
    #uncomment everything below this line if you don't want to save the animations
    fps = 20
    #choose the video converter (ffmpeg)
    plt.rcParams['animation.ffmpeg_path'] = get_ffmpeg_path()
    #create the writer (this one stores each frame as temporary pngs and then stitches them
    #together at the end)
    writer = animation.FFMpegFileWriter(fps=fps)
    #make sure the output directory exists, if not make it
    if not os.path.exists(os.path.join(cwd, 'results')):
        os.mkdir(os.path.join(cwd, 'results'))
    #create the animation file name
    output_path = os.path.join(cwd, 'results', '{}_months_{}_fps.mp4'.format(key, fps))
    #save the animation
    anim.save(output_path, writer=writer)
    #stop the animation
    anim.event_source.stop()
    #delete the animation object (need to do this to close the animated plot)
    del anim
    #close the figure
    plt.close()
