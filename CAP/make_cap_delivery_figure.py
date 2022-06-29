import sys
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from shapely.geometry import Point, Polygon, LineString
import load_data as ld
from sanker import Sanker
from osgeo import gdal
from apng import APNG

if len(sys.argv) > 1:
  if int(sys.argv[1]) == 1:
    plot_type = 'monthly'
  else:
    plot_type = 'annual'
  plot_gif = False
  if len(sys.argv) > 2:
    try:
      plot_gif = bool(int(sys.argv[2]))
    except:
      plot_gif = False
else:
  plot_type = 'annual'
  plot_gif = False
##############################################
###MAP PARAMETERS
##############################################
#map projection
projection_string = 'EPSG:3857'
##data to load satellite raster background
raster_name_pt1 = 'LC09_CU_'
raster_name_pt2 = '_02_SR'
raster_band_list = ['B4', 'B3', 'B2']
raster_id_list = {}
raster_id_list['006012'] = ['20220610_20220615','20220518_20220523']
raster_id_list['006013'] = ['20220610_20220615','20220603_20220608']
raster_id_list['006014'] = ['20220603_20220608','20220610_20220615']
raster_id_list['007012'] = ['20220612_20220617','20220603_20220608']
raster_id_list['007013'] = ['20220612_20220617','20220603_20220608']
raster_id_list['007014'] = ['20220612_20220617','20220603_20220608']
raster_id_list['007015'] = ['20220612_20220617',]
raster_id_list['008013'] = ['20220612_20220617',]
raster_id_list['008014'] = ['20220612_20220617',]
#plotting range (time)
year_start = 2018
year_end = 2023
num_years = year_end - year_start
#plotting range (space) - epsg coordinates
plot_xlim = [-12670000.0, -12321000.0]
plot_ylim = [3755000.0, 4030000.0]
yrange = plot_ylim[1] - plot_ylim[0]
xrange = plot_xlim[1] - plot_xlim[0]
sankey_colors = {}
sankey_colors['Municipal'] = 'darkorange'
sankey_colors['Irrigation'] = 'forestgreen'
sankey_colors['Tribal'] = 'maroon'
sankey_colors['Recharge'] = 'darkturquoise'
map_colors = {}
map_colors['Municipal'] = 'goldenrod'
map_colors['Irrigation'] = 'teal'
map_colors['Tribal'] = 'indianred'
map_colors['Recharge'] = 'darkturquoise'

##############################################
###GIS DATA
##############################################
#contractor names - bridge between CAP delivery data names and GIS shapefile names
key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge, total_deliveries = ld.load_key_dicts(year_start, year_end)
#load contractor shapefiles
epsg_use = 3857
cws = gpd.read_file('CWS_Service_Area/CWS_Service_Area.shp')#M&I
irrdst = gpd.read_file('Irrigation_District/Irrigation_District.shp')#Agriculture
tribal = gpd.read_file('ALRIS_amerind/ALRIS_amerind.shp')#Tribal Lands
amas = gpd.read_file('AMA_and_INA/AMA_and_INA.shp')#Active Management Areas
gwstorage = gpd.read_file('Underground_Storage_Facilities/Underground_Storage_Facilities.shp')#GW recharge basins
canals = gpd.read_file('Cen_AZ_Proj/Cen_AZ_Proj.shp')#shape of distribution canal
cws = cws[pd.notna(cws['geometry'])]
cws = cws.to_crs(epsg = epsg_use)
irrdst = irrdst.to_crs(epsg = epsg_use)
tribal = tribal.to_crs(epsg = epsg_use)
amas = amas.to_crs(epsg = epsg_use)
gwstorage = gwstorage.to_crs(epsg = epsg_use)
basin_facilities = gpd.sjoin(gwstorage, amas, how = 'inner', op = 'within')#link gw storage facilities to active management areas
#shape of distribution path (CAP canal)
canals = canals.to_crs(epsg = epsg_use)

##############################################
###DELIVERY/FLOW DATA
##############################################
#online data is read by get_cap_delivery_data.py
#read in CAP contractor deliveries from .csv file 
district_deliveries = pd.read_csv('deliveries_by_district.csv', index_col = 0)
district_deliveries.index = pd.to_datetime(district_deliveries.index)
district_deliveries['year_num'] = pd.DatetimeIndex(district_deliveries.index).year
annual_deliveries = pd.DataFrame()
for x in district_deliveries:
  annual_deliveries[x] = district_deliveries.groupby(['year_num'])[x].aggregate(sum)
if plot_type == 'annual':
  deliveries_use = annual_deliveries
  scale_factor = 15000000.0
else:
  deliveries_use = district_deliveries
  scale_factor = 1750000.0

##############################################
###CREATE MAPS
##############################################
monthnum_datetime = 1
yearnum_datetime = 2018
off_ramp_buffer = 50000#how much distance (map units) between the delivery point/centroid and where the delivery breaks from the main path
for dt_index, del_row in deliveries_use.iterrows():  
  #Make new sankey map for each timestep (either monthly/annual)  
  sankey_plot = Sanker()
  #load raster files as map background - comment out to reduce file size
  #sankey_plot.load_batch_raster('az_satellite_photos/', raster_id_list, raster_name_pt1, raster_name_pt2, raster_band_list, projection_string, max_bright = (6000.0, 15000.0), use_gamma = 0.8, upscale_factor = 0.1)
  #initialize the shape of the distribution 'path' - divided into individual 'segments'
  sankey_plot.create_sankey_path(canals, two_d_smoothing_reach = [-12650000.0, 12550000.0], outlier_values = 1000000.0)  
  canal_vert_range = np.max(sankey_plot.sankey_path_y) - np.min(sankey_plot.sankey_path_y)#total verticle height of the distribution path, in map units
  #location of all the delivery points (contractor shapefile centroids)
  district_centroids = ld.create_district_centroids(district_deliveries, key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge, sankey_plot.sankey_path_x, off_ramp_buffer, epsg_use)
  geometry = [Point(xy) for xy in zip(district_centroids['x_loc'], district_centroids['y_loc'])]
  crs = {'init': 'epsg:4326'}
  district_centroids = gpd.GeoDataFrame(district_centroids, crs = crs, geometry = geometry)
  
  #find total deliveries in the current timestep (i.e., initial flow on the distribution path)  
  total_section_flow = 0.0
  for dict_use, shape_use, column_use in zip([key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge], [cws, irrdst, tribal, basin_facilities], ['CWS_NAME', 'LONG_NAME', 'NAME', 'BASIN_NAME']):
    for x in dict_use:
      counter = 0
      for y in dict_use[x]:
        this_district = shape_use[shape_use[column_use] == y]
        if len(this_district) > 0.0:
          counter = 1
      if counter == 1:
        total_section_flow += del_row[x]
  
  old_inverse_slope = 0.0
  #plotting path colors
  #loop through each 'segment' of the distribution path
  for counter in range(0, len(sankey_plot.sankey_path_y) - 1):
    #find districts that are delivered at this segment of the distribution path
    remaining_districts = district_centroids[district_centroids['sankey_loc'] == counter]
    #find the total deliveries made at this segment of the distribution path in this timestep
    total_deliveries, sorted_deliveries = sankey_plot.custom_delivery_path(remaining_districts, del_row, counter, total_section_flow)
    #calculate starting width, in map units
    total_width = canal_vert_range * total_section_flow / scale_factor
    #plot current segment of the distribution channel 
    #the current segment starts with a width equal to the remaining water from last segment, and ends with a width equal to the remaining water less this path deliveries
    old_inverse_slope_prev = sankey_plot.add_sankey_custom_path(total_width, total_deliveries/total_section_flow, old_inverse_slope, counter)
    
    #plot diversion from the main canal stem - pt. 1
    xy_edge, xy_middle = sankey_plot.add_sankey_custom_diversion(total_width, total_deliveries/total_section_flow, old_inverse_slope, counter)
    bottom_y = min([xy_edge[1], xy_middle[1]])
    dist_y = max([xy_edge[1], xy_middle[1]]) - min([xy_edge[1], xy_middle[1]])
    #plot path to each diversion location - pt. 2
    for index_sorted, row_sorted in sorted_deliveries.iterrows():
      delivery_width = canal_vert_range * row_sorted['deliveries'] / scale_factor
      top_y = bottom_y + delivery_width * dist_y / (canal_vert_range * total_deliveries / scale_factor)
      sankey_plot.add_sankey_custom_diversion2(xy_edge[0], [bottom_y, top_y], delivery_width, row_sorted['x_loc'], row_sorted['y_loc'], sankey_colors[row_sorted['type']])
      bottom_y = top_y * 1.0
    #calculate new total flow and path direction of main distribution path
    total_section_flow -= total_deliveries
    old_inverse_slope = old_inverse_slope_prev * 1.0

  #plot district outlines
  for x in key_dict:
    for y in key_dict[x]:
      this_district = cws[cws['CWS_NAME'] == y]
      if len(this_district) != 0:
        this_district.plot(ax = sankey_plot.ax, facecolor = map_colors['Municipal'], edgecolor = 'black', alpha = 0.4, zorder = 1, linewidth = 0.5)
  for x in key_dict_irr:
    for y in key_dict_irr[x]:
      this_district = irrdst[irrdst['LONG_NAME'] == y]
      if len(this_district) != 0:
        this_district.plot(ax = sankey_plot.ax, facecolor = map_colors['Irrigation'], edgecolor = 'black', alpha = 0.4, zorder = 1, linewidth = 0.5)
  for x in key_dict_tribal:
    for y in key_dict_tribal[x]:
      this_district = tribal[tribal['NAME'] == y]
      if len(this_district) != 0:
        this_district.plot(ax = sankey_plot.ax, facecolor = map_colors['Tribal'], edgecolor = 'black', alpha = 0.4, zorder = 1, linewidth = 0.5)
  basin_facilities.plot(ax = sankey_plot.ax, facecolor =map_colors['Recharge'], edgecolor = 'darkturquoise', alpha = 0.4, zorder = 1, linewidth = 0.75)
  
  #plot outlines for legend/titles
  xl = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * 0.0425
  xr = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * 0.3725
  by = plot_ylim[1] - (plot_ylim[1] - plot_ylim[0]) * 0.46
  uy = plot_ylim[1] - (plot_ylim[1] - plot_ylim[0]) * 0.3325
  
  p2 = Polygon([(xl, by), (xr, by), (xr, uy), (xl, uy)])
  s2 = gpd.GeoSeries(p2)
  df2 = gpd.GeoDataFrame({'geometry': s2, 'df2':[1]})
  df2.crs = {'init' :'epsg:3857'}
  df2.plot(ax = sankey_plot.ax, facecolor = 'beige', edgecolor = 'indianred', linewidth = 3.0, alpha = 0.8)

  xlb = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * 0.0325
  xrb = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * 0.3825
  byb = plot_ylim[0] + (plot_ylim[1] - plot_ylim[0]) * 0.05
  uyb = plot_ylim[0] + (plot_ylim[1] - plot_ylim[0]) * 0.49
  
  p3 = Polygon([(xlb, byb), (xrb, byb), (xrb, uyb), (xlb, uyb)])
  s3 = gpd.GeoSeries(p3)
  df3 = gpd.GeoDataFrame({'geometry': s3, 'df3':[1]})
  df3.crs = {'init' :'epsg:3857'}
  df3.plot(ax = sankey_plot.ax, facecolor = 'beige', edgecolor = 'beige', linewidth = 3.0, alpha = 0.4)
  
  #plot district centroids + map labels
  label_font = 32
  title_font = 52
  subtitle_font = 46
  legend_font = 32
  district_plot_cutoff = 50000.0
  for index_cent, row_cent in district_centroids.iterrows():
    #plot centroids with a size equal to the average annual delivery
    this_district = district_centroids[district_centroids.index == index_cent]
    size_factor = row_cent['delivery'] / np.max(district_centroids['delivery'])
    this_district.plot(ax = sankey_plot.ax, facecolor = map_colors[row_cent['type']], edgecolor = 'black', linewidth = 2.5 * size_factor,  markersize = 300.0 * size_factor, zorder = 4)
    #label the largest contractors
    if row_cent['delivery'] > district_plot_cutoff:
      if row_cent['name'] == 'Phoenix' or row_cent['name'] == 'MSIDD':
        sankey_plot.add_map_label(row_cent['x_loc'] + xrange * 0.0025, row_cent['y_loc'] - yrange * 0.0125, row_cent['name'], fontsize = label_font, ha = 'left', va = 'top')
      else:
        sankey_plot.add_map_label(row_cent['x_loc'] + xrange * 0.0025, row_cent['y_loc'] + yrange * 0.0075, row_cent['name'], fontsize = label_font, ha = 'left', va = 'bottom')

  #plot map titles
  if plot_type == 'monthly':
    current_datetime = datetime(yearnum_datetime, monthnum_datetime, 1, 0, 0)
    print(current_datetime)
    monthnum_datetime += 1
    if monthnum_datetime == 13:
      monthnum_datetime = 1
      yearnum_datetime += 1
  else:
    print(str(dt_index))
  sankey_plot.add_map_label(plot_xlim[0] + xrange * 0.06, plot_ylim[1] - yrange * 0.39875, 'The Central\nArizona Project', fontsize = title_font, ha = 'left')
  if plot_type == 'annual':
    sankey_plot.add_map_label(plot_xlim[0] + xrange * 0.2125, plot_ylim[1] - yrange * 0.55, str(dt_index) + ' Deliveries', fontsize = subtitle_font)
  else:
    sankey_plot.add_map_label(plot_xlim[0] + xrange * 0.2125, plot_ylim[1] - yrange * 0.55, current_datetime.strftime('%b, %Y'), fontsize = subtitle_font)
  
  #plot map labels
  legend_points = pd.DataFrame([0.0,0.0,0.0,0.0,0.0,0.0,0.0])
  legend_x_loc = 0.0675
  legend_y_start = 0.615
  legend_spacing = 0.05
  facecolor_list = ['gray', 'gray', map_colors['Municipal'], map_colors['Irrigation'], map_colors['Tribal'], map_colors['Recharge'], 'steelblue']
  alpha_list = [1.0, 1.0, 0.8, 0.8, 0.8, 0.8, 0.9]
  total_size = [50000, 200000]
  
  #get the geometry of the map labels (plot as map objects so sizes are to scale with map)
  geometry = []
  legend_labels = []
  counter = 0
  for x in range(0, 2):
    geometry.append(Point(plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * legend_x_loc, plot_ylim[1] - (plot_ylim[1] - plot_ylim[0]) * (legend_y_start + x * legend_spacing)))
    legend_labels.append(str(int(total_size[counter]/1000.0)) + ' taf/year')
    counter += 1
  for x in ['Municipal Provider', 'Irrigation District', 'Tribal Lands', 'Recharge', 'Canal']:
    xlb = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * (legend_x_loc - 0.02)
    xrb = plot_xlim[0] + (plot_xlim[1] - plot_xlim[0]) * (legend_x_loc + 0.02)
    uyb =  plot_ylim[1] - (plot_ylim[1] - plot_ylim[0]) * (legend_y_start - 0.01 + counter * legend_spacing)
    byb =  plot_ylim[1] - (plot_ylim[1] - plot_ylim[0]) * (legend_y_start + 0.01 + counter * legend_spacing)
    geometry.append(Polygon([(xlb, byb), (xrb, byb), (xrb, uyb), (xlb, uyb)]))
    legend_labels.append(x)
    counter += 1
    
  #plot legend & legend labels
  legend_points = gpd.GeoDataFrame(legend_points, crs = canals.crs, geometry = geometry)
  counter = 0
  for index, row in legend_points.iterrows():
    this_legend = legend_points[legend_points.index == index]
    if facecolor_list[counter] == 'gray':
      this_legend.plot(ax = sankey_plot.ax, facecolor = facecolor_list[counter], edgecolor = 'black', linewidth = 2.5,  markersize = 300.0 * total_size[counter] /  np.max(district_centroids['delivery']))
    else:
      this_legend.plot(ax = sankey_plot.ax, facecolor = facecolor_list[counter], edgecolor = 'black', linewidth = 2.5,  alpha = alpha_list[counter])
    sankey_plot.add_map_label(plot_xlim[0] + xrange * 0.1225, plot_ylim[1] - (plot_ylim[1] - plot_ylim[0]) * (legend_y_start + counter * legend_spacing), legend_labels[counter], fontsize = legend_font, ha = 'left')
    counter += 1
    
  #plot map inset of the state
  inset_districts = pd.DataFrame()
  for index, row in district_centroids.iterrows():
    if row['name'] == 'Phoenix' or row['name'] == 'Tucson':
      this_district = district_centroids[district_centroids.index == index]
      inset_districts = pd.concat([inset_districts, this_district])
  crs = {'init': 'epsg:4326'}
  inset_districts = gpd.GeoDataFrame(inset_districts, crs = crs, geometry = inset_districts.geometry)
  shapefile_name = 'states/states.shp'
  sankey_plot.add_inset_figure(shapefile_name, canals, inset_districts, plot_xlim, plot_ylim, epsg_use)

  #save file
  if plot_type == 'annual':
    sankey_plot.format_sankey(plot_xlim, plot_ylim, 'sankey_figs/deliveries_' + str(dt_index))
  else:
    sankey_plot.format_sankey(plot_xlim, plot_ylim, 'sankey_figs/deliveries_' + current_datetime.strftime('%Y_%m'))
  del sankey_plot
  plt.close()

##############################################
###CREATE GIF
##############################################
#plot gif of monthly maps
if plot_gif and plot_type == 'monthly':
  figure_base = 'sankey_figs/deliveries_'
  for year_num in range(2018, 2022):
    print(year_num)
    image_list = []
    
    for month_num in range(1, 13):
      file_name = figure_base + str(year_num) + '_' + str(month_num).zfill(2) + '.png'
      image_list.append(file_name)
    kargs = { 'duration': 1.0 }
    APNG.from_files(image_list, delay=1000).save('sankey_figs/gif_' + str(year_num) + '.png')
