import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd
import scipy.interpolate as interp
from osgeo import gdal
import rasterio
from rasterio.plot import reshape_as_image
from skimage import exposure
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from shapely.geometry import Point, Polygon, LineString

class Sanker():

  def __init__(self, nr = 1, nc = 0):
    self.sub_rows = nr
    self.sub_cols = nc
    if self.sub_cols == 0:
      self.fig, self.ax = plt.subplots(self.sub_rows, figsize=(20,16))
      if self.sub_rows == 1:
        self.type = 'single'
        self.ax.grid(False)
      else:
        self.type = '1d'
    else: 
      self.fig, self.ax = plt.subplots(self.sub_rows, self.sub_cols, figsize=(25, 10))
      self.type = '2d'
    self.band_min = np.ones(3)*(-1.0)
    self.band_max = np.ones(3)*(-1.0)
    self.brightness_adj = 2.0

      
  def create_sankey_path(self, sankey_pathline, two_d_smoothing_reach = [100, 0], outlier_values = 0.0):
    counter = 0
    x_list = []
    y_list = []
    for index, row in sankey_pathline.iterrows():  
      line = row['geometry']
      line_list = list(line.coords)
      for i in range(1, len(line_list)):
        start_point = line_list[i-1]
        end_point = line_list[i]
        x_list.append(start_point[0])
        y_list.append(start_point[1])
        x_list.append(end_point[0])
        y_list.append(end_point[1])

    x_list_ar = np.asarray(x_list)
    y_list_ar = np.asarray(y_list)
    sorted_x = np.argsort(x_list_ar)
    values_x = x_list_ar[sorted_x]
    values_y = y_list_ar[sorted_x]
    single_direction_x = [values_x[0],]
    single_direction_y = [values_y[0],]
    current_y = values_y[0] * 1.0
    current_x = values_x[0] * 1.0
    for xx in range(1, len(values_x)):
      if values_x[xx] < two_d_smoothing_reach[0] or values_x[xx] > two_d_smoothing_reach[1]:
        if values_y[xx] < current_y and values_x[xx] > current_x:
          single_direction_x.append(values_x[xx])
          single_direction_y.append(values_y[xx])
          current_y = values_y[xx] * 1.0
          current_x = values_x[xx] * 1.0
      else:
        if values_x[xx] > current_x:
          single_direction_x.append(values_x[xx])
          single_direction_y.append(values_y[xx])
          current_y = values_y[xx] * 1.0
          current_x = values_x[xx] * 1.0
  
    values_x = np.asarray(single_direction_x)
    values_y = np.asarray(single_direction_y)
    sorted_x = np.argsort(values_x)
    x_range = np.max(values_x) - np.min(values_x)
    self.sankey_path_x = np.linspace(np.min(values_x), np.max(values_x) - x_range * 0.00, 20)
    yspline = interp.InterpolatedUnivariateSpline(values_x[sorted_x],values_y[sorted_x], k=3, ext='zeros')
    self.sankey_path_y = yspline(self.sankey_path_x)
    self.sankey_path_x = self.sankey_path_x[self.sankey_path_y>outlier_values]
    self.sankey_path_y = self.sankey_path_y[self.sankey_path_y>outlier_values]
    
  def custom_delivery_path(self, current_deliveries, del_row, counter, section_flow = 0.0):
    total_deliveries = 0.0
    ind_deliveries = []
    ind_x = []
    ind_y = []
    del_type = []
    for index, row in current_deliveries.iterrows():
      district_name = row['name']
      total_deliveries += del_row[district_name]
      ind_deliveries.append(del_row[district_name])
      ind_x.append(row['x_loc'])
      ind_y.append(row['y_loc'])
      del_type.append(row['type'])
      
    this_step_deliveries = pd.DataFrame()
    this_step_deliveries['deliveries'] = ind_deliveries
    this_step_deliveries['x_loc'] = ind_x
    this_step_deliveries['y_loc'] = ind_y
    this_step_deliveries['type'] = del_type
    sorted_deliveries = this_step_deliveries.sort_values(by=['y_loc'])
      
    return total_deliveries, sorted_deliveries
  
  def add_sankey_custom_diversion(self, total_width, deliveries, old_inverse_slope, counter):
    start_x = self.sankey_path_x[counter]
    end_x = self.sankey_path_x[counter+1]
    start_y = self.sankey_path_y[counter]
    end_y = self.sankey_path_y[counter+1]
    slope = (end_y - start_y) / (end_x - start_x)
    inverse_slope = -1.0/slope
    if counter == 0:
      old_inverse_slope = inverse_slope * 1.0
      
    if inverse_slope * old_inverse_slope < 0.0:
      box_x1 = start_x - (0.5 - deliveries) * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      box_x2 = start_x - 0.5 * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      box_y1 = start_y - (0.5 - deliveries) * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      box_y2 = start_y - 0.5 * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      x_range = [box_x2, box_x1]
      y_bottom = [box_y2, box_y1]
      y_top = [box_y2, box_y2]
    else:
      box_x1 = start_x + 0.5 * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      box_x2 = start_x + (0.5 - deliveries) * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      box_y1 = start_y + 0.5 * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      box_y2 = start_y + (0.5 - deliveries) * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      x_range = [box_x2, box_x1]
      y_bottom = [box_y2, box_y2]
      y_top = [box_y2, box_y1]
    
    self.ax.fill_between(x_range,y_bottom, y_top, facecolor = 'steelblue', edgecolor = 'steelblue', alpha = 0.9, linewidth = 0.05, zorder = 3)
    
    return [box_x1, box_y1], [box_x1, box_y2]
 
  def add_sankey_custom_diversion2(self, river_diversion_x, river_diversion_y, total_diversion, destination_centroid_x, destination_centroid_y, color_use):
    y = np.linspace(-10, 10, 100)
    y_plot = np.linspace(river_diversion_x, destination_centroid_x, 100)
        
    bottom_sig_height = destination_centroid_y - total_diversion * 0.5 -  min(river_diversion_y[0], river_diversion_y[1])
    top_sig_height = destination_centroid_y + total_diversion * 0.5 -  max(river_diversion_y[0], river_diversion_y[1])

    z_bottom = bottom_sig_height/(1 + np.exp(-1 * y)) + min(river_diversion_y[0], river_diversion_y[1])
    z_top = top_sig_height/(1 + np.exp(-1 * y)) + max(river_diversion_y[0], river_diversion_y[1])
    
    self.ax.fill_between(y_plot, z_bottom, z_top, facecolor=color_use, edgecolor = 'none', alpha = 1.0, zorder = 4)
 
  def add_sankey_custom_path(self, total_width, deliveries, old_inverse_slope, counter):
    start_x = self.sankey_path_x[counter]
    end_x = self.sankey_path_x[counter+1]
    start_y = self.sankey_path_y[counter]
    end_y = self.sankey_path_y[counter+1]
    slope = (end_y - start_y) / (end_x - start_x)
    inverse_slope = -1.0/slope

    if counter == 0:
      old_inverse_slope = inverse_slope * 1.0
      
    if old_inverse_slope < 0.0:
      if inverse_slope < 0.0:
        box_x1 = start_x + (0.5 - deliveries) * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_x2 = start_x - 0.5 * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y1 = start_y + (0.5 - deliveries) * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y2 = start_y - 0.5 * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      else:
        box_x1 = start_x + 0.5 * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_x2 = start_x - (0.5 - deliveries) * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y1 = start_y + 0.5 * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y2 = start_y - (0.5 - deliveries) * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      
    else:
      if inverse_slope > 0.0:
        box_x1 = start_x - 0.5 * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_x2 = start_x + (0.5 - deliveries) * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y1 = start_y - 0.5 * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y2 = start_y + (0.5 - deliveries) * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      else:
        box_x1 = start_x - (0.5 - deliveries) * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_x2 = start_x + 0.5 * total_width / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y1 = start_y - (0.5 - deliveries) * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
        box_y2 = start_y + 0.5 * total_width * old_inverse_slope / np.power(1 + np.power(old_inverse_slope, 2), 0.5)
      

    if slope > 0.0:
      box_x3 = end_x + (0.5 - deliveries/2.0) * total_width / np.power(1 + np.power(inverse_slope, 2), 0.5)
      box_x4 = end_x - (0.5 - deliveries/2.0) * total_width / np.power(1 + np.power(inverse_slope, 2), 0.5)
      box_y3 = end_y + (0.5 - deliveries/2.0) * total_width * inverse_slope / np.power(1 + np.power(inverse_slope, 2), 0.5)
      box_y4 = end_y - (0.5 - deliveries/2.0) * total_width * inverse_slope / np.power(1 + np.power(inverse_slope, 2), 0.5)
    else:
      box_x3 = end_x - (0.5 - deliveries/2.0) * total_width / np.power(1 + np.power(inverse_slope, 2), 0.5)
      box_x4 = end_x + (0.5 - deliveries/2.0) * total_width / np.power(1 + np.power(inverse_slope, 2), 0.5)
      box_y3 = end_y - (0.5 - deliveries/2.0) * total_width * inverse_slope / np.power(1 + np.power(inverse_slope, 2), 0.5)
      box_y4 = end_y + (0.5 - deliveries/2.0) * total_width * inverse_slope / np.power(1 + np.power(inverse_slope, 2), 0.5)
      
    box_length = np.asarray([box_x1, box_x2, box_x3, box_x4])
    sorted_length = np.sort(box_length)
    total_length_1 = sorted_length[1] - sorted_length[0]
    total_length_2 = sorted_length[2] - sorted_length[1]
    total_length_3 = sorted_length[3] - sorted_length[2]
    slope_top = (box_y4 - box_y2) / (box_x4 - box_x2)
    slope_bottom = (box_y3 - box_y1) / (box_x3 - box_x1)
    if box_x2 < box_x3:
      slope_use_2t = slope_top * 1.0
      slope_use_2b = slope_bottom * 1.0
    else:
      slope_use_2t = old_inverse_slope * 1.0
      slope_use_2b = inverse_slope * 1.0
    if old_inverse_slope < 0.0:
      box_y5 = box_y2 + old_inverse_slope * total_length_1
      box_y6 = box_y2 + slope_top * total_length_1
    else:
      box_y5 = box_y1 + slope_bottom * total_length_1 #bottom
      box_y6 = box_y1 + old_inverse_slope * total_length_1 #top
    box_y7 = box_y5 + slope_use_2b * total_length_2#top
    box_y8 = box_y6 + slope_use_2t * total_length_2#bottom
    box_y9 = box_y7 + slope_bottom * total_length_3
    box_y10 = box_y8 + inverse_slope * total_length_3

    if old_inverse_slope < 0.0:
      self.ax.fill_between([sorted_length[0], sorted_length[1]], [box_y2, box_y5], [box_y2, box_y6], facecolor = 'steelblue', edgecolor = 'steelblue', alpha = 0.9, linewidth = 0.3, zorder = 3)
    else:
      self.ax.fill_between([sorted_length[0], sorted_length[1]], [box_y1, box_y5], [box_y1, box_y6], facecolor = 'steelblue', edgecolor = 'steelblue', alpha = 0.9, linewidth = 0.3, zorder = 3)  
    self.ax.fill_between([sorted_length[1], sorted_length[2]], [box_y5, box_y7], [box_y6, box_y8], facecolor = 'steelblue', edgecolor = 'steelblue', alpha = 0.9, linewidth = 0.3, zorder = 3)
    if slope > 0.0:
      self.ax.fill_between([sorted_length[2], sorted_length[3]], [box_y7, box_y3], [box_y8, box_y3], facecolor = 'steelblue', edgecolor = 'steelblue', alpha = 0.9, linewidth = 0.3, zorder = 3)
    else:
      self.ax.fill_between([sorted_length[2], sorted_length[3]], [box_y7, box_y4], [box_y8, box_y4], facecolor = 'steelblue', edgecolor = 'steelblue', alpha = 0.9, linewidth = 0.3, zorder = 3)
  
    return inverse_slope

  def add_sankey_title(self, column, title, height):
    self.ax.text(column, height, title, fontsize = 20, weight = 'bold', fontname = 'Gill Sans MT',verticalalignment='center', horizontalalignment='center') 

  def format_sankey(self, xlim, ylim, fig_name):
    self.ax.set_xlim(xlim)
    self.ax.set_ylim(ylim)
    self.ax.axis('off')    
    plt.tight_layout()
    plt.savefig(fig_name + '.png', dpi = 300, bbox_inches = 'tight', pad_inches = 0.0)
                  
  def add_map_label(self, x_loc, y_loc, label_text, fontsize = 16, fontname = 'Constantia', va = 'center', ha = 'center', zorder = 20):
    
    self.ax.text(x_loc, y_loc, label_text, fontsize = fontsize, weight = 'bold', fontname = fontname,verticalalignment=va, horizontalalignment=ha, zorder = zorder)

  def load_batch_raster(self, project_folder, raster_id_list, raster_name_pt1, raster_name_pt2, raster_band_list, projection_string, max_bright = (100.0, 1500.0), use_gamma = 0.4, contrast_range = (0.5, 99.5), upscale_factor = 'none'):
    for grid_cell in raster_id_list:
      for date_range in raster_id_list[grid_cell]:
        file_folder = project_folder + raster_name_pt1 + grid_cell + '_' + date_range + raster_name_pt2 + '/'
        rgb_list = []
        for band_name in raster_band_list:
          rgb_list.append(raster_name_pt1 + grid_cell + '_' + date_range + raster_name_pt2 + '_' + band_name)

        self.load_sat_bands(file_folder, rgb_list, projection_string, max_bright = max_bright, use_gamma =use_gamma, contrast_range = contrast_range, upscale_factor = upscale_factor)
  
  def load_sat_bands(self, file_folder, rgb_list, projection_string, max_bright = (100.0, 1500.0), use_gamma = 0.4, contrast_range = (0.5, 99.5), upscale_factor = 'none'):
    
    counter = 0
    for raster_name in rgb_list:
      ds = gdal.Open(file_folder + raster_name + '.tif')
      #geotransform - tells you where pixel 0 (x,y) is in the coordinate system and the scale of pixels to coordinates
      geotransform = ds.GetGeoTransform()
      ##clip the raster - note: this is done in UTM coordinates before projection
      ##so that you don't need to re-project entire contiguous US
      clip_name = raster_name + '_clipped.tif'
      #ds = gdal.Translate(clip_name, ds, projWin = [geotransform[0] + x_bound[0]*geotransform[1], geotransform[3] + y_bound[0]*geotransform[5], geotransform[0] + x_bound[1]*geotransform[1], geotransform[3] +y_bound[1]*geotransform[5]])
      output_name = file_folder + raster_name + '_projected.tif'
      ##re-project raster from UTM to LAT/LONG
      gdal.Warp(output_name, ds, dstSRS = projection_string)
      raster = rasterio.open(output_name)
      if upscale_factor == 'none':
        ind_rgb = raster.read()
      else:
        ind_rgb = raster.read(out_shape=(raster.count, int(raster.height * upscale_factor), int(raster.width * upscale_factor)), resampling=rasterio.enums.Resampling.bilinear)
      zero_mask = np.zeros(ind_rgb.shape)
      ones_mask = np.ones(ind_rgb.shape)
#      ind_rgb[ind_rgb <= 1.0] = ones_mask[ind_rgb <= 1.0]
      real_values = ind_rgb[~np.isnan(ind_rgb)]
      pLow = max_bright[0] 
      pHigh = max_bright[1] 
      ind_rgb = exposure.rescale_intensity(ind_rgb, in_range=(pLow,pHigh))
      if self.band_min[counter] == -1.0:
        self.band_min[counter], self.band_max[counter] = np.min(ind_rgb), np.max(ind_rgb)
      ind_rgb = (ind_rgb - self.band_min[counter])/(self.band_max[counter] - self.band_min[counter])
      ind_rgb = exposure.adjust_gamma(ind_rgb, gamma = use_gamma, gain = 1)
      if counter == 0:
        rgb_bands = np.zeros((4, ind_rgb.shape[1], ind_rgb.shape[2]))
      rgb_bands[counter,:,:] = ind_rgb
      
      counter = counter + 1
    
    true_value_mask = np.ones((1,ind_rgb.shape[1], ind_rgb.shape[2]))
    false_value_overlay = np.zeros((1, ind_rgb.shape[1], ind_rgb.shape[2]))
    false_value_mask = ind_rgb == 0.0
    true_value_mask[false_value_mask] = false_value_overlay[false_value_mask]
    rgb_bands[3,:,:] = true_value_mask
    image = reshape_as_image(rgb_bands)
   
    del rgb_bands
    del ind_rgb
    del true_value_mask
    del false_value_overlay
    del false_value_mask

    ##put bands into RGB order
    #image = rasterio.plot.reshape_as_image(b)
    ##give coordinates of image bands
    spatial_extent = rasterio.plot.plotting_extent(raster)
    ##plot image
    if self.type == 'single':
      self.ax.imshow(image, extent = spatial_extent)
    elif self.type == '1d':
      self.ax[nr].imshow(image, extent = spatial_extent)
    elif self.type == '2d':
      for nr in range(0, self.sub_rows):
        for nc in range(0, self.sub_cols):

          self.ax[nr][nc].imshow(image, extent = spatial_extent)


  def add_legend(self, legend_location, legend_handle, legend_properties, nr = 0, nc = 0):
    self.ax.legend(handles=legend_handle, loc=legend_location, prop=legend_properties, frameon = False)

  def add_inset_figure(self, shapefile_name, shape2, shape3, box_lim_x, box_lim_y, epsg_num, use_water = False, inset_lim_x = 0.0, inset_lim_y = 0.0):

    axins = inset_axes(self.ax, width = '30%', height = '30%', loc = 1, bbox_to_anchor=(0,0,1,1), bbox_transform=self.ax.transAxes)
    map_shape = gpd.read_file(shapefile_name)
    map_shape = map_shape.to_crs(epsg = epsg_num)
    map_shape_state = map_shape[map_shape['STATE_ABBR'] == 'AZ']
    map_shape_state.plot(ax = axins, facecolor = 'beige', edgecolor = 'black', linewidth = 1.0, alpha = 0.6, zorder = 2)
    try:
      shape2.plot(ax = axins, color = 'black', linewidth = 3.0, zorder = 5)
    except:
      pass
    try:
      shape3.plot(ax = axins, facecolor = 'slategray', edgecolor = 'black', linewidth = 1.0, markersize = 75.0, zorder = 10)
    except:
      pass

    p2 = Polygon([(box_lim_x[0], box_lim_y[0]), (box_lim_x[1], box_lim_y[0]), (box_lim_x[1], box_lim_y[1]), (box_lim_x[0], box_lim_y[1])])
    df1 = gpd.GeoDataFrame({'geometry': p2, 'df1':[1,1]})
    df1.crs = {'init' :'epsg:'+ str(epsg_num)}
    df1.plot(ax = axins, facecolor = 'none', edgecolor = 'crimson', alpha = 1.0, linewidth = 5.0, zorder = 3)
    if use_water:
      p3 = Polygon([(inset_lim_x[0], inset_lim_y[0]), (inset_lim_x[1], inset_lim_y[0]), (inset_lim_x[1], inset_lim_y[1]), (inset_lim_x[0], inset_lim_y[1])])
      df2 = gpd.GeoDataFrame({'geometry': p3, 'df2':[1,1]})
      df2.plot(ax = axins, facecolor = 'steelblue', edgecolor = 'steelblue', linewidth = 1.0, alpha = 1.0, zorder = 1)
      
    extent = map_shape_state.bounds
    extent = pd.DataFrame(extent)
    extent = extent.reset_index()
    index = 0
    x_extent = extent.loc[index, 'maxx'] - extent.loc[index, 'minx']
    y_extent = extent.loc[index, 'maxy'] - extent.loc[index, 'miny']
    xl = extent.loc[index, 'minx'] - x_extent/8.0
    xr = extent.loc[index, 'maxx'] + x_extent/8.0
    by = extent.loc[index, 'miny'] - y_extent/8.0
    uy = extent.loc[index, 'maxy'] + y_extent/8.0
    axins.axis('off')
    axins.set_xlim(xl, xr)
    axins.set_ylim(by, uy)
    axins.set_xticklabels('')
    axins.set_yticklabels('')      
