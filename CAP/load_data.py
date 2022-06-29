import pandas as pd
import numpy as np
import geopandas as gpd

def load_key_dicts(year_start, year_end):

  key_dict = {}
  key_dict['ASARCO'] = ['ASARCO-HAYDEN OPS', 'ASARCO, LLC']
  key_dict['Avondale'] = ['CITY OF AVONDALE',]
  key_dict['AZ State Land'] = []
  key_dict['AZWC, Coolidge'] = ['AZ WATER CO - COOLIDGE',]
  key_dict['AZWC, White Tank'] = ['AZ WATER CO - WHITE TANK',]
  key_dict['AZWC, Casa Grande'] = ['AZ WATER CO - WHITE TANK',]
  key_dict['AZWC, Superstition'] = ['AZ WATER CO - WHITE TANK',]
  #key_dict['AZWC, Casa Grande'] = ['AZ WATER CO - APACHE JUNCT','AZ WATER CO - BISBEE', 'AZ WATER CO - LAKESIDE','AZ WATER CO - MIAMI/CLAYPOOL','AZ WATER CO - ORACLE','AZ WATER CO - OVERGAARD','AZ WATER CO - PINETOP LAKES', 'AZ WATER CO - PINEWOOD','AZ WATER CO - SAN MANUEL','AZ WATER CO - SEDONA WATER SYSTEM','AZ WATER CO - SIERRA VISTA','AZ WATER CO - STANFIELD','AZ WATER CO - SUPERIOR SYS','AZ WATER CO - TIERRA GRANDE','AZ WATER CO - WINKELMAN']
  #key_dict['AZWC, Superstition'] = ['AZ WATER CO - APACHE JUNCT','AZ WATER CO - BISBEE', 'AZ WATER CO - LAKESIDE','AZ WATER CO - MIAMI/CLAYPOOL','AZ WATER CO - ORACLE','AZ WATER CO - OVERGAARD','AZ WATER CO - PINETOP LAKES', 'AZ WATER CO - PINEWOOD','AZ WATER CO - SAN MANUEL','AZ WATER CO - SEDONA WATER SYSTEM','AZ WATER CO - SIERRA VISTA','AZ WATER CO - STANFIELD','AZ WATER CO - SUPERIOR SYS','AZ WATER CO - TIERRA GRANDE','AZ WATER CO - WINKELMAN']
  key_dict['Buckeye'] = ['CITY OF BUCKEYE - VALENCIA TOWN DIVISION', 'CITY OF BUCKEYE SONORA-SUNDANCE', 'CITY OF BUCKEYE TARTESSO WATER SYSTEM', 'TOWN OF BUCKEYE']
  key_dict['CAGRD'] = []
  key_dict['Carefree WC'] = ['CAREFREE WATER COMPANY',]
  key_dict['Cave Creek'] = ['CAVE CREEK WATER COMPANY',]
  key_dict['Chandler'] = ['CITY OF CHANDLER',]
  key_dict['Chaparral City WC'] = ['CHAPARRAL WATER COMPANY',]
  key_dict['El Mirage'] = ['CITY OF EL MIRAGE',]
  key_dict['Eloy'] = ['CITY OF ELOY']
  key_dict['EPCOR, AF'] = ['EPCOR - AGUA FRIA',]
  key_dict['EPCOR, PV'] = ['EPCOR - PARADISE VALLEY/SCOTTSDALE',]
  key_dict['EPCOR, SC'] = ['EPCOR -SUN CITY',]
  key_dict['EPCOR, SCW'] = ['EPCOR - SUN CITY WEST',]
  key_dict['EPCOR, ST'] = ['EPCOR - SAN TAN ANTHEM',]
  key_dict['Florence'] = []
  key_dict['Freeport'] = ['FREEPORT MCMORAN - TOWN OF BAGDAD']
  key_dict['Freeport-Miami'] = []
  key_dict['FWID'] = ['FLOWING WELLS IRRIGATION DIST',]
  key_dict['Gilbert'] = ['TOWN OF GILBERT']
  key_dict['Glendale'] = ['CITY OF GLENDALE',]
  key_dict['Goodyear'] = ['CITY OF GOODYEAR',]
  key_dict['Greater Tonopah, Water Utility'] = []
  key_dict['Marana'] = ['MARANA MUNICIPAL - HARTMAN VISTAS', 'MARANA MUNICIPAL -AIRLINE LAMBERT', 'MARANA MUNICIPAL- LA PUERTA', 'MARANA MUNICIPAL- MARANA', 'MARANA MUNICIPAL-CORTARO', 'MARANA MUNICIPAL-PALO VERDE', 'MARANA-PICTURE ROCKS']
  key_dict['Maricopa Cty P&R'] = ['MARICOPA CONSOLIDATED DWID', 'MARICOPA MTN DWID 1', 'MARICOPA MTN DWID 2']
  key_dict['Mesa'] = ['CITY OF MESA',]
  key_dict['Metro DWID'] = ['METROPOLITAN DWID', 'METROPOLITAN DWID - DIABLO VILLAGE', 'METROPOLITAN DWID - E&T', 'METROPOLITAN DWID - LAZY B']
  key_dict['Oro Valley'] = ['ORO VALLEY WATER UTILITY', 'ORO VALLEY WATER-COUNTRYSIDE']
  key_dict['Peoria'] = ['CITY OF PEORIA']
  key_dict['Phoenix'] = ['CITY OF PHOENIX', 'CITY OF PHOENIX - ANTHEM']
  key_dict['Queen Creek'] = ['TOWN OF QUEEN CREEK',]
  key_dict['Resolution Copper'] = []
  key_dict['Rio Verde Utilities'] = ['RIO VERDE UTILITIES, INC.',]
  key_dict['Rosemont Copper'] = []
  key_dict['Scottsdale'] = ['CITY OF SCOTTSDALE',]
  key_dict['Spanish Trail WC'] = ['SPANISH TRAIL WATER COMPANY',]
  key_dict['SRP'] = []
  key_dict['Surprise'] = ['SURPRISE, CITY OF - MOUNTAIN VISTA',]
  key_dict['Tempe'] = ['CITY OF TEMPE WATER MGMT DIV',]
  key_dict['Tonto Hills DWID'] = ['TONTO HILLS DWID',]
  key_dict['Tucson'] = ['TUCSON, CITY OF', 'TUCSON WATER SILVERBELL WEST', 'TUCSON WATER-CATALINA', 'TUCSON WATER-CORONA', 'TUCSON WATER-DIAMOND BELL', 'TUCSON WATER-RANCHO DEL SOL LINDO', 'TUCSON WATER-SIERRITA FOOTHILLS', 'TUCSON WATER-SUNSET RANCH', 'TUCSON WATER-THUNDERHEAD RANCH', 'TUCSON WATER-VALLEY VIEW ACRES']
  key_dict['Vail WC'] = ['VAIL WATER COMPANY',]
  key_dict['WUCFD, Apache Junction'] = ['APACHE JUNCTION WATER FACILITIES DISTRICT',]
  key_dict['Coolidge'] = 'duplicate'
  key_dict['Tank'] = 'duplicate'
  key_dict['Superstition'] = 'duplicate'
  key_dict['Carefree'] = 'duplicate'
  key_dict['Cave'] = 'duplicate'
  key_dict['Chaparral'] = 'duplicate'
  key_dict['Mirage'] = 'duplicate'
  key_dict['AF'] = 'duplicate'
  key_dict['PV'] = 'duplicate'
  key_dict['SC'] = 'duplicate'
  key_dict['SCW'] = 'duplicate'
  key_dict['ST'] = 'duplicate'
  key_dict['Tonopah'] = 'duplicate'
  key_dict['Maricopa'] = 'duplicate'
  key_dict['Metro'] = 'duplicate'
  key_dict['Oro'] = 'duplicate'
  key_dict['Queen'] = 'duplicate'
  key_dict['Resolution'] = 'duplicate'
  key_dict['Verde'] = 'duplicate'
  key_dict['Rosemont'] = 'duplicate'
  key_dict['Spanish'] = 'duplicate'
  key_dict['Tonto'] = 'duplicate'
  key_dict['Vail'] = 'duplicate'
  key_dict['WUCFD,'] = 'duplicate'
  key_dict['Morenci,'] = 'duplicate'
  key_dict['Safford,'] = 'duplicate'

  key_dict_irr = {}
  key_dict_irr['BKW Farms'] = []
  key_dict_irr['CAIDD'] = ['Central Arizona Irrigation and Drainage District',]
  key_dict_irr['CHCID'] = ['Chandler Heights Citrus Irrigation District',]
  key_dict_irr['HIDD'] = ['Hohokam Irrigation District',]
  key_dict_irr['HVID'] = ['Harquahala Valley Irrigation District',]
  key_dict_irr['MSIDD'] = ['Maricopa - Stanfield Irrigation and Drainage Dist',]
  key_dict_irr['MWD'] = ['Maricopa Water District',]
  key_dict_irr['NMIDD'] = ['New Magma Irrigation and Drainage District',]
  key_dict_irr['QCID'] = ['Queen Creek Irrigation District',]
  key_dict_irr['SCIDD'] = ['San Carlos Irrigation and Drainage District',]
  key_dict_irr['RWCD'] = ['Roosevelt Water Conservation District',]
  key_dict_irr['SRP'] = []
  key_dict_irr['TID'] = ['Tonopah Irrigation District',]

  key_dict_tribal = {}
  key_dict_tribal['Ak-Chin'] = ['Maricopa (Ak-Chin) Reservation',]
  key_dict_tribal['GRIC'] = ['Gila River Reservation',]
  key_dict_tribal['Pascua Yaqui'] = ['Pascua Yaqui Reservation',]
  key_dict_tribal["Tohono O'odham - ST"] = ['Papago Reservation']
  key_dict_tribal["Tohono O'odham - SX"] = ['San Xavier Reservation',]
  key_dict_tribal['Yaqui'] = 'duplicate'
  key_dict_tribal['ST'] = 'duplicate'
  key_dict_tribal['SX'] = 'duplicate'
  key_dict_tribal['SCAT'] = 'duplicate'
  
  key_dict_recharge = {}
  key_dict_recharge['AWBA Phx AMA'] = ['PHOENIX AMA',]
  key_dict_recharge['AWBA Pinal AMA'] = ['PINAL AMA',]
  key_dict_recharge['AWBA Tucson AMA'] = ['TUCSON AMA',]
  key_dict_recharge['CAGRD'] = ['TUCSON AMA',]
  key_dict_recharge['USBR'] = ['TUCSON AMA',]
  key_dict_recharge['AWBA Interstate'] = ['TUCSON AMA',]
  key_dict_recharge['Phx'] = 'duplicate'
  key_dict_recharge['Pinal'] = 'duplicate'
  key_dict_recharge['Tucson'] = 'duplicate'
  key_dict_recharge['Interstate'] = 'duplicate'
   
  total_deliveries = {}
  for x in key_dict:
    if key_dict[x] == 'duplicate':
      skipthis = 1
    else:
      total_deliveries[x] = np.zeros(12 * (year_end - year_start))
  for x in key_dict_irr:
    if key_dict_irr[x] == 'duplicate':
      skipthis = 1
    else:
      total_deliveries[x] = np.zeros(12 * (year_end - year_start))
  for x in key_dict_tribal:
    if key_dict_tribal[x] == 'duplicate':
      skipthis = 1
    else:
      total_deliveries[x] = np.zeros(12 * (year_end - year_start))
  for x in key_dict_recharge:
    if key_dict_recharge[x] == 'duplicate':
      skipthis = 1
    else:
      total_deliveries[x] = np.zeros(12 * (year_end - year_start))


  return key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge, total_deliveries
  
def set_buffer_vals(year_num):
  buffer_vals = {}
  for group_type in ['Irrigation', 'Excess', 'Tribal', 'Exchange', 'Muni']:
    buffer_vals[group_type] = {}
  if year_num == 2022:
    buffer_vals['Tribal']['Ak-Chin'] = 1.0

  buffer_vals['Muni']['Chandler'] = 1.0
  buffer_vals['Muni']['Phoenix'] = 1.0
  buffer_vals['Muni']['Scottsdale'] = 1.0
  buffer_vals['Muni']['Mesa'] = 1.0
  if year_num == 2018:
    buffer_vals['Muni']['Tucson'] = 1.0
  
  if year_num == 2021:
    buffer_vals['Exchange']['Ak-Chin'] = 1.0
  if year_num == 2020:
    buffer_vals['Exchange']['CAGRD'] = 7.0
  buffer_vals['Exchange']['Chandler'] = 3.0
  if year_num == 2020:
    buffer_vals['Exchange']['Chandler'] = 4.0
  buffer_vals['Exchange']['Morenci'] = 2.0
  if year_num < 2020:
    buffer_vals['Exchange']['Safford'] = 1.0
  buffer_vals['Exchange']['Gilbert'] = 3.0
  buffer_vals['Exchange']['Glendale'] = 2.0
  buffer_vals['Exchange']['Mesa'] = 3.0
  buffer_vals['Exchange']['Yaqui'] = 2.0
  buffer_vals['Exchange']['SCAT'] = 1.0
  if year_num == 2018:
    buffer_vals['Exchange']['Peoria'] = 1.0
    buffer_vals['Exchange']['Tucson'] = 1.0
  buffer_vals['Exchange']['Phoenix'] = 3.0
  buffer_vals['Exchange']['Scottsdale'] = 3.0
  buffer_vals['Exchange']['Tempe'] = 2.0
  buffer_vals['Exchange']['SX'] = 1.0
  buffer_vals['Exchange']['ST'] = 1.0
  buffer_vals['Exchange']['Tonto'] = 1.0
  buffer_vals['Exchange']['WUCFD'] = 1.0
  return buffer_vals
  
def load_cap_deliveries(page, key_dict, tolerance, buffer_vals, find_left_column = False, below_word = 'none'):
  text = page.get_text('words')
  user_names = []
  user_locs_high = []
  user_locs_low = []
  user_values = {}
  id_low_ix = 1
  id_high_ix = 3
  for x in key_dict:
    user_values[x] = []
    user_values[x + '_left'] = []
  left_align = 999.999
  if find_left_column:
    for x in text:
      left_align = min(left_align, x[0])
  high_line = -999
  if below_word != 'none':
    for x in text:
      if x[4] == below_word:
        high_line = x[3]  
    
  for x in text:
    if find_left_column:
      if x[0] <= left_align + tolerance:
        use_line = True
      else:
        use_line = False
    else:
      use_line = True
    if x[3] < high_line:
      use_line = False
    if x[4] in key_dict and use_line:
      if x[4] in user_names:
        counter = 0
        for nm, hb, lb in zip(user_names, user_locs_high, user_locs_low):
          if nm == x[4] and hb > x[id_high_ix]:
            user_locs_high[counter] = x[id_high_ix]
            user_locs_low[counter] = x[id_low_ix]
          counter += 1
      else:
        user_names.append(x[4])
        user_locs_high.append(x[id_high_ix])
        user_locs_low.append(x[id_low_ix])            
      
  for x in text:
    for nm, hb, lb in zip(user_names, user_locs_high, user_locs_low):
      use_val = False
      if nm in buffer_vals:
        buffer_extra = buffer_vals[nm] * (hb - lb)
      else:
        buffer_extra = 0.0
      if x[id_high_ix] <= hb + tolerance + buffer_extra and x[id_low_ix] >= lb - tolerance - buffer_extra:
        use_val = True
        try:
          this_value = int(str(x[4]).replace(',',''))
        except:
          use_val = False
        if use_val:
          if key_dict[nm] == 'duplicate':
            for xxx in key_dict:
              if nm in xxx and nm != xxx:
                user_values[xxx].append(this_value)
                user_values[xxx + '_left'].append(x[0])
                break
          else:
            user_values[nm].append(this_value)
            user_values[nm + '_left'].append(x[0])

  return user_values
  
def find_monthly_deliveries(user_values, nm):
  location_vals = np.asarray(user_values[nm + '_left'])
  sorted_locations = np.argsort(location_vals)
  user_deliveries = np.asarray(user_values[nm])
  sorted_deliveries = user_deliveries[sorted_locations]
  return sorted_deliveries
  
def aggregate_deliveries(total_deliveries, sorted_deliveries, nm, buffer_vals_current, year_num, year_start):
  if nm in buffer_vals_current:
    counter1 = 0
    counter2 = 0
    for xx in range(0,len(sorted_deliveries)):
      total_deliveries[int(counter2 + (year_num - year_start) * 12)] += sorted_deliveries[xx]
      counter1 += 1
      if counter1 == int(len(sorted_deliveries) / 15):
        counter2 += 1
        counter1 = 0
        if counter2 == 12:
          break
  else:
    for xx in range(0, 12):
      total_deliveries[int(xx + (year_num - year_start) * 12)] += sorted_deliveries[xx]
  return total_deliveries

  
def create_district_centroids(district_deliveries, key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge, sankey_index, off_ramp_buffer, epsg_use):
  cws = gpd.read_file('CWS_Service_Area/CWS_Service_Area.shp')
  irrdst = gpd.read_file('Irrigation_District/Irrigation_District.shp')
  tribal = gpd.read_file('ALRIS_amerind/ALRIS_amerind.shp')
  canals = gpd.read_file('Cen_AZ_Proj/Cen_AZ_Proj.shp')
  amas = gpd.read_file('AMA_and_INA/AMA_and_INA.shp')
  gwstorage = gpd.read_file('Underground_Storage_Facilities/Underground_Storage_Facilities.shp')
  cws = cws[pd.notna(cws['geometry'])]
  cws = cws.to_crs(epsg = epsg_use)
  irrdst = irrdst.to_crs(epsg = epsg_use)
  tribal = tribal.to_crs(epsg = epsg_use)
  canals = canals.to_crs(epsg = epsg_use)
  amas = amas.to_crs(epsg = epsg_use)
  gwstorage = gwstorage.to_crs(epsg = epsg_use)
  basin_facilities = gpd.sjoin(gwstorage, amas, how = 'inner', op = 'within')
  average_deliveries = []
  centroid_points_x = []
  centroid_points_y = []
  index_locs = []
  district_name_list = []
  district_type = []
  for shapefile_data, key_dict_use, shapefile_name in zip([cws, irrdst, tribal, basin_facilities], [key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge], ['Municipal', 'Irrigation', 'Tribal', 'Recharge']):
    for district_name in key_dict_use:
      subset_list = key_dict_use[district_name]
      centroid_x = 0.0
      centroid_y = 0.0
      subset_counter = 0.0
      for ind_name in subset_list:
        if shapefile_name == 'Municipal':
          this_district = shapefile_data[shapefile_data['CWS_NAME'] == ind_name]          
        elif shapefile_name == 'Irrigation':
          this_district = shapefile_data[shapefile_data['LONG_NAME'] == ind_name]          
        elif shapefile_name == 'Tribal':
          this_district = shapefile_data[shapefile_data['NAME'] == ind_name]          
        elif shapefile_name == 'Recharge':
          this_district = shapefile_data[shapefile_data['BASIN_NAME'] == ind_name]      
        if len(this_district) > 0:
          if shapefile_name == 'Recharge':
            max_area = 0.0
            max_area_index = 0
            for index, row in this_district.iterrows():
              if row['Shape__A_1'] > max_area:
                max_area = row['Shape__A_1'] * 1.0
                max_area_index = index * 1
            cent_cords = this_district.loc[max_area_index, 'geometry'].centroid.coords
            cent_tuple = cent_cords[0]
            centroid_x += cent_tuple[0]
            centroid_y += cent_tuple[1]
            subset_counter += 1.0                
            
          else:            
            for index, row in this_district.iterrows():
              cent_cords = row.geometry.centroid.coords
              cent_tuple = cent_cords[0]
              centroid_x += cent_tuple[0]
              if district_name == 'Phoenix':
                centroid_y += cent_tuple[1] - 13750.0
              else:
                centroid_y += cent_tuple[1]
              subset_counter += 1.0                
      if subset_counter > 0:
        centroid_x = float(centroid_x / subset_counter)
        centroid_y = float(centroid_y / subset_counter)
        for index_loc in range(0, len(sankey_index)):
          if sankey_index[index_loc] > centroid_x - off_ramp_buffer:
            break
        total_deliveries = 12.0 * np.sum(district_deliveries[district_name]) / len(district_deliveries[district_name])
        average_deliveries.append(total_deliveries)
        district_name_list.append(district_name)
        centroid_points_x.append(centroid_x)
        centroid_points_y.append(centroid_y)
        index_locs.append(index_loc)
        district_type.append(shapefile_name)
  
  district_centroids = pd.DataFrame()
  district_centroids['name'] = district_name_list
  district_centroids['x_loc'] = centroid_points_x
  district_centroids['y_loc'] = centroid_points_y
  district_centroids['sankey_loc'] = index_locs
  district_centroids['type'] = district_type
  district_centroids['delivery'] = average_deliveries
  
  return district_centroids
  
  
