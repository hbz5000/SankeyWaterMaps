import sys
import pandas as pd
import numpy as np
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from datetime import datetime, timedelta
import requests
import fitz
import load_data as ld

##############################################
###FILE READ PARAMETERS
##############################################
if len(sys.argv) > 1:
  try:
    update_data = bool(int(sys.argv[1]))
  except:
    update_data = False
else:
  update_data = False
year_start = 2018
year_end = 2023
num_years = year_end - year_start
tolerance = 0.1#for reading PDF tables - what is the tolerance for 'same row'
#load district names
key_dict, key_dict_irr, key_dict_tribal, key_dict_recharge, total_deliveries = ld.load_key_dicts(year_start, year_end)

#annual loop - one year per pdf
for year_num in range(year_start, year_end):
  #some of the tables have multi-row entries, this has to be added manually
  buffer_vals = ld.set_buffer_vals(year_num)
  year_string = str(year_num)
  
  ##############################################
  ###DOWNLOAD DATA
  ##############################################
  #download delivery data pdfs from the cap website (this saves the pdf)
  #cap data is updated monthly, this doesn't need to be updated every time
  temp_filename = 'cap_deliveries_' + year_string + '.pdf'
  if update_data:
    print('read from web')
    url1 = 'https://library.cap-az.com/documents/departments/water-operations/year-to-date-by-contract-type-' + year_string + '.pdf'
    response = requests.get(url1)      
    my_raw_data = response.content
    with open(temp_filename, 'wb') as my_data:
      my_data.write(my_raw_data)
  ##############################################
  ###READ PDF DATA
  ##############################################
  #read pdf
  with fitz.open(temp_filename) as doc:
    counter = 0
    for page in doc:
      #first page of document has 'excess' deliveries
      #they go to irrigation & recharge customers
      if counter == 1:
        #irrigation customers
        buffer_vals_current = buffer_vals['Irrigation']
        user_values = ld.load_cap_deliveries(page, key_dict_irr, tolerance, buffer_vals_current)
        for nm in key_dict_irr:
          if len(user_values[nm]) > 0:
            #assign one or more row to the appropriate month
            sorted_deliveries = ld.find_monthly_deliveries(user_values, nm)
            total_deliveries[nm] = ld.aggregate_deliveries(total_deliveries[nm], sorted_deliveries, nm, buffer_vals_current, year_num, year_start)                        
        #recharge customers
        buffer_vals_current = buffer_vals['Excess']
        user_values = ld.load_cap_deliveries(page, key_dict_recharge, tolerance, buffer_vals_current)
        for nm in key_dict_recharge:
          if len(user_values[nm]) > 0:
            #assign one or more row to the appropriate month
            sorted_deliveries = ld.find_monthly_deliveries(user_values, nm)
            total_deliveries[nm] = ld.aggregate_deliveries(total_deliveries[nm], sorted_deliveries, nm, buffer_vals_current, year_num, year_start)            
      #second page of document has 'tribal' deliveries
      if counter == 2:
        #tribal customers
        buffer_vals_current = buffer_vals['Tribal']
        user_values = ld.load_cap_deliveries(page, key_dict_tribal, tolerance, buffer_vals_current)
        for nm in key_dict_tribal:
          if len(user_values[nm]) > 0:
            sorted_deliveries = ld.find_monthly_deliveries(user_values, nm)
            total_deliveries[nm] = ld.aggregate_deliveries(total_deliveries[nm], sorted_deliveries, nm, buffer_vals_current, year_num, year_start)            
      #second through fourth page of document has 'exchange' deliveries
      #exchanges go to both municipal and tribal users
      if counter >= 2 and counter <= 4:
        buffer_vals_current = buffer_vals['Exchange']
        #exchange deliveries only reads the name on the farthest left column
        user_values_a = ld.load_cap_deliveries(page, key_dict, tolerance, buffer_vals_current, find_left_column = True, below_word = 'Off-Res')
        user_values_b = ld.load_cap_deliveries(page, key_dict_tribal, tolerance, buffer_vals_current, find_left_column = True, below_word = 'Off-Res')
        #municipal users
        for nm in key_dict:
          if len(user_values_a[nm]) > 0:
            sorted_deliveries_a = ld.find_monthly_deliveries(user_values_a, nm)
            total_deliveries[nm] = ld.aggregate_deliveries(total_deliveries[nm], sorted_deliveries_a, nm, buffer_vals_current, year_num, year_start)    
        #tribal users            
        for nm in key_dict_tribal:
          if len(user_values_b[nm]) > 0:
            sorted_deliveries_b = ld.find_monthly_deliveries(user_values_b, nm)
            total_deliveries[nm] = ld.aggregate_deliveries(total_deliveries[nm], sorted_deliveries_b, nm, buffer_vals_current, year_num, year_start)          
      #page 5+ has municipal deliveries
      if counter >= 5:
        buffer_vals_current = buffer_vals['Muni']
        user_values = ld.load_cap_deliveries(page, key_dict, tolerance, buffer_vals_current)
        for nm in key_dict:
          if len(user_values[nm]) > 0:
            sorted_deliveries = ld.find_monthly_deliveries(user_values, nm)
            total_deliveries[nm] = ld.aggregate_deliveries(total_deliveries[nm], sorted_deliveries, nm, buffer_vals_current, year_num, year_start)                
      counter += 1    

##############################################
##WRITE DATA TO CSV
##############################################
datetime_index = []
year_num = 2018
month_num = 1
max_length = 0
for x in total_deliveries:
  max_length = max(max_length, len(total_deliveries[x]))
for x in range(0,max_length):
  datetime_index.append(datetime(year_num, month_num, 1, 0, 0))
  month_num += 1
  if month_num == 13:
    month_num = 1
    year_num += 1
total_deliveries = pd.DataFrame(total_deliveries, index = datetime_index) 
total_deliveries.index = datetime_index
total_deliveries.to_csv('deliveries_by_district.csv')

##############################################
##PLOT DATA
##############################################
delivery_check = True
if delivery_check:
  fig, ax = plt.subplots(5)  
  counter = 0    
  running_total = np.zeros(len(total_deliveries['Phoenix']))
  for xx in key_dict_irr:
    if key_dict_irr[xx] != 'duplicate':
      if counter == 0:
        baseline = np.zeros(len(total_deliveries[xx]))
        counter += 1
      ax[0].fill_between(total_deliveries.index, baseline, baseline + total_deliveries[xx])
      ax[4].fill_between(total_deliveries.index, baseline + running_total, baseline + total_deliveries[xx] + running_total, facecolor = 'forestgreen')
      baseline += total_deliveries[xx]
  counter = 0    
  running_total += baseline
  for xx in key_dict_tribal:
    if key_dict_tribal[xx] != 'duplicate':
      if counter == 0:
        baseline = np.zeros(len(total_deliveries[xx]))
        counter += 1
      ax[1].fill_between(total_deliveries.index, baseline, baseline + total_deliveries[xx])
      ax[4].fill_between(total_deliveries.index, baseline + running_total, baseline + total_deliveries[xx] + running_total, facecolor = 'indianred')
      baseline += total_deliveries[xx]
  counter = 0    
  running_total += baseline
  for xx in key_dict:
    if key_dict[xx] != 'duplicate' and xx != 'CAGRD':
      if counter == 0:
        baseline = np.zeros(len(total_deliveries[xx]))
        counter += 1
      ax[2].fill_between(total_deliveries.index, baseline, baseline + total_deliveries[xx])
      ax[4].fill_between(total_deliveries.index, baseline + running_total, baseline + total_deliveries[xx] + running_total, facecolor = 'goldenrod')
      baseline += total_deliveries[xx]
  counter = 0    
  running_total += baseline
  for xx in key_dict_recharge:
    if key_dict_recharge[xx] != 'duplicate':
      if counter == 0:
        baseline = np.zeros(len(total_deliveries[xx]))
        counter += 1
      ax[3].fill_between(total_deliveries.index, baseline, baseline + total_deliveries[xx])
      ax[4].fill_between(total_deliveries.index, baseline + running_total, baseline + total_deliveries[xx] + running_total, facecolor = 'steelblue')
      baseline += total_deliveries[xx]
  ax[0].set_ylabel('Irrigation')
  ax[1].set_ylabel('Tribal')
  ax[2].set_ylabel('Municipal')
  ax[3].set_ylabel('Recharge')
  ax[4].set_ylabel('Total')
  for axis_no in range(0, 5):   
    ax[axis_no].set_xlim([total_deliveries.index[0], datetime(2022, 5, 1, 0, 0)])
  plt.show()

