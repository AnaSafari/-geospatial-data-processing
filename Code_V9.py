

############################################################################
####### practical experience in some of the low level filtering, ###########
####### aggregation and processing of geospatial data, and       ###########
####### introduction to the principle of operationalization      ###########
############################################################################

############# Geolocation Filtering
from pyproj import Proj, transform
import mysql.connector
import numpy as np
import matplotlib.pyplot as plt

conn = mysql.connector.connect(host="crepe.usask.ca", database='SHED7', user="ans911", password="Change_this_396")
cursor = conn.cursor()  
query = """\
SELECT COUNT(*) FROM gps WHERE (lat >= 52.058367 AND lat <= 52.214608) 
AND (lon >= -106.7649138128 AND lon <= -106.52225318)"""
cursor.execute(query)
recordsInsideSaskatoon=cursor.fetchall()


query = """\
SELECT COUNT(*) FROM gps WHERE (lat >= 52.058367 AND lat <= 52.214608) 
AND (lon >= -106.7649138128 AND lon <= -106.52225318)"""
cursor.execute(query)
totalGPS_records=cursor.fetchall()

eliminatedGPS_records = totalGPS_records - recordsInsideSaskatoon

############# battery Filtering 

query = """\
SELECT COUNT(distinct user_id) from ((SELECT user_id, COUNT(user_id) As countRecordPerUser 
FROM battery group by user_id) AS T) where countRecordPerUser<80
"""
cursor.execute(query)
numberOfRecordsPerPersonLessThanThreshold=cursor.fetchall()

query = """\
SELECT COUNT(distinct user_id) from ((SELECT user_id, COUNT(user_id) As countRecordPerUser 
FROM battery group by user_id) AS T) where countRecordPerUser<50
"""
cursor.execute(query)
numberOfRecordsPerPersonLessThaThreshold=cursor.fetchall()

######## I also used tableau for this part and the sql code is demonstrated in the report 


############################################################################
############################### ANSWER TO Q.2 ##############################
#############################################################################
from pyproj import Proj, transform
import mysql.connector
import numpy as np
import matplotlib.pyplot as plt

conn = mysql.connector.connect(host="crepe.usask.ca", database='SHED7', user="ans911", password="Change_this_396")
cursor = conn.cursor()  
query = """\
SELECT lat, lon FROM gps WHERE (lat >= 52.058367 AND lat <= 52.214608) 
AND (lon >= -106.7649138128 AND lon <= -106.52225318)"""
cursor.execute(query)
numberOfPointInsideGeoBox=list(cursor.fetchall())

lat = [ seq[0] for seq in numberOfPointInsideGeoBox ]
lon = [ seq[1] for seq in numberOfPointInsideGeoBox ]

# LatLon with WGS84 datum used by GPS units and Google Earth
WGS84 = Proj(init='EPSG:4326') 
# UTM zone 13N
UTM13N = Proj(init='EPSG:32613')

# Projects from latitude, longitude to EPSG:32613
Easting, Northing = transform(WGS84, UTM13N, lon, lat)

# To create a 2D histogram of projected points and then to plot a heatmap 
heatmap, xedges, yedges = np.histogram2d(Easting, Northing, bins=(70, 70))
extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
fig, ax = plt.subplots()
cax = ax.imshow(heatmap.T, interpolation='nearest', origin='low', extent= extent)
ax.set_title('Heat Map of projected points to EPSG:32613')
plt.xlabel('Easting')
plt.ylabel('Northing')
cbar = plt.colorbar(cax, orientation='vertical')
cbar.ax.set_xticklabels(['Low', 'Medium', 'High'])
plt.show()

#############################################################################
############################################################################
############################### ANSWER TO Q.3 ##############################
#############################################################################
# Saving mysql query result in Panda Dataframe, adding Northing, 
# Easting and duty_cycle columns to it 
# and saving the dataframe in .csv file

import mysql.connector
import pandas as pd
from pyproj import Proj, transform


conn = mysql.connector.connect(host="crepe.usask.ca", database='SHED7', user="ans911", password="Change_this_396")
query = """\
SELECT user_id, record_time, lat, lon FROM gps WHERE (lat >= 52.058367 AND lat <= 52.214608) 
AND (lon >= -106.7649138128 AND lon <= -106.52225318) AND user_id IN (511,764,807)"""
df = pd.read_sql(query, con=conn) 

df_backuped = df.copy()

# to project gps point to EPSG:32613 
lat = df.lat.tolist()
lon = df.lon.tolist()
WGS84 = Proj(init='EPSG:4326') 
UTM13N = Proj(init='EPSG:32613')
easting, northing = transform(WGS84, UTM13N, lon, lat)


# to add the easting and northing columns to out dataframe
se = pd.Series(easting)
df['easting'] = se.values

se = pd.Series(northing)
df['northing'] = se.values

#################################################################################
########################### Finding Visit Frequency #############################
#################################################################################

from pyproj import Proj, transform
import mysql.connector
import numpy as np
import matplotlib.pyplot as plt

conn = mysql.connector.connect(host="crepe.usask.ca", database='SHED7', user="ans911", password="Change_this_396")
cursor = conn.cursor()  
query = """\
SELECT lat, lon FROM gps WHERE (lat >= 52.058367 AND lat <= 52.214608) 
AND (lon >= -106.7649138128 AND lon <= -106.52225318) AND user_id IN (511,764,807)"""
cursor.execute(query)
numberOfPointInsideGeoBox=list(cursor.fetchall())

lat = [ seq[0] for seq in numberOfPointInsideGeoBox ]
lon = [ seq[1] for seq in numberOfPointInsideGeoBox ]

# LatLon with WGS84 datum used by GPS units and Google Earth
WGS84 = Proj(init='EPSG:4326') 
# UTM zone 13N
UTM13N = Proj(init='EPSG:32613')

# Projects from latitude, longitude to EPSG:32613
Easting, Northing = transform(WGS84, UTM13N, lon, lat)

# To create a 2D histogram of projected points and then to plot a heatmap 
visitFrequency100G, xedges_100G, yedges_100G = np.histogram2d(Easting, Northing, bins=(140, 140))
visitFrequency400G, xedges_400G, yedges_400G = np.histogram2d(Easting, Northing, bins=(35, 35))
visitFrequency1600G, xedges_1600G, yedges_1600G = np.histogram2d(Easting, Northing, bins=(8, 8))

############# Transfering Array of Tuples to List 
visitFrequency100G_inList = visitFrequency100G.tolist()
visitFrequency400G_inList = visitFrequency400G.tolist()
visitFrequency1600G_inList = visitFrequency1600G.tolist()

flattened_visitFrequency100G = [y for x in visitFrequency100G for y in x]
flattened_visitFrequency400G = [y for x in visitFrequency400G for y in x]
flattened_visitFrequency1600G = [y for x in visitFrequency1600G for y in x]
###### Removing zero frequencies
flattened_visitFrequency100G_filtered = [x for x in flattened_visitFrequency100G if x != 0]
flattened_visitFrequency400G_filtered = [x for x in flattened_visitFrequency400G if x != 0]
flattened_visitFrequency1600G_filtered = [x for x in flattened_visitFrequency1600G if x != 0]

import matplotlib.pyplot as plt
 ############## plotting visit of frequency
plt.hist(flattened_visitFrequency100G_filtered, bins = 5000, histtype = 'step', stacked = True, fill = False ,alpha=0.5, label='Grid Size = 100')
plt.hist(flattened_visitFrequency400G_filtered, bins = 5000, histtype = 'step', stacked = True, fill = False ,alpha=0.5, label='Grid Size = 400')
plt.hist(flattened_visitFrequency1600G_filtered, bins = 5000, histtype = 'step', stacked = True, fill = False ,alpha=0.5, label='Grid Size = 1600')
plt.legend(loc='upper right')
plt.title('Frequency of visits for three random participants')
plt.xlabel('Vist counts')
plt.xscale('log')
plt.yscale('log')
plt.ylabel('Count')
plt.show()


#################################################################################
########################### Finding Dwell time and trip length ##################
#################################################################################


import datetime
import time
import numpy as np
import math

# 2016,07,07,12,6,1 is the minimum record date which exists in our gps query result and is set to be 
# the starting time of the study corresponding to duty_cycle = 0
duty_cycle = np.round(pd.Series(df.record_time-datetime.datetime(2016,7,7,12,6,1)).dt.total_seconds()/120)

# to convert to int from float
duty_cycle = duty_cycle.astype(int)

# to add the duty_cycle column to out dataframe
se = pd.Series(duty_cycle)
df['duty_cycle'] = se.values

# a trick to delete duplicate records for each user which fall in one duty-cycle
tuples = list(zip(df.user_id,df.duty_cycle))
se = pd.Series(tuples)
df['tuples'] = se.values
#df.drop_duplicates(['tuples'], keep='last')
df = df.drop_duplicates(['tuples'], keep='last')
del df['tuples']

df_befor_griding = df.copy()



# to find the index of each (100 by 100) bin that our data falls in
# the bins can be set to larger one by changing the values of delta_x and delta_y 
###########################################################################
##### 100
df = df_befor_griding.copy()
delta_x_100G = 100
delta_y_100G = 100
bin_i_100G = (np.floor((df.easting - np.min(df.easting))/delta_x_100G)).astype(int)
bin_j_100G = (np.floor((df.northing - np.min(df.northing))/delta_y_100G)).astype(int)

########## creating a column including a tuple of user and bin
df['user_bin_100G'] = pd.Series(list(zip(df.user_id,bin_i_100G, bin_j_100G))).values

tuples = list(zip(df.user_bin_100G,df.duty_cycle))
se = pd.Series(tuples)
df['tuples2'] = se.values
df = df.drop_duplicates(['tuples2'], keep = 'last')
del df['tuples2']


df = df.sort_values(['user_id','duty_cycle'], ascending=[True, True])

dwell_time=[]
dwell_bin=[]
distance=[]
dc_start = df.iloc[0]['duty_cycle']
bin_start = df.iloc[0]['user_bin_100G']
trip_bin_start = bin_start
count = 0
for i in df.duty_cycle:
	if (i) ==  df.shape[0]:
		break
	elif (df.iloc[i+1]['user_bin_100G'] == df.iloc[i]['user_bin_100G'] and df.iloc[i+1]['duty_cycle']- df.iloc[i]['duty_cycle'] <= 3):
		trip_bin_finish = bin_start
		count += 1
		dist = math.sqrt(pow((trip_bin_finish[1]-trip_bin_start[1]),2) + pow((trip_bin_finish[2]-trip_bin_start[2]),2))
		distance.append(dist)
	else:
		trip_bin_start = bin_start
		dwell_time.append(count)
		dwell_bin.append(bin_start)
		dc_start = df.iloc[i+1]['duty_cycle']
		bin_start = df.iloc[i+1]['user_bin_100G']
		count = 0

columns = [('DwellingBin', dwell_bin), ('DwellingTime', dwell_time)]
dwell_table = pd.DataFrame.from_items(columns)
DwellingDF_100G = dwell_table.loc[dwell_table['DwellingTime'] >= 3]
Distance100G = distance

###########################################################################
##### Grid size 400
df = df_befor_griding.copy()
delta_x_400G = 400
delta_y_400G = 400
bin_i_400G = (np.floor((df.easting - np.min(df.easting))/delta_x_400G)).astype(int)
bin_j_400G = (np.floor((df.northing - np.min(df.northing))/delta_y_400G)).astype(int)


df['user_bin_400G'] = pd.Series(list(zip(df.user_id,bin_i_400G, bin_j_400G))).values

tuples = list(zip(df.user_bin_400G,df.duty_cycle))
se = pd.Series(tuples)
df['tuples2'] = se.values
df = df.drop_duplicates(['tuples2'], keep = 'last')
del df['tuples2']

df = df.sort_values(['user_id','duty_cycle'], ascending=[True, True])

dwell_time=[]
dwell_bin=[]
distance=[]
dc_start = df.iloc[0]['duty_cycle']
bin_start = df.iloc[0]['user_bin_400G']
trip_bin_start = bin_start
count = 0
for i in df.duty_cycle:
	if (i) ==  df.shape[0]:
		break
	elif (df.iloc[i+1]['user_bin_400G'] == df.iloc[i]['user_bin_400G'] and df.iloc[i+1]['duty_cycle']- df.iloc[i]['duty_cycle'] <= 3):
		trip_bin_finish = bin_start
		count += 1
		dist = math.sqrt(pow((trip_bin_finish[1]-trip_bin_start[1]),2) + pow((trip_bin_finish[2]-trip_bin_start[2]),2))
		distance.append(dist)
	else:
		trip_bin_start = bin_start
		dwell_time.append(count)
		dwell_bin.append(bin_start)
		dc_start = df.iloc[i+1]['duty_cycle']
		bin_start = df.iloc[i+1]['user_bin_400G']
		count = 0

columns = [('DwellingBin', dwell_bin), ('DwellingTime', dwell_time)]
dwell_table = pd.DataFrame.from_items(columns)
DwellingDF_400G = dwell_table.loc[dwell_table['DwellingTime'] >= 3]
Distance400G = distance

#########################################################################
######## Grid size 1600

df = df_befor_griding.copy()
delta_x_1600G = 1600
delta_y_1600G = 1600
bin_i_1600G = (np.floor((df.easting - np.min(df.easting))/delta_x_1600G)).astype(int)
bin_j_1600G = (np.floor((df.northing - np.min(df.northing))/delta_y_1600G)).astype(int)


df['user_bin_1600G'] = pd.Series(list(zip(df.user_id,bin_i_1600G, bin_j_1600G))).values

tuples = list(zip(df.user_bin_1600G,df.duty_cycle))
se = pd.Series(tuples)
df['tuples2'] = se.values
df = df.drop_duplicates(['tuples2'], keep = 'last')
del df['tuples2']


df = df.sort_values(['user_id','duty_cycle'], ascending=[True, True])

dwell_time=[]
dwell_bin=[]
distance=[]
dc_start = df.iloc[0]['duty_cycle']
bin_start = df.iloc[0]['user_bin_1600G']
trip_bin_start = bin_start
count = 0
for i in df.duty_cycle:
	if (i) ==  df.shape[0]:
		break
	elif (df.iloc[i+1]['user_bin_1600G'] == df.iloc[i]['user_bin_1600G'] and df.iloc[i+1]['duty_cycle']- df.iloc[i]['duty_cycle'] <= 3):
		trip_bin_finish = bin_start
		count += 1
		dist = math.sqrt(pow((trip_bin_finish[1]-trip_bin_start[1]),2) + pow((trip_bin_finish[2]-trip_bin_start[2]),2))
		distance.append(dist)
	else:
		trip_bin_start = bin_start
		dwell_time.append(count)
		dwell_bin.append(bin_start)
		dc_start = df.iloc[i+1]['duty_cycle']
		bin_start = df.iloc[i+1]['user_bin_1600G']
		count = 0

columns = [('DwellingBin', dwell_bin), ('DwellingTime', dwell_time)]
dwell_table = pd.DataFrame.from_items(columns)
DwellingDF_1600G = dwell_table.loc[dwell_table['DwellingTime'] >= 3]
Distance1600G = distance

############################# plots


import matplotlib.pyplot as plt


plt.hist(DwellingDF_100G.DwellingTime, bins = 100, histtype='step', stacked=True, fill=False, alpha=0.5, label='Grid Size = 100')
plt.hist(DwellingDF_400G.DwellingTime, bins = 100, histtype='step', stacked=True, fill=False, alpha=0.5, label='Grid Size = 400')
plt.hist(DwellingDF_1600G.DwellingTime, bins = 100, histtype='step', stacked=True, fill=False, alpha=0.5, label='Grid Size = 1600')
plt.legend(loc='upper right')
plt.title('Dwell Time Histogram')
plt.xlabel('Dwell Time bin (in duty_cycle)')
plt.ylabel('Count')
plt.xscale('log')
plt.show()

Distance100G = [x for x in Distance100G if x!=0]
Distance400G = [x for x in Distance400G if x!=0]
Distance1600G = [x for x in Distance1600G if x!=0]
import matplotlib.pyplot as plt
plt.hist(Distance100G, bins = 100, histtype='step', stacked=True, fill=False, alpha=0.5, label='Grid Size = 100')
plt.hist(Distance400G, bins = 100, histtype='step', stacked=True, fill=False, alpha=0.5, label='Grid Size = 400')
plt.hist(Distance1600G, bins = 100, histtype='step', stacked=True, fill=False, alpha=0.5, label='Grid Size = 1600')
plt.legend(loc='upper right')
plt.title('Trip Length Histogram')
plt.xlabel('Trip Length')
plt.ylabel('Count')
plt.xscale('log')
plt.show()


#########################################################################################################
###################### Useful Websites ##################################################################
# http://www.uvm.edu/~jbagrow/dsv/heatmap_basemap.html
# http://all-geo.org/volcan01010/2012/11/change-coordinates-with-pyproj/
# https://peak5390.wordpress.com/2012/12/08/matplotlib-basemap-tutorial-installing-matplotlib-and-basemap/
# https://gist.github.com/teechap/9c066a9ab054cc322877
# https://ocefpaf.github.io/python4oceanographers/blog/2013/12/16/utm/
#http://stackoverflow.com/questions/2369492/generate-a-heatmap-in-matplotlib-using-a-scatter-data-set

#http://stackoverflow.com/questions/9706041/finding-index-of-an-item-closest-to-the-value-in-a-list-thats-not-entirely-sort
###########################################################################################################
