# Sankey Water Maps
This repository contains all code and data for Sankey Water Maps, a mapping tool that uses publicly available data to map the distribution of state/federal water projects  
### Existing Maps with Examples:  
Central Arizona Project (2018 - current)  
### Future Maps (no examples yet):  
Central Valley Project: Tehama-Colusa Canal  
Central Valley Project: Friant-Kern Canal  
Sacramento-San Joaquin Delta: CVP/SWP combined operations  

## Installation and setup  
1. Clone this repository to your local machine.  
2. Ensure that the following python packages are installed via pip or conda:  
    a. datetime  
	b. numpy  
	c. pandas  
	d. geopandas  
	e. matplotlib  
	e. requests  
	f. fitz (install with pip install pymupdf)  
	g. shapely  
	h. apng  
	i. scipy  
	j. rasterio  
	k. skimage  
	l. mpl_toolkits  
	
3. Example maps (in SankeyWaterMaps/CAP/sankey_figs/example_(year_num) are shown with background satellite images but image data is not stored on Github  
	a. if you don't want to make maps with satellite images, skip this step  
	b. Map background images are available to be downloaded here: https://drive.google.com/drive/folders/1U8uOTz3eLYN2WZXkn0SpghhIL8qhCtgX?usp=sharing  
	c. place folder az_satellite_photos (from google drive link) in SankeyWaterMaps/CAP working directory  
	d. un-comment line 115 in SankeyWaterMaps/CAP/make_cap_delivery_figure.py so it reads: sankey_plot.load_batch_raster('az_satellite_photos/', raster_id_list, raster_name_pt1, raster_name_pt2, raster_band_list, projection_string, max_bright = (6000.0, 15000.0), use_gamma = 0.8, upscale_factor = 0.1)  
4. To read delivery data from CAP, navigate to the working directory SankeyWaterMaps/CAP and run: python -W ignore get_cap_delivery_data.py  
	a. this will create a file 'deliveries_by_district.csv', which is monthly delivery data to all CAP contractors, Jan 2018 - May 2022  
	b. the default data is created from .pdf delivery data files in the SankeyWaterMaps/CAP folder  
	c. if you want to update data with new delivery files from the CAP website (updated monthly by CAP), navigate to SankeyWaterMaps/CAP an run: python -W ignore get_cap_delivery_data.py 1  
5. After creating the file 'deliveries_by_district.csv' (either from stored or online data), stay in the working directory SankeyWaterMaps/CAP and run: python -W ignore make_cap_delivery_figure.py  
    a. The script will default to creating maps of annual deliveries (2018 - 2022).  
    b. If you want to develop maps of monthly deliveries, run: python -W ignore make_cap_delivery_figure.py 1  
	c. If you want to make an animated PNG (like a .gif) map of monthly deliveries, run: python -W ignore make_cap_delivery_figure.py 1 1  
	
6. Maps will write to SankeyWaterMaps/CAP/sankey_figs  
	a. Animated PNGs can be viewed with Mozilla Firefox or Google Chrome  