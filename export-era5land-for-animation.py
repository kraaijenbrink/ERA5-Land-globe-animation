# Quick export of ERA5-Land data from GEE for animation
# Philip Kraaijenbrink, 20201217

import ee
from datetime import datetime as dt
from datetime import timedelta
ee.Initialize()

# get era5 data
e5l       = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY")
expcrs    = e5l.first().projection().crs().getInfo()
expscale  = e5l.first().projection().nominalScale().multiply(2).getInfo()
expextent = ee.Geometry.LinearRing([[-180, -90], [180, -90], [180, 90], [-180, 90], [-180, -90]], 'EPSG:4326', False)

# set list of variables to export
allbands     = e5l.first().bandNames().getInfo()
export_bands = ['dewpoint_temperature_2m','u_component_of_wind_10m', 'total_evaporation_hourly', 'soil_temperature_level_1', 'total_precipitation_hourly', 
                'soil_temperature_level_4', 'volumetric_soil_water_layer_1', 'skin_temperature', 'snowfall_hourly', 'surface_sensible_heat_flux_hourly','v_component_of_wind_10m',
                'surface_solar_radiation_downwards_hourly', 'leaf_area_index_low_vegetation', 'surface_latent_heat_flux_hourly', 'temperature_2m', 'potential_evaporation_hourly', 
                'surface_pressure','surface_net_solar_radiation_hourly']

delta = 0
days_per_var = 1
startdate = dt(2020,1,1)
for band in export_bands:

    # get time filter params
    sdate = startdate + timedelta(days=delta)
    edate = startdate + timedelta(days=delta+days_per_var)
    sdate_str = dt.strftime(sdate,'%Y-%m-%d')
    edate_str = dt.strftime(edate,'%Y-%m-%d')
    delta += days_per_var
    
    # get bandstack of the hourly data for current var
    expimg = e5l.filterDate(sdate_str, edate_str)\
                .select(band)\
                .toBands()

    # export
    task = ee.batch.Export.image.toDrive(
        image=expimg,
        region=expextent,
        scale=expscale,
        crs=expcrs,
        description='E5L_'+band,
        folder='era5land_animation',
        fileNamePrefix=sdate_str.replace('-','')+'-'+edate_str.replace('-','')+'_'+band
    )
    task.start()

    print(sdate_str.replace('-','')+'-'+edate_str.replace('-','')+'_'+band)
