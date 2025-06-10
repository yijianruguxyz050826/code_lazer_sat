def radar_to_ecef(radar_row):
    lon = radar_row["部署经度(°）"] * u.deg
    lat = radar_row["部署纬度(°）"] * u.deg
    alt = radar_row["部署高度(km)"] * u.km
    loc = EarthLocation(lon=lon, lat=lat, height=alt)
    return loc.get_itrs(obstime=start_time).cartesian.xyz.to(u.km)