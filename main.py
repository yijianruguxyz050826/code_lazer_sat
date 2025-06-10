import pandas as pd
import pulp
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
from astropy.coordinates import EarthLocation

## 加载数据
sensor_data = pd.read_excel(r"simData/sensorData.xlsx")
require_data = pd.read_excel(r"simData/requireData.xlsx")

## 参数设置
num_radars = len(sensor_data)
num_targets = len(require_data)

# 设置探测时间段：2021年10月14日 4:00:00 到 2021年10月15日 4:00:00 UTCG
start_time_str = "2021-10-14T04:00:00"
end_time_str = "2021-10-15T04:00:00"

# 使用 astropy 定义开始和结束时间（注意：Time 默认是 UTC 时间）
start_time = Time(start_time_str, format='isot', scale='utc')
end_time = Time(end_time_str, format='isot', scale='utc')

# 计算总时间片数量（按每分钟划分）
dt_per_minute = 60 * u.s  # 每分钟的时间间隔
total_seconds = (end_time - start_time) * 86400  # 转换为秒
time_slots = int(total_seconds.value // dt_per_minute.value)

## 相关函数定义
def create_orbit_from_elements(row):
    a = row["半长轴（km)"] * u.km
    ecc = row["偏心率"] * u.one
    inc = row["轨道倾角(°)"] * u.deg
    raan = row["升交点赤经(°)"] * u.deg
    argp = row["近地点幅角(°）"] * u.deg
    nu = row["平近点角(°)"] * u.deg

    return Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)

def radar_to_ecef(radar_row):
    lon = radar_row["部署经度(°）"] * u.deg
    lat = radar_row["部署纬度(°）"] * u.deg
    alt = radar_row["部署高度(km)"] * u.km
    loc = EarthLocation(lon=lon, lat=lat, height=alt)
    return loc.get_itrs(obstime=start_time).cartesian.xyz.to(u.km)

orbits = [create_orbit_from_elements(require_data.iloc[i]) for i in range(num_targets)]
radar_positions = [radar_to_ecef(sensor_data.iloc[i]) for i in range(num_radars)]
required_times = require_data["需要的观测时间(min)"].values
radar_capacities = sensor_data["最大探测目标数"].values

# 创建问题
prob = pulp.LpProblem("Space_Target_Detection_Planning", pulp.LpMaximize)