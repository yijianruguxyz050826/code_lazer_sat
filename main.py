import pandas as pd
import numpy as np
import pulp
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
from astropy.time import Time
import math

# 加载数据
sensor_data = pd.read_excel(r"simData/sensorData.xlsx")
require_data = pd.read_excel(r"simData/requireData.xlsx")

# 参数设置
num_radars = len(sensor_data)
num_targets = len(require_data)
time_slots = 1440  # 24小时 × 60分钟
start_time = Time.now()  # 假设当前时刻为起始时间
dt_per_minute = 60 * u.s  # 每分钟的时间间隔

def create_orbit_from_elements(row):
    a = row["半长轴（km)"] * u.km
    ecc = row["偏心率"] * u.one
    inc = row["轨道倾角(°)"] * u.deg
    raan = row["升交点赤经(°)"] * u.deg
    argp = row["近地点幅角(°）"] * u.deg
    nu = row["平近点角(°)"] * u.deg

    return Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)

from astropy.coordinates import EarthLocation

def radar_to_ecef(radar_row):
    lon = radar_row["部署经度(°）"] * u.deg
    lat = radar_row["部署纬度(°）"] * u.deg
    alt = radar_row["部署高度(km)"] * u.km
    loc = EarthLocation(lon=lon, lat=lat, height=alt)
    return loc.get_itrs(obstime=start_time).cartesian.xyz.to(u.km)

def is_visible(radar_idx, target_idx, time_slot):
    orbit = orbits[target_idx]

    # 使用时间偏移量（相对时间）
    time_of_flight = time_slot * dt_per_minute  # dt_per_minute = 60 * u.s

    try:
        r_target = orbit.propagate(time_of_flight).r.to(u.km)
    except Exception as e:
        print(f"轨道传播失败 at slot {time_slot}: {e}")
        return False

    r_radar = radar_positions[radar_idx]
    distance = np.linalg.norm((r_target - r_radar).value)

    max_distance = sensor_data.loc[radar_idx, "探测距离(km)"]
    return distance <= max_distance

orbits = [create_orbit_from_elements(require_data.iloc[i]) for i in range(num_targets)]
radar_positions = [radar_to_ecef(sensor_data.iloc[i]) for i in range(num_radars)]
required_times = require_data["需要的观测时间(min)"].values
radar_capacities = sensor_data["最大探测目标数"].values

# 创建问题
prob = pulp.LpProblem("Space_Target_Detection_Planning", pulp.LpMaximize)

# 定义变量（仅当可见时才创建）
x = {}

for r in range(num_radars):
    for s in range(num_targets):
        for t in range(time_slots):
            if is_visible(r, s, t):  # 只有可见时才允许分配
                x[(r, s, t)] = pulp.LpVariable(f"x_{r}_{s}_{t}", 0, 1, pulp.LpBinary)

# 添加目标函数：最大化总探测次数
prob += pulp.lpSum([x[r, s, t] for (r, s, t) in x.keys()]), "Total_Detections"

# 添加约束1：每个目标至少被探测 required_times[s] 分钟
for s in range(num_targets):
    prob += pulp.lpSum([x[r, s, t] for (r, s_, t) in x.keys() if s_ == s]) >= required_times[s], f"Target_{s}_Requirement"

# 添加约束2：每个雷达在每分钟最多探测 capacity 个目标
for r in range(num_radars):
    for t in range(time_slots):
        prob += pulp.lpSum([x[r, s, t] for (r_, s, t_) in x.keys() if r_ == r and t_ == t]) <= radar_capacities[r], f"Radar_{r}_Time_{t}_Capacity"
        
# 求解
prob.solve()

# 输出状态
print("Status:", pulp.LpStatus[prob.status])

# 输出变量（只显示非零变量）
for v in prob.variables():
    if v.varValue > 0:
        print(v.name, "=", v.varValue)