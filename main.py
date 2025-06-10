import pandas as pd
import pulp
from astropy.time import Time
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from astropy import units as u
from astropy.coordinates import EarthLocation

## 加载数据
print("数据加载开始...")
sensor_data = None
require_data = None
try:
    sensor_path = r"simData/sensorData.xlsx"
    require_path = r"simData/requireData.xlsx"
    sensor_data = pd.read_excel(sensor_path)
    require_data = pd.read_excel(require_path)
except FileNotFoundError as e:
    print(f"文件未找到: {e}")
except pd.errors.ParserError as e:
    print(f"Excel 文件解析失败: {e}")
except Exception as e:
    print(f"发生未知错误: {e}")
else:
    print("数据加载成功")
finally:
    print("数据加载结束。")

## 参数设置
print("参数设置开始...")
# 检查 sensor_data 和 require_data 是否已加载成功
if sensor_data is None or require_data is None:
    raise ValueError("数据未正确加载，请先完成数据读取。")

try:
    num_radars = len(sensor_data)
    num_targets = len(require_data)

    # 设置探测时间段：UTC 时间
    start_time_str = "2021-10-14T04:00:00"
    end_time_str = "2021-10-15T04:00:00"

    # 使用 astropy 定义时间
    start_time = Time(start_time_str, format='isot', scale='utc')
    end_time = Time(end_time_str, format='isot', scale='utc')

    # 计算总时间片数量（按每分钟划分）
    dt_per_minute = 60 * u.s  # 每分钟的时间间隔
    total_seconds = (end_time - start_time) * 86400  # 转换为秒
    time_slots = int(total_seconds.value // dt_per_minute.value)

    ## 相关函数定义
    def create_orbit_from_elements(row):
        try:
            a = row["半长轴（km)"] * u.km
            ecc = row["偏心率"] * u.one
            inc = row["轨道倾角(°)"] * u.deg
            raan = row["升交点赤经(°)"] * u.deg
            argp = row["近地点幅角(°）"] * u.deg
            nu = row["平近点角(°)"] * u.deg

            return Orbit.from_classical(Earth, a, ecc, inc, raan, argp, nu)
        except KeyError as e:
            raise KeyError(f"轨道参数列缺失: {e}")
        except Exception as e:
            raise ValueError(f"轨道参数转换失败: {e}")

    def radar_to_ecef(radar_row):
        try:
            lon = radar_row["部署经度(°）"] * u.deg
            lat = radar_row["部署纬度(°）"] * u.deg
            alt = radar_row["部署高度(km)"] * u.km
            loc = EarthLocation(lon=lon, lat=lat, height=alt)
            return loc.get_itrs(obstime=start_time).cartesian.xyz.to(u.km)
        except KeyError as e:
            raise KeyError(f"雷达部署信息列缺失: {e}")
        except Exception as e:
            raise ValueError(f"雷达位置转换失败: {e}")

    # 创建轨道和雷达位置列表
    orbits = [create_orbit_from_elements(require_data.iloc[i]) for i in range(num_targets)]
    radar_positions = [radar_to_ecef(sensor_data.iloc[i]) for i in range(num_radars)]

    # 提取观测时间和雷达容量
    required_times = require_data["需要的观测时间(min)"].values
    radar_capacities = sensor_data["最大探测目标数"].values

    # 创建优化问题
    prob = pulp.LpProblem("Space_Target_Detection_Planning", pulp.LpMaximize)

except KeyError as e:
    print(f"数据列缺失或名称错误: {e}")
except TypeError as e:
    print(f"数据类型错误或未正确加载: {e}")
except Exception as e:
    print(f"参数设置过程中发生错误: {e}")
else:
    print("参数设置成功")
finally:
    print("参数设置结束")