import pulp
import numpy as np
import pandas as pd
from deap import base, creator, tools

# 加载数据
sensor_data = pd.read_excel("sensorData.xlsx")
object_data = np.loadtxt("objectData.txt")  # TLE 数据
require_data = pd.read_excel("requireData.xlsx")

# 参数设置
num_radars = len(sensor_data)
num_targets = len(object_data)
time_slots = 1440  # 24小时 × 60分钟

# 定义变量
x = pulp.LpVariable.dicts("x", [(r, s, t) for r in range(num_radars)
                                 for s in range(num_targets)
                                 for t in range(time_slots)], 0, 1, pulp.LpBinary)

# 创建问题
prob = pulp.LpProblem("Space_Target_Detection_Planning", pulp.LpMaximize)

# 添加目标函数
prob += pulp.lpSum([x[r][s][t] for r in range(num_radars)
                    for s in range(num_targets)
                    for t in range(time_slots)])

# 添加约束
for s in range(num_targets):
    prob += pulp.lpSum([x[r][s][t] for r in range(num_radars)
                        for t in range(time_slots)]) >= require_data.loc[s, "所需探测次数"]

for r in range(num_radars):
    for t in range(time_slots):
        prob += pulp.lpSum([x[r][s][t] for s in range(num_targets)]) <= sensor_data.loc[r, "目标容量"]

# 求解
prob.solve()

# 输出结果
print("Status:", pulp.LpStatus[prob.status])
for v in prob.variables():
    if v.varValue > 0:
        print(v.name, "=", v.varValue)