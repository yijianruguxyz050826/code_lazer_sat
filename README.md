# 空间目标探测任务规划与评估任务书

## 相关库安装

请使用以下命令安装所需的 Python 库。将 `XXX` 替换为下方列出的库名称：

```bash
pip install -i https://pypi.tuna.tsinghua.edu.cn/simple XXX
```

需要安装的库包括：

- `numpy`
- `pandas`
- `pulp`
- `openpyxl`
- `deap`
- `poliastro`

> **注意：** 若在安装过程中遇到网络限制问题（如大量下载行为被拦截），建议更换网络环境或尝试其他镜像源。

---

## 任务场景描述

本任务旨在利用一组雷达设备对给定的空间目标进行探测任务规划与评估。具体包括：

- 使用多种类型的雷达进行探测，包括**相控阵雷达**和**机械跟踪雷达**。
- 雷达参数详见 `sensorData.xlsx` 文件，包含部署位置、探测能力等信息。
- 每个空间目标的轨道参数及其探测需求（如所需测站数量、最小探测次数、最短探测时间等）详见 `requireData.xlsx` 文件。
- 探测时间段为：**2021年10月14日 4时0分0秒 至 2021年10月15日 4时0分0秒（UTCG时间）**

### 输入数据

- 雷达配置信息
- 空间目标轨道数据
- 探测任务需求

### 输出结果

输出应包含每部雷达在指定时间段内对各目标的详细探测安排，包括：

- 探测起止时间
- 是否执行探测任务
- 探测顺序安排