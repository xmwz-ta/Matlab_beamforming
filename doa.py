import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# 跳过前面 5 行注释，读取真实数据，header=None 表示无标题行
data = pd.read_csv('datafinal1.csv', skiprows=5, header=None)

# 筛选第 0 列时间在 0~23ms 范围内的行
data[0] = pd.to_numeric(data[0], errors='coerce')  # 转换为数值
data = data[(data[0] >= 0) & (data[0] <= 23)]

# 重排列顺序（根据你之前指定的顺序）
X = data[[2, 3, 1, 4, 5, 6, 7, 8, 9]].to_numpy().T  # 转置：阵元数 × 时间

arrayNum, L = X.shape

# 参数设置
c = 1500  # 声速 (m/s)
fc = 1000  # 中心频率 (Hz)
lambda_ = c / fc  # 波长
d = 0.75  # 阵元间距
arrayPos = np.arange(arrayNum) * d  # 阵元位置

# 角度范围
angleIndex = np.arange(-90, 90.5, 0.5)
angleIndexLength = len(angleIndex)

# 构建权矢量 W
W = np.zeros((arrayNum, angleIndexLength), dtype=complex)
for i in range(angleIndexLength):
    W[:, i] = 1 / np.sqrt(arrayNum) * np.exp(-1j * (2 * np.pi / lambda_) * arrayPos * np.sin(np.radians(angleIndex[i])))

# 协方差矩阵估计
R_hat = (X @ X.conj().T) / L

# CBF计算
doa_cbf = np.zeros(angleIndexLength, dtype=complex)
for i in range(angleIndexLength):
    doa_cbf[i] = W[:, i].conj().T @ R_hat @ W[:, i]

# 转换为 dB
doa_cbf_db = 10 * np.log10(np.abs(doa_cbf))

# 绘图：直角坐标
plt.figure()
plt.plot(angleIndex, doa_cbf_db, 'b', linewidth=1.5)
plt.xlabel('Angle (deg)')
plt.ylabel('Power (dB)')
plt.title('DOA Estimation via CBF')
plt.grid(True)
plt.show()

# 绘图：极坐标图
plt.figure()
ax = plt.subplot(111, polar=True)
theta_rad = np.deg2rad(angleIndex)
ax.plot(theta_rad, doa_cbf_db, 'b', linewidth=1.5)
ax.set_title('DOA Estimation via CBF (Polar View)', va='bottom')

# 极坐标图设置
ax.set_theta_zero_location('N')   # 0° 在顶部
ax.set_theta_direction(-1)        # 顺时针
ax.set_thetamin(-90)
ax.set_thetamax(90)
ax.set_xticks(np.deg2rad(np.arange(-90, 91, 30)))
ax.set_rlim([np.max(doa_cbf_db) - 30, np.max(doa_cbf_db)])  # 限制动态范围

plt.show()
