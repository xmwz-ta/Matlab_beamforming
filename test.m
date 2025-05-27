%% 水声阵列信号处理 - DOA估计（10°目标信号仿真）
clear; clc; close all;

%% 参数设置
c = 1500;        % 水中声速(m/s)
fc = 1000;       % 信号中心频率(Hz)
lambda = c/fc;   % 波长
d = lambda/2;    % 阵元间距(λ/2)
arrayNum = 8;    % 阵元数
arrayPos = 0:d:(arrayNum-1)*d;  % 阵元位置

sampleRate = 480000;  % 采样率
ts = 1/sampleRate;    % 采样间隔
L = 2000;             % 采样点数
t = (0:L-1)*ts;       % 时间序列

theta_deg = 10;       % 入射角（单位：度）
theta_rad = deg2rad(theta_deg);

%% 构造信号（高斯调制正弦脉冲）
tau = 0.0005;                         % 脉冲宽度控制
t0 = 0.002;                           % 脉冲中心时刻
s = sin(2*pi*fc*t) .* exp(-((t - t0)/tau).^2);  % 发射信号

%% 构造阵列接收信号 X
X = zeros(arrayNum, L);
for m = 1:arrayNum
    delay = arrayPos(m) * sind(theta_deg) / c;  % 到达时间延迟
    s_delayed = sin(2*pi*fc*(t - delay)) .* exp(-((t - delay - t0)/tau).^2);
    X(m, :) = s_delayed;
end

%% 加入噪声（可选）
SNR_dB = 20;
X = awgn(X, SNR_dB, 'measured');

%% DOA估计设置
angleIndex = asin((-512:511)/512) * 180 / pi;  % 角度扫描
angleIndexLength = length(angleIndex);

%% 构造方向矢量矩阵 W
W = zeros(arrayNum, angleIndexLength);
for i = 1:angleIndexLength
    W(:,i) = exp(-1j*(2*pi/lambda)*arrayPos.'*sind(angleIndex(i)));
end

%% 计算协方差矩阵
R_hat = (X * X') / L;

%% CBF方法估计
doa_cbf = zeros(angleIndexLength,1);
for i = 1:angleIndexLength
    doa_cbf(i) = W(:,i)' * R_hat * W(:,i);
end

%% 结果绘图
figure;
plot(angleIndex, 10*log10(abs(doa_cbf)/max(abs(doa_cbf))), 'b', 'LineWidth', 1.5);
xlabel('Angle (deg)');
ylabel('Normalized Power (dB)');
title('DOA Estimation via CBF (Simulated 10° Signal)');
grid on;
