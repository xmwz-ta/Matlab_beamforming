%% 水声阵列信号处理
clear; clc; close all;

%% 参数设置
c = 1500;        % 水中声速(m/s)
fc = 1000;       % 信号中心频率(Hz)
d = 0.75;        % 阵元间距(m)，设置为半波长
lambda = c/fc;

%% 数据读取与预处理
% 读取CSV文件
data = readmatrix('datafinal1.csv'); % 
col1 = data(:,1);
data = data(col1>=0&col1<=23,:);

% 转置数据为阵元数×时间序列
%X = data(:,2:6)'; 

X = [data(:,3), data(:,4), data(:,2), data(:,5), data(:,6), data(:,7), data(:,8), data(:,9), data(:,10)]';

[arrayNum, L] = size(X);

arrayPos = 0:d:(arrayNum-1) * d; %阵元位置

angleIndex = -90:0.5:90;  % 以0.5°为步长扫过整个角度范围

angleIndexLength = length(angleIndex);

targetNum = 1;

%% 信号预处理
% As = zeros(arrayNum,targetNum);  % A(theta)
% S_signal = zeros(targetNum,L);   % s(n)
% for i=1:1:targetNum
%     As(:,i) = exp(-1j*(2*pi/lambda).*arrayPos.'*(sind(theta_T(i))));
%     S_signal(i,:) = exp(1j*2*pi*fc*timeIndex);
% end

%S_signal = As * S_signal;      % X(n) = A(theta) * S(n)

%定义权矢量矩阵 W
W = zeros(arrayNum,angleIndexLength);
for i=1:angleIndexLength
    W(:,i)=1/sqrt(arrayNum).*exp(-1j*(2*pi/lambda).*arrayPos.'*sind(angleIndex(i)));
end
 
%信号协方差矩阵的最大似然估计
R_hat = (X * X.' ) / L;
%R_hatInv = inv(R_hat);
%DOA Algorithms
 
%CBF方法
doa_cbf=zeros(angleIndexLength,1);
 
%P(w)=W^H R W
for i=1:1:angleIndexLength
    doa_cbf(i) = W(:,i)' * R_hat * W(:,i);
end


figure;
plot(angleIndex, 10*log10(doa_cbf), 'b', 'LineWidth', 1.5);
xlabel('Angle (deg)');
ylabel('Power (dB)');
title('DOA Estimation via CBF');
grid on;

% 取模并转为 dB
doa_cbf_db = 10 * log10(abs(doa_cbf));

% 极坐标绘图
figure;
polarplot(deg2rad(angleIndex), doa_cbf_db, 'b', 'LineWidth', 1.5);
title('DOA Estimation via CBF (Polar View)');

% 设置极坐标显示参数
rlim([max(doa_cbf_db)-30, max(doa_cbf_db)]);     % 显示主瓣附近的动态范围
thetalim([-90 90]);                              % 显示角度范围
thetaticks(-90:30:90);                           % 角度刻度
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');  


