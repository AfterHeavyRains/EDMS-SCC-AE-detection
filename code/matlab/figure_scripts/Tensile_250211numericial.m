clc; clear; close all;

% 采样频率 & 时间参数
fs = 50;  % 采样率 100 Hz
T = 20;    % 信号持续时间 20 秒
t = 0:1/fs:T-1/fs; % 时间向量

% 复杂调频信号 (参考论文公式)
% signal = 1.7 * sin(14 * pi * t + 3.6 * sin(0.4 * pi * t)) ...
%        + 1.8 * sin(20 * pi * t + 3.2 * sin(0.6 * pi * t));
phi1 = 0.1 * randn(size(t));
phi2 = 0.1 * randn(size(t));

signal = 1.6 * sin(13 * pi * t + 3.2 * sin(0.5 * pi * t) + phi1) ...
       + 1.9 * sin(21 * pi * t + 2.9 * sin(0.7 * pi * t) + phi2);


% 1. 绘制原始信号
figure('Color','w'); %

subplot(3, 1, 1);
plot(t, signal, 'b');
xlabel('Time (s)');
ylabel('Voltage (V)');
title('(a) Simulated AE Crack Signal');
%grid on;

% 2. 计算 CWT 并绘制时频图 (使用 contour 使背景白色)
subplot(3, 1, 2);
[cwt_coeffs, freq] = cwt(signal, 'bump', fs, 'VoicesPerOctave',14);
 % 使用 Morlet 小波计算 CWT
cwt_coeffs_smooth = imgaussfilt(abs(cwt_coeffs), 1.5);  % σ=1 适中
contour(t, freq, cwt_coeffs_smooth, 10, 'LineWidth', 1);
%  contour(t, freq, abs(cwt_coeffs), 15, 'LineWidth', 0.8); % 使用 contour 代替 imagesc
clim([0, max(abs(cwt_coeffs(:)))]);
  % 使颜色范围与论文一致
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([0,20]);
title('(b) Continuous Wavelet Transform (CWT)');
%grid on;
 colormap('hot'); % 设为白色背景
colorbar;

% 3. 计算 WSST (SWT) 并绘制时频图 (使用 contour)
subplot(3, 1, 3);
[sst_coeffs, freq_sst] = wsst(signal, fs, 'bump'); % 计算同步压缩小波变换
contour(t, freq_sst, abs(sst_coeffs), 'LineWidth', 0.8); % 使用 contour 绘制
axis xy;
xlabel('Time (s)');
ylabel('Frequency (Hz)');
title('(c) Synchrosqueezed Wavelet Transform (SWT)');
ylim([0,20]);
%grid on;
colormap('jet'); % 设为 jet 颜色，区分频率成分
colorbar;

% 统一字体为 Times New Roman
set(findall(gcf, '-property', 'FontName'), 'FontName', 'Times New Roman');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 11);

% 保存图像
exportgraphics(gcf, 'contour_cwt_swt.tiff', 'Resolution', 650);
disp('图像已保存为 contour_cwt_swt.tiff');
exportgraphics(gcf, 'fig2.png', 'Resolution', 300);
disp('fig2.png 已保存');
