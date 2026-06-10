clc; clear;
%%本段脚本用于自动识别ae事件 改进前，只有两个图
%% 读取信号
repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
file_path = fullfile(repoRoot, 'data', 'Data5.txt');
signal = load(file_path);
output_folder = fullfile(repoRoot, 'results', 'figures', 'check_run'); % 本地输出目录
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
output_file_name = 'rate.tiff'; % 替换为你的文件名
output_file_path = fullfile(output_folder, output_file_name); % 组合完整路径

fs = 2000000; % 采样频率 2 MHz
time = (0:length(signal)-1) / fs; % 时间轴 (秒)

%% 提取指定时间范围的信号（5ms - 25ms）
time_start = 0e-3; 
time_end = 25e-3; 
time_idx = (time >= time_start) & (time < time_end);

signal_selected = signal(time_idx);
time_selected = time(time_idx);

%% 滤波 & 小波去噪
filtered_signal = bandpass_filter(signal_selected, 100000, 700000, fs, 6);
denoised_signal = wavelet_denoise(filtered_signal, 'sym4', 5, 0.3);

%% 滑动窗口参数
window_length = 200e-6; 
step_size = 100e-6; 
window_samples = floor(window_length * fs);
step_samples = floor(step_size * fs);
num_windows = floor((length(denoised_signal) - window_samples) / step_samples) + 1;

%% 计算 TF 熵
tf_entropy = zeros(1, num_windows);
time_axis = zeros(1, num_windows);

for i = 1:num_windows
    start_idx = (i-1) * step_samples + 1;
    end_idx = start_idx + window_samples - 1;
    if end_idx > length(denoised_signal), break; end

    window_signal = denoised_signal(start_idx:end_idx);
    [sst, ~] = wsst(window_signal, fs, 'bump');
    power_spectrum = abs(sst).^2;
    power_spectrum = power_spectrum ./ sum(power_spectrum, 1);

    % 计算 Shannon 熵
    tf_entropy(i) = -sum(power_spectrum .* log2(power_spectrum + eps), 'all') + 1e-6;
    tf_entropy_norm = (tf_entropy - min(tf_entropy)) / (max(tf_entropy) - min(tf_entropy));
    % 记录窗口中心时间
    time_axis(i) = mean(time_selected(start_idx:end_idx));
end

%% 计算熵变化率（微分）
entropy_rate = diff(tf_entropy) ./ diff(time_axis);
time_rate = time_axis(2:end); % 时间对齐


%% 计算 AE 事件触发阈值
mean_rate = mean(entropy_rate);
std_rate = std(entropy_rate);
threshold = mean_rate - 3.5 * std_rate;  % 2.8设定

%% 检测 AE 事件（熵变化率低于阈值）
AE_events = find(entropy_rate < threshold);
AE_times = time_rate(AE_events) * 1000; % 转换为 ms

%% 绘制熵变化率及 AE 事件
figure;
commonFont = 'Times New Roman';
commonFS = 13;
subplot(2,1,1);
plot(time_axis * 1000, tf_entropy, 'b', 'LineWidth', 1.5);
hold on;
scatter(AE_times, tf_entropy(AE_events+1), 50, 'r', 'filled'); % 标记 AE 事件
xlabel('Time (ms)');
ylabel('Tsallis Entropy');
title('(a) Tsallis Entropy and AE Event Localization');
set(gca,'FontName',commonFont,'FontSize',commonFS,'LineWidth',0.8);
% grid on;
legend({'Tsallis Entropy', 'AE Events'}, 'Location', 'southeast');

subplot(2,1,2);
plot(time_rate * 1000, entropy_rate, 'k', 'LineWidth', 1.5);
hold on;
yline(threshold, 'r--', 'LineWidth', 2); % 画出阈值
scatter(AE_times, entropy_rate(AE_events), 50, 'r', 'filled'); % 标记 AE 事件
xlabel('Time (ms)');
ylabel('Entropy Rate');
set(gca,'FontName',commonFont,'FontSize',commonFS,'LineWidth',0.8);
title('(b) Entropy Rate Evolution during AE Process');
% grid on;
legend({'Entropy Rate', 'Adaptive Baseline', 'AE Events'}, 'Location', 'southeast');
% 保存图片
exportgraphics(gcf, output_file_path, 'Resolution', 650); % 设置分辨率为 300 DPI
disp(['图像已保存至: ', output_file_path]);
%% 输出 AE 事件信息
disp(['AE 事件发生在时间点（ms）: ', num2str(AE_times)]);
disp(['建议的 AE 触发阈值: ', num2str(threshold)]);
%% 保存 AE 事件时间列表（单位：毫秒）为 txt 文件
output_txt_file = fullfile(output_folder, 'AE_detected_times_ms.txt');
writematrix(AE_times', output_txt_file, 'Delimiter', 'tab');
disp(['AE 检测时间点（单位 ms）已保存至：', output_txt_file]);


%% 小波去噪函数
function denoised_signal = wavelet_denoise(signal, wavelet, level, threshold_factor)
    [C, L] = wavedec(signal, level, wavelet);
    denoised_coeffs = C;
    for i = 2:length(L)
        coeff_start = sum(L(1:i-1)) + 1;
        coeff_end = sum(L(1:i));
        if coeff_end > length(C), continue; end
        coeff = C(coeff_start:coeff_end);
        sigma = median(abs(coeff)) / 0.6745;
        threshold = sigma * sqrt(2 * log(length(coeff))) * threshold_factor;
        denoised_coeffs(coeff_start:coeff_end) = wthresh(coeff, 's', threshold);
    end
    denoised_signal = waverec(denoised_coeffs, L, wavelet);
end

%% 带通滤波函数
function filtered_signal = bandpass_filter(signal, low_cutoff, high_cutoff, fs, order)
    [b, a] = butter(order, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
    filtered_signal = filtfilt(b, a, signal);
end
