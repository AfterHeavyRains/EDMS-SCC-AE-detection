clc;
clear;

%% 读取信号
repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
file_path = fullfile(repoRoot, 'data', 'Data5.txt');
signal = load(file_path);
output_folder = fullfile(repoRoot, 'results', 'figures', 'check_run');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end
output_file_name = 'swt_tsallis_20250220_20250108_2mmQ235H_1000N_2_Data5_sym4_0.3.tiff';
output_file_path = fullfile(output_folder, output_file_name);

fs = 2000000; % 采样频率 2 MHz
time = (0:length(signal)-1) / fs;

%% 提取指定时间范围的信号（5ms - 25ms）
time_start = 0e-3;
time_end = 28e-3;
time_idx = (time >= time_start) & (time < time_end);
signal_selected = signal(time_idx);
time_selected = time(time_idx);

%% 滤波
filtered_signal = bandpass_filter(signal_selected, 100000, 700000, fs, 6);

%% 小波去噪
denoised_signal = wavelet_denoise(filtered_signal, 'sym4', 5, 0.3);

%% 滑动窗口参数
window_length = 200e-6;
step_size = 100e-6;
window_samples = floor(window_length * fs);
step_samples = floor(step_size * fs);
num_windows = floor((length(denoised_signal) - window_samples) / step_samples) + 1;

tf_entropy = zeros(1, num_windows);
tf_entropy_band1 = zeros(1, num_windows);
tf_entropy_band2 = zeros(1, num_windows);
tf_entropy_band3 = zeros(1, num_windows);
tf_entropy_sum = zeros(1, num_windows);
tf_error = zeros(1, num_windows);
time_axis = zeros(1, num_windows);

%% 滑动窗口计算 TF 熵
for i = 1:num_windows
    start_idx = (i-1) * step_samples + 1;
    end_idx = start_idx + window_samples - 1;
    if end_idx > length(denoised_signal)
        break;
    end
    
    window_signal = denoised_signal(start_idx:end_idx);
    [sst, freq_axis] = wsst(window_signal, fs, 'bump');
    power_spectrum = abs(sst).^2;
    power_spectrum = power_spectrum ./ sum(power_spectrum, 1);
    
    tf_entropy(i) = tsallis_entropy(power_spectrum, 2);
    
    % 分频段计算 WSST
    band1_signal = bandpass_filter(window_signal, 100000, 300000, fs, 6);
    band2_signal = bandpass_filter(window_signal, 300000, 500000, fs, 6);
    band3_signal = bandpass_filter(window_signal, 500000, 700000, fs, 6);
    
    [sst1, ~] = wsst(band1_signal, fs);
    [sst2, ~] = wsst(band2_signal, fs);
    [sst3, ~] = wsst(band3_signal, fs);
    
    power_band1 = abs(sst1).^2 ./ sum(abs(sst1).^2, 1);
    power_band2 = abs(sst2).^2 ./ sum(abs(sst2).^2, 1);
    power_band3 = abs(sst3).^2 ./ sum(abs(sst3).^2, 1);
    
    tf_entropy_band1(i) = tsallis_entropy(power_band1, 2);
    tf_entropy_band2(i) = tsallis_entropy(power_band2, 2);
    tf_entropy_band3(i) = tsallis_entropy(power_band3, 2);
    
    tf_entropy_sum(i) = tf_entropy_band1(i) + tf_entropy_band2(i) + tf_entropy_band3(i);
    tf_entropy_norm = (tf_entropy - min(tf_entropy)) / (max(tf_entropy) - min(tf_entropy) + eps);
    tf_entropy_sum_norm = (tf_entropy_sum - min(tf_entropy_sum)) / (max(tf_entropy_sum) - min(tf_entropy_sum) + eps);

    if tf_entropy_norm(i) ~= 0
        tf_error(i) = abs((tf_entropy_sum_norm(i) - tf_entropy_norm(i)) / tf_entropy_norm(i)) * 100;
%       if tf_entropy(i) ~= 0 
% tf_error(i) = abs((tf_entropy_sum(i) - tf_entropy(i)) / tf_entropy(i)) * 100;
    else
        tf_error(i) = NaN;
    end
    
    time_axis(i) = mean(time_selected(start_idx:end_idx));
end

%% 绘制图像
% figure;
% subplot(4,1,1);
% plot(time_selected * 1000, denoised_signal,  'Color', [0,0,255]/255,'LineWidth', 1.5);
% xlabel('Time (ms)');
% xlim([0 25]);
% ylabel('Voltage (V)');
% title('(a) Crack Signal (Denoised)');
% % grid on;
% 
% subplot(4,1,2);
% plot(time_axis * 1000, tf_entropy_norm, 'Color', [0,141,255]/255, 'LineWidth', 1.5);
% % plot(time_axis * 1000, tf_entropy, 'm', 'LineWidth', 1.5);
% xlabel('Time (ms)');
% xlim([0 25]);
% ylabel('TF Entropy');
% title('(b) Tsallis Entropy');
% % grid on;
% 
% subplot(4,1,3);
% plot(time_axis * 1000, tf_entropy_sum_norm, 'Color', [123,116,133]/255, 'LineWidth', 1.5);
% % plot(time_axis * 1000, tf_entropy_sum, 'c', 'LineWidth', 1.5);
% xlabel('Time (ms)');
% xlim([0 25]);
% ylabel('Tsallis Entropy Sum');
% title('(c) Summed Tsallis Entropy of Frequency Bands');
% % grid on;
% 
% subplot(4,1,4);
% plot(time_axis * 1000, tf_error, 'k', 'LineWidth', 1.5);
% xlabel('Time (ms)');
% xlim([0 25]);
% ylim([0 10]);
% ylabel('TF Error (%)');
% title('(d) Tsallis Entropy Error Between Summed Bands and Total)');
% % grid on;
figure('Position',[100,100,600,800]); % 竖排四图
set(gcf,'Color','w');

%% (a) Crack Signal
subplot(4,1,1);
plot(time_selected*1000, denoised_signal, 'Color',[0 0.2 0.8],'LineWidth',1.2);
xlim([0 25]);
ylim([-6e-3 6e-3]);
ylabel('Voltage (V)');
title('(a) Crack Signal (Denoised)','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',0.8);
box on;

%% (b) Tsallis Entropy
subplot(4,1,2);
plot(time_axis*1000, tf_entropy_norm, 'Color',[0.4 0.4 0.5],'LineWidth',1.3);
xlim([0 25]);
ylim([-0.1 1]);
ylabel('TF Entropy');
title('(b) Tsallis Entropy','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',0.8);
box on;

%% (c) Summed Tsallis Entropy of Frequency Bands
subplot(4,1,3);
plot(time_axis*1000, tf_entropy_sum_norm, 'b','LineWidth',1.3);
xlim([0 25]);
ylim([-0.1 1]);
ylabel('Entropy Sum');
title('(c) Summed Tsallis Entropy of Frequency Bands','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',0.8);
box on;

%% (d) Tsallis Entropy Error
subplot(4,1,4);
hold on;
% 灰带 ±2% 范围
% err_band = 2*ones(size(tf_error));
% fill([time_axis*1000, fliplr(time_axis*1000)], ...
%      [zeros(size(tf_error)), err_band], ...
%      [0.9 0.9 0.9], 'EdgeColor','none','FaceAlpha',0.4);
% 主线
plot(time_axis*1000, tf_error, 'k','LineWidth',1.3);
% 平均线
% yline(mean(tf_error),'--','Color',[0.5 0.5 0.5],'LineWidth',1);
xlim([0 25]);
ylim([0 10]);
xlabel('Time (ms)');
ylabel('TF Error (%)');
title('(d) Tsallis Entropy Error Between Summed Bands and Total','FontWeight','normal');
set(gca,'FontName','Times New Roman','FontSize',13,'LineWidth',0.8);
box on;

%% 全局调整
% 保持上下间距一致
subplotHandles = get(gcf,'Children');
for i = 1:length(subplotHandles)
    set(subplotHandles(i),'TickDir','out','TickLength',[0.01 0.01]);
end

exportgraphics(gcf, output_file_path, 'Resolution', 650);
disp(['图像已保存至: ', output_file_path]);

%% Tsallis 熵计算函数
function tsallis_H = tsallis_entropy(signal, q)
    prob = abs(signal) / sum(abs(signal));
    tsallis_H = (1 / (q - 1)) * (1 - sum(prob.^q));
end

%% 小波去噪函数
function denoised_signal = wavelet_denoise(signal, wavelet, level, threshold_factor)
    % 小波分解
    [C, L] = wavedec(signal, level, wavelet);
    
    denoised_coeffs = C;
    for i = 2:length(L)  % 跳过近似系数
        coeff_start = sum(L(1:i-1)) + 1;
        coeff_end = sum(L(1:i));

        if coeff_end > length(C)
            warning(['Skipping level ', num2str(i), ' due to index out of bounds.']);
            continue;
        end

        coeff = C(coeff_start:coeff_end);
        sigma = median(abs(coeff)) / 0.6745;
        threshold = sigma * sqrt(2 * log(length(coeff))) * threshold_factor;

        % 软阈值去噪
        denoised_coeffs(coeff_start:coeff_end) = wthresh(coeff, 's', threshold);
    end

    % 信号重构
    denoised_signal = waverec(denoised_coeffs, L, wavelet);
end

%% 带通滤波函数
function filtered_signal = bandpass_filter(signal, low_cutoff, high_cutoff, fs, order)
    [b, a] = butter(order, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
    filtered_signal = filtfilt(b, a, signal);
end
