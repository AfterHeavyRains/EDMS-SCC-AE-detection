clc; clear;

%% 读取信号
repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
file_path = fullfile(repoRoot, 'data', 'Data5.txt');
signal = load(file_path);

fs = 2000000; % 采样频率 2 MHz
time = (0:length(signal)-1) / fs; % 时间轴 (秒)

%% 提取指定时间范围的信号（0ms - 25ms）
time_start = 0e-3; % 0 ms
time_end = 25e-3; % 25 ms
time_idx = (time >= time_start) & (time < time_end);
signal_selected = signal(time_idx);
time_selected = time(time_idx);

%% 滤波 & 小波去噪
filtered_signal = bandpass_filter(signal_selected, 100000, 700000, fs, 6);
denoised_signal = wavelet_denoise(filtered_signal, 'sym4', 5, 0.3);

%% 时间窗口参数
window_lengths = [100e-6, 200e-6, 300e-6]; % 200µs, 400µs, 600µs
step_size = 100e-6; % 100 µs 步长

TFSWE_results = struct();

for win_idx = 1:length(window_lengths)
    window_length = window_lengths(win_idx);
    window_samples = floor(window_length * fs);
    step_samples = floor(step_size * fs);
    num_windows = floor((length(denoised_signal) - window_samples) / step_samples) + 1;
    
    tf_entropy = zeros(1, num_windows);
    tf_entropy_band1 = zeros(1, num_windows);
    tf_entropy_band2 = zeros(1, num_windows);
    tf_entropy_band3 = zeros(1, num_windows);
    tf_entropy_sum = zeros(1, num_windows);
    time_axis = zeros(1, num_windows);

    %% 计算 TFSWE 分布
    for i = 1:num_windows
        start_idx = (i-1) * step_samples + 1;
        end_idx = start_idx + window_samples - 1;
        if end_idx > length(denoised_signal), break; end

        window_signal = denoised_signal(start_idx:end_idx);
        [sst, ~] = wsst(window_signal, fs, 'bump');
        power_spectrum = abs(sst).^2;
        power_spectrum = power_spectrum ./ sum(power_spectrum, 1);
        tf_entropy(i) = -sum(power_spectrum .* log2(power_spectrum + eps), 'all') + 1e-6;

        % 三个频段的 WSST
        band1_signal = bandpass_filter(window_signal, 100000, 300000, fs, 6);
        band2_signal = bandpass_filter(window_signal, 300000, 500000, fs, 6);
        band3_signal = bandpass_filter(window_signal, 500000, 700000, fs, 6);
        
        [sst1, ~] = wsst(band1_signal, fs);
        [sst2, ~] = wsst(band2_signal, fs);
        [sst3, ~] = wsst(band3_signal, fs);

        % 计算每个频段的熵值
        power_band1 = abs(sst1).^2 ./ sum(abs(sst1).^2, 1);
        power_band2 = abs(sst2).^2 ./ sum(abs(sst2).^2, 1);
        power_band3 = abs(sst3).^2 ./ sum(abs(sst3).^2, 1);

        tf_entropy_band1(i) = -sum(power_band1 .* log2(power_band1 + eps), 'all');
        tf_entropy_band2(i) = -sum(power_band2 .* log2(power_band2 + eps), 'all');
        tf_entropy_band3(i) = -sum(power_band3 .* log2(power_band3 + eps), 'all');

        % 总和 TFSWE
        tf_entropy_sum(i) = tf_entropy_band1(i) + tf_entropy_band2(i) + tf_entropy_band3(i);
        time_axis(i) = mean(time_selected(start_idx:end_idx));
    end

    % 归一化处理
    tf_entropy_norm = (tf_entropy - min(tf_entropy)) / (max(tf_entropy) - min(tf_entropy));
    tf_entropy_sum_norm = (tf_entropy_sum - min(tf_entropy_sum)) / (max(tf_entropy_sum) - min(tf_entropy_sum));

    % 存储结果
    TFSWE_results(win_idx).time_axis = time_axis;
    TFSWE_results(win_idx).tf_entropy_norm = tf_entropy_norm;
    TFSWE_results(win_idx).tf_entropy_sum_norm = tf_entropy_sum_norm;
end

%% 绘制图像
fnt = 'Times New Roman';
fsz = 15;

fig = figure('Position', [100, 100, 680, 580], 'Color', 'w');

for win_idx = 1:length(window_lengths)
    ax = subplot(3, 1, win_idx);
    plot(TFSWE_results(win_idx).time_axis * 1000, TFSWE_results(win_idx).tf_entropy_norm, ...
        'b', 'LineWidth', 1.2);
    hold on;
    plot(TFSWE_results(win_idx).time_axis * 1000, TFSWE_results(win_idx).tf_entropy_sum_norm, ...
        'r', 'LineWidth', 1.2);
    hold off;
    xlabel('Time (ms)', 'FontName', fnt, 'FontSize', fsz);
    ylabel('Normalized TFSWE', 'FontName', fnt, 'FontSize', fsz);
    title(['(Window = ', num2str(window_lengths(win_idx)*1e6), ' \mus)'], ...
        'FontName', fnt, 'FontSize', fsz, 'FontWeight', 'normal');
    lg = legend({'Total TFSWE', 'Summed Band TFSWE'}, 'Location', 'southeast');
    set(lg, 'FontName', fnt, 'FontSize', fsz - 1);
    set(ax, 'FontName', fnt, 'FontSize', fsz, 'Box', 'on', ...
        'XGrid', 'off', 'YGrid', 'off', 'TickDir', 'out');
end

exportgraphics(fig, 'window.eps', 'ContentType', 'vector');
exportgraphics(fig, 'window.png', 'Resolution', 300);
fprintf('Saved: window.eps / window.png\n');

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
