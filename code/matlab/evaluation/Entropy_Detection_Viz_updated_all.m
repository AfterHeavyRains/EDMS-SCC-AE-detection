%% =============================================================
% AE Event Detection (Method 1: Multi-step Composite Detection) with Pseudo Labeling & Visualization
% Author: [Your Name]   Date: 2025-08-26
%
% 功能概览：
%   ✅ 多步骤复合检测（归一化熵 + 熵率辅助 + 能量过滤 + 突发性检测）
%   ✅ 分层伪标签（高置信/低置信/包络-only）
%   ✅ 三联可视化图：信号+包络 / 熵值 / 熵率
%% =============================================================

clc; clear;
tic
%% ---------------- 1. 基础设置 ----------------
repoRoot = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
% Raw 8 mm waveform folder, expects Data1.txt, Data2.txt, ...
% The full per-thickness waveform dataset is available from the corresponding author on request;
% place it here (or repoint data_folder) to reproduce the comparison figures.
data_folder   = fullfile(repoRoot, 'data', 'raw_8mm');
output_folder = fullfile(repoRoot, 'results', 'viz_run');
if ~exist(output_folder,'dir'), mkdir(output_folder); end

fs = 2e6;
time_start = 0e-3; 
time_end   = 1000e-3;
num_files  = 2;

AGR_all = zeros(num_files,1);
EOR_all = zeros(num_files,1);
AGR_orig_all = zeros(num_files,1);   % ⭐ 新增
EOR_orig_all = zeros(num_files,1);   % ⭐ 新增
dt_all_orig   = [];
dt_all_layered= [];
weak_quality_all = [];
counts_orig   = [0 0 0 0 0]; % HC, LC, MC, WC, EO
counts_layered= [0 0 0 0 0];

% 只挑部分样本可视化
% plot_samples = [1, 50, 100, 150, 200];
plot_samples = [2];
for sample_id = 1:num_files
    fprintf('\n================ 样本 %d =================\n', sample_id);

    % 原始文件路径
    file_path = fullfile(data_folder, sprintf('Data%d.txt', sample_id));
    signal = load(file_path);
    time   = (0:numel(signal)-1)/fs;
    
    % 预处理
%     signal_denoised = wavelet_denoise( bandpass_filter(signal,1e5,7e5,fs,6), 'sym4', 5, 0.3);
    signal_denoised = signal;

    % 如果需要，仍然保存处理结果（可选）
%     % 创建一个新的子文件夹保存降噪结果
%     denoised_folder = fullfile(output_folder, 'data_denoised');
%     if ~exist(denoised_folder, 'dir')
%         mkdir(denoised_folder);
%     end
    
    % 保存为和原始文件同名的 .mat
%     [~, fname, ~] = fileparts(file_path);  % 比如 Data1
%     denoised_file = fullfile(denoised_folder, [fname '.mat']);
%     save(denoised_file, 'signal_denoised', 'fs');


    %% ---------------- 2. 滑窗熵计算 + 熵率 ----------------
    win_len = 200e-6; step = 100e-6;
    win_pts = floor(win_len*fs); 
    step_pts = floor(step*fs);

    start_idx  = 1:step_pts:(numel(signal)-win_pts+1);
    center_idx = start_idx + floor(win_pts/2);
    time_axis  = time(center_idx);

    % 工具函数（单位转换）
    rate_dt  = step_pts / fs;                         
    ms2rate  = @(ms) max(1, round((ms*1e-3) / rate_dt)); 
    ms2samp  = @(ms) max(1, round(ms*1e-3 * fs));        

    % 熵计算
    W = numel(start_idx);
    q = 1.6; 
    tf_entropy = zeros(1,W);

    for i = 1:W
        seg = signal_denoised(start_idx(i):start_idx(i)+win_pts-1);
        [sst, ~] = wsst(seg, fs, 'bump');
        P = abs(sst).^2; P = P ./ max(sum(P,1), eps);
        tf_entropy(i) = (1 - sum(P(:).^q)) / (q - 1);
    end

    %% ---------------- 3. 多步骤复合检测 ----------------
    % (1) 归一化熵值
    tf_entropy_norm = (tf_entropy - min(tf_entropy)) / (max(tf_entropy) - min(tf_entropy));
    
    % (2) 熵变化率
    entropy_rate = diff(tf_entropy_norm) ./ diff(time_axis);
    time_rate    = time_axis(2:end);
    
    % (3) 熵阈值 + 熵率辅助
%     entropy_thresh = 0.38;
    entropy_thresh = 0.48;
    mean_rate = mean(entropy_rate); 
    std_rate  = std(entropy_rate);
    rate_thresh = mean_rate - 4.0 * std_rate;
    low_entropy_idx = find(tf_entropy_norm(2:end) < entropy_thresh);
    AE_events_combined = [];
    
    % ±1 ms 搜索
    wR = ms2rate(1.0);
    for idx = low_entropy_idx'
        left  = max(1, idx - wR);
        right = min(length(entropy_rate), idx + wR);
        if min(entropy_rate(left:right)) < rate_thresh
            AE_events_combined = [AE_events_combined, idx];
        end
    end
    AE_events_combined = unique(AE_events_combined(:));
    
    % (4) 突发性 + 能量过滤
    window_energy = ms2samp(1.0);
    energy_bg_global = movsum(signal_denoised.^2, window_energy);
    energy_thresh = 1.1 * median(energy_bg_global); %energy_global_mult
    burst_thresh  = 1.2;
    
    AE_events_filtered_burst = [];
    for k = 1:numel(AE_events_combined)
        rate_idx = AE_events_combined(k);
        samp_idx = center_idx(rate_idx);
        left  = max(1, samp_idx - window_energy);
        right = min(numel(signal_denoised), samp_idx + window_energy);
        win_sig = signal_denoised(left:right);
        energy_local = sum(win_sig.^2);
        burstiness   = kurtosis(abs(win_sig));
        if energy_local > energy_thresh && burstiness > burst_thresh
            AE_events_filtered_burst = [AE_events_filtered_burst; rate_idx];
        end
    end
    AE_events_combined = AE_events_filtered_burst;
    
    % (5) 动态背景能量过滤
    window_energy_small = ms2samp(0.05);
    window_background   = ms2samp(2.5);
    alpha = 3.0;    %alpha
    
    AE_events_filtered_dynamic = [];
    for k = 1:numel(AE_events_combined)
        rate_idx = AE_events_combined(k);
        samp_idx = center_idx(rate_idx);
        % 局部能量
        left  = max(1, samp_idx - window_energy_small);
        right = min(numel(signal_denoised), samp_idx + window_energy_small);
        energy_local = sum(signal_denoised(left:right).^2);
        % 背景能量
        left_bg  = max(1, samp_idx - window_background);
        right_bg = min(numel(signal_denoised), samp_idx + window_background);
        bg_seg   = signal_denoised(left_bg:right_bg);
        energy_bg_series = movsum(bg_seg.^2, window_energy_small);
        mean_bg = mean(energy_bg_series);
        std_bg  = std(energy_bg_series);
        local_energy_thresh = mean_bg + alpha * std_bg;
        if energy_local > local_energy_thresh
            AE_events_filtered_dynamic = [AE_events_filtered_dynamic; rate_idx];
        end
    end
    AE_events_combined = AE_events_filtered_dynamic;
    
    % (6) 时间转换 + 最小间隔过滤
    if isempty(AE_events_combined)
        fprintf('⚠️ 无检测事件通过能量/突发性判据，取前10个候选.\n');
        AE_events_combined = low_entropy_idx(1:min(10, length(low_entropy_idx)))';
    end
    
    AE_times_all = sort(time_axis(AE_events_combined)*1000);
    t_entropy = [];
    min_interval_ms = 5;
    for i = 1:numel(AE_times_all)
        if isempty(t_entropy) || (AE_times_all(i)-t_entropy(end) > min_interval_ms)
            t_entropy(end+1,1) = AE_times_all(i);
        end
    end

    %% ---------------- 4. 包络候选 (STA/LTA) ----------------
    stalta_evt_ms = stalta_detect(signal_denoised, fs, time, time_start, time_end);

    %% ---------------- 5. 分层伪标签 (原始匹配) ----------------
    tol_ms = 0.8;
    [HC, LC, EO, ~, ~, ~, ~, dt_signed] = match_pseudo_labels(t_entropy, stalta_evt_ms, tol_ms);
    fprintf('样本 %d: HC=%d, LC=%d, EO=%d\n', sample_id, length(HC), length(LC), length(EO));

        %% ---------------- 6. 微弱事件检测 ----------------
    fprintf('\n=== 启动微弱事件检测通道 ===\n');
    
    % (1) 放宽检测参数
%     entropy_thresh_weak = 0.50;  % 放宽熵阈值
    entropy_thresh_weak = 0.69; 
    mean_rate = mean(entropy_rate);
    std_rate  = std(entropy_rate);
    rate_thresh_weak = mean_rate - 0.2 * std_rate;
    burst_thresh_weak = 1.1;           % 降低突发性要求
    alpha_weak        = 2.5;           % 放宽能量要求
    
    % (2) 候选事件（熵 + 熵率）
    low_entropy_weak = find(tf_entropy_norm(2:end) < entropy_thresh_weak & ...
                           tf_entropy_norm(2:end) >= 0.38);   % 排除主通道已检测的
    AE_events_weak_candidates = [];
    wR = ms2rate(1.5);  % 稍大的搜索窗口
    
    for idx = low_entropy_weak'
        left  = max(1, idx - wR);
        right = min(length(entropy_rate), idx + wR);
        if min(entropy_rate(left:right)) < rate_thresh_weak
            AE_events_weak_candidates = [AE_events_weak_candidates, idx];
        end
    end
    
    if isempty(AE_events_weak_candidates)
        fprintf('⚠️ 无微弱候选，进一步放宽阈值...\n');
        entropy_thresh_weak = 0.55;
        low_entropy_weak = find(tf_entropy_norm(2:end) < entropy_thresh_weak & ...
                               tf_entropy_norm(2:end) >= entropy_thresh);
        for idx = low_entropy_weak'
            left  = max(1, idx - wR);
            right = min(length(entropy_rate), idx + wR);
            if min(entropy_rate(left:right)) < rate_thresh_weak
                AE_events_weak_candidates = [AE_events_weak_candidates, idx];
            end
        end
    end
    fprintf('微弱候选事件数量: %d\n', length(AE_events_weak_candidates));
    
    % (3) 候选质量评估
    t_weak = [];
    weak_quality_scores = [];
    
    % 主事件索引用于避免重复
    main_events_idx = [];
    for i = 1:length(t_entropy)
        [~, idx] = min(abs(time_axis*1000 - t_entropy(i)));
        main_events_idx = [main_events_idx, idx];
    end
    
    for k = 1:length(AE_events_weak_candidates)
        rate_idx = AE_events_weak_candidates(k);
    
        % 跳过与主事件过近的候选
        if any(abs(main_events_idx - rate_idx) <= 5), continue; end
    
        samp_idx   = center_idx(rate_idx);
        event_time = time_axis(rate_idx) * 1000;  % ms
    
        % (3.1) 能量 + 突发性检查
        window_energy_weak = ms2samp(0.1);
        left  = max(1, samp_idx - window_energy_weak);
        right = min(numel(signal_denoised), samp_idx + window_energy_weak);
        win_sig = signal_denoised(left:right);
        energy_local = sum(win_sig.^2);
        burstiness   = kurtosis(abs(win_sig));
    
        % 动态背景能量
        window_background_weak = ms2samp(2.0);
        left_bg  = max(1, samp_idx - window_background_weak);
        right_bg = min(numel(signal_denoised), samp_idx + window_background_weak);
        bg_seg   = signal_denoised(left_bg:right_bg);
        energy_bg_series = movsum(bg_seg.^2, window_energy_weak);
        mean_bg = mean(energy_bg_series); std_bg = std(energy_bg_series);
        local_energy_thresh = mean_bg + alpha_weak * std_bg;
    
        if energy_local <= local_energy_thresh || burstiness <= burst_thresh_weak
            continue;
        end
    
        % (3.2) 多维特征打分
        quality_features = zeros(1,6);
    
        % 特征1: 熵值
        quality_features(1) = max(0, 1 - (tf_entropy_norm(rate_idx)-0.38)/0.17);
        % 特征2: 熵率显著性
        rate_significance = abs(entropy_rate(rate_idx-1) - mean_rate) / (std_rate+eps);
        quality_features(2) = min(1, rate_significance/3);
        % 特征3: 能量比
        energy_ratio = energy_local / (mean_bg+eps);
        quality_features(3) = max(0, min(1, log10(energy_ratio)/1.5));
        % 特征4: 突发性
        quality_features(4) = max(0, min(1, (burstiness-1)/2));
        % 特征5: SNR
        if length(win_sig) >= 8
            signal_power = max(movsum(win_sig.^2, max(1, round(length(win_sig)/4))));
            noise_power  = median(movsum(win_sig.^2, max(1, round(length(win_sig)/4))));
            snr = signal_power / (noise_power+eps);
            quality_features(5) = max(0, min(1, log10(snr)/2));
        else
            quality_features(5) = 0.3;
        end
        % 特征6: 频域特征
        if length(win_sig) >= 32
            try
                [psd,f] = pwelch(win_sig, [], [], [], fs);
                ae_band = (f>=1e5) & (f<=6e5);
                noise_band = f<5e4;
                ae_power = sum(psd(ae_band));
                noise_power_freq = sum(psd(noise_band));
                spectral_ratio = ae_power/(noise_power_freq+eps);
                quality_features(6) = max(0, min(1, log10(spectral_ratio)/2));
            catch
                quality_features(6) = 0.3;
            end
        else
            quality_features(6) = 0.3;
        end
    
        % (3.3) 综合得分
        weights_weak = [0.25,0.20,0.15,0.15,0.15,0.10];
        overall_quality = sum(weights_weak .* quality_features);
    
        if overall_quality >= 0.4  %weak_quality_cut
            t_weak = [t_weak; event_time];
            weak_quality_scores = [weak_quality_scores; overall_quality];
        end
    end
    
    % (4) 最小间隔过滤
    if ~isempty(t_weak)
        [t_weak, sort_idx] = sort(t_weak);
        weak_quality_scores = weak_quality_scores(sort_idx);
        min_interval_weak = 3.0;  %min_interval_weak
        keep_mask = true(length(t_weak),1);
        for i = 2:length(t_weak)
            if t_weak(i)-t_weak(i-1) < min_interval_weak
                if weak_quality_scores(i) > weak_quality_scores(i-1)
                    keep_mask(i-1)=false;
                else
                    keep_mask(i)=false;
                end
            end
        end
        t_weak = t_weak(keep_mask);
        weak_quality_scores = weak_quality_scores(keep_mask);
    end
    
    % (5) 质量统计
    fprintf('微弱事件检测完成: %d个事件\n', length(t_weak));
    if ~isempty(weak_quality_scores)
        fprintf('质量得分范围: %.3f - %.3f (平均: %.3f)\n', ...
            min(weak_quality_scores), max(weak_quality_scores), mean(weak_quality_scores));
        fprintf('质量分级: 高(%d), 中(%d), 低(%d)\n', ...
            sum(weak_quality_scores>=0.6), ...
            sum(weak_quality_scores>=0.5 & weak_quality_scores<0.6), ...
            sum(weak_quality_scores<0.5));
    end

   %% ---------------- 7. 分层伪标签 ---------------- 
    if ~isempty(t_weak)
    
        % (1) 主通道匹配（HC/LC）
        used_env_main = false(length(stalta_evt_ms),1);
        for i = 1:length(HC)
            [~,closest_idx] = min(abs(stalta_evt_ms-HC(i)));
            if abs(stalta_evt_ms(closest_idx)-HC(i)) <= tol_ms
                used_env_main(closest_idx) = true;
            end
        end
        remaining_env = stalta_evt_ms(~used_env_main);
    
        % (2) 微弱通道匹配（MC / WC）
        MC=[]; WC=[];
        used_env_weak = false(length(remaining_env),1);
    
        for i = 1:length(t_weak)
            weak_time = t_weak(i); 
            weak_quality = weak_quality_scores(i);
    
            if ~isempty(remaining_env)
                [min_dist,closest_idx] = min(abs(remaining_env-weak_time));
                if min_dist<=tol_ms && ~used_env_weak(closest_idx)
                    % 匹配成功：weak+包络 → MC
                    MC = [MC; weak_time];
                    used_env_weak(closest_idx) = true;
                else
                    % 匹配失败：weak没包络 → WC（候选，不算事件）
                    WC = [WC; weak_time];
                end
            else
                % 没有剩余包络 → WC
                WC = [WC; weak_time];
            end
        end
    
        % (3) 最终分类
        HC_layered = HC; 
        LC_layered = LC; 
        MC_layered = MC; 
        WC_layered = WC;
        EO_layered = remaining_env(~used_env_weak);
    
    else
        fprintf('⚠️ 微弱通道无事件，使用原始标签\n');
        HC_layered=HC; MC_layered=[]; LC_layered=LC; WC_layered=[]; EO_layered=EO;
    end
    
    % ⭐ 新增：计算分层 Δt (HC + MC)
    dt_signed_layered = [];
    
    % HC
    for i = 1:numel(HC_layered)
        [dmin,j] = min(abs(stalta_evt_ms - HC_layered(i)));
        if dmin <= tol_ms
            dt_signed_layered(end+1,1) = HC_layered(i) - stalta_evt_ms(j);
        end
    end
    
    % MC
    for i = 1:numel(MC_layered)
        [dmin,j] = min(abs(stalta_evt_ms - MC_layered(i)));
        if dmin <= tol_ms
            dt_signed_layered(end+1,1) = MC_layered(i) - stalta_evt_ms(j);
        end
    end
    
    % 打印统计
    fprintf('样本 %d (分层): HC=%d, LC=%d, MC=%d, WC=%d, EO=%d | Δt分层=%d\n', ...
        sample_id, length(HC_layered), length(LC_layered), length(MC_layered), ...
        length(WC_layered), length(EO_layered), length(dt_signed_layered));

%% -------- 8. 导出与评估（精简版） --------
    sample_folder = fullfile(output_folder, sprintf('Data%d', sample_id));
    if ~exist(sample_folder,'dir'), mkdir(sample_folder); end
    
    % --- 原始方案（假定 HC/LC/EO 已存在） ---
    AGR_orig = numel(HC) / max(1, (numel(HC)+numel(LC)));
    EOR_orig = numel(EO) / max(1, (numel(HC)+numel(EO)));
    writematrix(HC(:), fullfile(sample_folder,'AE_eval_HC_orig.txt'));
    writematrix(LC(:), fullfile(sample_folder,'AE_eval_LC_orig.txt'));
    writematrix(EO(:), fullfile(sample_folder,'AE_eval_EO_orig.txt'));
    
    % --- 分层方案（若存在） ---
    has_layered = exist('HC_layered','var') && ~isempty(HC_layered) && ...
                  exist('MC_layered','var') && ~isempty(MC_layered) && ...
                  exist('LC_layered','var') && exist('EO_layered','var');
    
    if has_layered
        Final = sort([HC_layered(:); MC_layered(:)]);
        hi_conf = numel(HC_layered) + numel(MC_layered);
        tot_ent = hi_conf + numel(LC_layered);
    
        AGR_final = hi_conf / max(1, tot_ent);
        EOR_final = numel(EO_layered) / max(1, (hi_conf + numel(EO_layered)));
        Coverage  = numel(Final) / max(1, (numel(Final) + numel(EO_layered)));
    
        writematrix(Final,             fullfile(sample_folder,'AE_final_times_ms.txt'));
        writematrix(HC_layered(:),     fullfile(sample_folder,'AE_eval_HC.txt'));
        writematrix(MC_layered(:),     fullfile(sample_folder,'AE_eval_MC.txt'));
        writematrix(LC_layered(:),     fullfile(sample_folder,'AE_eval_LC.txt'));
        if exist('WC_layered','var'), writematrix(WC_layered(:), fullfile(sample_folder,'AE_eval_WC.txt')); end
        writematrix(EO_layered(:),     fullfile(sample_folder,'AE_eval_EO.txt'));
    
    else
        % 无分层就把原始的最终清单也存一下（按你原流程可换成 t_entropy）
        Final     = sort([HC(:)]);          % 最小替代：仅 HC 作为“发布清单”
        AGR_final = AGR_orig;               % 指标退回原始
        EOR_final = EOR_orig;
        Coverage  = NaN;
        writematrix(Final, fullfile(sample_folder,'AE_final_times_ms.txt'));
    end
    
    % --- 写入汇总数组 ---
    AGR_all(sample_id)      = AGR_final;
    EOR_all(sample_id)      = EOR_final;
    AGR_orig_all(sample_id) = AGR_orig;
    EOR_orig_all(sample_id) = EOR_orig;
    
    % （可选）类别计数（精简，不判空）
    counts_orig    = counts_orig    + [numel(HC), numel(LC), 0, 0, numel(EO)];
    if has_layered
        counts_layered = counts_layered + [numel(HC_layered), numel(LC_layered), ...
                                           numel(MC_layered), exist('WC_layered','var')*numel(WC_layered), ...
                                           numel(EO_layered)];
    end
    fprintf('样本 %d: 原始Δt数=%d, 分层Δt数=%d\n', ...
        sample_id, ...
        (exist('dt_signed','var') && ~isempty(dt_signed)) * length(dt_signed), ...
        (exist('dt_signed_layered','var') && ~isempty(dt_signed_layered)) * length(dt_signed_layered));
    fprintf('✅ 样本 %d 检测完成！结果已保存到 %s\n', sample_id, sample_folder);
    % ====== 汇总表：打印 + 追加保存到 CSV ======
    summary_csv = fullfile(output_folder,'AGR_EOR_summary.csv');
    
    % 事件计数（无分层时 MC/WC 记 0）
    nHC = numel(HC);
    nLC = numel(LC);
    nEO = numel(EO);
    nMC = (exist('MC_layered','var') && ~isempty(MC_layered)) * numel(MC_layered);
    nWC = (exist('WC_layered','var') && ~isempty(WC_layered)) * numel(WC_layered);
    nFinal = exist('Final','var') * numel(Final);
    if ~exist('Coverage','var'), Coverage = NaN; end  % 防止未定义
    
      % 首次写入加表头，否则追加
    new_file = ~exist(summary_csv,'file');
    
    % 确保目录存在
    [pth,~,~] = fileparts(summary_csv);
    if ~isempty(pth) && ~exist(pth,'dir'), mkdir(pth); end
    
    % 选择写入模式
    mode = 'wt';            % write text
    if ~new_file
        mode = 'at';        % append text
    end
    
    fid = fopen(summary_csv, mode);
    assert(fid ~= -1, '无法打开文件：%s', summary_csv);
    
    if new_file
        fprintf(fid, 'SampleID,HC,MC,LC,WC,EO,Final,AGR_orig,EOR_orig,AGR_layered,EOR_layered,Coverage\n');
    end
    
    fprintf(fid, '%d,%d,%d,%d,%d,%d,%d,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
            sample_id, nHC, nMC, nLC, nWC, nEO, nFinal, ...
            AGR_orig, EOR_orig, AGR_final, EOR_final, Coverage);
    
    fclose(fid);

    
    % 命令行一行汇总
    fprintf('Summary #%d | HC=%d MC=%d LC=%d WC=%d EO=%d | Final=%d | AGR_o=%.3f EOR_o=%.3f | AGR_l=%.3f EOR_l=%.3f | Cov=%.3f\n', ...
            sample_id, nHC, nMC, nLC, nWC, nEO, nFinal, ...
            AGR_orig, EOR_orig, AGR_final, EOR_final, Coverage);

    %% ---------------- 9. 可视化（只画选定样本） ----------------
    if ismember(sample_id, plot_samples)
        plot_entropy_triple(signal_denoised, fs, time, ...
            tf_entropy_norm, entropy_rate, time_axis, time_rate, ...
            HC_layered, LC_layered, EO_layered, rate_thresh, ...
            MC_layered, WC_layered, t_weak);
        saveas(gcf, fullfile(output_folder, sprintf('triple_sample_%d.png', sample_id)));
    end
end
%     %% Figure 2: AGR/EOR 分布
%     figure('Color','w','Position',[100 100 600 400]);
%     boxplot([AGR_all(:), EOR_all(:)], 'Labels', {'AGR','EOR'});
%     ylabel('Score');
%     title('Figure 2. AGR/EOR 分布 (200个样本)');
%     saveas(gcf, fullfile(output_folder, 'Figure2_AGR_EOR.png'));
%     figure;
% %     subplot(1,2,1); histogram(AGR_all,10); title('AGR 分布');
% %     subplot(1,2,2); histogram(EOR_all,10); title('EOR 分布');
% 
% 
%     %% Figure 3: Δt 分布 (改进版，概率密度风格)
%     figure('Color','w','Position',[100 100 600 400]);
%     bin_width = 0.01; % 每个柱子 0.01 ms
%     % 原始方案直方图
%     histogram(dt_all_orig, 'BinWidth',0.05,  'FaceColor',[0.6 0.6 0.6], 'FaceAlpha',0.5,'Normalization','pdf');
% % ⭐ 归一化为概率密度
%     hold on;
%     % 分层方案直方图
%    histogram(dt_all_layered, 'BinWidth',0.05, 'FaceColor',[0.2 0.6 0.6], 'FaceAlpha',0.6, 'Normalization','pdf');
% 
%     % 可选：叠加核密度估计曲线 (让分布更平滑)
% %     [f_orig,xi_orig] = ksdensity(dt_all_orig);
% %     plot(xi_orig,f_orig,'--','Color',[0.3 0.3 0.3],'LineWidth',1.5);
% %     
% %     [f_layered,xi_layered] = ksdensity(dt_all_layered);
% %     plot(xi_layered,f_layered,'-','Color',[0 0.5 0.5],'LineWidth',1.8);
%     
%     xlabel('\Delta t (ms)'); ylabel('Density');
%     legend({'Original','Layered'},'Location','northeast');
%     title('Figure 3. \Delta t Distribution (Entropy - Envelope, 200 samples)');
%     grid on;
%     
%     saveas(gcf, fullfile(output_folder, 'Figure3_Dt_Distribution_density.png'));
% 
%     %% Figure 4: Weak 事件质量分布
% 
%     figure('Color','w','Position',[100 100 600 400]);
%     hold on;
%     scatter(find(weak_quality_all<0.5), weak_quality_all(weak_quality_all<0.5), ...
%         40,[0.4 0.4 0.4],'filled','MarkerFaceAlpha',0.6,'DisplayName','Low');
%     scatter(find(weak_quality_all>=0.5 & weak_quality_all<0.6), ...
%         weak_quality_all(weak_quality_all>=0.5 & weak_quality_all<0.6), ...
%         40,[0.65 0.65 0.65],'filled','MarkerFaceAlpha',0.6,'DisplayName','Mid');
%     scatter(find(weak_quality_all>=0.6), weak_quality_all(weak_quality_all>=0.6), ...
%         40,[0.2 0.6 0.6],'filled','MarkerFaceAlpha',0.6,'DisplayName','High');
%     
%     yline(0.5,'--r','LineWidth',1.2);
%     yline(0.6,'--g','LineWidth',1.2);
%     
%     ylim([0.3 1.2]);
%     xlabel('Weak Event Index (all samples)'); ylabel('Quality Score');
%     title('Figure 4. Weak 事件质量分布 (200个样本)');
%     legend('Location','best'); grid on;
% 
%     figure('Color','w','Position',[100 100 400 400]);
%     boxchart(ones(size(weak_quality_all)), weak_quality_all,'MarkerStyle','.','BoxFaceColor',[0.2 0.6 0.9]);
%     yline(0.5,'--r'); yline(0.6,'--g');
%     ylim([0.3 1.2]);
%     ylabel('Quality Score');
%     title('Figure 4. Weak 事件质量分布 (200个样本)');
%     grid on;

    %% Figure 5: 事件类别分布变化
%     figure('Color','w','Position',[100 100 600 400]);
%     bar([counts_orig; counts_layered],'stacked'); grid on;
%     set(gca,'XTickLabel',{'原始方案','分层方案'});
%     ylabel('Total Count (200 samples)');
%     legend({'HC','LC','MC','WC','EO'},'Location','northoutside','Orientation','horizontal');
%     title('Figure 5. 事件类别分布变化 (200个样本累计)');
%     saveas(gcf, fullfile(output_folder, 'Figure5_Event_Category.png'));
%% Figure 2: AGR / EOR 分布 (增强版)
    figure('Color','w','Position',[100 100 600 400]);
    hold on;
    boxplot([AGR_all(:), EOR_all(:)], 'Labels', {'AGR','EOR'});
    scatter(ones(size(AGR_all)), AGR_all, 20, [0.3 0.6 0.8],'filled','MarkerFaceAlpha',0.4);
    scatter(2*ones(size(EOR_all)), EOR_all, 20, [0.8 0.4 0.4],'filled','MarkerFaceAlpha',0.4);
    ylabel('Score');
    title('Figure 2. AGR/EOR 分布 (200个样本)');
    grid on;
    saveas(gcf, fullfile(output_folder, 'Figure2_AGR_EOR.png'));

    %% Figure 3: Δt 分布 (用 boxchart + scatter 替代 violinplot)
    
    figure('Color','w','Position',[100 100 600 400]);
    
    % 绘制 boxchart
    group = [repmat({'Original'}, numel(dt_all_orig), 1); repmat({'Layered'}, numel(dt_all_layered), 1)];
    values = [dt_all_orig(:); dt_all_layered(:)];
    boxchart(categorical(group), values, 'BoxFaceColor',[0.6 0.6 0.9],'WhiskerLineColor','k');
    
    hold on;
    
    % 添加散点（模拟 violin 效果）
    g1 = randn(size(dt_all_orig))*0.05 + 1; % Original
    g2 = randn(size(dt_all_layered))*0.05 + 2; % Layered
    scatter(g1, dt_all_orig, 20, 'r', 'filled', 'MarkerFaceAlpha',0.3);
    scatter(g2, dt_all_layered, 20, 'b', 'filled', 'MarkerFaceAlpha',0.3);
    
    ylabel('\Delta t (ms)');
    title('Figure 3. 熵-包络 \Delta t 分布 (200个样本)');
    grid on;
    
    % 均值 ± 标准差
    mean_orig = mean(dt_all_orig); std_orig = std(dt_all_orig);
    mean_layered = mean(dt_all_layered); std_layered = std(dt_all_layered);
    
    text(1, mean_orig, sprintf('%.3f ± %.3f', mean_orig,std_orig), ...
        'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    
    text(2, mean_layered, sprintf('%.3f ± %.3f', mean_layered,std_layered), ...
        'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',10,'FontWeight','bold');
    
    saveas(gcf, fullfile(output_folder, 'Figure3_Dt_BoxScatter.png'));

    %% Figure 4: Weak 事件质量分布 (增强版)
    figure('Color','w','Position',[100 100 700 400]);
    
    subplot(1,2,1); % 散点分级
    hold on;
    scatter(find(weak_quality_all<0.5), weak_quality_all(weak_quality_all<0.5), ...
        36,'r','filled','MarkerFaceAlpha',0.6,'DisplayName','Low (<0.5)');
    scatter(find(weak_quality_all>=0.5 & weak_quality_all<0.6), ...
        weak_quality_all(weak_quality_all>=0.5 & weak_quality_all<0.6), ...
        36,'y','filled','MarkerFaceAlpha',0.6,'DisplayName','Mid (0.5-0.6)');
    scatter(find(weak_quality_all>=0.6), weak_quality_all(weak_quality_all>=0.6), ...
        36,'g','filled','MarkerFaceAlpha',0.6,'DisplayName','High (>=0.6)');
    yline(0.5,'--r'); yline(0.6,'--g');
    ylim([0.3 1.2]); xlabel('Weak Event Index'); ylabel('Quality Score');
    title('(a) Weak 事件分级散点'); legend('Location','best'); grid on;
    
    subplot(1,2,2); % 箱线图/分布
    boxplot(weak_quality_all, 'Labels', {'Weak Quality'});
    ylabel('Score'); ylim([0.3 1.2]);
    title('(b) 分布统计'); grid on;
    
    saveas(gcf, fullfile(output_folder, 'Figure4_Weak_Quality.png'));
    
%% Figure 5: 事件类别分布变化 (堆叠柱状图 + 百分比)
figure('Color','w','Position',[100 100 650 400]);

counts_matrix = [counts_orig; counts_layered]; % [2 x 5]

bar_data = counts_matrix' ./ sum(counts_matrix,2)'; % 转成百分比
hb = bar(bar_data,'stacked'); grid on;

set(gca,'XTickLabel',{'HC','LC','MC','WC','EO'});
ylabel('Proportion (%)');
legend({'原始方案','分层方案'},'Location','northoutside','Orientation','horizontal');
title('Figure 5. 事件类别比例变化 (200个样本累计)');

% 在柱子上标注百分比
for i = 1:size(bar_data,1)
    for j = 1:size(bar_data,2)
        if bar_data(i,j) > 0.05
            text(i, sum(bar_data(i,1:j))-bar_data(i,j)/2, ...
                sprintf('%.1f%%',bar_data(i,j)*100), ...
                'HorizontalAlignment','center','FontSize',8,'Color','w');
        end
    end
end

saveas(gcf, fullfile(output_folder, 'Figure5_Event_Category_Stacked.png'));
%% Figure 6: EOR Comparison (Original vs Layered)
    figure('Color','w','Position',[100 100 800 400]);
    
    plot(1:num_files, EOR_orig_all, '-o', ...
        'Color',[0.6 0.6 0.6],'LineWidth',1.2,'MarkerFaceColor',[0.6 0.6 0.6], ...
        'DisplayName','Original');
    hold on;
    
    plot(1:num_files, EOR_all, '-o', ...
        'Color',[0.2 0.6 0.9],'LineWidth',1.5,'MarkerFaceColor',[0.2 0.6 0.9], ...
        'DisplayName','Layered');
    
    xlabel('Sample ID');
    ylabel('EOR');
    legend('Location','best');
    title('Figure 6. EOR Comparison (200 samples)');
    grid on;
    
    saveas(gcf, fullfile(output_folder, 'Figure6_EOR_Comparison.png'));
    
    
    %% Figure 7: AGR Comparison (Original vs Layered)
    figure('Color','w','Position',[100 100 800 400]);
    
    plot(1:num_files, AGR_orig_all, '-o', ...
        'Color',[0.6 0.6 0.6],'LineWidth',1.2,'MarkerFaceColor',[0.6 0.6 0.6], ...
        'DisplayName','Original');
    hold on;
    
    plot(1:num_files, AGR_all, '-o', ...
        'Color',[0.2 0.6 0.9],'LineWidth',1.5,'MarkerFaceColor',[0.2 0.6 0.9], ...
        'DisplayName','Layered');
    
    xlabel('Sample ID');
    ylabel('AGR');
    legend('Location','best');
    title('Figure 7. AGR Comparison (200 samples)');
    grid on;
    
    saveas(gcf, fullfile(output_folder, 'Figure7_AGR_Comparison.png'));


    toc
    %% ---------------- 辅助函数区 ----------------
    % (A) prctile_clip       - 百分位裁剪
    % (B) match_pseudo_labels- 匹配熵事件与包络事件，生成 HC/LC/EO
    % (C) stalta_detect      - STA/LTA 包络候选检测
    % (D) eval_vs_env        - 熵检测 vs 包络评估 (Precision / Recall / F1)
    % (E) plot_entropy_triple- 三联可视化 (信号+包络 / 熵 / 熵率)
    % ---------- (A) 百分位裁剪 ----------
    function y = prctile_clip(x, lo, hi)
    % 将数据裁剪到 [lo, hi] 百分位之间
        plo = prctile(x, lo);
        phi = prctile(x, hi);
        y   = min(max(x, plo), phi);
    end
    
    
    % ---------- (B) 匹配伪标签 ----------
    function [HC, LC, EO, idx_HC_ent, idx_HC_env, idx_LC_ent, idx_EO_env, dt_signed] = ...
        match_pseudo_labels(t_entropy, t_env, tol_ms)
    % 匹配熵检测事件与包络检测事件（最近邻 + 一对一 + 容差）
    %
    % 输入：
    %   t_entropy : 熵检测事件时间（ms）
    %   t_env     : 包络检测事件时间（ms）
    %   tol_ms    : 容差 (ms)
    %
    % 输出：
    %   HC, LC, EO     : 高置信 / 低置信 / 仅包络
    %   idx_HC_ent/env : HC 对应两侧索引
    %   idx_LC_ent     : LC 对应熵索引
    %   idx_EO_env     : EO 对应包络索引
    %   dt_signed      : 熵-包络的时间差 (ms)
    
        t_entropy = t_entropy(:);
        t_env     = t_env(:);
        used_env  = false(numel(t_env),1);
    
        HC=[]; LC=[]; EO=[];
        idx_HC_ent=[]; idx_HC_env=[]; idx_LC_ent=[];
        dt_signed=[];
    
        for i = 1:numel(t_entropy)
            if isempty(t_env)
                LC(end+1,1)       = t_entropy(i);
                idx_LC_ent(end+1) = i;
                continue;
            end
    
            % 最近邻匹配
            [dmin,j] = min(abs(t_env - t_entropy(i)));
            if dmin <= tol_ms && ~used_env(j)
                % 命中 → HC
                HC(end+1,1)         = t_entropy(i);
                idx_HC_ent(end+1,1) = i;
                idx_HC_env(end+1,1) = j;
                dt_signed(end+1,1)  = (t_entropy(i) - t_env(j));
                used_env(j) = true;
            else
                % 未命中 → LC
                LC(end+1,1)       = t_entropy(i);
                idx_LC_ent(end+1) = i;
            end
        end
    
        % 未匹配 env → EO
        idx_EO_env = find(~used_env);
        EO = t_env(idx_EO_env);
    end
    
    
    % ---------- (C) STA/LTA 包络检测 ----------
    function t_ms = stalta_detect(sig_denoised, fs, time, time_start, time_end)
    % 生成包络候选事件（毫秒级）
    %
    % 输入：
    %   sig_denoised : 已去噪信号
    %   fs           : 采样率
    %   time         : 秒级时间轴
    %   time_start/end : 输出时段范围 (s)
    %
    % 输出：
    %   t_ms : 候选事件时间 (ms)
    
        % 参数
        sta_ms = 0.10;       % STA窗口
        lta_ms = 2.00;       % LTA窗口
        env_smooth_ms = 0.02;% 包络平滑
    
        sta_pts = max(1, round(sta_ms*1e-3*fs));
        lta_pts = max(sta_pts+1, round(lta_ms*1e-3*fs));
    
        env = abs(hilbert(sig_denoised));
        env = movmean(env, max(1, round(env_smooth_ms*1e-3*fs)));
    
        sta = movmean(env, sta_pts);
        lta = movmedian(env, lta_pts);
        ratio = sta ./ (lta + eps);
    
        % 阈值 = 高分位数 (自适应)
        ratio_bg = ratio(isfinite(ratio));
        thr = quantile(ratio_bg, 0.9979);
    
        above = ratio > thr;
        d = diff([0; above(:); 0]);
        seg_st = find(d==1); seg_ed = find(d==-1)-1;
    
        % 去毛刺
        min_dur_ms  = 0.04;
        min_dur_pts = max(1, round(min_dur_ms*1e-3*fs));
    
        t_list = [];
        for k = 1:numel(seg_st)
            seg = seg_st(k):seg_ed(k);
            if numel(seg) < min_dur_pts, continue; end
            [~, rel] = max(ratio(seg));
            idx = seg(1) + rel - 1;
            t_list(end+1) = time(idx)*1000; %#ok<AGROW>
        end
    
        % 合并过近事件
        min_sep_ms = 0.30;
        t_list = sort(t_list(:));
        keep = [];
        for i = 1:numel(t_list)
            if isempty(keep) || (t_list(i)-keep(end) >= min_sep_ms)
                keep(end+1) = t_list(i); %#ok<AGROW>
            end
        end
    
        t_ms = keep(:);
        if nargin >= 5
            mask = (t_ms >= time_start*1000) & (t_ms <= time_end*1000);
            t_ms = t_ms(mask);
        end
    end
    
    
    % ---------- (D) 熵检测 vs 包络评估 ----------
    function [P,R,F1,TP_times,FP_times,FN_times] = eval_vs_env(D_times_ms, Env_ms, tol_ms)
    % 将检测事件 D_times_ms 与参考 Env_ms 对比
    %
    % 输出：
    %   P,R,F1   : Precision / Recall / F1
    %   TP/FP/FN : 各类时间戳
    
        D = D_times_ms(:)'; Env = Env_ms(:)';
        used = false(size(Env)); TP=[]; FP=[];
    
        for i = 1:numel(D)
            [dmin,j] = min(abs(Env-D(i)));
            if dmin<=tol_ms && ~used(j)
                TP(end+1) = D(i); used(j) = true;
            else
                FP(end+1) = D(i);
            end
        end
    
        FN = Env(~used);
    
        P = numel(TP)/max(1,numel(TP)+numel(FP));
        R = numel(TP)/max(1,numel(TP)+numel(FN));
        F1= 2*P*R/max(1,P+R);
    
        TP_times=TP; FP_times=FP; FN_times=FN;
    end
    
    
    % ---------- (E) 三联可视化 ----------
    function plot_entropy_triple(sig_denoised, fs, time_s, entropy, rate,  time_axis, time_rate, HC, LC, EO, rate_thresh, MC, WC, TWEAK)
    % 可视化三联图：
    %   (a) 时域信号+包络+事件
    %   (b) 熵
    %   (c) 熵率
    %
    % 可选输入：
    %   MC（橙色◆）、WC（青色■）、TWEAK（灰色▲）
    
        if nargin < 11 || isempty(rate_thresh)
            m = median(rate(:));
            madv = 1.4826*median(abs(rate(:)-m));
            rate_thresh = m - 2*madv;
        end
        if nargin < 12 || isempty(MC),   MC=[];   end
        if nargin < 13 || isempty(WC),   WC=[];   end
        if nargin < 14 || isempty(TWEAK),TWEAK=[];end
    
        % 统一 ms 轴
        t_sig_ms  = time_s(:)*1000;
        t_ent_ms  = time_axis(:)*1000;
        t_rate_ms = time_rate(:)*1000;
    
        % 包络
        env = abs(hilbert(sig_denoised));
        env = movmean(env, max(1, round(0.02e-3*fs)));
    
        % 最近邻采样点
        get_y = @(t_ms, x_ms, y) interp1(x_ms, y, t_ms, 'nearest','extrap');
    
        figure('Color','w','Position',[100 100 1200 780]);
    
        % ===== (a) 时域 =====
        subplot(3,1,1); hold on;
        plot(t_sig_ms, sig_denoised,'b','DisplayName','Signal');
        plot(t_sig_ms, env,'k','LineWidth',1.1,'DisplayName','Envelope');
        scatter(HC,get_y(HC,t_sig_ms,env),36,'g','filled','DisplayName','HC');
        scatter(LC,get_y(LC,t_sig_ms,env),36,'r','filled','DisplayName','LC');
        scatter(EO,get_y(EO,t_sig_ms,env),36,'m','filled','DisplayName','EO');
        scatter(MC,get_y(MC,t_sig_ms,env),48,'d','filled','MarkerFaceColor',[1 .6 0],'DisplayName','MC');
        scatter(WC,get_y(WC,t_sig_ms,env),42,'s','filled','MarkerFaceColor',[0 .7 .9],'DisplayName','WC');
        scatter(TWEAK,get_y(TWEAK,t_sig_ms,env),36,'^','filled','MarkerFaceColor',[.6 .6 .6],'DisplayName','Weak(all)');
        title('(a) Time-domain Envelope & Events'); grid on;
        xlabel('Time (ms)'); ylabel('Amplitude');
        legend('Location','northeast');
    
        % ===== (b) 熵 =====
        subplot(3,1,2); hold on;
        plot(t_ent_ms, entropy,'b','LineWidth',1.2,'DisplayName','Entropy');
        scatter(HC,get_y(HC,t_ent_ms,entropy),40,'g','filled','DisplayName','HC');
        scatter(LC,get_y(LC,t_ent_ms,entropy),40,'r','filled','DisplayName','LC');
        scatter(EO,get_y(EO,t_ent_ms,entropy),40,'m','filled','DisplayName','EO');
        scatter(MC,get_y(MC,t_ent_ms,entropy),48,'d','filled','MarkerFaceColor',[1 .6 0],'DisplayName','MC');
        scatter(WC,get_y(WC,t_ent_ms,entropy),42,'s','filled','MarkerFaceColor',[0 .7 .9],'DisplayName','WC');
        scatter(TWEAK,get_y(TWEAK,t_ent_ms,entropy),36,'^','filled','MarkerFaceColor',[.6 .6 .6],'DisplayName','Weak(all)');
        title('(b) Normalized TF Entropy'); grid on;
        xlabel('Time (ms)'); ylabel('Entropy (norm)');
        legend('Location','northeast');
    
        % ===== (c) 熵率 =====
        subplot(3,1,3); hold on;
        plot(t_rate_ms, rate,'k','LineWidth',1.2,'DisplayName','Entropy Rate');
        yline(rate_thresh,'r--','LineWidth',1.4,'DisplayName','Threshold');
        scatter(HC,get_y(HC,t_rate_ms,rate),40,'g','filled','DisplayName','HC');
        scatter(LC,get_y(LC,t_rate_ms,rate),40,'r','filled','DisplayName','LC');
        scatter(EO,get_y(EO,t_rate_ms,rate),40,'m','filled','DisplayName','EO');
        scatter(MC,get_y(MC,t_rate_ms,rate),48,'d','filled','MarkerFaceColor',[1 .6 0],'DisplayName','MC');
        scatter(WC,get_y(WC,t_rate_ms,rate),42,'s','filled','MarkerFaceColor',[0 .7 .9],'DisplayName','WC');
        scatter(TWEAK,get_y(TWEAK,t_rate_ms,rate),36,'^','filled','MarkerFaceColor',[.6 .6 .6],'DisplayName','Weak(all)');
        title('(c) Entropy Rate'); grid on;
        xlabel('Time (ms)'); ylabel('Entropy Rate');
        legend('Location','northeast');
    end

%% **小波去噪函数**
function denoised_signal = wavelet_denoise(signal, wavelet, level, threshold_factor)
    % 小波分解
    [C, L] = wavedec(signal, level, wavelet);
    
    denoised_coeffs = C;

    % 处理每个细节系数
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