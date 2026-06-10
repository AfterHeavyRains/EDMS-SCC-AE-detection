%% =============================================================
% AE Event Detection (Method 1: Multi-step Composite Detection)
% with Pseudo Labeling (Layered) & Visualization   [4mm版]
% Author: (you)    Date: 2025-09-26
% =============================================================
%% 4mm的单个测试版，从Entropy_Detection_updated修改来的
clc; clear;

%% ---------------- 0) 路径与采样率 ----------------
file_path     = 'E:/Data/SCC_tensile/20250109_2000kHz_4mmQ235H_1000N_1/Data100.txt';
output_folder = 'E:/Code/SCC_matlab_code/results';
if ~exist(output_folder,'dir'), mkdir(output_folder); end
fs = 2e6;

% 画图/导出区间（秒）
time_start = 0e-3;
time_end   = 200e-3;
signal = load(file_path);
time   = (0:numel(signal)-1)/fs;

signal_denoised = wavelet_denoise(bandpass_filter(signal,1e5,7e5,fs,6),'sym4',5,0.3);

%% ---------------- 1) 滑窗熵 + 熵率 ----------------
% 4mm 常用 200us/100us；若想更敏：120us/60us
win_len = 200e-6; step = 100e-6;
[tf_entropy, time_axis] = compute_entropy_wsst(signal_denoised, fs, win_len, step);
[tf_entropy_norm, entropy_rate, time_rate] = post_entropy_rate(tf_entropy, time_axis);

% 实用工具（毫秒 <-> 索引）
rate_dt  = (time_axis(2)-time_axis(1));     % s
ms2rate  = @(ms) max(1, round((ms*1e-3)/rate_dt));
ms2samp  = @(ms) max(1, round(ms*1e-3*fs));

%% ---------------- 2) 复合检测（主通道） ----------------
% 用分位替代固定0.38，避免样本间归一化差异导致候选泛滥
qE = quantile(tf_entropy_norm, 0.1);          % 0.10 可试
params_main.entropy_thresh = qE;
params_main.wR_ms          = 1.0;    % 熵率邻域 ±wR（ms）
params_main.min_interval_ms= 8.0;    %8  熵候选最小间隔（ms）

% 能量/突发门（全局 + 动态）
params_main.energy_global_mult = 1.25;% 1.25 
params_main.burst_thresh       = 3;% 3.0
params_main.alpha_dyn          = 6;   % 6.0动态背景倍率
params_main.win_energy_ms      = 1.0;   % 全局能量窗
params_main.win_small_ms       = 0.08;  % 动态局部窗
params_main.win_bg_ms          = 2.5;   % 动态背景窗

[t_entropy] = detect_main_channel(signal_denoised, fs, time, time_axis, ...
    tf_entropy_norm, entropy_rate, params_main);

%% ---------------- 3) 包络候选 (STA/LTA, 自适应分位阈) ----------------
params_env.sta_ms         = 0.12;
params_env.lta_ms         = 2.00;
params_env.env_smooth_ms  = 0.02;
params_env.quantile_thr   = 0.9972;%0.9972
params_env.min_dur_ms     = 0.052;%0.052
params_env.min_sep_ms     = 0.45;%0.45

stalta_evt_ms = stalta_detect(signal_denoised, fs, time, time_start, time_end, params_env);

%% ---------------- 4) 原始匹配 → HC/LC/EO ----------------
tol_ms = 0.55;
[HC, LC, EO, ~, ~, ~, ~, dt_signed_orig] = match_pseudo_labels(t_entropy, stalta_evt_ms, tol_ms);

AGR_orig = numel(HC)/max(1, (numel(HC)+numel(LC)));
EOR_orig = numel(EO)/max(1, (numel(HC)+numel(EO)));

% === 轻诊断 ===
fprintf('\n[Diag] ENV触发数: %d\n', numel(stalta_evt_ms));
fprintf('[Diag] 主通道候选数(=HC+LC): %d\n', numel(HC)+numel(LC));
fprintf('[Diag] HC=%d, LC=%d, EO=%d, AGR=%.3f\n', numel(HC), numel(LC), numel(EO), AGR_orig);

% 看 Δt 是否有系统偏差（群时延/窗中心误差）
if ~isempty(dt_signed_orig)
    med_dt = median(dt_signed_orig); mad_dt = 1.4826*median(abs(dt_signed_orig - med_dt));
    fprintf('[Diag] Δt中位数=%.3f ms, MAD=%.3f ms\n', med_dt, mad_dt);
end
% === ENV 对齐补偿（若存在系统性偏移）===
if ~isempty(dt_signed_orig)
    med_dt = median(dt_signed_orig);
    if abs(med_dt) > 0.25       % 阈值可取 0.25~0.5 ms
        stalta_evt_ms = stalta_evt_ms + med_dt;   % 对齐 ENV
        [HC, LC, EO, ~, ~, ~, ~, dt_signed_orig] = ...
            match_pseudo_labels(t_entropy, stalta_evt_ms, tol_ms);
        AGR_orig = numel(HC)/max(1, (numel(HC)+numel(LC)));
        EOR_orig = numel(EO)/max(1, (numel(HC)+numel(EO)));
        fprintf('[Align] Δt中位数=%.3f ms，已对齐：AGR=%.3f, EOR=%.3f\n', med_dt, AGR_orig, EOR_orig);
    end
end

%% ---------------- 5) 微弱事件通道（Weak） ----------------
params_weak.entropy_lo       = params_main.entropy_thresh;       % 不与主通道重叠的下界
params_weak.entropy_hi       = min(0.9, params_main.entropy_thresh + 0.12);        % 首轮弱通道上界
params_weak.rate_k_weak      = 1.0;         % 熵率阈（更松）
params_weak.wR_ms            = 1.5;
params_weak.burst_thresh     = 2.0;
params_weak.alpha_dyn        = 3.5;
params_weak.win_small_ms     = 0.10;
params_weak.win_bg_ms        = 2.0;
params_weak.quality_cut      = 0.35;        % 综合质量阈
params_weak.min_sep_ms       = 3.0;         % 弱通道 NMS

[t_weak, weak_scores] = detect_weak_channel( ...
    signal_denoised, fs, time_axis, center_index(time, time_axis), ...
    tf_entropy_norm, entropy_rate, params_weak);

%% ---------------- 6) 分层伪标签（HC/LC + MC/WC + EO_layered） ----------------
[HC_l, LC_l, MC_l, WC_l, EO_l, dt_signed_layered] = layered_labeling( ...
    t_entropy, stalta_evt_ms, t_weak, weak_scores, tol_ms);

% 发布清单：Final = HC + MC（高置信）
Final = sort([HC_l; MC_l]);

% 指标（Layered）
high_conf = numel(HC_l) + numel(MC_l);
tot_ent   = high_conf + numel(LC_l);
Coverage  = numel(Final)/max(1, (numel(Final)+numel(EO_l)));
EOR_layer = numel(EO_l)/max(1, (high_conf+numel(EO_l)));
AGR_layer = high_conf/max(1, tot_ent);

%% ---------------- 7) 导出 ----------------
writematrix(Final,      fullfile(output_folder,'AE_final_times_ms.txt'));
writematrix(HC_l,       fullfile(output_folder,'AE_eval_HC.txt'));
writematrix(MC_l,       fullfile(output_folder,'AE_eval_MC.txt'));
writematrix(LC_l,       fullfile(output_folder,'AE_eval_LC.txt'));
writematrix(WC_l,       fullfile(output_folder,'AE_eval_WC.txt'));  % 仅候选
writematrix(EO_l,       fullfile(output_folder,'AE_eval_EO.txt'));
if ~isempty(weak_scores)
    weak_table = table(t_weak, weak_scores, 'VariableNames', {'Time_ms','QualityScore'});
    writetable(weak_table, fullfile(output_folder,'weak_events_quality.csv'));
end

%% ---------------- 8) 打印汇总 ----------------
fprintf('\n=== 方案对比 ===\n');
fprintf(' 原始: HC=%d, LC=%d, EO=%d | AGR=%.3f, EOR=%.3f\n', numel(HC), numel(LC), numel(EO), AGR_orig, EOR_orig);
fprintf(' 分层: HC=%d, MC=%d, LC=%d, EO=%d | AGR=%.3f, EOR=%.3f, Cov=%.3f\n', ...
    numel(HC_l), numel(MC_l), numel(LC_l), numel(EO_l), AGR_layer, EOR_layer, Coverage);

%% ---------------- 9) 可视化 ----------------
% 9.1 三联图（信号+包络 / 熵 / 熵率）
rate_thresh_vis = mean(entropy_rate) - 2.0*std(entropy_rate);  % 仅显示用
plot_entropy_triple(signal_denoised, fs, time, ...
    tf_entropy_norm, entropy_rate, time_axis, time_rate, ...
    HC_l, LC_l, EO_l, rate_thresh_vis, MC_l, WC_l, t_weak);

% 9.2 AGR/EOR 条形图（使用本次真实数值）
figure('Color','w','Position',[100 100 560 380]);
bar_data = [AGR_orig EOR_orig; AGR_layer EOR_layer];
bar(bar_data,'grouped'); grid on; ylim([0 1]);
set(gca,'XTickLabel',{'原始','分层'});
legend({'AGR','EOR'},'Location','northoutside','Orientation','horizontal');
title('AGR / EOR 对比');

% 9.3 Δt 直方图（原始 vs 分层）
%   - 原始：用 HC 的 dt_signed_orig
%   - 分层：对 Final 重新与 ENV 匹配，得到 dt_signed_layered
if ~isempty(dt_signed_orig) || ~isempty(dt_signed_layered)
    figure('Color','w','Position',[100 100 560 380]);
    hold on; grid on;
    if ~isempty(dt_signed_orig)
        histogram(dt_signed_orig, 'BinWidth',0.05,'FaceAlpha',0.6);
    end
    if ~isempty(dt_signed_layered)
        histogram(dt_signed_layered, 'BinWidth',0.05,'FaceAlpha',0.6);
    end
    xlabel('\Delta t = t_{entropy}-t_{env} (ms)'); ylabel('Count');
    legend({'原始 HC','分层 Final'},'Location','best');
    title('\Delta t 分布对比');
end

fprintf('完成。\n');
%% ====================== 参数轨迹日志表格 ======================
logfile = fullfile(output_folder, 'param_log.csv');

% 如果不存在，就写入表头
if ~exist(logfile, 'file')
    fid = fopen(logfile, 'w');
    fprintf(fid, ['Timestamp,qE,thrR_coeff,energy_global_mult,burst_thresh,alpha_dyn,',...
                  'quantile_thr,min_dur_ms,min_sep_ms,tol_ms,',...
                  'HC,LC,EO,AGR,EOR,Cov\n']);
    fclose(fid);
end

% 写入当前运行结果
fid = fopen(logfile, 'a');
fprintf(fid, '%s,%.4f,%.2f,%.2f,%.2f,%.2f,%.6f,%.3f,%.3f,%.3f,%d,%d,%d,%.3f,%.3f,%.3f\n', ...
    datestr(now, 'yyyy-mm-dd HH:MM:SS'), ...
    qE, 2.3, ...  % ← 这里填你实际的 thrR 系数（2.3 / 2.4 / 2.6）
    params_main.energy_global_mult, ...
    params_main.burst_thresh, ...
    params_main.alpha_dyn, ...
    params_env.quantile_thr, ...
    params_env.min_dur_ms, ...
    params_env.min_sep_ms, ...
    tol_ms, ...
    numel(HC), numel(LC), numel(EO), ...
    AGR_orig, EOR_orig, Coverage);
fclose(fid);

fprintf('\n[Log] 参数与结果已记录至 param_log.csv\n');

%% ========================== 局部函数 ==========================
% 放在脚本末尾，避免“脚本中的函数定义必须在结尾”的错误。

function [tf_entropy, time_axis] = compute_entropy_wsst(sig, fs, win_len, step)
    win_pts  = floor(win_len*fs);
    step_pts = floor(step*fs);
    start_idx  = 1:step_pts:(numel(sig)-win_pts+1);
    center_idx = start_idx + floor(win_pts/2);
    time_axis  = (center_idx-1)/fs;

    W = numel(start_idx);
    q = 1.6; tf_entropy = zeros(1,W);
    for i = 1:W
        seg = sig(start_idx(i):start_idx(i)+win_pts-1);
        [sst, ~] = wsst(seg, fs, 'bump');
        P = abs(sst).^2; P = P ./ max(sum(P,1), eps);
        tf_entropy(i) = (1 - sum(P(:).^q)) / (q - 1);
    end
end

function [ent_norm, ent_rate, time_rate] = post_entropy_rate(tf_entropy, time_axis)
    % Min-Max 归一化（与你当前稳定版本保持一致）
    ent_norm = (tf_entropy - min(tf_entropy)) / max(eps, (max(tf_entropy) - min(tf_entropy)));
    ent_s = movmedian(ent_norm, 5);
    ent_rate = diff(ent_s) ./ diff(time_axis);
    time_rate = time_axis(2:end);
end

function idx = center_index(time_s, time_axis)
    % 把熵中心点映射回原始采样索引（用于能量窗口）
    t_cent = time_axis(:);
    idx = max(1, min(numel(time_s), round(t_cent * (numel(time_s)-1)/time_s(end)) ));
end

function t_ms = detect_main_channel(sig, fs, time, time_axis, ent_norm, ent_rate, p)
    % 1) 熵 + 熵率双门
    thrE  = p.entropy_thresh;

    mR = median(ent_rate);
    madR = 1.4826*median(abs(ent_rate - mR));
    thrR = mR - 2.6*madR;   % 2.6 可试
    thrRateMag = 1.2*madR;   % <<< 新增：负斜率幅度门
% --- 局部显著性门：收回到 2.4/2.5 ---
idx_low = find(ent_norm(2:end) < thrE);
wProm = ms2rate_local(max(1.2, p.wR_ms), time_axis);
keepProm = false(size(idx_low));
for kk = 1:numel(idx_low)
    ii = idx_low(kk);
    L  = max(2, ii - wProm);
    R  = min(numel(ent_norm)-1, ii + wProm);
    neigh   = ent_norm(L:R);
    medLoc  = median(neigh);
    madLoc  = 1.4826*median(abs(neigh - medLoc));
    keepProm(kk) = (ent_norm(ii+1) < (medLoc - 2.5*madLoc));  % 2.5
end
idx_low = idx_low(keepProm);
fprintf('[Stage] after local-saliency: %d\n', numel(idx_low));

% --- 负熵率门：只用“邻域最小”这一条（AND/OR都不要了） ---
wR = ms2rate_local(p.wR_ms, time_axis);
cand = [];
for k = 1:numel(idx_low)
    idx = idx_low(k);
    L = max(1, idx - wR);
    R = min(numel(ent_rate), idx + wR);
    if min(ent_rate(L:R)) < thrR
        cand(end+1) = idx; %#ok<AGROW>
    end
end
cand = unique(cand(:));
fprintf('[Stage] after rate-gate: %d\n', numel(cand));
if isempty(cand), t_ms = []; return; end


    if isempty(cand), t_ms = []; return; end

    % 2) 全局能量 + 突发
    window_energy = ms2samp_local(p.win_energy_ms, fs);
    energy_bg_global = movsum(sig.^2, window_energy);
    Eglob = p.energy_global_mult * median(energy_bg_global);

    center_idx = nearest_sample_index(time, time_axis);
    keep = false(size(cand));
    for k = 1:numel(cand)
        samp_idx = center_idx(cand(k));
        L = max(1, samp_idx - window_energy);
        R = min(numel(sig), samp_idx + window_energy);
        win_sig = sig(L:R);
        elocal  = sum(win_sig.^2);
        burst   = kurtosis(abs(win_sig));
        keep(k) = (elocal > Eglob) && (burst > p.burst_thresh);
    end
    cand = cand(keep);
    if isempty(cand), t_ms = []; return; end

    % 3) 动态背景能量
    ws = ms2samp_local(p.win_small_ms, fs);
    wb = ms2samp_local(p.win_bg_ms, fs);
    keep2 = false(size(cand));
    for k = 1:numel(cand)
        samp_idx = center_idx(cand(k));
        L = max(1, samp_idx - ws);
        R = min(numel(sig), samp_idx + ws);
        elocal = sum(sig(L:R).^2);

        Lbg = max(1, samp_idx - wb);
        Rbg = min(numel(sig), samp_idx + wb);
        bg_seg = sig(Lbg:Rbg);
        ebg_series = movsum(bg_seg.^2, ws);
        m = mean(ebg_series); s = std(ebg_series);
        keep2(k) = (elocal > m + p.alpha_dyn*s);
    end
    cand = cand(keep2);
    if isempty(cand), t_ms = []; return; end

    % 4) 最小间隔
    AE_times_ms = sort(time_axis(cand)*1000);
    t_ms = [];
    for i = 1:numel(AE_times_ms)
        if isempty(t_ms) || AE_times_ms(i) - t_ms(end) > p.min_interval_ms
            t_ms(end+1,1) = AE_times_ms(i); %#ok<AGROW>
        end
    end
end

function [t_weak, scores] = detect_weak_channel(sig, fs, time_axis, center_idx, ent_norm, ent_rate, p)
%     muR = mean(ent_rate); sdR = std(ent_rate);
%     thrR_w = muR - p.rate_k_weak*sdR;
    mR = median(ent_rate);
    madR = 1.4826*median(abs(ent_rate - mR));
    thrR_w = mR - 1.6*madR;   % 弱通道更松
    idx_w = find(ent_norm(2:end) < p.entropy_hi & ent_norm(2:end) >= p.entropy_lo);
    wR = ms2rate_local(p.wR_ms, time_axis);

    cand_w = [];
    for idx = idx_w'
        L=max(1,idx-wR); R=min(numel(ent_rate), idx+wR);
        if min(ent_rate(L:R)) < thrR_w, cand_w = [cand_w idx]; end %#ok<AGROW>
    end
    cand_w = unique(cand_w(:));

    t_weak = []; scores = [];
    if isempty(cand_w), return; end

    ws = ms2samp_local(p.win_small_ms, fs);
    wb = ms2samp_local(p.win_bg_ms, fs);

    % 计算质量
    for ii = 1:numel(cand_w)
        rix = cand_w(ii);
        sidx = center_idx(rix);

        L = max(1, sidx-ws); R = min(numel(sig), sidx+ws);
        win_sig = sig(L:R);
        elocal  = sum(win_sig.^2);
        burst   = kurtosis(abs(win_sig));

        Lbg = max(1, sidx-wb); Rbg = min(numel(sig), sidx+wb);
        bg_seg = sig(Lbg:Rbg);
        ebg_series = movsum(bg_seg.^2, ws);
        mean_bg = mean(ebg_series); std_bg = std(ebg_series);

        if elocal <= mean_bg + p.alpha_dyn*std_bg || burst <= p.burst_thresh
            continue;
        end

        % 6维特征
        feat = zeros(1,6);
        feat(1) = max(0, 1 - (ent_norm(rix)-p.entropy_lo)/max(1e-6,(p.entropy_hi-p.entropy_lo)));
% 替换 feat(2) 这一行
        feat(2) = min(1, abs(ent_rate(max(1,rix-1)) - mR)/(3*(madR+eps)));

        energy_ratio = elocal/(mean_bg+eps);
        feat(3) = max(0, min(1, log10(energy_ratio)/1.5));
        feat(4) = max(0, min(1, (burst-1)/2));
        if numel(win_sig) >= 8
            signal_power = max(movsum(win_sig.^2, max(1, round(numel(win_sig)/4))));
            noise_power  = median(movsum(win_sig.^2, max(1, round(numel(win_sig)/4))));
            snr = signal_power/(noise_power+eps);
            feat(5) = max(0, min(1, log10(snr)/2));
        else
            feat(5) = 0.3;
        end
        % 频域粗特征
        if numel(win_sig) >= 32
            try
                [psd,f] = pwelch(win_sig, [], [], [], fs);
                ae_band = (f>=1e5) & (f<=6e5);
                noise_band = f<5e4;
                ae_power = sum(psd(ae_band));
                noise_power_freq = sum(psd(noise_band));
                spectral_ratio = ae_power/(noise_power_freq+eps);
                feat(6) = max(0, min(1, log10(spectral_ratio)/2));
            catch
                feat(6) = 0.3;
            end
        else
            feat(6) = 0.3;
        end

        w = [0.25,0.20,0.15,0.15,0.15,0.10];
        sc = sum(w.*feat);
        if sc >= p.quality_cut
            t_weak(end+1,1) = time_axis(rix)*1000; %#ok<AGROW>
            scores(end+1,1) = sc;                 %#ok<AGROW>
        end
    end

    % NMS（最小间隔）
    if ~isempty(t_weak)
        [t_weak, ord] = sort(t_weak); scores = scores(ord);
        keep = true(size(t_weak));
        for i=2:numel(t_weak)
            if t_weak(i)-t_weak(i-1) < p.min_sep_ms
                if scores(i) > scores(i-1), keep(i-1)=false; else, keep(i)=false; end
            end
        end
        t_weak = t_weak(keep); scores = scores(keep);
    end
end

function [HC_l, LC_l, MC_l, WC_l, EO_l, dt_signed_layered] = layered_labeling(t_entropy, env_ms, t_weak, weak_scores, tol_ms)
    % 先用主通道匹配出 HC/LC，剩余 env 给弱通道配 MC
    [HC, LC, EO, idx_HC_ent, idx_HC_env] = match_pseudo_labels(t_entropy, env_ms, tol_ms);
    used_env = false(numel(env_ms),1);
    used_env(idx_HC_env) = true;
    rem_env = env_ms(~used_env);

    MC=[]; WC=[]; used2=false(numel(rem_env),1);
    for i=1:numel(t_weak)
        if isempty(rem_env), WC(end+1,1)=t_weak(i); continue; end %#ok<AGROW>
        [dmin,j] = min(abs(rem_env - t_weak(i)));
        if dmin<=tol_ms && ~used2(j)
            MC(end+1,1)=t_weak(i); used2(j)=true;
        else
            WC(end+1,1)=t_weak(i);
        end
    end
    EO_layered = rem_env(~used2);

    HC_l = HC; LC_l = LC; MC_l = MC; WC_l = WC; EO_l = EO_layered;

    % 计算分层 Δt：Final 与 ENV 重新做一次匹配
    Final = sort([HC_l; MC_l]);
    [~,~,~,~,~,~,~, dt_signed_layered] = match_pseudo_labels(Final, env_ms, tol_ms);
end

function [HC, LC, EO, idx_HC_ent, idx_HC_env, idx_LC_ent, idx_EO_env, dt_signed] = ...
    match_pseudo_labels(t_entropy, t_env, tol_ms)
    t_entropy = t_entropy(:); t_env = t_env(:);
    used = false(numel(t_env),1);
    HC=[]; LC=[]; EO=[];
    idx_HC_ent=[]; idx_HC_env=[]; idx_LC_ent=[];
    dt_signed=[];

    for i=1:numel(t_entropy)
        if isempty(t_env)
            LC(end+1,1)=t_entropy(i); idx_LC_ent(end+1,1)=i; continue;
        end
        [dmin,j] = min(abs(t_env - t_entropy(i)));
        if dmin<=tol_ms && ~used(j)
            HC(end+1,1)=t_entropy(i);
            idx_HC_ent(end+1,1)=i; idx_HC_env(end+1,1)=j;
            dt_signed(end+1,1) = t_entropy(i) - t_env(j);
            used(j)=true;
        else
            LC(end+1,1)=t_entropy(i);
            idx_LC_ent(end+1,1)=i;
        end
    end
    idx_EO_env = find(~used);
    EO = t_env(idx_EO_env);
end

function t_ms = stalta_detect(sig, fs, time, t0, t1, p)
    sta = max(1, round(p.sta_ms*1e-3*fs));
    lta = max(sta+1, round(p.lta_ms*1e-3*fs));
    env = abs(hilbert(sig));
    env = movmean(env, max(1, round(p.env_smooth_ms*1e-3*fs)));
    S = movmean(env, sta);
    L = movmedian(env, lta);
    ratio = S ./ (L + eps);
    valid = ratio(isfinite(ratio));
    if isempty(valid), t_ms = []; return; end
    thr = quantile(valid, p.quantile_thr);

    above = ratio > thr;
    d = diff([0; above(:); 0]);
    st = find(d==1); ed = find(d==-1)-1;

    min_dur_pts = max(1, round(p.min_dur_ms*1e-3*fs));
    evt = [];
    for i=1:numel(st)
        seg = st(i):ed(i);
        if numel(seg) < min_dur_pts, continue; end
        [~,ix] = max(ratio(seg));
        evt(end+1) = st(i)+ix-1; %#ok<AGROW>
    end

    % 合并过近
    min_sep_pts = max(1, round(p.min_sep_ms*1e-3*fs));
    evt = sort(evt(:));
    keep = true(size(evt));
    for i=2:numel(evt)
        if (evt(i)-evt(i-1)) < min_sep_pts
            if ratio(evt(i)) >= ratio(evt(i-1)), keep(i-1)=false; else, keep(i)=false; end
        end
    end
    evt = evt(keep);

    t_ms = time(evt)*1000;
    mask = (t_ms >= t0*1000) & (t_ms <= t1*1000);
    t_ms = t_ms(mask);
end

function center_idx = nearest_sample_index(time, time_axis)
    % 将熵中心点（time_axis）映射到原始采样索引
    center_idx = round(interp1(time, 1:numel(time), time_axis, 'nearest','extrap'));
    center_idx = min(max(center_idx,1), numel(time));
end

function n = ms2rate_local(ms, time_axis)
    dt = time_axis(2)-time_axis(1);
    n  = max(1, round((ms*1e-3)/dt));
end

function n = ms2samp_local(ms, fs)
    n = max(1, round(ms*1e-3*fs));
end

function plot_entropy_triple(sig_denoised, fs, time_s, entropy, rate, ...
    time_axis, time_rate, HC, LC, EO, rate_thresh, MC, WC, TWEAK)

    if nargin < 11 || isempty(rate_thresh)
        m = median(rate(:)); madv = 1.4826*median(abs(rate(:)-m));
        rate_thresh = m - 2*madv;
    end
    if nargin < 12 || isempty(MC),   MC=[];   end
    if nargin < 13 || isempty(WC),   WC=[];   end
    if nargin < 14 || isempty(TWEAK),TWEAK=[];end

    t_sig_ms  = time_s(:)*1000;
    t_ent_ms  = time_axis(:)*1000;
    t_rate_ms = time_rate(:)*1000;

    env = abs(hilbert(sig_denoised));
    env = movmean(env, max(1, round(0.02e-3*fs)));

    get_y = @(t_ms, x_ms, y) interp1(x_ms, y, t_ms, 'nearest','extrap');

    figure('Color','w','Position',[100 100 1200 780]);

    % (a) 时域
    subplot(3,1,1); hold on; grid on;
    plot(t_sig_ms, sig_denoised,'b','DisplayName','Signal');
    plot(t_sig_ms, env,'k','LineWidth',1.1,'DisplayName','Envelope');
    if ~isempty(HC), scatter(HC,get_y(HC,t_sig_ms,env),36,'g','filled','DisplayName','HC'); end
    if ~isempty(LC), scatter(LC,get_y(LC,t_sig_ms,env),36,'r','filled','DisplayName','LC'); end
    if ~isempty(EO), scatter(EO,get_y(EO,t_sig_ms,env),36,'m','filled','DisplayName','EO'); end
    if ~isempty(MC), scatter(MC,get_y(MC,t_sig_ms,env),48,'d','filled','MarkerFaceColor',[1 .6 0],'DisplayName','MC'); end
    if ~isempty(WC), scatter(WC,get_y(WC,t_sig_ms,env),42,'s','filled','MarkerFaceColor',[0 .7 .9],'DisplayName','WC'); end
    if ~isempty(TWEAK), scatter(TWEAK,get_y(TWEAK,t_sig_ms,env),36,'^','filled','MarkerFaceColor',[.6 .6 .6],'DisplayName','Weak(all)'); end
    title('(a) Time-domain Envelope & Events'); xlabel('Time (ms)'); ylabel('Amp'); legend('Location','northeast');

    % (b) 熵
    subplot(3,1,2); hold on; grid on;
    plot(t_ent_ms, entropy,'b','LineWidth',1.2,'DisplayName','Entropy');
    if ~isempty(HC), scatter(HC,get_y(HC,t_ent_ms,entropy),40,'g','filled','DisplayName','HC'); end
    if ~isempty(LC), scatter(LC,get_y(LC,t_ent_ms,entropy),40,'r','filled','DisplayName','LC'); end
    if ~isempty(EO), scatter(EO,get_y(EO,t_ent_ms,entropy),40,'m','filled','DisplayName','EO'); end
    if ~isempty(MC), scatter(MC,get_y(MC,t_ent_ms,entropy),48,'d','filled','MarkerFaceColor',[1 .6 0],'DisplayName','MC'); end
    if ~isempty(WC), scatter(WC,get_y(WC,t_ent_ms,entropy),42,'s','filled','MarkerFaceColor',[0 .7 .9],'DisplayName','WC'); end
    if ~isempty(TWEAK), scatter(TWEAK,get_y(TWEAK,t_ent_ms,entropy),36,'^','filled','MarkerFaceColor',[.6 .6 .6],'DisplayName','Weak(all)'); end
    title('(b) Normalized TF Entropy'); xlabel('Time (ms)'); ylabel('Entropy (norm)'); legend('Location','northeast');

    % (c) 熵率
    subplot(3,1,3); hold on; grid on;
    plot(t_rate_ms, rate,'k','LineWidth',1.2,'DisplayName','Entropy Rate');
    yline(rate_thresh,'r--','LineWidth',1.4,'DisplayName','Threshold');
    if ~isempty(HC), scatter(HC,get_y(HC,t_rate_ms,rate),40,'g','filled','DisplayName','HC'); end
    if ~isempty(LC), scatter(LC,get_y(LC,t_rate_ms,rate),40,'r','filled','DisplayName','LC'); end
    if ~isempty(EO), scatter(EO,get_y(EO,t_rate_ms,rate),40,'m','filled','DisplayName','EO'); end
    if ~isempty(MC), scatter(MC,get_y(MC,t_rate_ms,rate),48,'d','filled','MarkerFaceColor',[1 .6 0],'DisplayName','MC'); end
    if ~isempty(WC), scatter(WC,get_y(WC,t_rate_ms,rate),42,'s','filled','MarkerFaceColor',[0 .7 .9],'DisplayName','WC'); end
    if ~isempty(TWEAK), scatter(TWEAK,get_y(TWEAK,t_rate_ms,rate),36,'^','filled','MarkerFaceColor',[.6 .6 .6],'DisplayName','Weak(all)'); end
    title('(c) Entropy Rate'); xlabel('Time (ms)'); ylabel('Entropy Rate'); legend('Location','northeast');
end
%% 带通滤波函数
function filtered_signal = bandpass_filter(signal, low_cutoff, high_cutoff, fs, order)
    [b, a] = butter(order, [low_cutoff, high_cutoff] / (fs / 2), 'bandpass');
    filtered_signal = filtfilt(b, a, signal);
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