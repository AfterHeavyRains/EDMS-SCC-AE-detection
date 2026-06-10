%% ======================== Ablation Study for Table 2 ========================
function run_ablation_study()
% eval_and_save_all_AblationStudy;
% 运行消融实验（Table 2）
% Liu Guangyu 2025-10-09

poolobj = gcp('nocreate');
delete(poolobj)
parpool(6)
rootBase = 'D:\Code\Entropy_Detection_Viz_updated_all';
mkData   = @(mm) fullfile(rootBase, sprintf('result%imm',mm), 'data_denoised');
dataDirs = struct('mm2', mkData(2), 'mm4', mkData(4), 'mm6', mkData(6), 'mm8', mkData(8));
thicks   = {'mm2','mm4','mm6','mm8'};
tlabel   = struct('mm2','2mm','mm4','4mm','mm6','6mm','mm8','8mm');

outDir   = fullfile(rootBase, 'eval_exports', 'ablation');
if ~exist(outDir,'dir'), mkdir(outDir); end
baseP = struct( ...
  'use_sst', false, ...
  'use_percentile_norm', true, 'entropy_pct_lo', 1, 'entropy_pct_hi', 99, 'rate_smooth_win', 3, ...
  'use_fixed_rate', true, 'fixed_rate_thr', -1700, ...        % 稍收紧
  'tol_ms', 1.0, 'min_interval_ms', 2.6, ...
  'use_energy_gates', true, 'energy_global_mult', 1.25, ...   % A/C中间值
  'burst_thresh', 1.35, ...                                   % C中较优值
  'use_dynamic_bg', true, 'alpha', 3.5, ...                   % 提高动态门系数
  'stalta_mode','ratio','stalta_thresh',2.5,'env_smooth_ms',0.02, ... % 介于2.3~3.25之间
  'snr_thresh', 9.6, ...
  'win_energy_pts', 200, 'bg_win_pts', 12000, ...
  'use_weak_scoring', true, ...
  'weak_quality_cut', 0.45, 'weak_min_sep_ms', 3.0, 'rate_k_weak', 0.4, ...
  'min_dur_ms', 0.04, 'min_sep_ms', 0.5, ...
  'pre_denoised', true, 'stalta_quantile', 0.995);
% ====== 通用参数 ======

baseP.use_fixed_rate   = true;
baseP.fixed_rate_thr   = -1700;
baseP.tol_ms           = 1.8;
baseP.rate_k           = 2.0;
baseP.rate_k_weak      = 0.40;
baseP.weak_quality_cut = 0.45;

bestO.ent = struct('mm2',0.540,'mm4',0.540,'mm6',0.540,'mm8',0.540);
bestL.delta = 0.24;

% ====== 实验组定义 ======
% ablations = {
%  'Full (Layered)',             @(p)p;                                  % baseline full
%  '- Entropy-rate (τ_R)', @(p)remove_entropy_rate_gate(p);
%  '- Tolerance match (Δt)', @(p)setfield(p,'tol_ms',0);  % 
%  '- Weak-event quality score', @(p)remove_weak_scoring(p);
%  '- Local fine-sweep', @(p)disable_fine_sweep(p);
%  'Original baseline', @(p)restore_original_mode(p);
% 
%  'Envelope-only (Hilbert/STA-LTA)', @(p)disable_entropy_branch(p);
%  'Entropy-only (no τ_R/no Δt)', @(p)disable_envelope_branch(p);
%  };
% ====== 实验组定义 ======
ablations = {
 'Full (Layered)',             @(p)p;                                  % baseline full
 '- Entropy-rate (τ_R)',       @(p)remove_entropy_rate_gate(p);         % 去掉熵率门控
 '- Dynamic Energy Gate (α)',  @(p)remove_dynamic_energy_gate(p);       % 去掉动态能量门控
 '- Weak-event quality score', @(p)remove_weak_scoring_soft(p);         % 去掉弱事件评分（宽容版）
 '- Tolerance match (Δt)',     @(p)reduce_tolerance(p);                 % 减小时间容差
 'Original baseline',          @(p)restore_original_mode(p);            % 单层原始基线
 'Envelope-only (Hilbert/STA-LTA)', @(p)disable_entropy_branch(p);     % 保留不变
 'Entropy-only (no env)',      @(p)disable_envelope_branch_soft(p);     % 宽容禁用包络
};


metrics = struct();

% 并行准备
if isempty(gcp('nocreate')), parpool; end

for ai = 1:numel(ablations)
    name = ablations{ai,1};
    fn   = ablations{ai,2};
    name = matlab.lang.makeValidName(name);
    fprintf('\n=== [%d/%d] %s ===\n', ai, numel(ablations), name);
    for t = 1:numel(thicks)
        tk = thicks{t};
        ddir = dataDirs.(tk);
        params = baseP;
        params.entropy_thresh       = bestO.ent.(tk);
        params.entropy_thresh_weak  = params.entropy_thresh + bestL.delta;
        params = fn(params);  % 应用禁用逻辑

%         [mAGRo, mAGRl, mEORo, mEORl, mCovO, mCovL, ~] = eval_on_set(ddir, 1:200, params);
%         metrics.(name).(tk) = [mAGRl, mEORl, mCovL]; % Layered为主，简化表格
        [AGRo, AGRl, EORo, EORl, CovO, CovL, ~] = eval_on_set_full(ddir, 1:200, params);
        % mean±std 统计

        metrics.(name).(tk).AGR_mean = mean(AGRl,'omitnan');
        metrics.(name).(tk).AGR_std  = std(AGRl,'omitnan');
        metrics.(name).(tk).EOR_mean = mean(EORl,'omitnan');
        metrics.(name).(tk).EOR_std  = std(EORl,'omitnan');
        metrics.(name).(tk).Cov_mean = mean(CovL,'omitnan');
        metrics.(name).(tk).Cov_std  = std(CovL,'omitnan');

        fprintf('[%s] AGR=%.3f  EOR=%.3f  Cov=%.3f\n', ...
    tlabel.(tk), metrics.(name).(tk).AGR_mean, ...
    metrics.(name).(tk).EOR_mean, metrics.(name).(tk).Cov_mean);

    end
end

% ====== 汇总输出 Table 2 (mean ± std) ======
outFile = fullfile(outDir, 'Table2_ablation.csv');
fid = fopen(outFile,'w');

% CSV 表头
fprintf(fid, ['Setting,', ...
    '2mm_AGR_mean,2mm_AGR_std,2mm_EOR_mean,2mm_EOR_std,2mm_Cov_mean,2mm_Cov_std,', ...
    '4mm_AGR_mean,4mm_AGR_std,4mm_EOR_mean,4mm_EOR_std,4mm_Cov_mean,4mm_Cov_std,', ...
    '6mm_AGR_mean,6mm_AGR_std,6mm_EOR_mean,6mm_EOR_std,6mm_Cov_mean,6mm_Cov_std,', ...
    '8mm_AGR_mean,8mm_AGR_std,8mm_EOR_mean,8mm_EOR_std,8mm_Cov_mean,8mm_Cov_std\n']);

names = fieldnames(metrics);
for i = 1:numel(names)
    nm = names{i};
    line = nm;
    for t = 1:numel(thicks)
        M = metrics.(nm).(thicks{t});
        line = sprintf('%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f', line, ...
            M.AGR_mean, M.AGR_std, M.EOR_mean, M.EOR_std, M.Cov_mean, M.Cov_std);
    end
    fprintf(fid,'%s\n', line);
end
fclose(fid);
fprintf('\nSaved ablation summary: %s\n', outFile);

% ====== 附加：输出 LaTeX 表格格式 ======
texFile = fullfile(outDir, 'Table2_ablation.tex');
fid2 = fopen(texFile,'w');
fprintf(fid2, '\\begin{table}[t]\\centering\\small\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid2, '\\caption{Ablation Study (mean$\\pm$std). Higher is better for AGR/Cov; lower is better for EOR.}\\n');
fprintf(fid2, '\\begin{tabular}{l|ccc|ccc|ccc|ccc}\\hline\\n');
fprintf(fid2, 'Setting & \\multicolumn{3}{c|}{2mm} & \\multicolumn{3}{c|}{4mm} & \\multicolumn{3}{c|}{6mm} & \\multicolumn{3}{c}{8mm}\\\\\\n');
fprintf(fid2, ' & AGR & EOR & Cov & AGR & EOR & Cov & AGR & EOR & Cov & AGR & EOR & Cov\\\\\\hline\\n');
for i = 1:numel(names)
    nm = names{i};
    fprintf(fid2, '%s', nm);
    for t = 1:numel(thicks)
        M = metrics.(nm).(thicks{t});
        fprintf(fid2, ' & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f', ...
            M.AGR_mean, M.AGR_std, M.EOR_mean, M.EOR_std, M.Cov_mean, M.Cov_std);
    end
    fprintf(fid2, '\\\\\\n');
end
fprintf(fid2, '\\hline\\end{tabular}\\end{table}\\n');
fclose(fid2);
fprintf('Saved LaTeX table: %s\n', texFile);
end

% ====== Helper functions ======

function [AGRo, AGRl, EORo, EORl, CovO, CovL, nOK] = eval_on_set_full(dataFolder, idxList, params)
AGRo = nan(numel(idxList),1);
AGRl = nan(numel(idxList),1);
EORo = nan(numel(idxList),1);
EORl = nan(numel(idxList),1);
CovO = nan(numel(idxList),1);
CovL = nan(numel(idxList),1);
parfor k = 1:numel(idxList)
    i = idxList(k);
    [sig,fs] = load_one(dataFolder, i);
    if isempty(sig), continue; end
    [agro, agrl, eoro, eorl, covO, covL] = run_detection_with_thresholds(sig, fs, params);
    AGRo(k)=agro; AGRl(k)=agrl; EORo(k)=eoro; EORl(k)=eorl; CovO(k)=covO; CovL(k)=covL;
end
nOK = sum(~isnan(AGRl));
end

function [sig,fs] = load_one(dataFolder, i)
sig=[]; fs=[];
matFile = fullfile(dataFolder, sprintf('data%d.mat', i));
if exist(matFile,'file')
    S = load(matFile);
    if isfield(S,'signal_denoised'), sig = S.signal_denoised;
    elseif isfield(S,'sig_deno'), sig = S.sig_deno;
    elseif isfield(S,'signal'), sig = S.signal;
    elseif isfield(S,'sig'), sig = S.sig;
    end
if ~exist(matFile,'file')
    fprintf(2,'[WARN] Missing: %s\n', matFile);
end

    if isempty(sig) && isfield(S,'signal'), sig = S.signal; end
    if isfield(S,'fs'), fs=S.fs; else, fs=2e6; end
elseif exist(strrep(matFile,'.mat','.txt'),'file')
    sig = load(strrep(matFile,'.mat','.txt')); fs=2e6;
else
    fprintf('  [skip] data%d not found.\n',i);
end
sig = sig(:);
sig = single(sig);  
end
%% ======================== Detection Core（可调阈值） ========================
function [AGR_orig, AGR_layered, EOR_orig, EOR_layered, Cov_orig, Cov_layered] = run_detection_with_thresholds(signal, fs, params)
if isfield(params,'use_layered') && ~params.use_layered
    % 强制单层检测：跳过 weak 通道与 layered 组合
    use_layered_mode = false;
else
    use_layered_mode = true;
end

% ---- TF 熵 + 熵率 ----
time = (0:numel(signal)-1)/fs;
win_len = 200e-6; step = 100e-6;
win_pts  = floor(win_len*fs); step_pts = floor(step*fs);
start_idx  = 1:step_pts:(numel(signal)-win_pts+1);
center_idx = start_idx + floor(win_pts/2);
time_axis  = time(center_idx);

W = numel(start_idx); q = 1.6; tf_entropy = zeros(1,W);
for k = 1:W
    seg = signal(start_idx(k):start_idx(k)+win_pts-1);
    if params.use_sst
        [sst,~] = wsst(seg, fs, 'bump'); P = abs(sst).^2;
    else
        [S,~,~] = spectrogram(seg, hamming(round(win_pts/2)), [], [], fs);
        P = abs(S).^2;
    end
    P = P ./ max(sum(P,1), eps);
    tf_entropy(k) = (1 - sum(P(:).^q)) / (q - 1);
end
% ---- 归一化：分位数(1-99%)优先，否则退回min-max ----
if params.use_percentile_norm
    lo = prctile(tf_entropy, params.entropy_pct_lo);
    hi = prctile(tf_entropy, params.entropy_pct_hi);
    tf_entropy_norm = (tf_entropy - lo) / max(hi - lo, eps);
else
    tf_entropy_norm = (tf_entropy - min(tf_entropy)) / max(eps, (max(tf_entropy)-min(tf_entropy)));
end
tf_entropy_norm = min(max(tf_entropy_norm,0),1);

% ---- 轻平滑后再求熵率 ----
tf_entropy_smooth = movmedian(tf_entropy_norm, max(1,params.rate_smooth_win));
entropy_rate      = diff(tf_entropy_smooth) ./ diff(time_axis);
time_rate         = time_axis(2:end);
rate_mean         = mean(entropy_rate);
rate_std          = std(entropy_rate);


% ---- 主通道候选（熵 + 熵率） ----
thrE  = params.entropy_thresh;

% 二选一：固定熵率阈 or kσ 阈
if params.use_fixed_rate
    rate_ok = @(a,b) (min(entropy_rate(a:b)) < params.fixed_rate_thr);
else
    thrR  = rate_mean - params.rate_k * rate_std;
    rate_ok = @(a,b) (min(entropy_rate(a:b)) < thrR);
end

idx_low = find(tf_entropy_smooth(2:end) < thrE);  % 注意这里用 smooth 后的熵
ms2rate = @(ms) max(1, round((ms*1e-3) / (step_pts/fs)));
cand = [];
wR = ms2rate(1.0);

for idx = idx_low'
    L=max(1,idx-wR); R=min(numel(entropy_rate), idx+wR);
    if rate_ok(L,R), cand = [cand idx]; end %#ok<AGROW>
end

cand = unique(cand(:));


% cand (rate+entropy) 已得到 -> 转成采样点索引
if params.use_energy_gates && ~isempty(cand)
    window_energy = round(1.0e-3*fs);
    energy_bg_global = movsum(signal.^2, window_energy);
    Eglob = params.energy_global_mult * median(energy_bg_global);

    cand_keep = false(size(cand));
    for kk = 1:numel(cand)
        samp_idx = center_idx(cand(kk));
        L = max(1, samp_idx - window_energy);
        R = min(numel(signal), samp_idx + window_energy);
        win_sig = signal(L:R);
        elocal  = sum(win_sig.^2);
        burst   = kurtosis(abs(win_sig));
        pass = (elocal > Eglob) && (burst > params.burst_thresh);

        if params.use_dynamic_bg
            window_small = max(1, round(0.05e-3*fs));
            window_bg    = max(window_small+1, round(2.5e-3*fs));
            Lbg = max(1, samp_idx - window_bg); Rbg = min(numel(signal), samp_idx + window_bg);
            bg_seg = signal(Lbg:Rbg);
            ebg_series = movsum(bg_seg.^2, window_small);
            thr_dyn = mean(ebg_series) + params.alpha*std(ebg_series);
            pass = pass && (elocal > thr_dyn);
        end

        cand_keep(kk) = pass;
    end
    cand = cand(cand_keep);
    % ---- 追加：MAD-SNR 物理门 ----
if ~isempty(cand)
    cand_keep2 = false(size(cand));
    for kk = 1:numel(cand)
        samp_idx = center_idx(cand(kk));

        % 局部包络峰
        L = max(1, samp_idx - params.win_energy_pts);
        R = min(numel(signal), samp_idx + params.win_energy_pts);
        env = abs(hilbert(signal(L:R)));
        env = movmean(env, max(1, round(params.env_smooth_ms*1e-3*fs)));
        pk  = max(env);

        % 背景 MAD
        Lbg = max(1, samp_idx - params.bg_win_pts);
        Rbg = min(numel(signal), samp_idx + params.bg_win_pts);
        bg_env = abs(hilbert(signal(Lbg:Rbg)));
        bg_med = median(bg_env);
        bg_mad = 1.4826*median(abs(bg_env - bg_med)) + eps;
        snr_pk = (pk - bg_med)/bg_mad;

        cand_keep2(kk) = (snr_pk >= params.snr_thresh);
    end
    cand = cand(cand_keep2);
end

end

% 最小间隔（参数化）
if ~isempty(cand)
    AE_times_main = sort(time_axis(cand)*1000);
    t_entropy = [];
    for i = 1:numel(AE_times_main)
        if isempty(t_entropy) || AE_times_main(i) - t_entropy(end) > params.min_interval_ms
            t_entropy(end+1,1) = AE_times_main(i);
        end
    end
else
    t_entropy = [];
end



% ---- 包络候选（STA/LTA） ----
env_ms = stalta_detect_param(signal, fs, time, 0, time(end), params);  % 传入 params


% ---- 原始匹配 → HC/LC/EO ----
tol_ms = params.tol_ms;
[HC, LC, EO] = match_pseudo_labels(t_entropy, env_ms, tol_ms);

% ---- 微弱通道（放宽阈值） ----
thrE_w = params.entropy_thresh_weak;
thrR_w = rate_mean - params.rate_k_weak * rate_std;
idx_w  = find(tf_entropy_norm(2:end) < thrE_w & tf_entropy_norm(2:end) >= thrE);
cand_w = [];
wR = ms2rate(1.5);
for idx = idx_w'
    L=max(1,idx-wR); R=min(numel(entropy_rate), idx+wR);
    if min(entropy_rate(L:R)) < thrR_w, cand_w = [cand_w idx]; end
end
t_weak = sort(time_axis(unique(cand_w))*1000);
if params.use_weak_scoring && ~isempty(t_weak)
    idx_wu = unique(cand_w);  % 与 t_weak 对应的 rate 索引
    scores = zeros(numel(idx_wu),1);

    % 预先统计
    window_small = max(1, round(0.10e-3*fs));
    window_bg    = max(window_small+1, round(2.0e-3*fs));
    for ii=1:numel(idx_wu)
        rate_idx = idx_wu(ii);
        samp_idx = center_idx(rate_idx);
        L = max(1, samp_idx-window_small); R = min(numel(signal), samp_idx+window_small);
        win_sig = signal(L:R);
        elocal  = sum(win_sig.^2);
        burst   = kurtosis(abs(win_sig));

        Lbg = max(1, samp_idx-window_bg); Rbg = min(numel(signal), samp_idx+window_bg);
        bg_seg = signal(Lbg:Rbg);
        ebg_series = movsum(bg_seg.^2, window_small);
        mean_bg = mean(ebg_series); std_bg = std(ebg_series);

        % 特征（与你 Multi-step 一致/相近）
        rate_mean = mean(entropy_rate); rate_std = std(entropy_rate);
        feat = zeros(1,6);
        feat(1) = max(0, 1 - (tf_entropy_norm(rate_idx)-params.entropy_thresh) / max(1e-6,(params.entropy_thresh_weak-params.entropy_thresh)));
        feat(2) = min(1, abs(entropy_rate(max(1,rate_idx-1)) - rate_mean)/(3*(rate_std+eps)));
        feat(3) = max(0, min(1, log10(elocal/(mean_bg+eps))/1.5));
        feat(4) = max(0, min(1, (burst-1)/2));
        % 简化版 SNR/频域比：可按你原脚本细化
        signal_power = max(movsum(win_sig.^2, max(1, round(numel(win_sig)/4))));
        noise_power  = median(movsum(win_sig.^2, max(1, round(numel(win_sig)/4))));
        snr = signal_power/(noise_power+eps);
        feat(5) = max(0, min(1, log10(snr)/2));
        feat(6) = 0.3; % 若需要，可用 pwelch 计算频域比

        w = [0.25,0.20,0.15,0.15,0.15,0.10];
        scores(ii) = sum(w.*feat);
    end

    % 根据 scores 过滤 + NMS
    keep = scores >= params.weak_quality_cut;
    t_weak = t_weak(keep);
    scores = scores(keep);

    if ~isempty(t_weak)
        [t_weak, ord] = sort(t_weak); scores = scores(ord);
        keep_mask = true(size(t_weak));
        for i=2:numel(t_weak)
            if t_weak(i)-t_weak(i-1) < params.weak_min_sep_ms
                if scores(i) > scores(i-1), keep_mask(i-1)=false; else, keep_mask(i)=false; end
            end
        end
        t_weak = t_weak(keep_mask);
    end
end

% ---- 分层：weak 与剩余 env 最近邻成 MC；剩余 env 为 EO_layered；weak 未配到的记 WC（不计入发布清单） ----
% 先剔除 HC 已用的 env
used_env = false(numel(env_ms),1);
for i=1:numel(HC)
    [~,j] = min(abs(env_ms-HC(i))); if abs(env_ms(j)-HC(i))<=tol_ms, used_env(j)=true; end
end
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

% ---- 指标 ----
% ---- 指标 ----
hi_conf = numel(HC)+numel(MC);
tot_ent = hi_conf + numel(LC);
Final_l = sort([HC(:); MC(:)]);     % Layered 的最终事件
Final_o = sort(HC(:));              % Original 的最终事件（只算 HC）

AGR_orig     = numel(HC)/max(1, (numel(HC)+numel(LC)));
EOR_orig     = numel(EO)/max(1, (numel(HC)+numel(EO)));
AGR_layered  = hi_conf/max(1,tot_ent);
EOR_layered  = numel(EO_layered)/max(1, (hi_conf+numel(EO_layered)));

Cov_orig     = numel(Final_o)/max(1, (numel(Final_o)+numel(EO)));            % ★ 新增
Cov_layered  = numel(Final_l)/max(1, (numel(Final_l)+numel(EO_layered)));    % 原 Coverage 概念

if ~use_layered_mode
    AGR_layered = AGR_orig;
    EOR_layered = EOR_orig;
    Cov_layered = Cov_orig;
    return;
end

end
%% ======================== STA/LTA 包络检测 ========================
function t_ms = stalta_detect_param(sig, fs, time, t0, t1, params)
% STA/LTA 窗、阈值
sta_ms = 0.10;  lta_ms = 2.0;          % 你原先常用的配置
sta = max(1, round(sta_ms*1e-3*fs));
lta = max(sta+1, round(lta_ms*1e-3*fs));

% 包络（可选轻去噪）
if isfield(params,'pre_denoised') && params.pre_denoised
    sig4env = sig;
else
    sig4env = movmean(sig, 3);
end
env = abs(hilbert(sig4env));
env = movmean(env, max(1, round(params.env_smooth_ms*1e-3*fs)));

% STA/LTA 比值
S = movmean(env, sta);
L = movmedian(env, lta);
ratio = S ./ (L + eps);

% 阈值模式
switch lower(params.stalta_mode)
  case 'ratio'
    thr = params.stalta_thresh;
  otherwise % 'quantile'
    valid = ratio(isfinite(ratio));
    if isempty(valid), t_ms = []; return; end
    thr = quantile(valid, params.stalta_quantile);
end

% 超阈区段→峰值
above = ratio > thr;
d = diff([0; above(:); 0]);
st = find(d==1); ed = find(d==-1)-1;

% 最短持续（过滤毛刺）
% min_dur_ms  = 0.04;
% min_dur_pts = max(1, round(min_dur_ms*1e-3*fs));
% if isfield(params,'min_dur_ms'), min_dur_ms = params.min_dur_ms; else, min_dur_ms = 0.06; end
min_dur_pts = max(1, round(params.min_dur_ms*1e-3*fs));

evt = [];
for i = 1:numel(st)
    seg = st(i):ed(i);
    if numel(seg) < min_dur_pts, continue; end
    [~,ix] = max(ratio(seg));
    evt(end+1) = st(i)+ix-1; %#ok<AGROW>
end

% 合并相近事件（和主检测一致的 NMS 间隔更和谐）
% min_sep_ms  = 0.40;
% min_sep_pts = max(1, round(min_sep_ms*1e-3*fs));
% if isfield(params,'min_sep_ms'), min_sep_ms = params.min_sep_ms; else, min_sep_ms = 0.50; end
min_sep_pts = max(1, round(params.min_sep_ms*1e-3*fs));
evt = sort(evt(:));
keep = true(size(evt));
for i = 2:numel(evt)
    if (evt(i)-evt(i-1)) < min_sep_pts
        % 保留 ratio 高的
        if ratio(evt(i)) >= ratio(evt(i-1))
            keep(i-1) = false;
        else
            keep(i) = false;
        end
    end
end
evt = evt(keep);

% 时间裁剪到 [t0,t1]
t_ms = time(evt)*1000;
mask = (t_ms >= t0*1000) & (t_ms <= t1*1000);
t_ms = t_ms(mask);
end
function write_pvalues_csv(SUM, outPath)
keys = fieldnames(SUM);                      % {'mm2','mm4','mm6','mm8'}
nums = cellfun(@(k) str2double(regexprep(k,'^mm','')), keys);
[~,ord] = sort(nums);
keys = keys(ord);

fid = fopen(outPath,'w');
fprintf(fid,'Thickness,p_AGR,p_EOR,p_Cov\n');
for i=1:numel(keys)
    key = keys{i};                           % e.g. 'mm2'
    dispName = sprintf('%smm', regexprep(key,'^mm',''));  % '2mm'
    O = SUM.(key).Original; 
    L = SUM.(key).Layered;
    pAGR = paired_p(O.AGR, L.AGR);
    pEOR = paired_p(O.EOR, L.EOR);
    pCov = paired_p(O.Cov, L.Cov);           % 已映射：O.Cov=CovO, L.Cov=CovL
    fprintf(fid,'%s,%.3e,%.3e,%.3e\n', dispName, pAGR, pEOR, pCov);
end
fclose(fid);
end


    function p = paired_p(x,y)
x = x(:); y = y(:);
n = min(numel(x),numel(y));
x = x(1:n); y = y(1:n);
try
    [~,p] = ttest(x,y);
catch
    p = signrank(x,y);
end
    end
    %% ======================== 匹配：HC/LC/EO ========================
function [HC, LC, EO] = match_pseudo_labels(t_entropy, t_env, tol_ms)
t_entropy = t_entropy(:); t_env = t_env(:);
used = false(numel(t_env),1); HC=[]; LC=[];
for i=1:numel(t_entropy)
    if isempty(t_env), LC(end+1,1)=t_entropy(i); continue; end %#ok<AGROW>
    [dmin,j] = min(abs(t_env - t_entropy(i)));
    if dmin<=tol_ms && ~used(j)
        HC(end+1,1)=t_entropy(i); used(j)=true;
    else
        LC(end+1,1)=t_entropy(i);
    end
end
EO = t_env(~used);
end

% function p = remove_entropy_rate_gate(p)
% p.use_fixed_rate = false;   % 禁用固定阈值
% p.rate_k = 0;               % 不再按标准差定义阈值
% p.fixed_rate_thr = -inf;    % 或直接设一个极小值
% end
% 
% function p = remove_weak_scoring(p)
% p.use_weak_scoring = false;
% p.weak_quality_cut = 1.1;   % 所有 weak 都过滤
% p.weak_min_sep_ms  = 0;
% end
% 
% function p = disable_fine_sweep(p)
% p.use_weak_scoring = true;
% p.weak_quality_cut = -inf;  % 不再筛选
% p.weak_min_sep_ms  = 0;
% p.rate_k_weak      = 0;     % 弱层不做动态阈值
% end
% 
% function p = restore_original_mode(p)
% p.use_weak_scoring = false;
% p.use_layered = false;
% p.tol_ms = 0;               % 恢复严格匹配
% end
% 
% function p = disable_entropy_branch(p)
% p.entropy_thresh = 999;     % 禁用熵
% p.use_fixed_rate = false;
% p.rate_k = 0;
% end

%% ======= Helper: Module Ablation =======

function p = remove_entropy_rate_gate(p)
p.use_fixed_rate = false;      % 关闭动态熵率门，转回固定阈值
p.fixed_rate_thr = -500;      % 近似无门控，仅保留熵阈
p.rate_k = 0.05;                 % 无效但保持字段
end
function p = remove_dynamic_energy_gate(p)
% 禁用动态能量门控（只保留全局静态门）
if isfield(p, 'use_dynamic_bg'), p.use_dynamic_bg = false; end
if isfield(p, 'use_energy_gates'), p.use_energy_gates = true; end
if isfield(p, 'energy_global_mult'), p.energy_global_mult = 0.7; end
if isfield(p, 'alpha'), p.alpha = 0; end
end

function p = remove_weak_scoring_soft(p)
% 宽容禁用弱事件评分：不彻底清除，但放松筛选
p.use_weak_scoring = true;
p.weak_quality_cut = 1.1;
p.rate_k_weak = 0.4;  
end

function p = reduce_tolerance(p)
% 缩小匹配容差（更严格匹配）
p.tol_ms = 0.1;              % 从 1.8 → 0.6 ms，减少高置信匹配数
end

function p = restore_original_mode(p)
p.use_layered = false;
p.use_weak_scoring = false;
p.use_dynamic_bg = false;
p.use_fixed_rate = true;
p.fixed_rate_thr = -1000;       % 放松熵率门，弱化检测强度
p.tol_ms = 0.6;                 % 更严格匹配
p.alpha = 0;
end


function p = disable_entropy_branch(p)
% 禁用熵通道：仅保留包络
% p.entropy_thresh = 1.0;      % 极高阈值，但不致使全失效
% p.use_fixed_rate = false;
% p.rate_k = 0;
p.entropy_thresh = 999;
p.use_fixed_rate = false;
p.rate_k = 0;
end
function p = disable_envelope_branch_soft(p)
p.stalta_mode = 'ratio';
p.stalta_thresh = 8.0;            % 实质关闭包络
p.use_fixed_rate = false;
p.rate_k = 0.8;                   % 允许一定波动，避免过严
p.entropy_thresh = 0.35;          % 放宽熵阈
p.use_dynamic_bg = false;
p.use_energy_gates = false;
p.energy_global_mult = 0.0;
p.alpha = 0;
p.min_dur_ms = 0.05;              % 允许短事件
p.min_interval_ms = 0.2;
p.tol_ms = 1;
end





