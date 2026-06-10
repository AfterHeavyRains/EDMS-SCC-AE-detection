%%全部改了并行的版本
% 共享主阈值重新搜索origin和layered最优 original优先
%  function optimize_and_sweep_2mm()
clear
clearvars -except DO_* rootPath dataFolder saveDir nSamples subsetN subsetIdx best_params
clc; if isempty(gcp('nocreate')), parpool; end
%% =============== Paths & basic config ===============
rootPath   = 'E:\Code\SCC_matlab_code\results\results\result2mm\result2mm';
dataFolder = fullfile(rootPath, 'data_denoised');   % data1.mat ... data200.mat
saveDir    = rootPath;
save_fig   = @(h,fn) exportgraphics(h, fullfile(saveDir,fn), 'Resolution', 300);
nSamples   = 200;

useSubset  = true;        % 优化阶段只用子集，提速
subsetN    = 30;
subsetIdx  = 1:min(subsetN,nSamples);

% 目标函数权重（两步都用）
w1=1.0;  % AGR 正向
w2=1.0;  % EOR 反向
w3=0.6;  % Coverage 正向
coverage_min = 0.60;      % 低于此覆盖率给惩罚（别太高，否则 Original 很难）
EOR_cap_O = 0.80;         % Original 允许的 EOR 上限（软惩罚）
EOR_cap_L = 0.60;         % Layered 允许的 EOR 上限（软惩罚）

% Layered 的固定开关（保留你现有选择）
fixed.use_weak_scoring   = true;
fixed.weak_quality_cut   = 0.4;
fixed.weak_min_sep_ms    = 3.0;
fixed.use_energy_gates   = true;
fixed.energy_global_mult = 1.1;
fixed.burst_thresh       = 1.2;
fixed.use_dynamic_bg     = true;
fixed.alpha              = 3.0;
% === defaults for gating & dynamic BG ===
if ~exist('best_params','var'), best_params = struct; end
if ~isfield(best_params,'use_energy_gates'),   best_params.use_energy_gates = true; end
if ~isfield(best_params,'energy_global_mult'), best_params.energy_global_mult = 1.1; end
if ~isfield(best_params,'burst_thresh'),       best_params.burst_thresh      = 1.2; end
if ~isfield(best_params,'use_dynamic_bg'),     best_params.use_dynamic_bg    = true; end
if ~isfield(best_params,'alpha'),              best_params.alpha             = 3.0; end

%% ===== Step-O (Original-first, maximize AGRo) =====
fprintf('\n[Step-O] Optimize Original baseline (maximize AGRo) on %d-sample subset...\n', numel(subsetIdx));

% —— 基线参数（保留你的门控；如果只追 AGRo，可适度放松 SNR 门）——
base = struct( ...
  'use_sst', false, ...
  'use_percentile_norm', true, 'entropy_pct_lo', 1, 'entropy_pct_hi', 99, 'rate_smooth_win', 3, ...
  'use_fixed_rate', true, 'fixed_rate_thr', -1900, ...
  'tol_ms', 1.0, 'min_interval_ms', 2.6, ...
  'use_energy_gates', true, 'energy_global_mult', 1.15, 'burst_thresh', 1.30, ...
  'use_dynamic_bg', true, 'alpha', 3.0, ...
  'stalta_mode','ratio','stalta_thresh',2.3,'env_smooth_ms',0.02, ...
  'snr_thresh', 9.0, ...              % 原 10 → 9；更容易保留“可配”的峰，方便 AGRo 上来
  'win_energy_pts', 200, 'bg_win_pts', 12000, ...
  'use_weak_scoring', false, ...
  'pre_denoised', true);

% —— 网格：主阈值更宽、固定熵率更负、更宽 tol —— 
ent_grid = 0.50:0.02:0.80;                    % ↑ 主阈值范围拉宽
frt_grid = [-1700 -1800 -1900 -2000 -2100 -2200 -2300];  % ↑ 更负更严
tol_grid = [0.8 1.0 1.2 1.4 1.6 1.8];         % ↑ 宽容差，尽量把 HC 匹上

% （可选）如果你想一起扫最小合并间隔，也可以打开：
% minsep_grid = [2.0 2.6 3.0];  % 变大会合并近邻，常能提高 AGRo；默认不用，稳住即可。

fprintf('[Step-O] grid sizes: ent=%d, frt=%d, tol=%d\n', numel(ent_grid), numel(frt_grid), numel(tol_grid));
ent_i = 0;

rowsO = [];
bestO = struct('score', -Inf);

for ent = ent_grid
  ent_i = ent_i + 1;
  fprintf('  [O] ent %.3f (%d/%d)...\n', ent, ent_i, numel(ent_grid));
  for frt = frt_grid
    for tol = tol_grid
      pO = base;
      pO.entropy_thresh = ent;
      pO.fixed_rate_thr = frt;
      pO.tol_ms         = tol;
      % 如果要扫 minsep_grid，则在此处 pO.min_interval_ms = x;

      [mAGRo, ~, mEORo, ~, mCov] = eval_on_set(dataFolder, subsetIdx, pO);
      if any(isnan([mAGRo mEORo mCov])), continue; end

      % ✅ 只看 AGRo；给极小权重的辅助项，避免完全异常解
      score = mAGRo + 0.03*mCov - 0.01*mEORo;

      rowsO(end+1,:) = [ent frt tol mAGRo mEORo mCov score]; %#ok<AGROW>

      if score > bestO.score
        bestO = struct('ent',ent,'frt',frt,'tol',tol, ...
          'AGRo',mAGRo,'EORo',mEORo,'Cov',mCov,'score',score);
      end
    end
  end
end

if isempty(rowsO)
  error('[Step-O] No feasible Original found. Try widening ent_grid (e.g., 0.46:0.02:0.86) or tol_grid up to 2.0, or set snr_thresh=8.5.');
end

To = array2table(rowsO,'VariableNames',{'ent','fixed_rate_thr','tol','AGRo','EORo','Cov','Score'});

% 打印 Top-10（按 AGRo 优先，再按 score）
[~,ord] = sortrows([To.AGRo, To.Score],[ -1, -1 ]); % AGRo 降序，其次 score
topN = min(10,height(To));
fprintf('[Step-O] Top-%d Original candidates (by AGRo):\n', topN);
for ii=1:topN
  r=To(ord(ii),:);
  fprintf('  #%d: ent=%.3f frt=%.0f tol=%.2f | AGRo=%.3f EORo=%.3f Cov=%.3f Score=%.3f\n', ...
    ii, r.ent, r.fixed_rate_thr, r.tol, r.AGRo, r.EORo, r.Cov, r.Score);
end
fprintf('[Step-O] Baseline*: ent=%.3f, frt=%.0f, tol=%.2f | AGRo=%.3f, EORo=%.3f, Cov=%.3f, Score=%.3f\n', ...
  bestO.ent, bestO.frt, bestO.tol, bestO.AGRo, bestO.EORo, bestO.Cov, bestO.score);

% 固定 Original baseline 的参数（供 Step-L 复用）
params_baseline = base;
params_baseline.entropy_thresh = bestO.ent;
params_baseline.fixed_rate_thr = bestO.frt;
params_baseline.tol_ms         = bestO.tol;



%% ===== Step-L (Layered on top of baseline ent/rate/tol) =====
fprintf('\n[Step-L] Tune Layered on fixed baseline (ent/tol/fixed-rate from Step-O)...\n');

% 在 baseline 上打开 weak，并微扫 delta / rk_w / quality_cut
weak_delta_grid = 0.16:0.02:0.24;     % ent_w = ent + delta
rk_w_grid       = [0.4 0.6 0.8];      % “弱熵率门”放宽
cut_grid        = [0.45 0.50 0.55];   % 先干净后扩量

rowsL = [];
bestL = struct('score', -Inf);
fprintf('[Step-L] grid sizes: delta=%d, rk_w=%d, cut=%d\n', numel(weak_delta_grid), numel(rk_w_grid), numel(cut_grid));
delta_i = 0;

for delta = weak_delta_grid
  delta_i = delta_i + 1;
  fprintf('  [L] delta %.3f (%d/%d)...\n', delta, delta_i, numel(weak_delta_grid));

  for rkw = rk_w_grid
    for cut = cut_grid
      pL = params_baseline;
      pL.use_weak_scoring     = true;
      pL.entropy_thresh_weak  = pL.entropy_thresh + delta;
      pL.rate_k_weak          = rkw;
      pL.weak_quality_cut     = cut;
      pL.weak_min_sep_ms      = 3.0;   % 保守NMS
      % 其余门控同 baseline

      [~, mAGRl, ~, mEORl, mCov] = eval_on_set(dataFolder, subsetIdx, pL);
      if any(isnan([mAGRl mEORl mCov])), continue; end

      % 简单的 Layered 评分：更看 AGRl，轻惩罚 EORl，覆盖加点分
      score = 1.4*mAGRl - 0.6*mEORl + 0.4*mCov;

      rowsL(end+1,:) = [delta rkw cut mAGRl mEORl mCov score]; %#ok<AGROW>

      if score > bestL.score
        bestL = struct('delta',delta,'rk_w',rkw,'cut',cut, ...
          'AGRl',mAGRl,'EORl',mEORl,'Cov',mCov,'score',score);
      end
    end
  end
end

if isempty(rowsL)
  error('[Step-L] No feasible Layered found. Try lowering weak_quality_cut or increasing delta to 0.26.');
end

Tl = array2table(rowsL,'VariableNames',{'delta','rk_w','cut','AGRl','EORl','Cov','Score'});
[~,ordL] = sort(Tl.Score,'descend');
topN = min(10,height(Tl));
fprintf('[Step-L] Top-%d Layered candidates on baseline:\n', topN);
for ii=1:topN
  r=Tl(ordL(ii),:);
  fprintf('  #%d: delta=%.3f rk_w=%.2f cut=%.2f | AGRl=%.3f EORl=%.3f Cov=%.3f Score=%.3f\n', ...
    ii, r.delta, r.rk_w, r.cut, r.AGRl, r.EORl, r.Cov, r.Score);
end
fprintf('[Step-L] Layered*: delta=%.3f (ent_w=%.3f), rk_w=%.2f, cut=%.2f | AGRl=%.3f, EORl=%.3f, Cov=%.3f, Score=%.3f\n', ...
  bestL.delta, (bestO.ent+bestL.delta), bestL.rk_w, bestL.cut, bestL.AGRl, bestL.EORl, bestL.Cov, bestL.score);

% —— 输出两个最终参数包（后续可用在全量200样本评估）——
params_O_final = params_baseline;  % 原样
params_L_final = params_baseline;  % 在原样基础上加弱通道
params_L_final.use_weak_scoring    = true;
params_L_final.entropy_thresh_weak = params_L_final.entropy_thresh + bestL.delta;
params_L_final.rate_k_weak         = bestL.rk_w;
params_L_final.weak_quality_cut    = bestL.cut;
params_L_final.weak_min_sep_ms     = 3.0;

%% ======================== Helper: 批量评估（子集） ========================
function [mAGRo, mAGRl, mEORo, mEORl, mCov, nOK] = eval_on_set(dataFolder, idxList, params)
AGRo=[]; AGRl=[]; EORo=[]; EORl=[]; Cov=[]; nOK=0;

AGRo_tmp = nan(numel(idxList),1);
AGRl_tmp = nan(numel(idxList),1);
EORo_tmp = nan(numel(idxList),1);
EORl_tmp = nan(numel(idxList),1);
Cov_tmp  = nan(numel(idxList),1);

parfor k = 1:numel(idxList)    % ← 关键
    i = idxList(k);
    [sig,fs] = load_one(dataFolder, i);
    if isempty(sig), continue; end
    [agro, agrl, eoro, eorl, cov] = run_detection_with_thresholds(sig, fs, params);
    AGRo_tmp(k)=agro; AGRl_tmp(k)=agrl; EORo_tmp(k)=eoro; EORl_tmp(k)=eorl; Cov_tmp(k)=cov;
end

AGRo = AGRo_tmp(~isnan(AGRo_tmp));
AGRl = AGRl_tmp(~isnan(AGRl_tmp));
EORo = EORo_tmp(~isnan(EORo_tmp));
EORl = EORl_tmp(~isnan(EORl_tmp));
Cov  = Cov_tmp(~isnan(Cov_tmp));
nOK  = numel(AGRo);

mAGRo = mean(AGRo,'omitnan'); mAGRl = mean(AGRl,'omitnan');
mEORo = mean(EORo,'omitnan'); mEORl = mean(EORl,'omitnan');
mCov  = mean(Cov,'omitnan');
end

%% ======================== Helper: 读取一个样本 ========================
function [sig,fs] = load_one(dataFolder, i)
sig=[]; fs=[];
matFile = fullfile(dataFolder, sprintf('data%d.mat', i));
if exist(matFile,'file')
    S = load(matFile);
    if isfield(S,'signal_denoised'), sig = S.signal_denoised; end
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
function [AGR_orig, AGR_layered, EOR_orig, EOR_layered, Coverage] = run_detection_with_thresholds(signal, fs, params)
if ~isfield(params,'use_sst'), params.use_sst=false; end
% ---- 新增：可选门控与弱通道打分的参数（默认关闭，兼容旧寻优） ----
if ~isfield(params,'use_energy_gates'), params.use_energy_gates = true; end
if ~isfield(params,'energy_global_mult'), params.energy_global_mult = 1.1; end   % 全局能量倍率
if ~isfield(params,'burst_thresh'),      params.burst_thresh      = 1.2; end    % 突发度门
if ~isfield(params,'use_dynamic_bg'),    params.use_dynamic_bg    = true; end  % 动态背景门
if ~isfield(params,'alpha'),             params.alpha             = 3.0; end    % 动态背景系数

if ~isfield(params,'use_weak_scoring'),  params.use_weak_scoring  = true; end  % 弱通道质量打分
if ~isfield(params,'weak_quality_cut'),  params.weak_quality_cut  = 0.4; end
if ~isfield(params,'weak_min_sep_ms'),   params.weak_min_sep_ms   = 3.0; end

% ---- 参数 ----
if ~isfield(params,'entropy_thresh'),      params.entropy_thresh = 0.38; end
if ~isfield(params,'rate_k'),              params.rate_k = 2.0; end
if ~isfield(params,'tol_ms'),              params.tol_ms = 0.3; end
if ~isfield(params,'entropy_thresh_weak'), params.entropy_thresh_weak = params.entropy_thresh + 0.12; end
if ~isfield(params,'rate_k_weak'),         params.rate_k_weak = 1.0; end
% ===== 新增：稳健归一化/平滑/固定熵率阈值/STA-LTA/SNR/合并间隔 =====
if ~isfield(params,'use_percentile_norm'), params.use_percentile_norm = true; end
if ~isfield(params,'entropy_pct_lo'),      params.entropy_pct_lo    = 1;    end
if ~isfield(params,'entropy_pct_hi'),      params.entropy_pct_hi    = 99;   end
if ~isfield(params,'rate_smooth_win'),     params.rate_smooth_win   = 3;    end

if ~isfield(params,'use_fixed_rate'),      params.use_fixed_rate    = false; end
if ~isfield(params,'fixed_rate_thr'),      params.fixed_rate_thr    = -1900; end  % 6mm推荐

% STA/LTA 模式：'ratio' 或 'quantile'
if ~isfield(params,'stalta_mode'),         params.stalta_mode       = 'quantile'; end
if ~isfield(params,'stalta_thresh'),       params.stalta_thresh     = 2.3;   end  % ratio 模式用
if ~isfield(params,'stalta_quantile'),     params.stalta_quantile   = 0.9979;end  % quantile 模式用
if ~isfield(params,'env_smooth_ms'),       params.env_smooth_ms     = 0.02;  end  % ms

% SNR 物理门（MAD）
if ~isfield(params,'snr_thresh'),          params.snr_thresh        = 10;    end
if ~isfield(params,'win_energy_pts'),      params.win_energy_pts    = 200;   end
if ~isfield(params,'bg_win_pts'),          params.bg_win_pts        = 12000; end

% 事件合并的最小时间间隔（ms）
if ~isfield(params,'min_interval_ms'),     params.min_interval_ms   = 2.6;   end

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
hi_conf = numel(HC)+numel(MC);
tot_ent = hi_conf + numel(LC);
Final   = sort([HC(:); MC(:)]);
AGR_orig     = numel(HC)/max(1, (numel(HC)+numel(LC)));
EOR_orig     = numel(EO)/max(1, (numel(HC)+numel(EO)));
AGR_layered  = hi_conf/max(1,tot_ent);
EOR_layered  = numel(EO_layered)/max(1, (hi_conf+numel(EO_layered)));
Coverage     = numel(Final)/max(1, (numel(Final)+numel(EO_layered)));
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

%% ======================== STA/LTA 包络检测 ========================
function t_ms = stalta_detect_param(sig, fs, time, t0, t1, params)
% 读取参数（带默认）
if ~isfield(params,'stalta_mode'),       params.stalta_mode = 'quantile'; end
if ~isfield(params,'stalta_thresh'),     params.stalta_thresh = 2.3; end
if ~isfield(params,'stalta_quantile'),   params.stalta_quantile = 0.9979; end
if ~isfield(params,'env_smooth_ms'),     params.env_smooth_ms = 0.02; end

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
min_dur_ms  = 0.04;
min_dur_pts = max(1, round(min_dur_ms*1e-3*fs));

evt = [];
for i = 1:numel(st)
    seg = st(i):ed(i);
    if numel(seg) < min_dur_pts, continue; end
    [~,ix] = max(ratio(seg));
    evt(end+1) = st(i)+ix-1; %#ok<AGROW>
end

% 合并相近事件（和主检测一致的 NMS 间隔更和谐）
min_sep_ms  = 0.30;
min_sep_pts = max(1, round(min_sep_ms*1e-3*fs));
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


