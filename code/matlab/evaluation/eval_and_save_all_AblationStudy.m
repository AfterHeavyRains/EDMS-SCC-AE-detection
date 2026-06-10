% function eval_and_save_all()
% 文件开头把签名改成：
function eval_and_save_all(mode)
% if nargin>=1 && strcmpi(mode,'sweep')
%     [TOP, ALL] = sweep_small();
%     assignin('base','TOP',TOP);  % 丢回工作区方便你查看
%     assignin('base','ALL',ALL);
%     return;
% end

% 一次评估四厚度（Original+Layered），保存 CSV/MAT/LaTeX，全程不做寻优
% 放到文件末尾的新函数：sweep_small() 用于小范围扫掠；在命令行调用即可。eval_and_save_all('sweep');
% 消融实验（Ablation Study），新增function run_ablation_study() ，在命令行调用run_ablation_study即可。


%% ===== 基本路径与目录拼接 =====
rootBase = 'E:\Code\SCC_matlab_code\results\Entropy_Detection_Viz_updated_all';
mkData   = @(mm) fullfile(rootBase, sprintf('result%imm',mm), sprintf('result%imm',mm), 'data_denoised');
dataDirs = struct('mm2', mkData(2), 'mm4', mkData(4), 'mm6', mkData(6), 'mm8', mkData(8));
thicks   = {'mm2','mm4','mm6','mm8'};
tlabel   = struct('mm2','2mm','mm4','4mm','mm6','6mm','mm8','8mm');

outDir = fullfile(rootBase, 'eval_exports');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ===== 数据规模 =====
% 若有缺失样本，直接把 idxList 改成存在的 id 列表
nSamples = 200;
idxList  = 1:nSamples;

%% ===== 你的“最佳参数”（来自日志）=====
% Original（基线）
bestO.ent = struct('mm2',0.520,'mm4',0.580,'mm6',0.520,'mm8',0.520);
bestO.tol = 1.80;
bestO.frt = -1700;                 % 固定熵率阈（对应 frt=-1700）

% Layered（分层）
bestL.delta = 0.240;               % entropy weak 增量
bestL.cut   = 0.45;                % 弱事件质量阈
% 说明：你当前 run_detection_with_thresholds 里没用 ent_w、rk_w，这里不强行注入

%% ===== 评估与落盘开关 =====
SAVE_LONG_CSV       = true;        % 统一长表
SAVE_PER_THICK_CSV  = true;        % 每厚度单独 CSV
SAVE_PER_THICK_MAT  = true;        % 每厚度 MAT（含原始向量与params）
WRITE_TABLE1_LATEX  = true;        % 自动生成 Table 1 的 LaTeX
WRITE_TABLE1_WIDECSV= true;        % 生成 mean/std 汇总 CSV

%% ===== 预备输出 =====
longCSV = fullfile(outDir, 'AGR_EOR_Cov_per_sample.csv');
if SAVE_LONG_CSV
    fid = fopen(longCSV,'w');
    fprintf(fid,'sample_id,thickness,method,AGR,EOR,Cov\n');
end

% 并行（可选）
if isempty(gcp('nocreate')), parpool; end

% 汇总容器（为了写 wide-CSV / LaTeX）
SUM = struct();  % SUM.(thickness).(method).{AGR,EOR,Cov} = vector

%% ===== 主循环：四厚度 =====
for t = 1:numel(thicks)
    tk   = thicks{t};
    Tnm  = tlabel.(tk);
    ddir = dataDirs.(tk);

    % ---- 构造参数（统一入口）----
    params = struct();
    params.use_fixed_rate    = true;
    params.fixed_rate_thr    = bestO.frt;
    params.entropy_thresh    = bestO.ent.(tk);
    params.tol_ms            = bestO.tol;

    % Layered 相关
    params.entropy_thresh_weak = params.entropy_thresh + bestL.delta;
    params.use_weak_scoring    = true;
    params.weak_quality_cut    = bestL.cut;


    params.rate_k_weak         = 0.40;   % ★ 必加，复现 rk_w
    params.weak_quality_cut    = bestL.cut;  % 0.45
    params.weak_min_sep_ms     = 3.0;    %（可保留默认，也可显式写上）

    % 其余门控/平滑/STA-LTA/SNR等都走你函数里的默认值

% ---- 逐样本评估 ----
AGRo = nan(numel(idxList),1); AGRl = AGRo;
EORo = AGRo; EORl = AGRo; 
CovO = AGRo; CovL = AGRo;             % ★ 新增两列
okIds = false(numel(idxList),1);

parfor k = 1:numel(idxList)
    i = idxList(k);
    [sig,fs] = load_one(ddir, i);
    if isempty(sig), continue; end
    try
        % ★ 现在接 6 个输出（第5个是 Cov_orig，第6个是 Cov_layered）
        [agro, agrl, eoro, eorl, covO, covL] = run_detection_with_thresholds(sig, fs, params);
        AGRo(k)=agro; AGRl(k)=agrl; EORo(k)=eoro; EORl(k)=eorl; 
        CovO(k)=covO; CovL(k)=covL;                              % ★ 赋值
        okIds(k)=true;
    catch ME
        fprintf(2,'[WARN] %s #%d failed: %s\n', Tnm, i, ME.message);
    end
end


    % ---- 写入长表 CSV ----
    if SAVE_LONG_CSV
        for k = find(okIds).'
            i = idxList(k);
            fprintf(fid,'%d,%s,%s,%.6f,%.6f,%.6f\n', i, Tnm, 'Original', AGRo(k), EORo(k), CovO(k));
            fprintf(fid,'%d,%s,%s,%.6f,%.6f,%.6f\n', i, Tnm, 'Layered',  AGRl(k), EORl(k), CovL(k));

        end
    end

    % ---- 每厚度分表 CSV/MAT ----
    if SAVE_PER_THICK_CSV
        TT = table( idxList(:), repmat({Tnm},numel(idxList),1), ...
                    AGRo, AGRl, EORo, EORl, CovO, CovL, okIds, ...
         'VariableNames',{'sample_id','thickness','AGRo','AGRl','EORo','EORl','CovO','CovL','valid'});

        writetable(TT, fullfile(outDir, sprintf('per_sample_%s.csv', Tnm)));
    end
    if SAVE_PER_THICK_MAT
        params_used = params; %#ok<NASGU>
        save(fullfile(outDir, sprintf('per_sample_%s.mat', Tnm)), ...
             'idxList','AGRo','AGRl','EORo','EORl','CovO','CovL','okIds','params_used','-v7.3');


    end

% ---- 汇总入内存（用于 wide-CSV 和 LaTeX）----
SUM.(tk).Original.AGR = AGRo(okIds);
SUM.(tk).Original.EOR = EORo(okIds);
SUM.(tk).Original.Cov = CovO(okIds);

SUM.(tk).Layered.AGR  = AGRl(okIds);
SUM.(tk).Layered.EOR  = EORl(okIds);
SUM.(tk).Layered.Cov  = CovL(okIds);

fprintf('[%s] n=%d | AGRo=%.3f AGRl=%.3f | EORo=%.3f EORl=%.3f | CovO=%.3f CovL=%.3f\n',...
    Tnm, sum(okIds), ...
    mean(AGRo(okIds),'omitnan'), mean(AGRl(okIds),'omitnan'), ...
    mean(EORo(okIds),'omitnan'), mean(EORl(okIds),'omitnan'), ...
    mean(CovO(okIds),'omitnan'), mean(CovL(okIds),'omitnan'));

end

if SAVE_LONG_CSV, fclose(fid); fprintf('Saved LONG CSV: %s\n', longCSV); end

%% ===== 生成 wide-CSV（mean±std）与 Table 1 LaTeX =====
if WRITE_TABLE1_WIDECSV || WRITE_TABLE1_LATEX
    [WIDE, latexStr] = build_table1_from_SUM(SUM);
    if WRITE_TABLE1_WIDECSV
        writetable(WIDE, fullfile(outDir, 'Table1_wide_summary.csv'));
    end
    if WRITE_TABLE1_LATEX
        texPath = fullfile(outDir, 'Table1_main.tex');
        fid2 = fopen(texPath,'w'); fprintf(fid2, '%s\n', latexStr); fclose(fid2);
        fprintf('Saved Table 1 LaTeX: %s\n', texPath);
    end
end

% 保存一次配置快照（便于复现实验）
cfg = struct();
cfg.time     = char(datetime('now'));
cfg.rootBase = rootBase;
cfg.dataDirs = dataDirs;
cfg.bestO    = bestO;
cfg.bestL    = bestL;
cfg.nSamples = nSamples;
cfg.idxList  = idxList;
save(fullfile(outDir,'config_snapshot.mat'),'-struct','cfg');

fprintf('All done. Outputs in: %s\n', outDir);
end

%% ====== 工具：从 SUM 生成 wide-CSV 和 LaTeX ======
function [WIDE, latexStr] = build_table1_from_SUM(SUM)
thicks = fieldnames(SUM);  % {'2mm','4mm','6mm','8mm'} 或类似
rows = [];
for i=1:numel(thicks)
    tk = thicks{i};
    O  = SUM.(tk).Original; 
    L  = SUM.(tk).Layered;

    mAGRo = mean(O.AGR,'omitnan'); sAGRo = std(O.AGR,'omitnan');
    mAGRl = mean(L.AGR,'omitnan'); sAGRl = std(L.AGR,'omitnan');

    mEORo = mean(O.EOR,'omitnan'); sEORo = std(O.EOR,'omitnan');
    mEORl = mean(L.EOR,'omitnan'); sEORl = std(L.EOR,'omitnan');

    % 注意：如果你希望 Original/Layered 各有自己的 Coverage，
    % SUM 里要分别存 CovO/CovL。若目前只有一份 Cov，则下行改成对应字段。
    mCovo = mean(O.Cov,'omitnan'); sCovo = std(O.Cov,'omitnan');
    mCovl = mean(L.Cov,'omitnan'); sCovl = std(L.Cov,'omitnan');

    rows = [rows; {tk, ...
        mAGRo, sAGRo, mAGRl, sAGRl, ...
        mEORo, sEORo, mEORl, sEORl, ...
        mCovo, sCovo, mCovl, sCovl}]; %#ok<AGROW>
end

WIDE = cell2table(rows, 'VariableNames', { ...
    'Thickness', ...
    'AGRo_mean','AGRo_std','AGRl_mean','AGRl_std', ...
    'EORo_mean','EORo_std','EORl_mean','EORl_std', ...
    'Covo_mean','Covo_std','Covl_mean','Covl_std'});

% ========== 直接在本函数内生成 LaTeX ==========
% 统一显示成 “2mm/4mm/6mm/8mm”
dispTk = @(s) sprintf('%smm', regexp(s,'\d+','match','once'));

buf = strings(0,1);
buf(end+1) = "\begin{table}[t]";
buf(end+1) = "\centering";
buf(end+1) = "\caption{Overall performance across thicknesses (mean$\pm$std, per-sample). Higher is better for AGR and Coverage; lower is better for EOR.}";
buf(end+1) = "\renewcommand{\arraystretch}{1.2}";
buf(end+1) = "\begin{tabular}{c|cc|cc|cc}";
buf(end+1) = "\hline";
buf(end+1) = "\multirow{2}{*}{Thickness} & \multicolumn{2}{c|}{AGR$\uparrow$} & \multicolumn{2}{c|}{EOR$\downarrow$} & \multicolumn{2}{c}{Coverage$\uparrow$} \\";
buf(end+1) = " & Orig. & Layered & Orig. & Layered & Orig. & Layered \\";
buf(end+1) = "\hline";

for i=1:height(WIDE)
    tk   = dispTk(WIDE.Thickness{i});
    line = sprintf('%s & %.3f$\\pm$%.3f & \\textbf{%.3f}$\\pm$%.3f & %.3f$\\pm$%.3f & \\textbf{%.3f}$\\pm$%.3f & %.3f$\\pm$%.3f & \\textbf{%.3f}$\\pm$%.3f \\\\', ...
        tk, ...
        WIDE.AGRo_mean(i), WIDE.AGRo_std(i), WIDE.AGRl_mean(i), WIDE.AGRl_std(i), ...
        WIDE.EORo_mean(i), WIDE.EORo_std(i), WIDE.EORl_mean(i), WIDE.EORl_std(i), ...
        WIDE.Covo_mean(i), WIDE.Covo_std(i), WIDE.Covl_mean(i), WIDE.Covl_std(i));
    buf(end+1) = string(line);
end

buf(end+1) = "\hline";
buf(end+1) = "\end{tabular}";
buf(end+1) = "\label{tab:main}";
buf(end+1) = "\end{table}";

latexStr = strjoin(buf, newline);
end


%% ======================== Helper: 批量评估（子集） ========================
function [mAGRo, mAGRl, mEORo, mEORl, mCovO, mCovL, nOK] = eval_on_set(dataFolder, idxList, params)
% 统一按“逐样本 → *_tmp 收集 → 聚合”的模式
AGRo_tmp = nan(numel(idxList),1);
AGRl_tmp = nan(numel(idxList),1);
EORo_tmp = nan(numel(idxList),1);
EORl_tmp = nan(numel(idxList),1);
CovO_tmp = nan(numel(idxList),1);   % ★
CovL_tmp = nan(numel(idxList),1);   % ★

parfor k = 1:numel(idxList)
    i = idxList(k);
    [sig,fs] = load_one(dataFolder, i);
    if isempty(sig), continue; end
    % ★ 接收 6 个输出（见下一个修正版函数）
    [agro, agrl, eoro, eorl, covO, covL] = run_detection_with_thresholds(sig, fs, params);
    AGRo_tmp(k)=agro; AGRl_tmp(k)=agrl; EORo_tmp(k)=eoro; EORl_tmp(k)=eorl; 
    CovO_tmp(k)=covO; CovL_tmp(k)=covL;                                   % ★
end

% 有效样本筛除 NaN
mask = ~isnan(AGRo_tmp);
AGRo = AGRo_tmp(mask);  AGRl = AGRl_tmp(mask);
EORo = EORo_tmp(mask);  EORl = EORl_tmp(mask);
CovO = CovO_tmp(mask);  CovL = CovL_tmp(mask);
nOK  = sum(mask);

% 聚合
mAGRo = mean(AGRo,'omitnan'); 
mAGRl = mean(AGRl,'omitnan');
mEORo = mean(EORo,'omitnan'); 
mEORl = mean(EORl,'omitnan');
mCovO = mean(CovO,'omitnan');   % ★
mCovL = mean(CovL,'omitnan');   % ★
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
function [AGR_orig, AGR_layered, EOR_orig, EOR_layered, Cov_orig, Cov_layered] = run_detection_with_thresholds(signal, fs, params)
if ~isfield(params,'use_sst'), params.use_sst=false; end
% ---- 新增：可选门控与弱通道打分的参数（默认关闭，兼容旧寻优） ----
if ~isfield(params,'use_energy_gates'), params.use_energy_gates = true; end
if ~isfield(params,'energy_global_mult'), params.energy_global_mult = 1.18; end   % 1.15 全局能量倍率
if ~isfield(params,'burst_thresh'),      params.burst_thresh      = 1.3; end    % 突发度门
if ~isfield(params,'use_dynamic_bg'),    params.use_dynamic_bg    = true; end  % 动态背景门
if ~isfield(params,'alpha'),             params.alpha             = 3.0; end    % 动态背景系数

if ~isfield(params,'use_weak_scoring'),  params.use_weak_scoring  = true; end  % 弱通道质量打分
if ~isfield(params,'weak_quality_cut'),  params.weak_quality_cut  = 0.4; end
if ~isfield(params,'weak_min_sep_ms'),   params.weak_min_sep_ms   = 3.0; end
% params.use_weak_scoring     = true;
% params.entropy_thresh_weak  = params.entropy_thresh + bestL.delta; % delta=0.240
% params.rate_k_weak          = 0.40;    % ← rk_w
% params.weak_quality_cut     = 0.45;    % ← cut
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
if ~isfield(params,'stalta_mode'),         params.stalta_mode       = 'ratio'; end
if ~isfield(params,'stalta_thresh'),       params.stalta_thresh     = 3.25;   end  % 3.2 ratio 模式用
if ~isfield(params,'stalta_quantile'),     params.stalta_quantile   = 0.9985;end  % quantile 模式用
if ~isfield(params,'env_smooth_ms'),       params.env_smooth_ms     = 0.03;  end  % ms
if ~isfield(params,'min_dur_ms'),          params.min_dur_ms = 0.06; end
if ~isfield(params,'min_sep_ms'),          params.min_sep_ms = 0.60; end %0.5
% SNR 物理门（MAD）
if ~isfield(params,'snr_thresh'),          params.snr_thresh        = 9.0;    end %8.8
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

% % 兼容旧签名：如果只要5个输出，就把 Coverage 当作 Layered 的
% if nargout<=5
%     Coverage = Cov_layered;
% else
%     Coverage = Cov_orig;     % 第5个输出作为 Orig 覆盖率
%     varargout{1} = Cov_layered;  % 第6个输出作为 Layered 覆盖率
% end
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

function [TOP, ALL] = sweep_small()
% 小范围参数扫掠（围绕当前较优配置）
% - 仅扫 2/6/8mm（4mm 暂不参与）
% - 先用子集样本快速评估，再对前K组合做全量复评
% - 结果自动落 CSV

%% 路径与数据源（复用你主程序里的设置）
rootBase = 'E:\Code\SCC_matlab_code\results\Entropy_Detection_Viz_updated_all';
mkData   = @(mm) fullfile(rootBase, sprintf('result%imm',mm), sprintf('result%imm',mm), 'data_denoised');
dataDirs = struct('mm2', mkData(2), 'mm6', mkData(6), 'mm8', mkData(8)); % ← 不含 4mm
thOrder  = {'mm2','mm6','mm8'};
tlabel   = struct('mm2','2mm','mm6','6mm','mm8','8mm');

outDir = fullfile(rootBase, 'eval_exports', 'sweep');
if ~exist(outDir,'dir'), mkdir(outDir); end

% 基线（与 eval_and_save_all 同步）
bestO.ent = struct('mm2',0.520,'mm6',0.520,'mm8',0.520);
bestO.tol = 1.80;
bestO.frt = -1700;
bestL.delta = 0.240;
bestL.cut   = 0.45;

% —— 搜索网格（小范围微调，可按需增减）——
G.tol_ms            = 1.80;
G.stalta_thresh     = [3.20, 3.25, 3.30];
G.min_sep_ms        = [0.50, 0.60];
G.min_dur_ms        = [0.060, 0.065];
G.env_smooth_ms     = [0.030, 0.035];
G.energy_global_mult= [1.15, 1.18];
G.snr_thresh        = [8.8, 9.0];
% Layered 相关（目前表现很好，先少动；若想再涨 AGRl 可稍放开 cut）
G.weak_quality_cut  = [0.45];   % 可改为 [0.45, 0.44]
G.rate_k_weak       = [0.40];   % 固定

% —— 快速搜索参数 —— 
Nsub     = 80;     % 子集样本数（快速评估）
seed     = 42;     % 固定随机子集
topK     = 12;     % 选出前 K 个组合进入全量复评

% ====== 生成组合 ======
grid = allcomb(G.tol_ms, G.stalta_thresh, G.min_sep_ms, G.min_dur_ms, ...
               G.env_smooth_ms, G.energy_global_mult, G.snr_thresh, ...
               G.weak_quality_cut, G.rate_k_weak);
nC = size(grid,1);
fprintf('Sweep combos: %d\n', nC);

if isempty(gcp('nocreate')), parpool; end

% ====== 快速评估（子集）======
rng(seed);
quickRows = [];
for ci = 1:nC
    % 打包参数
    p = base_params();                         % 默认值（见下方子函数）
    p.use_fixed_rate     = true;
    p.fixed_rate_thr     = bestO.frt;
    p.tol_ms             = grid(ci,1);
    p.stalta_thresh      = grid(ci,2);
    p.min_sep_ms         = grid(ci,3);
    p.min_dur_ms         = grid(ci,4);
    p.env_smooth_ms      = grid(ci,5);
    p.energy_global_mult = grid(ci,6);
    p.snr_thresh         = grid(ci,7);
    p.use_weak_scoring   = true;
    p.weak_quality_cut   = grid(ci,8);
    p.rate_k_weak        = grid(ci,9);

    % Layered 阈值按厚度设，其他共用
    metrics = struct(); % 汇总每厚度
    for t = 1:numel(thOrder)
        tk   = thOrder{t};
        ddir = dataDirs.(tk);

        p.entropy_thresh       = bestO.ent.(tk);
        p.entropy_thresh_weak  = p.entropy_thresh + bestL.delta;

        % 取子集样本
        nAll = 200;
        idx  = randperm(nAll, Nsub);

        [mAGRo, mAGRl, mEORo, mEORl, mCovO, mCovL, nOK] = eval_on_set(ddir, idx, p);
        metrics.(tk) = [mAGRo, mAGRl, mEORo, mEORl, mCovO, mCovL, nOK];
    end

    % —— 目标函数（按你的偏好）——
    % 重点：AGRo/AGRl ↑、EORo/EORl ↓，Coverage 略加分
    sc = 0;
    for t = 1:numel(thOrder)
        M = metrics.(thOrder{t});
        AGRo=M(1); AGRl=M(2); EORo=M(3); EORl=M(4); CovO=M(5); CovL=M(6);
        sc = sc ...
            + 0.40*AGRo - 0.30*EORo ...
            + 0.60*AGRl - 0.40*EORl ...
            + 0.10*((CovO+CovL)/2);
    end

    quickRows = [quickRows; ...
        ci, grid(ci,:), sc, ...
        metrics.mm2, metrics.mm6, metrics.mm8]; %#ok<AGROW>
end

% 整理快速评估结果
Vnames = {'cid','tol_ms','stalta_thresh','min_sep_ms','min_dur_ms','env_smooth_ms','energy_mult','snr_thresh','weak_cut','rk_weak','score', ...
          'm2_AGRo','m2_AGRl','m2_EORo','m2_EORl','m2_CovO','m2_CovL','m2_n', ...
          'm6_AGRo','m6_AGRl','m6_EORo','m6_EORl','m6_CovO','m6_CovL','m6_n', ...
          'm8_AGRo','m8_AGRl','m8_EORo','m8_EORl','m8_CovO','m8_CovL','m8_n'};
QUICK = array2table(quickRows,'VariableNames',Vnames);
QUICK = sortrows(QUICK,'score','descend');

% 保存快速评估
writetable(QUICK, fullfile(outDir,'sweep_quick_subset.csv'));
fprintf('Saved quick sweep: %s\n', fullfile(outDir,'sweep_quick_subset.csv'));

% ====== 全量复评（Top K）======
K = min(topK, height(QUICK));
fullRows = [];
for r = 1:K
    row = QUICK(r,:);
    % 参数
    p = base_params();
    p.use_fixed_rate     = true;
    p.fixed_rate_thr     = bestO.frt;
    p.tol_ms             = row.tol_ms;
    p.stalta_thresh      = row.stalta_thresh;
    p.min_sep_ms         = row.min_sep_ms;
    p.min_dur_ms         = row.min_dur_ms;
    p.env_smooth_ms      = row.env_smooth_ms;
    p.energy_global_mult = row.energy_mult;
    p.snr_thresh         = row.snr_thresh;
    p.use_weak_scoring   = true;
    p.weak_quality_cut   = row.weak_cut;
    p.rate_k_weak        = row.rk_weak;

    % 全量评
    met = [];
    for t = 1:numel(thOrder)
        tk   = thOrder{t};
        ddir = dataDirs.(tk);

        p.entropy_thresh       = bestO.ent.(tk);
        p.entropy_thresh_weak  = p.entropy_thresh + bestL.delta;

        idx = 1:200; % 全量
        [mAGRo, mAGRl, mEORo, mEORl, mCovO, mCovL, nOK] = eval_on_set(ddir, idx, p);
        met = [met, mAGRo, mAGRl, mEORo, mEORl, mCovO, mCovL, nOK]; %#ok<AGROW>
    end

    % 重新算分
    sc = 0;
    for t = 1:numel(thOrder)
        base = 1 + (t-1)*7;  % 每厚度 7 列
        AGRo=met(base); AGRl=met(base+1); EORo=met(base+2); EORl=met(base+3); CovO=met(base+4); CovL=met(base+5);
        sc = sc ...
            + 0.40*AGRo - 0.30*EORo ...
            + 0.60*AGRl - 0.40*EORl ...
            + 0.10*((CovO+CovL)/2);
    end

    fullRows = [fullRows; table( ...
        row.cid, row.tol_ms, row.stalta_thresh, row.min_sep_ms, row.min_dur_ms, row.env_smooth_ms, row.energy_mult, row.snr_thresh, row.weak_cut, row.rk_weak, sc, ...
        'VariableNames', {'cid','tol_ms','stalta_thresh','min_sep_ms','min_dur_ms','env_smooth_ms','energy_mult','snr_thresh','weak_cut','rk_weak','score'}) ...
        array2table(met, 'VariableNames', { ...
        'm2_AGRo','m2_AGRl','m2_EORo','m2_EORl','m2_CovO','m2_CovL','m2_n', ...
        'm6_AGRo','m6_AGRl','m6_EORo','m6_EORl','m6_CovO','m6_CovL','m6_n', ...
        'm8_AGRo','m8_AGRl','m8_EORo','m8_EORl','m8_CovO','m8_CovL','m8_n'})]; %#ok<AGROW>
end

ALL = sortrows(fullRows,'score','descend');
writetable(ALL, fullfile(outDir,'sweep_full_topK.csv'));
fprintf('Saved full sweep (Top %d): %s\n', K, fullfile(outDir,'sweep_full_topK.csv'));

% 返回前 5 条方便查看
TOP = ALL(1:min(5,height(ALL)),:);
disp(TOP);
end

% ====== 小工具：默认参数底座（与你主函数保持一致）======
function p = base_params()
p = struct();
p.use_sst             = false;
p.use_energy_gates    = true;
p.energy_global_mult  = 1.15;
p.burst_thresh        = 1.3;
p.use_dynamic_bg      = true;
p.alpha               = 3.0;
p.use_weak_scoring    = true;
p.weak_quality_cut    = 0.45;
p.weak_min_sep_ms     = 3.0;

p.entropy_thresh      = 0.38;
p.rate_k              = 2.0;
p.tol_ms              = 0.30;
p.entropy_thresh_weak = p.entropy_thresh + 0.12;
p.rate_k_weak         = 1.0;

p.use_percentile_norm = true;
p.entropy_pct_lo      = 1;
p.entropy_pct_hi      = 99;
p.rate_smooth_win     = 3;

p.use_fixed_rate      = false;
p.fixed_rate_thr      = -1900;

p.stalta_mode         = 'ratio';
p.stalta_thresh       = 3.2;
p.stalta_quantile     = 0.9985;
p.env_smooth_ms       = 0.03;
p.min_dur_ms          = 0.06;
p.min_sep_ms          = 0.50;

p.snr_thresh          = 8.8;
p.win_energy_pts      = 200;
p.bg_win_pts          = 12000;

p.min_interval_ms     = 2.6;
end

% ====== 笛卡尔积（无依赖实现）======
function M = allcomb(varargin)
% 生成所有组合的矩阵，每列对应一个向量
args = varargin;
n = numel(args);
[args{1:n}] = ndgrid(args{:});
M = zeros(numel(args{1}), n);
for i=1:n
    M(:,i) = args{i}(:);
end
end

