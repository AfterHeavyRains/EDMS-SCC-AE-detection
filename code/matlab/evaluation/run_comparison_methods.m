%% ======================== Comparison Methods for New Table ========================
% 对比实验：STA/LTA、Hilbert、改进AIC、Spectral Kurtosis、Shannon熵、EDMS
% 统一以STA/LTA为参考基准，计算AGR/EOR/Coverage
% 修复：AIC向量化、每10样本打印进度、去掉不必要的parpool
% Liu Guangyu 2026

function run_comparison_methods()

%% ===== 路径配置 =====
rootBase = 'PATH\TO\full_dataset';   % <-- set to your local full-dataset root (see README)
mkData   = @(mm) fullfile(rootBase, sprintf('result%imm',mm), sprintf('result%imm',mm), 'data_denoised');
dataDirs = struct('mm2', mkData(2), 'mm4', mkData(4), 'mm6', mkData(6), 'mm8', mkData(8));
thicks   = {'mm2','mm4','mm6','mm8'};
tlabel   = struct('mm2','2mm','mm4','4mm','mm6','6mm','mm8','8mm');
outDir   = fullfile(rootBase, 'eval_exports', 'comparison');
if ~exist(outDir,'dir'), mkdir(outDir); end

%% ===== STA/LTA基准参数（与原框架一致）=====
sp.stalta_mode   = 'ratio';
sp.stalta_thresh = 2.4;
sp.env_smooth_ms = 0.05;
sp.min_dur_ms    = 0.07;
sp.min_sep_ms    = 0.5;
sp.pre_denoised  = true;
sp.tol_ms        = 1.0;
sp.env_bg_ms     = 1.2;
sp.hilbert_k     = 4.0;
sp.hilbert_snr   = 3.0;
sp.aic_env_k     = 2.5;
sp.aic_snr       = 3.2;
sp.aic_contrast  = 0.08;
sp.sk_k          = 4.4;
sp.sk_snr        = 3.4;
sp.kg_k          = 4.6;
sp.kg_snr        = 3.6;
sp.shannon_k     = 1.8;
sp.shannon_rate_k= 1.6;
sp.shannon_snr   = 2.8;

%% ===== 对比方法定义（6行×2列 cell）=====
methods = {
    'STA/LTA (Baseline)',            @(s,f) detect_stalta(s,f,sp);
    'Hilbert Fixed-Thr',             @(s,f) detect_hilbert(s,f,sp);
    'Improved AIC (Hou 2025)',       @(s,f) detect_aic(s,f,sp);
    'Spectral Kurtosis',             @(s,f) detect_sk(s,f,sp);
    'Kurtogram',                     @(s,f) detect_kurtogram(s,f,sp);
    'Shannon Entropy (multi-scale)', @(s,f) detect_shannon(s,f,sp);
    'EDMS (Proposed)',               [];
};
nMethods = size(methods,1);
nSamples = 200;

%% ===== 读入EDMS已有结果（避免重跑）=====
edms_res = struct();
for t = 1:numel(thicks)
    tk  = thicks{t};
    Tnm = tlabel.(tk);
    csv = fullfile(rootBase,'eval_exports',sprintf('per_sample_%s.csv',Tnm));
    if exist(csv,'file')
        T = readtable(csv);
        edms_res.(tk).AGR   = T.AGRl;
        edms_res.(tk).EOR   = T.EORl;
        edms_res.(tk).Cov   = T.CovL;
        edms_res.(tk).valid = T.valid;
        fprintf('[EDMS] Loaded %s: %d samples\n', Tnm, height(T));
    else
        fprintf('[WARN] EDMS CSV not found: %s\n', csv);
        edms_res.(tk) = [];
    end
end

%% ===== 预分配结构 =====
metrics = struct();
for t = 1:numel(thicks)
    tk = thicks{t};
    for m = 1:nMethods
        mn = matlab.lang.makeValidName(methods{m,1});
        metrics.(mn).(tk).AGR = nan(nSamples,1);
        metrics.(mn).(tk).EOR = nan(nSamples,1);
        metrics.(mn).(tk).Cov = nan(nSamples,1);
    end
end

%% ===== 主循环：四厚度 =====
for t = 1:numel(thicks)
    tk   = thicks{t};
    Tnm  = tlabel.(tk);
    ddir = dataDirs.(tk);
    fprintf('\n======== %s ========\n', Tnm);

    t0 = tic;
    for k = 1:nSamples

        [sig, fs] = load_one(ddir, k);
        if isempty(sig), continue; end

        % 参考事件集（STA/LTA）
        t_ref = detect_stalta(sig, fs, sp);

        % 逐方法检测（跳过最后一个EDMS）
        for m = 1:(nMethods-1)
            mn = matlab.lang.makeValidName(methods{m,1});
            try
                t_det = methods{m,2}(sig, fs);
                [agr,eor,cov] = compute_metrics(t_det, t_ref, sp.tol_ms);
                metrics.(mn).(tk).AGR(k) = agr;
                metrics.(mn).(tk).EOR(k) = eor;
                metrics.(mn).(tk).Cov(k) = cov;
            catch ME
                fprintf('[WARN] %s #%d [%s]: %s\n', Tnm, k, methods{m,1}, ME.message);
            end
        end

        % 进度打印：第1、每10个、最后一个
        if k==1 || mod(k,10)==0 || k==nSamples
            el  = toc(t0);
            eta = el/k*(nSamples-k);
            % 打印当前各方法最新均值（只用已完成的非NaN样本）
            line = sprintf('  [%s] %3d/%d  %.0fs elapsed  ETA~%.0fs |', Tnm, k, nSamples, el, eta);
            for m = 1:(nMethods-1)
                mn  = matlab.lang.makeValidName(methods{m,1});
                agr_now = metrics.(mn).(tk).AGR(1:k);
                eor_now = metrics.(mn).(tk).EOR(1:k);
                line = sprintf('%s %s:AGR=%.2f/EOR=%.2f', line, ...
                    strtrim(methods{m,1}(1:min(8,end))), ...
                    mean(agr_now,'omitnan'), mean(eor_now,'omitnan'));
            end
            fprintf('%s\n', line);
        end
    end

    % 填入EDMS结果
    edmsN = matlab.lang.makeValidName('EDMS (Proposed)');
    if ~isempty(edms_res.(tk))
        metrics.(edmsN).(tk).AGR = edms_res.(tk).AGR;
        metrics.(edmsN).(tk).EOR = edms_res.(tk).EOR;
        metrics.(edmsN).(tk).Cov = edms_res.(tk).Cov;
    end

    % 厚度汇总打印
    fprintf('\n--- %s Summary ---\n', Tnm);
    for m = 1:nMethods
        mn = matlab.lang.makeValidName(methods{m,1});
        a  = metrics.(mn).(tk).AGR;
        e  = metrics.(mn).(tk).EOR;
        c  = metrics.(mn).(tk).Cov;
        fprintf('  %-35s AGR=%.3f+-%.3f  EOR=%.3f+-%.3f  Cov=%.3f+-%.3f\n', ...
            methods{m,1}, mean(a,'omitnan'),std(a,'omitnan'), ...
            mean(e,'omitnan'),std(e,'omitnan'), ...
            mean(c,'omitnan'),std(c,'omitnan'));
    end
end

%% ===== 输出表格 =====
write_table(metrics, methods, thicks, tlabel, outDir);
fprintf('\nDone. Outputs: %s\n', outDir);
end

%% ====================================================================
%%  方法1：STA/LTA（参考基准）
%% ====================================================================
function t_ms = detect_stalta(sig, fs, p)
time = (0:numel(sig)-1)/fs;
sig  = double(sig(:));
sta  = max(1, round(0.10e-3*fs));
lta  = max(sta+1, round(2.0e-3*fs));
env  = abs(hilbert(sig));
env  = movmean(env, max(1,round(p.env_smooth_ms*1e-3*fs)));
S    = movmean(env, sta);
L    = movmedian(env, lta);
ratio = S./(L+eps);
above = ratio > p.stalta_thresh;
d = diff([0;above(:);0]);
st = find(d==1); ed = find(d==-1)-1;
mdp = max(1,round(p.min_dur_ms*1e-3*fs));
msp = max(1,round(p.min_sep_ms*1e-3*fs));
evt = [];
for i=1:numel(st)
    seg=st(i):ed(i);
    if numel(seg)<mdp, continue; end
    [~,ix]=max(ratio(seg)); evt(end+1)=st(i)+ix-1; %#ok
end
evt = nms(evt(:), ratio, msp);
t_ms = time(evt)*1000;
end

%% ====================================================================
%%  方法2：Hilbert包络 + 固定阈值
%% ====================================================================
function t_ms = detect_hilbert(sig, fs, p)
time = (0:numel(sig)-1)/fs;
sig  = double(sig(:));
env  = abs(hilbert(sig));
env  = movmean(env, max(1,round(p.env_smooth_ms*1e-3*fs)));
thr  = robust_high_thresh(env, p.hilbert_k);
above = env > thr;
d = diff([0;above(:);0]);
st = find(d==1); ed = find(d==-1)-1;
mdp = max(1,round(p.min_dur_ms*1e-3*fs));
msp = max(1,round(p.min_sep_ms*1e-3*fs));
evt = [];
for i=1:numel(st)
    seg=st(i):ed(i);
    if numel(seg)<mdp, continue; end
    [pk,ix]=max(env(seg));
    idx = st(i)+ix-1;
    if local_peak_snr(env, idx, fs, p.env_bg_ms) >= p.hilbert_snr && pk >= thr
        evt(end+1)=idx; %#ok
    end
end
evt = nms(evt(:), env, msp);
t_ms = time(evt)*1000;
end

%% ====================================================================
%%  方法3：改进AIC（向量化累积方差，快速版）
%%  参考：Hou et al. MSSP 2025
%% ====================================================================
function t_ms = detect_aic(sig, fs, p)
time = (0:numel(sig)-1)/fs;
sig  = double(sig(:));
N    = numel(sig);

% 分析窗参数（缩短窗长减少计算量）
win_ms   = 3.0;
step_ms  = 0.5;
win_pts  = round(win_ms*1e-3*fs);
step_pts = round(step_ms*1e-3*fs);

% 包络预筛选阈值
env_raw = abs(hilbert(sig));
env_sm  = movmean(env_raw, round(0.5e-3*fs));
bg_thr  = robust_high_thresh(env_sm, p.aic_env_k);

starts = 1:step_pts:(N-win_pts+1);
evt_cand = [];

for i = 1:numel(starts)
    s = starts(i);
    e = s + win_pts - 1;

    % 预筛选：包络峰值须超阈
    if max(env_sm(s:e)) < bg_thr, continue; end

    seg = sig(s:e);
    n   = numel(seg);

    %% ---- 向量化AIC（核心加速）----
    % 利用累积均值/累积平方和，O(n)计算所有分割点的方差
    cs  = cumsum(seg);            % 累积和
    cs2 = cumsum(seg.^2);        % 累积平方和

    k_vec = (2:n-2)';            % 分割点索引向量

    % 前段 [1..k] 的方差
    mu1  = cs(k_vec) ./ k_vec;
    var1 = cs2(k_vec)./k_vec - mu1.^2;
    var1 = max(var1, eps);

    % 后段 [k+1..n] 的方差
    n2   = n - k_vec;
    sum2 = cs(n) - cs(k_vec);
    sq2  = cs2(n) - cs2(k_vec);
    mu2  = sum2 ./ n2;
    var2 = sq2./n2 - mu2.^2;
    var2 = max(var2, eps);

    aic_vec = k_vec.*log(var1) + (n2-1).*log(var2);

    [aic_min, idx] = min(aic_vec);
    kmin = k_vec(idx);
    edge_guard = max(3, round(0.08*n));
    if kmin <= edge_guard || kmin >= (n-edge_guard)
        continue;
    end
    aic_med = median(aic_vec);
    aic_scale = 1.4826*median(abs(aic_vec - aic_med)) + eps;
    if (aic_med - aic_min) < p.aic_contrast * max(aic_scale, 1)
        continue;
    end

    % 换算全局索引
    idx_glb = s + kmin - 1;
    if local_peak_snr(env_sm, idx_glb, fs, p.env_bg_ms) >= p.aic_snr
        evt_cand(end+1) = idx_glb; %#ok
    end
end

if isempty(evt_cand), t_ms = []; return; end

msp = round(p.min_sep_ms*1e-3*fs);
evt_cand = nms(evt_cand(:), env_sm, msp);
t_ms = time(evt_cand)*1000;
end

%% ====================================================================
%%  方法4：Spectral Kurtosis
%% ====================================================================
function t_ms = detect_sk(sig, fs, p)
time = (0:numel(sig)-1)/fs;
sig  = double(sig(:));
N    = numel(sig);

% 谱峭度：STFT功率谱沿时间轴的峭度
win_sk = 256; hop_sk = 64;
[S,F,~] = spectrogram(sig, hamming(win_sk), win_sk-hop_sk, win_sk, fs);
P = abs(S).^2;

nF = size(P,1);
sk_vec = zeros(nF,1);
for fi = 1:nF
    row = P(fi,:);
    mu2 = mean(row.^2);
    if mu2 > 0
        sk_vec(fi) = mean(row.^4)/(mu2^2) - 2;
    end
end

% 高SK频带
sk_thr = mean(sk_vec) + 1.5*std(sk_vec);
hf = F(sk_vec > sk_thr);
if isempty(hf)
    [~,ord] = sort(sk_vec,'descend');
    hf = F(ord(1:max(1,round(0.2*nF))));
end

f_lo = max(200e3, min(hf));
f_hi = min(700e3, max(hf));
if f_lo >= f_hi, f_lo=200e3; f_hi=700e3; end
Wp = [f_lo, f_hi]/(fs/2);
Wp = max(0.01, min(0.99, Wp));

try
    [b,a] = butter(4, Wp, 'bandpass');
    sf = filtfilt(b,a,sig);
catch
    sf = sig;
end

env = abs(hilbert(sf));
env = movmean(env, max(1,round(0.2e-3*fs)));
thr = robust_high_thresh(env, p.sk_k);

above = env > thr;
d = diff([0;above(:);0]);
st = find(d==1); ed = find(d==-1)-1;
mdp = max(1,round(p.min_dur_ms*1e-3*fs));
msp = max(1,round(p.min_sep_ms*1e-3*fs));

evt = [];
for i=1:numel(st)
    seg=st(i):ed(i);
    if numel(seg)<mdp, continue; end
    [pk,ix]=max(env(seg));
    idx = st(i)+ix-1;
    if local_peak_snr(env, idx, fs, p.env_bg_ms) >= p.sk_snr && pk >= thr
        evt(end+1)=idx; %#ok
    end
end
evt = nms(evt(:), env, msp);
t_ms = time(evt)*1000;
end

%% ====================================================================
%%  方法5：Kurtogram风格子带搜索
%% ====================================================================
function t_ms = detect_kurtogram(sig, fs, p)
time = (0:numel(sig)-1)/fs;
sig  = double(sig(:));

fc_list = [250e3, 350e3, 450e3, 550e3, 650e3];
bw_list = [80e3, 120e3, 180e3, 240e3];
best_score = -inf;
best_env = [];

for ic = 1:numel(fc_list)
    for ib = 1:numel(bw_list)
        fc = fc_list(ic);
        bw = bw_list(ib);
        flo = max(80e3, fc - bw/2);
        fhi = min(850e3, fc + bw/2);
        if flo >= fhi || fhi >= 0.98*(fs/2)
            continue;
        end
        try
            [b,a] = butter(4, [flo, fhi]/(fs/2), 'bandpass');
            sf = filtfilt(b, a, sig);
        catch
            continue;
        end
        env_tmp = abs(hilbert(sf));
        env_tmp = movmean(env_tmp, max(1, round(0.15e-3*fs)));
        score = kurtosis(env_tmp) + 0.15*log(var(env_tmp)+eps);
        if isfinite(score) && score > best_score
            best_score = score;
            best_env = env_tmp;
        end
    end
end

if isempty(best_env)
    t_ms = [];
    return;
end

env = best_env;
thr = robust_high_thresh(env, p.kg_k);
above = env > thr;
d = diff([0; above(:); 0]);
st = find(d==1);
ed = find(d==-1)-1;
mdp = max(1,round(p.min_dur_ms*1e-3*fs));
msp = max(1,round(p.min_sep_ms*1e-3*fs));

evt = [];
for i = 1:numel(st)
    seg = st(i):ed(i);
    if numel(seg) < mdp, continue; end
    [pk,ix] = max(env(seg));
    idx = st(i)+ix-1;
    if local_peak_snr(env, idx, fs, p.env_bg_ms) >= p.kg_snr && pk >= thr
        evt(end+1) = idx; %#ok
    end
end

evt = nms(evt(:), env, msp);
t_ms = time(evt)*1000;
end

%% ====================================================================
%%  方法6：多尺度Shannon熵（entropy-driven对比）
%% ====================================================================
function t_ms = detect_shannon(sig, fs, p)
time   = (0:numel(sig)-1)/fs;
sig    = double(sig(:));
win_pts  = floor(200e-6*fs);
step_pts = floor(100e-6*fs);
starts   = 1:step_pts:(numel(sig)-win_pts+1);
centers  = starts + floor(win_pts/2);
time_ax  = time(centers);
W = numel(starts);
H = zeros(1,W);

spec_win = hamming(floor(win_pts/4));
for k = 1:W
    seg = sig(starts(k):starts(k)+win_pts-1);
    [Sk,~] = spectrogram(seg, spec_win, [], [], fs);
    Pc = mean(abs(Sk).^2, 2);
    Pc = Pc/(sum(Pc)+eps);
    H(k) = -sum(Pc.*log(Pc+eps));
end

Hn = (H - min(H))/(max(H)-min(H)+eps);
Hs = movmedian(Hn, 5);
Hr = diff(Hs);
rate_med = median(Hr);
rate_mad = 1.4826*median(abs(Hr - rate_med)) + eps;
thr_lo = robust_low_thresh(Hs, p.shannon_k);
thr_rate = rate_med - p.shannon_rate_k * rate_mad;

is_below = Hs < thr_lo;
d = diff([0, is_below, 0]);
st = find(d==1); ed2 = find(d==-1)-1;
msp_s = p.min_sep_ms*1e-3;
env = abs(hilbert(sig));
env = movmean(env, max(1,round(0.2e-3*fs)));

evt_t = [];
evt_score = [];
for i=1:numel(st)
    seg_i = st(i):ed2(i);
    [hmin,ix] = min(Hs(seg_i));
    idx_h = st(i)+ix-1;
    rL = max(1, idx_h-2);
    rR = min(numel(Hr), idx_h+2);
    if min(Hr(rL:rR)) > thr_rate
        continue;
    end
    samp_idx = centers(idx_h);
    pk_snr = local_peak_snr(env, samp_idx, fs, p.env_bg_ms);
    if pk_snr < p.shannon_snr
        continue;
    end
    evt_t(end+1) = time_ax(idx_h); %#ok
    evt_score(end+1) = (thr_lo - hmin) + 0.15*max(0, pk_snr - p.shannon_snr); %#ok
end

% 时间域NMS
if isempty(evt_t)
    t_ms = [];
    return;
end
[evt_t, ord] = sort(evt_t);
evt_score = evt_score(ord);
keep  = true(size(evt_t));
for i=2:numel(evt_t)
    if (evt_t(i)-evt_t(i-1)) < msp_s
        if evt_score(i) >= evt_score(i-1)
            keep(i-1) = false;
        else
            keep(i) = false;
        end
    end
end
t_ms = evt_t(keep)*1000;
end

%% ====================================================================
%%  NMS辅助（保留幅值大的）
%% ====================================================================
function evt = nms(evt, amp, min_sep_pts)
if isempty(evt), return; end
evt = sort(evt(:));
keep = true(size(evt));
for i = 2:numel(evt)
    if (evt(i)-evt(i-1)) < min_sep_pts
        if amp(evt(i)) >= amp(evt(i-1))
            keep(i-1) = false;
        else
            keep(i) = false;
        end
    end
end
evt = evt(keep);
end

function thr = robust_high_thresh(x, k)
x = x(:);
medx = median(x,'omitnan');
madx = 1.4826*median(abs(x - medx),'omitnan') + eps;
thr = medx + k*madx;
end

function thr = robust_low_thresh(x, k)
x = x(:);
medx = median(x,'omitnan');
madx = 1.4826*median(abs(x - medx),'omitnan') + eps;
thr = medx - k*madx;
end

function snr_pk = local_peak_snr(env, idx, fs, bg_ms)
idx = max(1, min(numel(env), round(idx)));
bg_pts = max(3, round(bg_ms*1e-3*fs));
Lbg = max(1, idx - bg_pts);
Rbg = min(numel(env), idx + bg_pts);
bg = env(Lbg:Rbg);
bg_med = median(bg);
bg_mad = 1.4826*median(abs(bg - bg_med)) + eps;
snr_pk = (env(idx) - bg_med)/bg_mad;
end

%% ====================================================================
%%  评估指标（严格对应论文公式9）
%% ====================================================================
function [AGR, EOR, Cov] = compute_metrics(t_det, t_ref, tol_ms)
t_det = t_det(:); t_ref = t_ref(:);
if isempty(t_det) || isempty(t_ref)
    AGR=0; EOR=1; Cov=0; return;
end
used = false(numel(t_ref),1);
n_match = 0;
for i = 1:numel(t_det)
    [dmin,j] = min(abs(t_ref - t_det(i)));
    if dmin<=tol_ms && ~used(j)
        n_match=n_match+1; used(j)=true;
    end
end
AGR = n_match/max(1,numel(t_det));
EOR = (numel(t_det)-n_match)/max(1,numel(t_det));
Cov = n_match/max(1,numel(t_ref));
end

%% ====================================================================
%%  数据加载（复用原有逻辑）
%% ====================================================================
function [sig,fs] = load_one(dataFolder, i)
sig=[]; fs=[];
f = fullfile(dataFolder, sprintf('data%d.mat',i));
if ~exist(f,'file'), return; end
S = load(f);
flds = {'signal_denoised','sig_deno','signal','sig'};
for fi=1:numel(flds)
    if isfield(S,flds{fi}), sig=S.(flds{fi}); break; end
end
if isfield(S,'fs'), fs=S.fs; else, fs=2e6; end
sig=single(sig(:));
end

%% ====================================================================
%%  输出CSV + LaTeX
%% ====================================================================
function write_table(metrics, methods, thicks, tlabel, outDir)
nMethods = size(methods,1);
mKeys = cellfun(@(x) matlab.lang.makeValidName(x), methods(:,1), 'UniformOutput',false);

%% CSV
fid = fopen(fullfile(outDir,'Table_comparison.csv'),'w');
fprintf(fid,'Method,2mm_AGR,2mm_AGR_std,2mm_EOR,2mm_EOR_std,2mm_Cov,2mm_Cov_std,4mm_AGR,4mm_AGR_std,4mm_EOR,4mm_EOR_std,4mm_Cov,4mm_Cov_std,6mm_AGR,6mm_AGR_std,6mm_EOR,6mm_EOR_std,6mm_Cov,6mm_Cov_std,8mm_AGR,8mm_AGR_std,8mm_EOR,8mm_EOR_std,8mm_Cov,8mm_Cov_std\n');
for m=1:nMethods
    line = methods{m,1};
    for t=1:numel(thicks)
        tk=thicks{t};
        a=metrics.(mKeys{m}).(tk).AGR;
        e=metrics.(mKeys{m}).(tk).EOR;
        c=metrics.(mKeys{m}).(tk).Cov;
        line=sprintf('%s,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f',line,...
            mean(a,'omitnan'),std(a,'omitnan'),...
            mean(e,'omitnan'),std(e,'omitnan'),...
            mean(c,'omitnan'),std(c,'omitnan'));
    end
    fprintf(fid,'%s\n',line);
end
fclose(fid);
fprintf('Saved: %s\n', fullfile(outDir,'Table_comparison.csv'));

%% LaTeX
fid2 = fopen(fullfile(outDir,'Table_comparison.tex'),'w');
fprintf(fid2,'\\begin{table*}[t]\\centering\\small\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid2,'\\caption{Detection performance comparison (mean$\\pm$std). $\\uparrow$ higher is better; $\\downarrow$ lower is better.}\n');
fprintf(fid2,'\\label{tab:comparison}\n');
fprintf(fid2,'\\begin{tabular}{l|ccc|ccc|ccc|ccc}\\hline\n');
fprintf(fid2,'Method & \\multicolumn{3}{c|}{2\\,mm} & \\multicolumn{3}{c|}{4\\,mm} & \\multicolumn{3}{c|}{6\\,mm} & \\multicolumn{3}{c}{8\\,mm}\\\\\n');
fprintf(fid2,' & AGR$\\uparrow$ & EOR$\\downarrow$ & Cov$\\uparrow$ & AGR$\\uparrow$ & EOR$\\downarrow$ & Cov$\\uparrow$ & AGR$\\uparrow$ & EOR$\\downarrow$ & Cov$\\uparrow$ & AGR$\\uparrow$ & EOR$\\downarrow$ & Cov$\\uparrow$\\\\\\hline\n');
for m=1:nMethods
    isEDMS = contains(methods{m,1},'EDMS');
    if isEDMS, fprintf(fid2,'\\hline\n'); end
    nm_tex = strrep(methods{m,1},'_','\_');
    if isEDMS
        fprintf(fid2,'\\textbf{%s}', nm_tex);
    else
        fprintf(fid2,'%s', nm_tex);
    end
    for t=1:numel(thicks)
        tk=thicks{t};
        a=metrics.(mKeys{m}).(tk).AGR;
        e=metrics.(mKeys{m}).(tk).EOR;
        c=metrics.(mKeys{m}).(tk).Cov;
        am=mean(a,'omitnan'); as_=std(a,'omitnan');
        em=mean(e,'omitnan'); es_=std(e,'omitnan');
        cm=mean(c,'omitnan'); cs_=std(c,'omitnan');
        if isEDMS
            fprintf(fid2,' & \\textbf{%.3f}$\\pm$%.3f & \\textbf{%.3f}$\\pm$%.3f & \\textbf{%.3f}$\\pm$%.3f',am,as_,em,es_,cm,cs_);
        else
            fprintf(fid2,' & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f & %.3f$\\pm$%.3f',am,as_,em,es_,cm,cs_);
        end
    end
    fprintf(fid2,'\\\\\n');
end
fprintf(fid2,'\\hline\n\\end{tabular}\n\\end{table*}\n');
fclose(fid2);
fprintf('Saved: %s\n', fullfile(outDir,'Table_comparison.tex'));
end
