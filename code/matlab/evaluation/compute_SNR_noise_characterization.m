%% =========================================================================
%  compute_SNR_noise_characterization.m
%  SNR estimation and low-SNR performance characterization for SCC-AE data
%
%  Why this script does NOT use peak/RMS as the reported SNR:
%    In a short AE search window, the maximum of band-limited noise is often
%    3-5 times its RMS. A peak/RMS definition therefore reports about
%    10-14 dB even when the AE waveform is visually buried. This is useful as
%    a diagnostic only, not as an honest low-SNR metric.
%
%  Reported SNR metrics:
%    - SNR_total_rms_db:
%        10*log10(P_event_window / P_noise)
%      This is near 0 dB when the event window is indistinguishable from
%      noise, and can be slightly negative due to finite-sample variation.
%
%    - SNR_excess_db:
%        10*log10((P_event_window - P_noise) / P_noise)
%      This estimates signal-only power. If P_event_window <= P_noise, the
%      excess signal power is not resolvable and the value is stored as NaN.
%
%  Outputs:
%    - SNR_summary.csv
%    - SNR_per_file_metrics.csv
%    - SNR_degradation_by_bin.csv
%    - snr_events_<thickness>.mat
%    - SNR_distribution.png / .eps
%    - SNR_degradation_curves.png / .eps
%% =========================================================================
clc; clear; close all;

%% ---- 1. Configuration --------------------------------------------------
fs             = 2e6;    % sampling rate (Hz)
bandpass_low   = 1e5;    % 100 kHz
bandpass_high  = 7e5;    % 700 kHz
bandpass_order = 6;
event_half_ms  = 0.5;    % +/-0.5 ms around event time for event energy
noise_buf_ms   = 7;      % exclude +/-7 ms around events from noise estimate

configs = {
    'PATH/TO/raw/20250108_2000kHz_2mmQ235H_1000N_1/', ...
        'PATH/TO/full_dataset/result2mm/result2mm/', ...
        '2 mm', '2mm';
    'PATH/TO/raw/20250109_2000kHz_4mmQ235H_1000N_1/', ...
        'PATH/TO/full_dataset/result4mm/result4mm/', ...
        '4 mm', '4mm';
    'PATH/TO/raw/20250110_2000kHz_6mmQ235H_1000N_1/', ...
        'PATH/TO/full_dataset/result6mm/result6mm/', ...
        '6 mm', '6mm';
    'PATH/TO/raw/20250116_2000kHz_8mmq235H_1000N_2/', ...
        'PATH/TO/full_dataset/result8mm/result8mm/', ...
        '8 mm', '8mm';
};

num_files  = 200;
root_dir   = 'PATH/TO/full_dataset/';   % <-- set to your local full-dataset root (see README)
output_dir = fullfile(root_dir, 'SNR_analysis');
eval_dir   = fullfile(root_dir, 'eval_exports');
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

event_half_pts = round(event_half_ms * 1e-3 * fs);
noise_buf_pts  = round(noise_buf_ms  * 1e-3 * fs);

%% ---- 2. Pre-compute bandpass filter -----------------------------------
[b_bp, a_bp] = butter(bandpass_order, [bandpass_low, bandpass_high] / (fs/2), 'bandpass');

%% ---- 3. Main loop ------------------------------------------------------
summary_rows = [];
per_file_rows = [];
snr_collections = cell(size(configs, 1), 1);
labels = configs(:, 3);

for ti = 1:size(configs, 1)
    data_folder   = configs{ti, 1};
    result_folder = configs{ti, 2};
    label         = configs{ti, 3};
    key           = configs{ti, 4};
    fprintf('\n===== Processing %s =====\n', label);

    snr_total_rms_all = [];
    snr_excess_all    = [];
    snr_peak_all      = [];
    event_file_id_all = [];
    noise_rms_files   = nan(num_files, 1);
    n_events_files    = nan(num_files, 1);
    file_snr_median   = nan(num_files, 1);
    file_snr_mean     = nan(num_files, 1);
    file_excess_med   = nan(num_files, 1);
    file_buried_frac  = nan(num_files, 1);

    for fid = 1:num_files
        raw_path    = fullfile(data_folder,   sprintf('Data%d.txt', fid));
        result_path = fullfile(result_folder, sprintf('Data%d', fid), 'AE_final_times_ms.txt');

        if ~exist(raw_path, 'file')
            continue;
        end

        sig_raw = load(raw_path);
        sig_raw = sig_raw(:);
        N = numel(sig_raw);

        % Remove slow offset before filtering. filtfilt gives zero-phase
        % band-limited data for both noise and event windows.
        sig_raw = sig_raw - median(sig_raw);
        sig_bp = filtfilt(b_bp, a_bp, sig_raw);

        if exist(result_path, 'file')
            ev_ms = load(result_path);
            ev_ms = ev_ms(:);
        else
            ev_ms = [];
        end
        n_ev = numel(ev_ms);

        noise_mask = true(N, 1);
        for ei = 1:n_ev
            c = round(ev_ms(ei) * 1e-3 * fs);
            lo = max(1, c - noise_buf_pts);
            hi = min(N, c + noise_buf_pts);
            noise_mask(lo:hi) = false;
        end

        noise_samples = sig_bp(noise_mask);
        if numel(noise_samples) < 1000
            noise_samples = sig_bp(1:min(N, round(0.1*fs)));
        end
        noise_power = mean(noise_samples.^2);
        noise_rms = sqrt(noise_power);
        if noise_power <= eps
            continue;
        end

        noise_rms_files(fid) = noise_rms;
        n_events_files(fid) = n_ev;

        file_snr_vals = nan(n_ev, 1);
        file_excess_vals = nan(n_ev, 1);
        for ei = 1:n_ev
            c = round(ev_ms(ei) * 1e-3 * fs);
            lo = max(1, c - event_half_pts);
            hi = min(N, c + event_half_pts);
            if hi - lo < 10
                continue;
            end

            event_samples = sig_bp(lo:hi);
            event_power = mean(event_samples.^2);

            snr_total = 10 * log10(event_power / noise_power);
            if event_power > noise_power
                snr_excess = 10 * log10((event_power - noise_power) / noise_power);
            else
                snr_excess = NaN;
            end

            % Diagnostic only. This reproduces why buried events can appear
            % to have about 11 dB SNR under a peak/RMS definition.
            snr_peak = 20 * log10(max(abs(event_samples)) / noise_rms);

            snr_total_rms_all(end+1, 1) = snr_total; %#ok<SAGROW>
            snr_excess_all(end+1, 1)    = snr_excess; %#ok<SAGROW>
            snr_peak_all(end+1, 1)      = snr_peak; %#ok<SAGROW>
            event_file_id_all(end+1, 1) = fid; %#ok<SAGROW>
            file_snr_vals(ei) = snr_total;
            file_excess_vals(ei) = snr_excess;
        end

        file_snr_median(fid) = median(file_snr_vals, 'omitnan');
        file_snr_mean(fid) = mean(file_snr_vals, 'omitnan');
        file_excess_med(fid) = median(file_excess_vals, 'omitnan');
        file_buried_frac(fid) = mean(isnan(file_excess_vals));

        if mod(fid, 20) == 0
            fprintf('  File %d/%d | noise_RMS=%.5f V | n_events=%d | median SNR_total=%.2f dB\n', ...
                fid, num_files, noise_rms, n_ev, file_snr_median(fid));
        end
    end

    fprintf('\n  [%s] Events analyzed: %d\n', label, numel(snr_total_rms_all));
    fprintf('  SNR_total_rms (dB): mean=%.2f  std=%.2f  median=%.2f  min=%.2f  max=%.2f\n', ...
        mean(snr_total_rms_all, 'omitnan'), std(snr_total_rms_all, 'omitnan'), ...
        median(snr_total_rms_all, 'omitnan'), min(snr_total_rms_all), max(snr_total_rms_all));
    fprintf('  SNR_peak diagnostic (dB): mean=%.2f  median=%.2f\n', ...
        mean(snr_peak_all, 'omitnan'), median(snr_peak_all, 'omitnan'));
    fprintf('  Unresolved excess-power fraction: %.1f%%\n', 100 * mean(isnan(snr_excess_all)));

    save(fullfile(output_dir, sprintf('snr_events_%s.mat', key)), ...
        'snr_total_rms_all', 'snr_excess_all', 'snr_peak_all', ...
        'event_file_id_all', 'noise_rms_files', 'n_events_files', ...
        'file_snr_median', 'file_snr_mean', 'file_excess_med', 'file_buried_frac');

    snr_collections{ti} = snr_total_rms_all(:);
    summary_rows(end+1, :) = [ti, numel(snr_total_rms_all), ...
        mean(noise_rms_files, 'omitnan'), std(noise_rms_files, 'omitnan'), ...
        mean(snr_total_rms_all, 'omitnan'), std(snr_total_rms_all, 'omitnan'), ...
        median(snr_total_rms_all, 'omitnan'), min(snr_total_rms_all), max(snr_total_rms_all), ...
        mean(snr_peak_all, 'omitnan'), median(snr_peak_all, 'omitnan'), ...
        100 * mean(isnan(snr_excess_all))]; %#ok<SAGROW>

    per_file_rows = [per_file_rows; ...
        [repmat(ti, num_files, 1), (1:num_files)', file_snr_median, file_snr_mean, ...
         file_excess_med, file_buried_frac, noise_rms_files, n_events_files]]; %#ok<AGROW>
end

%% ---- 4. Save SNR summaries --------------------------------------------
summary_names = {'ThicknessIdx','N_Events','NoiseFloor_mean_V','NoiseFloor_std_V', ...
    'SNR_total_rms_mean_dB','SNR_total_rms_std_dB','SNR_total_rms_median_dB', ...
    'SNR_total_rms_min_dB','SNR_total_rms_max_dB', ...
    'SNR_peak_mean_dB_diagnostic','SNR_peak_median_dB_diagnostic', ...
    'Unresolved_excess_fraction_pct'};
T = array2table(summary_rows, 'VariableNames', summary_names);
T.Thickness = labels;
writetable(T, fullfile(output_dir, 'SNR_summary.csv'));

PF = array2table(per_file_rows, 'VariableNames', ...
    {'ThicknessIdx','SampleID','SNR_total_median_dB','SNR_total_mean_dB', ...
     'SNR_excess_median_dB','BuriedFraction','NoiseFloor_RMS_V','N_Events'});
PF.Thickness = labels(PF.ThicknessIdx);
writetable(PF, fullfile(output_dir, 'SNR_per_file_metrics.csv'));
fprintf('\nSummary saved to %s\n', fullfile(output_dir, 'SNR_summary.csv'));

%% ---- 5. Link SNR with existing AGR/EOR metrics -------------------------
deg_rows = [];
bin_edges = [-Inf, -0.5, 0, 0.5, 1, 2, Inf];
bin_labels = {'<-0.5', '-0.5 to 0', '0 to 0.5', '0.5 to 1', '1 to 2', '>2'};

for ti = 1:size(configs, 1)
    key = configs{ti, 4};
    eval_path = fullfile(eval_dir, sprintf('per_sample_%s_full_eval_withCovOL.csv', key));
    if ~exist(eval_path, 'file')
        eval_path = fullfile(eval_dir, sprintf('per_sample_%s.csv', key));
    end
    if ~exist(eval_path, 'file')
        continue;
    end

    M = readtable(eval_path);
    S = PF(PF.ThicknessIdx == ti, :);
    J = innerjoin(S, M, 'Keys', 'SampleID');
    for bi = 1:numel(bin_labels)
        in_bin = J.SNR_total_median_dB > bin_edges(bi) & J.SNR_total_median_dB <= bin_edges(bi+1);
        if ~any(in_bin)
            continue;
        end

        AGR = pick_metric(J, in_bin, {'AGRl','AGR_layered','AGR','Final'});
        EOR = pick_metric(J, in_bin, {'EORl','EOR_layered','EOR'});
        COV = pick_metric(J, in_bin, {'Cov','Coverage','CovL'});
        deg_rows(end+1, :) = [ti, bi, sum(in_bin), ...
            mean(J.SNR_total_median_dB(in_bin), 'omitnan'), ...
            mean(AGR, 'omitnan'), std(AGR, 'omitnan'), ...
            mean(EOR, 'omitnan'), std(EOR, 'omitnan'), ...
            mean(COV, 'omitnan'), std(COV, 'omitnan')]; %#ok<SAGROW>
    end
end

if ~isempty(deg_rows)
    D = array2table(deg_rows, 'VariableNames', ...
        {'ThicknessIdx','SNRBinIdx','N_Samples','SNR_total_median_mean_dB', ...
         'AGR_mean','AGR_std','EOR_mean','EOR_std','Coverage_mean','Coverage_std'});
    D.Thickness = labels(D.ThicknessIdx);
    D.SNRBin = bin_labels(D.SNRBinIdx)';
    writetable(D, fullfile(output_dir, 'SNR_degradation_by_bin.csv'));
end

%% ---- 6. Plot SNR distribution -----------------------------------------
set(0, 'DefaultFigureVisible', 'off');
colors = [0.22 0.47 0.72;
          0.30 0.69 0.29;
          0.89 0.10 0.11;
          0.60 0.40 0.12];

fig = figure('Position', [100 100 560 420], 'Color', 'w', 'Visible', 'off');
ax = axes(fig);
hold(ax, 'on');
for ti = 1:4
    d = snr_collections{ti};
    if isempty(d), continue; end
    c = colors(ti,:);
    q1 = quantile(d,0.25); q3 = quantile(d,0.75);
    med = median(d); iqr_v = q3 - q1;
    wlo = max(min(d), q1 - 1.5*iqr_v);
    whi = min(max(d), q3 + 1.5*iqr_v);
    bw = 0.4;

    rectangle(ax,'Position',[ti-bw/2,q1,bw,q3-q1], ...
        'EdgeColor',c,'LineWidth',1.5,'FaceColor',[c 0.12]);
    plot(ax,[ti-bw/2,ti+bw/2],[med,med],'-','Color',c,'LineWidth',2);
    plot(ax,[ti,ti],[wlo,q1],'-','Color',c,'LineWidth',1.2);
    plot(ax,[ti,ti],[q3,whi],'-','Color',c,'LineWidth',1.2);
    plot(ax,[ti-bw/4,ti+bw/4],[wlo,wlo],'-','Color',c,'LineWidth',1.2);
    plot(ax,[ti-bw/4,ti+bw/4],[whi,whi],'-','Color',c,'LineWidth',1.2);

    nsub = min(350,numel(d));
    idx = randperm(numel(d),nsub);
    jitter = 0.12*(rand(nsub,1)-0.5);
    scatter(ax,ti+jitter,d(idx),5,c,'filled','MarkerFaceAlpha',0.22);
    plot(ax,ti,mean(d,'omitnan'),'d','MarkerSize',7, ...
        'MarkerEdgeColor',[0.2 0.2 0.2],'MarkerFaceColor',[0.2 0.2 0.2]);
    text(ax,ti,whi+0.12,sprintf('%.2f',mean(d,'omitnan')), ...
        'HorizontalAlignment','center','FontSize',9,'Color',c);
end
plot(ax,[0.5,4.5],[0,0],'k--','LineWidth',0.8);
set(ax,'XTick',1:4,'XTickLabel',labels,'FontSize',11,'FontName','Times New Roman', ...
    'XLim',[0.35 4.65],'TickDir','out','Box','off');
ylabel(ax,'SNR_{total,rms} (dB)','FontSize',12,'FontName','Times New Roman');
xlabel(ax,'Specimen thickness','FontSize',12,'FontName','Times New Roman');
grid(ax,'off');
exportgraphics(fig, fullfile(output_dir,'SNR_distribution.png'), 'Resolution', 200);
exportgraphics(fig, fullfile(output_dir,'SNR_distribution.eps'), 'ContentType', 'vector');
close(fig);

%% ---- 7. Plot degradation curves ----------------------------------------
if exist('D', 'var') && ~isempty(D)
    fig = figure('Position', [100 100 720 360], 'Color', 'w', 'Visible', 'off');
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile; hold(ax1, 'on');
    ax2 = nexttile; hold(ax2, 'on');
    for ti = 1:4
        R = D(D.ThicknessIdx == ti, :);
        if isempty(R), continue; end
        [x, ord] = sort(R.SNR_total_median_mean_dB);
        plot(ax1, x, R.AGR_mean(ord), '-o', 'Color', colors(ti,:), ...
            'LineWidth', 1.3, 'MarkerFaceColor', colors(ti,:));
        plot(ax2, x, R.EOR_mean(ord), '-o', 'Color', colors(ti,:), ...
            'LineWidth', 1.3, 'MarkerFaceColor', colors(ti,:));
    end
    format_degradation_axis(ax1, 'AGR');
    format_degradation_axis(ax2, 'EOR');
    legend(ax2, labels, 'Location', 'best', 'Box', 'off');
    exportgraphics(fig, fullfile(output_dir,'SNR_degradation_curves.png'), 'Resolution', 200);
    exportgraphics(fig, fullfile(output_dir,'SNR_degradation_curves.eps'), 'ContentType', 'vector');
    close(fig);
end

fprintf('Figures saved to %s\n', output_dir);

%% ---- 8. Noise spectral characterization --------------------------------
fprintf('\n---- Noise Spectral Characterization ----\n');
for ti = 1:size(configs,1)
    data_folder = configs{ti,1};
    label = configs{ti,3};
    sfm_vals = [];
    for fid = [1 10 50 100 150 200]
        raw_path = fullfile(data_folder, sprintf('Data%d.txt', fid));
        if ~exist(raw_path,'file'), continue; end
        sig_raw = load(raw_path);
        sig_raw = sig_raw(:) - median(sig_raw(:));
        sig_f = filtfilt(b_bp, a_bp, sig_raw);
        noise_seg = sig_f(1:min(numel(sig_f), round(0.02*fs)));
        [pxx,~] = periodogram(noise_seg,[],[],fs);
        pxx = pxx + eps;
        sfm_vals(end+1) = exp(mean(log(pxx))) / mean(pxx); %#ok<SAGROW>
    end
    fprintf('  [%s] Spectral Flatness: mean=%.4f (0=tonal, 1=white)\n', ...
        label, mean(sfm_vals, 'omitnan'));
end

fprintf('\nDone. Results in: %s\n', output_dir);

%% ---- Local helper functions --------------------------------------------
function vals = pick_metric(T, rows, candidates)
vals = nan(sum(rows), 1);
for ci = 1:numel(candidates)
    name = candidates{ci};
    if ismember(name, T.Properties.VariableNames)
        vals = T{rows, name};
        return;
    end
end
end

function format_degradation_axis(ax, yname)
set(ax, 'FontName', 'Times New Roman', 'FontSize', 11, ...
    'LineWidth', 0.8, 'TickDir', 'out', 'Box', 'off');
xlabel(ax, 'Median SNR_{total,rms} per sample (dB)', ...
    'FontName', 'Times New Roman', 'FontSize', 12);
ylabel(ax, yname, 'FontName', 'Times New Roman', 'FontSize', 12);
ylim(ax, [0 1]);
grid(ax, 'on');
end
