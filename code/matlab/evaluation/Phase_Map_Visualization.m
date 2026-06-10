% 这段代码生成的是论文中的图8，热力分布图。
% %% === Global Entropy–Rate Phase Map Visualization (Enhanced Color Scheme) ===
% clear; clc;
% load('phase_data_mm4.mat');   % ← 载入 all_H, all_R, all_cls
% 
% figure('Color','w','Position',[300 200 850 620]);
% 
% %% --- 配色方案（清新平衡型） ---
% C = struct(...
%     'HC', [0.12 0.55 0.30], ...   % 深绿（高置信）
%     'MC', [0.95 0.75 0.15], ...   % 金黄（中置信）
%     'LC', [0.25 0.45 0.85], ...   % 蓝色（低置信）
%     'BG', [0.55 0.55 0.55]);      % 背景灰
% 
% %% --- 背景柔雾点云 ---
% scatter(all_H, all_R, 6, C.BG, 'filled', ...
%     'MarkerFaceAlpha', 0.025, 'MarkerEdgeAlpha', 0.01);
% hold on;
% 
% %% --- 防堆叠轻微抖动 ---
% jitter = @(n,scale) (randn(n,1)*scale);
% 
% %% --- 半透明柔点绘制函数 ---
% plot_soft = @(mask, col) scatter(...
%     all_H(mask) + jitter(sum(mask), 0.002), ...
%     all_R(mask) + jitter(sum(mask), 40), ...
%     28, ...
%     'MarkerFaceColor', col, ...
%     'MarkerEdgeColor', [0.15 0.15 0.15], ...
%     'MarkerEdgeAlpha', 0.25, ...
%     'MarkerFaceAlpha', 0.6);
% 
% %% --- 绘制三类事件 ---
% plot_soft(strcmp(all_cls,'HC'), C.HC);
% plot_soft(strcmp(all_cls,'MC'), C.MC);
% plot_soft(strcmp(all_cls,'LC'), C.LC);
% 
% %% --- 坐标轴与网格 ---
% xlabel('Normalized Entropy (H_t)', 'FontName','Times New Roman','FontSize',12);
% ylabel('Entropy Rate (dH/dt)', 'FontName','Times New Roman','FontSize',12);
% grid on; box on;
% ax = gca;
% ax.GridLineStyle = ':';
% ax.LineWidth = 0.8;
% ax.XLim = [0 1];
% ax.YLim = [-4000 4000];
% ax.FontName = 'Times New Roman';
% ax.FontSize = 11;
% 
% %% --- 背景对角虚线（淡化） ---
% % plot([0 1], [4000 -4000], '-', 'Color', [0.7 0.7 0.7 0.12], 'LineWidth', 2.5);
% % plot([0 1], [-4000 4000], '-', 'Color', [0.7 0.7 0.7 0.12], 'LineWidth', 2.5);
% % yline(0, '-', 'Color', [0.6 0.6 0.6 0.15], 'LineWidth', 2);
% 
% %% --- 图例（框内右上角） ---
% lgd = legend({'Background','HC','MC','LC'}, ...
%     'Box','off','Location','northeast');
% lgd.Title.String = 'Event Category';
% set(lgd, 'FontSize',10, 'Color','none', 'FontName','Times New Roman');
% 
% % 轻微上移图例，避免与数据点重叠
% lgd.Position(2) = lgd.Position(2) + 0.02;  
% uistack(lgd,'top');  % 确保图例在最上层
% 
% %% --- 导出高分辨率图像 ---
% exportgraphics(gcf, 'phase_softscatter_mm4_colorEnhanced_inLegend.png', 'Resolution', 600);
% fprintf('✅ Saved: phase_softscatter_mm4_colorEnhanced_inLegend.png\n');
clear; clc;
load('phase_data_mm4.mat');

set(groot, 'DefaultAxesFontName', 'Times New Roman', ...
    'DefaultAxesFontSize', 17, ...
    'DefaultTextFontName', 'Times New Roman', ...
    'DefaultTextFontSize', 17, ...
    'DefaultAxesLineWidth', 0.9, ...
    'DefaultLineLineWidth', 0.8);

edgesH = 0:0.02:1.0;
edgesR = linspace(-4000,4000,120);
cats = {'HC','MC','LC'};
titles = {'HC Density Map','MC Density Map','LC Density Map'};

figure('Color','w','Position',[300 150 960 600]);

% === 创建均匀分布的tile布局 ===
tiledlayout(1,3, 'Padding','compact', 'TileSpacing','compact');

for i = 1:3
    mask = strcmp(all_cls, cats{i});
    N = histcounts2(all_H(mask), all_R(mask), edgesH, edgesR);

    % 平滑
    if strcmp(cats{i},'LC'), N = imgaussfilt(N, 1.0);
    else, N = imgaussfilt(N, 0.6);
    end

    % 归一化+对数增强
    vmax = prctile(N(:),95);
    N = log10(N./vmax + 1e-3);
    N(N<-3) = -3;  % 防止极暗值

    nexttile(i);
    imagesc(edgesH(1:end-1), edgesR(1:end-1), N'); axis xy;
    colormap(parula);

    % 只在最后一个tile添加colorbar
    if i == 3
        cb = colorbar;
        cb.Layout.Tile = 'east';  % 放在整个tile右侧，不占用子图空间
        cb.Label.String = 'log-scaled relative density';
        cb.Label.FontSize = 15;
    end

    % 动态clim
    if strcmp(cats{i},'HC'), clim([-2.8 0]);
    elseif strcmp(cats{i},'MC'), clim([-2.5 0]);
    else, clim([-2.2 0]);
    end

    hold on;
    contour(edgesH(1:end-1), edgesR(1:end-1), N', 5, 'k', 'LineWidth',0.45);

    xlabel('$H_t$', 'Interpreter','latex', 'FontSize',18);
    ylabel('$dH/dt$', 'Interpreter','latex', 'FontSize',18);
    title(['(' char('a'+i-1) ') ' titles{i}], 'FontSize',18);
    xlim([0 1]); ylim([-4000 4000]);
    box on;
end

%% === 导出高分辨率图像 ===
exportgraphics(gcf,'EntropyRate_MultiView_mm4_3panel.tiff','Resolution',650);
fprintf('Saved: EntropyRate_MultiView_mm4_3panel.tiff\n');

% sgtitle('Entropy–Entropy-Rate Multi-View Representation','FontName','Times New Roman','FontSize',14);

%% === 导出高分辨率图像 ===
exportgraphics(gcf,'EntropyRate_MultiView_mm4.png','Resolution',600);
fprintf('Saved: EntropyRate_MultiView_mm4.png\n');
