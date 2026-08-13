function make_figure7_local_sweep_crossing_schematic
% MAKE_FIGURE7_LOCAL_SWEEP_CROSSING_SCHEMATIC
%
% Figure 7:
% Schematic of the local sigma-sweep corrector and threshold crossing
% selection at fixed b_pred.
%
% Intended for Chapter 3.3.

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure7_local_sweep_crossing_schematic.png');

%% ---------------- USER SETTINGS ----------------
r_thr = 0.35;

sigma_pred = 0.5120;
sigma_grid = linspace(0.5075, 0.5150, 26);

% synthetic feature profile shaped like a sharp transition
y = 0.58 ./ (1 + exp((sigma_grid - 0.5109)/0.00018)) ...
    + 0.02*exp(-((sigma_grid-0.5134)/0.00035).^2);

% downward threshold crossing
above = (y >= r_thr);
cross_idx = find(above(1:end-1) & ~above(2:end), 1, 'first');

if isempty(cross_idx)
    error('No crossing found in synthetic data.');
end

sigma_cross = 0.5*(sigma_grid(cross_idx) + sigma_grid(cross_idx+1));

% optional neighboring crossing window marker
track_left  = sigma_cross - 0.00055;
track_right = sigma_cross + 0.00055;

%% ---------------- FIGURE ----------------
fig = figure('Color','w','Position',[100 100 950 700]);
ax = axes(fig); hold(ax,'on');

% feature curve
plot(ax, sigma_grid, y, 'o-', 'LineWidth', 2, 'MarkerSize', 6, ...
    'Color', [0.0 0.35 0.8], 'MarkerFaceColor', 'w');

% threshold
yline(ax, r_thr, 'r--', 'LineWidth', 2, ...
    'Label', '$r_{\mathrm{thr}}$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');

% predictor sigma
xline(ax, sigma_pred, ':', 'LineWidth', 2, 'Color', [0.2 0.2 0.2], ...
    'Label', '$\sigma_{\mathrm{pred}}$', ...
    'Interpreter','latex', ...
    'LabelVerticalAlignment','top', ...
    'LabelHorizontalAlignment','center');

% crossing
xline(ax, sigma_cross, '--', 'LineWidth', 2.0, 'Color', [0.75 0 0], ...
    'Label', '$\sigma_{\mathrm{corr}}$', ...
    'Interpreter','latex', ...
    'LabelVerticalAlignment','middle', ...
    'LabelHorizontalAlignment','center');

% crossing bracket / tracking window
plot(ax, [track_left track_right], [0.045 0.045], 'k-', 'LineWidth', 2.0);
plot(ax, [track_left track_left],   [0.038 0.052], 'k-', 'LineWidth', 2.0);
plot(ax, [track_right track_right], [0.038 0.052], 'k-', 'LineWidth', 2.0);
text(mean([track_left track_right]) - 0.0005, 0.065, 'tracking / hysteresis window', ...
    'FontSize', 14);

% highlight crossing points
plot(ax, sigma_grid(cross_idx),   y(cross_idx),   'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');
plot(ax, sigma_grid(cross_idx+1), y(cross_idx+1), 'ko', 'MarkerSize', 8, 'MarkerFaceColor','k');



%% ---------------- STYLE ----------------
xlabel(ax, '\sigma', 'Interpreter','tex', 'FontSize',17);
ylabel(ax, 'feature value', 'FontSize',17);

set(ax, 'FontSize', 14);
box(ax, 'on');
grid(ax, 'on');

xlim(ax, [min(sigma_grid)-0.0003, max(sigma_grid)+0.0003]);
ylim(ax, [0, 0.65]);

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 7 to %s\n', png_name);
end
end
