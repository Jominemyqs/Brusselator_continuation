function plot_combined_curve
% PLOT_COMBINED_CURVE_ONLY
%
% Combined curve-only plot for spiral and target transition curves.
% No background tiles.
%
% Expected files (adjust names if needed):
%   Spiral:
%     right_up.mat, right_down.mat, left_up.mat, left_down.mat,
%     middle_right.mat, middle_left.mat
%
%   Target:
%     target_right_up.mat, target_right_down.mat, target_left_up.mat, target_left_down.mat,
%     target_middle_right.mat, target_middle_left.mat

clear; clc; close all;

%% ---------------- AXES WINDOW ----------------
sigma_min = 0.28;
sigma_max = 0.60;
b_min     = 9.20;
b_max     = 16.00;

lw = 2.4;
ms = 6.0;

save_png = false;
png_name = figure_output_path('combined_curve_only.png');

%% ---------------- LOAD SPIRAL ----------------
RU  = load_branch('right_up.mat',          'up');
RD  = load_branch('right_down.mat',        'down');
LU  = load_branch('left_up.mat',           'up');
LD  = load_branch('left_down.mat',         'down');
MR  = load_branch('middle_right.mat',      'mid');
ML  = load_branch('middle_left.mat',       'mid');

%% ---------------- LOAD TARGET ----------------
TRU = load_branch('target_right_up.mat',   'up');
TRD = load_branch('target_right_down.mat', 'down');
TLU = load_branch('target_left_up.mat',    'up');
TLD = load_branch('target_left_down.mat',  'down');
TMR = load_branch('target_middle_right.mat','mid');
TML = load_branch('target_middle_left.mat', 'mid');

%% ---------------- PLOT ----------------
figure('Color','w','Position',[100 80 900 750]);
hold on;

% ---------- Spiral family ----------
plot(LD(1,:), LD(2,:), 'ko-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Spiral Left-Down');
plot(LU(1,:), LU(2,:), 'ko-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'HandleVisibility', 'off');

plot(RD(1,:), RD(2,:), 'b^-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Spiral Right-Down');
plot(RU(1,:), RU(2,:), 'b^-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'HandleVisibility', 'off');

plot(ML(1,:), ML(2,:), 'ks--', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Spiral Middle-Left');
plot(MR(1,:), MR(2,:), 'bd--', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Spiral Middle-Right');

% ---------- Target family ----------
plot(TLD(1,:), TLD(2,:), 'm s:', 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor','w', 'DisplayName', 'Target Left-Down');
plot(TLU(1,:), TLU(2,:), 'm s:', 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor','w', 'HandleVisibility', 'off');

plot(TRD(1,:), TRD(2,:), 'r d:', 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor','w', 'DisplayName', 'Target Right-Down');
plot(TRU(1,:), TRU(2,:), 'r d:', 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor','w', 'HandleVisibility', 'off');

plot(TML(1,:), TML(2,:), 'm x-.', 'LineWidth', lw, 'MarkerSize', ms+1, ...
    'DisplayName', 'Target Middle-Left');
plot(TMR(1,:), TMR(2,:), 'r x-.', 'LineWidth', lw, 'MarkerSize', ms+1, ...
    'DisplayName', 'Target Middle-Right');

%% ---------------- AXES / STYLE ----------------
xlim([sigma_min, sigma_max]);
ylim([b_min, b_max]);

xlabel('\sigma', 'Interpreter', 'tex');
ylabel('b', 'Interpreter', 'tex');
set(gca, 'FontSize', 15);
box on;
grid on;
axis square;


if save_png
    exportgraphics(gcf, png_name, 'Resolution', 300);
    fprintf('Saved figure to %s\n', png_name);
end

end

%% ========================================================================
function P = load_branch(fname, direction_tag)
S = load(fname);

P = [];
if isfield(S,'p_hist_up') && strcmp(direction_tag,'up')
    P = S.p_hist_up;
end
if isfield(S,'p_hist_dn') && strcmp(direction_tag,'down')
    P = S.p_hist_dn;
end
if isempty(P) && isfield(S,'p_hist_mid') && strcmp(direction_tag,'mid')
    P = S.p_hist_mid;
end
if isempty(P) && isfield(S,'p_history')
    P = S.p_history;
end
if isempty(P) && isfield(S,'p_hist')
    P = S.p_hist;
end

if isempty(P)
    error('Could not find branch history in %s', fname);
end

if size(P,2) >= 2 && norm(P(:,1)-P(:,2)) < 1e-12
    P = P(:,2:end);
end
end
