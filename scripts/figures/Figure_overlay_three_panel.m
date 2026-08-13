function Figure_overlay_three_panel
% FIGURE_OVERLAY_THREE_PANEL
%
% (a) Figure10.jpg
% (b) overlaid spiral/target curves, two colors only
% (c) Figure14.jpg

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
imgA = 'Figure10.jpg';
imgC = 'Figure14.jpg';

save_png = true;
png_name = figure_output_path('Figure_overlay_three_panel.png');

sigma_min = 0.28;
sigma_max = 0.60;
b_min     = 9.20;
b_max     = 16.00;

c_spiral = [0.05 0.25 0.80];
c_target = [0.85 0.15 0.15];

lw_main = 1.8;
ms_main = 4.0;

use_middle_spiral_left  = true;
use_middle_spiral_right = true;
use_middle_target_left  = true;
use_middle_target_right = true;

%% ---------------- LOAD CURVES ----------------
RU = load_branch_if_exists('right_up.mat',   'up');
RD = load_branch_if_exists('right_down.mat', 'down');
LU = load_branch_if_exists('left_up.mat',    'up');
LD = load_branch_if_exists('left_down.mat',  'down');
MR = load_branch_if_exists('middle_right.mat','mid');
ML = load_branch_if_exists('middle_left.mat', 'mid');

TRU = load_branch_if_exists('target_right_up.mat',   'up');
TRD = load_branch_if_exists('target_right_down.mat', 'down');
TLU = load_branch_if_exists('target_left_up.mat',    'up');
TLD = load_branch_if_exists('target_left_down.mat',  'down');
TMR = load_branch_if_exists('target_middle_right.mat','mid');
TML = load_branch_if_exists('target_middle_left.mat', 'mid');

%% ---------------- FIGURE LAYOUT ----------------
fig = figure('Color','w','Position',[80 80 1500 480]);
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

%% ================= PANEL (a) =================
ax1 = nexttile;
if isfile(imgA)
    A = imread(imgA);
    imshow(A, 'Parent', ax1, 'Border', 'tight');
else
    axis(ax1,'off');
    text(ax1,0.5,0.5,'Figure10.jpg not found','HorizontalAlignment','center');
end
axis(ax1,'off');
title(ax1,'(a)','FontWeight','normal','FontSize',15);

%% ================= PANEL (b) =================
ax2 = nexttile; hold(ax2,'on');

% Spiral family
plot_branch(ax2, LD, c_spiral, '-',  'o', lw_main, ms_main);
plot_branch(ax2, LU, c_spiral, '-',  'o', lw_main, ms_main);
plot_branch(ax2, RD, c_spiral, '-',  'o', lw_main, ms_main);
plot_branch(ax2, RU, c_spiral, '-',  'o', lw_main, ms_main);

if use_middle_spiral_left
    plot_branch(ax2, ML, c_spiral, '--', 'o', lw_main, ms_main);
end
if use_middle_spiral_right
    plot_branch(ax2, MR, c_spiral, '--', 'o', lw_main, ms_main);
end

% Target family
plot_branch(ax2, TLD, c_target, ':',  'd', lw_main, ms_main);
plot_branch(ax2, TLU, c_target, ':',  'd', lw_main, ms_main);
plot_branch(ax2, TRD, c_target, ':',  'd', lw_main, ms_main);
plot_branch(ax2, TRU, c_target, ':',  'd', lw_main, ms_main);

if use_middle_target_left
    plot_branch(ax2, TML, c_target, '-.', 'd', lw_main, ms_main);
end
if use_middle_target_right
    plot_branch(ax2, TMR, c_target, '-.', 'd', lw_main, ms_main);
end

xlim(ax2,[sigma_min, sigma_max]);
ylim(ax2,[b_min, b_max]);
xlabel(ax2,'\sigma','Interpreter','tex');
ylabel(ax2,'b','Interpreter','tex');
set(ax2,'FontSize',13);
box(ax2,'on');
grid(ax2,'on');
axis(ax2,'square');
title(ax2,'(b)','FontWeight','normal','FontSize',15);

%% ================= PANEL (c) =================
ax3 = nexttile;
if isfile(imgC)
    C = imread(imgC);
    imshow(C, 'Parent', ax3, 'Border', 'tight');
else
    axis(ax3,'off');
    text(ax3,0.5,0.5,'Figure14.jpg not found','HorizontalAlignment','center');
end
axis(ax3,'off');
title(ax3,'(c)','FontWeight','normal','FontSize',15);

%% ---------------- FORCE SAME PANEL BOX SIZES ----------------
drawnow;
pos2 = ax2.Position;
ax1.Position = pos2;
ax3.Position = pos2;

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved figure to %s\n', png_name);
end

end

%% ========================================================================
function P = load_branch_if_exists(fname, direction_tag)
if ~isfile(fname)
    P = [];
    return;
end

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
    warning('Could not find branch history in %s', fname);
    return;
end

if size(P,2) >= 2 && norm(P(:,1)-P(:,2)) < 1e-12
    P = P(:,2:end);
end
end

function plot_branch(ax, P, color_in, line_style, marker_style, lw, ms)
if isempty(P)
    return;
end

plot(ax, P(1,:), P(2,:), ...
    'LineStyle', line_style, ...
    'Color', color_in, ...
    'LineWidth', lw, ...
    'Marker', marker_style, ...
    'MarkerSize', ms, ...
    'MarkerIndices', 1:4:size(P,2), ...
    'MarkerFaceColor', 'w');
end
