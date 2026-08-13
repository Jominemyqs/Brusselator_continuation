function Figure20_spiral_curve_only(middle_dir)
% FIGURE20_SPIRAL_CURVE_ONLY
%
% Curve-only plot for the spiral-seed continuation results.
% Includes:
%   left_up, left_down, right_up, right_down,
%   middle_left, middle_right
%
% Expected files:
%   left_up.mat
%   left_down.mat
%   right_up.mat
%   right_down.mat
%   middle_left.mat
%   middle_right.mat
%
% Optional input:
%   middle_dir  directory containing newly computed middle_left.mat and
%               middle_right.mat. Side branches remain loaded from the
%               project root.

clc; close all;

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
if nargin < 1 || isempty(middle_dir)
    middle_dir = project_root;
end
if ~isfolder(middle_dir)
    error('Middle-results directory does not exist: %s', middle_dir);
end

%% ---------------- USER SETTINGS ----------------
sigma_min = 0.28;
sigma_max = 0.60;
b_min     = 8.20;
b_max     = 16.00;

lw = 2.4;
ms = 6.0;

save_png = false;
png_name = figure_output_path('Figure20_spiral_curve_only.png');

%% ---------------- LOAD BRANCHES ----------------
LU = load_branch(fullfile(project_root, 'left_up.mat'),    'up');
LD = load_branch(fullfile(project_root, 'left_down.mat'),  'down');
RU = load_branch(fullfile(project_root, 'right_up.mat'),   'up');
RD = load_branch(fullfile(project_root, 'right_down.mat'), 'down');
ML = load_branch(fullfile(middle_dir, 'middle_left.mat'),  'mid');
MR = load_branch(fullfile(middle_dir, 'middle_right.mat'), 'mid');

%% ---------------- PLOT ----------------
figure('Color','w','Position',[100 80 900 760]);
hold on;

% left family
plot(LD(1,:), LD(2,:), 'ko-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Left-Down');
plot(LU(1,:), LU(2,:), 'ko-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Left-Up');

% right family
plot(RD(1,:), RD(2,:), 'b^-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Right-Down');
plot(RU(1,:), RU(2,:), 'b^-', 'LineWidth', lw, 'MarkerSize', ms, ...
    'DisplayName', 'Right-Up');

% middle branches
plot(ML(1,:), ML(2,:), 'ks--', 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor','w', 'DisplayName', 'Middle-Left');
plot(MR(1,:), MR(2,:), 'bd--', 'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor','w', 'DisplayName', 'Middle-Right');

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
