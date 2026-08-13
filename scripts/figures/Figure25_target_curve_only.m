function Figure25_target_curve_only(middle_dir)
% FIGURE25_TARGET_CURVE_ONLY
%
% Curve-only plot for the target-seed continuation results.
% Left family:
%   target_left_up, target_left_down, target_middle_left
% Right family:
%   target_right_up, target_right_down, target_middle_right
%
% Optional input:
%   middle_dir  directory containing newly computed
%               target_middle_left.mat and target_middle_right.mat.

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
sigma_max = 0.55;
b_min     = 9.20;
b_max     = 16.00;

lw = 2.4;
ms = 6.0;

save_png = false;
png_name = figure_output_path('Figure25_target_curve_only.png');

%% ---------------- LOAD BRANCHES ----------------
TLU = load_branch(fullfile(project_root, 'target_left_up.mat'),      'up');
TLD = load_branch(fullfile(project_root, 'target_left_down.mat'),    'down');
TRU = load_branch(fullfile(project_root, 'target_right_up.mat'),     'up');
TRD = load_branch(fullfile(project_root, 'target_right_down.mat'),   'down');
TML = load_branch(fullfile(middle_dir, 'target_middle_left.mat'),    'mid');
TMR = load_branch(fullfile(middle_dir, 'target_middle_right.mat'),   'mid');

%% ---------------- PLOT ----------------
figure('Color','w','Position',[100 80 900 760]);
hold on;

% left family: all same color/style
plot(TLD(1,:), TLD(2,:), 'ks-', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Left-Down');

plot(TLU(1,:), TLU(2,:), 'ks-', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Left-Up');

plot(TML(1,:), TML(2,:), 'ks-', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Middle-Left');

% right family: all same color/style
plot(TRD(1,:), TRD(2,:), 'bd-', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Right-Down');

plot(TRU(1,:), TRU(2,:), 'bd-', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Right-Up');

plot(TMR(1,:), TMR(2,:), 'bd-', ...
    'LineWidth', lw, 'MarkerSize', ms, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Middle-Right');

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
