function Figure21_spiral_background
% FIGURE21_SPIRAL_BACKGROUND
%
% Coarse late-time spacetime background with spiral continuation curves.
% Includes:
%   left_up, left_down, right_up, right_down,
%   middle_left, middle_right
%
% Expected files:
%   right_up.mat
%   right_down.mat
%   left_up.mat
%   left_down.mat
%   middle_left.mat
%   middle_right.mat
%   ic_0.45_10.mat
%   solve_brusselator_1d.m

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
sigma_min = 0.28;
sigma_max = 0.55;
b_min     = 9.2;
b_max     = 16.00;

agrid = linspace(sigma_min, sigma_max, 8);
bgrid = linspace(b_min, b_max, 11);

T_relax = 360;
T_obs   = 360;
n_keep  = 120;

tile_alpha = 0.50;

lw_branch = 2.2;
ms_branch = 5.5;

save_png = false;
png_name = figure_output_path('Figure21_spiral_background.png');

%% ---------------- LOAD SEED ----------------
Sseed = load('ic_0.45_10.mat');
if isfield(Sseed,'ic')
    ic_seed = Sseed.ic;
else
    fn = fieldnames(Sseed);
    ic_seed = Sseed.(fn{1});
end

%% ---------------- LOAD BRANCHES ----------------
RU = load_branch('right_up.mat',     'up');
RD = load_branch('right_down.mat',   'down');
LU = load_branch('left_up.mat',      'up');
LD = load_branch('left_down.mat',    'down');
ML = load_branch('middle_left.mat',  'mid');
MR = load_branch('middle_right.mat', 'mid');

%% ---------------- BUILD TILE GRID ----------------
figure('Color','w');
hold on;

da = agrid(2) - agrid(1);
db = bgrid(2) - bgrid(1);

fprintf('Generating %d x %d = %d tiles...\n', ...
    numel(agrid), numel(bgrid), numel(agrid)*numel(bgrid));

for ib = 1:numel(bgrid)
    b0 = bgrid(ib);

    for ia = 1:numel(agrid)
        sigma0 = agrid(ia);

        fprintf('  tile at (sigma,b) = (%.4f, %.4f)\n', sigma0, b0);

        Z = make_tile(ic_seed, sigma0, b0, T_relax, T_obs, n_keep);

        X = linspace(sigma0 - 0.5*da, sigma0 + 0.5*da, size(Z,2));
        Y = linspace(b0     - 0.5*db, b0     + 0.5*db, size(Z,1));

        imagesc(X, Y, Z, 'AlphaData', tile_alpha);
        hold on;
    end
end

%% ---------------- COLORMAP / AXES ----------------
try
    C = brewermap([], 'RdYlBu');
    colormap(C);
catch
    colormap(parula(256));
end

caxis([-1.5, 1.5]);
set(gca, 'YDir', 'normal');
box on;
grid on;
xlabel('\sigma', 'Interpreter', 'tex');
ylabel('b',      'Interpreter', 'tex');
set(gca, 'FontSize', 15);

xlim([sigma_min - 0.5*da, sigma_max + 0.5*da]);
ylim([b_min - 0.5*db,     b_max     + 0.5*db]);

axis square;

%% ---------------- OVERLAY CONTINUATION CURVES ----------------
% Left branches
plot(LD(1,:), LD(2,:), 'ko:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Left-Down');
plot(LU(1,:), LU(2,:), 'ko:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Left-Up');

% Right branches
plot(RD(1,:), RD(2,:), 'b^-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Right-Down');
plot(RU(1,:), RU(2,:), 'b^-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Right-Up');

% Middle branches
plot(ML(1,:), ML(2,:), 'ks--', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Middle-Left');
plot(MR(1,:), MR(2,:), 'bd--', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Middle-Right');


%% ---------------- SAVE ----------------
if save_png
    exportgraphics(gcf, png_name, 'Resolution', 300);
    fprintf('Saved figure to %s\n', png_name);
end

end

%% ========================================================================
function P = load_branch(fname, direction_tag)
S = load(fname);

P = [];
if isfield(S,'p_hist_up')   && strcmp(direction_tag,'up')
    P = S.p_hist_up;
end
if isfield(S,'p_hist_dn')   && strcmp(direction_tag,'down')
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

function Z = make_tile(ic_seed, sigma, b, T_relax, T_obs, n_keep)
par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic_seed, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel,  par, T_obs,   0);

if isempty(V)
    Z = zeros(n_keep, size(ic_seed,1));
    return;
end

n_keep = min(n_keep, size(V,1));
Vlate  = V(end-n_keep+1:end, :);

if max(Vlate(:)) - min(Vlate(:)) > 1e-8
    Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    Z = zeros(size(Vlate));
end
end
