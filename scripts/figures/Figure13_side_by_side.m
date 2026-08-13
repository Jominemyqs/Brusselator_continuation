function Figure13_side_by_side
% FIGURE13_SIDE_BY_SIDE
%
% Side-by-side Figure-13-style panels:
%   LEFT  = spiral continuation + spacetime tiles
%   RIGHT = target continuation + spacetime tiles
%
% Common axes are enforced so the two panels can be compared directly.
%
% Required files:
%   Spiral:
%       right_up.mat, right_down.mat, left_up.mat, left_down.mat, ic_0.45_10.mat
%   Target:
%       target_right_up.mat, target_right_down.mat, target_left_up.mat, target_left_down.mat
%       half_target_seed.mat
%   Solver:
%       solve_brusselator_1d.m

clear; clc; close all;

%% ---------------- COMMON DISPLAY WINDOW ----------------
sigma_min = 0.28;
sigma_max = 0.60;
b_min     = 9.20;
b_max     = 16.00;

% Keep same coarse grid idea
agrid = linspace(sigma_min, sigma_max, 8);
bgrid = linspace(b_min, b_max, 11);

da = agrid(2) - agrid(1);
db = bgrid(2) - bgrid(1);

tile_alpha = 0.50;

lw_branch = 2.2;
ms_branch = 5.8;

save_png = false;
png_name = figure_output_path('Figure13_side_by_side.png');

%% ---------------- TILE SETTINGS ----------------
% Spiral
T_relax_sp = 360;
T_obs_sp   = 360;
n_keep_sp  = 120;

% Target
T_relax_tg = 240;
T_obs_tg   = 240;
n_keep_tg  = 120;

%% ---------------- LOAD SEEDS ----------------
% Spiral seed
Sseed_sp = load('ic_0.45_10.mat');
if isfield(Sseed_sp,'ic')
    ic_seed_sp = Sseed_sp.ic;
else
    fn = fieldnames(Sseed_sp);
    ic_seed_sp = Sseed_sp.(fn{1});
end

% Target seed
Sseed_tg = load('half_target_seed.mat');
if isfield(Sseed_tg,'ic')
    ic_seed_tg = Sseed_tg.ic;
elseif isfield(Sseed_tg,'ic_end')
    ic_seed_tg = Sseed_tg.ic_end;
else
    fn = fieldnames(Sseed_tg);
    ic_seed_tg = Sseed_tg.(fn{1});
end

%% ---------------- LOAD BRANCHES ----------------
% Spiral
RU = load_branch('right_up.mat',   'up');
RD = load_branch('right_down.mat', 'down');
LU = load_branch('left_up.mat',    'up');
LD = load_branch('left_down.mat',  'down');

% Target
TRU = load_branch('target_right_up.mat',   'up');
TRD = load_branch('target_right_down.mat', 'down');
TLU = load_branch('target_left_up.mat',    'up');
TLD = load_branch('target_left_down.mat',  'down');

%% ---------------- FIGURE / LAYOUT ----------------
figure('Color','w','Position',[100 80 1500 700]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

%% =======================================================================
% LEFT PANEL: SPIRAL
ax1 = nexttile(tl,1);
hold(ax1,'on');

fprintf('Generating spiral tiles...\n');
for ib = 1:numel(bgrid)
    b0 = bgrid(ib);
    for ia = 1:numel(agrid)
        sigma0 = agrid(ia);
        fprintf('  spiral tile at (sigma,b) = (%.4f, %.4f)\n', sigma0, b0);

        Z = make_tile(ic_seed_sp, sigma0, b0, T_relax_sp, T_obs_sp, n_keep_sp);

        X = linspace(sigma0 - 0.5*da, sigma0 + 0.5*da, size(Z,2));
        Y = linspace(b0     - 0.5*db, b0     + 0.5*db, size(Z,1));

        imagesc(ax1, X, Y, Z, 'AlphaData', tile_alpha);
    end
end

% Spiral colormap
try
    C1 = brewermap([], 'RdYlBu');
    colormap(ax1, C1);
catch
    colormap(ax1, parula(256));
end
caxis(ax1, [-1.5, 1.5]);
set(ax1, 'YDir', 'normal');
box(ax1, 'on');
grid(ax1, 'on');
axis(ax1, 'square');
xlim(ax1, [sigma_min - 0.5*da, sigma_max + 0.5*da]);
ylim(ax1, [b_min - 0.5*db,     b_max     + 0.5*db]);
xlabel(ax1, '\sigma', 'Interpreter', 'tex');
ylabel(ax1, 'b',      'Interpreter', 'tex');
title(ax1, 'Spiral');
set(ax1, 'FontSize', 15);

% Overlay spiral branches
plot(ax1, LD(1,:), LD(2,:), 'ko:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Left-Down');
plot(ax1, LU(1,:), LU(2,:), 'ko:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Left-Up');

plot(ax1, RD(1,:), RD(2,:), 'b^-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Right-Down');
plot(ax1, RU(1,:), RU(2,:), 'b^-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Right-Up');

legend(ax1, 'Location', 'southoutside', 'NumColumns', 2);

%% =======================================================================
% RIGHT PANEL: TARGET
ax2 = nexttile(tl,2);
hold(ax2,'on');

fprintf('Generating target tiles...\n');
for ib = 1:numel(bgrid)
    b0 = bgrid(ib);
    for ia = 1:numel(agrid)
        sigma0 = agrid(ia);
        fprintf('  target tile at (sigma,b) = (%.4f, %.4f)\n', sigma0, b0);

        Z = make_tile(ic_seed_tg, sigma0, b0, T_relax_tg, T_obs_tg, n_keep_tg);

        X = linspace(sigma0 - 0.5*da, sigma0 + 0.5*da, size(Z,2));
        Y = linspace(b0     - 0.5*db, b0     + 0.5*db, size(Z,1));

        imagesc(ax2, X, Y, Z, 'AlphaData', tile_alpha);
    end
end

% Target colormap
colormap(ax2, turbo(256));
caxis(ax2, [-1.5, 1.5]);
set(ax2, 'YDir', 'normal');
box(ax2, 'on');
grid(ax2, 'on');
axis(ax2, 'square');
xlim(ax2, [sigma_min - 0.5*da, sigma_max + 0.5*da]);
ylim(ax2, [b_min - 0.5*db,     b_max     + 0.5*db]);
xlabel(ax2, '\sigma', 'Interpreter', 'tex');
ylabel(ax2, 'b',      'Interpreter', 'tex');
title(ax2, 'Target');
set(ax2, 'FontSize', 15);

% Overlay target branches
plot(ax2, TLD(1,:), TLD(2,:), 'ks:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Target Left-Down');
plot(ax2, TLU(1,:), TLU(2,:), 'ks:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Target Left-Up');

plot(ax2, TRD(1,:), TRD(2,:), 'bd-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Target Right-Down');
plot(ax2, TRU(1,:), TRU(2,:), 'bd-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'MarkerFaceColor', 'w', 'DisplayName', 'Target Right-Up');

legend(ax2, 'Location', 'southoutside', 'NumColumns', 2);

%% ---------------- OVERALL TITLE ----------------
title(tl, 'Spiral and target continuation curves with Figure-13-style spacetime tiles', ...
    'FontSize', 17, 'FontWeight', 'normal');

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
if isfield(S,'p_hist_up') && strcmp(direction_tag,'up')
    P = S.p_hist_up;
end
if isfield(S,'p_hist_dn') && strcmp(direction_tag,'down')
    P = S.p_hist_dn;
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
