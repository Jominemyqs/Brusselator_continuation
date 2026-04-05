function Figure13_spiral_grid
% FIGURE13_SPIRAL_GRID
%
% Wenjun-style Figure 13 analogue for the spiral continuation curves:
%   - one late-time spacetime tile per sampled grid cell
%   - overlays the four continuation branches with stronger styling
%
% Assumes these files exist in current directory:
%   right_up.mat
%   right_down.mat
%   left_up.mat
%   left_down.mat
%   ic_0.45_10.mat
%   solve_brusselator_1d.m
%
% This first version does NOT specially trim the top/bottom regions.
% You can inspect the output first, then refine.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
% Parameter window
sigma_min = 0.28;
sigma_max = 0.60;
b_min     = 9.20;
b_max     = 13.00;

% Sample grid (coarse first pass)
agrid = linspace(sigma_min, sigma_max, 8);   % sigma centers
bgrid = linspace(b_min, b_max, 11);          % b centers

% Tile simulation settings
T_relax = 360;
T_obs   = 360;       % total observation window
n_keep  = 120;       % keep only late-time rows for the tile

tile_alpha = 0.50;

% Branch styling (closer to Wenjun)
lw_branch = 2.2;
ms_branch = 5.5;

% Optional: output png
save_png = false;
png_name = 'Figure13_spiral_grid.png';

%% ---------------- LOAD SEED ----------------
Sseed = load('ic_0.45_10.mat');
if isfield(Sseed,'ic')
    ic_seed = Sseed.ic;
else
    fn = fieldnames(Sseed);
    ic_seed = Sseed.(fn{1});
end

%% ---------------- LOAD BRANCHES ----------------
RU = load_branch('right_up.mat',   'up');
RD = load_branch('right_down.mat', 'down');
LU = load_branch('left_up.mat',    'up');
LD = load_branch('left_down.mat',  'down');

%% ---------------- BUILD TILE GRID ----------------
figure('Color','w');
hold on;

da = agrid(2) - agrid(1);
db = bgrid(2) - bgrid(1);

fprintf('Generating %d x %d = %d tiles...\n', numel(agrid), numel(bgrid), numel(agrid)*numel(bgrid));

for ib = 1:numel(bgrid)
    b0 = bgrid(ib);

    for ia = 1:numel(agrid)
        sigma0 = agrid(ia);

        fprintf('  tile at (sigma,b) = (%.4f, %.4f)\n', sigma0, b0);

        Z = make_tile(ic_seed, sigma0, b0, T_relax, T_obs, n_keep);

        % map tile into this parameter cell
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
title('Spiral continuation curves');
set(gca, 'FontSize', 15);

xlim([sigma_min - 0.5*da, sigma_max + 0.5*da]);
ylim([b_min - 0.5*db,     b_max     + 0.5*db]);

axis square;

%% ---------------- OVERLAY CONTINUATION CURVES ----------------
% Left branches: black circles dotted
plot(LD(1,:), LD(2,:), 'ko:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Left-Down');
plot(LU(1,:), LU(2,:), 'ko:', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Left-Up');

% Right branches: blue triangles dash-dot
plot(RD(1,:), RD(2,:), 'b^-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Right-Down');
plot(RU(1,:), RU(2,:), 'b^-.', 'LineWidth', lw_branch, 'MarkerSize', ms_branch, ...
    'DisplayName', 'Right-Up');

legend('Location', 'bestoutside');

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(gcf, png_name, 'Resolution', 300);
    fprintf('Saved figure to %s\n', png_name);
end

end

%% ========================================================================
function P = load_branch(fname, direction_tag)
% direction_tag is only for readability

S = load(fname);

P = [];
if isfield(S,'p_hist_up')   && strcmp(direction_tag,'up'),   P = S.p_hist_up; end
if isfield(S,'p_hist_dn')   && strcmp(direction_tag,'down'), P = S.p_hist_dn; end
if isempty(P) && isfield(S,'p_history'), P = S.p_history; end
if isempty(P) && isfield(S,'p_hist'),    P = S.p_hist;    end

if isempty(P)
    error('Could not find branch history in %s', fname);
end

% Remove duplicated initial point if present
if size(P,2) >= 2 && norm(P(:,1)-P(:,2)) < 1e-12
    P = P(:,2:end);
end
end

function Z = make_tile(ic_seed, sigma, b, T_relax, T_obs, n_keep)
% Simulate one parameter point and return a normalized late-time spacetime tile

par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic_seed, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel,  par, T_obs,   0);

if isempty(V)
    Z = zeros(n_keep, size(ic_seed,1));
    return;
end

% keep late-time portion only
n_keep = min(n_keep, size(V,1));
Vlate  = V(end-n_keep+1:end, :);

% normalize like Wenjun's script
if max(Vlate(:)) - min(Vlate(:)) > 1e-8
    Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    Z = zeros(size(Vlate));
end
end
