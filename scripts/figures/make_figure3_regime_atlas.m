function make_figure3_regime_atlas
% MAKE_FIGURE3_REGIME_ATLAS
%
% Figure 3:
% Background-only coarse qualitative regime atlas for the 1D Brusselator.
% No continuation curves are overlaid.
%
% The figure is intended for Chapter 2 (Pattern Regimes of Interest).
% Each tile shows a normalized late-time spacetime pattern at one
% parameter value (sigma,b).
%
% Requires:
%   - half_target_seed.mat
%   - solve_brusselator_1d.m
%
% Notes:
%   - Uses the target seed, since it gives a broad survey of the pattern
%     regimes of interest in the relevant region.
%   - If desired, you can later make a second version with ic_0.45_10.mat
%     and compare.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
% smaller window than the final Figure-13-style plots
sigma_min = 0.30;
sigma_max = 0.53;
b_min     = 9.80;
b_max     = 12.60;

% fewer grid cells than your Figure13 scripts
agrid = linspace(sigma_min, sigma_max, 5);   % sigma centers
bgrid = linspace(b_min, b_max, 5);           % b centers

T_relax = 240;
T_obs   = 240;
n_keep  = 100;

tile_alpha = 0.95;     % since no curves are overlaid, can keep tiles stronger

save_png = true;
png_name = figure_output_path('Figure3_regime_atlas.png');

%% ---------------- LOAD SEED ----------------
Sseed = load('half_target_seed.mat');
if isfield(Sseed,'ic')
    ic_seed = Sseed.ic;
elseif isfield(Sseed,'ic_end')
    ic_seed = Sseed.ic_end;
else
    fn = fieldnames(Sseed);
    ic_seed = Sseed.(fn{1});
end

%% ---------------- BUILD TILE GRID ----------------
figure('Color','w','Position',[100 100 950 800]);
hold on;

da = agrid(2) - agrid(1);
db = bgrid(2) - bgrid(1);

fprintf('Generating Figure 3 regime atlas: %d x %d = %d tiles\n', ...
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

%% ---------------- STYLE ----------------
colormap(parula(256));
caxis([-1.5, 1.5]);
cb = colorbar;
cb.Label.String = 'normalized activity';
cb.Label.FontSize = 14;

set(gca, 'YDir', 'normal');
box on;
grid on;

xlabel('\sigma', 'Interpreter','tex', 'FontSize',16);
ylabel('b',      'Interpreter','tex', 'FontSize',16);

set(gca, 'FontSize', 14);
xlim([sigma_min - 0.5*da, sigma_max + 0.5*da]);
ylim([b_min - 0.5*db,     b_max     + 0.5*db]);

axis square;

%% ---------------- OPTIONAL GRID LABELS ----------------
% Uncomment if you want tile-center markers:
% for ib = 1:numel(bgrid)
%     for ia = 1:numel(agrid)
%         plot(agrid(ia), bgrid(ib), 'k.', 'MarkerSize', 8);
%     end
% end

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(gcf, png_name, 'Resolution', 300);
    fprintf('Saved Figure 3 to %s\n', png_name);
end

end

%% ========================================================================
function Z = make_tile(ic_seed, sigma, b, T_relax, T_obs, n_keep)
% Simulate one parameter point and return a normalized late-time spacetime tile

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
