function figure_pattern_examples_paper
% FIGURE_PATTERN_EXAMPLES_PAPER
%
% Paper/thesis-ready figure showing representative 1D Brusselator
% late-time spacetime patterns:
%   (a) wave / stripe
%   (b) target or half-target
%   (c) spiral / source-defect
%
% Uses the same simulation pipeline as continuation diagnostics:
%   relax from seed -> observe -> plot late-time spacetime window
%
% Required files:
%   solve_brusselator_1d.m
%   ic_0.45_10.mat
%   half_target_seed.mat   (for half-target panel if desired)
%
% Output:
%   Figure1_pattern_examples.png
%
% NOTE:
%   If use_half_target_panel = true, then panel (b) is plotted on [0,L]
%   using half_target_seed.mat and should be described in the caption as
%   a half-target pattern (one half of the full target on [-L,L]).
%
%   If use_half_target_panel = false, then panel (b) is generated from the
%   full-domain seed ic_0.45_10.mat and should be described as a full
%   target pattern, provided the selected parameters indeed produce one.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
save_png = true;
png_name = figure_output_path('Figure1_pattern_examples.png');

% Common simulation settings
T_relax_main   = 360;
T_obs_main     = 360;
T_relax_target = 240;
T_obs_target   = 240;

n_keep = 120;   % late-time rows shown

% ------------------------------------------------
% Choose whether panel (b) is half-target on [0,L]
% or a full target on [-L,L].
use_half_target_panel = true;

% ------------------------------------------------
% Representative parameter values
%
% (a) wave / stripe
par_wave.sigma = 0.300;
par_wave.b     = 10.800;

% (b) target / half-target
par_target.sigma = 0.430;
par_target.b     = 10.500;

% (c) spiral / source-defect
% Replace if you have a cleaner true spiral example.
par_spiral.sigma = 0.510;
par_spiral.b     = 10.000;

%% ---------------- LOAD SEEDS ----------------
% Full-domain / spiral seed
S1 = load('ic_0.45_10.mat');
if isfield(S1,'ic')
    ic_seed_main = S1.ic;
else
    fn = fieldnames(S1);
    ic_seed_main = S1.(fn{1});
end

% Half-target seed
if isfile('half_target_seed.mat')
    S2 = load('half_target_seed.mat');
    if isfield(S2,'ic')
        ic_seed_target = S2.ic;
    elseif isfield(S2,'ic_end')
        ic_seed_target = S2.ic_end;
    else
        fn = fieldnames(S2);
        ic_seed_target = S2.(fn{1});
    end
else
    warning('half_target_seed.mat not found; falling back to full-domain seed.');
    ic_seed_target = ic_seed_main;
    use_half_target_panel = false;
end

%% ---------------- GENERATE PATTERNS ----------------
% (a) wave / stripe from full-domain seed
Z_wave = make_spacetime_tile( ...
    ic_seed_main, par_wave, T_relax_main, T_obs_main, n_keep);

% (b) target / half-target
if use_half_target_panel
    Z_target = make_spacetime_tile( ...
        ic_seed_target, par_target, T_relax_target, T_obs_target, n_keep);
    panel_b_text = 'Target / half-target';
else
    Z_target = make_spacetime_tile( ...
        ic_seed_main, par_target, T_relax_main, T_obs_main, n_keep);
    panel_b_text = 'Target';
end

% (c) spiral / source-defect from full-domain seed
Z_spiral = make_spacetime_tile( ...
    ic_seed_main, par_spiral, T_relax_main, T_obs_main, n_keep);

%% ---------------- PLOT ----------------
fig = figure('Color','w','Position',[80 80 1500 500]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

% same color scaling for all panels
allvals = [Z_wave(:); Z_target(:); Z_spiral(:)];
cmax = max(abs(prctile(allvals,99)));
if cmax < 1e-8
    cmax = 1;
end

titlestr = { ...
    sprintf('(a) Wave / stripe  $(\\sigma,b)=(%.3f, %.3f)$', ...
        par_wave.sigma, par_wave.b), ...
    sprintf('(b) %s  $(\\sigma,b)=(%.3f, %.3f)$', ...
        panel_b_text, par_target.sigma, par_target.b), ...
    sprintf('(c) Spiral / source-defect  $(\\sigma,b)=(%.3f, %.3f)$', ...
        par_spiral.sigma, par_spiral.b) ...
    };

Zs = {Z_wave, Z_target, Z_spiral};

for j = 1:3
    ax = nexttile(tl,j);
    imagesc(ax, Zs{j});
    set(ax,'YDir','normal');
    axis(ax,'tight');
    axis(ax,'square');
    caxis(ax,[-cmax cmax]);

    % cleaner diverging map if available
    try
        C = brewermap([], 'RdYlBu');
        colormap(ax, C);
    catch
        colormap(ax, parula(256));
    end

    title(ax, titlestr{j}, ...
        'Interpreter','latex', ...
        'FontSize', 14, ...
        'FontWeight','normal');

    xlabel(ax, 'space index', 'FontSize', 13);

    if j == 1
        ylabel(ax, 'late-time index', 'FontSize', 13);
    else
        set(ax,'YTickLabel',[]);
    end

    set(ax, 'FontSize', 12, 'LineWidth', 0.8);
    box(ax,'on');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'normalized activity';
cb.Label.FontSize = 13;
cb.FontSize = 12;

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(fig, png_name, 'Resolution', 400);
    fprintf('Saved figure to %s\n', png_name);
end

%% ---------------- PRINT REMINDER ----------------
if use_half_target_panel
    fprintf(['\nReminder for caption/text:\n' ...
        'Panel (b) is a HALF-TARGET pattern on [0,L].\n' ...
        'In the caption, clarify that it represents one half of the full target\n' ...
        'pattern on [-L,L], with the source core at x=0.\n\n']);
else
    fprintf(['\nReminder for caption/text:\n' ...
        'Panel (b) is being generated from the full-domain seed.\n' ...
        'Please verify visually that it is indeed a clean FULL target pattern.\n\n']);
end

end

%% ========================================================================
function Z = make_spacetime_tile(ic_seed, par, T_relax, T_obs, n_keep)
% Relax to the attractor / pattern
[ic_rel, ~, ~] = solve_brusselator_1d(ic_seed, par, T_relax, 0);

% Observe spacetime
[~, ~, V] = solve_brusselator_1d(ic_rel, par, T_obs, 0);

if isempty(V)
    Z = zeros(n_keep, size(ic_seed,1));
    return;
end

% Keep only late-time part
n_keep = min(n_keep, size(V,1));
Vlate = V(end-n_keep+1:end, :);

% Normalize in same spirit as continuation/tile figures
if max(Vlate(:)) - min(Vlate(:)) > 1e-8
    Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    Z = zeros(size(Vlate));
end
end
