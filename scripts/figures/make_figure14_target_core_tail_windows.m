function make_figure14_target_core_tail_windows
% MAKE_FIGURE14_TARGET_CORE_TAIL_WINDOWS
%
% Figure 14:
% Target/source-defect spacetime plot with core and tail windows overlaid.
%
% Intended for Section 5.2.

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure14_target_core_tail_windows.png');

%% ---------------- USER SETTINGS ----------------
sigma = 0.430;
b     = 10.500;

T_relax = 240;
T_obs   = 240;
n_keep  = 120;

% Window definitions in spatial index coordinates
% Adjust to match your actual feature_evaluation_target.m choices
Nx_est = 256;   % only used to express intended placement
core_x1 = 1;
core_x2 = 70;

tail_x1 = 180;
tail_x2 = 256;

%% ---------------- LOAD TARGET SEED ----------------
S = load('half_target_seed.mat');
if isfield(S,'ic')
    ic0 = S.ic;
elseif isfield(S,'ic_end')
    ic0 = S.ic_end;
else
    fn = fieldnames(S);
    ic0 = S.(fn{1});
end

%% ---------------- SIMULATE ----------------
Vlate = get_late_spacetime(ic0, sigma, b, T_relax, T_obs, n_keep);
Z = normalize_tile(Vlate);
Nx = size(Z,2);

% rescale window endpoints if Nx differs from Nx_est
scale = Nx / Nx_est;
cx1 = max(1, round(core_x1 * scale));
cx2 = min(Nx, round(core_x2 * scale));
tx1 = max(1, round(tail_x1 * scale));
tx2 = min(Nx, round(tail_x2 * scale));

%% ---------------- PLOT ----------------
fig = figure('Color','w','Position',[100 80 950 720]);
ax = axes(fig);

imagesc(ax, Z);
set(ax,'YDir','normal');
axis(ax,'tight');
colormap(ax, parula(256));
cb = colorbar(ax);
cb.Label.String = 'normalized activity';
cb.Label.FontSize = 14;

hold(ax,'on');

% rectangles: [x y w h]
rectangle(ax, 'Position', [cx1, 1, cx2-cx1, size(Z,1)-1], ...
    'EdgeColor', 'w', 'LineWidth', 2.5, 'LineStyle', '-');
rectangle(ax, 'Position', [tx1, 1, tx2-tx1, size(Z,1)-1], ...
    'EdgeColor', [1 0.8 0], 'LineWidth', 2.5, 'LineStyle', '--');

text(ax, cx1+4, 10, 'core window', 'Color','w', 'FontSize',16, 'FontWeight','bold');
text(ax, tx1+4, 10, 'tail window', 'Color',[1 0.9 0], 'FontSize',16, 'FontWeight','bold');

hold(ax,'off');

xlabel(ax,'space index','FontSize',16);
ylabel(ax,'late-time index','FontSize',16);

set(ax,'FontSize',14);

if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 14 to %s\n', png_name);
end

end

%% ========================================================================
function Vlate = get_late_spacetime(ic0, sigma, b, T_relax, T_obs, n_keep)
par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

if isempty(V)
    Vlate = zeros(n_keep, size(ic0,1));
    return;
end

n_keep = min(n_keep, size(V,1));
Vlate = V(end-n_keep+1:end,:);
end

function Z = normalize_tile(Vlate)
if max(Vlate(:)) - min(Vlate(:)) > 1e-8
    Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    Z = zeros(size(Vlate));
end
end
