function make_figure10_spiral_symmetry_defect_construction
% MAKE_FIGURE10_SPIRAL_SYMMETRY_DEFECT_CONSTRUCTION
%
% Figure 10:
% Visual construction of the branch-aware spiral symmetry-defect feature.
%
% Top row: representative right-branch example using left-right reflection
% Bottom row: representative left-branch example using up-down reflection
%
% Columns:
%   (1) original V_late
%   (2) reflected V_late
%   (3) absolute difference
%
% Requires:
%   - ic_0.45_10.mat
%   - solve_brusselator_1d.m

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure10_spiral_symmetry_defect_construction.png');

%% ---------------- USER SETTINGS ----------------
n_keep  = 120;     % number of late-time frames to display
T_relax = 360;
T_obs   = 1200;    % use longer observation for cleaner late-time behavior

% Representative examples
right_ex.sigma = 0.505;
right_ex.b     = 10.000;

left_ex.sigma  = 0.330;
left_ex.b      = 10.000;

%% ---------------- LOAD SPIRAL SEED ----------------
S = load('ic_0.45_10.mat');
if isfield(S,'ic')
    ic0 = S.ic;
else
    fn = fieldnames(S);
    ic0 = S.(fn{1});
end

%% ---------------- SIMULATE ----------------
V_right = get_late_spacetime(ic0, right_ex.sigma, right_ex.b, T_relax, T_obs, n_keep);
V_left  = get_late_spacetime(ic0, left_ex.sigma,  left_ex.b,  T_relax, T_obs, n_keep);

% branch-aware reflections on RAW data (for actual metric)
V_right_ref = fliplr(V_right);
V_left_ref  = flipud(V_left);

D_right = abs(V_right - V_right_ref);
D_left  = abs(V_left  - V_left_ref);

% numerical defect values from raw data
rV_lr = norm(V_right - V_right_ref, 'fro') / max(norm(V_right,'fro'), eps);
rV_ud = norm(V_left  - V_left_ref,  'fro') / max(norm(V_left,'fro'),  eps);

% normalized copies for plotting only
Zr     = normalize_tile(V_right);
Zr_ref = normalize_tile(V_right_ref);

Zl     = normalize_tile(V_left);
Zl_ref = normalize_tile(V_left_ref);

%% ---------------- BUILD FIGURE ----------------
fig = figure('Color','w','Position',[60 60 1700 950]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% ---- Top row: right branch ----
ax1 = nexttile(tl);
imagesc(ax1, Zr);
set(ax1,'YDir','normal');
axis(ax1,'tight');
colorbar(ax1);
xlabel(ax1,'space index');
ylabel(ax1,'late-time index');
title(ax1, sprintf('(a) Right branch: $V_{\\mathrm{late}}$\n$(\\sigma,b)=(%.3f,%.3f)$', ...
    right_ex.sigma, right_ex.b), 'Interpreter','latex', 'FontSize',16);

ax2 = nexttile(tl);
imagesc(ax2, Zr_ref);
set(ax2,'YDir','normal');
axis(ax2,'tight');
colorbar(ax2);
xlabel(ax2,'space index');
ylabel(ax2,'late-time index');
title(ax2, '(b) $\mathrm{fliplr}(V_{\mathrm{late}})$', ...
    'Interpreter','latex', 'FontSize',16);

ax3 = nexttile(tl);
imagesc(ax3, D_right);
set(ax3,'YDir','normal');
axis(ax3,'tight');
colorbar(ax3);
xlabel(ax3,'space index');
ylabel(ax3,'late-time index');
title(ax3, sprintf('(c) $|V_{\\mathrm{late}}-\\mathrm{fliplr}(V_{\\mathrm{late}})|$\n$r_{V,lr}=%.3e$', ...
    rV_lr), 'Interpreter','latex', 'FontSize',16);

% ---- Bottom row: left branch ----
ax4 = nexttile(tl);
imagesc(ax4, Zl);
set(ax4,'YDir','normal');
axis(ax4,'tight');
colorbar(ax4);
xlabel(ax4,'space index');
ylabel(ax4,'late-time index');
title(ax4, sprintf('(d) Left branch: $V_{\\mathrm{late}}$\n$(\\sigma,b)=(%.3f,%.3f)$', ...
    left_ex.sigma, left_ex.b), 'Interpreter','latex', 'FontSize',16);

ax5 = nexttile(tl);
imagesc(ax5, Zl_ref);
set(ax5,'YDir','normal');
axis(ax5,'tight');
colorbar(ax5);
xlabel(ax5,'space index');
ylabel(ax5,'late-time index');
title(ax5, '(e) $\mathrm{flipud}(V_{\mathrm{late}})$', ...
    'Interpreter','latex', 'FontSize',16);

ax6 = nexttile(tl);
imagesc(ax6, D_left);
set(ax6,'YDir','normal');
axis(ax6,'tight');
colorbar(ax6);
xlabel(ax6,'space index');
ylabel(ax6,'late-time index');
title(ax6, sprintf('(f) $|V_{\\mathrm{late}}-\\mathrm{flipud}(V_{\\mathrm{late}})|$\n$r_{V,ud}=%.3e$', ...
    rV_ud), 'Interpreter','latex', 'FontSize',16);

% ---- Colormaps / style ----
for ax = [ax1 ax2 ax4 ax5]
    colormap(ax, parula(256));
    set(ax,'FontSize',12);
end
for ax = [ax3 ax6]
    colormap(ax, hot(256));
    set(ax,'FontSize',12);
end

if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 10 to %s\n', png_name);
end

fprintf('\nComputed defect values:\n');
fprintf('  Right branch (lr): rV_lr = %.6e at (sigma,b) = (%.3f, %.3f)\n', ...
    rV_lr, right_ex.sigma, right_ex.b);
fprintf('  Left  branch (ud): rV_ud = %.6e at (sigma,b) = (%.3f, %.3f)\n', ...
    rV_ud, left_ex.sigma, left_ex.b);

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
Vlate = V(end-n_keep+1:end, :);

end

function Z = normalize_tile(Vlate)
% z-score for display only
if max(Vlate(:)) - min(Vlate(:)) > 1e-12
    Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    Z = zeros(size(Vlate));
end
end
