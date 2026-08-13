function make_figure23_spiral_LCR_diagnostics
% MAKE_FIGURE23_SPIRAL_LCR_DIAGNOSTICS
%
% Figure 23:
% L/C/R spacetime diagnostics for a representative spiral continuation point.
%
% Default: right-branch example using the left-right spacetime symmetry metric.
%
% Requires:
%   - ic_0.45_10.mat
%   - solve_brusselator_1d.m

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure23_spiral_LCR_diagnostics.png');

%% ---------------- USER SETTINGS ----------------
% Representative branch point
b_fixed = 10.00;
sigmaC  = 0.5100;
eps_trip = 0.0010;

sigmaL = sigmaC - eps_trip;
sigmaR = sigmaC + eps_trip;

T_relax = 360;
T_obs   = 360;
n_keep  = 120;

branch_label = 'Right spiral branch';

%% ---------------- LOAD SPIRAL SEED ----------------
S = load('ic_0.45_10.mat');
if isfield(S,'ic')
    ic0 = S.ic;
else
    fn = fieldnames(S);
    ic0 = S.(fn{1});
end

%% ---------------- SIMULATE ----------------
dataL = sim_obs(ic0, sigmaL, b_fixed, T_relax, T_obs, n_keep);
dataC = sim_obs(ic0, sigmaC, b_fixed, T_relax, T_obs, n_keep);
dataR = sim_obs(ic0, sigmaR, b_fixed, T_relax, T_obs, n_keep);

% branch-aware metric for right branch
rfun = @(V) norm(V - fliplr(V), 'fro') / max(norm(V,'fro'), eps);

rL = rfun(dataL.Z);
rC = rfun(dataC.Z);
rR = rfun(dataR.Z);

%% ---------------- FIGURE ----------------
fig = figure('Color','w','Position',[80 90 1600 560]);
tl = tiledlayout(fig,1,3,'TileSpacing','compact','Padding','compact');

ax1 = nexttile(tl);
imagesc(ax1, dataL.Z); set(ax1,'YDir','normal'); axis(ax1,'tight');
colormap(ax1, parula(256)); colorbar(ax1);
xlabel(ax1,'space index','FontSize',15);
ylabel(ax1,'late-time index','FontSize',15);
title(ax1, sprintf('(a) L: $\\sigma=%.4f$, $r_{V,lr}=%.3e$', sigmaL, rL), ...
    'Interpreter','latex', 'FontSize',16);
set(ax1,'FontSize',13);

ax2 = nexttile(tl);
imagesc(ax2, dataC.Z); set(ax2,'YDir','normal'); axis(ax2,'tight');
colormap(ax2, parula(256)); colorbar(ax2);
xlabel(ax2,'space index','FontSize',15);
ylabel(ax2,'late-time index','FontSize',15);
title(ax2, sprintf('(b) C: $\\sigma=%.4f$, $r_{V,lr}=%.3e$', sigmaC, rC), ...
    'Interpreter','latex', 'FontSize',16);
set(ax2,'FontSize',13);

ax3 = nexttile(tl);
imagesc(ax3, dataR.Z); set(ax3,'YDir','normal'); axis(ax3,'tight');
colormap(ax3, parula(256)); colorbar(ax3);
xlabel(ax3,'space index','FontSize',15);
ylabel(ax3,'late-time index','FontSize',15);
title(ax3, sprintf('(c) R: $\\sigma=%.4f$, $r_{V,lr}=%.3e$', sigmaR, rR), ...
    'Interpreter','latex', 'FontSize',16);
set(ax3,'FontSize',13);

title(tl, sprintf('%s at fixed $b=%.2f$', branch_label, b_fixed), ...
    'Interpreter','latex', 'FontSize',18);

if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 23 to %s\n', png_name);
end

end

%% ========================================================================
function data = sim_obs(ic0, sigma, b, T_relax, T_obs, n_keep)
par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

if isempty(V)
    Vlate = zeros(n_keep, size(ic0,1));
else
    n_keep = min(n_keep, size(V,1));
    Vlate = V(end-n_keep+1:end,:);
end

data.V = Vlate;
if max(Vlate(:)) - min(Vlate(:)) > 1e-8
    data.Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    data.Z = zeros(size(Vlate));
end
end
