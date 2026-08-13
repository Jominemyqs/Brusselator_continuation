function make_figure24_spiral_sweep_validation
% MAKE_FIGURE24_SPIRAL_SWEEP_VALIDATION
%
% Figure 24:
% Fixed-b horizontal sweep validation for the spiral feature.
%
% Default: right-branch example using r_{V,lr}.
%
% Layout:
%   top: metric vs sigma
%   bottom: L/C/R spacetime plots
%
% Requires:
%   - ic_0.45_10.mat
%   - solve_brusselator_1d.m

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure24_spiral_sweep_validation.png');

%% ---------------- USER SETTINGS ----------------
b_fixed = 10.00;

sigma_start = 0.500;
sigma_end   = 0.520;
ds_sigma    = 0.001;

r_thr = 2.6e-3;

use_manual_triplet = true;
sigmaL = 0.5090;
sigmaC = 0.5100;
sigmaR = 0.5110;

T_relax = 360;
T_obs   = 360;
n_keep  = 120;

%% ---------------- LOAD SEED ----------------
S = load('ic_0.45_10.mat');
if isfield(S,'ic')
    ic0 = S.ic;
else
    fn = fieldnames(S);
    ic0 = S.(fn{1});
end

%% ---------------- SWEEP ----------------
sigmas = sigma_start:ds_sigma:sigma_end;
M = numel(sigmas);
rvals = zeros(1,M);

for j = 1:M
    Z = get_metric_tile(ic0, sigmas(j), b_fixed, T_relax, T_obs, n_keep);
    rvals(j) = norm(Z - fliplr(Z), 'fro') / max(norm(Z,'fro'), eps);
end

if ~use_manual_triplet
    sigmaC = first_downward_crossing(sigmas, rvals, r_thr);
    sigmaL = sigmaC - ds_sigma;
    sigmaR = sigmaC + ds_sigma;
end

dataL = get_metric_tile(ic0, sigmaL, b_fixed, T_relax, T_obs, n_keep);
dataC = get_metric_tile(ic0, sigmaC, b_fixed, T_relax, T_obs, n_keep);
dataR = get_metric_tile(ic0, sigmaR, b_fixed, T_relax, T_obs, n_keep);

rL = norm(dataL - fliplr(dataL), 'fro') / max(norm(dataL,'fro'), eps);
rC = norm(dataC - fliplr(dataC), 'fro') / max(norm(dataC,'fro'), eps);
rR = norm(dataR - fliplr(dataR), 'fro') / max(norm(dataR,'fro'), eps);

%% ---------------- FIGURE ----------------
fig = figure('Color','w','Position',[70 60 1300 950]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% top row spans all columns
ax1 = nexttile(tl, [1 3]);
plot(ax1, sigmas, rvals, 'o-', 'LineWidth',1.6, 'MarkerSize',5, ...
    'Color',[0 0.35 0.8], 'MarkerFaceColor','w');
hold(ax1,'on');
yline(ax1, r_thr, '--', 'LineWidth',1.3, 'Label','$r_{\mathrm{thr}}$', ...
    'Interpreter','latex', ...
    'LabelHorizontalAlignment','right', ...
    'LabelVerticalAlignment','bottom');

xline(ax1, sigmaL, ':', 'L', 'LineWidth',1.2);
xline(ax1, sigmaC, '--', 'C', 'LineWidth',1.2);
xline(ax1, sigmaR, ':', 'R', 'LineWidth',1.2);

plot(ax1, sigmaL, rL, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);
plot(ax1, sigmaC, rC, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);
plot(ax1, sigmaR, rR, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);

hold(ax1,'off');
grid(ax1,'on');
xlabel(ax1,'\sigma', 'Interpreter','tex', 'FontSize',15);
ylabel(ax1,'$r_{V,lr}$', 'Interpreter','latex', 'FontSize',15);
title(ax1, sprintf('(a) Spiral metric along a horizontal sweep at b = %.2f', b_fixed), ...
    'FontSize',17);
set(ax1,'FontSize',13);

% bottom row: L/C/R
ax2 = nexttile(tl);
imagesc(ax2, dataL); set(ax2,'YDir','normal'); axis(ax2,'tight');
colormap(ax2, parula(256)); colorbar(ax2);
xlabel(ax2,'space index','FontSize',14);
ylabel(ax2,'late-time index','FontSize',14);
title(ax2, sprintf('(b) L: $\\sigma=%.4f$', sigmaL), 'Interpreter','latex', 'FontSize',15);

ax3 = nexttile(tl);
imagesc(ax3, dataC); set(ax3,'YDir','normal'); axis(ax3,'tight');
colormap(ax3, parula(256)); colorbar(ax3);
xlabel(ax3,'space index','FontSize',14);
ylabel(ax3,'late-time index','FontSize',14);
title(ax3, sprintf('(c) C: $\\sigma=%.4f$', sigmaC), 'Interpreter','latex', 'FontSize',15);

ax4 = nexttile(tl);
imagesc(ax4, dataR); set(ax4,'YDir','normal'); axis(ax4,'tight');
colormap(ax4, parula(256)); colorbar(ax4);
xlabel(ax4,'space index','FontSize',14);
ylabel(ax4,'late-time index','FontSize',14);
title(ax4, sprintf('(d) R: $\\sigma=%.4f$', sigmaR), 'Interpreter','latex', 'FontSize',15);

if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 24 to %s\n', png_name);
end

end

%% ========================================================================
function Z = get_metric_tile(ic0, sigma, b, T_relax, T_obs, n_keep)
par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

if isempty(V)
    Vlate = zeros(n_keep, size(ic0,1));
else
    n_keep = min(n_keep, size(V,1));
    Vlate = V(end-n_keep+1:end,:);
end

if max(Vlate(:)) - min(Vlate(:)) > 1e-8
    Z = reshape(zscore(Vlate(:)), size(Vlate,1), size(Vlate,2));
else
    Z = zeros(size(Vlate));
end
end

function sigC = first_downward_crossing(sigmas, y, thr)
above = (y >= thr);
idx = find(above(1:end-1) & ~above(2:end), 1, 'first');
if isempty(idx)
    sigC = NaN;
else
    sigC = 0.5*(sigmas(idx) + sigmas(idx+1));
end
end
