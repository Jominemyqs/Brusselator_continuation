function make_figure5_target_transition_motivation
% MAKE_FIGURE5_TARGET_TRANSITION_MOTIVATION
%
% Figure 5:
% Fixed-b horizontal sweep showing that a sharp change in the target
% feature corresponds to a sharp qualitative change in spacetime pattern.
%
% Layout:
%   Top row: target feature vs sigma with threshold and L/C/R markers
%   Bottom row: L/C/R late-time spacetime plots
%
% Intended for Chapter 2.3: Why Transition Curves Matter
%
% Requires:
%   - half_target_seed.mat
%   - solve_brusselator_1d.m
%   - feature_evaluation_target.m

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
b_fixed = 10.02;

sigma_start = 0.305;
sigma_end   = 0.345;
ds_sigma    = 0.001;

T_relax = 240;
T_obs   = 240;

% threshold for target feature
r_thr = 8.0;

% use these specific L/C/R sigmas if you want exact control
use_manual_triplet = true;
sigmaL = 0.3245;
sigmaC = 0.3255;
sigmaR = 0.3265;

% if false, choose automatically from sweep
eps_trip = 1e-3;

save_png = true;
png_name = figure_output_path('Figure5_target_transition_motivation.png');

%% ---------------- LOAD TARGET SEED ----------------
Sseed = load('half_target_seed.mat');
if isfield(Sseed,'ic')
    ic_seed = Sseed.ic;
elseif isfield(Sseed,'ic_end')
    ic_seed = Sseed.ic_end;
else
    fn = fieldnames(Sseed);
    ic_seed = Sseed.(fn{1});
end

%% ---------------- SWEEP IN SIGMA ----------------
sig_grid = sigma_start:ds_sigma:sigma_end;
M = numel(sig_grid);
feat_vals = zeros(1,M);

fprintf('Running Figure 5 sweep at b = %.6f\n', b_fixed);

for j = 1:M
    feat_vals(j) = eval_target_feature(ic_seed, sig_grid(j), b_fixed, T_relax, T_obs);
end

%% ---------------- CHOOSE L/C/R ----------------
if use_manual_triplet
    % keep your known clean example
    % nothing to do
else
    [sigmaC, ~] = choose_sigma_up(sig_grid, feat_vals, r_thr);
    sigmaL = max(sig_grid(1), sigmaC - eps_trip);
    sigmaR = min(sig_grid(end), sigmaC + eps_trip);
end

featL = eval_target_feature(ic_seed, sigmaL, b_fixed, T_relax, T_obs);
featC = eval_target_feature(ic_seed, sigmaC, b_fixed, T_relax, T_obs);
featR = eval_target_feature(ic_seed, sigmaR, b_fixed, T_relax, T_obs);

%% ---------------- SIMULATE L/C/R ----------------
dataL = sim_obs(ic_seed, sigmaL, b_fixed, T_relax, T_obs);
dataC = sim_obs(ic_seed, sigmaC, b_fixed, T_relax, T_obs);
dataR = sim_obs(ic_seed, sigmaR, b_fixed, T_relax, T_obs);

%% ---------------- MAKE FIGURE ----------------
fig = figure('Color','w','Position',[80 80 1200 900]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing','compact', 'Padding','compact');

% ---- top: feature plot spanning all columns ----
ax1 = nexttile(tl, [1 3]);
plot(ax1, sig_grid, feat_vals, 'o-', 'LineWidth',1.5, 'MarkerSize',6);
hold(ax1, 'on');
yline(ax1, r_thr, '--', 'LineWidth',1.2, 'Label','r_{thr}', ...
    'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom');

xline(ax1, sigmaL, ':', 'L', 'LineWidth',1.2, ...
    'LabelVerticalAlignment','top', 'LabelHorizontalAlignment','center');
xline(ax1, sigmaC, '--', 'C', 'LineWidth',1.2, ...
    'LabelVerticalAlignment','top', 'LabelHorizontalAlignment','center');
xline(ax1, sigmaR, ':', 'R', 'LineWidth',1.2, ...
    'LabelVerticalAlignment','top', 'LabelHorizontalAlignment','center');

plot(ax1, sigmaL, featL, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);
plot(ax1, sigmaC, featC, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);
plot(ax1, sigmaR, featR, 'ko', 'MarkerFaceColor','w', 'MarkerSize',8);

hold(ax1, 'off');
grid(ax1, 'on');
xlabel(ax1, '\sigma', 'Interpreter','tex', 'FontSize',15);
ylabel(ax1, 'target feature', 'Interpreter','tex', 'FontSize',15);

set(ax1, 'FontSize',13);

% ---- bottom left: L ----
ax2 = nexttile(tl);
imagesc(ax2, dataL.V);
set(ax2, 'YDir','normal');
axis(ax2, 'tight');
colormap(ax2, parula);
colorbar(ax2);
xlabel(ax2, 'space index', 'FontSize',14);
ylabel(ax2, 'late-time index', 'FontSize',14);
set(ax2, 'FontSize',12);

% ---- bottom center: C ----
ax3 = nexttile(tl);
imagesc(ax3, dataC.V);
set(ax3, 'YDir','normal');
axis(ax3, 'tight');
colormap(ax3, parula);
colorbar(ax3);
xlabel(ax3, 'space index', 'FontSize',14);
ylabel(ax3, 'late-time index', 'FontSize',14);
set(ax3, 'FontSize',12);

% ---- bottom right: R ----
ax4 = nexttile(tl);
imagesc(ax4, dataR.V);
set(ax4, 'YDir','normal');
axis(ax4, 'tight');
colormap(ax4, parula);
colorbar(ax4);
xlabel(ax4, 'space index', 'FontSize',14);
ylabel(ax4, 'late-time index', 'FontSize',14);
set(ax4, 'FontSize',12);

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 5 to %s\n', png_name);
end

fprintf('\nFigure 5 summary:\n');
fprintf('  b      = %.6f\n', b_fixed);
fprintf('  sigmaL = %.6f, featL = %.6e\n', sigmaL, featL);
fprintf('  sigmaC = %.6f, featC = %.6e\n', sigmaC, featC);
fprintf('  sigmaR = %.6f, featR = %.6e\n', sigmaR, featR);

end

%% ========================================================================
function feat = eval_target_feature(ic0, sigma, b, T_relax, T_obs)
featpar = struct();
featpar.feature = 'target';
featpar.N = 1;

modelpar = struct();
modelpar.model   = 'Brusselator_1D';
modelpar.ic0     = ic0;
modelpar.a       = sigma;
modelpar.b       = b;
modelpar.T_relax = T_relax;
modelpar.T_obs   = T_obs;

feat = feature_evaluation_target(featpar, modelpar);
end

function data = sim_obs(ic0, sigma, b, T_relax, T_obs)
par = struct('sigma', sigma, 'b', b);

[ic_relaxed, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
[ic_end, ~, V]     = solve_brusselator_1d(ic_relaxed, par, T_obs, 0);

data.ic_end  = ic_end;
data.u_final = ic_end(:,2);
data.V       = V;
end

function [sigC, idx] = choose_sigma_up(sig_grid, y, r_thr)
above = (y >= r_thr);
cross_all = find(~above(1:end-1) & above(2:end));

if isempty(cross_all)
    [~, idx] = min(abs(y - r_thr));
    sigC = sig_grid(idx);
else
    idx = cross_all(1);
    sigC = 0.5*(sig_grid(idx) + sig_grid(idx+1));
end
end
