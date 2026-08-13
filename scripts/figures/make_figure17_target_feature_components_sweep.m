function make_figure17_target_feature_components_sweep
% MAKE_FIGURE17_TARGET_FEATURE_COMPONENTS_SWEEP
%
% Figure 17:
% Horizontal sweep showing S_c, T_r, and min{S_c,T_r}.
%
% Intended for Section 5.5.

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure17_target_feature_components_sweep.png');

%% ---------------- USER SETTINGS ----------------
b_fixed = 10.02;

sigma_start = 0.315;
sigma_end   = 0.332;
ds_sigma    = 0.0005;

T_relax = 240;
T_obs   = 240;

Nx_est = 256;
core_x1 = 1;
core_x2 = 70;
tail_x1 = 180;
tail_x2 = 256;

%% ---------------- LOAD SEED ----------------
S = load('half_target_seed.mat');
if isfield(S,'ic')
    ic0 = S.ic;
elseif isfield(S,'ic_end')
    ic0 = S.ic_end;
else
    fn = fieldnames(S);
    ic0 = S.(fn{1});
end

%% ---------------- SWEEP ----------------
sigmas = sigma_start:ds_sigma:sigma_end;
M = numel(sigmas);

Sc_vals = zeros(1,M);
Tr_vals = zeros(1,M);
Ft_vals = zeros(1,M);

for j = 1:M
    Vlate = get_late_spacetime(ic0, sigmas(j), b_fixed, T_relax, T_obs);

    Nx = size(Vlate,2);
    scale = Nx / Nx_est;

    cx1 = max(1, round(core_x1 * scale));
    cx2 = min(Nx, round(core_x2 * scale));
    tx1 = max(1, round(tail_x1 * scale));
    tx2 = min(Nx, round(tail_x2 * scale));

    Vcore = Vlate(:,cx1:cx2);
    Vtail = Vlate(:,tx1:tx2);

    Sc = mean(var(Vcore, 0, 2));  % spatial variance over core, averaged over time
    Tr = mean(var(Vtail, 0, 1));  % temporal variance over tail, averaged over space
    Ft = min(Sc, Tr);

    Sc_vals(j) = Sc;
    Tr_vals(j) = Tr;
    Ft_vals(j) = Ft;
end

%% ---------------- PLOT ----------------
fig = figure('Color','w','Position',[100 100 1000 700]);
ax = axes(fig); hold(ax,'on');

plot(ax, sigmas, Sc_vals, 'o-', 'LineWidth',1.8, 'MarkerSize',5, ...
    'DisplayName','$S_c$');
plot(ax, sigmas, Tr_vals, 's-', 'LineWidth',1.8, 'MarkerSize',5, ...
    'DisplayName','$T_r$');


hold(ax,'off');
grid(ax,'on');
xlabel(ax,'\sigma', 'Interpreter','tex', 'FontSize',16);
ylabel(ax,'feature component value', 'FontSize',16);

set(ax,'FontSize',14);

if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 17 to %s\n', png_name);
end

end

%% ========================================================================
function V = get_late_spacetime(ic0, sigma, b, T_relax, T_obs)
par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);
end
