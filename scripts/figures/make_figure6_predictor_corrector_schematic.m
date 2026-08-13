function make_figure6_predictor_corrector_schematic
% MAKE_FIGURE6_PREDICTOR_CORRECTOR_SCHEMATIC
%
% Figure 6:
% Schematic of the secant predictor + local horizontal sweep corrector
% in the (sigma,b)-plane.

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure6_predictor_corrector_schematic.png');

%% ---------------- SETUP ----------------
% Two previously accepted points
pkm1 = [0.332, 10.20];
pk   = [0.338, 10.55];

% predictor step length multiplier
alpha = 0.85;

% secant direction
d = pk - pkm1;
d = d / norm(d);

% predicted point
ppred = pk + alpha * d;

% local horizontal sweep at fixed b_pred
half_width = 0.018;
sigma_left  = ppred(1) - half_width;
sigma_right = ppred(1) + half_width;

%% ---------------- NOTIONAL TRANSITION CURVE ----------------
sig_curve = linspace(0.318, 0.368, 1000);
b_curve = 9.85 + 18*(sig_curve-0.318) - 160*(sig_curve-0.343).^2;

% corrected point = intersection of horizontal line b = b_pred with curve
[~, idx_corr] = min(abs(b_curve - ppred(2)));
pcorr = [sig_curve(idx_corr), b_curve(idx_corr)];

%% ---------------- DRAW FIGURE ----------------
fig = figure('Color','w','Position',[100 100 950 700]);
ax = axes(fig); hold(ax,'on');

% --- transition curve ---
plot(ax, sig_curve, b_curve, 'k-', 'LineWidth', 2.5);

% --- secant segment between previous points ---
plot(ax, [pkm1(1), pk(1)], [pkm1(2), pk(2)], '-', ...
    'Color', [0.2 0.4 0.8], 'LineWidth', 2.0);

% --- predictor arrow ---
quiver(ax, pk(1), pk(2), ppred(1)-pk(1), ppred(2)-pk(2), 0, ...
    'Color', [0.2 0.4 0.8], 'LineWidth', 2.2, 'MaxHeadSize', 0.8);

% --- horizontal sweep line ---
plot(ax, [sigma_left, sigma_right], [ppred(2), ppred(2)], 'r--', 'LineWidth', 2.0);

% --- dashed guide from predictor to corrected point ---
plot(ax, [ppred(1), pcorr(1)], [ppred(2), pcorr(2)], ':', ...
    'Color', [0.6 0 0], 'LineWidth', 1.8);

% --- points ---
plot(ax, pkm1(1), pkm1(2), 'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k');
plot(ax, pk(1),   pk(2),   'ko', 'MarkerSize', 9, 'MarkerFaceColor', 'k');
plot(ax, ppred(1),ppred(2),'o', ...
    'Color',[0.2 0.4 0.8], 'MarkerSize', 10, 'MarkerFaceColor', [0.2 0.4 0.8]);
plot(ax, pcorr(1),pcorr(2),'rs', ...
    'MarkerSize', 10, 'MarkerFaceColor', 'r');

% --- annotations ---
text(pkm1(1)-0.004, pkm1(2)-0.08, '$p_{k-1}$', ...
    'Interpreter','latex', 'FontSize',16);
text(pk(1)+0.001,   pk(2)-0.08,   '$p_k$', ...
    'Interpreter','latex', 'FontSize',16);
text(ppred(1)+0.001, ppred(2)+0.03, '$p_{\mathrm{pred}}$', ...
    'Interpreter','latex', 'FontSize',16, 'Color',[0.2 0.4 0.8]);
text(pcorr(1)+0.001, pcorr(2)-0.08, '$p_{k+1}$', ...
    'Interpreter','latex', 'FontSize',16, 'Color',[0.6 0 0]);

text((sigma_left+sigma_right)/2 - 0.010, ppred(2)+0.06, ...
    'local sweep at fixed $b_{\mathrm{pred}}$', ...
    'Interpreter','latex', 'FontSize',15, 'Color',[0.65 0 0]);

%% ---------------- STYLE ----------------
xlabel(ax, '\sigma', 'Interpreter','tex', 'FontSize',17);
ylabel(ax, 'b',      'Interpreter','tex', 'FontSize',17);

set(ax, 'FontSize', 14);
box(ax, 'on');
grid(ax, 'on');

xlim(ax, [0.318 0.368]);
ylim(ax, [9.95 11.25]);

axis(ax, 'square');

%% ---------------- SAVE ----------------
if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 6 to %s\n', png_name);
end
end
