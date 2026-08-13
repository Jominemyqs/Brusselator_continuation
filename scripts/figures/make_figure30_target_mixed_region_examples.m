function make_figure30_target_mixed_region_examples
% MAKE_FIGURE30_TARGET_MIXED_REGION_EXAMPLES
%
% Figure 30:
% Representative ambiguous / mixed localized-source patterns.
%
% This version includes only the first two examples:
%   (a) spiral-like / source-defect-like
%   (b) spiral-like / source-defect-like
%
% Intended for Section 6.6 (Failure modes and ambiguous regions).
%
% Requires:
%   - ic_0.45_10.mat
%   - solve_brusselator_1d.m

clear; clc; close all;

save_png = true;
png_name = figure_output_path('Figure30_ambiguous_cross_family_examples.png');

%% ---------------- USER SETTINGS ----------------
n_keep = 120;

% Case definitions
cases = struct([]);

cases(1).label     = '(a)';
cases(1).sigma     = 0.420;
cases(1).b         = 9.400;
cases(1).seed_type = 'spiral';
cases(1).family    = 'spiral-like';

cases(2).label     = '(b)';
cases(2).sigma     = 0.330;
cases(2).b         = 14.200;
cases(2).seed_type = 'spiral';
cases(2).family    = 'spiral-like';

show_metric_text = false;  % keep titles clean

T_relax_spiral = 360;
T_obs_spiral   = 360;

%% ---------------- LOAD SEED ----------------
S1 = load('ic_0.45_10.mat');
if isfield(S1,'ic')
    ic_spiral = S1.ic;
else
    fn = fieldnames(S1);
    ic_spiral = S1.(fn{1});
end

%% ---------------- BUILD FIGURE ----------------
fig = figure('Color','w','Position',[90 90 1150 520]);
tl = tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');

for k = 1:numel(cases)

    ic0 = ic_spiral;
    T_relax = T_relax_spiral;
    T_obs   = T_obs_spiral;

    [Vlate, metric_text] = get_spacetime_and_optional_text( ...
        ic0, cases(k).sigma, cases(k).b, T_relax, T_obs, n_keep, cases(k).family);

    ax = nexttile(tl);
    imagesc(ax, Vlate);
    set(ax,'YDir','normal');
    axis(ax,'tight');
    colormap(ax, parula(256));

    cb = colorbar(ax);
    if k == numel(cases)
        cb.Label.String = 'activity';
        cb.Label.FontSize = 14;
    end

    xlabel(ax,'space index','FontSize',15);
    if k == 1
        ylabel(ax,'late-time index','FontSize',15);
    end

    if show_metric_text
        title(ax, sprintf('%s %s\n$(\\sigma,b)=(%.3f,%.3f)$\n%s', ...
            cases(k).label, cases(k).family, cases(k).sigma, cases(k).b, metric_text), ...
            'Interpreter','latex', 'FontSize',15);
    else
        title(ax, sprintf('%s %s\n$(\\sigma,b)=(%.3f,%.3f)$', ...
            cases(k).label, cases(k).family, cases(k).sigma, cases(k).b), ...
            'Interpreter','latex', 'FontSize',16);
    end

    set(ax,'FontSize',13);
end



if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved Figure 30 to %s\n', png_name);
end

end

%% ========================================================================
function [Vlate, metric_text] = get_spacetime_and_optional_text(ic0, sigma, b, T_relax, T_obs, n_keep, family)
par = struct('sigma', sigma, 'b', b);

[ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
[~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

if isempty(V)
    Vlate = zeros(n_keep, size(ic0,1));
else
    n_keep = min(n_keep, size(V,1));
    Vlate = V(end-n_keep+1:end,:);
end

switch lower(family)
    case 'spiral-like'
        Z = Vlate;
        if max(Z(:)) - min(Z(:)) > 1e-8
            Z = reshape(zscore(Z(:)), size(Z,1), size(Z,2));
        else
            Z = zeros(size(Z));
        end
        rVspace = norm(Z - fliplr(Z), 'fro') / max(norm(Z,'fro'), eps);
        rVtime  = norm(Z - flipud(Z), 'fro') / max(norm(Z,'fro'), eps);
        metric_text = sprintf('$r_{V,lr}=%.2e$, $r_{V,ud}=%.2e$', rVspace, rVtime);

    otherwise
        metric_text = '';
end
end
