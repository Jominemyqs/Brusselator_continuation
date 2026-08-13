%% vertical_sweep_compare_metrics_target.m
% Vertical sweep in b for multiple fixed sigma values.
%
% For each fixed sigma, sweep over b and compare candidate metrics.
% Intended for diagnosing lower / middle target branches with U/C/D logic.
%
% Computes:
%   (1) target_feature   = feature_evaluation_target(featpar, modelpar)
%   (2) r_coh           = temporal coherence score (smaller = more regular)
%   (3) optionally save spacetime plots and final profiles
%
% Output per sigma:
%   - sweep_log.csv
%   - metrics_overlay.png
%   - target_feature_vs_b.png
%   - coherence_vs_b.png
%   - spacetime/profile pngs per b (if enabled)

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
sigma_list = [0.42, 0.44, 0.46, 0.47, 0.48];

b_start = 9.4;
b_end   = 9.6;
db      = 0.02;    % try 0.02 first, then 0.01 if needed

T_relax = 240;

% diagnostic observation after relaxation
do_spacetime = true;
do_profiles  = true;
T_diag       = 240;

% use only last fraction of V for spacetime-based metrics
use_last_window = true;
last_frac       = 0.25;

% IMPORTANT:
% false = independent b evaluations from same IC
% true  = warm-start along b
use_b_warmstart = false;

% target seed
ic_file = 'half_target_seed.mat';

% threshold line for reference
r_thr = 10;

timestamp = datestr(now,'yyyymmdd_HHMMSS');
rootdir = sprintf('vertSweep_target_compareMetrics_db%.0e_%s', db, timestamp);
if ~exist(rootdir,'dir'); mkdir(rootdir); end

%% ---------------- LOAD SEED IC ----------------
S = load(ic_file);
if isfield(S,'ic')
    ic_seed = S.ic;
elseif isfield(S,'ic_end')
    ic_seed = S.ic_end;
else
    fn = fieldnames(S);
    ic_seed = S.(fn{1});
end

%% ---------------- b GRID ----------------
bvals = b_start:db:b_end;
nb = numel(bvals);

fprintf('Running sweeps: %d sigma-values, %d b-values each\n', numel(sigma_list), nb);
fprintf('use_b_warmstart = %d\n', use_b_warmstart);

for ss = 1:numel(sigma_list)
    sigma_fixed = sigma_list(ss);

    outdir = fullfile(rootdir, sprintf('sigma_%.4f', sigma_fixed));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    % optional warm-start state along b
    ic0 = ic_seed;

    % columns:
    % b, target_feature, r_coh
    log_mat = nan(nb, 3);

    fprintf('\n=== sigma = %.6f ===\n', sigma_fixed);

    for i = 1:nb
        b = bvals(i);

        fprintf('[sigma=%.4f] %d/%d b=%.6f\n', sigma_fixed, i, nb, b);

        stepdir = fullfile(outdir, sprintf('b_%.6f', b));
        if (do_spacetime || do_profiles) && ~exist(stepdir,'dir')
            mkdir(stepdir);
        end

        par = struct('sigma', sigma_fixed, 'b', b);

        % ---- Relaxation ----
        ic_rel_start = ic_seed;
        if use_b_warmstart
            ic_rel_start = ic0;
        end

        [ic_rel, ~, ~] = solve_brusselator_1d(ic_rel_start, par, T_relax, 0);

        % ---- Diagnostic observation from relaxed state ----
        [ic_end, ~, V_diag] = solve_brusselator_1d(ic_rel, par, T_diag, 0);

        % ---- Target feature ----
        featpar = struct();
        featpar.feature = 'target';
        featpar.alpha = 0.01;
        featpar.avoid_steady = 0;
        featpar.N = 1;

        modelpar = struct();
        modelpar.model = 'Brusselator_1D';
        modelpar.a = sigma_fixed;
        modelpar.b = b;
        modelpar.ic0 = ic_rel_start;
        modelpar.T_relax = T_relax;
        modelpar.T_obs = T_diag;
        modelpar.ic_base = ic_seed;
        modelpar.ic_det = ic_seed;
        modelpar.crossing_direction = 'down';  % only for gate fallback if used internally

        target_feature = feature_evaluation_target(featpar, modelpar);

        % ---- coherence score ----
        r_coh = temporal_coherence_score(V_diag, use_last_window, last_frac);

        log_mat(i,:) = [b, target_feature, r_coh];

        % ---- Save spacetime ----
        if do_spacetime
            save_spacetime(stepdir, sprintf('target_sigma%.6f_b%.6f', sigma_fixed, b), ...
                V_diag, target_feature, r_coh);
        end

        % ---- Save final profile ----
        if do_profiles
            save_profile(stepdir, sprintf('target_sigma%.6f_b%.6f', sigma_fixed, b), ...
                ic_end(:,2), target_feature);
        end

        % ---- Warm-start update for next b (optional) ----
        if use_b_warmstart
            ic0 = ic_rel;
        end
    end

    % ---- Save CSV ----
    hdr = {'b','target_feature','r_coh'};
    writecell(hdr, fullfile(outdir,'sweep_log.csv'));
    writematrix(log_mat, fullfile(outdir,'sweep_log.csv'), 'WriteMode','append');

    % ---- Overlay plot ----
    fig = figure('Visible','off');
    yyaxis left
    plot(log_mat(:,1), log_mat(:,2), 'o-', 'LineWidth', 1.2, 'DisplayName','target feature'); hold on;
    yline(r_thr,'--','r_{thr}');
    ylabel('target feature');

    yyaxis right
    plot(log_mat(:,1), log_mat(:,3), 's-', 'LineWidth', 1.2, 'DisplayName','coherence score');
    ylabel('r_{coh}');

    grid on;
    xlabel('b');
    title(sprintf('\\sigma=%.4f | target feature and coherence vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'metrics_overlay.png'));
    close(fig);

    % ---- target feature plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,2), 'o-'); hold on; yline(r_thr,'--'); hold off;
    grid on;
    xlabel('b'); ylabel('target feature');
    title(sprintf('\\sigma=%.4f | target feature vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'target_feature_vs_b.png'));
    close(fig);

    % ---- coherence plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,3), 'o-');
    grid on;
    xlabel('b'); ylabel('r_{coh}');
    title(sprintf('\\sigma=%.4f | temporal coherence vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'coherence_vs_b.png'));
    close(fig);

    fprintf('Saved: %s\n', outdir);
end

fprintf('\nAll done. Root output folder: %s\n', rootdir);

%% ---------------- Helpers ----------------
function r = temporal_coherence_score(V, use_last, frac)
    if isempty(V)
        r = NaN;
        return;
    end

    if use_last
        T = size(V,1);
        t0 = max(1, floor((1-frac)*T));
        V = V(t0:end,:);
    end

    nT = size(V,1);
    if nT < 2
        r = NaN;
        return;
    end

    diffs = zeros(nT-1,1);
    for k = 2:nT
        diffs(k-1) = norm(V(k,:) - V(k-1,:)) / (norm(V(k,:)) + 1e-12);
    end
    r = mean(diffs);
end

function save_spacetime(stepdir, tag, V, featval, rcoh)
    fig = figure('Visible','off');
    imagesc(V); axis tight; set(gca,'YDir','normal');
    colormap(gca,'sky'); colorbar;
    xlabel('space'); ylabel('time');
    title(sprintf('%s | feature=%.3e | rcoh=%.3e', tag, featval, rcoh), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('spacetime_%s.png', tag)));
    close(fig);
end

function save_profile(stepdir, tag, u2, featval)
    fig = figure('Visible','off');
    plot(u2, 'LineWidth', 1.2); grid on;
    xlabel('grid index'); ylabel('u_2');
    title(sprintf('%s | feature=%.3e', tag), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('profile_%s.png', tag)));
    close(fig);
end