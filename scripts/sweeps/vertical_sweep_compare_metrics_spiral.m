%% vertical_sweep_compare_metrics_spiral.m
% Vertical sweep in b for multiple fixed sigma values.
%
% Computes and compares:
%   (1) rProf_rel = final-profile symmetry of u2 after T_relax
%   (2) rV_space  = spacetime symmetry using fliplr(V)
%   (3) rV_time   = spacetime symmetry using flipud(V)
%   (4) r_coh     = temporal coherence score
%
% Intended use:
%   diagnose U/C/D structure for spiral lower / middle branches.

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
sigma_list = [0.43, 0.45, 0.57, 0.51];

b_start = 9.42;
b_end   = 9.52;
db      = 0.02;      % then refine to 0.02 or 0.01 if needed

T_relax = 360;

% diagnostic observation after relaxation
do_spacetime = true;
do_profiles  = true;
T_diag       = 360;

% use only the last fraction of V for rV metrics
use_last_window = true;
last_frac       = 0.25;

% IMPORTANT:
% false = independent b evaluations from same IC
% true  = warm-start along b
use_b_warmstart = false;

% same seed IC as continuation
ic_file = 'ic_0.45_10.mat';

% threshold lines for reference
r_thr_right = 2.6e-3;
r_thr_left  = 0.5;

timestamp = datestr(now,'yyyymmdd_HHMMSS');
rootdir = sprintf('vertSweep_spiral_compareMetrics_db%.0e_%s', db, timestamp);
if ~exist(rootdir,'dir'); mkdir(rootdir); end

%% ---------------- LOAD SEED IC ----------------
S = load(ic_file);
if isfield(S,'ic')
    ic_seed = S.ic;
else
    fn = fieldnames(S);
    ic_seed = S.(fn{1});
end

%% ---------------- METRICS ----------------
% final profile symmetry in u2(x)
r_profile = @(u) norm(u - flipud(u)) / max(norm(u), eps);

% spacetime metrics
rV_space_full = @(V) norm(V - fliplr(V), 'fro') / max(norm(V,'fro'), eps);
rV_time_full  = @(V) norm(V - flipud(V), 'fro') / max(norm(V,'fro'), eps);

rV_space = @(V) windowed_metric(V, rV_space_full, use_last_window, last_frac);
rV_time  = @(V) windowed_metric(V, rV_time_full,  use_last_window, last_frac);

%% ---------------- b GRID ----------------
bvals = b_start:db:b_end;
n = numel(bvals);

fprintf('Running sweeps: %d sigma-values, %d b-values each\n', numel(sigma_list), n);
fprintf('use_b_warmstart = %d\n', use_b_warmstart);

for ss = 1:numel(sigma_list)
    sigma_fixed = sigma_list(ss);

    outdir = fullfile(rootdir, sprintf('sigma_%.4f', sigma_fixed));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    % optional warm-start state along b
    ic0 = ic_seed;

    % columns:
    % b, rProf_rel, rV_space, rV_time, r_coh
    log_mat = nan(n, 5);

    fprintf('\n=== sigma = %.6f ===\n', sigma_fixed);

    for i = 1:n
        b = bvals(i);
        par = struct('sigma', sigma_fixed, 'b', b);

        fprintf('[sigma=%.4f] %d/%d b=%.6f\n', sigma_fixed, i, n, b);

        stepdir = fullfile(outdir, sprintf('b_%.6f', b));
        if (do_spacetime || do_profiles) && ~exist(stepdir,'dir')
            mkdir(stepdir);
        end

        % ---- Relaxation ----
        ic_rel_start = ic_seed;
        if use_b_warmstart
            ic_rel_start = ic0;
        end

        [ic_rel, ~, ~] = solve_brusselator_1d(ic_rel_start, par, T_relax, 0);

        % profile metric after T_relax
        u2_rel = ic_rel(:,2);
        rProf_rel = r_profile(u2_rel);

        % ---- Diagnostic observation from relaxed state ----
        [ic_end, ~, V_diag] = solve_brusselator_1d(ic_rel, par, T_diag, 0);

        rVsp = rV_space(V_diag);
        rVtm = rV_time(V_diag);
        rcoh = temporal_coherence_score(V_diag, use_last_window, last_frac);

        log_mat(i,:) = [b, rProf_rel, rVsp, rVtm, rcoh];

        % ---- Save spacetime ----
        if do_spacetime
            save_spacetime(stepdir, sprintf('diag_sigma%.6f_b%.6f', sigma_fixed, b), ...
                V_diag, rProf_rel, rVsp, rVtm, rcoh);
        end

        % ---- Save final profile ----
        if do_profiles
            save_profile(stepdir, sprintf('diag_sigma%.6f_b%.6f', sigma_fixed, b), ...
                ic_end(:,2), rProf_rel);
        end

        % ---- Warm-start update for next b (optional) ----
        if use_b_warmstart
            ic0 = ic_rel;
        end
    end

    % ---- Save CSV ----
    hdr = {'b','rProf_rel','rV_space','rV_time','rcoh'};
    writecell(hdr, fullfile(outdir,'sweep_log.csv'));
    writematrix(log_mat, fullfile(outdir,'sweep_log.csv'), 'WriteMode','append');

    % ---- Overlay plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,2), 'o-', 'LineWidth', 1.2, 'DisplayName','rProf\_rel'); hold on;
    plot(log_mat(:,1), log_mat(:,3), 's-', 'LineWidth', 1.2, 'DisplayName','rV\_space');
    plot(log_mat(:,1), log_mat(:,4), 'd-', 'LineWidth', 1.2, 'DisplayName','rV\_time');
    yyaxis right
    plot(log_mat(:,1), log_mat(:,5), '^-', 'LineWidth', 1.2, 'DisplayName','rcoh');
    ylabel('rcoh');
    yyaxis left
    yline(r_thr_right,'--','r_{thr,right}');
    yline(r_thr_left,'--','r_{thr,left}');
    grid on;
    xlabel('b');
    ylabel('metric value');
    title(sprintf('\\sigma=%.4f | spiral metrics vs b', sigma_fixed));
    legend('Location','best');
    exportgraphics(fig, fullfile(outdir,'metrics_overlay.png'));
    close(fig);

    % ---- rProf plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,2), 'o-'); grid on;
    xlabel('b'); ylabel('rProf\_rel');
    title(sprintf('\\sigma=%.4f | rProf\\_rel vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'rProf_rel_vs_b.png'));
    close(fig);

    % ---- rV_space plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,3), 'o-'); hold on; yline(r_thr_right,'--'); hold off;
    grid on;
    xlabel('b'); ylabel('rV\_space');
    title(sprintf('\\sigma=%.4f | rV\\_space vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'rV_space_vs_b.png'));
    close(fig);

    % ---- rV_time plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,4), 'o-'); hold on; yline(r_thr_left,'--'); hold off;
    grid on;
    xlabel('b'); ylabel('rV\_time');
    title(sprintf('\\sigma=%.4f | rV\\_time vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'rV_time_vs_b.png'));
    close(fig);

    % ---- coherence plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,5), 'o-');
    grid on;
    xlabel('b'); ylabel('rcoh');
    title(sprintf('\\sigma=%.4f | temporal coherence vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'coherence_vs_b.png'));
    close(fig);

    fprintf('Saved: %s\n', outdir);
end

fprintf('\nAll done. Root output folder: %s\n', rootdir);

%% ---------------- Helpers ----------------
function r = windowed_metric(V, rfun, use_last, frac)
    if ~use_last
        r = rfun(V);
        return;
    end
    T = size(V,1);
    t0 = max(1, floor((1-frac)*T));
    Vuse = V(t0:end,:);
    r = rfun(Vuse);
end

function rcoh = temporal_coherence_score(V, use_last, frac)
    if isempty(V)
        rcoh = NaN;
        return;
    end

    if use_last
        T = size(V,1);
        t0 = max(1, floor((1-frac)*T));
        V = V(t0:end,:);
    end

    nT = size(V,1);
    if nT < 2
        rcoh = NaN;
        return;
    end

    diffs = zeros(nT-1,1);
    for k = 2:nT
        diffs(k-1) = norm(V(k,:) - V(k-1,:)) / (norm(V(k,:)) + 1e-12);
    end
    rcoh = mean(diffs);
end

function save_spacetime(stepdir, tag, V, rProf, rVsp, rVtm, rcoh)
    fig = figure('Visible','off');
    imagesc(V); axis tight; set(gca,'YDir','normal');
    colormap(gca,'sky'); colorbar;
    xlabel('space'); ylabel('time');
    title(sprintf('%s | rProf=%.2e | rVsp=%.2e | rVtm=%.2e | rcoh=%.2e', ...
        tag, rProf, rVsp, rVtm, rcoh), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('spacetime_%s.png', tag)));
    close(fig);
end

function save_profile(stepdir, tag, u2, rProf)
    fig = figure('Visible','off');
    plot(u2, 'LineWidth', 1.2); grid on;
    xlabel('grid index'); ylabel('u_2');
    title(sprintf('%s | rProf=%.2e', tag), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('profile_%s.png', tag)));
    close(fig);
end