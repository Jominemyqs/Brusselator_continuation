%% horizontal_sweep_compare_metrics.m
% Horizontal sweep in sigma for multiple fixed b values.
%
% Computes and compares:
%   (1) rProf_rel = symmetry of u2 after T_relax
%   (2) rV_space  = spacetime symmetry using fliplr(V)  [space reflection]
%   (3) rV_time   = spacetime symmetry using flipud(V)  [time reflection]
%
% Optional:
%   - cold-start or warm-start across sigma
%   - save spacetime plots
%   - save final u2 profiles
%
% Output per b:
%   - sweep_log.csv
%   - metrics_overlay.png
%   - rProf_rel_vs_sigma.png
%   - rV_space_vs_sigma.png
%   - rV_time_vs_sigma.png
%   - spacetime/profile pngs per sigma (if enabled)

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
b_list = [9.46];   % try one or several b values

sigma_start = 0.516;
sigma_end   = 0.520;
ds_sigma    = 0.002;             % try 1e-3 first, then 5e-4, then 2.5e-4 if needed

T_relax = 360;

% diagnostic observation after relaxation
do_spacetime = true;
do_profiles  = true;
T_diag       = 360;

% use only the last fraction of V for rV metrics
use_last_window = true;
last_frac       = 0.25;

% IMPORTANT:
% false = independent sigma evaluations from same IC
% true  = warm-start along sigma
use_sigma_warmstart = false;

% same seed IC as continuation
ic_file = 'ic_0.45_10.mat';

% threshold line for reference in plots
r_thr = 2.6e-3;

timestamp = datestr(now,'yyyymmdd_HHMMSS');
rootdir = sprintf('horizSweep_compareMetrics_ds%.0e_%s', ds_sigma, timestamp);
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

%% ---------------- SIGMA GRID ----------------
sigmas = sigma_start:ds_sigma:sigma_end;
n = numel(sigmas);

fprintf('Running sweeps: %d b-values, %d sigmas each\n', numel(b_list), n);
fprintf('use_sigma_warmstart = %d\n', use_sigma_warmstart);

for bb = 1:numel(b_list)
    b_fixed = b_list(bb);

    outdir = fullfile(rootdir, sprintf('b_%.4f', b_fixed));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    % optional warm-start state along sigma
    ic0 = ic_seed;

    % columns:
    % sigma, rProf_rel, rV_space, rV_time
    log_mat = nan(n, 4);

    fprintf('\n=== b = %.6f ===\n', b_fixed);

    for i = 1:n
        sigma = sigmas(i);
        par = struct('sigma', sigma, 'b', b_fixed);

        fprintf('[b=%.4f] %d/%d sigma=%.6f\n', b_fixed, i, n, sigma);

        stepdir = fullfile(outdir, sprintf('sigma_%.6f', sigma));
        if (do_spacetime || do_profiles) && ~exist(stepdir,'dir')
            mkdir(stepdir);
        end

        % ---- Relaxation ----
        ic_rel_start = ic_seed;
        if use_sigma_warmstart
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

        log_mat(i,:) = [sigma, rProf_rel, rVsp, rVtm];

        % ---- Save spacetime ----
        if do_spacetime
            save_spacetime(stepdir, sprintf('diag_sigma%.6f', sigma), V_diag, rVsp, rVtm);
        end

        % ---- Save final profile ----
        if do_profiles
            save_profile(stepdir, sprintf('diag_sigma%.6f', sigma), ic_end(:,2), rProf_rel);
        end

        % ---- Warm-start update for next sigma (optional) ----
        if use_sigma_warmstart
            ic0 = ic_rel;
        end
    end

    % ---- Save CSV ----
    hdr = {'sigma','rProf_rel','rV_space','rV_time'};
    writecell(hdr, fullfile(outdir,'sweep_log.csv'));
    writematrix(log_mat, fullfile(outdir,'sweep_log.csv'), 'WriteMode','append');

    % ---- Overlay plot (all metrics together) ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,2), 'o-', 'LineWidth', 1.2, 'DisplayName','rProf\_rel'); hold on;
    plot(log_mat(:,1), log_mat(:,3), 's-', 'LineWidth', 1.2, 'DisplayName','rV\_space');
    plot(log_mat(:,1), log_mat(:,4), 'd-', 'LineWidth', 1.2, 'DisplayName','rV\_time');
    yline(r_thr,'--','r_{thr}');
    grid on;
    xlabel('\sigma');
    ylabel('metric value');
    title(sprintf('b=%.4f | all metrics vs sigma', b_fixed));
    legend('Location','best');
    exportgraphics(fig, fullfile(outdir,'metrics_overlay.png'));
    close(fig);

    % ---- rProf plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,2), 'o-'); hold on; yline(r_thr,'--'); hold off;
    grid on;
    xlabel('\sigma'); ylabel('rProf\_rel');
    title(sprintf('b=%.4f | rProf\\_rel vs sigma', b_fixed));
    exportgraphics(fig, fullfile(outdir,'rProf_rel_vs_sigma.png'));
    close(fig);

    % ---- rV_space plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,3), 'o-'); hold on; yline(r_thr,'--'); hold off;
    grid on;
    xlabel('\sigma'); ylabel('rV\_space');
    title(sprintf('b=%.4f | rV\\_space vs sigma', b_fixed));
    exportgraphics(fig, fullfile(outdir,'rV_space_vs_sigma.png'));
    close(fig);

    % ---- rV_time plot ----
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,4), 'o-'); hold on; yline(r_thr,'--'); hold off;
    grid on;
    xlabel('\sigma'); ylabel('rV\_time');
    title(sprintf('b=%.4f | rV\\_time vs sigma', b_fixed));
    exportgraphics(fig, fullfile(outdir,'rV_time_vs_sigma.png'));
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

function save_spacetime(stepdir, tag, V, rVsp, rVtm)
    fig = figure('Visible','off');
    imagesc(V); axis tight; set(gca,'YDir','normal');
    colormap(gca,'sky'); colorbar;
    xlabel('space'); ylabel('time');
    title(sprintf('%s | rVspace=%.2e | rVtime=%.2e', tag, rVsp, rVtm), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('spacetime_%s.png', tag)));
    close(fig);
end

function save_profile(stepdir, tag, u2, rProf)
    fig = figure('Visible','off');
    plot(u2, 'LineWidth', 1.2); grid on;
    xlabel('grid index'); ylabel('u_2');
    title(sprintf('%s | rProf=%.2e', tag, rProf), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('profile_%s.png', tag)));
    close(fig);
end