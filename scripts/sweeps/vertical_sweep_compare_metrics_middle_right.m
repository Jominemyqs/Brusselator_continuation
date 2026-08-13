%% vertical_sweep_compare_metrics_middle_right.m
% Compare several candidate metrics for the MIDDLE-RIGHT branch
% by fixing sigma and sweeping in b (U/C/D logic).
%
% Supports two modes:
%   seed_type = 'target'
%   seed_type = 'spiral'
%
% For each fixed sigma:
%   - relax from a fixed seed
%   - observe spacetime V
%   - compute candidate metrics
%   - save spacetime plots
%   - save final profiles
%   - save csv + overlay plots
%
% Goal:
%   identify which metric transitions at the visually correct b-value,
%   so that middle-right continuation can use that metric instead of
%   drifting downward.

clear; clc; close all;

%% ================= USER SETTINGS =================
seed_type = 'target';   % 'target' or 'spiral'

% Suggested first tests:
% target: [0.42 0.44 0.46 0.48 0.50 0.51]
% spiral: [0.42 0.44 0.46 0.48 0.50 0.51]
sigma_list = [0.512, 0.515];

b_start = 9.40;
b_end   = 9.44;
db      = 0.01;

% independent evaluations are best for diagnostics
use_b_warmstart = false;

% simulation times
switch lower(seed_type)
    case 'target'
        T_relax = 240;
        T_diag  = 240;
        ic_file = 'half_target_seed.mat';
    case 'spiral'
        T_relax = 360;
        T_diag  = 360;
        ic_file = 'ic_0.45_10.mat';
    otherwise
        error('seed_type must be ''target'' or ''spiral''.');
end

% save diagnostics
do_spacetime = true;
do_profiles  = true;

% use last fraction of spacetime for V-based metrics
use_last_window = true;
last_frac       = 0.25;

timestamp = datestr(now,'yyyymmdd_HHMMSS');
rootdir = sprintf('vertSweep_middleRight_%s_db%.0e_%s', seed_type, db, timestamp);
if ~exist(rootdir,'dir'); mkdir(rootdir); end

fprintf('seed_type = %s\n', seed_type);
fprintf('sigma count = %d | b count = %d\n', numel(sigma_list), numel(b_start:db:b_end));
fprintf('use_b_warmstart = %d\n', use_b_warmstart);

%% ================= LOAD SEED =================
S = load(ic_file);
if isfield(S,'ic')
    ic_seed = S.ic;
elseif isfield(S,'ic_end')
    ic_seed = S.ic_end;
else
    fn = fieldnames(S);
    ic_seed = S.(fn{1});
end

bvals = b_start:db:b_end;
nb = numel(bvals);

%% ================= MAIN LOOP =================
for ss = 1:numel(sigma_list)
    sigma_fixed = sigma_list(ss);

    outdir = fullfile(rootdir, sprintf('sigma_%.4f', sigma_fixed));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    ic0 = ic_seed;

    switch lower(seed_type)
        case 'target'
            % columns:
            % b, target_feature, Sc, Tr, rcoh, gate_pass
            log_mat = nan(nb, 6);
            headers = {'b','target_feature','Sc','Tr','rcoh','gate_pass'};
        case 'spiral'
            % columns:
            % b, rProf_rel, rV_space, rV_time, rcoh, gate_pass
            log_mat = nan(nb, 6);
            headers = {'b','rProf_rel','rV_space','rV_time','rcoh','gate_pass'};
    end

    fprintf('\n=== sigma = %.6f ===\n', sigma_fixed);

    for i = 1:nb
        b = bvals(i);
        fprintf('[%s] sigma=%.4f | %d/%d | b=%.6f\n', seed_type, sigma_fixed, i, nb, b);

        stepdir = fullfile(outdir, sprintf('b_%.6f', b));
        if (do_spacetime || do_profiles) && ~exist(stepdir,'dir')
            mkdir(stepdir);
        end

        par = struct('sigma', sigma_fixed, 'b', b);

        % ----- relaxation -----
        ic_rel_start = ic_seed;
        if use_b_warmstart
            ic_rel_start = ic0;
        end

        [ic_rel, ~, ~] = solve_brusselator_1d(ic_rel_start, par, T_relax, 0);

        % ----- observation -----
        [ic_end, ~, V_diag] = solve_brusselator_1d(ic_rel, par, T_diag, 0);

        % ----- common metrics -----
        rcoh = temporal_coherence_score(V_diag);
        gate_pass = double(rcoh < 0.11);

        switch lower(seed_type)
            case 'target'
                [target_feature, Sc, Tr] = compute_target_metrics(V_diag);
                log_mat(i,:) = [b, target_feature, Sc, Tr, rcoh, gate_pass];

                if do_spacetime
                    save_spacetime_target(stepdir, sigma_fixed, b, V_diag, target_feature, Sc, Tr, rcoh);
                end
                if do_profiles
                    save_profile(stepdir, sprintf('target_sigma%.6f_b%.6f', sigma_fixed, b), ic_end(:,2), target_feature);
                end

            case 'spiral'
                [rProf_rel, rV_space, rV_time] = compute_spiral_metrics(ic_rel, V_diag, use_last_window, last_frac);
                log_mat(i,:) = [b, rProf_rel, rV_space, rV_time, rcoh, gate_pass];

                if do_spacetime
                    save_spacetime_spiral(stepdir, sigma_fixed, b, V_diag, rProf_rel, rV_space, rV_time, rcoh);
                end
                if do_profiles
                    save_profile(stepdir, sprintf('spiral_sigma%.6f_b%.6f', sigma_fixed, b), ic_end(:,2), rProf_rel);
                end
        end

        if use_b_warmstart
            ic0 = ic_rel;
        end
    end

    % ----- save csv -----
    writecell(headers, fullfile(outdir, 'sweep_log.csv'));
    writematrix(log_mat, fullfile(outdir, 'sweep_log.csv'), 'WriteMode', 'append');

    % ----- plots -----
    switch lower(seed_type)
        case 'target'
            make_target_plots(outdir, sigma_fixed, log_mat);
        case 'spiral'
            make_spiral_plots(outdir, sigma_fixed, log_mat);
    end

    fprintf('Saved: %s\n', outdir);
end

fprintf('\nAll done. Root folder: %s\n', rootdir);

%% ================= HELPERS =================

function [target_feature, Sc, Tr] = compute_target_metrics(V)
    if isempty(V)
        target_feature = NaN;
        Sc = NaN;
        Tr = NaN;
        return;
    end

    [nt, nx] = size(V);

    % last 40% in time, matching current target feature
    t0 = max(1, round(0.60 * nt));
    Vuse = V(t0:end, :);

    % core window
    j_core1 = max(1, round(0.40 * nx));
    j_core2 = min(nx, round(0.60 * nx));
    Vcore = Vuse(:, j_core1:j_core2);

    % tail window
    j_tail1 = max(1, round(0.80 * nx));
    j_tail2 = nx;
    Vtail = Vuse(:, j_tail1:j_tail2);

    Sc = mean(var(Vcore, 0, 2));
    Tr = mean(var(Vtail, 0, 1));
    target_feature = min(Sc, Tr);
end

function [rProf_rel, rV_space, rV_time] = compute_spiral_metrics(ic_rel, V, use_last_window, last_frac)
    u_final = ic_rel(:,2);
    denp = norm(u_final);
    if denp < eps
        rProf_rel = 0;
    else
        rProf_rel = norm(u_final - flipud(u_final)) / denp;
    end

    if isempty(V)
        rV_space = NaN;
        rV_time  = NaN;
        return;
    end

    if use_last_window
        T = size(V,1);
        t0 = max(1, floor((1-last_frac)*T));
        Vuse = V(t0:end,:);
    else
        Vuse = V;
    end

    den = norm(Vuse,'fro');
    if den < eps
        rV_space = 0;
        rV_time  = 0;
    else
        rV_space = norm(Vuse - fliplr(Vuse), 'fro') / den;
        rV_time  = norm(Vuse - flipud(Vuse), 'fro') / den;
    end
end

function r_coh = temporal_coherence_score(V)
    if isempty(V)
        r_coh = Inf;
        return;
    end

    nT = size(V,1);
    if nT < 2
        r_coh = Inf;
        return;
    end

    n_keep = min(80, nT);
    Vlate = V(end-n_keep+1:end, :);

    diffs = zeros(n_keep-1,1);
    for k = 2:n_keep
        diffs(k-1) = norm(Vlate(k,:) - Vlate(k-1,:)) / (norm(Vlate(k,:)) + 1e-12);
    end

    r_coh = mean(diffs);
end

function save_spacetime_target(stepdir, sigma, b, V, featval, Sc, Tr, rcoh)
    fig = figure('Visible','off');
    imagesc(V); axis tight; set(gca,'YDir','normal');
    colormap(gca,'sky'); colorbar;
    xlabel('space'); ylabel('time');
    title(sprintf('target sigma%.6f b%.6f | feat=%.3e | Sc=%.3e | Tr=%.3e | rcoh=%.3e', ...
        sigma, b, featval, Sc, Tr, rcoh), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('spacetime_target_sigma%.6f_b%.6f.png', sigma, b)));
    close(fig);
end

function save_spacetime_spiral(stepdir, sigma, b, V, rProf, rVsp, rVtm, rcoh)
    fig = figure('Visible','off');
    imagesc(V); axis tight; set(gca,'YDir','normal');
    colormap(gca,'sky'); colorbar;
    xlabel('space'); ylabel('time');
    title(sprintf('spiral sigma%.6f b%.6f | rProf=%.3e | rVsp=%.3e | rVtm=%.3e | rcoh=%.3e', ...
        sigma, b, rProf, rVsp, rVtm, rcoh), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('spacetime_spiral_sigma%.6f_b%.6f.png', sigma, b)));
    close(fig);
end

function save_profile(stepdir, tag, u2, metricval)
    fig = figure('Visible','off');
    plot(u2, 'LineWidth', 1.2); grid on;
    xlabel('grid index'); ylabel('u_2');
    title(sprintf('%s | metric=%.3e', tag, metricval), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('profile_%s.png', tag)));
    close(fig);
end

function make_target_plots(outdir, sigma_fixed, log_mat)
    b = log_mat(:,1);
    feat = log_mat(:,2);
    Sc = log_mat(:,3);
    Tr = log_mat(:,4);
    rcoh = log_mat(:,5);
    gate = log_mat(:,6);

    % overlay
    fig = figure('Visible','off');
    yyaxis left
    plot(b, feat, 'o-', 'LineWidth', 1.2, 'DisplayName','target feature'); hold on;
    plot(b, Sc,   's-', 'LineWidth', 1.2, 'DisplayName','Sc');
    plot(b, Tr,   'd-', 'LineWidth', 1.2, 'DisplayName','Tr');
    yline(10,'--','r_{thr}');
    ylabel('feature value');

    yyaxis right
    plot(b, rcoh, '^-', 'LineWidth', 1.2, 'DisplayName','rcoh');
    yline(0.11, '--', 'coh_{thr}');
    ylabel('rcoh');

    grid on;
    xlabel('b');
    title(sprintf('\\sigma=%.4f | target metrics vs b', sigma_fixed));
    legend('Location','best');
    exportgraphics(fig, fullfile(outdir,'metrics_overlay.png'));
    close(fig);

    % target feature
    fig = figure('Visible','off');
    plot(b, feat, 'o-'); hold on; yline(10,'--'); hold off;
    grid on; xlabel('b'); ylabel('target feature');
    title(sprintf('\\sigma=%.4f | target feature vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'target_feature_vs_b.png'));
    close(fig);

    % Sc
    fig = figure('Visible','off');
    plot(b, Sc, 'o-'); grid on;
    xlabel('b'); ylabel('Sc');
    title(sprintf('\\sigma=%.4f | Sc vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'Sc_vs_b.png'));
    close(fig);

    % Tr
    fig = figure('Visible','off');
    plot(b, Tr, 'o-'); grid on;
    xlabel('b'); ylabel('Tr');
    title(sprintf('\\sigma=%.4f | Tr vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'Tr_vs_b.png'));
    close(fig);

    % coherence
    fig = figure('Visible','off');
    plot(b, rcoh, 'o-'); hold on; yline(0.11,'--'); hold off;
    grid on; xlabel('b'); ylabel('rcoh');
    title(sprintf('\\sigma=%.4f | temporal coherence vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'coherence_vs_b.png'));
    close(fig);

    % gate pass
    fig = figure('Visible','off');
    stairs(b, gate, 'LineWidth', 1.2); ylim([-0.1 1.1]); grid on;
    xlabel('b'); ylabel('gate pass');
    title(sprintf('\\sigma=%.4f | gate pass vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'gate_pass_vs_b.png'));
    close(fig);
end

function make_spiral_plots(outdir, sigma_fixed, log_mat)
    b = log_mat(:,1);
    rProf = log_mat(:,2);
    rVsp  = log_mat(:,3);
    rVtm  = log_mat(:,4);
    rcoh  = log_mat(:,5);
    gate  = log_mat(:,6);

    fig = figure('Visible','off');
    yyaxis left
    plot(b, rProf, 'o-', 'LineWidth', 1.2, 'DisplayName','rProf\_rel'); hold on;
    plot(b, rVsp,  's-', 'LineWidth', 1.2, 'DisplayName','rV\_space');
    plot(b, rVtm,  'd-', 'LineWidth', 1.2, 'DisplayName','rV\_time');
    yline(2.6e-3,'--','r_{thr,right}');
    yline(0.5,'--','r_{thr,left}');
    ylabel('metric value');

    yyaxis right
    plot(b, rcoh, '^-', 'LineWidth', 1.2, 'DisplayName','rcoh');
    yline(0.11,'--','coh_{thr}');
    ylabel('rcoh');

    grid on;
    xlabel('b');
    title(sprintf('\\sigma=%.4f | spiral metrics vs b', sigma_fixed));
    legend('Location','best');
    exportgraphics(fig, fullfile(outdir,'metrics_overlay.png'));
    close(fig);

    fig = figure('Visible','off');
    plot(b, rProf, 'o-'); grid on;
    xlabel('b'); ylabel('rProf\_rel');
    title(sprintf('\\sigma=%.4f | rProf\\_rel vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'rProf_rel_vs_b.png'));
    close(fig);

    fig = figure('Visible','off');
    plot(b, rVsp, 'o-'); hold on; yline(2.6e-3,'--'); hold off;
    grid on; xlabel('b'); ylabel('rV\_space');
    title(sprintf('\\sigma=%.4f | rV\\_space vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'rV_space_vs_b.png'));
    close(fig);

    fig = figure('Visible','off');
    plot(b, rVtm, 'o-'); hold on; yline(0.5,'--'); hold off;
    grid on; xlabel('b'); ylabel('rV\_time');
    title(sprintf('\\sigma=%.4f | rV\\_time vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'rV_time_vs_b.png'));
    close(fig);

    fig = figure('Visible','off');
    plot(b, rcoh, 'o-'); hold on; yline(0.11,'--'); hold off;
    grid on; xlabel('b'); ylabel('rcoh');
    title(sprintf('\\sigma=%.4f | temporal coherence vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'coherence_vs_b.png'));
    close(fig);

    fig = figure('Visible','off');
    stairs(b, gate, 'LineWidth', 1.2); ylim([-0.1 1.1]); grid on;
    xlabel('b'); ylabel('gate pass');
    title(sprintf('\\sigma=%.4f | gate pass vs b', sigma_fixed));
    exportgraphics(fig, fullfile(outdir,'gate_pass_vs_b.png'));
    close(fig);
end