function post_diagnostics_fast(mat_path, outdir_base, diag_stride)
% POST_DIAGNOSTICS_FAST
%
% Fast branch audit for a continuation .mat file.
% Uses cold local sigma sweeps only, and is branch-aware by filename:
%
%   right_up.mat / right_down.mat:
%       metric = rV_lr
%       threshold = 2.6e-3
%       center rule = first downward crossing
%
%   left_up.mat / left_down.mat:
%       metric = rV_ud
%       threshold = 0.5
%       center rule = strongest upward crossing
%
% For every diag_stride-th branch point:
%   1) do ONE cold local sigma sweep at fixed b_k
%   2) compute metrics on that sweep:
%        - rU      = final-profile symmetry
%        - rV_ud   = spacetime asymmetry under flipud(V)
%        - rV_lr   = spacetime asymmetry under fliplr(V)
%   3) choose center sigma from the branch-aware metric
%   4) save L/C/R spacetime plots and final profiles
%   5) write a CSV log
%
% Usage:
%   post_diagnostics_fast('right_up.mat');
%   post_diagnostics_fast('right_up.mat', 'diag_outputs');
%   post_diagnostics_fast('right_up.mat', 'diag_outputs', 5);

    if nargin < 2 || isempty(outdir_base)
        outdir_base = pwd;
    end
    if nargin < 3 || isempty(diag_stride)
        diag_stride = 5;
    end

    if ~isfile(mat_path)
        error('Cannot find mat file: %s', mat_path);
    end

    S = load(mat_path);

    % ---- infer branch side from filename ----
    [~, mat_name, ~] = fileparts(mat_path);
    mat_name_lower = lower(mat_name);

    is_left_branch  = contains(mat_name_lower, 'left_up')  || contains(mat_name_lower, 'left_down');
    is_right_branch = contains(mat_name_lower, 'right_up') || contains(mat_name_lower, 'right_down');

    if ~is_left_branch && ~is_right_branch
        warning('Could not infer branch side from filename "%s". Defaulting to right-side settings.', mat_name);
        is_right_branch = true;
    end

    % ---- find branch history ----
    Praw = [];
    if isfield(S,'p_history'),      Praw = S.p_history; end
    if isempty(Praw) && isfield(S,'p_hist_up'), Praw = S.p_hist_up; end
    if isempty(Praw) && isfield(S,'p_hist_dn'), Praw = S.p_hist_dn; end
    if isempty(Praw) && isfield(S,'p_hist'),    Praw = S.p_hist; end

    if isempty(Praw)
        error('No p_history / p_hist_up / p_hist_dn / p_hist found in %s', mat_path);
    end

    % Often first column duplicates start
    P = Praw(:,2:end);
    nPts = size(P,2);

    % ---- cold-start seed ----
    ic_file = 'ic_0.45_10.mat';
    Tseed = load(ic_file);
    if isfield(Tseed,'ic')
        ic_seed = Tseed.ic;
    else
        fn = fieldnames(Tseed);
        ic_seed = Tseed.(fn{1});
    end

    % ---- parameters: matched to continuation ----
    par = struct();
    par.T_relax    = 360;
    par.T_obs      = 600;
    par.dsigma     = 2.5e-4;
    par.half_width = 0.008;
    par.eps_trip   = 5e-4;

    if is_right_branch
        par.use_metric         = 'rV_lr';
        par.rthr               = 2.6e-3;
        par.crossing_direction = 'down';
        par.pick_mode          = 'first_crossing';
        branch_label           = 'RIGHT';
    else
        par.use_metric         = 'rV_ud';
        par.rthr               = 0.5;
        par.crossing_direction = 'up';
        par.pick_mode          = 'strongest_crossing';
        branch_label           = 'LEFT';
    end

    % ---- output folder ----
    tag = datestr(now,'yyyymmdd_HHMMSS');
    outdir = fullfile(outdir_base, sprintf('diag_fast_%s_%s', mat_name, tag));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    % ---- metrics ----
    rU_fun    = @(u) norm(u - flipud(u)) / max(norm(u), eps);
    rV_ud_fun = @(V) norm(V - flipud(V), 'fro') / max(norm(V,'fro'), eps);
    rV_lr_fun = @(V) norm(V - fliplr(V), 'fro') / max(norm(V,'fro'), eps);

    idxs = 1:diag_stride:nPts;

    fprintf('\n=== POST DIAGNOSTICS FAST ===\n');
    fprintf('MAT: %s\n', mat_path);
    fprintf('Branch side inferred: %s\n', branch_label);
    fprintf('Points in MAT: %d\n', nPts);
    fprintf('Points checked (stride=%d): %d\n', diag_stride, numel(idxs));
    fprintf('Metric used: %s\n', par.use_metric);
    fprintf('Threshold: %.6g\n', par.rthr);
    fprintf('Crossing direction: %s\n', par.crossing_direction);
    fprintf('Pick mode: %s\n', par.pick_mode);
    fprintf('Output dir: %s\n\n', outdir);

    log_rows = {};

    for ii = 1:numel(idxs)
        k = idxs(ii);
        sigma_k = P(1,k);
        b_k     = P(2,k);

        stepdir = fullfile(outdir, sprintf('k%04d_sigma%.6f_b%.6f', k, sigma_k, b_k));
        if ~exist(stepdir,'dir'); mkdir(stepdir); end

        fprintf('[%d/%d] k=%d  (sigma,b)=(%.6f, %.6f)\n', ...
            ii, numel(idxs), k, sigma_k, b_k);

        sig_grid = (sigma_k - par.half_width):par.dsigma:(sigma_k + par.half_width);
        M = numel(sig_grid);

        % =========================
        % 1) COLD local sigma sweep
        % =========================
        cold = init_metric_struct(M);
        for j = 1:M
            data = sim_obs(ic_seed, sig_grid(j), b_k, par.T_relax, par.T_obs);
            cold.rU(j)    = rU_fun(data.u_final);
            cold.rV_ud(j) = rV_ud_fun(data.V);
            cold.rV_lr(j) = rV_lr_fun(data.V);
        end

        % choose metric
        metric_cold = get_metric(cold, par.use_metric);

        % choose center sigma
        [sigC, pick_info] = choose_sigma(sig_grid, metric_cold, par);
        sigL = max(sig_grid(1), sigC - par.eps_trip);
        sigR = min(sig_grid(end), sigC + par.eps_trip);

        fprintf('    chosen center sigC=%.6f using %s | %s\n', ...
            sigC, par.use_metric, pick_info.kind);

        % =========================
        % 2) L/C/R diagnostics
        % =========================
        dataL = sim_obs(ic_seed, sigL, b_k, par.T_relax, par.T_obs);
        dataC = sim_obs(ic_seed, sigC, b_k, par.T_relax, par.T_obs);
        dataR = sim_obs(ic_seed, sigR, b_k, par.T_relax, par.T_obs);

        rUL    = rU_fun(dataL.u_final);    rUC    = rU_fun(dataC.u_final);    rUR    = rU_fun(dataR.u_final);
        rVudL  = rV_ud_fun(dataL.V);       rVudC  = rV_ud_fun(dataC.V);       rVudR  = rV_ud_fun(dataR.V);
        rVlrL  = rV_lr_fun(dataL.V);       rVlrC  = rV_lr_fun(dataC.V);       rVlrR  = rV_lr_fun(dataR.V);

        % =========================
        % 3) minimal plots
        % =========================

        % --- chosen metric sweep only ---
        fig1 = figure('Visible','off');
        plot(sig_grid, metric_cold, 'o-', 'LineWidth',1.0);
        hold on;
        yline(par.rthr,'--','r_{thr}');
        xline(sigC,'--','C'); xline(sigL,':','L'); xline(sigR,':','R');
        hold off;
        grid on;
        xlabel('\sigma');
        ylabel(par.use_metric);
        title(sprintf('%s | k=%d | b=%.6f | metric=%s', branch_label, k, b_k, par.use_metric), ...
              'Interpreter','none');
        exportgraphics(fig1, fullfile(stepdir, sprintf('chosen_metric_k%04d.png',k)));
        close(fig1);

        % --- L/C/R spacetime ---
        fig2 = figure('Visible','off');
        tl2 = tiledlayout(fig2,1,3,'Padding','compact','TileSpacing','compact');

        nexttile;
        imagesc(dataL.V); axis tight; set(gca,'YDir','normal'); colorbar;
        title(sprintf('L: \\sigma=%.6f\nrU=%.2e\nrVud=%.2e\nrVlr=%.2e', ...
            sigL, rUL, rVudL, rVlrL), 'Interpreter','tex');
        xlabel('space'); ylabel('time');

        nexttile;
        imagesc(dataC.V); axis tight; set(gca,'YDir','normal'); colorbar;
        title(sprintf('C: \\sigma=%.6f\nrU=%.2e\nrVud=%.2e\nrVlr=%.2e', ...
            sigC, rUC, rVudC, rVlrC), 'Interpreter','tex');
        xlabel('space'); ylabel('time');

        nexttile;
        imagesc(dataR.V); axis tight; set(gca,'YDir','normal'); colorbar;
        title(sprintf('R: \\sigma=%.6f\nrU=%.2e\nrVud=%.2e\nrVlr=%.2e', ...
            sigR, rUR, rVudR, rVlrR), 'Interpreter','tex');
        xlabel('space'); ylabel('time');

        title(tl2, sprintf('%s | L/C/R spacetime | k=%d | b=%.6f', branch_label, k, b_k), ...
              'Interpreter','tex');
        exportgraphics(fig2, fullfile(stepdir, sprintf('LCR_spacetime_k%04d.png',k)));
        close(fig2);

        % --- L/C/R final u2 profiles ---
        fig3 = figure('Visible','off');
        plot(dataL.u_final,'LineWidth',1.2); hold on;
        plot(dataC.u_final,'LineWidth',1.2);
        plot(dataR.u_final,'LineWidth',1.2);
        grid on;
        legend({'L','C','R'},'Location','best');
        title(sprintf('%s | Final u2 profiles | k=%d | b=%.6f', branch_label, k, b_k), ...
              'Interpreter','tex');
        exportgraphics(fig3, fullfile(stepdir, sprintf('LCR_profiles_k%04d.png',k)));
        close(fig3);

        % --- save mats ---
        save(fullfile(stepdir,'diag_L.mat'),'dataL','sigL','rUL','rVudL','rVlrL');
        save(fullfile(stepdir,'diag_C.mat'),'dataC','sigC','rUC','rVudC','rVlrC');
        save(fullfile(stepdir,'diag_R.mat'),'dataR','sigR','rUR','rVudR','rVlrR');

        % --- log row ---
        log_rows(end+1,:) = {k, sigma_k, b_k, sigL, sigC, sigR, ...
            pick_info.kind, pick_info.idx, ...
            rUL, rUC, rUR, rVudL, rVudC, rVudR, rVlrL, rVlrC, rVlrR}; %#ok<SAGROW>
    end

    Tlog = cell2table(log_rows, 'VariableNames', ...
        {'k','sigma_k','b_k','sigma_L','sigma_C','sigma_R', ...
         'pick_kind','pick_idx', ...
         'rU_L','rU_C','rU_R', ...
         'rVud_L','rVud_C','rVud_R', ...
         'rVlr_L','rVlr_C','rVlr_R'});
    writetable(Tlog, fullfile(outdir,'diag_fast_log.csv'));

    fprintf('\nDone. Saved diagnostics to:\n%s\n', outdir);
end


% ========================= helpers =========================

function out = init_metric_struct(M)
    out = struct();
    out.rU    = zeros(1,M);
    out.rV_ud = zeros(1,M);
    out.rV_lr = zeros(1,M);
end

function y = get_metric(S, name)
    switch name
        case 'rU'
            y = S.rU;
        case 'rV_ud'
            y = S.rV_ud;
        case 'rV_lr'
            y = S.rV_lr;
        otherwise
            error('Unknown metric: %s', name);
    end
end

function [sigC, info] = choose_sigma(sig_grid, y, par)
    info = struct('kind','','idx',NaN);

    above = (y >= par.rthr);

    switch par.crossing_direction
        case 'down'
            % right side: first downward crossing
            cross_all = find( above(1:end-1) & ~above(2:end) );

            if isempty(cross_all)
                [~, idx] = min(abs(y - par.rthr));
                sigC = sig_grid(idx);
                info.kind = 'closest_to_threshold';
                info.idx  = idx;
            else
                cross_idx = cross_all(1);
                sigC = 0.5*(sig_grid(cross_idx) + sig_grid(cross_idx+1));
                info.kind = 'first_downward_crossing';
                info.idx  = cross_idx;
            end

        case 'up'
            % left side: strongest upward crossing
            cross_all = find( ~above(1:end-1) & above(2:end) );

            if isempty(cross_all)
                [~, idx] = min(abs(y - par.rthr));
                sigC = sig_grid(idx);
                info.kind = 'closest_to_threshold';
                info.idx  = idx;
            else
                dy = diff(y);
                jump_sizes = dy(cross_all);   % positive jumps for upward crossings
                good = jump_sizes > 0;

                if ~any(good)
                    [~, idx] = min(abs(y - par.rthr));
                    sigC = sig_grid(idx);
                    info.kind = 'closest_to_threshold';
                    info.idx  = idx;
                else
                    cross_keep = cross_all(good);
                    jump_keep  = jump_sizes(good);
                    [~, jj] = max(jump_keep);
                    cross_idx = cross_keep(jj);

                    sigC = 0.5*(sig_grid(cross_idx) + sig_grid(cross_idx+1));
                    info.kind = 'strongest_upward_crossing';
                    info.idx  = cross_idx;
                end
            end

        otherwise
            error('Unknown crossing_direction: %s', par.crossing_direction);
    end
end

function data = sim_obs(ic0, sigma, b, T_relax, T_obs)
    par = struct('sigma', sigma, 'b', b);

    [ic_relaxed, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
    [ic_end,     ~, V] = solve_brusselator_1d(ic_relaxed, par, T_obs, 0);

    data.ic_end  = ic_end;
    data.u_final = ic_end(:,2);
    data.V       = V;
end
