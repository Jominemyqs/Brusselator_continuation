function post_diagnostics_fast_target(mat_path, outdir_base, diag_stride)
% POST_DIAGNOSTICS_FAST_TARGET
%
% Fast branch audit for TARGET-pattern continuation.
%
% Branch-aware by filename:
%   target_right_up.mat / target_right_down.mat:
%       crossing = downward through r_thr
%
%   target_left_up.mat / target_left_down.mat:
%       crossing = upward through r_thr
%
% Uses:
%   - feature_evaluation_target.m
%   - cold local sigma sweeps only
%   - cold L/C/R diagnostics
%
% Usage:
%   post_diagnostics_fast_target('target_right_up.mat');
%   post_diagnostics_fast_target('target_right_up.mat', 'diag_outputs');
%   post_diagnostics_fast_target('target_right_up.mat', 'diag_outputs', 5);

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

    is_left_branch  = contains(mat_name_lower, 'target_left_up')  || contains(mat_name_lower, 'target_left_down');
    is_right_branch = contains(mat_name_lower, 'target_right_up') || contains(mat_name_lower, 'target_right_down');

    if ~is_left_branch && ~is_right_branch
        warning('Could not infer target branch side from filename "%s". Defaulting to right-side settings.', mat_name);
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

        % ---- load target seed ----
    seed_file = 'half_target_seed.mat';
    if ~isfile(seed_file)
        error('Could not find target seed file: %s', seed_file);
    end

    Tseed = load(seed_file);
    if isfield(Tseed,'ic')
        ic_seed = Tseed.ic;
    elseif isfield(Tseed,'ic_end')
        ic_seed = Tseed.ic_end;
    else
        fn = fieldnames(Tseed);
        ic_seed = Tseed.(fn{1});
    end
    % ---- parameters ----
    par = struct();
    par.T_relax    = 240;
    par.T_obs      = 240;
    par.dsigma     = 1e-3;
    par.half_width = 0.020;
    par.eps_trip   = 1e-3;
    par.rthr       = 8;

    if is_right_branch
        par.crossing_direction = 'down';
        branch_label = 'TARGET-RIGHT';
    else
        par.crossing_direction = 'up';
        branch_label = 'TARGET-LEFT';
    end

    % ---- output folder ----
    tag = datestr(now,'yyyymmdd_HHMMSS');
    outdir = fullfile(outdir_base, sprintf('diag_fast_target_%s_%s', mat_name, tag));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    idxs = 1:diag_stride:nPts;

    fprintf('\n=== POST DIAGNOSTICS FAST TARGET ===\n');
    fprintf('MAT: %s\n', mat_path);
    fprintf('Branch side inferred: %s\n', branch_label);
    fprintf('Points in MAT: %d\n', nPts);
    fprintf('Points checked (stride=%d): %d\n', diag_stride, numel(idxs));
    fprintf('Threshold: %.6g\n', par.rthr);
    fprintf('Crossing direction: %s\n', par.crossing_direction);
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

        feat_vals = zeros(1,M);

        % =========================
        % 1) COLD local sigma sweep
        % =========================
        for j = 1:M
            feat_vals(j) = eval_target_feature(ic_seed, sig_grid(j), b_k, par.T_relax, par.T_obs);
        end

        % choose center sigma
        [sigC, pick_info] = choose_sigma(sig_grid, feat_vals, par);
        sigL = max(sig_grid(1), sigC - par.eps_trip);
        sigR = min(sig_grid(end), sigC + par.eps_trip);

        fprintf('    chosen center sigC=%.6f | %s\n', sigC, pick_info.kind);

        % =========================
        % 2) L/C/R diagnostics
        % =========================
        dataL = sim_obs(ic_seed, sigL, b_k, par.T_relax, par.T_obs);
        dataC = sim_obs(ic_seed, sigC, b_k, par.T_relax, par.T_obs);
        dataR = sim_obs(ic_seed, sigR, b_k, par.T_relax, par.T_obs);

        featL = eval_target_feature(ic_seed, sigL, b_k, par.T_relax, par.T_obs);
        featC = eval_target_feature(ic_seed, sigC, b_k, par.T_relax, par.T_obs);
        featR = eval_target_feature(ic_seed, sigR, b_k, par.T_relax, par.T_obs);

        % optional extra profile metric for sanity
        rU_fun = @(u) norm(u - flipud(u)) / max(norm(u), eps);
        rUL = rU_fun(dataL.u_final);
        rUC = rU_fun(dataC.u_final);
        rUR = rU_fun(dataR.u_final);

        % =========================
        % 3) plots
        % =========================

        % --- target feature sweep ---
        fig1 = figure('Visible','off');
        plot(sig_grid, feat_vals, 'o-', 'LineWidth',1.0); hold on;
        yline(par.rthr,'--','r_{thr}');
        xline(sigC,'--','C'); xline(sigL,':','L'); xline(sigR,':','R');
        hold off;
        grid on;
        xlabel('\sigma');
        ylabel('target feature');
        title(sprintf('%s | k=%d | b=%.6f', branch_label, k, b_k), 'Interpreter','none');
        exportgraphics(fig1, fullfile(stepdir, sprintf('target_feature_k%04d.png',k)));
        close(fig1);

        % --- L/C/R spacetime ---
        fig2 = figure('Visible','off');
        tl2 = tiledlayout(fig2,1,3,'Padding','compact','TileSpacing','compact');

        nexttile;
        imagesc(dataL.V); axis tight; set(gca,'YDir','normal'); colorbar;
        title(sprintf('L: \\sigma=%.6f\nfeat=%.2e\nrU=%.2e', sigL, featL, rUL), 'Interpreter','tex');
        xlabel('space'); ylabel('time');

        nexttile;
        imagesc(dataC.V); axis tight; set(gca,'YDir','normal'); colorbar;
        title(sprintf('C: \\sigma=%.6f\nfeat=%.2e\nrU=%.2e', sigC, featC, rUC), 'Interpreter','tex');
        xlabel('space'); ylabel('time');

        nexttile;
        imagesc(dataR.V); axis tight; set(gca,'YDir','normal'); colorbar;
        title(sprintf('R: \\sigma=%.6f\nfeat=%.2e\nrU=%.2e', sigR, featR, rUR), 'Interpreter','tex');
        xlabel('space'); ylabel('time');

        title(tl2, sprintf('%s | L/C/R spacetime | k=%d | b=%.6f', branch_label, k, b_k), ...
            'Interpreter','tex');
        exportgraphics(fig2, fullfile(stepdir, sprintf('LCR_spacetime_k%04d.png',k)));
        close(fig2);

        % --- L/C/R final profiles ---
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
        save(fullfile(stepdir,'diag_L.mat'),'dataL','sigL','featL','rUL');
        save(fullfile(stepdir,'diag_C.mat'),'dataC','sigC','featC','rUC');
        save(fullfile(stepdir,'diag_R.mat'),'dataR','sigR','featR','rUR');

        % --- log row ---
        log_rows(end+1,:) = {k, sigma_k, b_k, sigL, sigC, sigR, ...
            pick_info.kind, pick_info.idx, featL, featC, featR, rUL, rUC, rUR}; %#ok<SAGROW>
    end

    Tlog = cell2table(log_rows, 'VariableNames', ...
        {'k','sigma_k','b_k','sigma_L','sigma_C','sigma_R', ...
         'pick_kind','pick_idx', ...
         'feat_L','feat_C','feat_R', ...
         'rU_L','rU_C','rU_R'});
    writetable(Tlog, fullfile(outdir,'diag_fast_target_log.csv'));

    fprintf('\nDone. Saved diagnostics to:\n%s\n', outdir);
end


% ========================= helpers =========================

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

function [sigC, info] = choose_sigma(sig_grid, y, par)
    info = struct('kind','','idx',NaN);

    above = (y >= par.rthr);

    switch par.crossing_direction
        case 'down'
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
            cross_all = find( ~above(1:end-1) & above(2:end) );

            if isempty(cross_all)
                [~, idx] = min(abs(y - par.rthr));
                sigC = sig_grid(idx);
                info.kind = 'closest_to_threshold';
                info.idx  = idx;
            else
                cross_idx = cross_all(1);
                sigC = 0.5*(sig_grid(cross_idx) + sig_grid(cross_idx+1));
                info.kind = 'first_upward_crossing';
                info.idx  = cross_idx;
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
