%% horizontal_sweep_target.m
% Small horizontal sweep for HALF-TARGET feature debugging.
%
% Uses feature_evaluation_target.m directly.
% Purpose:
%   - locate possible transition sigma values
%   - check whether the target feature matches spacetime plots
%
% Output per b:
%   - sweep_log.csv
%   - target_feature_vs_sigma.png
%   - spacetime pngs per sigma

clear; clc; close all;

%% ---------------- USER SETTINGS ----------------
b_list = [9.498];   % try one or several b values

sigma_start = 0.419;
sigma_end   = 0.425;
ds_sigma    = 0.001;          % start coarse, then refine later

T_relax = 240;
T_obs   = 240;

do_spacetime = true;

seed_file = fullfile('ic_half_target_0.45_10.00.mat');
if ~isfile(seed_file)
    seed_file = fullfile('half_target_seed.mat');
end

timestamp = datestr(now,'yyyymmdd_HHMMSS');
rootdir = sprintf('horizSweep_target_%s', timestamp);
if ~exist(rootdir,'dir'); mkdir(rootdir); end

%% ---------------- LOAD SEED ----------------
S = load(seed_file);
if isfield(S,'ic')
    ic_seed = S.ic;
elseif isfield(S,'ic_end')
    ic_seed = S.ic_end;
else
    fn = fieldnames(S);
    ic_seed = S.(fn{1});
end

%% ---------------- FEATURE HANDLE ----------------
featpar = struct();
featpar.feature = 'target';
featpar.N = 1;

sigmas = sigma_start:ds_sigma:sigma_end;
n = numel(sigmas);

fprintf('Running target sweeps: %d b-values, %d sigmas each\n', numel(b_list), n);

for bb = 1:numel(b_list)
    b_fixed = b_list(bb);

    outdir = fullfile(rootdir, sprintf('b_%.4f', b_fixed));
    if ~exist(outdir,'dir'); mkdir(outdir); end

    log_mat = nan(n, 2); % sigma, feature

    fprintf('\n=== b = %.6f ===\n', b_fixed);

    for i = 1:n
        sigma = sigmas(i);

        fprintf('[b=%.4f] %d/%d sigma=%.6f\n', b_fixed, i, n, sigma);

        stepdir = fullfile(outdir, sprintf('sigma_%.6f', sigma));
        if do_spacetime && ~exist(stepdir,'dir'); mkdir(stepdir); end

        modelpar = struct();
        modelpar.model   = 'Brusselator_1D';
        modelpar.ic0     = ic_seed;
        modelpar.a       = sigma;
        modelpar.b       = b_fixed;
        modelpar.T_relax = T_relax;
        modelpar.T_obs   = T_obs;

        % feature
        feat_val = feature_evaluation_target(featpar, modelpar);
        log_mat(i,:) = [sigma, feat_val];

        % spacetime diagnostics
        if do_spacetime
            par = struct('sigma', sigma, 'b', b_fixed);
            [ic_rel, ~, ~] = solve_brusselator_1d(ic_seed, par, T_relax, 0);
            [~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

            save_spacetime(stepdir, sprintf('target_sigma%.6f', sigma), V, feat_val);
        end
    end

    % save csv
    hdr = {'sigma','target_feature'};
    writecell(hdr, fullfile(outdir,'sweep_log.csv'));
    writematrix(log_mat, fullfile(outdir,'sweep_log.csv'), 'WriteMode','append');

    % plot feature vs sigma
    fig = figure('Visible','off');
    plot(log_mat(:,1), log_mat(:,2), 'o-'); grid on;
    xlabel('\sigma');
    ylabel('target feature');
    title(sprintf('b=%.4f | target feature vs sigma', b_fixed));
    exportgraphics(fig, fullfile(outdir,'target_feature_vs_sigma.png'));
    close(fig);

    fprintf('Saved: %s\n', outdir);
end

fprintf('\nAll done. Root output folder: %s\n', rootdir);

%% ---------------- Helpers ----------------
function save_spacetime(stepdir, tag, V, feat_val)
    fig = figure('Visible','off');
    imagesc(V); axis tight; set(gca,'YDir','normal');
    colormap(gca,'sky'); colorbar;
    xlabel('space'); ylabel('time');
    title(sprintf('%s | feature=%.3e', tag, feat_val), 'Interpreter','tex');
    exportgraphics(fig, fullfile(stepdir, sprintf('spacetime_%s.png', tag)));
    close(fig);
end