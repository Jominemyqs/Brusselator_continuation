%% Main_spiral.m
% Runs continuation along the spiral/stripe transition curve in both directions.
% Designed for Slurm array jobs via environment variable TASK.
%
% TASK convention:
%   0 = RIGHT-UP
%   1 = RIGHT-DOWN
%   2 = LEFT-UP
%   3 = LEFT-DOWN
%   4 = MIDDLE-RIGHT
%   5 = MIDDLE-LEFT
clear; clc; close all;
rng(5);

project_root = fileparts(mfilename('fullpath'));

task = str2double(getenv("TASK"));
if isnan(task); task = 0; end

% ---- Define the six branches ----
sigma_mid = 0.4200;
b_mid     = 9.4980;
ds_mid    = 0.0010;

switch task
    case 0  % RIGHT-UP
        p0 = [0.51025, 10.00];
        p1 = [0.50950, 10.02];

    case 1  % RIGHT-DOWN
        p0 = [0.51025, 10.00];
        p1 = [0.51175,  9.98];

    case 2  % LEFT-UP
        p0 = [0.3204, 10.0000];
        p1 = [0.3202, 10.0200];

    case 3  % LEFT-DOWN
        p0 = [0.3204, 10.0000];
        p1 = [0.32055,  9.9800];

    case 4  % MIDDLE-RIGHT
        % Direction seed only. The first accepted point is obtained by the
        % vertical b-corrector and the approximate seed is not reported.
        p0 = [sigma_mid,          b_mid];
        p1 = [sigma_mid + ds_mid, b_mid];

    case 5  % MIDDLE-LEFT
        p0 = [sigma_mid,          b_mid];
        p1 = [sigma_mid - ds_mid, b_mid];

    otherwise
        error("Unknown TASK=%d", task);
end

% secant direction
d = (p1 - p0); d = d(:) / norm(d);
p_start = p1(:);

run_up   = ismember(task, [0 2]);
run_dn   = ismember(task, [1 3]);
run_mid  = ismember(task, [4 5]);

is_right_branch = ismember(task, [0 1 4]);
is_left_branch  = ismember(task, [2 3 5]);

% ---------------- FEATURE HANDLE ----------------
if is_right_branch
    feature_handle = @feature_evaluation_right;
elseif is_left_branch
    feature_handle = @feature_evaluation_left;
else
    error('Could not determine branch side from TASK=%d', task);
end

obj_handle = @objective_evaluation;

%% ---------------- MODEL PARAMETERS ----------------
modelpar           = struct;
modelpar.model     = 'Brusselator_1D';
modelpar.init      = 'random';
modelpar.sets      = 'solution';
modelpar.xdim      = 1;

modelpar.init_db       = 0.02;
modelpar.init_da       = 0.02;
modelpar.branch        = 'spiral';
modelpar.save_IC       = true;
modelpar.perturb_IC    = false;
modelpar.find_max_time = false;

% Relax / observation times used inside feature_evaluation
modelpar.T_relax   = 360;
modelpar.T_obs     = 360;
modelpar.last_frac = 0.25;

% ---- Branch-side detector settings ----
if is_right_branch
    modelpar.r_thr              = 2.6e-3;
    modelpar.crossing_direction = 'down';   % high -> low
    modelpar.detector_side      = 'right';
elseif is_left_branch
    modelpar.r_thr              = 0.5;
    modelpar.crossing_direction = 'up';     % low -> high
    modelpar.detector_side      = 'left';
end

% --- Local sweep parameters ---
modelpar.dsigma_local          = 1e-3;
modelpar.local_halfwidth_steps = 10;
modelpar.expand_factors        = [1 2 4];

% Hysteresis / smoothing
modelpar.use_hysteresis   = true;
modelpar.drop_smooth_alpha = 0.35;
modelpar.track_window      = 8 * modelpar.dsigma_local;

% fixed threshold crossing strength
modelpar.strength_floor = 2e-4;

% hard cap
modelpar.eval_cap = 80;

% Middle branches run primarily left/right in sigma, so their corrector
% fixes sigma and searches vertically in b. These values recover the
% settings stored in the genuine middle-left continuation output.
modelpar.corrector_mode = 'horizontal_side';
if run_mid
    modelpar.corrector_mode          = 'vertical_middle';
    modelpar.db_local                = 0.01;
    modelpar.local_halfwidth_steps_b = 6;
    modelpar.expand_factors_b        = [1 2 4];
    modelpar.track_window_b          = 8 * modelpar.db_local;
    modelpar.eval_cap                = 120;
    modelpar.invalid_feature_value   = -1;

    if is_right_branch
        % The desired middle-right transition is the downward crossing
        % closest to the preceding b-value. A second, stronger upward
        % crossing can occur roughly 0.1 lower in b.
        modelpar.crossing_direction_b = 'down';
        modelpar.use_hysteresis       = true;
    else
        % Restored from the genuine historical middle-left metadata.
        modelpar.crossing_direction_b = 'up';
        modelpar.use_hysteresis       = false;
    end
end

% ---- Seed IC ----
S = load('ic_0.45_10.mat');
if isfield(S,'ic')
    ic_seed = S.ic;
else
    fn = fieldnames(S);
    ic_seed = S.(fn{1});
end

% Warm-start base for continuation
modelpar.ic_base = ic_seed;

% Fixed detection IC for local sweeps (NO sigma warm-start)
modelpar.ic_det  = ic_seed;

%% ---------------- FEATURE PARAMETERS ----------------
featpar              = struct;
featpar.feature      = 'spiral';
featpar.alpha        = 0.01;
featpar.avoid_steady = 0;
featpar.N            = 1;

%% ---------------- CONTINUATION PARAMETERS ----------------
contpar = struct;

contpar.step_size           = 0.01;
contpar.max_step_arc_length = 0.01;
contpar.min_step_arc_length = 0.01;

contpar.max_sim   = 3;
contpar.stopping  = -Inf;
contpar.b_bounds  = [8.0, 14.5];
contpar.sigma_bounds = [0.345879, Inf];
contpar.max_steps = 600;
contpar.dir_smooth_alpha = 0.15;
contpar.L = inf;

if run_mid
    contpar.step_size           = 0.003;
    contpar.max_step_arc_length = 0.003;
    contpar.min_step_arc_length = 0.003;
    contpar.sigma_bounds        = [0.28, 0.58];
    contpar.include_start       = false;
end

%% ---------------- OUTPUT DIR ----------------
tag = datestr(now,'yyyymmdd_HHMMSS');
output_root = fullfile(project_root, 'outputs', 'continuation');
if ~exist(output_root,'dir'); mkdir(output_root); end
outdir = fullfile(output_root, ...
    sprintf('full_continuation_spiral_task%d_%s', task, tag));
if ~exist(outdir,'dir'); mkdir(outdir); end

%% ---------------- RUN UPWARD CONTINUATION ----------------
if run_up
    start_up = struct;
    start_up.point  = p_start;
    start_up.normal = d;

    fprintf('\n=== Continuation UP starting at (%.6f, %.6f) ===\n', ...
        start_up.point(1), start_up.point(2));
    fprintf('Detector side: %s | r_thr = %.6g | crossing_direction = %s\n', ...
        modelpar.detector_side, modelpar.r_thr, modelpar.crossing_direction);

    [p_hist_up, counts_up, L_hist_up, p_all_up, metric_all_up] = continuation( ...
        contpar, featpar, modelpar, feature_handle, obj_handle, start_up, 1);

    save(fullfile(outdir,'right_up.mat'), ...
        'p_hist_up','counts_up','L_hist_up','p_all_up','metric_all_up', ...
        'start_up','contpar','featpar','modelpar','p0','p1','d');
end

%% ---------------- RUN DOWNWARD CONTINUATION ----------------
if run_dn
    start_dn = struct;
    start_dn.point  = p_start;
    start_dn.normal = d;

    fprintf('\n=== Continuation DOWN starting at (%.6f, %.6f) ===\n', ...
        start_dn.point(1), start_dn.point(2));
    fprintf('Detector side: %s | r_thr = %.6g | crossing_direction = %s\n', ...
        modelpar.detector_side, modelpar.r_thr, modelpar.crossing_direction);

    [p_hist_dn, counts_dn, L_hist_dn, p_all_dn, metric_all_dn] = continuation( ...
        contpar, featpar, modelpar, feature_handle, obj_handle, start_dn, 1);

    save(fullfile(outdir,'right_down.mat'), ...
        'p_hist_dn','counts_dn','L_hist_dn','p_all_dn','metric_all_dn', ...
        'start_dn','contpar','featpar','modelpar','p0','p1','d');
end

%% ---------------- RUN MIDDLE RE-SEEDED CONTINUATION ----------------
if run_mid
    start_mid = struct;
    start_mid.point  = p_start;
    start_mid.normal = d;

    fprintf('\n=== Continuation MIDDLE starting at (%.6f, %.6f) ===\n', ...
        start_mid.point(1), start_mid.point(2));
    fprintf('Detector side: %s | r_thr = %.6g | crossing_direction = %s\n', ...
        modelpar.detector_side, modelpar.r_thr, modelpar.crossing_direction);

    [p_hist_mid, counts_mid, L_hist_mid, p_all_mid, metric_all_mid] = continuation( ...
        contpar, featpar, modelpar, feature_handle, obj_handle, start_mid, 1);

    switch task
        case 4
            save_name = 'middle_right.mat';
        case 5
            save_name = 'middle_left.mat';
        otherwise
            save_name = 'middle_branch.mat';
    end

    save(fullfile(outdir, save_name), ...
        'p_hist_mid','counts_mid','L_hist_mid','p_all_mid','metric_all_mid', ...
        'start_mid','contpar','featpar','modelpar','p0','p1','d', ...
        'sigma_mid','b_mid','ds_mid');
end

fprintf('\nDone. Saved to: %s\n', outdir);
