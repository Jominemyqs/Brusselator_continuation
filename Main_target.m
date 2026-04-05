%% Main_target.m
% Runs continuation along the TARGET-pattern transition curves in both directions.
% Designed for Slurm array jobs via environment variable TASK.
%
% TASK:
%   0 = RIGHT-UP
%   1 = RIGHT-DOWN
%   2 = LEFT-UP
%   3 = LEFT-DOWN
%   4 = MIDDLE-RIGHT
%   5 = MIDDLE-LEFT
% Requires:
%   continuation.m
%   bisect_interval.m
%   objective_evaluation.m
%   feature_evaluation_target.m
%   solve_brusselator_1d.m
%
% Seed file expected:
%   out_half_target/ic_half_target_0.45_10.00.mat
% or fallback:
%   out_half_target/half_target_seed.mat

clear; clc; close all;
rng(5);

addpath(genpath(pwd));

task = str2double(getenv("TASK"));
if isnan(task); task = 0; end

%% ---------------- DEFINE FOUR TARGET BRANCHES ----------------
% Based on your refined horizontal sweeps:
%   left transition  ~ sigma = 0.326
%   right transition ~ sigma = 0.511
%
% Since the transition location is nearly unchanged between b=9.98,10,10.02,
% the initial tangent is essentially vertical.

%% ---------------- DEFINE SIX TARGET BRANCHES ----------------
sigma_mid = 0.4200;
b_mid     = 9.4980;
ds_mid    = 0.0010;

switch task
    case 0  % RIGHT-UP
        p0 = [0.51100, 10.0000];
        p1 = [0.51125, 10.0200];

    case 1  % RIGHT-DOWN
        p0 = [0.51100, 10.0000];
        p1 = [0.51100,  9.9800];

    case 2  % LEFT-UP
        p0 = [0.32600, 10.0000];
        p1 = [0.32600, 10.0200];

    case 3  % LEFT-DOWN
        p0 = [0.32600, 10.0000];
        p1 = [0.32600,  9.9800];

    case 4  % MIDDLE-RIGHT
        p0 = [0.42, 9.498];
        p1 = [0.419, 9.498];

    case 5  % MIDDLE-LEFT
        p0 = [0.42, 9.498];
        p1 = [0.419, 9.498];

    otherwise
        error("Unknown TASK=%d", task);
end

% Compute tangent direction from secant
d = (p1 - p0);
d = d(:) / norm(d);
p_start = p1(:);

run_up   = ismember(task, [0 2]);
run_dn   = ismember(task, [1 3]);
run_mid  = ismember(task, [4 5]);

%% ---------------- FEATURE / OBJECTIVE ----------------
feature_handle = @feature_evaluation_target;
obj_handle     = @objective_evaluation;

%% ---------------- MODEL PARAMETERS ----------------
modelpar           = struct;
modelpar.model     = 'Brusselator_1D';
modelpar.init      = 'random';
modelpar.sets      = 'solution';
modelpar.xdim      = 1;

modelpar.init_db       = 0.02;
modelpar.init_da       = 0.02;
modelpar.branch        = 'target';
modelpar.save_IC       = true;
modelpar.perturb_IC    = false;
modelpar.find_max_time = false;

% Times used by target feature
modelpar.T_relax = 240;
modelpar.T_obs   = 240;

% Threshold from your horizontal sweeps
modelpar.r_thr = 8;

% Local sweep parameters
modelpar.dsigma_local          = 1e-3;
modelpar.local_halfwidth_steps = 10;
modelpar.expand_factors        = [1 2 4];

% Hysteresis / smoothing
modelpar.track_window      = 8 * modelpar.dsigma_local;
modelpar.use_hysteresis    = true;
modelpar.drop_smooth_alpha = 0.35;
modelpar.eval_cap          = 80;
modelpar.strength_floor    = 2e-4;

% Crossing direction depends on branch side
if ismember(task, [0 1 4])   % right + middle-right
    modelpar.crossing_direction = 'down';
else                         % left + middle-left
    modelpar.crossing_direction = 'up';
end

%% ---------------- LOAD HALF-TARGET SEED ----------------
seed_file = 'half_target_seed.mat';
if ~isfile(seed_file)
    error('Could not find target seed file: %s', seed_file);
end

S = load(seed_file);
if isfield(S,'ic')
    ic_seed = S.ic;
elseif isfield(S,'ic_end')
    ic_seed = S.ic_end;
else
    fn = fieldnames(S);
    ic_seed = S.(fn{1});
end

modelpar.ic_base = ic_seed;
modelpar.ic_det  = ic_seed;



%% ---------------- FEATURE PARAMETERS ----------------
featpar              = struct;
featpar.feature      = 'target';
featpar.alpha        = 0.01;
featpar.avoid_steady = 0;
featpar.N            = 1;

%% ---------------- CONTINUATION PARAMETERS ----------------
contpar = struct;

% First-pass conservative settings
contpar.step_size           = 0.01;
contpar.max_step_arc_length = 0.01;
contpar.min_step_arc_length = 0.01;

contpar.max_sim   = 3;
contpar.stopping  = -Inf;

% Figure-13-like target range
contpar.b_bounds  = [7, 16];

contpar.max_steps = 600;
contpar.dir_smooth_alpha = 0.15;
contpar.L = inf;

%% ---------------- OUTPUT DIR ----------------
tag = datestr(now,'yyyymmdd_HHMMSS');
outdir = sprintf('full_continuation_target_task%d_%s', task, tag);
if ~exist(outdir,'dir'); mkdir(outdir); end

%% ---------------- RUN UPWARD CONTINUATION ----------------
if run_up
    start_up = struct;
    start_up.point  = p_start;
    start_up.normal = d;

    fprintf('\n=== TARGET continuation UP starting at (%.6f, %.6f) ===\n', ...
        start_up.point(1), start_up.point(2));

    [p_hist_up, counts_up, L_hist_up, p_all_up, metric_all_up] = continuation( ...
        contpar, featpar, modelpar, feature_handle, obj_handle, start_up, 1);

    switch task
        case 0
            save_name = 'target_right_up.mat';
        case 2
            save_name = 'target_left_up.mat';
        otherwise
            save_name = 'target_up.mat';
    end

    save(fullfile(outdir, save_name), ...
        'p_hist_up','counts_up','L_hist_up','p_all_up','metric_all_up', ...
        'start_up','contpar','featpar','modelpar','p0','p1','d','seed_file');
end

%% ---------------- RUN DOWNWARD CONTINUATION ----------------
if run_dn
    start_dn = struct;
    start_dn.point  = p_start;
    start_dn.normal = d;

    fprintf('\n=== TARGET continuation DOWN starting at (%.6f, %.6f) ===\n', ...
        start_dn.point(1), start_dn.point(2));

    [p_hist_dn, counts_dn, L_hist_dn, p_all_dn, metric_all_dn] = continuation( ...
        contpar, featpar, modelpar, feature_handle, obj_handle, start_dn, 1);

    switch task
        case 1
            save_name = 'target_right_down.mat';
        case 3
            save_name = 'target_left_down.mat';
        otherwise
            save_name = 'target_down.mat';
    end

    save(fullfile(outdir, save_name), ...
        'p_hist_dn','counts_dn','L_hist_dn','p_all_dn','metric_all_dn', ...
        'start_dn','contpar','featpar','modelpar','p0','p1','d','seed_file');
end

%% ---------------- RUN MIDDLE RE-SEEDED CONTINUATION ----------------
if run_mid
    start_mid = struct;
    start_mid.point  = p_start;
    start_mid.normal = d;

    fprintf('\n=== TARGET continuation MIDDLE starting at (%.6f, %.6f) ===\n', ...
        start_mid.point(1), start_mid.point(2));

    [p_hist_mid, counts_mid, L_hist_mid, p_all_mid, metric_all_mid] = continuation( ...
        contpar, featpar, modelpar, feature_handle, obj_handle, start_mid, 1);

    switch task
        case 4
            save_name = 'target_middle_right.mat';
        case 5
            save_name = 'target_middle_left.mat';
        otherwise
            save_name = 'target_middle.mat';
    end

    save(fullfile(outdir, save_name), ...
        'p_hist_mid','counts_mid','L_hist_mid','p_all_mid','metric_all_mid', ...
        'start_mid','contpar','featpar','modelpar','p0','p1','d', ...
        'seed_file','sigma_mid','b_mid','ds_mid');
end

fprintf('\nDone. Saved to: %s\n', outdir);

exit
