function [p_current, counter, sim_p_history, sim_metric_history, val_left, val_right, aux] = ...
    bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count, aux_in)
%BISECT_INTERVAL  Threshold-crossing corrector with adaptive window widening (NO sigma warm-start).
%
% Right side:
%   - uses feature_evaluation_right
%   - downward crossing of r = r_thr
%
% Left side:
%   - uses feature_evaluation_left
%   - upward crossing of r = r_thr
%
% Hysteresis:
%   choose crossing nearest prev_drop, tie-break by largest jump strength.

if nargin < 9 || isempty(aux_in), aux_in = struct(); end
aux = aux_in;

counter = 0;
sim_p_history = [];
sim_metric_history = [];
val_left = NaN;
val_right = NaN;

if ~isfield(aux,'no_cross_count') || isempty(aux.no_cross_count)
    aux.no_cross_count = 0;
end

% sanity
if ~strcmp(modelpar.model,'Brusselator_1D')
    error('bisect_interval: specialized for Brusselator_1D.');
end

if ~(strcmp(featpar.feature,'spiral') || strcmp(featpar.feature,'target'))
    error('bisect_interval: featpar.feature must be ''spiral'' or ''target''.');
end

if ~isfield(modelpar,'ic_base') || isempty(modelpar.ic_base)
    error('bisect_interval: modelpar.ic_base missing.');
end

% predictor midpoint
p_pred = 0.5*(p_c1 + p_c2);
sigma_pred = p_pred(1);
b_fix      = p_pred(2);

% sweep resolution
ds = 2.5e-4;
if isfield(modelpar,'dsigma_local') && ~isempty(modelpar.dsigma_local)
    ds = modelpar.dsigma_local;
end

base_hw = 14;
if isfield(modelpar,'local_halfwidth_steps') && ~isempty(modelpar.local_halfwidth_steps)
    base_hw = modelpar.local_halfwidth_steps;
end

% threshold
r_thr = 2e-3;
if isfield(modelpar,'r_thr') && ~isempty(modelpar.r_thr)
    r_thr = modelpar.r_thr;
end

% crossing direction
crossing_direction = 'down';
if isfield(modelpar,'crossing_direction') && ~isempty(modelpar.crossing_direction)
    crossing_direction = modelpar.crossing_direction;
end

% hysteresis window
track_window = 8*ds;
if isfield(modelpar,'track_window') && ~isempty(modelpar.track_window)
    track_window = modelpar.track_window;
end

use_hysteresis = true;
if isfield(modelpar,'use_hysteresis')
    use_hysteresis = logical(modelpar.use_hysteresis);
end

alpha = 0.35;
if isfield(modelpar,'drop_smooth_alpha') && ~isempty(modelpar.drop_smooth_alpha)
    alpha = modelpar.drop_smooth_alpha;
end

% widening plan
expand_factors = [1 2 4];
if isfield(modelpar,'expand_factors') && ~isempty(modelpar.expand_factors)
    expand_factors = modelpar.expand_factors(:).';
end

% hard cap
eval_cap = 80;
if isfield(modelpar,'eval_cap') && ~isempty(modelpar.eval_cap)
    eval_cap = modelpar.eval_cap;
else
    if ~isempty(max_count) && max_count > 0
        eval_cap = min(120, max(60, 20 * max_count));
    end
end
eval_cap = max(10, round(eval_cap));

% early quit setting
strength_floor = 2e-4;
if isfield(modelpar,'strength_floor') && ~isempty(modelpar.strength_floor)
    strength_floor = modelpar.strength_floor;
end

% fixed IC for ALL sigma evaluations in this local sweep
if isfield(modelpar,'ic_det') && ~isempty(modelpar.ic_det)
    ic0_fixed = modelpar.ic_det;
else
    ic0_fixed = modelpar.ic_base;
end

% initialize debug fields
aux.last_b_fix        = b_fix;
aux.last_sigma_pred   = sigma_pred;
aux.last_sigma_window = [NaN, NaN];
aux.last_rLR          = [NaN, NaN];
aux.last_r_thr        = r_thr;
aux.last_hw_used      = NaN;
aux.last_eval_cap     = eval_cap;
aux.last_counter      = 0;

best_found = false;
best = struct('mids',[],'strength',[],'hw',[],'mode',"");

have_prev = isfield(aux_in,'prev_drop') && ~isempty(aux_in.prev_drop);

for ff = 1:numel(expand_factors)
    hw = max(1, round(base_hw * expand_factors(ff)));
    sigmas = sigma_pred + (-hw:hw)*ds;
    M = numel(sigmas);

    rvals = zeros(1,M);
    k_last = 0;

    for k = 1:M
        if counter >= eval_cap
            break;
        end

        rvals(k) = eval_r_from_ic(ic0_fixed, modelpar, featpar, feature_handle, sigmas(k), b_fix);
        counter = counter + 1;
        k_last = k;
    end

    sigmas = sigmas(1:k_last);
    rvals  = rvals(1:k_last);
    M      = k_last;

    if M >= 1
        sim_p_history      = [sigmas; b_fix*ones(1,M)];
        sim_metric_history = rvals;
        val_left  = rvals(1);
        val_right = rvals(end);

        aux.last_sigma_window = [sigmas(1), sigmas(end)];
        aux.last_rLR          = [val_left, val_right];
        aux.last_hw_used      = hw;
        aux.last_counter      = counter;
    end

    if M < 2
        continue;
    end

    % threshold crossing candidates
    g = rvals - r_thr;
    cross_all = find(g(1:end-1).*g(2:end) <= 0);

    if isempty(cross_all)
        continue;
    end

    % keep only crossings with the requested direction
    switch crossing_direction
        case 'down'
            keep = g(cross_all) > 0 & g(cross_all+1) <= 0;
        case 'up'
            keep = g(cross_all) < 0 & g(cross_all+1) >= 0;
        otherwise
            error('bisect_interval: unknown crossing_direction = %s', crossing_direction);
    end

    cross = cross_all(keep);

    if isempty(cross)
        continue;
    end

    dr = diff(rvals);
    mids = 0.5*(sigmas(cross) + sigmas(cross+1));
    strength = abs(dr(cross));

    if use_hysteresis && have_prev
        prev = aux_in.prev_drop;
        near = find(abs(mids - prev) <= track_window);

        if ~isempty(near)
            [smax, jj] = max(strength(near));
            sigma_candidate = mids(near(jj));

            best_found = true;
            best.mids = sigma_candidate;
            best.strength = smax;
            best.hw = hw;
            best.mode = "hysteresis_near_prev_threshold_crossing";

            if smax >= strength_floor
                break;
            end
        end
    end

    [smax_all, jj_all] = max(strength);
    sigma_candidate = mids(jj_all);

    if ~best_found || (smax_all > best.strength)
        best_found = true;
        best.mids = sigma_candidate;
        best.strength = smax_all;
        best.hw = hw;

        if use_hysteresis && have_prev
            best.mode = "fallback_strongest_threshold_crossing";
        else
            best.mode = "init_strongest_threshold_crossing";
        end
    end

    if counter >= eval_cap
        break;
    end
end

% ---------------- output ----------------
if ~best_found
    sigma_star = sigma_pred;
    p_current  = [sigma_star; b_fix];

    aux.prev_drop = sigma_star;
    aux.used_mode = "no_crossing_stick_pred";
    aux.no_cross_count = aux.no_cross_count + 1;
    return
end

aux.no_cross_count = 0;

sigma_star = best.mids;
aux.used_mode = best.mode;

if have_prev
    sigma_star = (1-alpha)*aux_in.prev_drop + alpha*sigma_star;
    aux.used_mode = aux.used_mode + "_smoothed";
end

p_current = [sigma_star; b_fix];
aux.prev_drop = sigma_star;

end

function r = eval_r_from_ic(ic0, modelpar, featpar, feature_handle, sigma, b)
mp = modelpar;
mp.a = sigma;
mp.b = b;
mp.ic0 = ic0;
r = feature_handle(featpar, mp);
if iscell(r); r = r{1}; end
end

