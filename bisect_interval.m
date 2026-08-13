function [p_current, counter, sim_p_history, sim_metric_history, val_left, val_right, aux] = ...
    bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count, aux_in) %#ok<INUSD>
%BISECT_INTERVAL  Local threshold-crossing corrector with adaptive widening.
%
% Two correction geometries are supported:
%   horizontal_side (default): fix b and sweep sigma.
%   vertical_middle:           fix sigma and sweep b.
%
% Every evaluation in one corrector sweep uses the same detection seed.
% The routine selects a crossing with the requested orientation, optionally
% tracks the crossing closest to the previous one, and returns the midpoint
% of the bracketing samples. The historical function name is retained so
% existing drivers continue to work.

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

if ~strcmp(modelpar.model,'Brusselator_1D')
    error('bisect_interval: specialized for Brusselator_1D.');
end
if ~(strcmp(featpar.feature,'spiral') || strcmp(featpar.feature,'target'))
    error('bisect_interval: featpar.feature must be ''spiral'' or ''target''.');
end
if ~isfield(modelpar,'ic_base') || isempty(modelpar.ic_base)
    error('bisect_interval: modelpar.ic_base missing.');
end

% p_c1 and p_c2 are symmetric around the predictor.
p_pred = 0.5*(p_c1 + p_c2);
sigma_pred = p_pred(1);
b_pred = p_pred(2);

corrector_mode = 'horizontal_side';
if isfield(modelpar,'corrector_mode') && ~isempty(modelpar.corrector_mode)
    corrector_mode = char(modelpar.corrector_mode);
end

switch lower(corrector_mode)
    case {'horizontal_side','horizontal','sigma'}
        sweep_axis = 'sigma';
        sweep_center = sigma_pred;
        fixed_value = b_pred;

        coordinate_step = 2.5e-4;
        if isfield(modelpar,'dsigma_local') && ~isempty(modelpar.dsigma_local)
            coordinate_step = modelpar.dsigma_local;
        end

        base_hw = 14;
        if isfield(modelpar,'local_halfwidth_steps') && ~isempty(modelpar.local_halfwidth_steps)
            base_hw = modelpar.local_halfwidth_steps;
        end

        expand_factors = [1 2 4];
        if isfield(modelpar,'expand_factors') && ~isempty(modelpar.expand_factors)
            expand_factors = modelpar.expand_factors(:).';
        end

        crossing_direction = get_char_field(modelpar, 'crossing_direction', 'down');
        track_window = get_numeric_field(modelpar, 'track_window', 8*coordinate_step);

    case {'vertical_middle','vertical','b'}
        sweep_axis = 'b';
        sweep_center = b_pred;
        fixed_value = sigma_pred;

        coordinate_step = 0.01;
        if isfield(modelpar,'db_local') && ~isempty(modelpar.db_local)
            coordinate_step = modelpar.db_local;
        end

        base_hw = 6;
        if isfield(modelpar,'local_halfwidth_steps_b') && ~isempty(modelpar.local_halfwidth_steps_b)
            base_hw = modelpar.local_halfwidth_steps_b;
        end

        expand_factors = [1 2 4];
        if isfield(modelpar,'expand_factors_b') && ~isempty(modelpar.expand_factors_b)
            expand_factors = modelpar.expand_factors_b(:).';
        elseif isfield(modelpar,'expand_factors') && ~isempty(modelpar.expand_factors)
            expand_factors = modelpar.expand_factors(:).';
        end

        crossing_direction = get_char_field(modelpar, 'crossing_direction_b', ...
            get_char_field(modelpar, 'crossing_direction', 'up'));
        track_window = get_numeric_field(modelpar, 'track_window_b', 8*coordinate_step);

    otherwise
        error('bisect_interval: unknown corrector_mode = %s', corrector_mode);
end

if ~isscalar(coordinate_step) || ~isfinite(coordinate_step) || coordinate_step <= 0
    error('bisect_interval: local sweep spacing must be a positive finite scalar.');
end

r_thr = get_numeric_field(modelpar, 'r_thr', 2e-3);
use_hysteresis = logical(get_numeric_field(modelpar, 'use_hysteresis', true));
alpha = get_numeric_field(modelpar, 'drop_smooth_alpha', 0.35);
strength_floor = get_numeric_field(modelpar, 'strength_floor', 2e-4);

eval_cap = 80;
if isfield(modelpar,'eval_cap') && ~isempty(modelpar.eval_cap)
    eval_cap = modelpar.eval_cap;
elseif ~isempty(max_count) && max_count > 0
    eval_cap = min(120, max(60, 20 * max_count));
end
eval_cap = max(10, round(eval_cap));

% Fixed IC for every point in this local sweep: no path dependence within
% the corrector grid.
if isfield(modelpar,'ic_det') && ~isempty(modelpar.ic_det)
    ic0_fixed = modelpar.ic_det;
else
    ic0_fixed = modelpar.ic_base;
end

aux.last_corrector_mode = corrector_mode;
aux.last_sweep_axis = sweep_axis;
aux.last_fixed_value = fixed_value;
aux.last_sigma_pred = sigma_pred;
aux.last_b_pred = b_pred;
aux.last_sigma_window = [NaN, NaN];
aux.last_b_window = [NaN, NaN];
aux.last_coordinate_window = [NaN, NaN];
aux.last_rLR = [NaN, NaN];
aux.last_r_thr = r_thr;
aux.last_hw_used = NaN;
aux.last_eval_cap = eval_cap;
aux.last_counter = 0;

best_found = false;
best = struct('coordinate',[],'strength',[],'hw',[],'mode',"");

have_prev = false;
previous_crossing = NaN;
if isfield(aux_in,'prev_crossing') && ~isempty(aux_in.prev_crossing)
    previous_crossing = aux_in.prev_crossing;
    have_prev = true;
elseif isfield(aux_in,'prev_drop') && ~isempty(aux_in.prev_drop)
    % Backward compatibility with previously saved continuation state.
    previous_crossing = aux_in.prev_drop;
    have_prev = true;
end

for ff = 1:numel(expand_factors)
    hw = max(1, round(base_hw * expand_factors(ff)));
    coordinates = sweep_center + (-hw:hw)*coordinate_step;
    M = numel(coordinates);

    rvals = zeros(1,M);
    k_last = 0;

    for k = 1:M
        if counter >= eval_cap
            break;
        end

        if strcmp(sweep_axis,'sigma')
            sigma = coordinates(k);
            b = fixed_value;
        else
            sigma = fixed_value;
            b = coordinates(k);
        end

        rvals(k) = eval_r_from_ic(ic0_fixed, modelpar, featpar, feature_handle, sigma, b);
        counter = counter + 1;
        k_last = k;
    end

    coordinates = coordinates(1:k_last);
    rvals = rvals(1:k_last);
    M = k_last;

    if M >= 1
        if strcmp(sweep_axis,'sigma')
            sim_p_history = [coordinates; fixed_value*ones(1,M)];
            aux.last_sigma_window = [coordinates(1), coordinates(end)];
        else
            sim_p_history = [fixed_value*ones(1,M); coordinates];
            aux.last_b_window = [coordinates(1), coordinates(end)];
        end
        sim_metric_history = rvals;
        val_left = rvals(1);
        val_right = rvals(end);

        aux.last_coordinate_window = [coordinates(1), coordinates(end)];
        aux.last_rLR = [val_left, val_right];
        aux.last_hw_used = hw;
        aux.last_counter = counter;
    end

    if M < 2
        continue;
    end

    g = rvals - r_thr;
    cross_all = find(g(1:end-1).*g(2:end) <= 0);
    if isempty(cross_all)
        continue;
    end

    switch lower(crossing_direction)
        case 'down'
            keep = g(cross_all) > 0 & g(cross_all+1) <= 0;
        case 'up'
            keep = g(cross_all) < 0 & g(cross_all+1) >= 0;
        otherwise
            error('bisect_interval: unknown crossing direction = %s', crossing_direction);
    end

    cross = cross_all(keep);
    if isempty(cross)
        continue;
    end

    dr = diff(rvals);
    mids = 0.5*(coordinates(cross) + coordinates(cross+1));
    strength = abs(dr(cross));

    if use_hysteresis && have_prev
        near = find(abs(mids - previous_crossing) <= track_window);
        if ~isempty(near)
            [smax, jj] = max(strength(near));
            best_found = true;
            best.coordinate = mids(near(jj));
            best.strength = smax;
            best.hw = hw;
            best.mode = "hysteresis_near_previous_threshold_crossing";

            if smax >= strength_floor
                break;
            end
        end
    end

    [smax_all, jj_all] = max(strength);
    coordinate_candidate = mids(jj_all);
    if ~best_found || (smax_all > best.strength)
        best_found = true;
        best.coordinate = coordinate_candidate;
        best.strength = smax_all;
        best.hw = hw;
        if use_hysteresis && have_prev
            best.mode = "fallback_strongest_threshold_crossing";
        else
            best.mode = "strongest_threshold_crossing";
        end
    end

    if counter >= eval_cap
        break;
    end
end

if ~best_found
    p_current = p_pred;
    aux.prev_crossing = sweep_center;
    aux.prev_drop = sweep_center;
    aux.used_mode = "no_crossing_stick_to_predictor_" + string(sweep_axis);
    aux.no_cross_count = aux.no_cross_count + 1;
    return
end

aux.no_cross_count = 0;
coordinate_star = best.coordinate;
aux.used_mode = best.mode + "_" + string(sweep_axis);

if have_prev
    coordinate_star = (1-alpha)*previous_crossing + alpha*coordinate_star;
    aux.used_mode = aux.used_mode + "_smoothed";
end

if strcmp(sweep_axis,'sigma')
    p_current = [coordinate_star; fixed_value];
else
    p_current = [fixed_value; coordinate_star];
end
aux.prev_crossing = coordinate_star;
aux.prev_drop = coordinate_star;

end


function r = eval_r_from_ic(ic0, modelpar, featpar, feature_handle, sigma, b)
mp = modelpar;
mp.a = sigma;
mp.b = b;
mp.ic0 = ic0;
r = feature_handle(featpar, mp);
if iscell(r), r = r{1}; end
if ~isscalar(r) || ~isfinite(r)
    error('bisect_interval: feature evaluation must return one finite scalar.');
end
end


function out = get_numeric_field(S, name, default_value)
out = default_value;
if isfield(S,name) && ~isempty(S.(name))
    out = S.(name);
end
end


function out = get_char_field(S, name, default_value)
out = default_value;
if isfield(S,name) && ~isempty(S.(name))
    out = char(S.(name));
end
end
