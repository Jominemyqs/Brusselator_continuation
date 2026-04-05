function [p_history, counts, L_history, p_history_all, metric_history_all] = continuation( ...
    contpar, featpar, modelpar, feature_handle, obj_handle, start, viz )
%CONTINUATION  Predictor-corrector continuation with warm-start and drop-tracking corrector.
%
% Requires bisect_interval signature:
%   [p_current, counter, sim_p_history, sim_metric_history, val_left, val_right, aux] = ...
%       bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count, aux_in);

if nargin < 7, viz = 1; end

% ------------------- Warm-start setup (Brusselator_1D + spiral feature) -------------------
use_warmstart = isfield(modelpar,'model') && strcmp(modelpar.model,'Brusselator_1D') && ...
                isfield(featpar,'feature') && strcmp(featpar.feature,'spiral');

if use_warmstart
    if ~isfield(modelpar,'ic_base') || isempty(modelpar.ic_base)
        S = load('ic_0.45_10.mat');
        if isfield(S,'ic'), modelpar.ic_base = S.ic;
        else, fn = fieldnames(S); modelpar.ic_base = S.(fn{1});
        end
    end
    if ~isfield(modelpar,'T_relax') || isempty(modelpar.T_relax)
        modelpar.T_relax = 360;
    end
end
% -----------------------------------------------------------------------------------------

% ------------------- Direction smoothing (kills zig-zag) ---------------------------------
dir_smooth_alpha = 0.25;
if isfield(contpar,'dir_smooth_alpha') && ~isempty(contpar.dir_smooth_alpha)
    dir_smooth_alpha = contpar.dir_smooth_alpha;
end
% -----------------------------------------------------------------------------------------

% Normalize initial tangent direction
start.normal = start.normal(:);
start.normal = start.normal / norm(start.normal);

p_current = start.point(:);
p_prior   = start.point(:);

% Store accepted points as columns (include start twice for legacy compatibility)
p_history = [p_prior, p_current];

direction      = start.normal / norm(start.normal);
direction_norm = [-direction(2); direction(1)];

L_current = 0;
i = 0;

max_count = contpar.max_sim;
min_step_arc_length = contpar.min_step_arc_length;
max_step_arc_length = contpar.max_step_arc_length;
step_size = contpar.step_size;

p_history_all      = [];
metric_history_all = [];
L_history          = [];

% aux state passed to bisect_interval
aux = struct();
aux.prev_drop = p_current(1);
aux.used_mode = "";
aux.no_cross_count = 0;

while L_current < contpar.L
    i = i + 1;
    fprintf('step %d\n', i);

    p_prior = p_current;

    % Predictor
    p_pred = p_prior + step_size * direction;

    % Bracket for corrector (endpoints along normal)
    p_c1 = p_pred + 0.5 * step_size * direction_norm;
    p_c2 = p_pred - 0.5 * step_size * direction_norm;
    dn   = 0.5 * step_size * direction_norm; %#ok<NASGU>

    % Inject warm-start IC into modelpar for feature evaluations
    if use_warmstart
        modelpar.ic0 = modelpar.ic_base; %#ok<NASGU>
    end

    % Corrector
    try
        [p_current, counts(i), p_history_new, metric_history_new, val_c1, val_c2, aux] = ...
            bisect_interval(modelpar, featpar, feature_handle, obj_handle, p_c1, p_c2, dn, max_count, aux);
    catch ME
        error('continuation: bisect_interval failed at step %d: %s', i, ME.message);
    end

    p_current = p_current(:);

    % ---------- ALWAYS print a compact step summary ----------
    used_mode = "";
    if isfield(aux,'used_mode'), used_mode = string(aux.used_mode); end

    % safe-get debug fields
    last_hw = NaN; last_ctr = NaN; last_win = [NaN NaN]; last_rLR = [NaN NaN]; last_thr = NaN;
    if isfield(aux,'last_hw_used'),      last_hw  = aux.last_hw_used; end
    if isfield(aux,'last_counter'),      last_ctr = aux.last_counter; end
    if isfield(aux,'last_sigma_window'), last_win = aux.last_sigma_window; end
    if isfield(aux,'last_rLR'),          last_rLR = aux.last_rLR; end
    if isfield(aux,'last_r_thr'),        last_thr = aux.last_r_thr; end

    fprintf('  pred=(%.6f, %.6f)  corr=(%.6f, %.6f)  dir=(%.4f, %.4f)\n', ...
        p_pred(1), p_pred(2), p_current(1), p_current(2), direction(1), direction(2));
    fprintf('  mode=%s | prev_drop=%.6f | sweep_hw=%g | evals=%g | rLR=[%.2e, %.2e] thr=%.2e win=[%.6f, %.6f]\n', ...
        used_mode, aux.prev_drop, last_hw, last_ctr, last_rLR(1), last_rLR(2), last_thr, last_win(1), last_win(2));
    % --------------------------------------------------------

    % STOP if repeated no-cross
    if isfield(aux,'no_cross_count') && aux.no_cross_count >= 3
        fprintf('STOP: no threshold crossing for %d consecutive steps.\n', aux.no_cross_count);
        break;
    end

    % Stop if outside target b-range
    if isfield(contpar,'b_bounds') && ~isempty(contpar.b_bounds)
        if p_current(2) < contpar.b_bounds(1) || p_current(2) > contpar.b_bounds(2)
            fprintf('STOP: b_bounds hit (b=%.6f outside [%.6f, %.6f])\n', ...
                p_current(2), contpar.b_bounds(1), contpar.b_bounds(2));
            break;
        end
    end

    % Append simulation histories
    if ~isempty(p_history_new)
        if size(p_history_new,1) ~= 2 && size(p_history_new,2) == 2
            p_history_new = p_history_new.';
        end
        if size(p_history_new,1) == 2
            p_history_all = [p_history_all, p_history_new];
        end
    end
    if ~isempty(metric_history_new)
        metric_history_new = metric_history_new(:).';
        metric_history_all = [metric_history_all, metric_history_new];
    end

    % Save accepted point
    p_history(:, end+1) = p_current; %#ok<AGROW>
    L_history(end+1)    = step_size; %#ok<AGROW>
    L_current           = L_current + step_size;

    % Warm-start update to accepted point
    if use_warmstart
        try
            par_tmp = struct('sigma', p_current(1), 'b', p_current(2));
            [ic_new, ~, ~] = solve_brusselator_1d(modelpar.ic_base, par_tmp, modelpar.T_relax, 0);
            modelpar.ic_base = ic_new;
        catch ME
            warning('Warm-start update failed at step %d (sigma=%.6f, b=%.6f): %s', ...
                i, p_current(1), p_current(2), ME.message);
        end
    end

    % Step size adaptation (your heuristic)
    if counts(i) <= 2
        step_size = max([ step_size / sqrt(sqrt(2)), min_step_arc_length ]);
        fprintf('step_size decreased to %.10g\n', step_size);
    end
    if counts(i) >= max_count
        step_size = min([ sqrt(sqrt(2)) * step_size, max_step_arc_length ]);
        fprintf('step_size increased to %.10g\n', step_size);
    end

    % Direction update (EMA smooth)
    direction_raw = p_current - p_prior;
    if norm(direction_raw) < 1e-14
        warning('Direction update nearly zero at step %d; keeping previous direction.', i);
        direction_raw = direction;
    end
    direction_raw = direction_raw / norm(direction_raw);

    direction = (1-dir_smooth_alpha)*direction + dir_smooth_alpha*direction_raw;
    direction = direction / norm(direction);
    direction_norm = [-direction(2); direction(1)];

    % Visualization
    if viz
        figure(19); clf;
        plot(p_history(1,:), p_history(2,:), 'ko-'); hold on;
        if ~isempty(p_history_all) && size(p_history_all,1) == 2 && ~isempty(metric_history_all)
            m = metric_history_all(:);
            n = size(p_history_all,2);
            nn = min(length(m), n);
            scatter(p_history_all(1,1:nn), p_history_all(2,1:nn), 30, m(1:nn), 'filled');
        end
        hold off;
        xlabel('\sigma'); ylabel('b'); title(sprintf('i=%d', i));
        drawnow; pause(0.05);
    end
end
end
