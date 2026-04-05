function features = feature_evaluation_target(featpar, modelpar)
% FEATURE_EVALUATION_TARGET
%
% First-pass feature for HALF-TARGET / source-defect continuation on [0,L].
%
% Input:
%   modelpar.ic0      : seed / warm-start initial condition (N x 3)
%   modelpar.a        : sigma
%   modelpar.b        : b
%   modelpar.T_relax  : relaxation time
%   modelpar.T_obs    : observation time
%
% Output:
%   scalar feature measuring coexistence of:
%     (1) spatial structure in the interior/core
%     (2) temporal oscillation in the far-right tail
%
% Suggested interpretation:
%   larger value = stronger half-target/source-defect structure
%   smaller value = more pure stripe / wave / non-target state
%
% Current scalar choice:
%   features = min(Sc, Tr)
% where
%   Sc = mean spatial variance in core window over time
%   Tr = mean temporal variance in tail window over space

    % ---------------- IC ----------------
    if isfield(modelpar,'ic0') && ~isempty(modelpar.ic0)
        ic0 = modelpar.ic0;
        else
        seed_file = 'half_target_seed.mat';
        if ~isfile(seed_file)
            error('feature_evaluation_target: cannot find target seed file %s.', seed_file);
        end

        S = load(seed_file);
        if isfield(S,'ic')
            ic0 = S.ic;
        elseif isfield(S,'ic_end')
            ic0 = S.ic_end;
        else
            fn = fieldnames(S);
            ic0 = S.(fn{1});
        end
    end

    % ---------------- parameters ----------------
    par = struct;
    par.sigma = modelpar.a;
    par.b     = modelpar.b;

    if isfield(modelpar,'T_relax') && ~isempty(modelpar.T_relax)
        T_relax = modelpar.T_relax;
    else
        T_relax = 240;
    end

    if isfield(modelpar,'T_obs') && ~isempty(modelpar.T_obs)
        T_obs = modelpar.T_obs;
    else
        T_obs = 240;
    end

    % ---------------- simulate: relax then observe ----------------
    [ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
    [~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

    if isempty(V)
        features = 0;
        return
    end

    % ---------------- windows ----------------
    % V is assumed [time x space]
    [nt, nx] = size(V);

    % use last 40% of time samples (same spirit as your defect_score)
    t0 = max(1, round(0.60 * nt));
    Vuse = V(t0:end, :);

    % core/interior window: center-ish portion of the half domain
    j_core1 = max(1, round(0.40 * nx));
    j_core2 = min(nx, round(0.60 * nx));
    Vcore = Vuse(:, j_core1:j_core2);

    % tail window: far right tail
    j_tail1 = max(1, round(0.80 * nx));
    j_tail2 = nx;
    Vtail = Vuse(:, j_tail1:j_tail2);

    % ---------------- feature pieces ----------------
    % Sc = mean spatial variance over time in the core
    % (measures stripe/Turing-like structure in the interior)
    Sc = mean(var(Vcore, 0, 2));

    % Tr = mean temporal variance over space in the tail
    % (measures Hopf-like oscillation in the tail)
    Tr = mean(var(Vtail, 0, 1));

    % scalar target feature: requires both to be present
    features = min(Sc, Tr);

    % optional normalization / squashing could be added later if needed
end

