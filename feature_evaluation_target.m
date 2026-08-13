function features = feature_evaluation_target(featpar, modelpar)
% FEATURE_EVALUATION_TARGET
%
% Original target feature with added temporal-coherence gate.
%
% Original scalar choice:
%   features = min(Sc, Tr)
% where
%   Sc = mean spatial variance in core window over time
%   Tr = mean temporal variance in tail window over space
%
% Added gate:
%   if temporal coherence fails, force feature to the non-target side.

features = NaN;   % always initialize output

if strcmp(featpar.feature,'target')

    % ---------------- IC ----------------
    if isfield(modelpar,'ic0') && ~isempty(modelpar.ic0)
        ic0 = modelpar.ic0;
    else
        % fallback target seed
        seed_file = fullfile('out_half_target','ic_half_target_0.45_10.00.mat');
        if ~isfile(seed_file)
            seed_file = fullfile('out_half_target','half_target_seed.mat');
        end
        if ~isfile(seed_file)
            error('feature_evaluation_target: cannot find target seed file.');
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

    coh_thr = 0.11;
    if isfield(modelpar,'coh_thr') && ~isempty(modelpar.coh_thr)
        coh_thr = modelpar.coh_thr;
    end

    % ---------------- simulate: relax then observe ----------------
    [ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
    [~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

    if isempty(V)
        features = 0;
        return
    end

    % ---------------- regularity gate ----------------
    if ~check_temporal_coherence(V, coh_thr)
        if isfield(modelpar,'invalid_feature_value') && ~isempty(modelpar.invalid_feature_value)
            features = modelpar.invalid_feature_value;
        % Historical side-branch behavior follows below.
        elseif isfield(modelpar,'crossing_direction') && strcmp(modelpar.crossing_direction,'down')
            features = -1;    % right branches
        else
            features = 1e6;   % left branches
        end
        return
    end

    % ---------------- windows ----------------
    % V is [time x space]
    [nt, nx] = size(V);

    % use last 40% of time samples
    t0 = max(1, round(0.60 * nt));
    Vuse = V(t0:end, :);

    % core/interior window
    j_core1 = max(1, round(0.40 * nx));
    j_core2 = min(nx, round(0.60 * nx));
    Vcore = Vuse(:, j_core1:j_core2);

    % tail window
    j_tail1 = max(1, round(0.80 * nx));
    j_tail2 = nx;
    Vtail = Vuse(:, j_tail1:j_tail2);

    % ---------------- original feature pieces ----------------
    Sc = mean(var(Vcore, 0, 2));
    Tr = mean(var(Vtail, 0, 1));

    features = min(Sc, Tr);
    return
end

% ---------------- otherwise: keep original fallback ----------------
pattern = cell(1,featpar.N);
for i = 1:featpar.N
    pattern{i} = model(modelpar); %#ok<NASGU>
end

features = 0;
end
