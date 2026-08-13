function features = feature_evaluation_left(featpar, modelpar)
%FEATURE_EVALUATION_LEFT
% Left-side spiral feature:
%   rV_time = ||V - flipud(V)||_F / ||V||_F
% computed after relaxing from ic0 and observing V(t,x).
%
% Gated by temporal coherence.

features = NaN;   % always initialize output

if strcmp(featpar.feature,'spiral')

    % ---------------- IC: prefer warm-start if provided ----------------
    if isfield(modelpar,'ic0') && ~isempty(modelpar.ic0)
        ic0 = modelpar.ic0;
    else
        S = load('ic_0.45_10.mat');
        if isfield(S,'ic')
            ic0 = S.ic;
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
        T_relax = 360;
    end

    if isfield(modelpar,'T_obs') && ~isempty(modelpar.T_obs)
        T_obs = modelpar.T_obs;
    else
        T_obs = 1200;
    end

    if isfield(modelpar,'use_last_window') && ~isempty(modelpar.use_last_window)
        use_last_window = logical(modelpar.use_last_window);
    else
        use_last_window = true;
    end

    if isfield(modelpar,'last_frac') && ~isempty(modelpar.last_frac)
        last_frac = modelpar.last_frac;
    else
        last_frac = 0.25;
    end

    coh_thr = 0.11;
    if isfield(modelpar,'coh_thr') && ~isempty(modelpar.coh_thr)
        coh_thr = modelpar.coh_thr;
    end

    % ---------------- simulate: relax then observe ----------------
    [ic_rel, ~, ~] = solve_brusselator_1d(ic0, par, T_relax, 0);
    [~, ~, V]      = solve_brusselator_1d(ic_rel, par, T_obs, 0);

    % fallback if V missing
    if isempty(V)
        u_final = ic_rel(:,2);
        denom = norm(u_final);
        if denom < eps
            features = 0;
        else
            features = norm(u_final - flipud(u_final)) / denom;
        end
        return
    end

    % ---------------- regularity gate ----------------
    if ~check_temporal_coherence(V, coh_thr)
        if isfield(modelpar,'invalid_feature_value') && ~isempty(modelpar.invalid_feature_value)
            features = modelpar.invalid_feature_value;
        else
            % Preserve the historical side-branch behavior unless the
            % caller explicitly chooses how irregular states are encoded.
            features = 1e6;
        end
        return
    end

    % ---------------- windowing ----------------
    if use_last_window
        T = size(V,1);
        t0 = max(1, floor((1-last_frac)*T));
        Vuse = V(t0:end,:);
    else
        Vuse = V;
    end

    denom = norm(Vuse,'fro');
    if denom < eps
        features = 0;
    else
        % LEFT detector = time-flip metric
        features = norm(Vuse - flipud(Vuse), 'fro') / denom;
    end

    return
end

% ---------------- otherwise: keep your original logic ----------------
pattern = cell(1,featpar.N);
for i = 1:featpar.N
    pattern{i} = model(modelpar); %#ok<NASGU>
end

features = 0;
end
