function r_coh = temporal_coherence_score(V)

% Use last part of spacetime
nT = size(V,1);
if nT < 2
    r_coh = Inf;
    return;
end

n_keep = min(80, nT);
Vlate = V(end-n_keep+1:end, :);

% compare adjacent frames
diffs = zeros(n_keep-1,1);

for k = 2:n_keep
    diffs(k-1) = norm(Vlate(k,:) - Vlate(k-1,:)) / (norm(Vlate(k,:)) + 1e-12);
end

r_coh = mean(diffs);
end