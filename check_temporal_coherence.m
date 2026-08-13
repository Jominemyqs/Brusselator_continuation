function is_regular = check_temporal_coherence(V, thr_coh)

if nargin < 2 || isempty(thr_coh)
    thr_coh = 0.11;
end

r_coh = temporal_coherence_score(V);
is_regular = (r_coh < thr_coh);
end