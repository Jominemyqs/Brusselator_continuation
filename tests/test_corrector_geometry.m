function test_corrector_geometry
%TEST_CORRECTOR_GEOMETRY  Fast tests that do not run the Brusselator PDE.

featpar = struct('feature','target','alpha',0.01,'avoid_steady',0,'N',1);
modelpar = struct;
modelpar.model = 'Brusselator_1D';
modelpar.ic_base = 1;
modelpar.ic_det = 1;
modelpar.use_hysteresis = false;
modelpar.drop_smooth_alpha = 1;
modelpar.expand_factors = 1;
modelpar.eval_cap = 40;

p_pred = [1; 2];
p_c1 = p_pred + [0.1; 0.1];
p_c2 = p_pred - [0.1; 0.1];

% Default side corrector: b is fixed and sigma is swept.
side = modelpar;
side.r_thr = 5.15;
side.crossing_direction = 'up';
side.dsigma_local = 0.1;
side.local_halfwidth_steps = 3;
aux = struct('prev_crossing',1.15,'no_cross_count',0);
[p_side, ~, hist_side, ~, ~, ~, aux_side] = bisect_interval( ...
    side, featpar, @linear_feature, @objective_evaluation, ...
    p_c1, p_c2, [], 3, aux);

assert(abs(p_side(1) - 1.15) < 1e-12);
assert(abs(p_side(2) - 2.00) < 1e-12);
assert(all(abs(hist_side(2,:) - 2.00) < 1e-12));
assert(strcmp(aux_side.last_sweep_axis,'sigma'));

% Middle corrector: sigma is fixed and b is swept.
middle = modelpar;
middle.corrector_mode = 'vertical_middle';
middle.r_thr = 5.10;
middle.crossing_direction_b = 'up';
middle.db_local = 0.1;
middle.local_halfwidth_steps_b = 3;
middle.expand_factors_b = 1;
aux = struct('prev_crossing',2.05,'no_cross_count',0);
[p_middle, ~, hist_middle, ~, ~, ~, aux_middle] = bisect_interval( ...
    middle, featpar, @linear_feature, @objective_evaluation, ...
    p_c1, p_c2, [], 3, aux);

assert(abs(p_middle(1) - 1.00) < 1e-12);
assert(abs(p_middle(2) - 2.05) < 1e-12);
assert(all(abs(hist_middle(1,:) - 1.00) < 1e-12));
assert(strcmp(aux_middle.last_sweep_axis,'b'));

% Middle-right orientation: a decreasing feature must use a downward
% crossing in increasing b.
middle_right = middle;
middle_right.r_thr = 4.90;
middle_right.crossing_direction_b = 'down';
middle_right.use_hysteresis = true;
aux = struct('prev_crossing',2.05,'no_cross_count',0);
[p_middle_right, ~, ~, ~, ~, ~, aux_middle_right] = bisect_interval( ...
    middle_right, featpar, @decreasing_linear_feature, ...
    @objective_evaluation, p_c1, p_c2, [], 3, aux);

assert(abs(p_middle_right(1) - 1.00) < 1e-12);
assert(abs(p_middle_right(2) - 2.05) < 1e-12);
assert(contains(aux_middle_right.used_mode, 'threshold_crossing_b'));

% Integration smoke test: continuation advances in sigma while the
% vertical corrector keeps points close to sigma + 2*b = 5.1.
middle.db_local = 0.01;
middle.local_halfwidth_steps_b = 4;
middle.track_window_b = 0.08;
middle.drop_smooth_alpha = 1;

contpar = struct;
contpar.step_size = 0.05;
contpar.min_step_arc_length = 0.05;
contpar.max_step_arc_length = 0.05;
contpar.max_sim = 3;
contpar.L = Inf;
contpar.max_steps = 3;
contpar.b_bounds = [1.5, 2.5];
contpar.sigma_bounds = [0.5, 1.5];
contpar.dir_smooth_alpha = 0.15;
contpar.include_start = false;

start = struct;
start.point = [1; 2.05];
start.normal = [2; -1] / sqrt(5);

[P, counts] = continuation(contpar, featpar, middle, ...
    @linear_feature, @objective_evaluation, start, 0);

assert(size(P,2) == 3);
assert(all(diff(P(1,:)) > 0));
assert(all(counts > 0));
assert(max(abs(P(1,:) + 2*P(2,:) - middle.r_thr)) <= middle.db_local + 1e-12);

fprintf('test_corrector_geometry: PASS\n');
end


function r = linear_feature(~, modelpar)
r = modelpar.a + 2*modelpar.b;
end


function r = decreasing_linear_feature(~, modelpar)
r = 10 - modelpar.a - 2*modelpar.b;
end
