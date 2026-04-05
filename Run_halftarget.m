%% Run_targethalf.m — make ONE half-target (source-defect) on [0,L]
% Produces a single wedge (target) on the half domain and saves artifacts.

clear; clc; close all; rng(1);

% ---------------- paths/output ----------------
outdir = 'out_half_target';
if ~exist(outdir,'dir'), mkdir(outdir); end

% ---------------- parameters (Björn’s target band) ----------------
par.sigma = 0.45;   % inside mixed Hopf–Turing band
par.b     = 10.00;  % important: b=10 for source-defect tails
T         = 240;    % integration time (raise to 300 if needed)
plot_flag = 1;      % 1 = space–time plot

% ---------------- load grid / initial condition ----------------
% Expect ic.txt as either:
%   (i) 3 cols: [x, u0, v0]  OR
%  (ii) 1 col : [x] only, in which case we construct u0,v0 below.
D = load('ic.txt');
if size(D,2) == 3
    x  = D(:,1); u0 = D(:,2); v0 = D(:,3);
else
    x  = D(:);
    % base steady state for Brusselator (De Wit params, consistent with solver)
    a  = 2.5;
    b0 = par.b;
    u0 = a * ones(size(x));
    v0 = (b0/a) * ones(size(x));

    % add a localized bump near the left boundary to create a half-target
    Lh = x(end) - x(1);
    w  = 0.04*Lh;                     % width of bump (0.02–0.06 works)
    amp = 0.6;                         % amplitude (tune 0.4–0.8 if needed)
    u0  = u0 + amp * exp(-((x - x(1)).^2)/(2*w^2));

    % tiny noise helps the Hopf tails
    u0 = u0 + 0.01*randn(size(x));
    v0 = v0 + 0.01*randn(size(x));
end

% enforce column vectors and assemble IC (N×3)
x  = x(:); u0 = u0(:); v0 = v0(:);
if ~(numel(x)==numel(u0) && numel(x)==numel(v0))
    error('Length mismatch: x=%d, u0=%d, v0=%d', numel(x), numel(u0), numel(v0));
end
ic = [x, u0, v0];

% ---------------- run solver (handle 2- or 3-output variants) ----------------
S = []; V = [];
try
    [ic_end, S, V] = solve_brusselator_1d(ic, par, T, plot_flag);
catch
    [ic_end, V]    = solve_brusselator_1d(ic, par, T, plot_flag);
end

% ---------------- save artifacts ----------------
png_name = sprintf('st_half_target_sigma-%.3f_b-%.2f.png', par.sigma, par.b);
if exist('S','var') && ~isempty(S)
    save(fullfile(outdir,'half_target_seed.mat'), 'ic_end','par','S','V');
else
    save(fullfile(outdir,'half_target_seed.mat'), 'ic_end','par','V');
end
exportgraphics(gcf, fullfile(outdir, png_name));

% ---------------- quick sanity check (optional) ----------------
if ~isempty(V)
    [Sc,Tails,ok] = defect_score(V);
    fprintf('Center spatial var=%.2e, Tail temporal var=%.2e → half-target=%d\n', Sc, Tails, ok);
end

fprintf('Saved seed + plot in %s\n', outdir);

% ===== helper: crude detector (keeps you honest) =====
function [Sc,Tails,tf] = defect_score(V)
    nt = size(V,1); nx = size(V,2); t0 = round(0.6*nt);
    Vc = V(t0:end, round(0.40*nx):round(0.60*nx));     % center (Turing core)
    Vr = V(t0:end, round(0.80*nx):nx);                 % right tail (Hopf)
    Sc    = mean(var(Vc,0,2));                         % spatial variance over time
    Tr    = mean(var(Vr,0,1));                         % temporal variance in tail
    Tails = Tr;
    tf = (Sc > 1e-3) && (Tails > 1e-3);
end

