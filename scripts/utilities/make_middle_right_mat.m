function make_middle_right_mat
% MAKE_MIDDLE_RIGHT_MAT
%
% Creates middle_right.mat from manually chosen points.
% Uses p_hist_mid so it behaves like your other middle-branch files.
%
% The curve here is synthetic / manually constructed for plotting.

clear; clc;

%% manual + interpolated points
pts = [ ...
    0.420, 9.515000;
    0.425, 9.513250;
    0.430, 9.511500;
    0.435, 9.509750;
    0.440, 9.508000;
    0.445, 9.504750;
    0.450, 9.501500;
    0.455, 9.498250;
    0.460, 9.495000;
    0.465, 9.492500;
    0.470, 9.490000;
    0.475, 9.487500;
    0.480, 9.485000;
    0.485, 9.480000;
    0.490, 9.475000;
    0.495, 9.470000;
    0.500, 9.465000;
    0.505, 9.440000;
    0.510, 9.415000;
    0.513, 9.418000;
    0.517, 9.455000
];

outfile = 'middle_right.mat';

%% build struct
S = struct();

% accepted continuation-style history: 2 x N
S.p_hist_mid = pts.';

% also store dense points the same way
S.p_all_mid = pts.';

% placeholder metric since only coordinates are known
S.metric_all_mid = nan(1, size(pts,1));

% simple cumulative arclength
dp = diff(pts,1,1);
seglen = sqrt(sum(dp.^2,2));
S.L_hist_mid = [0; cumsum(seglen)].';

% start points
S.p0 = pts(1,:);
S.p1 = pts(2,:);

% optional placeholders / metadata
S.counts_mid = 1:(size(pts,1)+1);
S.note = 'Synthetic middle-right curve created from manually supplied points with piecewise-linear interpolation in sigma.';
S.sigma_mid = pts(1,1);
S.b_mid = pts(1,2);

%% save
save(outfile, '-struct', 'S');
fprintf('Saved %s\n', outfile);

%% print
disp('Points saved to middle_right.mat:');
disp(pts);

end