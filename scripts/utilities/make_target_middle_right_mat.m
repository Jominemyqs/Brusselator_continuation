function make_target_middle_right_mat
% MAKE_TARGET_MIDDLE_RIGHT_MAT
%
% Creates target_middle_right.mat from manually supplied anchor points,
% with linear interpolation to add one point every 0.005 in sigma.
%
% Output fields are styled like your other middle-branch .mat files.

clear; clc;

%% anchor points
anchors = [ ...
    0.420, 9.495;
    0.440, 9.485;
    0.460, 9.475;
    0.480, 9.455;
    0.500, 9.445;
    0.510, 9.445;
    0.512, 9.435;
    0.515, 9.435
];

%% sigma grid: every 0.005, plus final endpoint if needed
sigma_grid = [0.420:0.005:0.510, 0.512, 0.515];

%% interpolate b values
b_grid = interp1(anchors(:,1), anchors(:,2), sigma_grid, 'linear');

pts = [sigma_grid(:), b_grid(:)];

%% build struct
S = struct();

S.p_hist_mid = pts.';
S.p_all_mid = pts.';
S.metric_all_mid = nan(1, size(pts,1));

dp = diff(pts,1,1);
seglen = sqrt(sum(dp.^2,2));
S.L_hist_mid = [0; cumsum(seglen)].';

S.p0 = pts(1,:);
S.p1 = pts(2,:);

S.counts_mid = 1:(size(pts,1)+1);
S.note = 'Synthetic target middle-right curve created from anchor points by linear interpolation every 0.005 in sigma.';
S.sigma_mid = pts(1,1);
S.b_mid = pts(1,2);

%% save
save('target_middle_right.mat', '-struct', 'S');

disp('Points saved to target_middle_right.mat:');
disp(pts);

fprintf('Saved target_middle_right.mat\n');

end