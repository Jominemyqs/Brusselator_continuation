function append_points_to_right_up
% APPEND_POINTS_TO_RIGHT_UP
%
% Appends two manually chosen points to the end of right_up.mat.
%
% This is fine for plotting / figure construction, but note that the new
% points are manually added and are not original continuation output.

clear; clc;

%% user settings
infile  = 'right_up.mat';
outfile = 'right_up_extended.mat';

new_pts = [ ...
    0.3450, 13.7100;
    0.3350, 13.7400
    ];

%% load
if ~isfile(infile)
    error('Cannot find file: %s', infile);
end

S = load(infile);

% find history field
hist_name = '';
if isfield(S,'p_hist_up')
    hist_name = 'p_hist_up';
elseif isfield(S,'p_hist_dn')
    hist_name = 'p_hist_dn';
elseif isfield(S,'p_hist_mid')
    hist_name = 'p_hist_mid';
elseif isfield(S,'p_hist')
    hist_name = 'p_hist';
else
    error('Could not find continuation history field in %s.', infile);
end

P = S.(hist_name);

% make sure 2 x N
if size(P,1) ~= 2 && size(P,2) == 2
    P = P.';
end
if size(P,1) ~= 2
    error('Unexpected continuation array size: %dx%d', size(P,1), size(P,2));
end

%% append
P_new = [P, new_pts.'];

fprintf('Original number of points: %d\n', size(P,2));
fprintf('New number of points:      %d\n', size(P_new,2));
fprintf('Old last point: sigma = %.10f, b = %.10f\n', P(1,end), P(2,end));
fprintf('New last point: sigma = %.10f, b = %.10f\n', P_new(1,end), P_new(2,end));

S.(hist_name) = P_new;

% recompute matching L_hist if present
if strcmp(hist_name,'p_hist_up')
    Lname = 'L_hist_up';
elseif strcmp(hist_name,'p_hist_dn')
    Lname = 'L_hist_dn';
elseif strcmp(hist_name,'p_hist_mid')
    Lname = 'L_hist_mid';
else
    Lname = '';
end

if ~isempty(Lname)
    pts_all = P_new.';
    dp = diff(pts_all,1,1);
    seglen = sqrt(sum(dp.^2,2));
    S.(Lname) = [0; cumsum(seglen)].';
end

%% optional note
S.manual_append_note = 'Two final points were manually appended for plotting: (0.345,13.71), (0.335,13.74).';

%% save
save(outfile, '-struct', 'S');
fprintf('Saved extended file: %s\n', outfile);

end