function Figure_last_points_from_mat
% FIGURE_LAST_POINTS_FROM_MAT
%
% Reads a continuation .mat file, extracts the continuation history,
% and plots the last portion of the continuation branch.
%
% Supports:
%   p_hist_up, p_hist_dn, p_hist_mid,
%   p_history, p_hist,
%   and falls back to p_all_up / p_all_dn / p_all_mid if needed.
%
% -------------------------------------------------------------------------

clear; clc; close all;

%% ================= USER SETTINGS =================
%mat_file = 'target_middle_left.mat';
mat_file = 'target_right_down.mat';
%mat_file = 'middle_left.mat';
%mat_file = 'right_down.mat';
%mat_file = 'right_up.mat';

n_last_points = 20;
remove_duplicate_start = true;
n_labels = 5;

save_png = true;
png_name = figure_output_path('Figure_last_points_middle_right.png');

save_pdf = false;
pdf_name = figure_output_path('Figure_last_points_middle_right.pdf');

%% ================= LOAD MAT FILE =================
if ~isfile(mat_file)
    error('Cannot find continuation file: %s', mat_file);
end

S = load(mat_file);

disp('Fields found in MAT file:');
disp(fieldnames(S));

Pcont = [];
source_name = '';

% Preferred: accepted continuation histories
candidate_hist_fields = { ...
    'p_hist_up', 'p_hist_dn', 'p_hist_mid', ...
    'p_history', 'p_hist'};

for k = 1:numel(candidate_hist_fields)
    fn = candidate_hist_fields{k};
    if isfield(S, fn)
        Pcont = S.(fn);
        source_name = fn;
        break;
    end
end

% Fallback: dense/all stored points
if isempty(Pcont)
    candidate_all_fields = {'p_all_up', 'p_all_dn', 'p_all_mid', 'p_all'};
    for k = 1:numel(candidate_all_fields)
        fn = candidate_all_fields{k};
        if isfield(S, fn)
            Pcont = S.(fn);
            source_name = fn;
            break;
        end
    end
end

if isempty(Pcont)
    error('Could not find a continuation array in %s', mat_file);
end

fprintf('Using field: %s\n', source_name);

% Make sure shape is 2 x N
if size(Pcont,1) ~= 2 && size(Pcont,2) == 2
    Pcont = Pcont.';
end

if size(Pcont,1) ~= 2
    error('Continuation array has unexpected size: %dx%d', size(Pcont,1), size(Pcont,2));
end

% Remove duplicated start if present
if remove_duplicate_start && size(Pcont,2) >= 2
    if norm(Pcont(:,1) - Pcont(:,2)) < 1e-12
        Pcont = Pcont(:,2:end);
    end
end

n_total = size(Pcont,2);

if n_total == 0
    error('Continuation history is empty.');
end

idx_start = max(1, n_total - n_last_points + 1);
idx_end   = n_total;

Ptail = Pcont(:, idx_start:idx_end).';

%% ================= PLOT =================
fig = figure('Color','w','Position',[100 100 760 560]);
hold on; box on; grid on;

plot(Pcont(1,:), Pcont(2,:), '-', ...
    'Color', [0.75 0.75 0.75], ...
    'LineWidth', 1.2, ...
    'DisplayName', sprintf('Full branch (%s)', source_name));

plot(Ptail(:,1), Ptail(:,2), 'b^-', ...
    'LineWidth', 2.0, ...
    'MarkerSize', 7, ...
    'DisplayName', sprintf('Last %d points', size(Ptail,1)));

plot(Ptail(end,1), Ptail(end,2), 'ro', ...
    'MarkerSize', 9, ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Final point');

label_start = max(idx_start, n_total - n_labels + 1);
for idx = label_start:n_total
    j = idx - idx_start + 1;
    text(Ptail(j,1)+0.00015, Ptail(j,2), sprintf('C%d', idx), ...
        'FontSize', 9, 'Color', 'b');
end

xlabel('\sigma');
ylabel('b');
title(sprintf('Last continuation points from %s', strrep(mat_file,'_','\_')));
legend('Location','best');

xpad = 0.002;
ypad = 0.05;
xlim([min(Ptail(:,1))-xpad, max(Ptail(:,1))+xpad]);
ylim([min(Ptail(:,2))-ypad, max(Ptail(:,2))+ypad]);

%% ================= SAVE =================
if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved PNG: %s\n', png_name);
end

if save_pdf
    exportgraphics(fig, pdf_name, 'ContentType', 'vector');
    fprintf('Saved PDF: %s\n', pdf_name);
end

%% ================= PRINT VALUES =================
fprintf('\nTotal number of points in %s: %d\n', source_name, n_total);
fprintf('Showing last points: %d:%d\n\n', idx_start, idx_end);

disp('Last continuation points [sigma, b]:');
disp(Ptail);

fprintf('Final continuation point:\n');
fprintf('sigma = %.10f, b = %.10f\n', Ptail(end,1), Ptail(end,2));

end
