function MAKE_FIGURE19_TARGET_THRESHOLD_FROM_REPRESENTATIVE_SWEEPS
% MAKE_FIGURE19_TARGET_THRESHOLD_FROM_REPRESENTATIVE_SWEEPS
%
% Three-panel version:
%   b = 9.98, 10.00, 10.02
%   no legend box
%   no sigma_vis text
%   no yellow marker

clear; clc; close all;

%% ================= USER INPUT =================
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
sweep_data_dir = fullfile(project_root, 'data', 'sweep_logs');
files = { ...
    fullfile(sweep_data_dir, 'sweep_log_98.csv'), ...
    fullfile(sweep_data_dir, 'sweep_log_100.csv'), ...
    fullfile(sweep_data_dir, 'sweep_log_102.csv') ...
    };

bvals = [9.98, 10.00, 10.02];

sigma_vis = [0.5110, 0.5150, 0.5100];

thr_list = [6, 8, 10];

save_png = true;
png_name = figure_output_path('Figure19_target_threshold_selection_3panels.png');

save_pdf = false;
pdf_name = figure_output_path('Figure19_target_threshold_selection_3panels.pdf');

%% ================= LOAD DATA =================
n = numel(files);
if numel(bvals) ~= n || numel(sigma_vis) ~= n
    error('files, bvals, and sigma_vis must have the same length.');
end

S = cell(1,n);

for i = 1:n
    if ~isfile(files{i})
        error('Cannot find file: %s', files{i});
    end

    T = readtable(files{i});

    if width(T) < 2
        error('File %s must contain at least two columns: sigma and feature.', files{i});
    end

    sigma = T{:,1};
    feat  = T{:,2};

    sigma = sigma(:);
    feat  = feat(:);

    [sigma, idx] = sort(sigma);
    feat = feat(idx);

    S{i}.sigma = sigma;
    S{i}.feat  = feat;
end

%% ================= PLOT =================
fig = figure('Color','w','Position',[100 100 1350 420]);
tl = tiledlayout(1,n,'TileSpacing','compact','Padding','compact');

for i = 1:n
    nexttile; hold on;

    sigma = S{i}.sigma;
    feat  = S{i}.feat;

    plot(sigma, feat, 'o-', 'LineWidth', 1.5, 'MarkerSize', 5);

    yline(thr_list(1), '--', '6',  'LineWidth', 1.1, 'LabelVerticalAlignment','bottom');
    yline(thr_list(2), '--', '8',  'LineWidth', 1.1, 'LabelVerticalAlignment','bottom');
    yline(thr_list(3), '--', '10', 'LineWidth', 1.1, 'LabelVerticalAlignment','bottom');

    xline(sigma_vis(i), ':', 'LineWidth', 1.5);

    grid on; box on;
    xlabel('\sigma');
    ylabel('target feature');
    title(sprintf('b = %.2f', bvals(i)));
end



%% ================= SAVE =================
if save_png
    exportgraphics(fig, png_name, 'Resolution', 300);
    fprintf('Saved PNG: %s\n', png_name);
end

if save_pdf
    exportgraphics(fig, pdf_name, 'ContentType','vector');
    fprintf('Saved PDF: %s\n', pdf_name);
end

end
