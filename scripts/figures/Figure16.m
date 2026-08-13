function Figure16
% FIGURE16
%
% Compares sweep-derived transition points with a selected continuation
% segment for the target right-up branch.
%
% Here we compare the sweep-derived points near b = 10.2--10.6 with
% continuation points 10:30 from the computed branch.

clear; clc; close all;

%% ================= USER SETTINGS =================

% Sweep-derived points from horizontal sweeps
Psweep = [ ...
    0.507, 10.2;
    0.504, 10.3;
    0.501, 10.4;
    0.499, 10.5;
    0.497, 10.6
    ];

% Continuation file
mat_file = 'target_right_up.mat';

% Continuation point range to plot
cont_idx_start = 10;
cont_idx_end   = 30;

% Remove duplicated start if present
remove_duplicate_start = true;

% Save output
save_png = true;
png_name = figure_output_path('Figure16_target_sweep_vs_continuation.png');

save_pdf = false;
pdf_name = figure_output_path('Figure16_target_sweep_vs_continuation.pdf');

%% ================= LOAD CONTINUATION =================
if ~isfile(mat_file)
    error('Cannot find continuation file: %s', mat_file);
end

S = load(mat_file);

Pcont = [];
if isfield(S,'p_hist_up')
    Pcont = S.p_hist_up;
elseif isfield(S,'p_history')
    Pcont = S.p_history;
elseif isfield(S,'p_hist')
    Pcont = S.p_hist;
end

if isempty(Pcont)
    error('Could not find continuation history in %s', mat_file);
end

% remove duplicated start if present
if remove_duplicate_start && size(Pcont,2) >= 2
    if norm(Pcont(:,1) - Pcont(:,2)) < 1e-12
        Pcont = Pcont(:,2:end);
    end
end

% clamp index range
n_total = size(Pcont,2);
cont_idx_start = max(1, cont_idx_start);
cont_idx_end   = min(cont_idx_end, n_total);

if cont_idx_start > cont_idx_end
    error('Invalid continuation index range.');
end

Pcont_snip = Pcont(:, cont_idx_start:cont_idx_end).';

%% ================= PLOT =================
fig = figure('Color','w','Position',[100 100 760 560]);
hold on; box on; grid on;

% sweep-derived points
plot(Psweep(:,1), Psweep(:,2), 'ko', ...
    'MarkerSize', 8, ...
    'LineWidth', 1.4, ...
    'MarkerFaceColor', 'y', ...
    'DisplayName', 'Sweep-derived transition points');

plot(Psweep(:,1), Psweep(:,2), 'k--', ...
    'LineWidth', 1.0, ...
    'HandleVisibility', 'off');

% continuation segment
plot(Pcont_snip(:,1), Pcont_snip(:,2), 'b^-', ...
    'LineWidth', 2.0, ...
    'MarkerSize', 7, ...
    'DisplayName', sprintf('Continuation points'));

% annotate only selected continuation points
label_idx = [10 15 20 25 30];
for idx = label_idx
    if idx >= cont_idx_start && idx <= cont_idx_end
        j = idx - cont_idx_start + 1;
        text(Pcont_snip(j,1)+0.00015, Pcont_snip(j,2), sprintf('C%d', idx), ...
            'FontSize', 9, 'Color', 'b');
    end
end

xlabel('\sigma');
ylabel('b');
legend('Location','best');

% zoom window
all_sigma = [Psweep(:,1); Pcont_snip(:,1)];
all_b     = [Psweep(:,2); Pcont_snip(:,2)];

xpad = 0.002;
ypad = 0.03;

xlim([min(all_sigma)-xpad, max(all_sigma)+xpad]);
ylim([min(all_b)-ypad,     max(all_b)+ypad]);

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
disp('Sweep-derived points:');
disp(Psweep);

fprintf('Continuation points used: %d:%d\n', cont_idx_start, cont_idx_end);
disp(Pcont_snip);

end
