function output_path = figure_output_path(filename)
%FIGURE_OUTPUT_PATH Return the standard location for generated figures.

script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
output_dir = fullfile(project_root, 'figures', 'generated');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

output_path = fullfile(output_dir, filename);
end
