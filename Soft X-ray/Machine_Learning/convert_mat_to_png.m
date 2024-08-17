function convert_mat_to_png(input_dir,output_dir, z_grid, r_grid, name)
    % Create output directory if it doesn't exist
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    

    % Get list of all .mat files in the directory
    mat_files = dir(fullfile(input_dir, '*.mat'));
    % Check if there are any .mat files
    if isempty(mat_files)
        error('No .mat files found in the specified directory.');
    end

    [mesh_z,mesh_r] = meshgrid(z_grid,r_grid);

    % Loop over each .mat file and process it
    for k = 1:length(mat_files)
        % Load the .mat file
        mat_filename = fullfile(input_dir, mat_files(k).name);
        data = load(mat_filename);
        if isfield(data, name)
            
            X = data.(name);

            fig = figure('Visible', 'off');
            % ax = axes(fig);
            imagesc(X);%c=colorbar(ax, 'Ticks',[0,20,40]);
            axis off;
            %xlabel(ax, 'Z [Pixels]');ylabel(ax, 'R [Pixels]');

            colormap('turbo');
            %c.Label.String = 'Intensity [a.u.]';c.FontSize = 18;
            output_file = sprintf('image_%04d.png', k);
            exportgraphics(fig, fullfile(output_dir, output_file));
            close(gcf);
        else
            warning('Variable %s not found in file: %s', name, mat_files(k).n);
        end
    end
    
    disp('Conversion from .mat to .png completed.');
end